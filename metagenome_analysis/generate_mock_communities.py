import sys
import os
import argparse
import random
import subprocess

def parse_arguments():

	"""Describe and parse input arguments"""
	descr='\nCreate mock communities based on specified assembly summary files and origin species genomes\n'
	parser=argparse.ArgumentParser(description=descr.replace("'", ''), formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--asm', help='Assembly summary file containing the genomes to create mock communities from', \
		required=True)
	parser.add_argument('--n_genomes', help='number of genomes to randomly sample for each species', default=1, type=int)
	parser.add_argument('--n_samples', help='number of mock communities to create', default=1, type=int)
	parser.add_argument('--out_dir', help='path to output directory', required=True)
	parser.add_argument('--genus_only', help='Generates n_samples communities per origin containing only genomes from the respective genus', action='store_true')
	parser.add_argument('--origins', help='path to file containing sampled origin genomes', required='--genus_only' not in sys.argv)

	args=parser.parse_args()
	return args

def download_genomes(accession_file, assembly_summary):

	urls=[l.split('\t')[19] for l in assembly_summary if l.split('\t')[0] in [x.rstrip('\n') for x in open(accession_file, 'r')]]

	if not os.path.exists(f'{os.path.dirname(accession_file)}/genomes'):
		os.mkdir(f'{os.path.dirname(accession_file)}/genomes')

	#download 
	for url in urls:
		download=f'wget -r -nH --cut-dirs=7 {url}/{os.path.basename(url)}_genomic.fna.gz -P {os.path.dirname(accession_file)}/genomes'
		subprocess.call(download, shell=True)

	#unzip
	unzip=f'gunzip {os.path.dirname(accession_file)}/genomes/*'
	subprocess.call(unzip, shell=True)

def genus_communities(assembly_summary):

	origins=['Providencia rettgeri', 'Providencia stuartii', 'Proteus terrae', 'Atlantibacter hermanii', 'Morganella morganii', 'Aeromonas media',\
	  	 'Citrobacter amalonaticus', 'Pseudochrobactrum asaccharolyticum', 'Kluyvera cryocrescens']
	exclude_specs=['Providencia thailandensis', 'genomosp', 'uncultured', 'unclassified']
	exclude_asms=['GCA_018138925.1', 'GCA_000981805.1', 'GCA_001559075.2', 'GCA_002591155.1', 'GCA_016647575.1', 'GCA_009931275.1', 'GCA_002749905.1']
	
	for o in origins:

		#filter assembly summary
		origin_genera=[l for l in assembly_summary if o.split(' ')[0] in l.split('\t')[7] and not ' '.join(l.split('\t')[7].split(' ')[0:2])==o]
		ori_gen_filtered=[l for l in origin_genera if not any(spec in ' '.join(l.split('\t')[7].split(' ')[0:2]) for spec in exclude_specs) and not\
						 any(acc in l.split('\t')[0] for acc in exclude_asms)]
		
		#Use assembly summary list with only origin genera (- those contained in the exclusion lists) to call sample_asms()
		genus_coms=sample_asms(ori_gen_filtered)
		
		write_output(genus_coms, assembly_summary, f'{o.replace(" ", "_").lower()}_')

def sample_asms(assembly_summary):

	"""Go through all species in the assembly summary file and sample a number of each.
	Repeat X times."""
	lines=[l for l in assembly_summary if not l.startswith('#') if not 'reference genome' in l \
	and not 'representative genome' in l and not l.split('\t')[7].split(' ')[1]=='sp.']

	#Create dict containing {species:[accession_numbers]}
	specs={}

	for l in lines:
		taxon=' '.join(l.split('\t')[7].split(' ')[:2])
		if not taxon in specs:
			specs[taxon]=[]
			specs[taxon].append(l.split('\t')[0])
		else:
			specs[taxon].append(l.split('\t')[0])
	
	#Now randomly sample n_samples x n_genomes per species
	communities=[[] for i in range(int(args.n_samples))]
	for i in range(int(args.n_samples)):
		for k, v in specs.items():
			if len(v)>=args.n_genomes:
				communities[i].extend(random.sample(v, args.n_genomes))
			else:
				communities[i].extend(random.sample(v, len(v)))
	return communities

def write_output(communities, assembly_summary, name=''):

	n='\n'
	if not os.path.exists(f'{args.out_dir.rstrip("/")}/communities'):
		os.mkdir(f'{args.out_dir.rstrip("/")}/communities')

	if args.genus_only==False:
		origin_file=args.origins
		olines=[l for l in open(origin_file)]
		with open(f'{origin_file}_unique_headers.fna', 'w') as of:
			for l in olines:
				if l.startswith('>'):
					seqid=l.split('|')[0]
					new_seqid=f'{seqid}X'
					of.write(f'{l.replace(seqid, new_seqid)}')
				else:
					of.write(l)

	for i, com in enumerate(communities):
		if not os.path.exists(f'{args.out_dir.rstrip("/")}/communities/{name}community{i}'):
			os.mkdir(f'{args.out_dir.rstrip("/")}/communities/{name}community{i}')

		accession_file=f'{args.out_dir.rstrip("/")}/communities/{name}community{i}/accessions.txt'
		with open(accession_file, 'w') as f:
			for c in com:
				f.write(f'{c}{n}')				 

		#Download genomes
		download_genomes(accession_file, assembly_summary)		

		concatenate=f'cat {args.out_dir.rstrip("/")}/communities/{name}community{i}/genomes/GC* > {args.out_dir.rstrip("/")}/communities/{name}community{i}/genomes/{name}community{i}_contigs.fna'	
		subprocess.call(concatenate, shell=True)
		
		#Rewrite contig headers into kraken2 format
		rewrite=f'python /home/stefan/postdoc/scripts/wgs_origins/kraken2_scripts/convert_headers_to_kraken2.py --fasta {args.out_dir.rstrip("/")}/communities/{name}community{i}/genomes/{name}community{i}_contigs.fna --asm {args.asm} --genomes all'
		subprocess.call(rewrite, shell=True)

		#Append origin file with unique headers to
		if args.genus_only==False:
			append=f'cat {origin_file}_unique_headers.fna >> {args.out_dir.rstrip("/")}/communities/{name}community{i}/genomes/{name}community{i}_contigs.fna.kraken2_headers.fna' 
			subprocess.call(append, shell=True)
		
		if not os.path.exists(f'{args.out_dir.rstrip("/")}/communities/{name}community{i}/reads'):
			os.mkdir(f'{args.out_dir.rstrip("/")}/communities/{name}community{i}/reads')

		#generate reads separately for origin genomes and closely related genomes to control coverage 
		generate_reads=f'art_illumina -ss HS25 -f 2 -i {args.out_dir.rstrip("/")}/communities/{name}community{i}/genomes/{name}community{i}_contigs.fna.kraken2_headers.fna -o {args.out_dir.rstrip("/")}/communities/{name}community{i}/reads/{name}community{i}_contigs_kraken2_headers -l 150 -m 350 --sdev 35 --paired -na'
		subprocess.call(generate_reads, shell=True)

def main():

	global args
	args=parse_arguments()
	assembly_summary=[l for l in open(args.asm, 'r') if not l.startswith('#')]

	if args.genus_only==False:
		communities=sample_asms(assembly_summary)
		write_output(communities, assembly_summary)
	
	else:
		genus_communities(assembly_summary)

if __name__=='__main__':
	main()				
