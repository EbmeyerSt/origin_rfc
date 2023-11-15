import argparse
import sys
import os
import random

def parse_args():
	
	"""Create and parse input arguments"""

	descr="""\nRewrite headers from ncbi assembly into kraken2 format.\n
		Extract subset of origin genomes to exclude from test database"""
	parser=argparse.ArgumentParser(description=descr.replace("'", ''), \
	formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--fasta', help='Path to input fasta file to rewrite headers for', required=True)
	parser.add_argument('--asm', help='Path to ncbi assembly summary file', \
	required=True)
	parser.add_argument('--sample', help='Randomly sample 1-3 assemblies per taxon to exclude', action='store_true')
	parser.add_argument('--genomes', help='State of genome assemblies to include', choices=['complete', 'reference', 'all'], default='all')
	args=parser.parse_args()

	return args

def select_subset(taxon_dict, target_dir):

	"""Create list of genomes to not include into kraken22 test db"""

	tab='\t'
	n='\n'
	#Write metadata to output file

	if args.sample==True:
		with open(f'{target_dir}/convert_headers_subsets_info.txt', 'w') as f:

			if args.genomes=='complete':
				f.write(f'Only complete genomes and reference genomes included in output.{n}')
			elif args.genomes=='reference':
				f.write(f'Reference and representative genomes included in output.{n}')
			elif args.genomes=='all':
				f.write(f'All genomes included in output.{n}')

			exclude=[]
			too_few=[]
			all_accs=[]
			for t in {' '.join(v.split(' ')[0:2]) for k, v in taxon_dict.items()}:
				accs=[k for k, v in taxon_dict.items() if v==t]
				f.write(f'{t}{tab}{len(accs)}{n}')
				
				if 10>len(accs)>1:
					exclude.extend(random.sample(accs, 1))
				elif len(accs)>=10:		
					exclude.extend(random.sample(accs, 3))
				elif len(accs)<=1:
					too_few.append(t)

				all_accs.extend(accs)
				print(f'{t} {len(accs)}')

			include=[a for a in all_accs if not a in exclude]
			print('Assemblies to exclude selected...')

			f.write(f'{n}Included: {len(include)}{n}')
			for i in include:
				f.write(f'{i}{n}')
			f.write(f'{n}')

			f.write(f'{n}Excluded: {len(exclude)}{n}')
			for e in exclude:
				f.write(f'{e}{n}')
			f.write(f'{n}')

			f.write(f'{n}Insufficient number of assemblies for testing:{n}')
			for fw in too_few:
				f.write(f'{fw}{n}')
	else:
		exclude='None'
		with open(f'{target_dir}/convert_headers_info.txt', 'w') as f:

			if args.genomes=='complete':
				f.write(f'Only complete genomes and reference genomes included in output.{n}')
			elif args.genomes=='reference':
				f.write(f'Reference and representative genomes included in output.{n}')
			elif args.genomes=='all':
				f.write(f'All genomes included in output.{n}')

			all_accs=[]
			for t in {' '.join(v.split(' ')[0:2]) for k, v in taxon_dict.items()}:
				accs=[k for k, v in taxon_dict.items() if v==t]
				all_accs.extend(accs)

				f.write(f'{t}{tab}{len(accs)}{n}')
			f.write(f'{n}')

			f.write(f'N genomes: {len(all_accs)}{n}')
			
			for a in all_accs:
				f.write(f'{a}{n}')

	print('Subsets selected...')

	return exclude

def map_contigs(target_dir, taxid_dict):

	"""Parse assembly files asnd map contig accessions to 
	assembly accession"""

	files=[f'{target_dir}/{f}' for f in os.listdir(target_dir) \
		if f.startswith('GCA') or f.startswith('GCF')]

	mapped={l:'_'.join(os.path.basename(f).split('_')[:2]) for f in files for l in open(f, 'r') \
		if l.startswith('>') and '_'.join(os.path.basename(f).split('_')[:2]) \
		in taxid_dict.keys()}

	print('Mapped contigs to accessions...')
	return mapped

def parse_summary(assembly_summary, target_dir):

	"""Parse assembly summary file to get dict with species:taxid.
	Exclude assemblies that do not have complete genome status
	due to possible contamination issues"""
	
	#filter assembly summary with accesions of files in directory - not all genomes
	#present in the assembly summary carry the respective gene, and those that do not carry it ar not in files!
	file_accs=['_'.join(os.path.basename(f).split('_')[0:2]) for f in os.listdir(target_dir) \
		if f.startswith('GCA') or f.startswith('GCF')]

	assemblies=[l for l in open(assembly_summary, 'r') if not l.startswith('#') and l.split('\t')[0] \
		in file_accs]
	
	#This list contains accessions of origin genomes/taxa to exclude because of unclear taxonomy
	specs_to_exclude=['Providencia thailandens', 'genomosp', 'uncultured', 'unclassified']
	asms_to_exclude=['GCA_018138925.1', 'GCA_000981805.1', 'GCA_001559075.2', 'GCA_002591155.1', 'GCA_016647575.1', 'GCA_009931275.1', 'GCA_002749905.1']

	asms=[l for l in assemblies if not any(x.lower() in l.split('\t')[7].lower() for x in specs_to_exclude)]
	asms=[l for l in asms if not any(y==l.split('\t')[0] for y in asms_to_exclude)]

	if args.genomes=='complete':
		lines=[l for l in asms if l.split('\t')[11]=='Complete Genome' or l.split('\t')[4]=='reference genome' or l.split('\t')[4]=='representative genome']

	elif args.genomes=='reference':
		lines=[l for l in asms if l.split('\t')[4]=='representative genome' \
			or l.split('\t')[4]=='reference genome']

	elif args.genomes=='all':
		lines=asms

	taxid_dict={line.split('\t')[0]:line.split('\t')[6] for line \
		 in lines}

	taxon_dict={line.split('\t')[0]:' '.join(line.split('\t')[7].split(' ')[:2]) for line \
		 in lines}

	print('Parsed assembly summary file for taxids...')
	return taxid_dict, taxon_dict

def yield_fasta(fasta):

	"""Generator function yielding one header and sequence at a time"""

	first=True
	for line in open(fasta, 'r'):
		if line.startswith('>'):
			if first==True:
				header=line
				seq=''
				first=False			

			else:
				yield header, seq
				header=line
				seq=''
		else:
			seq+=line
	yield header, seq

def rewrite_headers(fasta, taxid_dict, taxon_dict, mapped, target_dir):

	"""Read fasta file to dict and rewrite headers to kraken2
	   format"""
	n='\n'

	print('Rewriting headers...')

	if not args.sample==True:
		select_subset(taxon_dict, target_dir)
		with open(f'{fasta}.kraken2_headers.fna', 'w') as f:
			for i, (header, seq) in enumerate(yield_fasta(fasta)):
				if header in mapped.keys():
					f.write(f'>seq{str(i)}|kraken:taxid|{taxid_dict[mapped[header]]} {taxon_dict[mapped[header]]}{n}'+seq.replace(n, '')+n)	

	
	else:
		print('Filtering contigs...')
		exclude=select_subset(taxon_dict, target_dir)
		with open(f'{fasta}.kraken2_headers.fna', 'w') as f:
			for i, (header, seq) in enumerate(yield_fasta(fasta)):
				if header in mapped.keys():
					if not mapped[header] in exclude:
						f.write(f'>seq{str(i)}|kraken:taxid|{taxid_dict[mapped[header]]} {taxon_dict[mapped[header]]}{n}'+seq.replace(n, '')+n)	
			
		with open(f'{fasta}.kraken2_headers_excluded.fna', 'w') as f:
			for i, (header, seq) in enumerate(yield_fasta(fasta)):
				if header in mapped.keys():
					if mapped[header] in exclude:
						f.write(f'>seq{str(i)}|kraken:taxid|{taxid_dict[mapped[header]]} {taxon_dict[mapped[header]]}{n}'+seq.replace(n, '')+n)	


	print('All done!')

def main():
	
	global args
	args=parse_args()	

	assembly_summary=args.asm #Assembly summary file
	fasta=args.fasta #sequences to rename in fasta format
	target_dir=os.path.dirname(fasta) #Path to dir containing single assembly files

	taxid_dict, taxon_dict=parse_summary(assembly_summary, target_dir)
	mapped=map_contigs(target_dir, taxid_dict)
	rewrite_headers(fasta, taxid_dict, taxon_dict, mapped, target_dir)
	
if __name__=='__main__':
	main()
