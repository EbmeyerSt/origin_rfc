import sys, os
import subprocess
import argparse
import math
from argparse import RawTextHelpFormatter
from pygenomeviz import GenomeViz
from matplotlib.patches import Patch

def parse_arguments():
	man_description='Create origin description figure from selected genome extracts.'
	parser=argparse.ArgumentParser(description=man_description.replace("'", ""), formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d', help='Path to directory containing the necessary files')
	parser.add_argument('-g', help='File containing genview accessions of selected genome extracts in the order that they shall be visualized in.')
	parser.add_argument('--annotations', help='Add annotations to genes in the plot', action='store_true')
	args=parser.parse_args()

	return args

def extract_metadata(path, file):

	#Extract sequences from '_contexts.fna' file
	seq_file=[f'{path}/{file}' for file in os.listdir(f'{path}') if file.endswith('_contexts.fna')]

	seqs={}
	for line in open(seq_file[0], 'r'):
		if line.startswith('>'):
			header=line
			seq=''
		else:
			seq+=line
			seqs[header]=seq

	seq_lens={key.split('__')[-2]:len(value) for key, value in seqs.items()}

	#Extract coding regions from 'visualization_meta.csv'
	lines=[line.split('\t') for line in open(f'{path}/annotation_meta.csv', 'r') if len(line.split('\t'))<2 or \
			len(line.split('\t'))>4]

	cds={}

	for i, line in enumerate(lines):

		if line[0].startswith('>'):

			id=line[0].split('__')[-3]
			cds[id]={}
			cds[id]['organism']=line[0].split('__')[-2]
			cds[id]['accession']=line[0].split('__')[-1].rstrip('\n')
			cds[id]['cds']={}
		else:

			cds[id]['cds'][i]={}
			cds[id]['cds'][i]['annotation']=line[0]
			cds[id]['cds'][i]['start']=int(line[1])
			cds[id]['cds'][i]['stop']=int(line[2])
			cds[id]['cds'][i]['seq']=line[4]
			cds[id]['seq_len']=seq_lens[id]

			#The target gene has no frame in the annotation file, but will always be in + orientation
			#therefore set to one
			if line[4].rstrip('\n')=='target':
				cds[id]['cds'][i]['frame']='1'
				cds[id]['cds'][i]['category']=line[4].rstrip('\n')

			else:
				if line[3]=='-':
					cds[id]['cds'][i]['frame']='-1'
				else:
					cds[id]['cds'][i]['frame']='1'

				try:
					cds[id]['cds'][i]['category']=line[5].rstrip('\n')
				except: 
					cds[id]['cds'][i]['category']='misc'
	
	#Extract info for specified genomes
	specs=[line.rstrip('\n') for line in open(f'{file}', 'r')]
	
	spec_cds={key:value for key, value in cds.items() if key in specs}
	spec_seqs={key:value for key, value in seqs.items() if key.split('__')[-2] in specs}

	return spec_cds, spec_seqs

def visualize(spec_cds, path, file, aln_dict):

	#Get order of genomes from accession input file and add alignments in that order
	genomes=[line.rstrip('\n') for line in open(file, 'r')]
	newline='\n'

	color_dict={'integron':'gold', 'transposase':'darkorchid', 'target':'firebrick', \
			'mobile':'limegreen', 'phage':'coral', 'hypothetical protein':\
			'darkgray', 'misc': 'cornflowerblue', 'hypothetical': 'darkgray',\
			'resistance':'salmon'}

	gv=GenomeViz(align_type='center', tick_style='bar', fig_track_height=1.1, feature_track_ratio=0.4, \
			tick_labelsize=10)

	for genome in genomes:

		#Add genome
		newline='\n'
		track=gv.add_feature_track(f'{spec_cds[genome]["organism"]}{newline}{spec_cds[genome]["accession"]}',\
				int(spec_cds[genome]["seq_len"]), labelsize=15)

		#Add cds and annotations
		for key2, value2 in spec_cds[genome]['cds'].items():
			if args.annotations==True:
				feature=track.add_feature(spec_cds[genome]["cds"][key2]['start'], spec_cds[genome]["cds"][key2]['stop'], \
						int(spec_cds[genome]["cds"][key2]['frame']), arrow_shaft_ratio=1.0, \
						label=spec_cds[genome]['cds'][key2]['annotation'], \
						facecolor=color_dict[spec_cds[genome]['cds'][key2]['category']], \
						labelvpos='top', labelhpos='center', labelrotation=20, linewidth=1, \
						labelsize=9)

			else:

				feature=track.add_feature(spec_cds[genome]["cds"][key2]['start'], spec_cds[genome]["cds"][key2]['stop'], \
						int(spec_cds[genome]["cds"][key2]['frame']), arrow_shaft_ratio=1.0, \
						facecolor=color_dict[spec_cds[genome]['cds'][key2]['category']], \
						 linewidth=1)

	#Get minimum % identity between alns
	ids=[]
	for key, value in aln_dict.items():
		for key2, value2 in aln_dict[key]['alns'].items():
			ids.append(float(aln_dict[key]['alns'][key2]['pident']))
	min_id=math.floor(min(ids))

	linkcolor='orange'
	#Add alignments
	for i, genome in enumerate(genomes):
		if i<len(genomes)-1:
			pair=f'{genomes[i]} {genomes[i+1]}'
			print(aln_dict[pair])
			print(pair)
			print(i)
			for key2, value2 in aln_dict[pair]['alns'].items():
				print(int(aln_dict[pair]['alns'][key2]['qstart']), int(aln_dict[pair]['alns'][key2]['qend']))
				link=gv.add_link((aln_dict[pair]['qname'],int(aln_dict[pair]['alns'][key2]['qstart']), \
						int(aln_dict[pair]['alns'][key2]['qend'])), (aln_dict[pair]['sname'],\
						int(aln_dict[pair]['alns'][key2]['sstart']), \
						int(aln_dict[pair]['alns'][key2]['send'])),size_ratio=1.0,\
						v=float(aln_dict[pair]['alns'][key2]['pident']), vmin=min_id, \
						normal_color=linkcolor)
	fig=gv.plotfig()

	handles=[
		Patch(color='gold', label='Integron associated'),
		Patch(color='darkorchid', label='Transposase'),
		Patch(color='firebrick', label='Target gene'),
		Patch(color='limegreen', label='Plasmid associated'),
		Patch(color='darkgray', label='Hypothetical protein'),
		Patch(color='cornflowerblue', label='Miscellaneous'),
		Patch(color='coral', label='Phage associated'),
		Patch(color='salmon', label='Resistance associated')
			]

	legend=fig.legend(handles=handles, loc='lower left')

	gv.set_colorbar(fig, vmin=min_id, bar_height=0.2, bar_colors=[linkcolor], \
			bar_label=f'% nucleotide{newline}identity', bar_labelsize=9, alpha=0.6)

	fig.savefig(f"{path}/origin_figure.png", dpi=300, pad_inches=0.1)

def annotate_cds(path, spec_cds):

	#Write cds to file
	with open(f'{path}/cds.fna', 'w') as f:
		for key, value in spec_cds.items():
			for key2, value2 in spec_cds[key]['cds'].items():
				if spec_cds[key]['cds'][key2]['category']!='target':
					f.write('>'+str(key2)+'\n'+spec_cds[key]['cds'][key2]['seq']+'\n')
	
	ncbi_prot='/storage/shared/databases/ncbi_protein/bacterial/ncbi_protein_bacterial_c95.dmnd'
	isfinder='/storage/shared/databases/ncbi_protein/bacterial/ISFinder.dmnd'
	#blast cds against ncbi protein database
	if not os.path.exists(f'{path}/cds_ncbi.csv'):
		prot_cmd=f'diamond blastx -p 40 -d {ncbi_prot} -q {path}/cds.fna -o {path}/cds_ncbi.csv --id 90 --more-sensitive --subject-cover 0.9 --max-target-seqs 1 -f 6 qseqid stitle pident'
		subprocess.call(prot_cmd, shell=True)

	#blast against ISFinder
	if not os.path.exists(f'{path}/cds_isfinder.csv'):
		isf_cmd=f'diamond blastx -p 40 -d {isfinder} -q {path}/cds.fna -o {path}/cds_isfinder.csv --id 90 --more-sensitive --subject-cover 0.9 --max-target-seqs 1 -f 6 qseqid stitle pident'
		subprocess.call(isf_cmd, shell=True)

	#Read results into dicts
	ncbi_dict={int(line.split('\t')[0]):' '.join(line.split('\t')[1].split(' ')[1:]).split(' [')[0] for line in open(f'{path}/cds_ncbi.csv', 'r')}
	isf_dict={int(line.split('\t')[0]):line.split('\t')[1] for line in open(f'{path}/cds_isfinder.csv', 'r')}

	#Now overwrite the old annotations with the new in spec_cds
	#Also update the categories based on the new annotations
	
	for key, value in spec_cds.items():
		for key2, value2 in spec_cds[key]['cds'].items():
			if key2 in ncbi_dict:
				spec_cds[key]['cds'][key2]['annotation']=ncbi_dict[key2].rstrip(' (plasmid)')

				if len(ncbi_dict[key2])>40 and '/' in ncbi_dict[key2]:
					spec_cds[key]['cds'][key2]['annotation']=ncbi_dict[key2].split('/')[0]

				if len(ncbi_dict[key2].split(' '))>3 and '(' in ' '.join(ncbi_dict[key2].split(' ')[-2:]):
					spec_cds[key]['cds'][key2]['annotation']=ncbi_dict[key2].split('(')[0]

			if key2 in isf_dict:
				spec_cds[key]['cds'][key2]['annotation']=isf_dict[key2]


			if key2 in ncbi_dict:
				if spec_cds[key]['cds'][key2]['category']=='resistance' and not any(i in ncbi_dict[key2]\
						.lower() for i in ['beta-lactam', 'aminoglycoside', \
						'carbapenem', 'tetracycline', 'macrolide', 'penicillin', 'quinolone',\
						'streptomycin', 'multidrug', 'tetr', 'stra', 'strb']):

					spec_cds[key]['cds'][key2]['category']='misc'

				if spec_cds[key]['cds'][key2]['category']=='transposase' and not any(i in ncbi_dict[key2]\
						.lower() for i in ['iscr', 'transpos', 'tnp', 'insertion']):

					spec_cds[key]['cds'][key2]['category']='misc'

			if 'hypothetical protein' in spec_cds[key]['cds'][key2]['annotation'].lower():
				spec_cds[key]['cds'][key2]['annotation']=''
				spec_cds[key]['cds'][key2]['category']='hypothetical protein'

			if any(i in spec_cds[key]['cds'][key2]['annotation'].lower() for i in ['beta-lactam', 'aminoglycoside', \
					'carbapenem', 'tetracycline', 'macrolide', 'penicillin', 'quinolon',\
					'streptomycin', 'multidrug', 'tetr', 'stra']):

				spec_cds[key]['cds'][key2]['category']='resistance'


			if  spec_cds[key]['cds'][key2]['category']=='hypothetical' and not 'hypothetical' in \
			ncbi_dict[key2]:

				spec_cds[key]['cds'][key2]['category']='misc'

			if any(keyword in spec_cds[key]['cds'][key2]['annotation'].lower() for keyword in \
					['integrase', 'inti', 'xerc', 'xerd']):
				spec_cds[key]['cds'][key2]['category']='integron'

			if any(keyword in spec_cds[key]['cds'][key2]['annotation'].lower() for keyword in \
					['mobiliza', 'moba', 'mobb', 'mobc', 'mobl', 'plasmid', 'relaxase', \
					'conjugation', 'secretion', 'conjugal transfer']):
				spec_cds[key]['cds'][key2]['category']='mobile'

			if any(keyword in spec_cds[key]['cds'][key2]['annotation'].lower() for keyword in \
					['phage']):
				spec_cds[key]['cds'][key2]['category']='phage'

			if any(keyword in spec_cds[key]['cds'][key2]['annotation'].lower() for keyword in \
					['iscr', 'transpos', 'tnp', 'insertion']) or \
					spec_cds[key]['cds'][key2]['annotation'].startswith('IS') or \
					spec_cds[key]['cds'][key2]['annotation'].startswith('Tn'):
				spec_cds[key]['cds'][key2]['category']='transposase'


	return spec_cds

def create_alignments(path, file, spec_seqs, spec_cds):

	#Write selected seqs to file
	with open(f'{path}/selected_seqs.fna', 'w') as f:
		for key, value in spec_seqs.items():
			f.write(key+value)

	#Create blastdb from file
	createdb=f'makeblastdb -in {path}/selected_seqs.fna -dbtype "nucl" -out {path}/selected_seqs.blastdb'
	subprocess.call(createdb,shell=True)

	#Blast file against itself
	blastn=f'blastn -query {path}/selected_seqs.fna -perc_identity 60 -db {path}/selected_seqs.blastdb -strand both -out {path}/selected_seqs.blastout -num_threads 48 -task blastn -outfmt "6 qseqid sseqid pident length qstart qend sstart send score"'
	subprocess.call(blastn, shell=True)

	#Filter alignments
	alns=[line for line in open(f'{path}/selected_seqs.blastout', 'r') if not line.split('\t')[0]\
			==line.split('\t')[1] and int(line.split('\t')[3])>100]

	aln_dict={}

	tab="\t"
	newline="\n"

	for i, line in enumerate(alns):
		if not f'{line.split(tab)[0].split("__")[-2]} {line.split(tab)[1].split("__")[-2]}' in aln_dict:
			key=f'{line.split(tab)[0].split("__")[-2]} {line.split(tab)[1].split("__")[-2]}'
			qid=line.split('\t')[0].split('__')[-2]
			sid=line.split('\t')[1].split('__')[-2]

			aln_dict[key]={}
			aln_dict[key]['qname']=f'{spec_cds[qid]["organism"]}{newline}{spec_cds[qid]["accession"]}'
			aln_dict[key]['sname']=f'{spec_cds[sid]["organism"]}{newline}{spec_cds[sid]["accession"]}'
			aln_dict[key]['alns']={}

		aln_dict[key]['alns'][i]={}
		aln_dict[key]['alns'][i]['qstart']=line.split('\t')[4]
		aln_dict[key]['alns'][i]['qend']=line.split('\t')[5]
		aln_dict[key]['alns'][i]['sstart']=line.split('\t')[6]
		aln_dict[key]['alns'][i]['send']=line.split('\t')[7]
		aln_dict[key]['alns'][i]['pident']=line.split('\t')[2]

	return aln_dict





def main():
	
	global args
	args=parse_arguments()
	path=args.d
	file=args.g

	spec_cds, spec_seqs=extract_metadata(path, file)
	
	aln_dict=create_alignments(path, file, spec_seqs, spec_cds)
	spec_cds=annotate_cds(path, spec_cds)
	visualize(spec_cds, path, file, aln_dict)

if __name__=='__main__':
	main()




