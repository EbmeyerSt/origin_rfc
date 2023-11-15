import sys, os
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import ete3
import sqlite3

""" Create phylogenies of marker genes from assemblies"""

def parse_args():

	"""Create input arguemnts for the script"""

	man_description='Create phylogenies of marker genes from assemblies.'
	parser=argparse.ArgumentParser(description=man_description.replace("'", ""), formatter_class=RawTextHelpFormatter)
	parser.add_argument('-g', help='Path to folder containing genomes', required=True)
	parser.add_argument('-m', help='Path to file containing marker genes', required=True)
	parser.add_argument('-gv', help='Path to folder containing genview_analysis folders for resistance genes', required=True)
	parser.add_argument('-id', help='minimum id of genes in GV database to reference sequence', default=80)
	parser.add_argument('--rgenes', help='resistance genes to look for',  nargs='+', required=True)
	parser.add_argument('-o', help='output directory', required=True)
	parser.add_argument('--taxon_headers', help='Display taxon on phylogeny label', action='store_true')
	parser.add_argument('--visualize', help='Use ete3 to visualize tree and generate figure', action='store_true')
	parser.add_argument('--long_seqs_only', help='Show ARGs only if contexts are longer than 15kbp', action='store_true')
	parser.add_argument('--color_gradient', help='Color node background based on perc identity towards reference ARG', action='store_true')
	args=parser.parse_args()

	return args

def fasta_to_dict(file, modify_header=0):

	"""Read fasta file into dictionary"""
	
	seqs={}

	if modify_header==0:
		for line in open(file, 'r'):
			if line.startswith('>'):
				header=line
				seq=''
			else:
				seq+=line
				seqs[header]=seq
	else:

		for line in open(file, 'r'):
			if line.startswith('>'):
				header=line.rstrip('\n')+'__'+'_'.join\
				(os.path.basename(file).split('_')[0:2])+'\n'
				seq=''
			else:
				seq+=line
				seqs[header]=seq
	return seqs


def rewrite_asms():

	"""Append assembly accession to every contig in the assembly
	and write to file"""


	files=[f'{args.g}{file}' for file in os.listdir(args.g) if file.endswith\
		('_genomic.fna')]

	#Read taxonomic classification to file
	taxa={line.split('\t')[0]:' '.join(line.split('\t')[7].split(' ')[0:2]) for line in open(f'{args.gv}assembly_summary.txt', 'r') \
	if not line.startswith('#')}

	#Read to dict
	seqs={}

	for file in files:
		seqs.update(fasta_to_dict(file, modify_header=1))
	
	with open(f'{args.o}all_asms.fna', 'w') as f:
		for key, value in seqs.items():
			f.write(key+value)

	all_asms=f'{args.o}all_asms.fna'
	return (all_asms, taxa)

def find_markers(genomes_fasta, marker_fasta):

	""" Use diamond to identify marker genes in the respective genomes"""

	markers=[line.split(' ')[1].rstrip('\n') for line in open(marker_fasta, 'r') if \
			line.startswith('>')]

	#Create output directory if not exists
	if not os.path.exists(args.o):
		os.makedir(args.o)

	#Create db
	mkdb=f"diamond makedb --in {marker_fasta} -d {marker_fasta}.dmnd"
	subprocess.call(mkdb, shell=True)

	#blastx
	blastx=f"diamond blastx -p 40 -d {marker_fasta}.dmnd -q {genomes_fasta} --id 70 --subject-cover 70 --max-target-seqs 1 --ultra-sensitive -o {args.o}{os.path.basename(genomes_fasta)}.blastout.csv -f 6 qtitle qseqid stitle pident length qseq"
	subprocess.call(blastx, shell=True)

	outfile=f'{args.o}{os.path.basename(genomes_fasta)}.blastout.csv'
	return (markers, outfile)

def find_args(markers, outfile, taxa):

	"""Extract 1 hit per assembly. Check whether any arg was identified
	in the assembly, if yes append to header. Write seq headers as
	ASM_rgene1_rgenex"""

	#If genomes are downloaded manually, there is a risk that newly uploaded genomes, that are not part of the
	#genbank database, are included in the analysis, which will throw an error. Therefore, parse the blast
	#output file and remove all assemblies not found in the genview assembly summary

	#Get list of assembly accessions from genview assembly summary file
	gv_asms={line.split('\t')[0] for line in open(f'{args.gv}assembly_summary.txt', 'r')}

	#parse lines from outfile and discared the ones that come from asms not in the gv assembly_summary
	keep_lines=[line for line in open(outfile, 'r') if line.split('\t')[0].split('__')[1] in gv_asms]

	seqs_dict={}
	new_seqs={}
	new_keys={}
	marker_files_anno=[]

	#Read sequences for marker into dict
	#Later: Add statistics and counts to print
	for marker in markers:
		if not args.taxon_headers==True:
			seqs_dict[marker]={line.split('__')[-1].split('\t')[0]:line.split('\t')[-1] \
					for line in keep_lines if marker in line}
		else:
			seqs_dict[marker]={line.split('__')[-1].split('\t')[0]+'_'+taxa[line.split('__')[-1].split('\t')[0]].replace(' ', '_'):line.split('\t')[-1] \
					for line in keep_lines if marker in line}

		new_seqs[marker]={}
		new_keys[marker]={}
		nl='\n'
		#Go to arg folder, and extract information from annotation file
		for r in args.rgenes:

			if args.long_seqs_only==True:
				contexts=fasta_to_dict(f'{args.gv}{r}_{args.id}_analysis/{r}_contexts.fna')
				long_contexts=[key.split('__')[-2] for key, value in contexts.items() if \
				len(value)>=15000]

				positives={line.split('__')[-1].rstrip('\n') for line in open(\
				f'{args.gv}{r}_{args.id}_analysis/annotation_meta.csv', 'r') if line.startswith('>')\
				and line.split('__')[-3] in long_contexts}

			else:
					
				positives={line.split('__')[-1].rstrip('\n') for line in open(\
				f'{args.gv}{r}_{args.id}_analysis/annotation_meta.csv', 'r') if line.startswith('>')}

			for key, value in seqs_dict[marker].items():
				accession='_'.join(key.split('_')[0:2]).rstrip('\n')

				if accession in positives and not accession \
				in new_keys[marker]:

					new_key=f'{key.rstrip(nl)}|{r}'
					new_keys[marker][accession]=new_key			

				elif accession in positives and accession in \
				new_keys[marker]:
					new_key=f'{new_keys[marker][accession]}|{r}'
					new_keys[marker][accession]=new_key			

		for key, value in seqs_dict[marker].items():

			accession='_'.join(key.split('_')[0:2]).rstrip('\n')
			if not accession in new_keys[marker]:
					new_keys[marker][accession]=key


		#Write annotated headers and respective marker gene sequences to file
		with open(f'{args.o}{marker}_seqs_annotated.fna', 'w') as f:
			for key, value in seqs_dict[marker].items():
				accession='_'.join(key.split('_')[0:2]).rstrip('\n')
				f.write('>'+new_keys[marker][accession]+'\n'+value.rstrip('\n')+'\n')
		marker_files_anno.append(f'{args.o}{marker}_seqs_annotated.fna')
	
	return marker_files_anno

def align(marker_files_anno):

	"""Create multiple sequence alignment and phylogeny"""

	for file in marker_files_anno:
		mafft=f'mafft --auto --reorder --thread 40 {file} > {file}.aln'
		subprocess.call(mafft, shell=True)

		phylo=f'FastTree -gtr -nt < {file}.aln > {file}.tree'
		subprocess.call(phylo, shell=True)

	return f'{file}.tree'

def visualize(tree):

	"""Use ete3 to create a circular tree with branch label backgrounds colored by resistance gene content"""

	t = ete3.Tree(tree)

	colors=['steelblue', 'darkcyan', 'lightskyblue', 'dodgerblue', 'cornflowerblue', 'slateblue', 'mediumslateblue', 'plum', \
			'thistle', 'lightsteelblue', 'mediumpurple', 'violet']

	if args.color_gradient==True:

		#Connect to sqlite database
		connection=sqlite3.connect(f'{args.gv}genview_database.db')
		cursor=connection.cursor()

		query=f"""
		SELECT args.id, args.perc_id, genomes.assembly FROM args \
		INNER JOIN genomes ON genomes.id=args.genome_id \
		WHERE LOWER(arg_name)='{args.rgenes[0].lower()}' \
		"""

		print('Querying genview database...')
		cursor.execute(query)
		results=cursor.fetchall()

		#For isolates with >1 copy of the respective gene, show the one with lowest id towards
		#the reference (Since the reference should be the mobile gene and we want to compare the 
		#putatively chromosomal ones)

		perc_ids={}

		for result in results:
			if result[2] in perc_ids:
				if float(result[1])<float(perc_ids[result[2]]):
					perc_ids[result[2]]=result[1]
			else:
				perc_ids[result[2]]=result[1]

		#Create list to check perc_id ranges
		ranges=[]

		#Go through leaves and color according to percent identity	
		for leaf in t.iter_leaves():
			asm_acc='_'.join(leaf.name.split('_')[0:2])
		
			if '|' in leaf.name and asm_acc in perc_ids:
				if 70<float(perc_ids[asm_acc])<=80:
					leaf.img_style['bgcolor']='gold'
					ranges.append('70-80')
				elif 80<float(perc_ids[asm_acc])<=90:
					leaf.img_style['bgcolor']='yellowgreen'
					ranges.append('80-90')
				elif 90<float(perc_ids[asm_acc])<=100:
					leaf.img_style['bgcolor']='forestgreen'
					ranges.append('90-100')
				elif 60<float(perc_ids[asm_acc])<=70:
					leaf.img_style['bgcolor']='khaki'
					ranges.append('60-70')
				elif 50<float(perc_ids[asm_acc])<=60:
					leaf.img_style['bgcolor']='peachpuff'
					ranges.append('50-60')

			leaf.img_style['size']=0			

	else:
		combs=[]
		#Get list of unique resistance gene combinations
		for leaf in t.iter_leaves():
			if '|' in leaf.name:
				combs.append('|'.join(leaf.name.split('|')[1:]))

		#Assign unique combinations to colors
		leaf_colors={}
		for i, comb in enumerate(set(combs)):
			leaf_colors[comb]=colors[i]

		#Set background colors based on gene combinations
		for leaf in t.iter_leaves():
			if '|' in leaf.name:
				leaf.img_style['bgcolor']=leaf_colors['|'.join(leaf.name.split('|')[1:])]
	 
			leaf.img_style['size']=0			
		
	#Set all branch lengths equal	
	for node in t.traverse():
		node.img_style['hz_line_width']=3
		node.img_style['vt_line_width']=3
		node.img_style['size']=0
		node.dist=4

	ts=ete3.TreeStyle()
	ts.mode = 'c'
	
	#Remove blue circles from nodes
	t.convert_to_ultrametric()

	#Add legend
	color_schemes={
			'50-60':'peachpuff',
			'60-70':'khaki',
			'70-80':'gold',
			'80-90':'yellowgreen',
			'90-100':'forestgreen'
			}

	color_legend={key:value for key, value in color_schemes.items() if key in ranges}

	legend_face=ete3.TextFace('% nucleotide identity to reference:', fsize=30)
	legend_face.margin_right=2
	legend_face.margin_left=40
	legend_face.margin_top=10
	legend_face.margin_bottom=5
	ts.legend.add_face(legend_face, column=0)

	for i, (attribute, color) in enumerate(color_legend.items(), start=1):

		attribute_face = ete3.TextFace(attribute, fgcolor='black', fsize=30)
		attribute_face.margin_left=10
		attribute_face.margin_top=9
		color_face = ete3.CircleFace(28, color)
		color_face.margin_left=16
		color_face.margin_bottom=5
		ts.legend.add_face(attribute_face, column=i)
		ts.legend.add_face(color_face, column=i)

	ts.legend_position=3
	ts.legend.scale=1.5
	ts.show_scale=False

	t.render(f'{tree.replace(".fna.tree", ".tree.svg")}', tree_style=ts)
	t.render(f'{tree.replace(".fna.tree", ".tree.png")}', tree_style=ts)

def main():

	global args
	args=parse_args()
	
	if len(args.rgenes)>1 and args.color_gradient==True:
		print('Color gradient works only with a single gene!')
		sys.exit()

	#if tree file already exists, just redraw the tree
	if not os.path.exists(f'{args.o}rpoB_seqs_annotated.tree.svg'):
		all_asms, taxa=rewrite_asms()
		markers, outfile = find_markers(all_asms, args.m)
		marker_files_anno=find_args(markers, outfile, taxa)
		tree=align(marker_files_anno)

	else: 
		tree=f'{args.o}rpoB_seqs_annotated.fna.tree'

	visualize(tree)

if __name__=='__main__':
	main()




