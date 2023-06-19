#!/usr/local/env python3.7

import sys, os, subprocess, multiprocessing, seaborn, argparse, shutil, time
from argparse import RawTextHelpFormatter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Define arguments and parse them
def parse_arguments():

	man_description="%r\n\nIdentify ORFs, calculate gANI and visualize gANI between selected genomes\
			\n%r" % ('_'*80, '_'*80)
	parser=argparse.ArgumentParser(description=man_description.replace("'", ""), formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d', help='path to folder containing genomes as .fna file', required=True)
	parser.add_argument('-g', help='path to .txt file containing one genome accession per line', required=False)
	parser.add_argument('-o', help='path to output directory', required=True)
	parser.add_argument('-asm', help='path to Genbank assembly summary file', required=True)

	args=parser.parse_args()
	return args

def download_genomes():

	"""Download genomes of specified taxa"""

	#Read file with genome accessions
	to_download=[line.rstrip('\n').lower() for line in open(args.g, 'r')]
	lines=[line for line in open(args.asm, 'r') if not line.startswith('#')]
	urls=[line.split('\t')[19] for line in lines if line.split('\t')[0].lower() \
	in to_download or line.split('\t')[17].lower() in to_download]

	for url in urls:
		try:
			if not os.path.exists(args.d.rstrip('/')+'/'+os.path.basename(url+'_genomic.fna.gz')) \
				or not os.path.exists(args.d.rstrip('/')+'/'+'_'.join(os.path.basename(url)\
				.split('_')+'.fna')[0:2]):
				download=f'wget -r -nH --cut-dirs=7 {url+"/"+os.path.basename(url+"_genomic.fna.gz")} -P {args.d}'
				subprocess.call(download, shell=True)
				time.sleep(1)
	
		except Exception as e:
			print(f'EXCEPTION: {e}')

	#unzip
	unzip=f'gunzip {args.d.rstrip("/")+"/*.gz"}'
	subprocess.call(unzip, shell=True)

	#Rename
	files=[(args.d.rstrip('/')+'/'+file,args.d.rstrip("/")+"/"+'_'.join(file.split('_')[0:2])+'.fna') for file \
	in os.listdir(args.d) if file.endswith('_genomic.fna')]

	for file in files:
		rename=f'mv {file[0]} {file[1]}'
		subprocess.call(rename, shell=True)	

#Multiprocessing function
def multiprocess(function, processes, *items):
	
	"""Multiprocess input function"""

	#Create queue from files to process, add poison pill for each process
	queue=multiprocessing.Queue()

	for element in items[0]:
		queue.put(element)

	proc_list=[]
	p=processes

	for i in range(p):
		queue.put('STOP')

	#Start processes
	started_procs=0
	for i in range(p):

		started_procs+=1

		if len(items)==1:
			proc_list.append(multiprocessing.Process(target=function, args=(queue,)))

		elif len(items)>1:
			proc_list.append(multiprocessing.Process(target=function, args=(queue, items[1])))

		proc_list[-1].start()
	
	#Join finished processes
	finished_procs=0
	for p in proc_list:
		finished_procs+=1
		p.join()
		print(f'{finished_procs} of {started_procs} processes finished...')

def predict_ORFs(queue):

	"""Run prodigal on specified genomes to predict ORFs"""

	element=queue.get()
	
	if element=='STOP':
		return

	while True:

		prodigal=f'prodigal -i {element} -d {element.replace(".fna", "_orfs.fna")} -p single'
		if not os.path.exists(element.replace('.fna', '_orfs.fna')):
			subprocess.call(prodigal, shell=True)
		
		element=queue.get()
		if element=='STOP':
			return

def calculate_gani(queue):
	
	"""Calculate gANI based on ORF predictions"""

	#Here, element has to be a tuple - each genome combination must be present
	element=queue.get()
	if element=='STOP':
		return

	while True:
		
		#Create temporary output dir - needed to multiprocess ANIcalculator, otherwise the processes overwrite temporary output
		if not os.path.exists(element[0]+"_vs_"+os.path.basename(element[1])+".gani"):
			tmp_dir=element[0]+"_vs_"+os.path.basename(element[1])+"_tmp"
			os.mkdir(tmp_dir)
			gani=f'/storage/stefan/downloads/ANIcalculator_v1/ANIcalculator -genome1fna {element[0]} \
			-genome2fna {element[1]} -outfile {element[0]+"_vs_"+os.path.basename(element[1])+".gani"} -outdir {tmp_dir} \
			-ignoreList /storage/stefan/downloads/ANIcalculator_v1/ANIignore'
			if not os.path.exists(element[0]+"_vs_"+os.path.basename(element[1])+".gani"):
				print(gani)
				subprocess.call(gani, shell=True)

			#Remove tmp output dir
			shutil.rmtree(tmp_dir)

		#Get next element from queue
		element=queue.get()
		if element=='STOP':
			return

def create_ani_dict():

	"""Create dictionary containing gANIs for all combinations of assemblies """

	#create file containing all ani comparisons
	concatenate=f'cat {args.o.rstrip("/")+"/"+"*.gani"} > {args.o.rstrip("/")+"/gani_all.txt"}'
	subprocess.call(concatenate, shell=True)

	all_anis={line for line in open(args.o.rstrip("/")+"/gani_all.txt", 'r') if not line.startswith('GENOME1\tGENOME2')}

	#Create dict with accession and species names
	spec_dict={line.split('\t')[0]:line.split('\t')[7] for line in open(args.asm, 'r') if not line.startswith('#')}

	#Now create dict from concatenated file
	ani_dict={}
	asms1=[line.split('\t')[0] for line in all_anis]
	asms2=[line.split('\t')[1] for line in all_anis]
	asms1.extend(asms2)
	asms=set(asms1)

	for line in asms:
		split_line=line.split('\t')
		acc='_'.join(split_line[0].split('_')[:2])
		ani_dict[acc]={}
		ani_dict[acc]['organism']=' '.join(spec_dict[acc].split(' ')[0:2])
		ani_dict[acc]['anis']={}
		ani_dict[acc]['af']={}

		#add ani and af to self as 100%
		ani_dict[acc]['anis'][acc]=int(100)
		ani_dict[acc]['af'][acc]=int(100)

	for line in all_anis:
		split_line=line.split('\t')
		acc='_'.join(split_line[0].split('_')[:2])
		acc2='_'.join(split_line[1].split('_')[:2])
		ani_dict[acc]['anis'][acc2]=int(float(split_line[2]))
		ani_dict[acc]['af'][acc2]=int(float(split_line[4]))

	return ani_dict
	


def create_ani_matrix(ani_dict):
	
	"""Create gANI matrix from genomes"""

	#Create list of all present accessions and species as tuples - those become column names in df
	#(after sorting)
	column_list_spec=[(key, ani_dict[key]['organism']) for key, value in ani_dict.items()]

	#Sort list of tuples by species name
	column_list_spec.sort(key=lambda x: x[1])

	#CReate list containing only the accessions
	column_list=[element[0]+' '+element[1] for element in column_list_spec]

	#Now create row dictionary
	row_dict={}
	for genome in column_list:
		ani_list=[]
		for genome2 in column_list:
			try:
				ani_list.append(ani_dict[genome.split(' ')[0]]['anis'][genome2.split(' ')[0]])
			except: 
				ani_list.append(int(0))

		row_dict[genome]=ani_list

	#Now create matrix from column list and row_dict
	matrix=pd.DataFrame.from_dict(row_dict, orient='index', columns=column_list)
	matrix=matrix.astype(float)
	print(matrix)
	#plot as heatmap

	mask=np.triu(np.ones_like(matrix, dtype=bool))
	mask=mask[1:, :-1]
	matrix_2=matrix.iloc[1:,:-1].copy()
	
	#Set figure parameters according to number of genomes
	if len(column_list) in range(10, 30):
		plt.subplots(figsize=(25,20))
	#	hm=seaborn.heatmap(matrix_2, cmap='RdYlGn', annot=True, fmt='.1f', mask=mask)
		hm=seaborn.clustermap(matrix_2, annot=True, cmap='RdYlGn', annot_kws={'fontsize':7}, fmt='g')
		hm.fig.subplots_adjust(right=0.65)
		hm.ax_cbar.set_position((0.9, 0.06, .015, .15))
	#	hm.set_xticklabels(hm.get_xticklabels(), rotation=90)

	elif len(column_list) in range(2,10):
		plt.subplots(figsize=(5,4))
		seaborn.set(font_scale=0.5)
		hm=seaborn.heatmap(matrix_2, annot=True, cmap='RdYlGn', annot_kws={'fontsize':7}, fmt='g')
		hm.fig.subplots_adjust(right=0.65)
		hm.ax_cbar.set_position((0.9, 0.06, .015, .15))
		hm.set_yticklabels(hm.get_yticklabels(), fontsize=7)

	plt.tight_layout()
	plt.savefig(args.o.rstrip('/')+'/gani_heatmap.svg')

def main():

	try:
		download_genomes()
	
	except Exception as e:
		print(f'EXCEPTION:{e}')
		print('Probably -g was not specified, try again if you have a list of accessions!')

	#Collect genome files and predict ORFs
	genomes=[args.d.rstrip('/')+'/'+file for file in os.listdir(args.d) if file.endswith('.fna') \
	and not file.endswith('_orfs.fna')]
	multiprocess(predict_ORFs, 30, genomes)
	print('Orfs predicted...')

	#Collect ORF files and calculate gANI
	orfs=[args.o.rstrip('/')+'/'+file for file in os.listdir(args.o) if file.endswith('_orfs.fna')]

	#Create list that contains all combination of files as tuple
	orf_comb=[(element_1, element_2) for element_2 in orfs for element_1 in orfs if not element_1==element_2]

	#Run
	multiprocess(calculate_gani,multiprocessing.cpu_count(),orf_comb)
	print('Pairwise ganis calculated...')

	#Create ANI dictionary
	gani_dict=create_ani_dict()

	#Create pandas dataframe from gani_dict and plot as heatmap
	create_ani_matrix(gani_dict)
	print('gani heatmap created...')
	#cleanup()

def cleanup():

	fnas=[args.d.rstrip('/')+'/'+file for file in os.listdir(args.d) if file.endswith('.fna')]

	for file in fnas:
		os.remove(file)

	print('Assembly files removed, finshed!')


if __name__=='__main__':

	args=parse_arguments()
	main()
