import sys
import os
import argparse
import subprocess
import pandas as pd

""" Run specified mock communities against specified kraken2 database and extract classification metrics
    for new origin species"""

def parse_args():

	"""Describe and parse input arguments"""
	descr='\nRun metagenome against kraken2 database and extract classification metrics for origin species\n'
	parser=argparse.ArgumentParser(description=descr.replace("'", ''), formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--kraken_db', help='path to kraken2 database for read classification', required=True)
	parser.add_argument('--out_dir', help='path to output directory', required=True)
	parser.add_argument('--threads', help='number of threads to run kraken2 on', default=20)
	parser.add_argument('--analyse_only', help='Skip running kraken2, only analyse output files', action='store_true')
	#The assembly_summary file should only contain the origin species, since we do not care about the rest
	parser.add_argument('--asms', help='assembly summary file to map taxid to species', required=True)	
	parser.add_argument('--analyse_multiple', help='path to directory containing multiple communities to analyse')	
	parser.add_argument('--reads', help='path to read files (forward and reverse) in fasta or fastq format', required='--analyse_multiple' not in sys.argv, nargs='+')

	args=parser.parse_args()

	return args	

def run_kraken2():

	"""Run kraken2 on specified files"""

	if args.analyse_only==False and args.analyse_multiple==False:
		kraken2=f'kraken2 --db {args.kraken_db} --paired {args.reads[0]} {args.reads[1]} --output {args.out_dir.rstrip("/")}/kraken2_out.txt --confidence 0.3 --threads {args.threads} --report {args.out_dir.rstrip("/")}/kraken2_rep.txt --report-zero-counts --use-names'
		subprocess.call(kraken2, shell=True)
		outfile=f'{args.out_dir.rstrip("/")}/kraken2_rep.txt'

	elif args.analyse_multiple!=False and not args.analyse_only==True:

		#go through all mock community directories
		dirs=[f'{args.analyse_multiple.rstrip("/")}/{d}' for d in os.listdir(args.analyse_multiple)]
		
		for d in dirs:

			r1=[f'{d.rstrip("/")}/reads/{file}' for file in os.listdir(f'{d.rstrip("/")}/reads') if file.endswith('_headers1.fq')]
			r2=[f'{d.rstrip("/")}/reads/{file}' for file in os.listdir(f'{d.rstrip("/")}/reads') if file.endswith('_headers2.fq')]

			kraken2=f'kraken2 --db {args.kraken_db} --paired {r1[0]} {r2[0]} --output {d}/kraken2_out.txt --confidence 0.3 --threads {args.threads} --report-zero-counts --report {d}/kraken2_rep.txt --use-names'
			subprocess.call(kraken2, shell=True)
			outfile='None'
			
	else:
		outfile=f'{args.out_dir.rstrip("/")}/kraken2_rep.txt'

	return outfile


def calculate_metrics(outfile):

	"""Go through kraken2 classifications and calculate precision and recall for each origin species.
	write the results to output file"""

	try:
		total_species_classified=sum([int(l.split('\t')[2]) for l in open(outfile, 'r') if l.split('\t')[3]=='S'])
		FP=[l.split('\t')[2] for l in open(outfile, 'r') if l.split('\t')[3]=='S' and l.split('\t')[5].rstrip('\n').lower().lstrip(' ')==' '.join(os.path.dirname(outfile)\
			.split('/')[-1].split('_')[:2]).lower()]
        
		if len(FP)<1:
			FP=0
			FPR=0
		else:
			FP=FP[0]
			FPR=float(FP)/total_species_classified
		
	except Exception as e:
		print(f'{e}')	
		total_species_classified=None
		FP=None
		FPR=None

	return (os.path.dirname(outfile).split('/')[-1].split('_')[-1], ' '.join(os.path.dirname(outfile).split('/')[-1].split('_')[:2]), FP, total_species_classified, FPR)

def group_analysis():

	"""Go through several directories with mock communities, calculate metrics and save result to dataframe"""
	
	data={'Community':[],'Species':[], 'FP':[], 'total_species_reads':[], 'FPR':[]}

	outfiles=[f'{args.analyse_multiple.rstrip("/")}/{d}/kraken2_rep.txt' for d in \
	os.listdir(args.analyse_multiple) if os.path.isdir(f'{args.analyse_multiple.rstrip("/")}/{d}/')]

	for o in outfiles:

		result=calculate_metrics(o)
		data['Community'].append(result[0])
		data['Species'].append(result[1])
		data['FP'].append(result[2])
		data['total_species_reads'].append(result[3])
		data['FPR'].append(result[4])

	df=pd.DataFrame(data)
	df.to_csv(f'{args.out_dir}/community_analysis.csv',index=False, sep=',')	

def main():

	global args
	args=parse_args()
	outfile=run_kraken2()
	if args.analyse_multiple!=False:
		group_analysis()
	else:
		metrics=calculate_metrics(outfile)

	print('Kraken2 classification metrics calculated!')

if __name__=='__main__':
	main()
