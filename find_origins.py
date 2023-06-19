#!/local/usr/env python3.8

import sys, os, subprocess, argparse, sqlite3
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
import random
import math
import logging

def parse_arguments():
	descr='Identify potential gene origins based on alignments'
	parser=argparse.ArgumentParser(description=descr.replace("'", ''), formatter_class=RawTextHelpFormatter)
	parser.add_argument('-gene', help='name of gene to analyze, "all" to analyze all genes', required =True)
	parser.add_argument('-db', help='path to genview database used to extract sequences', required=True)
	parser.add_argument('-d', '--directory', help='absolute path to directory containing directories with context files', required=True)
	parser.add_argument('-c', '--cluster', help='clustering threshold for cd-hit', default=0.95)
	parser.add_argument('-mc', '--more_centroids', help='Select more sequences in addition to centroids', action='store_true')
	args=parser.parse_args()
	
	return args

def extract_files():

	if args.gene=='all':
		context_files=[f'{root}/{file}' for root, dirs, files in os.walk(args.directory)\
			for file in files if file.endswith('_contexts.fna')]

	else:

		context_files=[f'{root}/{file}' for root, dirs, files in os.walk(args.directory)\
			for file in files if file.endswith('_contexts.fna') and\
			args.gene.lower() in file.lower()]

	if len(context_files)>0:
		if not os.path.exists(f'{os.path.abspath(args.directory)}/context_files'):
			os.mkdir(f'{os.path.abspath(args.directory)}/context_files')

		for file in context_files:

			c_file=file

			try:

				#Read centroid seqs to dict
				all_seqs={}
				for line in open(file, 'r'):
					if line.startswith('>'):
						header=line
						seq=''
						gid=line.split('__')[-2]
						all_seqs[gid]={}
						all_seqs[gid]['header']=header
						all_seqs[gid]['genus']=line.split('__')\
							[-1].split(' ')[0]
					else:
						seq+=line
						all_seqs[gid]['seq']=seq
						all_seqs[gid]['seq_len']=len(seq)

				#Cluster seqs to reduce redundancy
				file=cluster(file)

				#Extract additional sequences if less than 4 genera are present
				#in cluster
				if args.more_centroids==True:
					print(f'Including more centroids...')
					#Read clusters into dict
					clusters={}
					for line in open(f'{file}.clstr', 'r'):
						if line.startswith('>C'):
							clust=line.rstrip('\n')
							clusters[clust]={}
							clusters[clust]['members']=[]

						elif '*' in line:
							clusters[clust]['centroid']=line\
								.split('__')[-2]

							clusters[clust]['members'].append(line.split('__')[-2])
						else:
							clusters[clust]['members'].append(line.split('__')[-2])

					print(f'Clusters: {len(clusters.keys())}')

					#Now go through each cluster and assess taxonomic diversity
					all_centroid_seqs=[]
					#genera to drop:uncultured,eubacteriaceae, unidentified, remove square brackets, type-E, stuff ending on eae, Proteobacteria, stuff ending on ales, IncQ, alpha, unverified:, Gammaproteobacteria, type-D, bacterium, '', unverified_asmbly, 'gamma', '', unknown, ''Betaproteobacteria, 'unverified_org', SAR86, type-c'

					trash_taxa=['uncultured', 'unidentified', 'type-e', 'type-d', 'proteobacteria', 'gammaproteobacteria', 'bacterium', 'unverified_asmbly', 'gamma', 'unknown', 'betaproteobacteria', 'unverified_org', 'sar86', 'type-c']

					for key, value in clusters.items():
						genera={all_seqs[i]['genus'] for i in\
							clusters[key]['members'] if not \
							all_seqs[i]['genus'].lower() in trash_taxa \
							and not all_seqs[i]['genus'].rstrip('\n')\
							.endswith('eae') and not \
							all_seqs[i]['genus'].rstrip('\n').endswith('ales')}
						avg_seq_len=np.mean([all_seqs[i]['seq_len'] for i in \
							clusters[key]['members']])

						#If average sequence length is >17000 and taxonomic
						#diversity low (<=2), select more sequences from the 
						#cluster to include in the centroid file
						if len(genera) <= 2 and avg_seq_len >= 17000 and \
						len(clusters[key]['members'])>1:

							if len(clusters[key]['members'])<7:
								all_centroid_seqs.extend(\
								clusters[key]['members'])

							else:
								all_centroid_seqs.extend(\
								random.sample(clusters[key]\
								['members'], 7))
						else:
							all_centroid_seqs.append(clusters[key]\
									['centroid'])

					#Write to new centroid file
					with open(file, 'w') as f:
						for i in set(all_centroid_seqs):
							f.write(all_seqs[i]['header']+all_seqs[i]['seq'])

				#Read centroid seqs to dict
				context_seqs={}
				for line in open(file, 'r'):
					if line.startswith('>'):
						header=line.rstrip('\n').lstrip('>').split(' ')[0]
						seq=''
					else:
						seq+=line
						context_seqs[header]=seq

				#Check how many sequences are present for a gene - if >5000, preprocess
				num_genomes=len([line for line in open(file, 'r') if line.startswith('>')])
				iterations=0

				print(f'Number of genomes before filtering: {num_genomes}')
				if num_genomes > 5000:
					print(f'Trying to reduce the number of sequences ({num_genomes}) in {file}')
					file=preprocess(file)
					num_genomes=len([line for line in open(file, 'r') if line.startswith('>')])
					print(f'number of genomes after filtering: {num_genomes}')

					if num_genomes > 5000:

						print(f'too many sequences ({num_genomes}) in {file} to handle, splitting file...')
						split_files=split_file(file)

						#Now append the new files to the end of the list
						context_files.extend(split_files)

						continue

				#Use blastn to align seqs to one another
				file=blastn(f'{file}')

				#Check if blast output file is >0, otherwise abort
				if not os.path.exists(f'{file}') or os.path.getsize(f'{file}')==0:
					print('Blast output size is 0 or does not exist, aborting...')
					continue		
					

				#1. Filter out overlapping alignments
				if not os.path.exists(f'{file}.filtered'):
					print('Filtering alignments...')
					file, filtered_lines=filter_blastaln(f'{file}')
					
					print('Calculating alignment length...')
					#Identify groups of seqs that form long alignments with one another
					#calculate length of alignments
					aln_dict={}
					for line in filtered_lines:
						#If the sequence pair is already in the dictionary, skip the occassion where
						#query is subject (i.e seq1 seq2 in dict, skip seq2 seq1 alns)
						if line[1]+' '+line[0] in aln_dict:
							continue	

						if not ' '.join(line[0:2]) in aln_dict:
							aln_dict[' '.join(line[0:2])]={}
							aln_dict[' '.join(line[0:2])]['aln_len']=int(line[4])
							aln_dict[' '.join(line[0:2])]['qlen']=int(line[3])

						else:
							aln_dict[' '.join(line[0:2])]['aln_len']+=int(float(line[4]))
						
					#for key, value in aln_dict.items():
						#print(f'{key} aln_len:{aln_dict[key]["aln_len"]} qlen:{aln_dict[key]["qlen"]}')

					#Write intermediate file 
					if not os.path.exists(f'{file}.alnblocks.txt'):
						print('Creating alignment blocks...')
						#Now check which seqs have long alignments to other seqs and create 'alignment blocks'
						block_list=[]

						#Create dict containing 1 for all seqs that already are in an alignment block
						all_blockmembers={e:0 for e in [i for key in list(aln_dict.keys()) for i in\
						key.split(' ')]}

						outcasts=[]
						for key, value in aln_dict.items():

							#Check if item is already in block list
							el1=key.split(' ')[0] in (element for sublist in block_list for element in sublist)
							el2=key.split(' ')[1] in (element for sublist in block_list for element in sublist)

							#Check if alignment is sufficiently long, sufficient % of query 
							#is covered AND seqs are not part of any alignment block yet
							#if so, initialize new block

						#OLD, percent cutoff if aln_dict[key]['qlen']>=15000 and aln_dict[key]['aln_len']/float(aln_dict[key]['qlen'])>=13000 \
							if aln_dict[key]['qlen']>=15000 and aln_dict[key]['aln_len']>=13000 \
							and not (el1==True or el2==True):

								#Create a list with all seq-alns that have a sufficient
								#ly long alignment with the key sequence - append this 
								#these to the block list
								sublist=[key2 for key2, value2 in aln_dict.items() if any(subkey in key2 for subkey in key.split(' ')) and aln_dict[key2]['qlen']>=15000 and aln_dict[key2]['aln_len']>=13000]
								single_list=[]
								for element in sublist:
									single_list.extend([element.split(' ')[0], \
									element.split(' ')[1]])

								sublist=[e for e in list(set(single_list)) if not \
									all_blockmembers[e]==1]
								#Check if any member in sublist is already
								#part of an alignment block. if so, remove
								

								for x in sublist:
									all_blockmembers[x]=1

								block_list.append(list(sublist))

						#Save entries which did not qualify for alignment block in list 
						all_seq=[key.split(' ')[0:2] for key, value in aln_dict.items()]
						all_flat=[i for sublist in all_seq for i in sublist]

						#Create own block for species without long alignment that are not yet
						#in any other block
						outcasts={i for i in all_flat if all_blockmembers[i]==0 and not i.split('__')\
								[-1] in [tax.split('__')[-1] for tax in [taxkey for taxkey in all_blockmembers.keys() if all_blockmembers[taxkey]==1]]}

						#Write output file containing seq accession and block number		
						blocknum=0
						with open(f'{file}.alnblocks.txt', 'w') as blockfile:
							for block in set(tuple(subblock) for subblock in block_list):
								blocknum+=1
								for header in block:
									blockfile.write(f'Block{blocknum} {header}'+'\n')
							if len(outcasts)>0:
								for header in outcasts:
									blockfile.write(f'BlockO {header}'+'\n')

				if not os.path.exists(f'{file}.alnblocks.is_out.csv'):
					print('Searching for MGEs')
					#now extract seqs in aln_blocks and search them for MGEs

					with open(f'{file}.alnblocks.fna', 'w') as outfile:
						for line in open(f'{file}.alnblocks.txt', 'r'):
							outfile.write(">"+line.split(' ')[1].rstrip('\n')\
							+'\n'+context_seqs[line.split(' ')[1].rstrip('\n')]+'\n')

					#Search MGEs
					mge_file=search_mges(f'{file}.alnblocks.fna')

					#Now add the information about MGEs to the alignment block file
					blocks=[line for line in open(f'{file}.alnblocks.txt', 'r')]
					mges=[line.split('\t')[0] for line in open(f'{mge_file}')]
					
					#Append sequences shorter than 10000bp to mobile list
					for key, value in context_seqs.items():
						if len(value)<10000:
							mges.append(key.split(' ')[0].lstrip('>'))
					
					mges=set(mges)

					with open(f'{file}.alnblocks_mges.txt', 'w') as outfile:
						for line in blocks:
							if line.split(' ')[1].rstrip('\n') in mges:
								outfile.write(line.rstrip('\n')+' MGE\n')
							else:
								outfile.write(line.rstrip('\n')+' absent\n')

				if not os.path.exists(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file).split(".fna")[0]}_summary.txt'):
					#Read results from block file into dict
					block_lines=[line for line in open(f'{file}.alnblocks_mges.txt', 'r')]	
					block_dict={}
					for line in block_lines:
						if not line.split(' ')[0] in block_dict:
							block_dict[line.split(' ')[0]]={}
						block_dict[line.split(' ')[0]][line.split(' ')[1]]={}
						if line.split(' ')[2].rstrip('\n')=='MGE':
							block_dict[line.split(' ')[0]][line.split(' ')[1]]['mge']='present'
						else:
							block_dict[line.split(' ')[0]][line.split(' ')[1]]['mge']='absent'
							
					#Integrate database in formation into block_dict - which species, which % similarity to the reference
					connection=sqlite3.connect(args.db)
					cursor=connection.cursor()
					print('Querying database for additional information...')

					id_dict={}
					ids=[key2.split('__')[-2] for key, value in block_dict.items() for key2, value2 in block_dict[key].items()]
					for element in ids:
						query="""
						SELECT \
						args.id, \
						args.perc_id, \
						genomes.organism \
						FROM args \
						INNER JOIN genomes ON genomes.id=args.genome_id \
						WHERE args.id = ?
						"""

						cursor.execute(query, (int(element),))
						results=cursor.fetchall()

						for result in results:
							id_dict[int(result[0])]={}
							id_dict[int(result[0])]['perc_id']=result[1]
							id_dict[int(result[0])]['taxon']=result[2]


					#Now write results to summary file
					with open(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file).split(".fna")[0]}_summary.txt', 'w') as outfile:
						outfile.write('_'*50+'\n\n')
						outfile.write(f'Target gene (family): {os.path.basename(file).split("_")[0]}'+'\n')
						outfile.write('-'*50+'\n')
						outfile.write(f'Identified {len(block_dict)} blocks with alignments > 10kb'+'\n\n')
						for key, value in block_dict.items():
							block_ids=[key2.split('__')[-2] for key2, value2 in block_dict[key].items()]

							block_taxa={key2.split('__')[-1] for key2, value2 in block_dict[key].items() if not key2.split('__')[-1].lower() in trash_taxa and not (key2.split('__')[-1].endswith('ales') or key2.split('__')[-1].endswith('eae'))}
							block_species={' '.join(id_dict[int(key2.split('__')[-2])]['taxon'].split(' ')[0:2]) for key2, value2 in block_dict[key].items()}
							mge_taxa=[block_dict[key][key2] for key2, value2 in block_dict[key].items()\
								if block_dict[key][key2]['mge']=='present']
							id_range=[float(id_dict[int(key2.split('__')[-2])]['perc_id']) for key2, value2 in block_dict[key].items()]
							mean_seq_len=round(np.mean([len(context_seqs[con_key]) for id in block_ids for con_key, con_val in context_seqs.items() if int(id)==int(con_key.split('__')[-2])]),2)


							outfile.write(f'{key} contains {len(block_dict[key])} seqs from {len(block_taxa)} genera and {len(block_species)} species'+'\n')
							outfile.write(f'{len(mge_taxa)} ({round(((len(mge_taxa)/float(len(block_dict[key])))*100),1)}%) of the seqs contain at least one MGE/Integrase'+'\n')
							outfile.write(f'Percent identity range to reference: {min(id_range)}-{max(id_range)}'+'\n')
							outfile.write(f'Genera: {",".join(block_taxa)}'+'\n')
							outfile.write(f'Species: {",".join(block_species)}'+'\n')
							outfile.write(f'Mean sequence length: {mean_seq_len}'+'\n\n')

						outfile.write('_'*50+'\n\n')
						print(f'Analysis finished for {c_file}')

			except Exception as e:
				print(f'EXCEPTION while analyzing {c_file}')
				logger.exception(str(e))

	else:
		print('No files matching the specified gene found, exiting...')
		sys.exit()

def filter_blastaln(blast_alnfile):

	linelist=[line.split('\t') for line in open(blast_alnfile, 'r') if not line.split('\t')[1]==line.split('\t')[0] and int(line.split('\t')[4])>300]

	alns={}
	#Go through every line in alignment
	for line in linelist:
		key=' '.join(line[0:2])
		unique_pos=0
		#Append qseq sseq pair as key in alignment dict if it does not already exist
		if not key in alns.keys():
			alns[key]={}
		#If no alns are present for the pair, append the first one automatically
		if len(alns[key].keys())==0:
			alns[key]['aln1']=line
		else:
			#If alns are present, go through all and check whether the two alns overlap. if yes, discard the shorter one
			for key2, value2 in alns[key].items():
		    
			#Check if new alignment is enclosed by larger aln - if yes, discard
				if int(line[6])>=int(value2[6]) and int(line[7])<=int(value2[7]):
					pass
		            
			#Check if new alignment encloses smaller alignment - if yes, replace shorter aln with this one
				elif int(line[6])<=int(value2[6]) and int(line[7])>=int(value2[7]):
					alns[key][key2]=line
		            
			#If none of both is the case, increase unique_pos by 1. If it is increased for all keys, append to dict
				else:
					unique_pos+=1
		        
		if unique_pos==len(alns[key].keys()):    
			alns[key]['aln'+str(len(alns[key].keys())+1)]=line

	aln_list=sorted([alns[key][key2] for key, value in alns.items() for key2, value2 in alns[key].items()],key=lambda x: (x[0],x[1]))

	#Write filtered lines to outfile
	with open(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(blast_alnfile)}.filtered', 'w') as outfile:
		for line in aln_list:
			outfile.write('\t'.join(line))

	filename=f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(blast_alnfile)}.filtered'
	return filename, aln_list

def cluster(file):

	gene_name=os.path.basename(file).split('_')[0]
	print(f'Analyzing {gene_name}...')

	#Check number if seqs in context file exceeds 5000 - if so, do mobility analysis first,
	#then select a number of mobile elements to keep at random 
	
	if not os.path.exists(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.centroids'):
		print('Clustering sequences...')
		#Cluster sequences in file with 95% id to remove redundancy
		cluster_command=f'cd-hit-est -T 0 -M 0 -c {args.cluster} -s 0.7 -n 9 -d 50 -i {file} -o {os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.centroids'
		try:
			subprocess.call(cluster_command, shell=True)
			filename=f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.centroids'

		except Exception as e:
			print('Clustering of sequences failed, aborting...')
			print(f'{str(e)}')
	else:
		filename=f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.centroids'

	return filename

def blastn(file):

	if not os.path.exists(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.centroids.blastout'):

		#Create blast database to enable multithreading
		createdb=f'makeblastdb -in {file} -dbtype "nucl" -out {file}.blastdb'
		subprocess.call(createdb, shell=True)

		print('Comparing sequences...')
		blastn=f'blastn -query {os.path.abspath(args.directory)}/context_files/{os.path.basename(file)} -db {file}.blastdb -perc_id 70 -strand both -task blastn -out {os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.blastout -max_target_seqs 5000 -num_threads 48 -outfmt "6 qseqid sseqid pident qlen length score qstart qend sstart send"'
	
		#Remove blast database files
		db_files=[file for file in os.listdir(os.path.dirname(file)) if file.endswith('.blastdb')]
		for db_file in db_files:
			os.remove(db_file)

		try:
			subprocess.call(blastn, shell=True)
			filename=f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.blastout'

		except Exception as e:
			print('Could not run blastn, aborting...')
			print(f'{str(e)}')
			return
	return filename

def search_mges(file):

	if not os.path.exists(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out.csv'):
		#Search the seqs for MGEs, integrons, other ARGs
		diamond=f'diamond blastx -p 30 -d /storage/stefan/databases/mobileOG/mobileog-db_v1.6_manual_nophage.dmnd -q {file} -o {os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out.csv --id 90 --ultra-sensitive --top 100 --masking 0 --subject-cover 90 -f 6 qseqid sseqid stitle pident qstart qend qlen slen length score qseq qframe qtitle'
		try:
			subprocess.call(diamond, shell=True)

		except Exception as e:
			print('Error running diamond, aborting...')
			print(f'{str(e)}')
			return

	print('Filtering hits...')
	#Sort by queryID and start position
	dmnd_rawout=pd.read_csv(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out.csv', sep='\t', names=['qseqid', 'sseqid', 'stitle', 'pident', 'qstart', 'qend', 'qlen', 'slen', 'length', 'score', 'qseq', 'qframe', 'qtitle'])

	dmnd_rawout.sort_values(by=['qseqid', 'qstart']).to_csv(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out_annotated.csv', sep='\t', index=False, header=False)

	#Assign same start and end position if several gene variants hit the same locus to be able to group by exact position and filter best hit
	with open(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out_annotated_rep.csv', 'w') as rawout:
		line_number=0
		for line in open(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out_annotated.csv', 'r'):
			line_number+=1
			if line_number==1:
				if not line.split('\t')[-2].startswith('-'):
					qseqid=line.split('\t')[0]
					seqstart=line.split('\t')[4]
					seqend=line.split('\t')[5]
					rawout.write(line)
				else: 

					qseqid=line.split('\t')[0]
					seqstart=line.split('\t')[5]
					seqend=line.split('\t')[4]
					rawout.write(line)
			
			else:	#Do this one time for forward, another time for reverse strand (for reverse strand, switch seqstart and seqend)
				if not line.split('\t')[-2].startswith('-'):
					if line.split('\t')[0]==qseqid:
						if int(line.split('\t')[4]) in range(int(seqstart)-400, int(seqstart)+400) and int(line.split('\t')[5]) in range(int(seqend)-400, int(seqend)+400):
							newline=line.split('\t')
							newline[4]=seqstart
							newline[5]=seqend
							rawout.write('\t'.join(newline))
						else:

							seqstart=line.split('\t')[4]
							seqend=line.split('\t')[5]
							rawout.write(line)
					else:

						qseqid=line.split('\t')[0]
						seqstart=line.split('\t')[4]
						seqend=line.split('\t')[5]
						rawout.write(line)
				else:

					if line.split('\t')[0]==qseqid:
						if int(line.split('\t')[5]) in range(int(seqstart)-400, int(seqstart)+400) and int(line.split('\t')[4]) in range(int(seqend)-400, int(seqend)+400):
							newline=line.split('\t')
							newline[5]=seqstart
							newline[4]=seqend
							rawout.write('\t'.join(newline))
						else:

							seqstart=line.split('\t')[5]
							seqend=line.split('\t')[4]
							rawout.write(line)
					else:

						qseqid=line.split('\t')[0]
						seqstart=line.split('\t')[5]
						seqend=line.split('\t')[4]
						rawout.write(line)


	#Select best hit
	dmnd_hits=pd.read_csv(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out_annotated_rep.csv', sep='\t', names=['qseqid', 'sseqid', 'stitle', 'pident', 'qstart', 'qend', 'qlen', 'slen', 'length', 'score', 'qseq', 'qframe', 'qtitle'])
	best_hits=dmnd_hits.loc[dmnd_hits.groupby(['qseqid', 'qstart', 'qend'])['score'].idxmax()]
	best_hits.to_csv(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file).rstrip(".fna")}.is_out_annotated.csv', sep='\t', index=False, header=False)

	hit_file=f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file).rstrip(".fna")}.is_out_annotated.csv'

	#remove intermediate files
	os.remove(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out_annotated_rep.csv')
	#os.remove(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.is_out_annotated.csv')
	print('Temporary files removed')
	return hit_file

def preprocess(file):

	print('Prefiltering sequences with MGEs...')
	#Read into dict
	seqs={}
	for line in open(file, 'r'):
		if line.startswith('>'):
			header=line
			seq=''
		else:
			seq+=line
			seqs[header]=seq
	
	#search file for MGEs
	mge_file=search_mges(file)

	#Get list of all seqs containing an MGE
	mge_contexts=[line.split('\t')[0] for line in open(mge_file, 'r')]

	#Append all sequences that are shorter than 6000bp to the mobile list
	for key, value in seqs.items():
		if len(value)<10000:
			mge_contexts.append(key.lstrip('>').split(' ')[0])

	mge_contexts=set(mge_contexts)

	#From these, randomly select 500, or as many as available if < 500
	if len(mge_contexts)>=300:
		sample=random.sample(mge_contexts, 300)

	elif len(mge_contexts)<300:
		#If less than 500 have mges, just return
		return file

	#Now go through dict and filter out all MGE containing seqs that hav not been randomly sampled
	with open(f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.filtered', 'w') as outfile:
		for key, value in seqs.items():
			if not key.split(' ')[0].lstrip('>') in mge_contexts:
				outfile.write(key+value)
			else:
				if not key.split(' ')[0].lstrip('>') in sample:
					pass
				else:
					outfile.write(key+value)

	filename=f'{os.path.abspath(args.directory)}/context_files/{os.path.basename(file)}.filtered'
	return filename

def split_file(file):

	#Read file into dict
	seqs={}
	for line in open(file, 'r'):
		if line.startswith('>'):
			header=line
			seq=''
		else:
			seq+=line
			seqs[header]=seq

	#Based on the number of seqs, write the sequences to files such that no more than 5000 seqs are in each file
	file_num=math.ceil(len(seqs)/5000)

	ind=0
	keys=list(seqs.keys())

	split_files=[]
	for i in range(file_num):
		print(f'file num: {file_num}, i:{i}, ind:{ind}, end: {(len(keys)/file_num)+ind}')
		with open(file.replace('_context', '-'+str(i)+'_context'), 'w') as f:
			if not i==file_num-1:
				subkeys=keys[ind:math.ceil(len(keys)/file_num)+ind]
				for subk in subkeys:
					f.write(subk+seqs[subk])
				ind+=math.ceil(len(keys)/file_num)

			else:
				subkeys=keys[ind:]
				for subk in subkeys:
					f.write(subk+seqs[subk])
				ind+=math.ceil(len(keys)/file_num)
		f.close()
		split_files.append(file.replace('_context', '-'+str(i)+'_context'))			
	
	return split_files


def main():

	global args
	args=parse_arguments()
	
	global logger
	logger=logging.getLogger()

	extract_files()

if __name__=='__main__':
	main()

