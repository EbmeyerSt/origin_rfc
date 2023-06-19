import sys, os
import subprocess
import pandas as pd

def sum_to_tabular(sum_files):

	rows=[]
	for file in sum_files:
		line_rows=[]
		mobile_blocks=0
		chromosomal_blocks=0
		#parse lines for info contained in blocks
		for line in open(file, 'r'):
			print(line)
			if line.startswith('Block'):

				block=[]

				block.append(f'{os.path.basename(file).split("_summary")[0]}')
				block.append(line.split(' ')[0])
				block.append(line.split('contains ')[1].split(' seqs')[0])
				block.append(line.split(' genera')[0].split('from ')[1])
				block.append(line.split(' species')[0].split('and ')[1])
			if '%)' in line:
				block.append(line.split('(')[1].split('%)')[0])
				if float(line.split('(')[1].split('%)')[0])>50:
					mobile_blocks+=1
				else:
					chromosomal_blocks+=1

			if line.startswith('Percent'):
				block.append(line.split(': ')[1].rstrip('\n'))
			if line.startswith('Genera'):
				block.append(line.split('Genera: ')[1].rstrip('\n'))
			if line.startswith('Species'):
				block.append(line.split('Species: ')[1].rstrip('\n'))

				if 'plasmid' in line:
					block.append(1)
				else:
					block.append(0)

			if line.startswith('Mean'):
				block.append(line.split('length: ')[1].rstrip('\n'))

				line_rows.append(block)

		for row in line_rows:
			row.extend([mobile_blocks, chromosomal_blocks])

		for row in line_rows:
			rows.append(row)

	return rows

def main():

	target_dir=sys.argv[1] #Directory containing summary files from find_origins_v5.py
	sum_files=[f'{os.path.abspath(target_dir)}/{file}' for file in os.listdir(target_dir)\
			if file.endswith('summary.txt')]

	blocks=sum_to_tabular(sum_files)
	df=pd.DataFrame(blocks, columns=['Gene', 'Block', 'Num_seqs','Num_genera', 'Num_species', 'Perc_mobile','Id_range','Genera', 'Species','Plasmid_in_block','Mean_seq_len', 'Mob_blocks', 'Chrom_blocks'])
	df.to_csv(f'{os.path.abspath(target_dir)}/summary_data.csv', index=False, sep=',')

if __name__=='__main__':
	main()
