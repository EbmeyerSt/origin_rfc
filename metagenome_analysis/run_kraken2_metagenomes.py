import sys
import os
import logging
import subprocess

"""Run kraken2 on selected metagenome samples"""

target_dir=sys.argv[1] #/storage/stefan/HMP/paired_reads
krakendb=sys.argv[2] #/storage/stefan/new_origins_kraken2db/origin_db_test3
outdir=sys.argv[3] #directory to write output to (/storage/stefan/new_origins_august23/metagenomes/metagenome_results)

logging.basicConfig(filename=f'{outdir.rstrip("/")}/run_kraken2.log', filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

read_pairs=[(f'{target_dir.rstrip("/")}/{f}', f'{target_dir.rstrip("/")}/{f.replace("1.fastq", "2.fastq")}') for f in os.listdir(target_dir) if f.endswith('1.fastq')]

#Check that read pairs actually are existing files
for p in read_pairs:
	if os.path.isfile(p[0]):
		pass
	if os.path.isfile(p[1]):
		pass
	if not os.path.isfile(p[0]):
		logging.warning(f'{p[0]} is not an existing file!')
	if not os.path.isfile(p[1]):
		logging.warning(f'{p[1]} is not an existing file!')

for p in read_pairs:
	try:
		run_kraken2=(f'nice -n 5 kraken2 --db {krakendb} --paired {p[0]} {p[1]} --output {outdir.rstrip("/")}/{os.path.basename(p[0]).replace("1.fastq", "kraken2_out.txt")}'
		f' --confidence 0.3 --threads 50 --report-zero-counts --report {outdir.rstrip("/")}/{os.path.basename(p[0]).replace("1.fastq", "kraken2_report.txt")} --use-names')
		subprocess.call(run_kraken2, shell=True)
	except Exception as e:
		logging.exception(f'Exception: {e}')
