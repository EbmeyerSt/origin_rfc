import sys, os
from fasta_to_dict import fasta_to_dict

""" Filter out contigs belonging to origin species from kraken2 bacteria standard library file"""

def read_fasta(file):

	""" Generator function that yields each header-sequence pair as it iterates over the file"""

	first=True
	for line in open(file):
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

	
def main():

	lib=sys.argv[1] #path to kraken2 library fasta file

	origins=['Proteus terrae', 'Providencia rettgeri', 'Providencia stuartii', 'Providencia thailandensis', 'Pseudochrobactrum', \
		'Ochrobactrum', 'Aeromonas media', 'Atlantibacter hermannii', 'Morganella morganii', \
		'Citrobacter amalonaticus', 'Kluyvera cryocrescens']
	
	print('Determining number of sequences...')
	num_seqs=[l for l in open(lib, 'r') if l.startswith('>')]
	
	print('filtering...')
	with open(f'{lib}.origins_filtered.fna', 'w') as f:
		for i, (header, seq) in enumerate(read_fasta(lib)):
			print(f'Processing seq {i}/{len(num_seqs)}...')
			if not any(o.lower() in header.lower() for o in origins):
				f.write(header+seq)

	print('Done!')

if __name__=='__main__':
	main()
