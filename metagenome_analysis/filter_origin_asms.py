import sys, os
from fasta_to_dict import fasta_to_dict

""" Filter out contigs shorter than 5000bp or having 'plasmid' in the header"""

def filter(contigs):

    conts=fasta_to_dict(contigs)
    filtered={k:v for k, v in conts.items() if len(v) > 5000 and not 'plasmid' in k.lower()}

    with open(f'{contigs}.filtered.fna', 'w') as f:
        for k, v in filtered.items():
            f.write(k+v)

def main():

    conts=sys.argv[1] #Path to fasta file with origin contigs
    filter(conts)

if __name__=='__main__':
    main()