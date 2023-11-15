import sys, os

"""Parse assembly summary file to extract accessions for origin genera"""

def parse_file(file, origins):

    lines=[l for l in open(file) if not l.startswith('#') and any(o.split(' ')[0]\
            in l.split('\t')[7] for o in origins) and not any(o in l.split('\t')[7] for o in origins)]
    
    with open(f'{file}.origin_genera_refs.txt', 'w') as f:
        for l in lines:
            f.write(l.split('\t')[0]+'\n')

def main():

    file=sys.argv[1] #assembly_summary.txt
    origins=['Proteus terrae', 'Providencia rettgeri', 'Providencia stuartii', 'Aeromonas media', 'Atlantibacter hermannii', 'Morganella morganii', \
		'Citrobacter amalonaticus', 'Kluyvera cryocrescens', 'Pseudochrobactrum asaccharolyticum']
    
    parse_file(file, origins)

if __name__=='__main__':
    main()
