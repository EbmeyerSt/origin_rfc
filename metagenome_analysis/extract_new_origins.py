"""Use specified genview database to extract all origin species containig the respective gene with >=90% 
sequence similarity"""

import sys, os
import sqlite3

def extract(db_path):

	"""Extract assembly accessions of origin species carrying the respective gene """

	con=sqlite3.connect(f'{db_path}')
	cursor=con.cursor()

	origins=[('Providencia stuartii', 'tetb'), ('Providencia stuartii', 'catiii'), ('Providencia rettgeri', 'tet59'), ('Proteus terrae', 'teth'), ('Atlantibacter hermannii', 'cati'), ('Atlantibacter hermannii', 'hera'),('Kluyvera cryocrescens', 'kluc'), ('Aeromonas media', 'oxa-427'), ('Morganella morganii','tetd'), ('Morganella morganii', 'catii'), ('Citrobacter amalonaticus', 'cdia'), ('Pseudochrobactrum', 'aac6-ian')]

	asms=[]

	for i in origins:
		if i[1]=='hera' or i[1]=='kluc':
			query=f"""
			SELECT assembly FROM genomes
			INNER JOIN args ON genomes.id=args.genome_id
			WHERE LOWER(arg_name) LIKE \"{i[1].lower()+"%"}\"
			AND LOWER(genomes.organism) LIKE \"{i[0].lower()+"%"}\"
			AND args.perc_id >= 90;
			"""
		else:
			query=f"""
			SELECT assembly FROM genomes
			INNER JOIN args ON genomes.id=args.genome_id
			WHERE LOWER(arg_name) = \"{i[1].lower()}\"
			AND LOWER(genomes.organism) LIKE \"{i[0].lower()+"%"}\"
			AND args.perc_id >=90;
			"""
		print(query)
		cursor.execute(query)
		results=cursor.fetchall()

		asms.extend(results)

	with open(f'{os.path.dirname(db_path)}/origin_accessions.txt', 'w') as f:
		for i in asms:
			f.write(i[0]+'\n')


def main():

	db_path=sys.argv[1]
	extract(db_path)


if __name__=='__main__':
	main()
