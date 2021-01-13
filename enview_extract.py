#!/usr/local/env python2.7
#extract one or several genes from sqlite context_db

#In v7.5, add 'extract all' option and genus/species filter, remove '-flank' parameter and replace with relative path to file

import sys, os, sqlite3, subprocess, argparse
from argparse import RawTextHelpFormatter

def parse_arguments():
	descr='%r\n\nExtract specified genes and flanking regions from db\
	\n%r' % ('_'*80, '_'*80)
	parser=argparse.ArgumentParser(description=descr.replace("'", ""), \
	formatter_class=RawTextHelpFormatter)
	parser.add_argument('-genes', help='genes to extract. If list: "qnrA" "qnrB" ...'\
	, nargs='+')
	parser.add_argument('-o', help='path to output directory', required=True)
	parser.add_argument('-db', help='path to sqlite db file in original directory', required=True)
	parser.add_argument('-id', help='percent identity threshold for genes to extract', required=True)
	parser.add_argument('-taxa', help='list of genera and/or species to extract.\n\
	By default all taxa are extracted ', default=False, nargs='+')

	args=parser.parse_args()

	return args

def extract(args):

	#Connect to db
	connection=sqlite3.connect(args.db)
	cursor=connection.cursor()
	
	#Loop through all genes to extract
	for gene in args.genes:
		
		#Join tables based on genome ids
		query=f"""SELECT \
		args.arg_name, \
		args.id, \
		genomes.organism \
		FROM args \
		INNER JOIN genomes \
		ON args.genome_id=genomes.id \
		WHERE args.arg_name \
		LIKE \'{gene+"%"}\' \
		AND args.perc_id >= {float(args.id)} \
		"""

		if args.taxa!=False:
			i=0
			for taxon in args.taxa:
				i+=1
				if len(args.taxa)>1:
					if i==1:
						query+='AND ('
						query+=f'genomes.organism LIKE \'{taxon+"%"}\' '
					elif 1<i<len(args.taxa):
						query+=f'OR genomes.organism LIKE \'{taxon+"%"}\' '
					else:
						query+=f'OR genomes.organism LIKE \'{taxon+"%"}\') '


		cursor.execute(query)

		results=cursor.fetchall()

		#write create result directory if not exists
		if not os.path.exists(args.o.rstrip('/')+'/'+gene.lower()+'_'+str(args.id)+'_analysis'):
			os.makedirs(args.o.rstrip('/')+'/'+gene.lower()+'_'+str(args.id)+'_analysis')
		
		#Using the id, go back to the file containing the flanking regions and extract the ones matching the ids 
		flank_dict={}
		for line in open(os.path.dirname(args.db).rstrip('/')+'/all_flanks.csv_tmp', 'r'):
			flank_dict[line.split('\t')[0]]=line.split('\t')[1]

		arg_ids=[result[1] for result in results]

		#write results to output directory 
		with open(args.o.rstrip('/')+'/'+gene.lower()+'_'+str(args.id)+'_analysis/'+gene+'_contexts.fna', 'w') as outfile:	
			for result in results:
				outfile.write('>'+result[0]+'__'+str(result[1])+'__'+result[2]\
				+'\n'+flank_dict[str(result[1])]+'\n')
		

def main():
	 args=parse_arguments()
	 extract(args)

if __name__ =='__main__':
	main()
