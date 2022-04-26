import pandas as pd
import timeit
import numpy as np
from Bio import SeqIO
from shutil import copyfile
import os, sys, argparse, subprocess, multiprocessing, re, time, json, shutil, sqlite3
from argparse import RawTextHelpFormatter
from multiprocessing import Manager
from collections import defaultdict
from Bio.Seq import Seq
from os.path import expanduser

"""Tutorial command:
genview-makedb -d /path/to/output/directory -db /path/to/reference/PER.fna -p 10 -id 80 --taxa rheinheimera --assemblies

Alternatively, if you already know which genomes you want to compare, you can specify their GenBank accessions via --accessions, e.g:

genview-makedb -d /path/to/output/directory -db /path/to/reference/PER.fna -p 10 -id 80 --accessions /path/to/accessions.txt --assemblies
"""


def parse_arguments():
	man_description='Creates sqlite3 database with genetic environment from genomes containing the provided reference gene(s).'
	parser=argparse.ArgumentParser(description=man_description.replace("'", ""), formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d', '--target_directory', help='path to output directory', required=True)
	parser.add_argument('-db', '--database', help='fasta/multifasta file containing amino acid sequences of translated genes to be annotated', required=True)
	parser.add_argument('-p', '--processes', help='number of cores to run the script on', type=int, default=int(multiprocessing.cpu_count()/2))
	parser.add_argument('-id', '--identity', help='identity cutoff for hits to be saved to the database (e.g 80 for 80%% cutoff)', type=float, default=90)
	parser.add_argument('-scov', '--subject_coverage', help='minimum coverage for a hit to be saved to db (e.g 80 for 80%% cutoff)', type=float, default=90)
	parser.add_argument('--update', help='update an existing genview database with new genomes', action='store_true', default='False')
	parser.add_argument('--uniprot_db', help='Path to uniprotKB database', required=False, default='False')
	parser.add_argument('--uniprot_cutoff', help='%% identity threshold for annotating orfs aurrounding the target sequence, default 60', default=60)
	parser.add_argument('--taxa', help='taxon/taxa names to download genomes for - use "all" do download all available genomes, cannot be specified at the same time as --accessions', nargs='+', default='False')
	parser.add_argument('--assemblies', help='Search NCBI Assembly database ', action='store_true', default='False')
	parser.add_argument('--plasmids', help='Search NCBI Refseq plasmid database', action='store_true', default='False')
	parser.add_argument('--local', help='path to local genomes', default='False')
	parser.add_argument('--save_tmps', help='keep temporary files', action='store_true', default='False')
	parser.add_argument('--accessions', help='csv file containing one genome accession number per row, cannot be specied at the same time as --taxa', default='False')
	parser.add_argument('--integron_finder', help=argparse.SUPPRESS, default='False')
	parser.add_argument('--flanking_length', help='Max length of flanking regions to annotate', type=int, default=10000)
	parser.add_argument('--long-reads', help=argparse.SUPPRESS, action='store_true', default='False')
	parser.add_argument('--kraken2', help='Path to kraken2 database. Uses kraken2 to classify metagenomic long-reads.', default='False')
	parser.add_argument('--log', help='Write log file for debugging', action='store_true', default='False')
	parser.add_argument('--clean', help='Erase files from previous genview runs from target directory', action='store_true', default='False')
	args=parser.parse_args()

	return args

def download_uniprot():

	#Check if gdown is installed
	if args.uniprot_db=='False':
	
		print('\nPath to uniprot database (to annotate the target genes genetic environment) "uniprotKB.dmnd" not specified.\nIf you have run GEnView previously, please specify the full file path using --uniprot_db.\nElse, the database will be downloaded now!\n')
		down=input('Download uniprotKB database? (y/n) ')
		if down=='y':
			pass
		elif down=='n':
			sys.exit()
		else:
			print('Invalid input, please answer y or n')
			sys.exit()

		try:
			import gdown

		except ImportError:
			install=input('Package gdown is not installed - install now?(y/n)')
			if install=='y':
				subprocess.call('pip install gdown', shell=True)
			else:
				print('Package gdown is not installed. Use "pip install gdown" to install\
				or download database file anually from https://drive.google.com/uc?id=1VY70ab47Pu2fodYKKm1_fbJ83NGT1dNc')
				sys.exit()

		#download modified uniprot DB
		download_uniprot='gdown https://drive.google.com/uc?id=1VY70ab47Pu2fodYKKm1_fbJ83NGT1dNc -O %s'\
		% args.target_directory.rstrip('/')+'/uniprotKBjan2019.fna.gz'
		if not os.path.exists(args.target_directory.rstrip('/')+'/uniprotKBjan2019.fna.gz') and not os.path.exists(args.target_directory.rstrip('/')+'/uniprotKBjan2019.fna'):
			print('Decompressing database file...\n')
			subprocess.call(download_uniprot, shell=True)
			unzip=f'gunzip {args.target_directory.rstrip("/")+"/uniprotKBjan2019.fna.gz"}'
			subprocess.call(unzip, shell=True)

		#Check if download was succesfull
		if not os.path.exists(args.target_directory.rstrip('/')+'/uniprotKBjan2019.fna'):
			print('\nDownload of the uniprot database failed, possibly due to exceeded bandwith limit on storage.\nPlease try again in a few hours or download the database manually from "https://drive.google.com/uc?id=1VY70ab47Pu2fodYKKm1_fbJ83NGT1dNc".\nTo convert the file into a diamond database, unzip it using "gunzip filename.gz" and then run "diamond makedb --in databasefile.fa -d yourdbname.dmnd".\nSpecify the path to the so created database using --uniprot_db.\nIf possible, use a version of the database created by a previous run.\n')
			sys.exit()

		#Transform to diamond database
		if not os.path.exists(args.target_directory.rstrip('/')+'/uniprotKB.dmnd'):

			dmnd=f'diamond makedb --in {args.target_directory.rstrip("/")+"/uniprotKBjan2019.fna"} -d {args.target_directory.rstrip("/")+"/uniprotKB.dmnd"}'
			subprocess.call(dmnd, shell=True)

		if args.log==True:
			if os.path.exists(args.target_directory.rstrip("/")+"/uniprotKBjan2019.fna"):
				log_lines.append('UniprotKB downloaded...\n')
			else:
				log_lines.append('UniprotKB not found...FAILED\n')
			write_log()			

def reformat():

	#Check that db contains amino acid seqs
	db_lines=[line for line in open(args.database, 'r') if not line.startswith('>')]
	letter_set={letter for line in db_lines for letter in line}

	if len(letter_set)<12:
		print(f'{args.database} does not look like it contains amino acid sequences, please provide an amino acid sequence database instead.\nexiting...')
		sys.exit()


	if not os.path.exists(args.database.replace(args.database.split('.')[-1], 'dmnd')):
		print('Converting fasta to .dmnd database...')
		#Read in multifasta

		seqs={}
		for line in open(args.database, 'r'):
			if line.startswith('>'):
				header=line.rstrip('\n').lstrip('>')
				seq=''
			else:
				seq+=line
				seqs[header]=seq

		#Reformat into compatible header
		with open(args.database+'_reformatted.fna', 'w') as outfile:
			i=0
			for key, value in seqs.items():
				i+=1
				new_format=f'>gb|{key+str(i)}|NOARO{i}|{key} [unknown]\n{value}\n'
				outfile.write(new_format)
		outfile.close()

		#Transform to diamond database
		diamond=f'diamond makedb --in {args.database+"_reformatted.fna"} -d {args.database.replace(args.database.split(".")[-1], "dmnd")}'
		subprocess.call(diamond, shell=True)

		args.db_new=args.database.replace(args.database.split(".")[-1], "dmnd")
		args.database=args.db_new
	else:

		args.db_new=args.database.replace(args.database.split(".")[-1], "dmnd")
		args.database=args.db_new

	if args.log==True:
		if os.path.exists(args.database.replace(args.database.split(".")[-1], "dmnd")):
			log_lines.append('diamond database found...\n')
		else:
			log_lines.append('diamond database not found...FAILED\n')

		write_log()		

	return args.database

def reverse_complement(seq):

	#Calculate reverse complement

	try:
		sequence=Seq(seq)
		rev_seq=sequence.reverse_complement()
		reversedseq=True
	except:
		reversedseq=False
	
	if args.log==True:
		if reversedseq==True:
			log_lines.append('Sequences reversed...\n')
		else:
			log_lines.append('Exception while reversing sequences...FAILED\n')
		write_log()

	return rev_seq



def download_new(queue):

	if args.update==True:
		target_dir=f'{os.path.abspath(args.target_directory)}/update_tmp/genomes'
	else:

		target_dir=f'{os.path.abspath(args.target_directory)}/genomes'

	#download newly published genomes
	url=queue.get()
	if url=='STOP':
		return

	while True:

		if not os.path.isfile(target_dir+'/'+os.path.basename(url+'_genomic.fna.gz')) \
		and not os.path.isfile(target_dir+'/'+os.path.basename(url+'_genomic.fna')):

			try:
			
				download_fa='wget -r -nH --cut-dirs=7 %s -P %s' % \
				(url+'/'+os.path.basename(url+'_genomic.fna.gz'), target_dir)

				subprocess.call(download_fa, shell=True)
				time.sleep(1)
				
				#unzip
				unzip='gunzip %s' % target_dir+'/'+os.path.basename(url+'_genomic.fna.gz')
				subprocess.call(unzip, shell=True)

				#Give permissions
				permission='chmod a+r+w %s' % target_dir+'/'+os.path.basename(url+'_genomic.fna')
				subprocess.call(permission, shell=True)

			except Exception as e:
				print('EXCEPTION: %s' % e)

		
		elif os.path.isfile(target_dir+'/'+os.path.basename(url+'_genomic.fna')):
			print('%s already downloaded!' % url)

		url=queue.get()
		if url=='STOP':
			return

def split_fasta():

	if args.update==True:
		#Split concatenated fasta file into several smaller ones
		file=args.target_directory.rstrip('/')+'/'+'update_tmp/all_assemblies.fna'
		target_dir=args.target_directory.rstrip('/')+'/'+'update_tmp'
	else:
		file=args.target_directory.rstrip('/')+'/all_assemblies.fna'
		target_dir=args.target_directory.rstrip('/')

	line_num=0
	for line in open(file, 'r'):
		line_num+=1

	multiplicator=1
	iterated_lines=0

	for line in open(file, 'r'):
		iterated_lines+=1

		if line.startswith('>'):
			previous_complete=True
		else:
			previous_complete=False

		if iterated_lines==1:
			outfile=open(target_dir+'/all_assemblies_'+str(multiplicator)+'.fna', 'w')

		if iterated_lines<=(line_num/split_num)*multiplicator:	
			outfile.write(line)

		elif iterated_lines>=(line_num/split_num)*multiplicator and previous_complete==False:
			outfile.write(line)
		else:
			outfile.close()
			multiplicator+=1
			outfile=open(target_dir+'/all_assemblies_'+str(multiplicator)+'.fna', 'w')
			outfile.write(line)
	outfile.close()

	if args.log==True:
		asm_splits=[file for file in os.listdir(target_dir) if file.startswith('all_assemblies_')\
			and file.endswith('.fna')]
		
		if len(asm_splits)==split_num:
			log_lines.append('Files split...\n')
		else:
			log_lines.append('Files split incorrectly...FAILED(?)\n')
		write_log()

def concatenate_and_split():

	if args.update==True:
		target_dict=args.target_directory.rstrip('/')+'/'+'update_tmp'
	else:
		target_dict=args.target_directory.rstrip('/')

	if not os.path.isfile(target_dict.rstrip('/')+'/all_assemblies.fna'):
		print('Concatenating novel assemblies...')
		#Collect new genome fasta files
		if args.local=='False':
			new_genomes=[content[0].rstrip('/')+'/'+element for content \
			in os.walk(target_dict.rstrip('/')+'/genomes') \
			for element in content[2] if element.endswith('_genomic.fna')]

			#Write content of all fasta files to one file, append assembly name to respective contig
			with open(target_dict+'/all_assemblies.fna', 'w') as outfile:
				for file in new_genomes:
					for line in open(file, 'r'):
						if line.startswith('>'):
							outfile.write(line.rstrip('\n')+'__'+\
							file.split('/')[-1].split('_')[0]+'_'+\
							file.split('/')[-1].split('_')[1]+'\n')
						else:
							outfile.write(line)
					outfile.write('\n')
			outfile.close()
		else:
			new_genomes=[f'{os.path.abspath(args.local)}/{file}' for file in os.listdir(args.local)]

			#Write content of all fasta files to one file, append assembly name to respective contig
			with open(target_dict+'/all_assemblies.fna', 'w') as outfile:
				gvread_id=0
				for file in new_genomes:
					for line in open(file, 'r'):
						if line.startswith('>'):
							gvread_id+=1
							if ' ' in line:
								split_line=line.split(' ')
								split_line[0]=str(split_line[0])+f'.{gvread_id}'
								new_line=' '.join(split_line)	
							else:
								new_line=line.replace('\n', f'.{gvread_id}'+'\n')
							outfile.write(new_line)
						else:
							outfile.write(line)
					outfile.write('\n')
			outfile.close()

	if args.log==True:
		if os.path.exists(target_dict+'/all_assemblies.fna') and os.path.getsize\
			(target_dict+'/all_assemblies.fna')>0:
			log_lines.append('all_assemblies.fna created...\n')
		else:
			
			log_lines.append('all_assemblies.fna does not exist or is empty...FAILED\n')
		write_log()

	#Use kraken2 to classify long reads
	if args.local!='False':
		if args.kraken2!='False':
			kraken2_classify(target_dict+'/all_assemblies.fna')	

			if args.log==True:
				if os.path.exists(os.path.abspath(f'{args.target_directory}/kraken2/all_assemblies.kraken2.out'))\
				and os.path.exists(os.path.abspath(f'{args.target_directory}/kraken2/all_assemblies.kraken2.out')):	
					log_lines.append('Kraken2 run sucessfully...\n')
				else:
					log_lines.append('Kraken2 run unsuccessful...FAILED\n')
				write_log()

	#Split fasta file into smaller files
	print('Splitting into %d files...' % split_num)
	split_fasta()

	#Create .csv version of each file
	#Get list of all target files
	fna_files=[target_dict+'/'+file \
	for file in os.listdir(target_dict)\
	 if file.endswith('.fna') and file.startswith('all_assemblies_')]

	#Multiprocess
	multiprocess(convert_fa_to_csv, args.processes, fna_files)
	print('All files converted to .csv')

	csv_files=[target_dict+'/'+file \
	for file in os.listdir(target_dict)\
	 if file.endswith('.csv') and file.startswith('all_assemblies_')]

	if args.log==True:
		if len(csv_files)>0:
			log_lines.append('.fna files converted to .csv...\n')
		else:
			log_lines.append('.fna files not properly converted to .csv...FAILED\n')
		write_log()

	#concatenate files back into number of split files specified
	#through number of genomes to download 01-04-2022, this seems redundant
	#reconcatenate(fna_files)
	#reconcatenate(csv_files)


def reconcatenate(file_list):

	file_ending={file.split('.')[-1] for file in file_list}
	index=0
	new_files=[]
	print(f'Reconcatenating .{list(file_ending)[0]} files into {split_num} file(s)')
	for num in range(split_num):

		all_files=False
		new_index=index+int(len(file_list)/split_num)+\
		(len(file_list)%split_num>0)

		if new_index<=len(file_list):
			sub_list=file_list[index:new_index]
			index=new_index

		else:
			all_files=True
			sub_list=file_list[index:]

		cat_command='cat '
		for file in sub_list:
			cat_command=cat_command+f'{file} '
		cat_command=cat_command+f'> {os.path.dirname(file).rstrip("/")}/all_assemblies_{num}.{list(file_ending)[0]}.reconcat'
		subprocess.call(cat_command, shell=True)
		new_files.append(f'{os.path.dirname(file).rstrip("/")}/all_assemblies_{num}.{list(file_ending)[0]}.reconcat')
		
	#remove previous files
	for file in file_list:
		os.remove(file)
	
	#Rename new files to match old names
	for file in new_files:
		os.rename(file, file.replace(list(file_ending)[0]+'.reconcat', list(file_ending)[0]))

def add_taxonomy_lineages(summary_files):

	#Create dict with accession and taxid from assembly summary
	#identify latest assembly/plasmid summary file
	newest=sorted(summary_files)[-1]

	if args.update==True and len(summary_files)>1:
		previous=sorted(summary_files)[-2]

	#Process assembly summary
	if 'assembly' in newest:
		if not args.update==True:
			tax_dict={line.split('\t')[0]:line.split('\t')[5] for line in open(newest, 'r') if not line.startswith('#')}
		else:
			if len(summary_files)>1:
				old_lines=[line for line in open(previous, 'r')]
				tax_dict={line.split('\t')[0]:line.split('\t')[5] for line in open(newest, 'r') if not line.startswith('#') and not line in previous}
			else:
				tax_dict={line.split('\t')[0]:line.split('\t')[5] for line in open(newest, 'r') if not line.startswith('#')}

	#Process plasmid summary
	elif 'plasmid' in newest:
		if not args.update==True:
			tax_dict={line.split('\t')[0]:line.split('\t')[1] for line in open(newest, 'r') if not line.startswith('#')}
		else:
			if len(summary_files)>1:
				old_lines=[line for line in open(previous, 'r')]
				tax_dict={line.split('\t')[0]:line.split('\t')[1] for line in open(newest, 'r') if not line.startswith('#') and not line in previous}
			else:
				tax_dict={line.split('\t')[0]:line.split('\t')[1] for line in open(newest, 'r') if not line.startswith('#')}

	#Process local summary
	elif 'local_summary' in newest:
		# check if summary file contains any information
		with open(newest) as f:
			lines = f.readlines()
			# If we don't have any information, assign to cellular organism
			if lines[0].count('\t')<2:
				taxid = 131567
			else:
				taxid = False
		if not args.update==True:
			if not taxid:
				tax_dict={line.split('\t')[0]:line.split('\t')[1] for line in open(newest, 'r') if not line.startswith('#')}
			else:
				tax_dict={line.split('\t')[0]:taxid for line in open(newest, 'r') if not line.startswith('#')}
		else:
			if len(summary_files)>1:
				old_lines=[line for line in open(previous, 'r')]
				if not taxid:
					tax_dict={line.split('\t')[0]:line.split('\t')[1] for line in open(newest, 'r') if not line.startswith('#') and not line in previous}
				else:
					tax_dict={line.split('\t')[0]:taxid for line in open(newest, 'r') if not line.startswith('#') and not line in previous}
			else:
				tax_dict={line.split('\t')[0]:line.split('\t')[1] for line in open(newest, 'r') if not line.startswith('#')}

	#Create lineage dict
	lin_dict={}
	for key, value in tax_dict.items():
		lin_dict[key]={}
		lin_dict[key]['tax_id']=value
		lin_dict[key]['lineage']=''

	#Read in taxonomy db files
	nodes=open(args.target_directory.rstrip("/")+"/taxonomy/nodes.dmp", 'r')
	names=open(args.target_directory.rstrip("/")+"/taxonomy/names.dmp", 'r')

	node_dict={line.split('\t|\t')[0]:line.split('\t|\t')[1] for line in nodes}
	sci_names={line.split('\t')[0]:line.split('\t|\t')[1] for line in names \
	if 'scientific name' in line.split('\t|\t')[3]}

	#Create lineages for every node
	lineages=Manager().dict()
	for node, parent in node_dict.items():
		lineage=[]
		old_node=node
		lineage.append(node)
		while str(node)!='1':
			lineage.append(node_dict[node])
			node=node_dict[node]

		if '2' in lineage:
			dicti={}
			dicti[old_node]={}
			dicti[old_node]['lineage']=lineage
			dicti[old_node]['name']=sci_names[old_node]
			lineages.update(dicti)

	new_lineages=Manager().dict()
	lineage_list=[[key, value['lineage'], value['name']] for key, value in lineages.items()]

	multiprocess(assign_lineage_names, args.processes, lineage_list, new_lineages)

	if args.log==True:
		if len(new_lineages)>0:
			log_lines.append('Lineages assigned...\n')
		else:
			log_lines.append('Lineages could not be assigned...FAILED\n')
		write_log()

	return new_lineages


def assign_lineage_names(queue, new_lineages):

	line=queue.get()
	if line=='STOP':
		return
	
	#Create list of scientific names
	names=open(args.target_directory.rstrip("/")+"/taxonomy/names.dmp", 'r')

	sci_names={line.split('\t')[0]:line.split('\t|\t')[1] for line in names \
	if 'scientific name' in line.split('\t|\t')[3]}

	while True:
		tax={}
		tax[line[0]]={}
		named_lineage=[]
		for id in line[1]:
			try:
				sci_name=sci_names[id]
				named_lineage.append(sci_name)
			except:
				sci_name='unassigned'
				named_lineage.append(sci_name)

		try:
			tax[line[0]]['sci_lineage']=named_lineage
			tax[line[0]]['lineage']=line[1]
			tax[line[0]]['name']=line[2]
			new_lineages.update(tax)

		except Exception as e:
			print(f'ERROR: {e}, line: {line}')

		line=queue.get()
		if line=='STOP':
			return




def merge_files():

	print('Start merging...')
	#Create a dictionary of matching files, so that they can be merged
	file_dict={}
	for filename in os.listdir(args.target_directory.rstrip('/')+'/update_tmp'):
		file_dict[args.target_directory.rstrip('/')+'/update_tmp/'+filename]=''
	
		if os.path.isfile(args.target_directory.rstrip('/')+'/update_tmp/'+filename.replace('/update_tmp', '')):
			file_dict[args.target_directory.rstrip('/')+'/update_tmp/'+filename]=\
			args.target_directory.rstrip('/')+'/'+filename

	#now merge files into database directory
	for key, value in file_dict.items():
		if args.log==True:
			log_lines.append('Merging files post-update...\n')

		if not value=='':
			merge_command='cat %s >> %s' % (key, value)
			subprocess.call(merge_command, shell=True)

			if args.log==True:
				log_lines.append(f'{merge_command}...\n')

	#Merge assembly and plasmid summary files, then delete unmerged
	#Assembly summary files:
	if args.assemblies==True:
		asm_files=[f'{os.path.abspath(args.target_directory)}/{file}' for file in os.listdir(args.target_directory)\
			 if file.startswith('assembly_summary')]
		asm_sums=[line for file in asm_files for line in open(file, 'r')]

		with open(args.target_directory.rstrip('/')+'/assembly_summary.txt', 'w') as asm_file:
			for line in asm_sums:
				asm_file.write(line)

		for file in [file for file in os.listdir(args.target_directory) if file.startswith('assembly_summary')\
			and len(file.split('.'))>2]:
			os.remove(args.target_directory.rstrip('/')+'/'+file)

	#Plasmid summary files:
	if args.plasmids==True:
		plas_files=[f'{os.path.abspath(args.target_directory)}/{file}' for file in \
		os.listdir(args.target_directory) if file.startswith('plasmid_summary')]
		plas_sums=[line for file in plas_files for line in open(file, 'r')]

		with open(args.target_directory.rstrip('/')+'/plasmid_summary.txt.0', 'w') as plas_file:
			for line in plas_sums:
				plas_file.write(line)

		for file in [file for file in os.listdir(args.target_directory) if file.startswith('plasmid_summary')\
			and str(file.split('.')[2])!=str(0)]:
			os.remove(args.target_directory.rstrip('/')+'/'+file)

	#local summary files:
	if args.local!='False':
		loc_files=[f'{os.path.abspath(args.target_directory)}/{file}' for file in os.listdir(args.target_directory)\
			 if file.startswith('local_summary')]
		loc_sums=[line for file in loc_files for line in open(file, 'r')]

		with open(args.target_directory.rstrip('/')+'/local_summary.txt.0', 'w') as loc_file:
			for line in loc_sums:
				loc_file.write(line)

		for file in [file for file in os.listdir(args.target_directory) if file.startswith('local_summary')\
			and str(file.split('.')[2])!=str(0)]:
			os.remove(args.target_directory.rstrip('/')+'/'+file)

	if args.log==True:
		log_lines.append('Summary files merged')
		write_log()

	#Delete temporary directory
	shutil.rmtree(args.target_directory.rstrip('/')+'/update_tmp')
	print('''Old and new files merged, temporary update directory removed!\nDatabase updated!''')

def convert_fa_to_csv(queue):

	fa_file=queue.get()
	if fa_file=='STOP':
		return

	while True:

		#Read seqs into dict
		print('Reading in %s...' % fa_file)
		seq_dict={}

		for line in open(fa_file, 'r'):
			if line.startswith('>'):
				header=line.rstrip('\n')
				seq=''
			else:
				seq+=line.rstrip('\n')
				seq_dict[header]=seq
		
		print('Writing %s to .csv...' % fa_file)
		#write to .csv format
		with open(fa_file.replace('.fna', '.csv'), 'w') as outfile:
			for key, value in seq_dict.items():
				outfile.write(key+'\t'+value+'\n')
		outfile.close()

		seq_dict.clear()

		#Grab next file from the queue
		fa_file=queue.get()
		if fa_file=='STOP':
			return

def process_plasmids(queue, plas_dict):

	plasmid_file=queue.get()
	if plasmid_file=='STOP':
		return

	while True:
		print(f'Processing {plasmid_file}...')
		tmp_dict={}
		for line in open(args.target_directory.rstrip('/')+'/plasmids_tmp/'+plasmid_file, 'r'):
			if line.startswith('>'):
				header=line
				seq=''
			else:
				seq+=line
				tmp_dict[header]=seq

		plas_dict.update(tmp_dict)
		
		plasmid_file=queue.get()
		if plasmid_file=='STOP':
			return

				
def download_plasmids():

	#Create novel temporary directory for plasmids
	if not os.path.exists(args.target_directory.rstrip('/')+'/plasmids_tmp'):
		os.mkdir(args.target_directory.rstrip('/')+'/plasmids_tmp')

		print('downloading plasmid taxonomy...')
		download_tax='wget -r -P %s -nd -l 1 -A .catalog.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/' \
		% (args.target_directory.rstrip('/')+'/plasmids_tmp')
		subprocess.call(download_tax, shell=True)


		print('downloading plasmid sequences...')
		download_file2='wget -r -P %s -nd -A .genomic.fna.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/' \
		% (args.target_directory.rstrip('/')+'/plasmids_tmp')
		subprocess.call(download_file2, shell=True)

		#Unzip 
		unzip='gunzip %s*' % (args.target_directory.rstrip('/')+'/plasmids_tmp/')
		subprocess.call(unzip, shell=True)

	if args.log==True:
		if len(os.listdir(args.target_directory.rstrip('/')+'/plasmids_tmp'))>0:
			log_lines.append('Plasmids downloaded...\n')
		else:
			log_lines.append('No plasmids downloaded...FAILED\n')
		write_log()

	print('Processing plasmids, this may take some time...')

	#Multiprocess reading in of plasmids
	#Read into dict
	plasmid_files=[file for file in os.listdir(args.target_directory.rstrip('/')+\
	'/plasmids_tmp') if file.endswith('.genomic.fna')]

	plas_procs=args.processes

	plas_dict=Manager().dict()
	multiprocess(process_plasmids, plas_procs, plasmid_files, plas_dict)
	plasmid_dict=dict(plas_dict)

	if args.log==True:
		if len(plasmid_dict)>0:
			log_lines.append('Plasmid sequences read into memory...\n')
		else:
			log_lines.append('Plasmid sequences could not be read into memory...FAILED\n')
		write_log()

	#Determine number of plasmid_summary files
	sum_files=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
	if file.startswith('plasmid_summary')]

	if len(sum_files)==1:
		sum_num=1
	else:
		sum_num=len(sum_files)

	if args.taxa!='False':
		print('Writing plasmid summary file based on specified taxa...')
		with open(args.target_directory.rstrip('/')+'/plasmid_summary.txt'+'.'+str(sum_num), 'w') as outfile:
			try:
				for key, value in plasmid_dict.items():
					if any(taxon.lower() in key.lower() for taxon in args.taxa):
						outfile.write(str(key.split(' ')[0].lstrip('>'))+'\t'+' '.join(key.split(' ')[1:3])+'\n')
			except:
				
					outfile.write('')
					print('Something went wrong when writing plasmid summary...')
		outfile.close()

	if args.accessions!='False':
		
		print('Writing plasmid summary file from accession list...')
		with open(args.target_directory.rstrip('/')+'/plasmid_summary.txt'+'.'+str(sum_num), 'w') as outfile:
			try:
				for key, value in plasmid_dict.items():
					accessions=[line.rstrip('\n') for line in open(args.accessions, 'r')]
					if any(acc.lower() in key.lower() for acc in accessions):
						outfile.write(str(key.split(' ')[0].lstrip('>'))+'\t'+' '.join(key.split(' ')[1:3])+'\n')
			except:
					outfile.write('')
					print('Something went wrong when writing plasmid summary...')
		outfile.close()
	
	#add lineage to plasmid entries
	print('adding lineages to plasmids...')

	sum_files=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
	if file.startswith('plasmid_summary')]
	new_lineages=add_taxonomy_lineages(sum_files)

	#Go through refseq catalogue to extract all bacterial plasmids
	refseq_cat=[file for file in os.listdir(args.target_directory.rstrip('/')+'/plasmids_tmp') \
	if file.endswith('.catalog')]

	bacterial_plasmids=[line for line in open(args.target_directory.rstrip('/')+'/plasmids_tmp/'\
	+refseq_cat[0], 'r') if 'plasmid' in line.split('\t')[3]]

	#Read to dict: {accession:tax_id}
	refseq_dict={}
	for line in bacterial_plasmids:
		refseq_dict[line.split('\t')[2]]=line.split('\t')[0]
		
	#add taxon to plasmid entries
	new_plasmid_dict={}

	for key, value in plasmid_dict.items():
		seq=value
		taxon=' '.join(key.split(' ')[1:3])
		key2=key.lstrip('>').split(' ')[0]
		new_plasmid_dict[key2]={}
		new_plasmid_dict[key2]['seq']=value
		new_plasmid_dict[key2]['taxon']=' '.join(key.split(' ')[1:3])
		new_plasmid_dict[key2]['tax_id']=refseq_dict[key2]
		try:
			new_plasmid_dict[key2]['lineage']=new_lineages[refseq_dict[key2]]['sci_lineage']
		except:

			new_plasmid_dict[key2]['lineage']='unassigned'

	if args.log==True:
		if len(new_plasmid_dict)>0:
			log_lines.append('Lineages assigned to plasmids...\n')
		else:
			log_lines.append('Lineages could not be assigned to plasmids...FAILED\n')
		write_log()

	if args.update==True:
		
		#Connect to database and fetch species present in previous database version
		try:
			connection=sqlite3.connect(f'{os.path.abspath(args.target_directory)}/genview_database.db')
			if args.log==True:
				log_lines.append('Successfully connected to GEnView database...\n')
				write_log()

		except Exception as e:
			print(e)
			print('Could not connect to genview_database.db, exiting...\n')
			sys.exit()
		
		cursor=connection.cursor()

		query="""SELECT organism from genomes;"""
		cursor.execute(query)
		old_specs=cursor.fetchall()
		old_spec_set={spec[0] for spec in old_specs}

		if args.taxa!='False':
			for spec in old_spec_set:
				args.taxa.append(spec)

		if args.accessions=='False':
			print('comparing summary files...')
			#compare old and new assembly_summary files, get list of genomes that are in new but not in old
			summary_files=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
			if file.startswith('plasmid_summary')]

			newest=sorted(summary_files)[-1]

			if len(summary_files)>1:
				previous=sorted(summary_files)[-2]

				previous_plas_accs=[line.split('\t')[0] for line in open(previous, 'r')]
				new_plas_accs=[line.split('\t')[0] for line in open(newest, 'r') if not line.split('\t')[0] \
				in previous_plas_accs]
			else:
				
				new_plas_accs=[line.split('\t')[0] for line in open(newest, 'r')]

			print('%d novel plasmids identified!' % len(new_plas_accs))

		else:
			
			summary_files=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
			if file.startswith('plasmid_summary')]
			newest=sorted(summary_files)[-1]
			new_plas_accs=[line.split('\t')[0] for line in open(newest, 'r') if not line.startswith('#')]

		if not args.accessions=='False':
			genome_accessions=[line.rstrip('\n') for line in open(args.accessions)]
			hits=0
			for key, value in new_plasmid_dict.items():
				#Filter by taxonomy:
				if key in genome_accessions and key in new_plas_accs:
					hits+=1
					with open(args.target_directory.rstrip('/')+'/update_tmp/genomes/'\
					+key+'_genomic.fna', 'w') as outfile:
						outfile.write('>'+key+' plasmid\n'+value['seq'])
					outfile.close()

			if args.log==True:
				if hits>0:
					log_lines.append('Plasmids written to file...\n')
				else:
					log_lines.append('No new plasmids within accession list...\n')
				write_log()

		else:

			if not args.taxa[0]=='all':
				hits=0
				for key, value in new_plasmid_dict.items():
					#Filter by taxonomy:
					if any(taxon.lower() in ' '.join(value['lineage']).lower() for taxon\
					in args.taxa) and key in new_plas_accs:

						hits+=1
						with open(args.target_directory.rstrip('/')+'/update_tmp/genomes/'\
						+key+'_genomic.fna', 'w') as outfile:
							outfile.write('>'+key+' plasmid\n'+value['seq'])
						outfile.close()

				if args.log==True:
					if hits>0:
						log_lines.append('Plasmids written to file...\n')
					else:
						log_lines.append('No new plasmids found for specified taxa...\n')
					write_log()
				
				if hits>=1:
					print('Plasmids written to file!')
				else:
					print('No plasmids found for the searched taxa!')
					if not args.assemblies==True:
						sys.exit()
			else:
				hits=0
				for key, value in new_plasmid_dict.items():
					if key in new_plas_accs:
						hits+=1
						with open(args.target_directory.rstrip('/')+'/update_tmp/genomes/'\
						+key+'_genomic.fna', 'w') as outfile:
							outfile.write('>'+key+' plasmid\n'+value['seq'])
						outfile.close()					

				if args.log==True:
					if hits>0:
						log_lines.append('Plasmids written to file...\n')
					else:
						log_lines.append('No new plasmids for specified taxa...\n')
					write_log()

	else:

		
		if not args.accessions=='False':
			genome_accessions=[line.rstrip('\n') for line in open(args.accessions)]
			hits=0
			for key, value in new_plasmid_dict.items():
				#Filter by taxonomy:
				if key in genome_accessions and key in new_plas_accs:
					hits+=1
					with open(args.target_directory.rstrip('/')+'/update_tmp/genomes/'\
					+key+'_genomic.fna', 'w') as outfile:
						outfile.write('>'+key+' plasmid\n'+value['seq'])
					outfile.close()
				
			if args.log==True:
				if hits>0:
					log_lines.append('Plasmids written to file...\n')
				else:
					log_lines.append('No new plasmids for specified taxa...\n')
				write_log()

		else:
			if not args.taxa[0]=='all':
				hits=0
				for key, value in new_plasmid_dict.items():
					#Filter by taxonomy:
					if any(taxon.lower() in ' '.join(value['lineage']).lower() for taxon\
					in args.taxa):

						hits+=1
						with open(args.target_directory.rstrip('/')+'/genomes/'\
						+key+'_genomic.fna', 'w') as outfile:
							outfile.write('>'+key+' plasmid\n'+value['seq'])
						outfile.close()
					
				if args.log==True:
					if hits>0:
						log_lines.append('Plasmids written to file...\n')
					else:
						log_lines.append('No new plasmids for specified taxa...\n')
					write_log()

				
				if hits>=1:
					print('Plasmids written to file!')
				else:
					print('No plasmids found for the searched taxa!')
					if not args.assemblies==True:
						sys.exit()
			else:
				for key, value in new_plasmid_dict.items():
					with open(args.target_directory.rstrip('/')+'/genomes/'\
					+key+'_genomic.fna', 'w') as outfile:
						outfile.write('>'+key+' plasmid\n'+value['seq'])
					outfile.close()

				if args.log==True:
					if len(new_plasmid_dict)>0:
						log_lines.append('Plasmids written to file...\n')
					else:
						log_lines.append('No new plasmids for specified taxa...\n')
					write_log()
	#Remove plasmids_tmp
	shutil.rmtree(args.target_directory.rstrip('/')+'/plasmids_tmp')

def update():

	#Check if genview database exists
	if not os.path.exists(f'{os.path.abspath(args.target_directory)}/genview_database.db'):
		print('No genview database found in the target directory, exiting...\n')
		sys.exit()

	download_uniprot()

	print('Update: creating temporary download directory...')
	#Create temporary directory for downloading novel assemblies
	if not os.path.exists(args.target_directory.rstrip('/')+'/'+'update_tmp'):
		os.mkdir(args.target_directory.rstrip('/')+'/'+'update_tmp')


	if args.log==True:
		if os.path.isdir(args.target_directory.rstrip('/')+'/'+'update_tmp'):
			log_lines.append('Update directory created...\n')
		else:
			log_lines.append('Update directory not found...FAILED\n')
		write_log()

	#Download assembly summary files
	if args.assemblies==True:

		print('Update: downloading new assembly summary...')
		download_file='wget -P %s ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt' \
		% (args.target_directory)
		#Disable for debugging
		subprocess.call(download_file, shell=True)

		#compare old and new assembly_summary files, get list of genomes that are in new but not in old
		summary_files=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
		if file.startswith('assembly_summary.')]

		newest=sorted(summary_files)[-1]
		
		if args.log==True:
			if os.path.exists(os.path.abspath(args.target_directory)+'/assembly_summary.txt'):
				log_lines.append('Assembly summary downloaded...\n')
			else:
				log_lines.append('Assembly summary not found...FAILED\n')
			write_log()
		
		#In case a database is updated with assemblies for the first time:
		if len(summary_files)>1:
			previous=sorted(summary_files)[-2]
		
		#Rewrite new assembly summary file such that it only contains genome accessions from specified species
		if not args.taxa=='False':
			if not args.taxa[0]=='all':
				tax_asms=[]
				taxon_lines=[line for line in open(newest, 'r') if not line.startswith('#')]

				for line in taxon_lines:
					if any(taxon.lower() in line.lower() for taxon in args.taxa):
						tax_asms.append(line)

				with open(newest, 'w') as f:
					for line in tax_asms:
						f.write(line)
				f.close()

		if not args.accessions=='False':
			accessions=[line.rstrip('\n') for line in open(args.accessions, 'r')]
			asm_sum_lines=[line for line in open(newest, 'r') if any(acc in line for acc in accessions)]

			with open(newest, 'w') as f:
				for line in asm_sum_lines:
					f.write(line)
			f.close()

		if args.log==True:
			if os.path.exists(newest):
				log_lines.append('New assembly summary rewritten...\n')
			else:
				log_lines.append('No new assembly summary file found...FAILED\n')
			write_log()

		if args.accessions=='False':
			if len(summary_files)>1:
				print('comparing assembly summary files...')
				previous_asms=[line.split('\t')[0] for line in open(previous, 'r') if not line.startswith('#')]
				new_genome_urls=[line.split('\t')[19] for line in open(newest, 'r') if not line.split('\t')[0] \
				in previous_asms and not line.startswith('#')]
				
				print('%d new assemblies found!' % len(new_genome_urls)) 
			else:
				new_genome_urls=[line.split('\t')[19] for line in open(newest, 'r')]
		else:
			if len(summary_files)>1:

				previous_asms=[line.split('\t')[0] for line in open(previous, 'r') if not line.startswith('#')]
				new_genome_urls=[line.split('\t')[19] for line in open(newest, 'r') if not line.startswith('#') and line.split('\t')[0] in [line.rstrip('\n') for line in open(args.accessions, 'r')] and not line.split('\t')[0] in previous_asms]

			else:

				new_genome_urls=[line.split('\t')[19] for line in open(newest, 'r') if not line.startswith('#') and line.split('\t')[0] in args.accessions]
	
		if args.log==True:
			if len(new_genome_urls)>=1:
				log_lines.append('New genomes identified...\n')
			else:
				log_lines.append('No new genomes identified...FAILED(?)\n')
			write_log()

		#Create dictionary with assembly accession as key, lineages as value
		new_lineages=add_taxonomy_lineages(summary_files)

		asm_dict={}
		excepts=0
		new_lines=[line for line in open(newest,'r') if not line.startswith('#')]
		newlines_inc=[line for line in new_lines if line.split('\t')[19] in new_genome_urls]
		for line in newlines_inc:
			
			asm_dict[line.split('\t')[0]]={}
			asm_dict[line.split('\t')[0]]['tax_id']=line.split('\t')[5]
			asm_dict[line.split('\t')[0]]['url']=line.split('\t')[19]

			if args.accessions=='False':
				
				print('Assigning lineages...')
				try:
					asm_dict[line.split('\t')[0]]['sci_lineage']=new_lineages[line.split('\t')[5]]['sci_lineage']
					asm_dict[line.split('\t')[0]]['lineage']=new_lineages[line.split('\t')[5]]['lineage']
				except:
					excepts+=1
					asm_dict[line.split('\t')[0]]['sci_lineage']='unassigned'
					asm_dict[line.split('\t')[0]]['lineage']='unassigned'

				print(f'lineages assigned, {excepts} of {len(asm_dict)} could not be assigned')

				#Connect to database and fetch species present in previous database version
				connection=sqlite3.connect(args.target_directory.rstrip('/')+'/genview_database.db')
				cursor=connection.cursor()

				query="""SELECT organism from genomes;"""
				cursor.execute(query)
				old_specs=cursor.fetchall()
				old_spec_set={spec[0] for spec in old_specs}

				for spec in old_spec_set:
					args.taxa.append(spec)

		#Now multiprocess the download of these new genomes into temporary folder
		#disable for now to make sure not all genomes have to be downloaded again
		#TODO add exact line lineage in summary file

		if not args.accessions=='False':
			genome_accessions=[line.rstrip('\n').lower() for line in open(args.accessions, 'r')]
			genome_urls=[asm_dict[key]['url'] for key, value in asm_dict.items() if key.lower() in genome_accessions]
		
		else:
			if not args.taxa[0]=='all':

				genome_urls=[asm_dict[key]['url'] for key, value in asm_dict.items() if any(taxon.lower() in \
				' '.join(asm_dict[key]['sci_lineage']).lower() for taxon in args.taxa)]
			else:
				
				genome_urls=[asm_dict[key]['url'] for key, value in asm_dict.items()]
		
		if args.log==True:
			if len(new_genome_urls)>=1:
				log_lines.append('New genomes identified...\n')
			else:
				log_lines.append('No new genomes identified...FAILED(?)\n')
			write_log()
		
		if not os.path.exists(f'{os.path.abspath(args.target_directory)}/update_tmp/genomes'):
			os.mkdir(f'{os.path.abspath(args.target_directory)}/update_tmp/genomes')
	
		multiprocess(download_new, args.processes, genome_urls)

	elif args.local!='False':

		#Create local summary file
		create_simple_summary_file()
		
		#find local assembly files
		local_asms=[f'{os.path.abspath(args.target_directory)}/{file}' for file in \
		os.listdir(args.target_directory) if file.startswith('local_summary')]
		
		if len(local_asms)>1:
		
			newest=sorted([line for line in open(local_asms[-1], 'r')])
			previous=sorted([line for line in open(local_asms[-2], 'r')])

			#Set 'genome_urls' to number of new genomes in specified local path
			genome_urls=[line for line in newest if not line in previous]
		else:
			genome_urls=[line for line in open(local_asms[0], 'r')]

	if args.plasmids==True:	

		if not os.path.exists(f'{os.path.abspath(args.target_directory)}/update_tmp/genomes'):
			os.mkdir(f'{os.path.abspath(args.target_directory)}/update_tmp/genomes')

		#Download new plasmids
		download_plasmids()
		if 'genome_urls' in locals():
			if args.taxa!='False':
				genome_urls.extend([line for line in open(f'{os.path.abspath(args.target_directory)}/plasmid_summary.txt.0') if any(taxon for taxon in args.taxa)])
			else:
				genome_urls.extend([line for line in open(f'{os.path.abspath(args.target_directory)}/plasmid_summary.txt.0')])
		else:
			if args.taxa!='False':
				genome_urls=[line for line in open(f'{os.path.abspath(args.target_directory)}/plasmid_summary.txt.0') if any(taxon for taxon in args.taxa)]
			else:
				genome_urls=[line for line in open(f'{os.path.abspath(args.target_directory)}/plasmid_summary.txt.0')]

	if len(genome_urls)==0:
		print('No new genomes identified, exiting...')
		#Remove newest summary file(s) and in case of update, update directory
		if args.update==True:
			shutil.rmtree(f'{os.path.abspath(args.target_directory)}/update_tmp')
		
			if args.assemblies!='False':
				del_asms=[f'{os.path.abspath(args.target_directory)}/{file}' for file in os.listdir\
				(args.target_directory) if 'assembly_summary' in file]
				os.remove(sorted(del_asms)[-1])	

			if args.plasmids!='False':
				del_plas=[f'{os.path.abspath(args.target_directory)}/{file}' for file in os.listdir\
				(args.target_directory) if 'plasmid_summary' in file]
				os.remove(sorted(del_plas)[-1])	

			if args.local!='False':
				del_loc=[f'{os.path.abspath(args.target_directory)}/{file}' for file in os.listdir\
				(args.target_directory) if 'local_summary' in file]
				os.remove(sorted(del_loc)[-1])	
		sys.exit()
		
	#make split_num available in other functions	
	global split_num

	if 0<=len(genome_urls)<1000:	
		split_num=2
	elif 1000<len(genome_urls)<10000:	
		split_num=5
	elif 10000<len(genome_urls)<100000:	
		split_num=20
	else:
		split_num=200

	#concatenate new genomes to fasta and csv, then split into several files
	#(corresponding to number of split files for creating the original db)
	concatenate_and_split()

	print('Update: collecting fasta files...')
	#Create list of files containing assemblies in fasta format
	fa_files=[args.target_directory.rstrip('/')+'/'+'update_tmp/'+fa_file for fa_file \
	in os.listdir(args.target_directory.rstrip('/')+'/update_tmp')\
	if fa_file.startswith('all_assemblies_')\
	and fa_file.endswith('.fna') and not fa_file.startswith('flanking') \
	and not '_orfs' in fa_file]
	
	#make sure diamond database format is correct:
	args.database=reformat()

	#Define number of processes and empty list for processes
	multiprocess(annotate, args.processes, fa_files)
	#Do the same again for the 'create_database' function

	if args.log==True:

		#move diamond.log from current wd to target directory
		if not os.path.exists(f'{args.target_directory}/update_tmp/diamond.log'):
			mv_command=f'mv {os.getcwd()}/diamond.log {os.path.abspath(args.target_directory)}/update_tmp'
		if os.path.exists(f'{os.path.abspath(args.target_directory)}/diamond.log') and not f'{os.path.abspath(args.target_directory)}'==os.getcwd():
			mv_command=f'cat {os.getcwd()}/diamond.log {os.path.abspath(args.target_directory)}/update_tmp/diamond.log > {os.path.abspath(args.target_directory)}/update_tmp/diamond.log'
		subprocess.call(mv_command, shell=True)

		dmnd_log=[line for line in open(f'{os.path.abspath(args.target_directory)}/update_tmp/diamond.log', 'r') \
		if 'queries aligned' in line]
		if len(dmnd_log)==len(fa_files):
			log_lines.append(f'diamond blastx run successfull...'+'\n')
		else:
			log_lines.append(f'diamond blastx run unsuccessfull...FAILED'+'\n')
		write_log()

	if args.log==True:
		#Check presence and size of annotation files
		anno_log=[file for file in os.listdir(args.target_directory.rstrip('/')+'/update_tmp')\
			if file.endswith('_annotated.csv')]
		if len(anno_log)>0:
			log_lines.append('Annotation files generated...\n')
		else:
			log_lines.append('No annotation files could be found...FAILED\n')
		write_log()

	#Create list with files that have no corresponding flanking region file yet
	no_flank_files=[element for element in fa_files if not os.path.exists(element+'_flanking_regions.csv')]

	if len(no_flank_files)>0:

		multiprocess(create_db, args.processes, no_flank_files)

	print('Update: collecting flanking regions...')
	flank_files=[args.target_directory.rstrip('/')+'/update_tmp/'+fa_file for fa_file in \
	os.listdir(args.target_directory.rstrip('/')+'/update_tmp')\
	if fa_file.startswith('all_assemblies_')\
	and fa_file.endswith('.fna_flanking_regions.csv')]


	if args.log==True:
		#Check presence and size of annotation files
		if len(flank_files)>0:
			log_lines.append('flanking region files generated...\n')
		else:
			log_lines.append('No flanking region files could be found...FAILED\n')
		write_log()

	#summarize all flanking regions into one temporary file
	lines=[line for file in flank_files for line in open(file, 'r')]


	#To update, all IDs of genes already in the database have to be fetched
	previous_ids=[int(line.split('\t')[-1]) for line in open(args.target_directory\
	.rstrip('/')+'/all_annos.fna_tmp')]

	print('Update: Writing temporary files...')
	with open(args.target_directory.rstrip('/')+'/update_tmp/all_flanks.csv_tmp', 'w') as outfile:
		#Now next ID is assigned to last gene id before update+1
		i=int(sorted(previous_ids)[-1])
		for line in sorted(lines):
			#append ID for identification of env genes
			asm_name=line.split('\t')[-4]
			if asm_name.endswith('.fna'):
				asm_name=asm_name.replace('.fna', '.1')

			i+=1

			#Write in csv format to outfile
			outfile.write(str(i)+'\t'+line.split('\t')[-6]+'\n')
	outfile.close()

	#after assigning id, split the file into specified number of smaller fasta files again
	splitted_files=[file for file in os.listdir(args.target_directory.rstrip('/')+'/update_tmp/')\
	if file.startswith('all_assemblies_') and file.endswith('.csv') and len(file.split('.'))==2]
	split_number=len(splitted_files)
	csv_lines=[line for line in open(args.target_directory.rstrip('/')+'/update_tmp/all_flanks.csv_tmp', 'r')]

	line_count=0
	file_count=1
	outfile=open(args.target_directory.rstrip('/')+'/update_tmp/'+\
	'flanking_regions_'+str(file_count)+'.fna', 'w')
	for line in csv_lines:
		line_count+=1
		if line_count<=round((len(csv_lines)/split_number)+0.5, 0):
			outfile.write('>'+line.split('\t')[0]+'\n'+line.split('\t')[-1])
		else:
			file_count+=1
			line_count=0
			outfile.close()
			outfile=open(args.target_directory.rstrip('/')+'/update_tmp/'+\
			'flanking_regions_'+str(file_count)+'.fna', 'w')
			outfile.write('>'+line.split('\t')[0]+'\n'+line.split('\t')[-1])
	outfile.close()

	#Create list with flanking regions file paths
	flanking_fa=[args.target_directory.rstrip('/')+'/update_tmp/'+file for file in os.listdir(args.target_directory\
	.rstrip('/')+'/update_tmp') \
	if file.startswith('flanking_regions') and file.endswith('.fna') and not '_orfs' in file]

	#Create a multiprocessing queue from the flanking region files
	multiprocess(run_prodigal, args.processes, flanking_fa)

	if args.log==True:
		orf_fna_log=[file for file in os.listdir(args.target_directory.rstrip('/')+'/update_tmp')\
		if file.endswith('_orfs.fna')]
		
		orf_gff_log=[file for file in os.listdir(args.target_directory.rstrip('/')+'/update_tmp')\
		if file.endswith('_orfs.gff')]
		if len(orf_fna_log)>0 and len(orf_gff_log)>0:
			log_lines.append('Orfs predicted...\n')
		else:
			log_lines.append('No orf files present, something might have gone wrong when running prodigal...FAILED\n')
		write_log()	

	#concatenate all orf files into one for clustering
	concat_command='cat %s/flanking_regions_*_orfs.fna > %s/all_orfs.fna' % \
	(args.target_directory.rstrip('/')+'/update_tmp', args.target_directory.rstrip('/')+'/update_tmp')
	subprocess.call(concat_command, shell=True)

	#Call clustering function which returns cluster dict
	clust_dict=cluster_orfs(args.target_directory.rstrip('/')+'/update_tmp'+'/all_orfs.fna', \
	args.target_directory.rstrip('/')+'/update_tmp')

	#Split resulting centroid file into several for annotation
	print('Update: splitting orf file for annotation...')
	centr_dict={}
	for line in open(args.target_directory.rstrip('/')+'/update_tmp'+'/orfs_clustered.fna'):
		if line.startswith('>'):
			header=line
			seq=''
		else:
			seq+=line
			centr_dict[header]=seq

	outfiles=split_num

	key_count=0
	file_count=1
	outfile=open(args.target_directory.rstrip('/')+'/update_tmp'+'/'+\
	'split_'+str(file_count)+'_orfs.fna', 'w')
	for key, value in centr_dict.items():
		key_count+=1
		if key_count<=round((len(centr_dict)/split_num)+0.5, 0):
			outfile.write(key+value)
		else:
			file_count+=1
			key_count=0
			outfile.close()
			outfile=open(args.target_directory.rstrip('/')+'/update_tmp'+'/'+\
			'split_'+str(file_count)+'_orfs.fna', 'w')
			outfile.write(key+value)
	outfile.close()
	
	if args.log==True:
		orf_split_log=[file for file in os.listdir(args.target_directory.rstrip('/')+'/update_tmp')\
		if 'split_' in file and file.endswith('_orfs.fna')]
		if len(orf_split_log)>0:
			log_lines.append('Orf files split for annotation...\n')
		else:
			log_lines.append('Orf files could not be split for annotation...FAILED\n')	
		write_log()


	#Now annotate orfs
	orf_files=[args.target_directory.rstrip('/')+'/update_tmp'+'/'+file for file in \
	os.listdir(args.target_directory.rstrip('/')+'/update_tmp') if file.endswith('_orfs.fna') \
	and file.startswith('split_')]

	multiprocess(annotate_orfs, args.processes, orf_files)

	if args.log==True:

		#move diamond.log from current wd to target directory
		if not os.path.exists(f'{args.target_directory}/update_tmp/diamond.log'):
			mv_command=f'mv {os.getcwd()}/diamond.log {os.path.abspath(args.target_directory)}/update_tmp'
		if os.path.exists(f'{os.path.abspath(args.target_directory)}/diamond.log') and not f'{os.path.abspath(args.target_directory)}'==os.getcwd():
			mv_command=f'cat {os.getcwd()}/diamond.log {os.path.abspath(args.target_directory)}/update_tmp/diamond.log > {os.path.abspath(args.target_directory)}/update_tmp/diamond.log'
		subprocess.call(mv_command, shell=True)

		if os.path.exists(f'{os.getcwd()}/diamond.log'):
			os.remove(f'{os.getcwd()}/diamond.log')

		#now rename diamond log
		if os.path.exists(f'{os.path.abspath(args.target_directory)}/dmnd.log'):
			mv2=f'mv {os.path.abspath(args.target_directory)}/dmnd.log {os.path.abspath(args.target_directory)}/diamond.log'
			subprocess.call(mv2, shell=True)

		dmnd_log=[line for line in open(f'{os.path.abspath(args.target_directory)}/update_tmp/diamond.log', 'r') \
		if 'queries aligned' in line]
		if not args.update==True:
			if len(dmnd_log)/2==len(fa_files):
				log_lines.append(f'diamond blastp run successfull...'+'\n')
			else:
				log_lines.append(f'diamond blastp run unsuccessfull...FAILED'+'\n')
		else:
			
			if len(dmnd_log)>1:
				log_lines.append(f'diamond blastp run successfull...'+'\n')
			else:
				log_lines.append(f'diamond blastp run unsuccessfull...FAILED'+'\n')
		write_log()

	#Create a temporary summary file containing the line id
	with open(args.target_directory.rstrip('/')+'/update_tmp/all_annos.fna_tmp', 'w') as outfile:
		i=int(sorted(previous_ids)[-1])
		for line in sorted(lines):
			i+=1
			outfile.write(line.rstrip('\n')+'\t'+str(i)+'\n')
	outfile.close()

	all_orfs_gff=[args.target_directory.rstrip('/')+'/update_tmp/'+file for file in os.listdir(args.target_directory\
	.rstrip('/')+'/update_tmp') if file.endswith('_orfs.gff')]	

	all_orf_annos=[args.target_directory.rstrip('/')+'/update_tmp/'+file for file in os.listdir(args.target_directory\
	.rstrip('/')+'/update_tmp') if file.endswith('_orfs_annotated.csv')]

	all_IS_annos=[args.target_directory.rstrip('/')+'/update_tmp/'+file for file in os.listdir(args.target_directory) \
	if file.endswith('_orfs_ISannotated.csv')]

	env_dict=create_env_dict(all_orfs_gff, all_orf_annos, clust_dict, all_IS_annos)

	to_sql_db(env_dict, args.target_directory.rstrip('/')+'/update_tmp/all_annos.fna_tmp', \
	args.target_directory)

	#run integron finder and integrate into database
	#Create list with flanking region files
	flanks=[args.target_directory.rstrip('/')+'/update_tmp/'+file for file in os.listdir(\
	args.target_directory.rstrip('/')+'/update_tmp')\
	if file.startswith('flanking_regions_') and \
	file.endswith('.fna') and not '_orfs' in file]

	#Create temporary directory
	if args.integron_finder=='True':
		if not os.path.exists(args.target_directory.rstrip('/')+'/update_tmp/integrons_tmp'):
			os.mkdir(args.target_directory.rstrip('/')+'/update_tmp/integrons_tmp')

		dir=args.target_directory.rstrip('/')+'/update_tmp/integrons_tmp'

		multiprocess(annotate_integrons, args.processes, flanks, dir)	
		integrons_to_db()

	#Create arg/tnp association table
	transposon_table()
	merge_files()

	if not args.save_tmps==True:
		remove_tmps()

#Annotate sequences. Queue is a list like object containing 'tasks' (in this case assembly names) for each process to grab, emptying the queue in the process
def annotate(queue):

	#Grab an item from the queue. If item is 'STOP', end the process
	fa_file=queue.get()
	if fa_file=='STOP':
		return

	ending=fa_file.split('.')[-1]
	while True:
		try:
			if not os.path.exists(fa_file.replace(ending, '_annotated.csv')):
				print('Annotating %s...' % os.path.basename(fa_file))
				#Use diamond to annotate genes with local thresholds. Accept only one hit per target
			if args.processes*split_num>multiprocessing.cpu_count()*2:	
				threads=int((multiprocessing.cpu_count()*2)/args.processes)
			else:
				threads=args.processes

			if args.long_reads==True:
				frameshift='-F 1'
			else:
				frameshift=''

			diamond_blast=f'diamond blastx -p {threads} --quiet -d {args.database} -q {fa_file} -o {fa_file.replace(ending, "_annotated.csv")} --id {args.identity} --more-sensitive --quiet --top 100 --masking 0 {frameshift} --subject-cover {args.subject_coverage} --log -f 6 qseqid sseqid stitle pident qstart qend qlen slen length score qseq qframe qtitle'
			subprocess.call(diamond_blast, shell=True)

			#Sort by queryID and start position
			dmnd_rawout=pd.read_csv(fa_file.replace(ending, '_annotated.csv'), sep='\t', \
			names=['qseqid', 'sseqid', 'stitle', 'pident', 'qstart', 'qend', 'qlen', 'slen', 'length', 'score', 'qseq', 'qframe', 'qtitle'])
			dmnd_rawout.sort_values(by=['qseqid', 'qstart']).to_csv(fa_file.replace(ending, '_annotated.csv'), sep='\t', \
			index=False, header=False)
			#Assign same start and end position if several gene variants hit the same locus to be able to group by exact position and filter best hit
			with open(fa_file.replace(ending, '_annotated_rep.csv'), 'w') as rawout:
				line_number=0
				for line in open(fa_file.replace(ending, '_annotated.csv'), 'r'):
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
								if int(line.split('\t')[4]) in range(int(seqstart)-100, int(seqstart)+100) and int(line.split('\t')[5]) in range(int(seqend)-100, int(seqend)+100):
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
								if int(line.split('\t')[5]) in range(int(seqstart)-100, int(seqstart)+100) and int(line.split('\t')[4]) in range(int(seqend)-100, int(seqend)+100):
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
			dmnd_hits=pd.read_csv(fa_file.replace(ending, '_annotated_rep.csv'), sep='\t', names=['qseqid', 'sseqid', 'stitle', 'pident', 'qstart', 'qend', 'qlen', 'slen', 'length', 'score', 'qseq', 'qframe', 'qtitle'])
			best_hits=dmnd_hits.loc[dmnd_hits.groupby(['qseqid', 'qstart', 'qend'])['score'].idxmax()]
			best_hits.to_csv(fa_file.replace(ending, '_annotated.csv'), sep='\t', index=False, header=False)

					
		except Exception as e:
			print('EXCEPTION: %s, line %r' % (e, sys.exc_info()[-1].tb_lineno))
		fa_file=queue.get()
		if fa_file=='STOP':
			return

#Create database, which is a list of dictionaries containing sequence and metainformation about a genome
def create_db(queue):

	#Grab items from queue, which contains paths to assembly directories
	fa_file=queue.get()
	if fa_file=='STOP':
		return
	print('processing %s' % fa_file)

	if args.assemblies==True:
		#Parse assembly summary files for organism information
		org_dict={}

		#Identify most recent assembly summary file
		summary_files=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
		if file.startswith('assembly_summary')]

		newest=sorted(summary_files)[-1]

		for line in open(newest, 'r'):
			if not line.startswith('#'):
				org_dict[line.split('\t')[0]]=line.split('\t')[7]

		tax_lines=[line for line in open(newest, 'r') if not line.startswith('#')]

	if args.plasmids==True:

		if not 'org_dict' in locals():
			org_dict={}
		#Identify most recent plasmid_summary file
		summary_files_plas=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
		if file.startswith('plasmid_summary')]

		newest_plas=sorted(summary_files_plas)[-1]

		for line in open(newest_plas, 'r'):
			if len(line)>1:
				org_dict[line.split('\t')[0]]=line.split('\t')[1].rstrip('\n')+' plasmid'	

	if args.local!='False':

		if not 'org_dict' in locals():
			org_dict={}

		#Identify most recent local_summary file
		summary_files_cust=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
		if file.startswith('local_summary')]

		newest_cust=sorted(summary_files_cust)[-1]

		for line in open(newest_cust, 'r'):
			if len(line)>1:
				if args.long_reads=='False':
					#Check how much information we have got in the local summary
					num_infos = len(line.split('\t'))
					print('num infos {}'.format(num_infos))
					if num_infos>3:
						info = line.split('\t')[2].rstrip('\n') + ' ' + line.split('\t')[3].rstrip('\n')
					elif num_infos > 1:
						info = line.split('\t')[2].rstrip('\n')
					else:
						info = 'Unknown'
					org_dict[line.split('\t')[0]]=info
				else:
					info=line.rstrip('\n')
					org_dict[line]=info
					

		tax_lines_cust=[line for line in open(newest_cust, 'r') if not line.startswith('#')]

	while True:
		
		ending=fa_file.split('.')[-1]

		#Create set of contig accessions that contain annotated gene
		arg_assemblies={line.split('\t')[0] for line in open(fa_file.replace(ending, '_annotated.csv'), 'r')}

		seq_dict={}
		#Parse genome .csv files
		for line in open(fa_file.replace(ending, 'csv'), 'r'):
			acc = line.split('\t')[0].split(' ')[0].lstrip('>')
			if acc in arg_assemblies:
				#Create a dictionary containing just the contigs with annotated genes
				seq_dict[acc]=line.split('\t')[1].rstrip('\n')

		#Go through annotated files and extract info
					
		copy_list=[]
		nrows=0
		matches=0
		updated_lines=[]

		#Exclude genomes that have no annotated genes through breaking the
		#loop upon detection of empty annotation file
		if len(open(fa_file.replace(ending, '_annotated.csv'), 'r')
.readlines())<1:
			break
		
		print('Parsing annotations...')
		for line in open(fa_file.replace(ending, '_annotated.csv'), 'r'):
			nrows+=1
			gene_name=line.split('\t')[2].split('|')[3].split('[')[0].strip()+'__'+str(nrows)
			card_id_copy=line.split('\t')[1].split('|')[1]+'__'+str(nrows)
			

		#Extract upstream and downstream contexts of annotated genes. Identify right contig
		#based on diamond output that was saved to the genome
			if line.split('\t')[0] in seq_dict:
				contig_acc=line.split('\t')[0]
				frame=line.split('\t')[11]
				align_start=line.split('\t')[4]
				align_end=line.split('\t')[5]
				genome_acc=line.split('\t')[12].split('__')[-1].rstrip('\n')
				flank_length = args.flanking_length

				if genome_acc.endswith('.fna'):
					genome_acc=genome_acc.replace('.fna', '.1')

				if frame.startswith('-'):

					#Determine upstr and downstr length for - strand
					down_len=len(seq_dict[contig_acc][:int(align_end)])
					up_len=len(seq_dict[contig_acc][int(align_start):])

					#now extract seqs
					#1.) Both upstream and downstream are longer than 10kb
					if down_len>=flank_length and up_len>=flank_length:
						flanking_seq=\
						seq_dict[contig_acc][int(align_end)-flank_length:\
						int(align_start)+flank_length]
						new_uplen=flank_length
						new_downlen=flank_length
					#2.) Only upstream is longer than 10kb
					elif down_len<flank_length and up_len>=flank_length:
						flanking_seq=\
						seq_dict[contig_acc][:int(align_start)+flank_length]
						new_uplen=flank_length
						new_downlen=down_len
					#3.) Only downstream is longer than 10kb
					elif down_len>=flank_length and up_len<flank_length:
						flanking_seq=\
						seq_dict[contig_acc][int(align_end)-flank_length:]
						new_uplen=up_len
						new_downlen=flank_length
					#4.) Both are shorter than 10kb
					else:
						flanking_seq=\
						seq_dict[contig_acc]
						new_uplen=up_len
						new_downlen=down_len

					#Calculate reverse complement for - strand seqs here, so that all
					#genes are oriented the same way for the calculation
					flanking_seq_rev=reverse_complement(str(flanking_seq))
					flanking_seq=str(flanking_seq_rev)
					print('Sequence reversed!')

				#Now do the same thing for plus strand seqs
				else:

					#Determine upstr and downstr length for + strand
					down_len=len(seq_dict[contig_acc][int(align_end):])
					up_len=len(seq_dict[contig_acc][:int(align_start)])

					#now extract seqs
					#1.) Both upstream and downstream are longer than 10kb
					if down_len>=flank_length and up_len>=flank_length:
						flanking_seq=\
						seq_dict[contig_acc][int(align_start)-flank_length:\
						int(align_end)+flank_length]
						new_uplen=flank_length
						new_downlen=flank_length
					#2.) Only upstream is longer than 10kb
					elif down_len<flank_length and up_len>=flank_length:
						flanking_seq=\
						seq_dict[contig_acc][int(align_start)-flank_length:]
						new_uplen=flank_length
						new_downlen=down_len
					#3.) Only downstream is longer than 10kb
					elif down_len>=flank_length and up_len<flank_length:
						flanking_seq=\
						seq_dict[contig_acc][:int(align_end)+flank_length]
						new_uplen=up_len
						new_downlen=flank_length
					#4.) Both are shorter than 10kb
					else:
						flanking_seq=\
						seq_dict[contig_acc]
						new_uplen=up_len
						new_downlen=down_len

				try:	
					taxon=org_dict[genome_acc]

				except KeyError as e:

					taxon='unknown'
					if args.assemblies==True:
						for tax_line in tax_lines:
							if genome_acc in tax_line:
								print('Identical assembly found!')
								taxon=org_dict[tax_line.split('\t')[0]]
					
						print('Something went wrong when updating an entry - %s' % str(e))

					elif args.local!='False':
						for tax_line in tax_lines_cust:
							if genome_acc in tax_line:
								print('Identical assembly found!')
								taxon=org_dict[tax_line.split('\t')[0]]

				updated_lines.append(line.replace(line.split('\t')[10], flanking_seq+'\t'+str(new_uplen)+\
			'\t'+str(new_downlen)+'\t'+taxon))	

		#Write all flanking regions to .csv file
		print('Writing flanking regions to file...')
		with open(fa_file+'_flanking_regions.csv', 'w') as outfile:
			for line in updated_lines:
				outfile.write(line)
		outfile.close()

		#After finishing the current queue item, grab a new one until 'STOP' is encountered
		fa_file=queue.get()
		if fa_file=='STOP':
			return

def erase_previous():

	print('Cleaning up target directory...')
	#Go through all assembly directories and delete the results of previous analyses
	for root, dirs, files in os.walk(args.target_directory):
		for file in os.listdir(root):
			if file.startswith('all_') or '_summary.txt' in file or file.startswith('orfs_clustered.')\
			or file.endswith('_database.db') or file.\
			endswith('_reformatted.fna') or file.endswith('.log') or file.startswith('split_') \
			or file.startswith('flanking_regions'):
				os.remove(root+'/'+file)

		if any(element==os.path.basename(root) for element in ['genomes', 'taxonomy', 'kraken2']):
			shutil.rmtree(root)

	print('Assembly directories cleaned up!')

		
def run_prodigal(queue):
	
	#Grab items from queue, which contains paths to flanking file files
	fa_file=queue.get()
	if fa_file=='STOP':
		return

	while True:
		
		#Run prodigal
		if not os.path.exists(fa_file.replace('.fna', '_orfs.gff')):
			prodigal_call='prodigal -i %s -a %s -o %s -f gff -p meta' % \
			(fa_file, fa_file.replace('.fna', '_orfs.fna'), fa_file.replace('.fna', '_orfs.gff'))
			subprocess.call(prodigal_call, shell=True)

			#Rewrite headers in _orfs.fna
			orfs={}
			for line in open(fa_file.replace('.fna', '_orfs.fna'), 'r'):
				if line.startswith('>'):
					header=line.split(' ')[0]
					seq=''
				else:
					seq+=line
					orfs[header]=seq

			with open(fa_file.replace('.fna', '_orfs.fna'), 'w') as outfile:
				for key, value in orfs.items():
					outfile.write(key+'\n'+value)
			outfile.close()


		fa_file=queue.get()
		if fa_file=='STOP':
			return

def cluster_orfs(orf_file, target_dir):

	#Take a file containing all env gene seqs and cluster it. After clustering, create dict as follows:
	if not os.path.isfile(target_dir.rstrip('/')+'/orfs_clustered.fna'):
		cluster_command='cd-hit -T 0 -M 100000 -c 0.95 -s 0.95 -n 5 -i %s -o %s' % (orf_file, \
		target_dir.rstrip('/')+'/orfs_clustered.fna')

		subprocess.call(cluster_command, shell=True)

	#Create dictionary that contains each gene id as key and the representativ ecentroid as value
	clust_dict={}
	clust_lines=[line.rstrip('\n') for line in open(target_dir\
	.rstrip('/')+'/orfs_clustered.fna.clstr', 'r') if not line.startswith('>')]

	for line in clust_lines:
		if line.endswith('*'):
			centroid=line.split('>')[1].split('...')[0]

		clust_dict[line.split('>')[1].split('...')[0]]={}
		try:
			clust_dict[line.split('>')[1].split('...')[0]]['centroid']=centroid
		except:

			clust_dict[line.split('>')[1].split('...')[0]]['centroid']=line.split('>')[1].split('...')[0]

	if args.log==True:
		if os.path.exists(target_dir.rstrip('/')+'/orfs_clustered.fna.clstr') and \
			os.path.getsize(target_dir.rstrip('/')+'/orfs_clustered.fna.clstr')>0:
				log_lines.append('ORFs clustered...\n')
		else:
			log_lines.append("ORF cluster file wasn't found or is empty...FAILED\n")
		write_log()
	
	return clust_dict

def annotate_orfs(queue):

	fa_file=queue.get()
	if fa_file=='STOP':
		return
	
	while True:

		#Annotate orfs
		if args.uniprot_db=='False':
			uniprot_db_path=args.target_directory.rstrip('/')+'/uniprotKB.dmnd'
		else:
			if args.uniprot_db.endswith('.dmnd'):
				uniprot_db_path=args.uniprot_db
			else:
				print('Please provide --uniprot_db in .dmnd format!\nAborting...')
				sys.exit()

		if not os.path.exists(fa_file.replace('_orfs.fna', '_orfs_annotated.csv')):

			if args.processes*split_num>multiprocessing.cpu_count()*2:	
				threads=int((multiprocessing.cpu_count()*2)/args.processes)
			else:
				threads=args.processes

			diamond_call=f'diamond blastp -p {threads} -d {uniprot_db_path} -q {fa_file} -o {fa_file.replace("_orfs.fna", "_orfs_annotated.csv")} --id {args.uniprot_cutoff} --more-sensitive --max-target-seqs 1 --masking 0 --subject-cover 60 --log -f 6 qseqid sseqid stitle pident qstart qend qlen slen length qframe qtitle'
			subprocess.call(diamond_call, shell=True)
		fa_file=queue.get()
		if fa_file=='STOP':
			return

def create_env_dict(all_orfs_gff, all_orf_annos, clust_dict, all_IS_annos):

	#Create list containing all orfs from all flanking regions
	orfs=[line for file in all_orfs_gff for line in open(file, 'r') if not line.startswith('#')]

	#Create list containg all orf annotation information
	anno_orfs=[line for file in all_orf_annos for line in open(file, 'r')]

	#Create list containing all orf annotation information with IS
	anno_orfs_IS=[line for file in all_IS_annos for line in open(file, 'r')]

	#Create dictionary containing orf id and name	
	anno_dict={}
	for line in anno_orfs:
		anno_dict[line.split('\t')[0]]=' '.join(line.split('\t')[1].split('|')[2]\
		.split('_OS=')[0].split('_')[2:])
	
	#Now overwrite these annotations with 'high quality' IS annotations
	for line in anno_orfs_IS:
		anno_dict[line.split('\t')[0]]=line.split('\t')[1].split('|')[3].split(' [')[0]

	#add name key to clust_dict
	print('Assigning names to orfs...')
	for key, value in clust_dict.items():
		try:
			clust_dict[key]['name']=anno_dict[clust_dict[key]['centroid']]
		except:
			clust_dict[key]['name']='hypothetical protein'

	#Create dictionary from all orfs
	env_dict={}
	
	for line in orfs:
		gene_id=line.split('\t')[0]
		env_gene_id=line.split('\t')[8].split(';')[0].split('_')[1]
		new_name=clust_dict[str(gene_id)+'_'+str(env_gene_id)]['name']

		if not gene_id in env_dict:
			env_dict[gene_id]={}

		env_dict[gene_id][env_gene_id]={}
		env_dict[gene_id][env_gene_id]['env_name']=new_name
		env_dict[gene_id][env_gene_id]['env_start']=line.split('\t')[3]
		env_dict[gene_id][env_gene_id]['env_stop']=line.split('\t')[4]
		env_dict[gene_id][env_gene_id]['env_strand']=line.split('\t')[6]

	return env_dict

def annotate_integrons(queue, dir):

	flank_file=queue.get()
	if flank_file=='STOP':
		return

	while True:

		#Run integron finder
		if not os.path.exists(dir.rstrip('/')+'/all_integrons.txt'):
			print('annotating integrons...')
			search_int='integron_finder %s --cpu 30 --quiet --linear --outdir %s' % \
			(flank_file, dir)
			subprocess.call(search_int, shell=True)

			flank_file=queue.get()
			if flank_file=='STOP':
				return

		else:
			print('integron annotations already exist!')
			return

def integrons_to_db():

	if not args.update==True:
		int_folder=args.target_directory.rstrip('/')+'/integrons_tmp/'
	else:
		int_folder=args.target_directory.rstrip('/')+'/update_tmp/integrons_tmp/'

	#Get all integron and summary files from subfolders
	int_files=[root+'/'+file for root, dirs, files in os.walk(int_folder) \
	for file in files if file.endswith('.integrons')]

	sum_files=[root+'/'+file for root, dirs, files in os.walk(int_folder) \
	for file in files if file.endswith('.summary')]

	#summarize results and parse results into dictionary
	with open(int_folder+'all_integrons.txt', 'w') as outfile:
		for file in int_files:
			for line in open(file, 'r'):
				outfile.write(line.rstrip('\n')+'\n')
	outfile.close()

	with open(int_folder+'all_summaries.txt', 'w') as outfile:
		for file in sum_files:
			for line in open(file, 'r'):
				outfile.write(line.rstrip('\n')+'\n')
	outfile.close()


	int_lines=[line for line in open(int_folder+'all_integrons.txt', 'r') \
	if not line.startswith('ID') and not line.startswith('#')]

	#Parse results into dictionary
	integron_dict={}
	print('collecting integron data...')
	for line in int_lines:
		gene_acc=line.split('\t')[1]
		int_acc=line.split('\t')[0]
		el_id=line.split('\t')[2]

		if not gene_acc in integron_dict:
			integron_dict[gene_acc]={}
		if not int_acc in integron_dict[gene_acc]:
			integron_dict[gene_acc][int_acc]={}
		integron_dict[gene_acc][int_acc][el_id]={}
		integron_dict[gene_acc][int_acc][el_id]['start']=int(line.split('\t')[3])
		integron_dict[gene_acc][int_acc][el_id]['stop']=int(line.split('\t')[4])
		integron_dict[gene_acc][int_acc][el_id]['name']=line.split('\t')[8]
		integron_dict[gene_acc][int_acc][el_id]['strand']=line.split('\t')[5]

	#Find start and end position of integron, see if complete integron
	for key, value in integron_dict.items():
		for key2, value2 in value.items():
			start_list=[]
			stop_list=[]
			name_list=[]
			for key3, value3 in value2.items():
				start_list.append(value3['start'])
				stop_list.append(value3['stop'])
				name_list.append(value3['name'])
			integron_dict[key][key2]['int_start']=min(start_list)
			integron_dict[key][key2]['int_stop']=max(stop_list)
			if 'attC' in name_list and 'intI' in name_list:
				integron_dict[key][key2]['int_complete']='complete'
			else:
				integron_dict[key][key2]['int_complete']='incomplete'

	#Save results to database
	connection=sqlite3.connect(args.target_directory.rstrip('/')+'/genview_database.db')
	cursor=connection.cursor()

	if not args.update==True:
		#Create integron table
		cursor.execute('DROP TABLE IF EXISTS integrons')

		create_table="""
		CREATE TABLE integrons(
		arg_id INTEGER NOT NULL,
		int_id INTEGER NOT NULL PRIMARY KEY,
		int_start INTEGER NOT NULL,
		int_stop INTEGER NOT NULL,
		int_complete VARCHAR(100) NOT NULL,
		FOREIGN KEY (arg_id) REFERENCES args(id) ON DELETE CASCADE
		);
		"""
		cursor.execute(create_table)

		cursor.execute('DROP TABLE IF EXISTS int_elements')

		int_elements="""
		CREATE TABLE int_elements(
		id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
		int_id INTEGER NOT NULL,
		el_id INTEGER NOT NULL,
		el_start INTEGER NOT NULL,
		el_stop INTEGER NOT NULL,
		el_name VARCHAR(500),
		el_strand VARCHAR(4),
		FOREIGN KEY (int_id) REFERENCES integrons(int_id) ON DELETE CASCADE
		);
		""" 
		cursor.execute(int_elements)

		#Insert data into tables
		integron_id=0

	else:
		cursor.execute('SELECT MAX(int_id) FROM integrons')
		results=cursor.fetchall()
		integron_id=int(results[0][0])

	for key, value in integron_dict.items():
		for key2, value2 in value.items():
			integron_id+=1
			#into integron table
			cursor.execute('INSERT INTO integrons(arg_id, int_id, int_start, int_stop, \
			int_complete) VALUES(?,?,?,?,?)', (key, integron_id, value2['int_start'], \
			value2['int_stop'], value2['int_complete']))

			for key3, value3 in value2.items():
				if not key3.startswith('int_'):

					if value3['strand'].startswith('-'):
						strand='-'
					else:
						strand='+'

					#into integron element table
					cursor.execute('INSERT INTO int_elements(int_id, el_id, el_start, el_stop, \
					el_name, el_strand) VALUES(?,?,?,?,?,?)', (integron_id, key3, value3['start'], \
					value3['stop'], value3['name'], strand))
	connection.commit()

	print('Integron annotations saved to database!')

def multiprocess(target_func, processes, *items):
	
	#Create queue from files to process
	#NOTE: just one extra item, files have to be at first position
	queue=multiprocessing.Queue()

	for element in items[0]:
		queue.put(element)

	proc_list=[]
	p=processes

	for i in range(p):
		queue.put('STOP')

	started_procs=0
	for i in range(p):
		started_procs+=1

		if len(items)==1:
			proc_list.append(multiprocessing.Process(target=target_func, args=(queue,)))
		elif len(items)>1:
			proc_list.append(multiprocessing.Process(target=target_func, args=(queue, items[1])))
		elif len(items)>2:
			proc_list.append(multiprocessing.Process(target=target_func, args=(queue, items[1], items[2])))
		
		proc_list[-1].start()

	finished_procs=0
	for p in proc_list:
		finished_procs+=1
		p.join()
		print(f'{finished_procs} of {started_procs} processes finished...')


def to_sql_db(env_dict, all_anno, target_directory):

	#Create a dict containing all neccessary metainfo
	sum_dict={}
	unique_genomes={line.split('\t')[-2].split('__')[-1] for line in open(all_anno, 'r')}

	
	names=open(args.target_directory.rstrip("/")+"/taxonomy/names.dmp", 'r')

	#Create dictonary with organism names and tax_ids
	sci_names={line.split('\t|\t')[1]:line.split('\t|\t')[0] for line in names \
	if 'scientific name' in line.split('\t|\t')[3]}

	for genome in unique_genomes:

		if genome.endswith('.fna'):
			genome=genome.replace('.fna', '.1')
		sum_dict[genome]={}
		sum_dict[genome]['annotated genes']={}
	
	for line in open(all_anno, 'r'):

		genome_acc=line.split('\t')[-2].split('__')[-1]
		if genome_acc.endswith('.fna'):
			genome_acc=genome_acc.replace('.fna', '.1')
		prot_id=line.split('\t')[1].split('|')[1]\
		+'__'+line.split('\t')[-1].rstrip('\n')

		if args.kraken2=='False':
			sum_dict[genome_acc]['organism']=\
			line.split('\t')[-4]
		else:
			try:
				if len(line.split('\t')[-2].split(' '))>=2:
					sum_dict[genome_acc]['organism']=\
					' '.join(line.split('\t')[-2].split(' ')[1:3])
				else:
					sum_dict[genome_acc]['organism']=\
					line.split('\t')[-2].split(' ')[1]

			except Exception as e:
				print(f'{e}')
				sum_dict[genome_acc]['organism']=\
				'unclassified'
		try:
			if not 'plasmid' in sum_dict[genome_acc]['organism']:
				sum_dict[genome_acc]['taxon']=sci_names[sum_dict[genome_acc]['organism']]
			else:

				sum_dict[genome_acc]['taxon']=sci_names[' '.join(sum_dict[genome_acc]['organism']\
				.split(' ')[0:2])]
		except:
			sum_dict[genome_acc]['taxon']='unclassified'


		sum_dict[genome_acc]['annotated genes'][prot_id]={}
		sum_dict[genome_acc]['annotated genes'][prot_id]['name']=\
		line.split('\t')[1].split('|')[-1]
		sum_dict[genome_acc]['annotated genes'][prot_id]['card_id']=\
		prot_id.split('__')[0]
		sum_dict[genome_acc]['annotated genes'][prot_id]['perc_id']=\
		line.split('\t')[3]
		sum_dict[genome_acc]['annotated genes'][prot_id]['frame']=\
		line.split('\t')[-3]
		sum_dict[genome_acc]['annotated genes'][prot_id]['contig_acc']=\
		line.split('\t')[0]
		sum_dict[genome_acc]['annotated genes'][prot_id]['down_len']=\
		line.split('\t')[-5]
		sum_dict[genome_acc]['annotated genes'][prot_id]['up_len']=\
		line.split('\t')[-6]
		sum_dict[genome_acc]['annotated genes'][prot_id]['aln_len']=\
		int(line.split('\t')[8])*3
	


	#Create an sql database to save information on flanking genes, seqs etc in. Delete fasta files after
	#Create database file
	connection=sqlite3.connect(target_directory.rstrip('/')+'/genview_database.db')


	#Create a cursor to perform sql commands
	cursor=connection.cursor()

	if not args.update==True:
		#Create Tables - drop pre-existing tables
		cursor.execute('DROP TABLE IF EXISTS X')

		#Create lineage tables		
		#TODO: disable this, this comes now at the beginning of main function

		if args.assemblies==True:
			summary_files=[args.target_directory.rstrip('/')+'/'+file for \
			file in os.listdir(args.target_directory) \
			if file.startswith('assembly_summary')]

			new_lineages=add_taxonomy_lineages(summary_files)

		if args.plasmids==True:

			summary_files=[args.target_directory.rstrip('/')+'/'+file for \
			file in os.listdir(args.target_directory) \
			if file.startswith('plasmid_summary')]

			new_lineages=add_taxonomy_lineages(summary_files)

		if args.local!='False':

			summary_files=[args.target_directory.rstrip('/')+'/'+file for \
			file in os.listdir(args.target_directory) \
			if file.startswith('local_summary')]

			new_lineages=add_taxonomy_lineages(summary_files)


		create_genomes_table="""
		CREATE TABLE IF NOT EXISTS genomes(
			id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
			assembly VARCHAR(500) UNIQUE NOT NULL,
			organism VARCHAR(500) NOT NULL,
			taxon_id INTEGER NOT NULL,
			FOREIGN KEY (taxon_id) REFERENCES lineages(taxon)
			);	
		"""

		cursor.execute(create_genomes_table)

		create_arg_table="""
		CREATE TABLE IF NOT EXISTS args(
			id INTEGER PRIMARY KEY NOT NULL UNIQUE,
			arg_name VARCHAR(100),
			perc_id FLOAT NOT NULL,
			genome_id INTEGER NOT NULL,
			contig_acc VARCHAR(100) NOT NULL,
			arg_len INTEGER NOT NULL,
			uplen INTEGER NOT NULL,
			downlen INTEGER NOT NULL,
			arg_start INTEGER NOT NULL,
			arg_end INTEGER NOT NULL,
			frame INTEGER NOT NULL,
			FOREIGN KEY (genome_id) REFERENCES genomes(id) ON DELETE CASCADE
			);
		"""

		cursor.execute(create_arg_table)

		create_env_table="""
		CREATE TABLE IF NOT EXISTS envs(
			id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
			env_name VARCHAR(100) NOT NULL,
			arg_id INTEGER NOT NULL,
			env_start INTEGER NOT NULL,
			env_end INTEGER NOT NULL,
			env_strand VARCHAR(5) NOT NULL,
			FOREIGN KEY (arg_id) REFERENCES args(id) ON DELETE CASCADE
			);
		"""

		cursor.execute(create_env_table)

	#insert data
	if not args.update==True:
		genome_id=0
	else:
		cursor.execute('SELECT MAX(id) FROM genomes')
		results=cursor.fetchall()
		genome_id=int(results[0][0])
		pre_update_count=int(results[0][0])

	#Insert genomes
	for key, value in sum_dict.items():
		genome_id+=1
		try:
			if args.kraken2!='False':
				cursor.execute('INSERT INTO genomes(assembly, organism, taxon_id) VALUES(?,?,?);',\
				 (key.split(' ')[0], sum_dict[key]['organism'], sum_dict[key]['taxon']))	
				print(key, sum_dict[key]['organism'], sum_dict[key]['taxon'])
			else:

				cursor.execute('INSERT INTO genomes(assembly, organism, taxon_id) VALUES(?,?,?);',\
				 (key.split(' ')[0], sum_dict[key]['organism'], sum_dict[key]['taxon']))	
				

		except Exception as e:
			print(f'{key} not inserted in database - {e}')
			continue
			

		#Insert annotated genes
		for key2, value2 in value['annotated genes'].items():


			#change frame to plus, as all sequences where the ARG is on - strand 
			#are reverse complemented in create_db() - so all args should be encoded on plus strand!
			#This was changed in v9.4
			arg_start=value2['up_len']
			arg_end=int(value2['up_len'])+int(value2['aln_len'])

			if value2['frame'].startswith('-'):

				arg_start=value2['down_len']
				arg_end=int(value2['down_len'])+int(value2['aln_len'])
				strand='+'
			else:
				strand='+'

			cursor.execute('INSERT INTO args(arg_name, genome_id, contig_acc, perc_id, \
			arg_len, frame, uplen, downlen, arg_start, arg_end, id) VALUES(?,?,?,?,?,?,?,?,?,?,?);',\
			(value2['name'], genome_id, value2['contig_acc'], value2['perc_id'], value2['aln_len'],\
			strand, value2['up_len'], value2['down_len'], arg_start, arg_end,\
			key2.split('__')[-1]))

			try:
				#Insert flanking regions
				for key3, value3 in env_dict[key2.split('__')[-1]].items():
					cursor.execute('INSERT INTO envs(env_name, arg_id, env_start, env_end, env_strand) \
							VALUES(?,?,?,?,?);', (value3['env_name'], key2.split('__')[-1], \
							value3['env_start'], value3['env_stop'], value3['env_strand']))
			except KeyError as e:
				print('No genetic environment found for %s - %s' % (key2, e))
	connection.commit()
	print('Genes and flanking regions saved to database!')

	
	cursor.execute('SELECT MAX(id) FROM genomes')
	results=cursor.fetchall()
	post_update_count=int(results[0][0])

	if args.log==True:
		if args.update==True:
			if pre_update_count<post_update_count:
				log_lines.append(f'{post_update_count-pre_update_count} genomes added to database...'+'\n')
			else:
				log_lines.append('No new genomes in database after update...FAILED(?)\n')
		else:
			if post_update_count>0:
				log_lines.append(f'{post_update_count} genomes added to database...'+'\n')
			else:
				log_lines.append('No genomes added to database...FAILED\n')
		write_log()	

def check_lr_headers():
	
	#Check weather headers of long reads are in right format
	lr_files=[file for file in os.listdir(args.local) if (file.endswith('.fna') or file.endswith('fasta')\
	 or file.endswith('.fa'))]

	lr_headers=[line for file in lr_files for line in open(f'{os.path.abspath(args.local)}/{file}', 'r') if line.startswith('>')]
	if len(lr_headers)==0:
		print('No headers detected for your long read files, make sure your header lines start with ">"')
		sys.exit()
	else:
		#Count numbers of tabs and spaces in header
		tabs=lr_headers[1].count('\t')
		spaces=lr_headers[1].count(' ')

		if tabs>0:
			print('\\t is not allowed in read headers. Please make sure your header format is ">readid\\n" or ">readid organism\\n", t.ex ">SRR1234.1\\n" or ">SRR1234.1 Staphylococcus aureus\\n" ')
			sys.exit()
	
		if spaces>2:
			print('Multiple spaces detected in local fasta headers. Please make sure your header format is ">readid\\n" or ">readid organism\\n", t.ex ">SRR1234.1\\n" or ">SRR1234.1 Staphylococcus aureus\\n"')
			sys.exit()


def main():

	global args
	args=parse_arguments()

	if not os.path.exists(os.path.abspath(args.database)):
		print('\nPath specified with -db does not exist, please provide a valid path!\n')
		sys.exit()

	if not (args.database.endswith('.fna') or args.database.endswith('.fa') or args.database.endswith('.fasta') or args.\
	database.endswith('.faa')):
		print('\nIs your -db file in fasta format? Please provide a file ending with .fna, .fa, .faa or .fasta.\n')
		sys.exit()

	if not os.path.exists(args.target_directory):
		print('Target directory does not exist, creating...')
		os.mkdir(args.target_directory)

	if args.clean==True and args.update==True:
		print('Output from previous run is required for update, --clean and --update cannot be specified at the same time!')
		sys.exit()

	if args.clean==True and not args.update==True:
		erase_previous()

	if args.log==True:
		global log_lines
		log_lines=[]		
		#Write command to log file
		log_lines.append('_'*20+'\n'*2)

		command_string='genview-makedb '
		for key, value in vars(args).items():
			if not value=='False':
				if not value==True:
					if not key=='taxa':						
						command_string+=f'--{key} {value} '
					else:
						command_string+=f"--{key} "+" ".join(f"'{taxon}'" for taxon in value)+' '
				else:
					command_string+=f'--{key} '

		log_lines.append(command_string+'\n'+'_'*20+'\n')
		write_log()

	if args.long_reads!='False':
		print('\nLong read function is deactivated atm. For corrected long reads, use --local only.')
		sys.exit()
	
	if args.local!='False':
		if not os.path.isdir(args.local):
			print('--local should provide the path to directory containing local genome files!')
			sys.exit()

		if os.path.abspath(args.local)==os.path.abspath(args.target_directory):
			print('\nTarget directory and local genome directory are the same - Please create a separate directory for your local sequences and specify it using --local')
			sys.exit()

		#Check that provided files are fasta files
		loc_files=[file for file in os.listdir(args.local)]

		non_fasta=False
		for file in loc_files:
			if file.endswith('.fna') or file.endswith('.fasta') or file.endswith('.fa') or \
			file.endswith('.faa') or file.endswith('.ffn'):
				pass
			else:
				non_fasta=True

		if non_fasta==True:
			print('\nFile endings indicate that there are non-fasta files in the path specified under --local.\nPlease provide fasta files only in this directory, ending on ".fna", ".fa" or ".fasta"\n') 		
			sys.exit()

	if not args.update==True:
		#Check if target directory contains files from a previous run
		targetdir_files=[file for file in os.listdir(args.target_directory)]
		if ('all_annos.fna_tmp' or 'all_flanks.csv_tmp' or 'all_orfs.fna' or 'assembly_summary.txt' or \
		'genview_database.db' or 'genview.log' or 'orfs_clustered.fna' or 'plasmid_summary.txt.0' or 'all_assemblies_0.fna') in \
		targetdir_files:
			print('\nFiles from previous genview run detected in target directory!\nIf you want to update an existing database, use --update.\nIf you want to start a new run on a previously used target directory, use --clean!\n')
			sys.exit()

	if args.log==True and args.update==True:
		
		#get log file from previous run and add content to log lines
		if os.path.exists(f'{os.path.abspath(args.target_directory)}/genview.log'):
			old_log=[line for line in open(f'{os.path.abspath(args.target_directory)}/genview.log')]
			log_lines.extend(old_log)
			old_log.append('\n'+'_'*20+'\n')

		#now delete old log files
		if os.path.exists(f'{os.path.abspath(args.target_directory)}/genview.log'):
			os.remove(f'{os.path.abspath(args.target_directory)}/genview.log')

		if os.path.exists(f'{os.path.abspath(args.target_directory)}/diamond.log'):
			os.remove(f'{os.path.abspath(args.target_directory)}/diamond.log')

	if args.taxa=='False' and args.accessions=='False' and args.local=='False':
		print('No genomes specified for analysis, please specify either "--accessions" or "" --taxa')
		sys.exit()


	if args.kraken2!='False':
		db_files=[file for file in os.listdir(args.kraken2) if file.endswith('.k2d')]
		if len(db_files)<1:
			print('No kraken2 database files detected in specified path, if you want to use kraken2 you need to build a kraken2 database first!')
			sys.exit()
	
	if args.long_reads==True and args.local=='False':
		print('If --long_reads is specified. --local also needs to be specified!')
		sys.exit()

	if args.long_reads==True:
		check_lr_headers()

	if (args.assemblies==True or args.plasmids==True) and (args.long_reads==True or args.local!='False'):
		print('--local and --assemblies/plasmids cannot be specified at the same time!\n\
		Create a database with either local data or assemblies/plasmids first and then update it.')
		sys.exit()

	if args.accessions!='False' and not (args.assemblies==True or args.plasmids==True):
		print('\nIf accessions are provided using --accession, please specify whether they are assembly accessions, plasmid accessions or both using --assemblies and/or --plasmids.')
		sys.exit()

	if args.taxa!='False' and args.accessions!='False':
		print('\n--taxa cannot be specified at the same time as --accessions, please choose only one option.\n\
		Use the update option to combine local data and data from NCBI.')
		sys.exit()

	if args.uniprot_db!='False' and not args.uniprot_db.endswith('.dmnd'):
		print('\n--uniprot_db must specify path to diamond database ending on .dmnd!\n')
		sys.exit()
	
	if args.update==True:
		print('Updating genview database...')
		update()		
		sys.exit()

	download_uniprot()

	#Disable for debugging
	args.db=reformat()

	#Create output directory if not exists
	if not os.path.exists(args.target_directory):
		os.mkdir(args.target_directory)

	
	#Create directory to store sequence data
	if args.update==True:
		#Set target directory to temporary update directory
		if not os.path.exists(args.target_directory.rstrip('/')+'/'+'update_tmp/genomes'):
			os.mkdir(args.target_directory.rstrip('/')+'/'+'update_tmp/genomes')
	else:

		if not os.path.exists(args.target_directory.rstrip('/')+'/genomes'):
			os.mkdir(args.target_directory.rstrip('/')+'/genomes')

	#Download NCBI taxonomy files for lineage assignment
	if not os.path.exists(args.target_directory.rstrip('/')+'/taxonomy'):
		os.mkdir(args.target_directory.rstrip('/')+'/taxonomy')
	
	if not os.path.exists(args.target_directory.rstrip('/')+'/taxonomy/taxdump.tar.gz'):

		print('Downloading taxonomy data...\n')

		download=f'wget -P {args.target_directory.rstrip("/")+"/taxonomy/"} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
		subprocess.call(download, shell=True)

		unzip=f'tar -C {args.target_directory.rstrip("/")+"/taxonomy"} -xvzf {args.target_directory.rstrip("/")+"/taxonomy/taxdump.tar.gz"}'
		subprocess.call(unzip, shell=True)


	#Download assembly summary files
	if args.assemblies==True:
		print('Fetching assemblies and adding lineage information...')
		if not os.path.exists(args.target_directory.rstrip('/')+'/assembly_summary.txt'):
			print('Downloading assembly summary...')
			download_file='wget -P %s ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt' \
			% (args.target_directory)
			subprocess.call(download_file, shell=True)
	
		#Create assembly file only containing the data on specified taxa
		if args.taxa!='False':
			if not args.taxa[0]=='all':
				tax_asms=[]
				taxon_lines=[line for line in open(args.target_directory.rstrip('/')+'/assembly_summary.txt', 'r')\
				if not line.startswith('#')]

				for line in taxon_lines:
					if any(taxon.lower() in line.lower() for taxon in args.taxa):
						tax_asms.append(line)

			with open(args.target_directory.rstrip('/')+'/assembly_summary.txt', 'w') as f:
				for line in tax_asms:
					f.write(line)
			f.close()
		
		if args.log==True:
			if os.path.exists(os.path.abspath(args.target_directory)+'/assembly_summary.txt') and \
			os.path.getsize(os.path.abspath(args.target_directory)+'/assembly_summary.txt')>0:
				log_lines.append('Assembly summary downloaded and rewritten...\n')
			else:
				log_lines.append('Assembly summary not found or empty...FAILED\n')
			write_log()
	
		if not args.accessions=='False':
			accessions=[line.rstrip('\n') for line in open(args.accessions, 'r')]
			asm_sum_lines=[line for line in open(args.target_directory.rstrip('/')\
			+'/assembly_summary.txt', 'r') if any(acc in line for acc in accessions)]

			with open(args.target_directory.rstrip('/')+'/assembly_summary.txt', 'w') as f:
				for line in asm_sum_lines:
					f.write(line)
			f.close()

		#Parse summary file for genomes to download
		asm_summary=[line for line in open(args.target_directory\
		.rstrip('/')+'/assembly_summary.txt', 'r') if not line.startswith('#')] 

		#Create dictionary with assembly accession as key, lineages as value
		summary_files=[args.target_directory.rstrip('/')+'/'+file for \
		file in os.listdir(args.target_directory) \
		if file.startswith('assembly_summary')]

		new_lineages=add_taxonomy_lineages(summary_files)

		asm_dict={}
		excepts=0
		for line in asm_summary:
			
			asm_dict[line.split('\t')[0]]={}
			asm_dict[line.split('\t')[0]]['tax_id']=line.split('\t')[5]
			asm_dict[line.split('\t')[0]]['url']=line.split('\t')[19]
			try:
				asm_dict[line.split('\t')[0]]['sci_lineage']=new_lineages[line.split('\t')[5]]['sci_lineage']
				asm_dict[line.split('\t')[0]]['lineage']=new_lineages[line.split('\t')[5]]['lineage']
			except:
				excepts+=1
				asm_dict[line.split('\t')[0]]['sci_lineage']='unassigned'
				asm_dict[line.split('\t')[0]]['lineage']='unassigned'

		print(f'lineages assigned, {excepts} of {len(asm_dict)} could not be assigned')
		if excepts>100:
			print('Maybe it\'s time to update the database?')


		#Now multiprocess the download of these new genomes into temporary folder
		#disable for now to make sure not all genomes have to be downloaded again
		#TODO add exact line lineage in summary file
		if not args.accessions=='False':
			genome_accessions=[line.rstrip('\n').lower() for line in open(args.accessions, 'r')]
			genome_urls=[asm_dict[key]['url'] for key, value in asm_dict.items() if key.lower() in genome_accessions]

		else:
			if not args.taxa[0]=='all':
				genome_urls=[asm_dict[key]['url'] for key, value in asm_dict.items() if any(taxon.lower() in \
				' '.join(asm_dict[key]['sci_lineage']).lower() for taxon in args.taxa)]
			else:
				
				genome_urls=[asm_dict[key]['url'] for key, value in asm_dict.items()]
		
		if args.log==True:
			if len(genome_urls)>=1:
				log_lines.append(f'{len(genome_urls)} target genomes identified...'+'\n')
			else:
				log_lines.append(f'No target genomes identified...FAILED'+'\n')
		
		if len(genome_urls)>=1:

			if not os.path.exists(f'{os.path.abspath(args.target_directory)}/genomes'):
				os.mkdir(f'{os.path.abspath(args.target_directory)}/genomes')

			multiprocess(download_new, args.processes, genome_urls)
		else:
			print('No genomes found that fit the provided input, exiting...')

	if args.plasmids==True:
		print('Fetching plasmids...')
		download_plasmids()

		if 'genome_urls' in locals():
			genome_urls.extend([line for line in open(f'{os.path.abspath(args.target_directory)}/plasmid_summary.txt.0', "r")])
		else:
			genome_urls=[line for line in open(f'{os.path.abspath(args.target_directory)}/plasmid_summary.txt.0', 'r')]

	if args.local!='False':
		genome_urls=os.listdir(args.local)

	#If no local summary file exist, create one with the genome IDs
	if args.local!='False' and not os.path.isfile(args.target_directory.rstrip('/') + '/local_summary.txt.0'):
		create_simple_summary_file()

	#make split_num available in other functions	
	global split_num

	if 1<=len(genome_urls)<1000:	
		split_num=2
	elif 1000<len(genome_urls)<10000:	
		split_num=5
	elif 10000<len(genome_urls)<100000:	
		split_num=20
	else:
		split_num=200
	#Determine whether split files are already present
	split_files=[file for file in os.listdir(args.target_directory) if file\
	.startswith('all_assemblies_') and file.endswith('.csv') and len(file.split('.'))==2]

	if not os.path.isfile(args.target_directory.rstrip('/')+'/all_assemblies.fna') \
	or len(split_files)==0:
		concatenate_and_split()	
	
	print('collecting fasta files for annotation...')
	#Create list of files containing assemblies in fasta format
	fa_files=[args.target_directory.rstrip('/')+'/'+fa_file for fa_file \
	in os.listdir(args.target_directory)\
	if fa_file.startswith('all_assemblies_')\
	and fa_file.endswith('.fna') and not fa_file.startswith('flanking') \
	and not '_orfs' in fa_file]
	
	#Delete 'all_assemblies.fna' at this stage to free up storage space
	all_assembliesfna=[args.target_directory.rstrip('/')+'/'+fa_file for fa_file \
	in os.listdir(args.target_directory) if fa_file==('all_assemblies.fna')]
	os.remove(all_assembliesfna[0])

	#Define number of processes and empty list for processes
	multiprocess(annotate, args.processes, fa_files)

	if args.log==True:

		#move diamond.log from current wd to target directory
		if not os.path.exists(f'{args.target_directory}/diamond.log'):
			mv_command=f'mv {os.getcwd()}/diamond.log {os.path.abspath(args.target_directory)}'
		if os.path.exists(f'{os.path.abspath(args.target_directory)}/diamond.log') and not f'{os.path.abspath(args.target_directory)}'==os.getcwd():
			mv_command=f'cat {os.getcwd()}/diamond.log {os.path.abspath(args.target_directory)}/diamond.log > {os.path.abspath(args.target_directory)}/diamond.log'
		subprocess.call(mv_command, shell=True)

		dmnd_log=[line for line in open(f'{os.path.abspath(args.target_directory)}/diamond.log', 'r') \
		if 'queries aligned' in line]
		if len(dmnd_log)==len(fa_files):
			log_lines.append(f'diamond blastx run successfull...'+'\n')
		else:
			log_lines.append(f'diamond blastx run unsuccessfull...FAILED'+'\n')
		write_log()
	#Do the same again for the 'create_database' function

	if args.log==True:
		#Check presence and size of annotation file
		anno_log=[file for file in os.listdir(args.target_directory.rstrip('/'))\
			if file.endswith('_annotated.csv')]
		if len(anno_log)>0:
			log_lines.append('Annotation files generated...\n')
		else:
			log_lines.append('No annotation files could be found...FAILED\n')
		write_log()

	#Create list with files that have no corresponding flanking region file yet
	no_flank_files=[element for element in fa_files if not os.path.exists(element+'_flanking_regions.csv')]

	if len(no_flank_files)>0:
		multiprocess(create_db, args.processes, no_flank_files)

	print('collecting flanking regions...')
	flank_files=[args.target_directory.rstrip('/')+'/'+fa_file for fa_file in os.listdir(args.target_directory)\
	if fa_file.startswith('all_assemblies_')\
	and fa_file.endswith('.fna_flanking_regions.csv')]

	if args.log==True:
		if len(flank_files)>0:
			log_lines.append('flanking region files generated...\n')
		else:
			log_lines.append('No flanking region files found...FAILED\n')
		write_log()	

	#summarize all flanking regions into one temporary file
	lines=[line for file in flank_files for line in open(file, 'r')]

	#If no flanking regions were found, exit
	if len(lines)==0:
		print('No database matches found, exiting...')
		sys.exit()

	print('Writing temporary files...')
	with open(args.target_directory.rstrip('/')+'/all_flanks.csv_tmp', 'w') as outfile:
		i=0
		for line in sorted(lines):
			#append ID for identification of env genes
			asm_name=line.split('\t')[-4]
			if asm_name.endswith('.fna'):
				asm_name=asm_name.replace('.fna', '.1')

			i+=1

			#Write in csv format to outfile
			outfile.write(str(i)+'\t'+line.split('\t')[-6]+'\n')
	outfile.close()

	#after assigning id, split the file into specified number of smaller fasta files again
	outfiles=args.processes
	csv_lines=[line for line in open(args.target_directory.rstrip('/')+'/all_flanks.csv_tmp', 'r')]

	line_count=0
	file_count=1
	outfile=open(args.target_directory.rstrip('/')+'/'+\
	'flanking_regions_'+str(file_count)+'.fna', 'w')
	for line in csv_lines:
		line_count+=1
		if line_count<=round((len(csv_lines)/args.processes)+0.5, 0):
			outfile.write('>'+line.split('\t')[0]+'\n'+line.split('\t')[-1])
		else:
			file_count+=1
			line_count=0
			outfile.close()
			outfile=open(args.target_directory.rstrip('/')+'/'+\
			'flanking_regions_'+str(file_count)+'.fna', 'w')
			outfile.write('>'+line.split('\t')[0]+'\n'+line.split('\t')[-1])
	outfile.close()

	#Create list with flanking regions file paths
	flanking_fa=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
	if file.startswith('flanking_regions') and file.endswith('.fna') and not '_orfs' in file]

	#Create a multiprocessing queue from the flanking region files
	multiprocess(run_prodigal, args.processes, flanking_fa)

	if args.log==True:
		orf_fna_log=[file for file in os.listdir(args.target_directory.rstrip('/'))\
		if file.endswith('_orfs.fna')]

		orf_gff_log=[file for file in os.listdir(args.target_directory.rstrip('/'))\
		if file.endswith('_orfs.gff')]
	
		if len(orf_fna_log)>0 and len(orf_gff_log)>0:
			log_lines.append('ORFs predicted...\n')
		else:
			log_lines.append('No ORF files present, something might have gone wrong when running prodigal...\n')
		write_log()	

	#concatenate all orf files into one for clustering
	concat_command='cat %s/flanking_regions_*_orfs.fna > %s/all_orfs.fna' % \
	(args.target_directory.rstrip('/'), args.target_directory.rstrip('/'))
	subprocess.call(concat_command, shell=True)

	#Call clustering function which returns cluster dict
	clust_dict=cluster_orfs(args.target_directory.rstrip('/')+'/all_orfs.fna', \
	args.target_directory)

	#Split resulting centroid file into several for annotation
	print('splitting orf file for annotation...')
	centr_dict={}
	for line in open(args.target_directory.rstrip('/')+'/orfs_clustered.fna'):
		if line.startswith('>'):
			header=line
			seq=''
		else:
			seq+=line
			centr_dict[header]=seq

	outfiles=split_num

	key_count=0
	file_count=1
	outfile=open(args.target_directory.rstrip('/')+'/'+\
	'split_'+str(file_count)+'_orfs.fna', 'w')
	for key, value in centr_dict.items():
		key_count+=1
		#Make number of splits dependent on split number, not process number
		if key_count<=round((len(centr_dict)/split_num)+0.5, 0):
			outfile.write(key+value)
		else:
			outfile.close()
			file_count+=1
			key_count=0
			outfile=open(args.target_directory.rstrip('/')+'/'+\
			'split_'+str(file_count)+'_orfs.fna', 'w')
			outfile.write(key+value)

	#memory needs to be written to file here
	outfile.close()

	if args.log==True:
		orf_split_log=[file for file in os.listdir(args.target_directory.rstrip('/'))\
		if 'split_' in file and file.endswith('_orfs.fna')]
		if len(orf_split_log)>0:
			log_lines.append('ORF files split for annotation...\n')
		else:
			log_lines.append('ORF files not split for annotation...FAILED\n')	

	#Do the exact same thing for 'annotate_orfs' function
	orf_files=[args.target_directory.rstrip('/')+'/'+file for file in \
	os.listdir(args.target_directory) if file.endswith('_orfs.fna') \
	and file.startswith('split_')]

	multiprocess(annotate_orfs, args.processes, orf_files)

	if args.log==True:
		#move diamond.log from current wd to target directory
		if not os.path.exists(f'{os.path.abspath(args.target_directory)}/diamond.log'):
			mv_command=f'mv {os.getcwd()}/diamond.log {os.path.abspath(args.target_directory)}'
		
		if os.path.exists(f'{os.path.abspath(args.target_directory)}/diamond.log') and not f'{os.path.abspath(args.target_directory)}'==os.getcwd():
			mv_command=f'cat {os.getcwd()}/diamond.log {os.path.abspath(args.target_directory)}/diamond.log > {os.path.abspath(args.target_directory)}/dmnd.log'
		subprocess.call(mv_command, shell=True)

		if os.path.exists(f'{os.getcwd()}/diamond.log'):
			os.remove(f'{os.getcwd()}/diamond.log')

		#now rename diamond log
		mv2=f'mv {os.path.abspath(args.target_directory)}/dmnd.log {os.path.abspath(args.target_directory)}/diamond.log'
		subprocess.call(mv2, shell=True)

		dmnd_log=[line for line in open(f'{os.path.abspath(args.target_directory)}/diamond.log', 'r') \
		if 'queries aligned' in line]
		if len(dmnd_log)/2==len(fa_files):
			log_lines.append(f'diamond blastp run successfull...'+'\n')
		else:
			log_lines.append(f'diamond blastp run unsuccessfull...FAILED'+'\n')
		write_log()

	#Create a temporary summary file containing the line id
	with open(args.target_directory.rstrip('/')+'/all_annos.fna_tmp', 'w') as outfile:
		i=0
		for line in sorted(lines):
			i+=1
			outfile.write(line.rstrip('\n')+'\t'+str(i)+'\n')
	outfile.close()

	print('summarizing files...')
	all_orfs_gff=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
	if file.endswith('_orfs.gff')]	

	all_orf_annos=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
	if file.endswith('_orfs_annotated.csv')]

	all_IS_annos=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory) \
	if file.endswith('_orfs_ISannotated.csv')]

	print('combining genome features...')
	#Now create a function that creates the env dict!
	env_dict=create_env_dict(all_orfs_gff, all_orf_annos, clust_dict, all_IS_annos)

	to_sql_db(env_dict, args.target_directory.rstrip('/')+'/all_annos.fna_tmp', args.target_directory)

	#Create a shared dict for multiple processes to write to

	#Create list with flanking region files
	flanks=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(\
	args.target_directory) if file.startswith('flanking_regions_') and \
	file.endswith('.fna') and not '_orfs' in file]

	if args.integron_finder=='True':
		#Create temporary directory
		if not os.path.exists(args.target_directory.rstrip('/')+'/integrons_tmp'):
			os.mkdir(args.target_directory.rstrip('/')+'/integrons_tmp')

		dir=args.target_directory.rstrip('/')+'/integrons_tmp'

		multiprocess(annotate_integrons, args.processes, flanks, dir)	
		integrons_to_db()

	#Create arg/tnp association table
	transposon_table()

	if not args.save_tmps==True:
		remove_tmps()

def remove_tmps():

	#Remove temporary files to minimize storage use
	#Create lists with files to delete
	asms=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory)\
	if file.startswith('all_assemblies_')]

	splits=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory)\
	if file.startswith('split_')]

	flanks=[args.target_directory.rstrip('/')+'/'+file for file in os.listdir(args.target_directory)\
	if file.startswith('flanking_')]

	tmps=[asms, splits, flanks]

	for tmp_list in tmps:
		for element in tmp_list:
			os.remove(element)
	
	#Also remove downloaded genomes
	if os.path.exists(f'{os.path.abspath(args.target_directory)}/genomes'):
		shutil.rmtree(f'{os.path.abspath(args.target_directory)}/genomes')

	print('temporary files removed!')

def transposon_table():

	#Create table that contains args having a transposon 3000kbp upstream/downstream
	#of them and associated transposon annotation
	connection=sqlite3.connect(args.target_directory.rstrip('/')+'/genview_database.db')
	cursor=connection.cursor()

	find_tnp_args="""
	SELECT args.id, arg_name, envs.id, env_name FROM args
	INNER JOIN envs ON args.id=envs.arg_id
	WHERE (env_name LIKE 'IS%' OR env_name LIKE '%ranspos%'
	OR env_name LIKE 'insertion sequen%'
	OR env_name LIKE 'ISCR%'
	OR env_name LIKE 'tnp%'
	OR env_name LIKE 'Tnp%')
	AND (env_name NOT LIKE 'Isocho%')
	AND (env_name NOT LIKE 'Istb%')
	AND (env_name NOT LIKE 'Isoleu')
	AND ((env_start>=uplen-4000 AND env_end<=uplen)
	OR (env_start>=uplen+arg_len AND env_end<=uplen+arg_len+4000))
	"""

	#Execute search statement
	cursor.execute(find_tnp_args)
	results=cursor.fetchall()

	cursor.execute('DROP TABLE IF EXISTS arg_tnps')

	#Creae new table from results
	tnp_arg_table="""
	CREATE TABLE IF NOT EXISTS arg_tnps(
	id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
	arg_id INTEGER NOT NULL,
	arg_name VARCHAR(100) NOT NULL,
	env_id INTEGER NOT NULL,
	tnp_name VARCHAR(200) NOT NULL,
	FOREIGN KEY (arg_id) REFERENCES args(id)
	FOREIGN KEY (env_id) REFERENCES envs(id)
	);
	"""

	cursor.execute(tnp_arg_table)

	#Insert data into table
	for element in results:
		cursor.execute('INSERT INTO arg_tnps(arg_id, arg_name, env_id, tnp_name) VALUES(?,?,?,?)', \
				(element[0], element[1], element[2], element[3]))

	connection.commit()

def create_simple_summary_file():
	
	#Write file names instead of assembly headers to get number of genomes/metagenomes
	loc_files=[file for file in os.listdir(args.local)]	

	#Determine number of local summary files present and write file
	local_sums=[file for file in os.listdir(args.target_directory) if file.startswith('local_summary')]

	with open(f'{args.target_directory}/local_summary.txt.{str(len(local_sums))}', 'w') as outfile:
		for file in loc_files:
			outfile.write(file.split('.')[0]+'\n')	
	outfile.close()


def kraken2_classify(infile):

	#Enable updating
	if args.update==True:
		target_dir=f'{os.path.abspath(args.target_directory)}/update_tmp'
	else:
		target_dir=args.target_directory

	#Create kraken2 output directory
	if not os.path.exists(target_dir.rstrip('/')+'/kraken2'):
		os.mkdir(target_dir.rstrip('/')+'/kraken2')

	#Classify reads with kraken2
	if not os.path.exists(f'{os.path.abspath(target_dir)}/kraken2/{os.path.basename(infile).replace(".fna", ".kraken2.out")}'):
		kraken2=f'kraken2 --db {args.kraken2} --confidence 0.05 --use-names --output {os.path.abspath(target_dir)}/kraken2/{os.path.basename(infile).replace(".fna", ".kraken2.out")} --report {os.path.abspath(target_dir)}/kraken2/{os.path.basename(infile).replace(".fna", ".kraken2.rep")} --threads {args.processes} {infile}'
		subprocess.call(kraken2, shell=True)

	#Read in read file
	reads={}
	for line in open(infile, 'r'):
		if line.startswith('>'):
			header=line.lstrip('>').split(' ')[0].rstrip('\n')
			seq=''
		else:
			seq+=line
			reads[header]=seq

	print(f'number of seqs:{len(reads)}')
	#Write species/genus assignment to file
	#parse output file from kraken	
	class_dict={line.split("\t")[1]+' '+line.split("\t")[2].split("(taxi")[0].rstrip(' ')+'\n':reads\
		[line.split('\t')[1]] for line in \
		open(f'{os.path.abspath(target_dir)}/kraken2/{os.path.basename(infile).replace(".fna", ".kraken2.out")}')}

	#Write species assignment and read number to file
	with open(f'{infile}', 'w') as f:
		for key, value in class_dict.items():
			f.write('>'+key+value)

def write_log():

	with open(f'{os.path.abspath(args.target_directory)}/genview.log', 'w') as outfile:
		for line in log_lines:
			outfile.write(line)
	outfile.close()


if __name__=='__main__':

	main()
