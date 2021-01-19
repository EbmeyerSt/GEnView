#!/usr/local/env python3.7

import sys, os, subprocess, argparse, ete3, sqlite3
from argparse import RawTextHelpFormatter
from ete3 import SeqMotifFace, TreeStyle, add_face_to_node, Tree, NodeStyle, TextFace, COLOR_SCHEMES, SVG_COLORS, random_color

def parse_arguments():
	descr='%r\n\nVisualize annotate genes and genetic environments\n%r' % ('_'*80, '_'*80)
	parser=argparse.ArgumentParser(description=descr.replace("'", ''), formatter_class=RawTextHelpFormatter)
	parser.add_argument('-db', help='sqlite3 db containing annotations', required=True)
	parser.add_argument('-o', help='target directory', required=True)
	parser.add_argument('--force', help='Force new alignment and phylogeny', action='store_true')
	parser.add_argument('--compressed', help='Compress number of displayed sequences', action='store_true')
	parser.add_argument('--all', help='Create visualizations for all gene analyses in directory', action='store_true')
	args=parser.parse_args()

	return args



def read_db(context_file):

	print('Reading db...')
	#Todo: Add cluster size to tree 
	#Create list of gene ids from the flanking region file
	if args.compressed==True:
		arg_ids=[line.split('__')[1] for line in open(context_file[0]+'.centroids', 'r') \
			if line.startswith('>')]

		clust_sizes={}
		for line in open(context_file[0]+'.clusters', 'r'):
			if line.startswith('C'):
				clust_sizes[line.split('\t')[8].split('__')[1]]={}
				clust_sizes[line.split('\t')[8].split('__')[1]]=line.split('\t')[2]
	else:
		arg_ids=[line.split('__')[1] for line in open(context_file[0], 'r') \
			if line.startswith('>')]

	#Read flanking regions into dict
	print('Reading flanking regions to dict...')
	seq_dict={}
	for line in open(context_file[0], 'r'):
		if line.startswith('>'):
			header=line.split('__')[1]
			seq=''
		else:
			seq+=line
			seq_dict[header]=seq

	#Connect to sqlite3 database
	connection=sqlite3.connect(args.db)
	cursor=connection.cursor()

	#Create dictionary to save information in
	gene_dict={}

	#Extract env_genes, perc_id, start, end, organism, genome and flanklen from db
	print('querying database...')
	for id in arg_ids:
		query="""
		SELECT \
		args.arg_name, \
		args.id, \
		args.perc_id, \
		args.arg_len, \
		args.uplen, \
		args.downlen,\
		args.frame, \
		genomes.id, \
		genomes.organism, \
		genomes.assembly,
		envs.env_name, \
		envs.id, \
		envs.env_start, \
		envs.env_end, \
		envs.env_strand \
		FROM args \
		INNER JOIN genomes ON genomes.id = args.genome_id \
		INNER JOIN envs ON envs.arg_id = args.id \
		WHERE args.id = ?
		"""
		
		cursor.execute(query, (int(id),))
		results=cursor.fetchall()

		#Save results to dict
		gene_dict[id]={}
		
		for result in results:
			#To avoid reassignment of same keys over and over, set if clause
			if not 'name' in gene_dict[id]:
				gene_dict[id]['name']=result[0]
				gene_dict[id]['perc_id']=result[2]
				gene_dict[id]['length']=result[3]
				gene_dict[id]['uplen']=result[4]
				gene_dict[id]['downlen']=result[5]
				gene_dict[id]['frame']=result[6]
				gene_dict[id]['organism']=result[8]
				gene_dict[id]['assembly']=result[9]

				if args.compressed==True:
					gene_dict[id]['cluster_size']=clust_sizes[id]
				gene_dict[id]['env_genes']={}

			gene_dict[id]['env_genes'][result[11]]={}
			gene_dict[id]['env_genes'][result[11]]['env_name']=result[10]
			gene_dict[id]['env_genes'][result[11]]['env_start']=result[-3]
			gene_dict[id]['env_genes'][result[11]]['env_stop']=result[-2]
			gene_dict[id]['env_genes'][result[11]]['env_strand']=result[-1]
	
	#Drop keys that are empty
	copy_dict={}
	for key, value in gene_dict.items():
		if len(gene_dict[key])>0:
			copy_dict[key]=value
	
	gene_dict=copy_dict

	print('Extracting sequences...')
	#add gene seqs to gene_dict
	for key, value in gene_dict.items():
		if str(gene_dict[key]['frame']).startswith('-'):
			gene_dict[key]['seq']=seq_dict[key][gene_dict[key]['downlen']:gene_dict[key]['downlen']\
			+gene_dict[key]['length']+1]
			gene_dict[key]['stop']=gene_dict[key]['downlen']
			gene_dict[key]['start']=gene_dict[key]['downlen']+gene_dict[key]['length']
		else:
			gene_dict[key]['seq']=seq_dict[key][gene_dict[key]['uplen']:gene_dict[key]['uplen']\
			+gene_dict[key]['length']+1]
			gene_dict[key]['start']=gene_dict[key]['uplen']
			gene_dict[key]['stop']=gene_dict[key]['uplen']+gene_dict[key]['length']

		for key2, value2 in value['env_genes'].items():

			value2['seq']=seq_dict[key][value2['env_start']:value2['env_stop']+1]

	#Write annotation metadata
	with open(args.o.rstrip('/')+'/'+'annotation_meta.csv', 'w') as outfile:
		
		#To sort genes in right order in annotation metafile, create list of gene starts and stops
		for key, value in gene_dict.items():

			outfile.write('_'*20+'\n'+gene_dict[key]['name']+'__'+str(key)+'__'+gene_dict[key]['organism']+\
			'__'+gene_dict[key]['assembly']+'\n'+'-'*20+'\n'
			)
			position_list=[]
			lines=[]
			position_list.append(gene_dict[key]['start']+gene_dict[key]['stop'])
			for key2, value2 in value['env_genes'].items():
				position_list.append(value2['env_start']+value2['env_stop'])

			sorted_pos=sorted(position_list)

			#Now write genes to the list, in the same order as they are in the list
			lines.append(gene_dict[key]['name']+'\t'+str(gene_dict[key]['start'])+'\t'+\
			str(gene_dict[key]['stop'])+'\t'+gene_dict[key]['seq']+'\n')

			for key2, value2 in value['env_genes'].items():
				lines.append(value2['env_name']+'\t'+str(value2['env_start'])+'\t'+\
				str(value2['env_stop'])+'\t'+value2['env_strand']+'\t'+value2['seq']+'\n')

			for element in sorted_pos:
				for line in lines:
					if int(element)==int(line.split('\t')[1])+\
					int(line.split('\t')[2]):
						outfile.write(line)

	return gene_dict

def cluster_seqs(context_file):

	#temporarily rewrite context file to remove spaces from names
	with open(context_file[0]+'_tmp', 'w') as outfile:
		for line in open(context_file[0], 'r'):
			if line.startswith('>'):
				outfile.write(line.replace(' ', '_'))
			else:
				outfile.write(line)

	#Use usearch to cluster the sequences at 95%
	if not os.path.exists(context_file[0]+'.sorted'):
		sort='usearch -sortbylength %s -fastaout %s' % (context_file[0]+'_tmp', context_file[0]+'.sorted')
		subprocess.call(sort, shell=True)

	if not os.path.exists(context_file[0]+'.centroids'):
		cluster='usearch -cluster_smallmem %s -id 0.95 -centroids %s -uc %s' % \
		(context_file[0]+'.sorted', context_file[0]+'.centroids', context_file[0]+'.clusters')
		subprocess.call(cluster, shell=True)

	#remove tmp_file
	os.remove(context_file[0]+'_tmp')

def cluster_profiles(profiles):
	
	#Create two lists - identical profiles and clusters (which will be both a list of lists)
	identical_profiles=[]
	clusters=[]

	#sort profiles by length
	profiles.sort(key=len)

	print('Identifying identical profiles...')
	#Compare each profile to each other profile to find the identical ones
	for i in range(len(profiles)):
		duplicates=[]

		duplicates.append(profiles[i])
		#compare all profiles to find the identical ones
		for profile in profiles:
			if not profile==profiles[i]:
				if sorted(profiles[i][1])==sorted(profile[1]):
					duplicates.append(profile)

		identical_profiles.append(duplicates)

	unique_profiles=[profiles[0] for profiles in identical_profiles]
	return unique_profiles




def align(context_file, unique_profiles):

	#If compressed==True, use centroid file for this
	if args.compressed==True:
		context_file[0]=context_file[0]+'.centroids'

	#Read in gene contexts
	seq_dict={}
	for line in open(context_file[0], 'r'):
		if line.startswith('>'):
			header=line.lstrip('>').rstrip('\n')
			seq=''
		else:
			seq+=line
			seq_dict[header]=seq

	#Get unique gene IDs from unique profiles
	unique_ids=[profile[0] for profile in unique_profiles]

	#Write unique context to extra file
	with open(context_file[0].replace('.fna', '.unique.fna'), 'w') as outfile:
		for key, value in seq_dict.items():
			if str(key.split('__')[1]) in unique_ids:
				outfile.write('>'+key.replace(' ', '_')+'\n'+value+'\n')

	#Create mafft alignment of unique sequences
	if not os.path.exists(context_file[0].replace('.fna', '.unique.aln')) or args.force == True:
	
		align='mafft --auto --reorder --thread 48 %s > %s' \
		% (context_file[0].replace('.fna', '.unique.fna'), context_file[0]\
		.replace('.fna', '.unique.aln'))

		subprocess.call(align, shell=True)
	
	#Create Phylogeny
	if not os.path.exists(context_file[0].replace('.fna', '.unique.tree')) or args.force == True:
		
		phylogeny='FastTree -gtr -nt < %s > %s' % \
		(context_file[0].replace('.fna', '.unique.aln'), \
		context_file[0].replace('.fna', '.unique.tree'))

		subprocess.call(phylogeny, shell=True)

def visualize_phylogeny(gene_dict, context_file):
	
	#Read in tree and assign additional information to each leaf
	t=Tree(context_file[0].replace('.fna', '.unique.tree'))

	for node in t.traverse():
		if node.is_leaf():
			id=node.name.split('__')[1]
			node.add_features(organism=gene_dict[id]['organism'])
			node.add_features(assembly=gene_dict[id]['assembly'])
			node.add_features(pident=gene_dict[id]['perc_id'])
			if args.compressed==True:
				node.add_features(cluster_size=gene_dict[id]['cluster_size'])

	#Create dictionary to append motifs to
	motif_dict={}

	#Create keyword lists to set gene color
	tnps=['iscr', 'transpos', 'tnp', 'insertion']
	ints=['inti', 'integrase', 'xerc', 'xerd']
	mobiles=['secretion', 'mobiliza', 'moba', 'mobb', 'mobc', 'mobl', 'plasmid', 'relaxase',\
		'conjugation', 'type iv']
	res=['lactam', 'aminoglyco', 'fluoroquinolo', 'tetracyclin', 'macrolid', 'carbapenem']

	print('decorating the tree...')
	#Create motifs for each gene associated with a leaf
	for leaf in t.traverse():
		if leaf.is_leaf():
		
			#traverse through environment genes for the respective sequence
			for key, value in gene_dict.items():
				motifs=[]

				#Assign start and end position for annotated gene
				gene_start=int(gene_dict[key]['start']/5)
				gene_end=int(gene_dict[key]['stop']/5)

				#Sort such that the greater number is end and smaller is start
				if gene_start>gene_end:

					gene_end=int(gene_dict[key]['start']/5)
					gene_start=int(gene_dict[key]['stop']/5)
				
				#Append motif for annotated gene
				gene_motif=[gene_start, gene_end,'()', \
				2, 10, 'red', 'red', 'arial|8|black|'+str(gene_dict[key]['name'])]

				if not str(gene_dict[key]['frame']).startswith('-'):
					ori_motif=[gene_end, gene_end+10, '>', 2, 10, 'red', 'red', None]

				else:
					ori_motif=[gene_start-10, gene_start, '<', 2, 10, \
					'red', 'red', None]

				motifs.extend([gene_motif, ori_motif])

				for key2, value2 in value['env_genes'].items():

					#Set color, default is orange
					color='orange'

					if any(keyword in value2['env_name'].lower() for keyword in tnps):
						color='violet'
					if any(keyword in value2['env_name'].lower() for keyword in ints):
						color='yellow'
					if any(keyword in value2['env_name'].lower() for keyword in mobiles):
						color='green'
					if any(keyword in value2['env_name'].lower() for keyword in res):
						color='DodgerBlue'
					if 'hypothetical' in value2['env_name']:
						color='grey'

					#Create motif for one env gene at a time and append to motif list
					motif=[int(value2['env_start']/5), int(value2['env_stop']/5), '()', 2, 10, color, color, \
					'arial|8|black|'+str(value2['env_name'])]

					#Set condition: If env gene != annotated gene, append motif
					arg_pos={i for i in range(int(gene_motif[0]), int(gene_motif[1]))}
					env_pos={i for i in range(int(motif[0]), int(motif[1]))}
					
					#Calculate overlap percentage between annotated gene and env gene
					total_overlap=float(len(arg_pos.intersection(env_pos)))
					overlap_perc=float(total_overlap/int(gene_dict[key]['length']/5))*100

					if overlap_perc<=70.0:
						motifs.append(motif)

						#Create additional motif to show gene orientation
						if value2['env_strand']=='+':
							ori_motif=[int(value2['env_stop']/5), int(value2['env_stop']/5+10), '>', 2, 10, \
							color, color, None]

						else:
							ori_motif=[int(value2['env_start']/5-10), int(value2['env_start']/5), '<', 2, 10, \
							color, color, None]

						motifs.append(ori_motif)

				#append motif lists to respective annotated gene in dict
				gene_dict[key]['motifs']=motifs

	#Set node style
	nst_plasmid=NodeStyle()
	nst_plasmid['bgcolor']='DarkSeaGreen'
	nst_other=NodeStyle()
	nst_other='AntiqueWhite'

	#Now annotate the tree with the motifs
	for node in t.traverse():
		if node.is_leaf():
			if 'plasmid' in node.organism:
				node.set_style(nst_plasmid)
			else:
				node.set_style(nst_other)
			
			seqFace=SeqMotifFace(seq=None, motifs=gene_dict[node.name.split('__')[1]]['motifs'], \
			seq_format='blank', gap_format='line')
			(t & node.name).add_face(seqFace, 1, 'aligned')

			#Create box showing gene percent id
			similarity=TextFace(node.pident, fsize=8)
			similarity.margin_top=2
			similarity.margin_bottom=2
			similarity.margin_left=2
			similarity.margin_right=2

			#Set box background color based on pident
			if node.pident<=90.0:
				similarity.background.color='DarkGoldenrod'
			elif 90.0<node.pident<=95.0:
				similarity.background.color='ForestGreen'
			elif 95.0<=node.pident:
				similarity.background.color='YellowGreen'

			node.add_face(similarity, column=2, position='aligned')
	
			#Create box showing cluster size
			if args.compressed==True:
				clust_box=TextFace(node.cluster_size, fsize=8)
				clust_box.margin_top=2
				clust_box.margin_bottom=2
				clust_box.margin_left=2
				clust_box.margin_right=2

				node.add_face(clust_box, column=3, position='aligned')

	#Return the annotated tree
	return t

def render_tree(tree, context_file):

	print('rendering tree...')
	print (tree)
	ts=TreeStyle()
	ts.tree_width=30
	if not args.compressed==True:
		tree.render(context_file[0].replace('fna', '_tree_annotated.pdf'), tree_style=ts)
	else:
		print(context_file[0].replace('.centroids', '_tree_compressed.pdf'))
		tree.render(context_file[0].replace('.centroids', '_tree_compressed.pdf'), tree_style=ts)
	print('tree rendered!')

def main():
	#Extract file containing flanking regions
	context_file=[args.o.rstrip('/')+'/'+file for file in os.listdir(args.o) if file.endswith('_contexts.fna')]

	if args.compressed==True:
		cluster_seqs(context_file)
	#Collect gene entries from the database
	gene_dict=read_db(context_file)

	#Create a list of profiles for each genome
	profiles=[(key,[value2['env_name'] for key2, value2 in value['env_genes'].items()]) for key, value\
	in gene_dict.items()]

	#Identify unique profiles
	if args.compressed==True:
		unique_profiles=cluster_profiles(profiles)
	else:
		unique_profiles=profiles

	#align unique profiles and create phylogeny
	align(context_file, unique_profiles)

	#Use ete3 to annotate and visualize the tree
	tree=visualize_phylogeny(gene_dict, context_file)

	#render tree
	render_tree(tree, context_file)

if __name__=='__main__':
	args=parse_arguments()
	if not args.all==True:
		main()
	else:
		outfiles=[element[0] for element in os.walk(args.o)]
		outdirs=outfiles[1:]
		for outdir in outdirs:
			args.o=outdir
			try:
				main()
			except:
				print(f'Something went wrong when trying to process {args.o}!')
