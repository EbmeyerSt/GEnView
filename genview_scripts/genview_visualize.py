import sys, os, subprocess, argparse, sqlite3, timeit, csv, math, re, multiprocessing
from argparse import RawTextHelpFormatter

"""Tutorial command:
	genview-visualize -gene PER-1 -db /path/to/genview_database.db -id 80 
"""

def parse_arguments():
	descr='\nExtract, visualize and annotate genes and genetic environments from genview database\n'
	parser=argparse.ArgumentParser(description=descr.replace("'", ''), formatter_class=RawTextHelpFormatter)
	parser.add_argument('-gene', help='name of gene/orf to extract and visualize', required=True)
	parser.add_argument('-db', help='genview database created by genview-create-db', required=True)
	parser.add_argument('-id', help='percent identity threshold for genes to extract', required=True)
	parser.add_argument('-nodes', help='should nodes be connected to genome with solid line (solid), connected by dashed line (dash) or no connection (none)', default='solid', type=str)
	parser.add_argument('-taxa', help='list of genera and/or species to extract\nBy default all taxa are extracted', default='False', nargs='+')
	parser.add_argument('--force', help='Force new alignment and phylogeny', action='store_true')
	parser.add_argument('--compressed', help='Compress number of displayed sequences, helpful with large number of identical sequences', action='store_true')
	parser.add_argument('--custom_colors', help='path to file containing RGB color codes for gene color customization', default='False')
	parser.add_argument('--log', help='Write log file', action='store_true')
	args=parser.parse_args()

	return args

def extract():

	#Connect to db
	connection=sqlite3.connect(args.db)
	cursor=connection.cursor()
	
	#Join tables based on genome ids
	query=f"""SELECT \
	args.arg_name, \
	args.id, \
	genomes.organism \
	FROM args \
	INNER JOIN genomes \
	ON args.genome_id=genomes.id \
	WHERE args.arg_name \
	LIKE \'{args.gene+"%"}\' \
	AND args.perc_id >= {float(args.id)} \
	"""

	if args.taxa!='False':
		i=0
		for taxon in args.taxa:
			print(taxon)
			i+=1
			if len(args.taxa)>1:
				if i==1 and args.taxa>1:
					query+='AND ('
					query+=f'genomes.organism LIKE \'{taxon+"%"}\' '
				elif 1<i<len(args.taxa):
					query+=f'OR genomes.organism LIKE \'{taxon+"%"}\' '
				else:
					query+=f'OR genomes.organism LIKE \'{taxon+"%"}\') '
			else:
				query+=f'AND genomes.organism LIKE \'{taxon+"%"}\' '


	cursor.execute(query)

	results=cursor.fetchall()

	#write create result directory if not exists
	if not os.path.exists(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis'):
		os.makedirs(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis')

	elif os.path.exists(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis') and not \
	args.force==True:
		print(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis'+'exists already!\nUse --force to overwrite previous results.\n')
		sys.exit()
	
	#Using the id, go back to the file containing the flanking regions and extract the ones matching the ids 
	flank_dict={}
	for line in open(os.path.dirname(args.db).rstrip('/')+'/all_flanks.csv_tmp', 'r'):
		flank_dict[line.split('\t')[0]]=line.split('\t')[1]

	arg_ids=[result[1] for result in results]

	#write results to output directory 
	with open(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis/'+args.gene+'_contexts.fna', 'w') as outfile:	
		for result in results:
			outfile.write('>'+result[0]+'__'+str(result[1])+'__'+result[2]\
			+'\n'+flank_dict[str(result[1])]+'\n')
		
	if args.log==True:
		if len([line for line in open(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis/'+args.gene+'_contexts.fna', 'w') if line.startswith('>')])>=2:
			log_lines.append('Target sequences extracted...\n')
		else:
			log_lines.append('No/too few target sequences extracted...FAILED\n')
		write_log()

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

			#Replace prodigal orf with target gene annotation
			if result[4]-200<result[-3]<result[4]+200 and result[4]\
			+result[3]-200<result[-2]<result[4]+result[3]+200:

				gene_dict[id]['env_genes'][result[11]]['env_name']=result[0]
			else:
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

	#Create keyword lists to set gene color
	tnps=['iscr', 'transpos', 'tnp', 'insertion']
	ints=['inti', 'integrase', 'xerc', 'xerd']
	mobiles=['secretion', 'mobiliza', 'moba', 'mobb', 'mobc', 'mobl', 'plasmid', 'relaxase',\
		'conjugation', 'type iv']
	res=['lactam', 'aminoglyco', 'fluoroquinolo', 'tetracyclin', 'macrolid', 'carbapenem']

	#Create list header for tmp visualization file
	vis_list = [['key', 'name', 'organism', 'assembly', 'gene', 'start', 'stop', 'direction', 'group', 'sequence', 'gvid']]

	with open(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis'.rstrip('/')+'/'+'annotation_meta.csv', 'w') as outfile:
		
		#To sort genes in right order in annotation metafile, create list of gene starts and stops
		for key, value in gene_dict.items():

			outfile.write('_'*20+'\n'+gene_dict[key]['name']+'__'+str(key)+'__'+gene_dict[key]['organism']+\
			'__'+gene_dict[key]['assembly']+'\n'+'-'*20+'\n')

			position_list=[]
			lines=[]
			position_list.append(gene_dict[key]['start']+gene_dict[key]['stop'])
			for key2, value2 in value['env_genes'].items():
				position_list.append(value2['env_start']+value2['env_stop'])

			sorted_pos=sorted(position_list)

			#Now write genes to the list, in the same order as they are in the list
			lines.append(gene_dict[key]['name']+'\t'+str(gene_dict[key]['start'])+'\t'+\
			str(gene_dict[key]['stop'])+'\t'+gene_dict[key]['seq']+'\ttarget\n')

			gvid = 1
			t_key = 1
			for key2, value2 in value['env_genes'].items():
				#add group (e.g transposon,integron, etc.)
				group_assigned=0
				if not args.custom_colors==True:
					if any(keyword.lower() in value2['env_name'].lower() for keyword in tnps):
						group='transposase'
						group_assigned=1
					if any(keyword.lower() in value2['env_name'].lower() for keyword in ints):
						group='integron'
						group_assigned=1
					if any(keyword.lower() in value2['env_name'].lower() for keyword in mobiles):
						group='mobile'
						group_assigned=1
					if any(keyword.lower() in value2['env_name'].lower() for keyword in res):
						group='resistance'
						group_assigned=1
					if 'hypothetical' in value2['env_name'].lower():
							group='hypothetical'
							group_assigned=1

					if not group_assigned==1:
						group='misc'
				
				else:
					cust_groups=[line.rstrip('\n') for line in open(args.custom_colors, 'r')]	
					for cust_group in cust_groups:
						if any(keyword.lower() in value2['env_name'].lower() for keyword in \
						cust_group.split('\t')[1]):
							group=cust_group.split('\t')[0]
							group_assigned=1

					if not group_assigned==1:
						group='misc'
					
				lines.append(value2['env_name']+'\t'+str(value2['env_start'])+'\t'+\
				str(value2['env_stop'])+'\t'+value2['env_strand']+'\t'+value2['seq']+'\t'+group+'\n')

				#Vis list
				vis_list.append([str(key), gene_dict[key]['name'], gene_dict[key]['organism'], gene_dict[key]['assembly'], value2['env_name'], value2['env_start'], value2['env_stop'], value2['env_strand'], group, value2['seq'], str(key) + '.' + str(gvid)])
				if t_key == key:
					gvid += 1
				else:
					gvid = 1
					t_key = key

			for element in sorted_pos:
				for line in lines:
					if int(element)==int(line.split('\t')[1])+\
					int(line.split('\t')[2]):
						outfile.write(line)
	
	#Create tmp file for visualization
	with open(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis'.rstrip('/')+'/'+'visualization_meta.csv', 'w') as outfile:
		write = csv.writer(outfile) 
		write.writerows(vis_list) 

	if args.log==True:
		if len(gene_dict)>=2:
			log_lines.append('Genetic environments extracted and visualization metadata created...\n')
		else:
			log_lines.append('Could not extract genetic environments/create visualization metadata...FAILED\n')
		write_log()

	return gene_dict


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

	if args.log==True:
		if len(unique_profiles)>=2:
			log_lines.append('Unique profiles identified...\n')
		else:
			log_lines.append('Unique profiles could not be identified...\n')
		write_log()

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
		print('Removing special characters in sequence headers...')
		for key, value in seq_dict.items():
			if str(key.split('__')[1]) in unique_ids:
				new_entry='>'+key.replace(' ', '_').replace('(', '').replace(')', '')\
				.replace('"','').replace("'", "").replace('*', '')\
				.replace('/', '').replace('|', '')+'\n'+value+'\n'
				outfile.write(new_entry)


	#Create mafft alignment of unique sequences
	if not os.path.exists(context_file[0].replace('.fna', '.unique.aln').replace('(', '\(').replace(')', '\)')) or args.force == True:
	
		align='mafft --auto --reorder --thread %r  %s > %s' \
		% (multiprocessing.cpu_count(), context_file[0].replace('.fna', '.unique.fna').replace('(', '\(').replace(')', '\)'), context_file[0]\
		.replace('.fna', '.unique.aln').replace('(', '\(').replace(')', '\)'))

		subprocess.call(align, shell=True)

	if args.log==True:
		if os.path.exists(context_file[0].replace('.fna', '.unique.aln').replace('(', '\(').replace(')', '\)')) and\
		 os.path.getsize(context_file[0].replace('.fna', '.unique.aln').replace('(', '\(').replace(')', '\)'))>0:
			log_lines.append('Sequences aligned...\n')
		else:
			log_lines.append('Sequence alignment failed...FAILED\n')
		write_log()
	
	#Create Phylogeny
	if not os.path.exists(context_file[0].replace('.fna', '.unique.tree').replace('(', '\(').replace(')', '\)')) or args.force == True:
		
		phylogeny='FastTree -gtr -nt < %s > %s' % \
		(context_file[0].replace('.fna', '.unique.aln').replace('(', '\(').replace(')', '\)'), \
		context_file[0].replace('.fna', '.unique.tree').replace('(', '\(').replace(')', '\)'))

		subprocess.call(phylogeny, shell=True)

		if args.log==True:
			if os.path.exists(context_file[0].replace('.fna', '.unique.tree').replace('(', '\(').replace(')', '\)')) and\
			 os.path.getsize(context_file[0].replace('.fna', '.unique.tree').replace('(', '\(').replace(')', '\)'))>0:
				log_lines.append('Phylogeny created...\n')
			else:
				log_lines.append('Phylogeny could not be created...FAILED\n')
			write_log()


###HTML VERSION###
#Create a text index file of the phylogenetic tree
def make_text(final_tree):
	if os.path.exists(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis' + "/index_tree.txt"):
		os.remove(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis' + "/index_tree.txt")
	f = open(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis' + "/index_tree.txt", "a")
	for row in final_tree:
		f.write(str(row))
		f.write('\n')
	f.close()

#From branch name get direct parent branch
def get_parent(child):
	t = list(child)
	f = 0
	for v in reversed(t):
		if v.isnumeric():
			f += 1
		else:
			f += 1
			break
	k = ''.join(t)
	k = k[:-f]
	return k

#Getting last number of branch name
def get_last_number(string):
	t = list(string)
	res = []
	f = 0
	for v in reversed(t):
		if v.isnumeric():
			res.insert(0, v)
		else:
			break
	k = ''.join(res)
	return k

#Create a series of indexes, where each row represents a horizontal line by 1's as lines and 0's as spaces.
def create_tree_index(tree_file_path):
	with open(tree_file_path) as tree:
		completed = False
		array = []
		final_tree = []
		#Extract first line of tree file which contains the tree information
		for line in tree: 
			#Remove all bootstrap values from data      
			i = 0
			line = list(line)
			while True:
				if line[i] == ')' and line[i+6] == ':':
					del line[i+1:i+6]
				i += 1
				if i == len(line)-6:
					break 
			line = ''.join(line)
			used_ids = []
			array.extend(line)
			current_id = '1'
			#Iterate through each character in tree file line
			i = 0
			while i < len(array):
				#"(" indicates a new branch and therefore creates a branch array
				if array[i] == ')' and array[i+1] == ';':
					completed = True
					break
				elif array[i] == "(":   
					#Adjust id's and add to ID list
					id_ext = 1
					while True:
						tmp_id = current_id + '.' + str(id_ext)
						if tmp_id in used_ids:
							id_ext += 1
						else:
							current_id = tmp_id
							break   
					if len(used_ids) == 0:
						final_tree.append([current_id])
					else:
						final_tree.append([current_id]) 
					used_ids.append(current_id) 
				#Skip bootsrtap values
				elif array[i] == ')' and array[i+6] == ':':
					i += 5
				elif array[i] == ":" and array[i-1] != ")":
					#Adjust id's and add to ID list
					id_ext = 1
					while True:
						tmp_id = current_id + '.' + str(id_ext)
						if tmp_id in used_ids:
							id_ext += 1
						else:
							break
					final_tree.append([tmp_id])
					used_ids.append(tmp_id)    
				#Marks end of branch and all sub branches so length of this branch must be added
				elif array[i] == ":" and array[i-1] == ")":              
					current_id = get_parent(current_id)
				i += 1
			if completed == True:
				break

	result = []
	for item in final_tree:
		l = len(item[0])/2+0.5
		result.append(l)
	res = max(result)
	#Organize tree branches
	p = 1
	new_tree = []
	while p <= res:
		for item in final_tree:
			if int((len(item[0])/2+0.5)) == p:
				parent = get_parent(item[0])
				try:
					n = new_tree.index([parent])    
				except:
					n = 0                            
				last_letter = int(get_last_number(item[0]))
				if last_letter == 1:
					new_tree.insert(n, item)
				elif last_letter == 2:
					new_tree.insert(n+1, item)
				else:
					last_sibling = parent + '.' + str(int(last_letter)-1)
					n = new_tree.index([last_sibling]) 
					new_tree.insert(n+1, item)
					new_tree.insert(n+1, ['-'])
		p += 1            

	#Create list of orderd organism names
	with open(tree_file_path) as tree:
		final_order = []
		for line in tree: 
			line = list(line)
			i = 0
			while i < (len(line)-1):
				if line[i] != '(' and line[i-1] == '(':
					name = ''
					while True:
						name += line[i]
						i += 1
						if line[i] == ':':
							final_order.append(int(name.split('__')[1]))
							break
				elif line[i-1] == ',' and line[i] != '(':
					name = ''
					while True:
						name += line[i]
						i += 1
						if line[i] == ':':
							final_order.append(int(name.split('__')[1]))
							break
				i += 1
	
	space = [0]*5
	lne = [1]*5
	#Assign distances to branches and nodes
	for item in new_tree:
		#Check for if branch or node
		again = 0
		for branch in new_tree:
			if branch[0].startswith(item[0]):
				again += 1
		if item[0] == '-': 
			item.extend(space)
		elif item[0] == '1.1':
			item.extend(space)
		elif again <= 1:
			#If node
			a = int((len(item[0])/2+0.5))-2
			pre = space*a
			item.extend(pre)
			item.extend(lne)
		else:
			#if branch
			a = int((len(item[0])/2+0.5))-2
			pre = space*a
			item.extend(pre)
			item.extend(lne)

	#Add white spacer before and after lists
	new_tree.insert(0, [0]*len(new_tree[0]))
	new_tree[0][0] = '-'
	new_tree.append([0]*len(new_tree[0]))
	new_tree[len(new_tree)-1][0] = '-'

	#Adding spacers between list elements where required
	i = 0 
	while i < len(new_tree[:-1]):
		parent_i = get_parent(new_tree[i][0])
		parent_n = get_parent(new_tree[i+1][0])
		if parent_i == parent_n:
			spacer = [0]*(len(new_tree[i]))
			spacer[0] = '-'
			new_tree.insert(n, spacer)
		i += 1
	long = len(max(new_tree, key=len)) + 1
	z = 0
	while z < len(new_tree):
		if len(new_tree[z]) <= long:
			val = long - len(new_tree[z])
			p = 1
			while p <= val:
				new_tree[z].append(0)
				p += 1
		z += 1
	i = 0
	while i < len(new_tree):
		new_tree[i].insert(1, 0)
		i += 1
	i = 0
	while i < len(new_tree):
		new_tree[i].insert(1, 0)
		i += 1

	#Making lengths of list names for each branch and node the same length
	lenghts = []
	for row in new_tree:
		lenghts.append(len(str(row[0])))
	m = max(lenghts)
	i = 0
	while i < len(new_tree):
		add = m - len(str(new_tree[i][0]))
		if add > 0:
			new_tree[i][0] = str(new_tree[i][0]) + '-'*add
		i += 1

	#Get max number of levels
	mx = 0
	for branch in new_tree[1:-1]:
		name = branch[0].strip('-')
		cur_level = name.count('.')
		if cur_level > mx:
			mx = cur_level
	#Center horizontal lines
	shift_index = []
	while mx >= 1:
		cur_row = 1
		for b in new_tree[1:-1]:
			name = b[0].strip('-')
			cur_level = name.count('.')
			if cur_level == mx:
				if len(name) > 0:
					#Intercept here to remove and replace 1's
					r_start = (cur_level-1)*5 + 3
					r_end = r_start + 4
					again = 0
					#See if b is branch or node
					pos_child = name + '.1'
					for branch in new_tree[1:-1]:
						if branch[0].startswith(pos_child):
								again += 1
					#IF branch get all direct children and count distance.
					if again > 0:
						ext = 1
						tmp_row = 1
						for branch in new_tree[1:-1]:
							child = name + '.' + str(ext)
							if branch[0].strip('-') == child:
								if ext == 1:
									if child in shift_index:
										shift = shift_index.index(child)
										first_row = shift_index[shift+2]
									else:
										first_row = tmp_row                                        
									ext += 1
								else:
									if child in shift_index:
										shift = shift_index.index(child)
										last_row = shift_index[shift+2]
									else:
										last_row = tmp_row
									ext += 1
							tmp_row += 1 
						#Remove existing 1's in range and replace in midpoint
						rel_mid = math.ceil(((last_row-first_row)/2))
						actual_mid_row = first_row+rel_mid                        
						if cur_row != actual_mid_row:
							shift_index.append(name)
							shift_index.append(cur_row)
							shift_index.append(actual_mid_row)
							i = r_start
							while i <= r_end:
								new_tree[cur_row][i] = 0
								i += 1                            
							i = r_start
							while i <= r_end:
								new_tree[actual_mid_row][i] = 1
								i += 1
						#Add vertical lines:  
						ver_line = r_end + 1
						k = first_row
						while k <= last_row:
							new_tree[k][ver_line] = 1
							k += 1
			cur_row += 1
		mx -= 1
			
	#Add stipple lines to nodes
	if args.nodes != 'none':
		if args.nodes == 'dash' or args.nodes == 'solid':
			i = 1
			while i < len(new_tree[:-1]):
				n = len(new_tree[0])-3
				while n >= 0:
					if new_tree[i][n-2] == 1:
						break
					else:
						new_tree[i][n] = 1
						new_tree[i][n-1] = 1
						if args.nodes == 'solid':
							new_tree[i][n-2] = 1
					if args.nodes == 'dash':
						n -= 4
					elif args.nodes == 'solid':
						n -= 1
				i += 2

	return [new_tree, final_order]


def create_tree_lines(tree_index, y_val, size, width=300):
   #Create string of phylogenetic tree HTML
	tree_string = ''
	x_sep = round(width/len(tree_index[0][1:]),2)
	#Calculate height ==== 
	height = len(tree_index)*y_val + y_val
	if size == "l_tree":
		tree_string += '<svg class="tree '+ size +'" width="'+ str(width) +'" height="'+ str(height) +'">'
	else:
		tree_string += '<svg class="hidden tree '+ size +'" width="'+ str(width) +'" height="'+ str(height) +'">'    
	#Get horizontal lines
	y = y_val + 1
	for row in tree_index[1:-1]:
		x1 = 0
		x2 = 0
		xval = 0
		i = 1
		while i < len(row):
			n = i + 1
			p = i - 1
			if row[i] == 1 and row[n] == 1 and row[p] == 0:
				x1 = xval
			elif row[i] == 1 and row[n] == 0 and row[p] == 1:
				x2 = xval
				tree_string += '<line class="hor" x1="'+ str(x1) +'" y1="'+ str(y) +'" x2="'+ str(x2) +'" y2="'+ str(y) +'" stroke="black" stroke-width="1"/>'
			i += 1
			xval += x_sep
		y += y_val

	#get vertical lines
	i = 1
	xval = 0
	while i < len(tree_index[0]):
		y1 = 0
		y2 = 0
		y = y_val + 1
		k = 1    
		while k < (len(tree_index)-1):
			if tree_index[k][i] == 1 and tree_index[k+1][i] == 1 and tree_index[k-1][i] == 0:
				y1 = y
			elif tree_index[k][i] == 1 and tree_index[k+1][i] == 0 and tree_index[k-1][i] == 1:
				y2 = y
				if y1 > 0 and y2 > 0:
					tree_string += '<line class="ver" x1="'+ str(xval) +'" y1="'+ str(y1) +'" x2="'+ str(xval) +'" y2="'+ str(y2) +'"/>'
			k += 1    
			y += y_val  
		xval += x_sep    
		i += 1
	tree_string += '</svg>'
	return tree_string

def create_html(tree_index, meta_file): 
	#Create Phylo Tree
	tree = tree_index[0]
	tree_string = ''
	tree_string += create_tree_lines(tree, 16, 'l_tree')
	tree_string += create_tree_lines(tree, 9, 'm_tree')
	tree_string += create_tree_lines(tree, 4, 's_tree')

	#Get max sequence length
	with open(meta_file, newline = '') as org_file:
		reader = list(csv.reader(org_file))[1:]
		values = []
		org = []
		key = int(reader[0][0])
		for row in reader:
			if key == int(row[0]):
				org.append([int(row[0]), int(row[6])])                
			else:
				key = int(row[0])
				values.append(org)
				org = []
				org.append([int(row[0]), int(row[6])]) 
		values.append(org)
	#Get max sequence length and max length of each organism
	vals = []
	final = []
	tmp_vals = []
	key = values[0][0][0]
	for item in values:
		for i in item:
			vals.append(i[1])
			if key == int(i[0]):
				tmp_vals.append(i[1])
			else:
				final.append([key, max(tmp_vals)])
				tmp_vals = []
				key = int(i[0])
				tmp_vals.append(i[1])
	final.append([int(i[0]), max(tmp_vals)])
	factor = 1000/max(vals)

	
	name_string = ''
	name_string += '<div class="names" style="grid-column-start: 2; grid-column-end: 3;grid-row-start: 1;grid-row-end: 2;"><div class="grid_container">'
	#Create names
	with open(meta_file, newline = '') as file:
		reader = list(csv.reader(file))[1:]
		start = 1
		stop = 5
		key = True
		for n in tree_index[1]:
			for row in reader:
				if key == True:
					if int(row[0]) == int(n):
						key = False
						name_string += '<div style="grid-column: 1 / 1;grid-row:' + str(start) + '/' + str(stop) + ';"><p>' + row[2].replace(' ', '_') + '</p></div>'
						start += 5
						stop += 5
			key = True
	name_string += ' </div></div>'

	tnps=['iscr', 'transpos', 'tnp', 'insertion']
	ints=['inti', 'integrase', 'xerc', 'xerd']
	mobiles=['secretion', 'mobiliza', 'moba', 'mobb', 'mobc', 'mobl', 'plasmid', 'relaxase', 'conjugation', 'type iv']
	res=['lactam', 'aminoglyco', 'fluoroquinolo', 'tetracyclin', 'macrolid', 'carbapenem']

	string = ""
	string += '<div class="sequences" style="grid-column-start: 3; grid-column-end: 4;grid-row-start: 1;grid-row-end: 2;">'


	#Create genetic environment
	with open(meta_file, newline = '') as file:
		reader = list(csv.reader(file))[1:]
		key = 0
		id = 0
		x = False
		for n in tree_index[1]:
			for row in reader:
				if int(row[0]) == int(n):
					if row[0] != key:
						for i in final:
							if int(i[0]) == int(row[0]):
								high = i[1]
								break
						if x == True:
							string += '</div>'
						string += '<div class="grid_container">'
						string += '<div id="' + str(id) + '_line" class="line" style="grid-column: 1 / ' + str(math.ceil(int(high)*factor)) + ';grid-row: 2 / 4;background-color: rgba(0, 0, 0);opacity: 0.2;"></div>'
						with open(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis'.rstrip('/')+'/../'+'all_flanks.csv_tmp') as file:	
							seq_reader = csv.reader(file, delimiter='\t')
							for row_n in seq_reader:
								if int(row_n[0]) == int(row[0]):
									string += '<div id="'+ str(id) +'_line_info" class="hidden info_box"><button class="exit">X</button><p><strong>GVID:</strong>' + row[0] + '</p><p><strong>Organism:</strong>' + row[2] + '</p><p><strong>Accession:</strong> ' + row[3] + ' </p><textarea id="'+ str(id) +'_gene_sequence">' + '&#62' + row[2].replace(' ', '_') + '_' + row[3] +'&#13;&#10;'+ row_n[1] +'</textarea><button id="'+ str(id) +'" class="copy btn";">Copy</button></div>'           
					indx = 4
					if args.custom_colors=='False':	
						if args.gene.lower() in row[4].lower():
							color='red'
							color = 'rgb(201, 0, 0)'
						elif any(keyword in row[4].lower() for keyword in tnps):
							color='purple'
							color = 'rgb(201, 50, 255)'
						elif any(keyword in row[4].lower() for keyword in ints):
							color='yellow'
							color = 'rgb(253, 228, 0)'
						elif any(keyword in row[4].lower() for keyword in mobiles):
							color='green'
							color = 'rgba(0, 255, 13)'
						elif any(keyword in row[4].lower() for keyword in res):
							indx = 4
							color='DodgerBlue'
							color = 'rgba(0, 183, 255)'
						elif 'hypothetical' in row[4].lower():
							color='grey'
							color = 'rgb(153, 153, 153)'
						else:
							color='orange'
							color = 'rgb(244, 153, 50)'
					else:
						
						cust_groups=[line.rstrip('\n') for line in open(args.custom_colors, 'r')]	
						color='rgb(244, 153, 50)'
						for cust_group in cust_groups:
							if any(keyword.lower() in row[4].lower() for keyword in \
							cust_group.split('\t')[1].split(',')):
								color=cust_group.split('\t')[2]

					start = math.ceil(int(row[5])*factor)
					end = math.ceil(int(row[6])*factor)

					if row[7] == '+':
						string += '<div id="'+ str(id) +'_gene" class="R all" style="grid-column: '+ str(start) +' / '+ str(end) +';grid-row: 1 / 5;background-color: '+ color +';margin-right: 15px;opacity: 1.0;"><p class="gene_annot" style="font-size:8px;">' + row[indx] + '</p></div>'
						string += '<div class="AR right" id="'+ str(id) +'_gene_arrow" style="grid-column: '+ str(end) +' / '+ str(end) +';grid-row: 1 / 5;--my-color-var: '+ color +';opacity: 1.0;"></div>'
					elif row[7] == '-':
						string += '<div class="AL left" id="'+ str(id) +'_gene_arrow" style="grid-column: '+ str(start) +' / '+ str(start) +';grid-row: 1 / 5;--my-color-var: '+ color +';opacity: 1.0;"></div>'
						string += '<div id="'+ str(id) +'_gene" class="L all" style="grid-column: '+ str(start) +' / '+ str(end) +';grid-row: 1 / 5;background-color: '+ color +';margin-left: 15px;opacity: 1.0;"><p class="gene_annot" style="font-size:8px;">' + row[indx] + '</p></div>'

					string += '<div id="'+ str(id) +'_gene_info" class="hidden info_box"><button class="exit">X</button><p><strong>GVID:</strong>' + row[0] + '</p><p><strong>Organism:</strong>' + row[2] + '</p><p><strong>Gene:</strong> '+ row[4] +'</p><p><strong>Position:</strong> '+ str(row[5]) + ' - ' + str(row[6]) + ' (' + row[7] + ')' + '</p><textarea id="'+ str(id) +'_gene_sequence">' + '&#62' + row[2].replace(' ', '_') + '_GVID.' + str(row[10]) +'&#13;&#10;'+ row[9] +'</textarea><button id="'+ str(id) +'" class="copy btn";">Copy</button></div>'           

					x = True
					key = row[0]
					id += 1   
		string += '</div>'  

	html = """
	<html>
	<head>
		<title>GenView Output</title>
		<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
		<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.4.0/jspdf.umd.min.js"></script>
		<script src="https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.3.2/html2canvas.min.js"></script>
	</head>
	<style>
		body,html{margin:0;height:100%}.space{position:relative;height:50px}.navigation{top:0;border:none;outline:0;position:fixed;width:100%;height:50px;background-color:#838383;margin-bottom:60px;display:inline-block;overflow:hidden;z-index:100;box-shadow:0 5px 5px rgb(131,131,131,.5)}.options{position:fixed;top:-350;background-color:rgb(94,94,94,.8);width:25%;padding:10px;z-index:1000;right:0}.options input{display:block;width:250px}label{color:#fff}.active{background-color:#46b5c9!important}.navigation button{position:relative;width:100px;background-color:#838383;display:inline-block;border:none;outline:0;transition:.3s;color:#fff;cursor:pointer;height:100%}.navigation button:hover{background-color:#46b5c9;transition:.3s}.navigation button:active{background-color:#838383;transition:.3s}line{stroke:#000;stroke-width:1}.phylo_div{height:100%}.phylo_div div{display:inline-block;height:10px;padding:0;margin:0 auto}.large_grid_container{display:grid;height:100%;grid-template-columns:300px auto}.grid_container{margin-top:2px;display:grid;grid-row-gap:2px;grid-template-columns:repeat(1000,1fr)}.grid_container p{margin-top:6px;margin-bottom:0;padding:none}.all{transition:.3s;text-align:center;margin:none;z-index:2;transition:.3s;height:30px;overflow:hidden;text-overflow:ellipsis;cursor:pointer}.all p{color:#000;position:relative;overflow:hidden;white-space:nowrap;text-overflow:clip}.line{z-index:1;transition:.3s;cursor:pointer;height:5px}.left{position:relative;outline:0;width:-webkit-calc(100% - 15px);width:-moz-calc(100% - 15px);width:calc(100% - 15px);z-index:2;transition:.3s}.left:after{content:'';position:absolute;top:0;border-top:15px solid transparent;border-right:15px solid var(--my-color-var);border-bottom:15px solid transparent;width:0}.left_s{position:relative;outline:0;width:-webkit-calc(100% - 15px);width:-moz-calc(100% - 15px);width:calc(100% - 15px);z-index:2;transition:.3s}.left_s:after{content:'';position:absolute;top:0;border-top:3px solid transparent;border-right:6px solid var(--my-color-var);border-bottom:3px solid transparent;width:0}.left_m{position:relative;outline:0;width:-webkit-calc(100% - 15px);width:-moz-calc(100% - 15px);width:calc(100% - 15px);z-index:2;transition:.3s}.left_m:after{content:'';position:absolute;top:0;border-top:8px solid transparent;border-right:8px solid var(--my-color-var);border-bottom:8px solid transparent;width:0}.right{position:relative;outline:0;width:-webkit-calc(100% - 15px);width:-moz-calc(100% - 15px);width:calc(100% - 15px);left:-15px;z-index:2;transition:.3s}.right:after{content:'';position:absolute;top:0;border-top:15px solid transparent;border-left:15px solid var(--my-color-var);border-bottom:15px solid transparent;width:0}.right_s{position:relative;outline:0;width:-webkit-calc(100% - 15px);width:-moz-calc(100% - 15px);width:calc(100% - 15px);left:-6px;z-index:2;transition:.3s}.right_s:after{content:'';position:absolute;top:0;border-top:3px solid transparent;border-left:6px solid var(--my-color-var);border-bottom:3px solid transparent;width:0}.right_m{position:relative;outline:0;width:-webkit-calc(100% - 15px);width:-moz-calc(100% - 15px);width:calc(100% - 15px);left:-8px;z-index:2;transition:.3s}.right_m:after{content:'';position:absolute;top:0;border-top:8px solid transparent;border-left:8px solid var(--my-color-var);border-bottom:8px solid transparent;width:0}.hidden{display:none}.info_box{position:fixed;top:0;right:-350px;width:300px;height:100%;text-align:center;z-index:100;background-color:rgba(194,192,192,.9)}.info_box p{margin-right:auto;margin-left:auto;text-align:left;width:80%;overflow:visible;word-break:break-all}.info_box textarea{margin-right:auto;margin-left:auto;text-align:left;width:80%;height:300px;overflow:visible;resize:vertical}.btn{position:relative;width:100px;height:30px;background-color:#838383;display:inline-block;border:none;outline:0;margin-top:15px;transition:.3s;color:#fff;cursor:pointer}.btn:hover{background-color:#46b5c9;transition:.3s}.btn:active{background-color:#838383;transition:.3s}.exit{position:absolute;height:25px;width:25px;outline:0;top:10px;right:10px;border:none;border-radius:100%;background-color:#838383;text-align:center;transition:.3s;color:#fff;cursor:pointer}.exit:hover{background-color:red;transition:.3s}.exit:active{background-color:#860000;transition:.3s}#statusmessage{position:fixed;top:0;padding:0 5px;line-height:25px;background-color:#68d1f1;margin-top:-25px;font-weight:700;right:0;text-align:center;z-index:500;width:300px}.names{overflow:hidden;padding-left:5px;padding-right:5px;margin-top:7px;width:200px}.names p{padding:0 auto;margin:0 auto;height:28px;font-size:12px}
	</style>
	<script>
		function reveal(){$(document).on("click",".all,.line",function(){$(".options").css({top:-300}),$(".options").removeClass("open");var e="#"+$(this).attr("id")+"_info";$(e).hasClass("hidden")?($(".info_box").addClass("hidden"),$(e).removeClass("hidden"),$(e).css({right:"0"})):$(".info_box").css({right:"-350"}).addClass("hidden")})}function copy_to_clipboard(e){document.getElementById(e).select(),document.execCommand("copy")}function saveAs(e,s){var n=document.createElement("a");"string"==typeof n.download?(n.href=e,n.download=s,document.body.appendChild(n),n.click(),document.body.removeChild(n)):window.open(e)}$(document).on("change","#gene_annot",function(){$(this).is(":checked")?$(".gene_annot").removeClass("hidden"):$(".gene_annot").addClass("hidden")}),$(document).on("change","#node_annot",function(){$(this).is(":checked")?$(".names").removeClass("hidden"):$(".names").addClass("hidden")}),$(document).on("mouseenter",".all",function(){var e="#"+$(this).attr("id")+"_arrow";$(this).css({opacity:1}),$(this).css({"z-index":10}),$(e).css({"z-index":10}),$(e).css({opacity:1})}),$(document).on("mouseenter",".line",function(){$(this).css({opacity:1}),$(this).css({"z-index":10})}),$(document).on("mouseleave",".all",function(){var e="#"+$(this).attr("id")+"_arrow",s=$("#opacity_slider").val();val=s/100,$(this).css({opacity:val}),$(this).css({"z-index":2}),$(e).css({"z-index":2}),$(e).css({opacity:val})}),$(document).on("mouseleave",".line",function(){var e=$("#line_opacity_slider").val();val=e/100,$(this).css({opacity:val}),$(this).css({"z-index":1})}),reveal(),$(document).on("click",".copy",function(){copy_to_clipboard($(this).attr("id")+"_gene_sequence"),$("#statusmessage").text("Sequence copied to clipboard!").animate({"margin-top":0},200),setTimeout(function(){$("#statusmessage").animate({"margin-top":-25},200)},3e3)}),$(document).on("click",".exit",function(){$(this).closest(".info_box").animate({right:"-350"},"slow").addClass("hidden")}),$(document).on("change","#opacity_slider",function(){var e=$(this).val();val=e/100,$(".left, .right, .left_m, .right_m, .left_s, .right_s, .all").css({opacity:val})}),$(document).on("change","#line_opacity_slider",function(){var e=$(this).val();val=e/100,$(".line").css({opacity:val})}),$(document).on("click","#small",function(){$(".size").removeClass("active"),$(this).addClass("active"),$(".R").css({marginRight:"6px"}),$(".L").css({marginLeft:"6px"}),$(".all").css({height:"6px"}),$(".line").css({height:"2px"}),$(".AR").removeClass("right"),$(".AR").removeClass("right_m"),$(".AR").removeClass("right_s"),$(".AR").addClass("right_s"),$(".AL").removeClass("left"),$(".AL").removeClass("left_m"),$(".AL").removeClass("left_s"),$(".AL").addClass("left_s"),$(".tree").addClass("hidden"),$(".s_tree").removeClass("hidden"),$(".names").addClass("hidden"),$("#node_annot").prop("checked",!1),$("#node_annot").prop("disabled",!0)}),$(document).on("click","#medium",function(){$(".size").removeClass("active"),$(this).addClass("active"),$(".R").css({marginRight:"8px"}),$(".L").css({marginLeft:"8px"}),$(".all").css({height:"16px"}),$(".line").css({height:"5px"}),$(".AL").removeClass("left"),$(".AL").removeClass("left_m"),$(".AL").removeClass("left_s"),$(".AL").addClass("left_m"),$(".AR").removeClass("right"),$(".AR").removeClass("right_m"),$(".AR").removeClass("right_s"),$(".AR").addClass("right_m"),$(".tree").addClass("hidden"),$(".m_tree").removeClass("hidden"),$(".names").css({marginTop:"3px"}),$(".names p").css({fontSize:"8px"}),$(".names p").css({height:"14px"}),$(".names").removeClass("hidden"),$("#node_annot").prop("checked",!0),$("#node_annot").prop("disabled",!1)}),$(document).on("click","#large",function(){$(".size").removeClass("active"),$(this).addClass("active"),$(".R").css({marginRight:"15px"}),$(".L").css({marginLeft:"15px"}),$(".all").css({height:"30px"}),$(".line").css({height:"5px"}),$(".AL").removeClass("left"),$(".AL").removeClass("left_m"),$(".AL").removeClass("left_s"),$(".AL").addClass("left"),$(".AR").removeClass("right"),$(".AR").removeClass("right_m"),$(".AR").removeClass("right_s"),$(".AR").addClass("right"),$(".tree").addClass("hidden"),$(".l_tree").removeClass("hidden"),$(".names p").css({fontSize:"12px"}),$(".names p").css({height:"28px"}),$(".names").css({marginTop:"7px"}),$(".names").removeClass("hidden"),$("#node_annot").prop("checked",!0),$("#node_annot").prop("disabled",!1)}),$(document).on("change","#scale",function(){var e=$(this).val();0==e?$(".sequences").css({width:"auto"}):$(".sequences").css({width:e+"px"})}),$(document).on("change","#node_width_slider",function(){var e=$(this).val();0==e?$(".names").css({width:"auto"}):$(".names").css({width:e+"px"})}),$(document).ready(function(){var e=$(window).width(),s=$("#node_width_slider").val();$(".names").hasClass("hidden")&&(s=0);var n=e-300-s-10;$(".sequences").css({width:n}),$("#scale").val(n),$("#node_annot").prop("checked",!0),$("#gene_annot").prop("checked",!0)}),$(document).on("click","#scale_reset",function(){var e=$(window).width(),s=$("#node_width_slider").val();$(".names").hasClass("hidden")&&(s=0);var n=e-300-s-10;$(".sequences").css({width:n}),$("#scale").val(n)}),$(document).on("click","#options",function(){$(".options").hasClass("open")?($(".options").css({top:-300}),$(".options").removeClass("open")):($(".options").css({top:50}),$(".options").addClass("open"))}),$(document).on("click","#create_image",function(){var e=$(".sequences").width()+$(".names").width()+$(".phylo").width()+30;html2canvas(document.querySelector("#capture"),{width:e}).then(function(e){console.log(e),saveAs(e.toDataURL(),"GEnView_Image.svg")})});
	</script>
	<body>
        <div class="navigation">
        <button id="small" class="size">Small</button>
        <button id="medium" class="size">Medium</button>
        <button id="large" class="size active">Large</button>
        <button id="options" class="size" style="float: right;">Options</button>
    </div>
    <div class="options">
        <label for="opacity_slider">Gene Opacity</label>
        <input id="opacity_slider" type="range" min="0" max="100" value="100">
        <label for="line_opacity_slider">Sequence Opacity</label>
        <input id="line_opacity_slider" type="range" min="0" max="100" value="20">
        <label for="scale">Environment Scale</label>
        <input id="scale" type="range" min="50" max="10000">
        <button id="scale_reset">Reset Scale</button>
		<label for="node_annot">Node Annotation</label>
        <input id="node_annot" type="checkbox" checked>
		<label for="node_width_slider">Node Width</label>
        <input id="node_width_slider" type="range" min="10" max="500" value="200">
		<label for="gene_annot">Gene Annotation</label>
		<input id="gene_annot" type="checkbox" checked>
		<button id="create_image">Create Image</button>
		<!--<button id="create_pdf">Create PDF</button>-->
    </div>
	<div class="space"></div>
	<div id="statusmessage"></div>
	<div id="capture" class="large_grid_container">
		<div class="phylo" style="grid-column-start: 1; grid-column-end: 2;grid-row-start: 1;grid-row-end: 2;">
			<div class="phylo_div">"""        
	html2 = """
			</div>  
		</div>
	"""
	html3 = """
	</div>
	</body>
	</html>"""

	output = html + tree_string + html2 + name_string + string + html3
	return output

def write_output(output, output_dir):
	if os.path.exists(output_dir +"/interactive_visualization.html"):
		os.remove(output_dir +"/interactive_visualization.html")
	f = open(output_dir +"/interactive_visualization.html", "a")
	f.write(output)
	f.write('\n')
	f.close()

	if args.log==True:
		if os.path.exists(output_dir +"/interactive_visualization.html") and os.path.getsize(output_dir +"/interactive_visualization.html")>0:
			log_lines.append('Interactive visualization created!\n')
		else:
			log_lines.append('Visualization could not be created...FAILED\n')
		write_log()

def check_color_format():
	
	cust_lines=[line.rstrip('\n') for line in open(args.custom_colors, 'r')]
	splits=cust_lines[0].split('\t')

	if not len(splits)==3:
		print('\n'+f'Incorrect file format detected for {args.custom_colors}.'+'\nFormat should be class\\tkeyword1,keyword2\\trgbcolor\\n\n')	
		sys.exit()

def main():

	global args
	args=parse_arguments()

	if args.log==True:
		global log_lines
		log_lines=[]

	if not args.custom_colors=='False':
		check_color_format()
	#Extract genes from db
	extract()
	#Extract file containing flanking regions
	context_file=[os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis'.rstrip('/')+'/'+file for file in os.listdir(os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis') if file.endswith('_contexts.fna')]

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

	###HTML OUTPUT
	#Create index phylogenetic tree f
	final_tree = create_tree_index(context_file[0].replace('.fna', '.unique.tree'))
	#Create html of phylogenetic tree and corresponding sequences
	output = create_html(final_tree, os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis'.rstrip('/')+'/'+'visualization_meta.csv')
	#Write HTML output file
	write_output(output, os.path.dirname(args.db).rstrip('/')+'/'+args.gene.lower()+'_'+str(args.id)+'_analysis')
	print('Visualization ready!')
	exit()

def write_log():
	
	with open(f'{os.path.abspath(args.target_directory)}/genview_viz.log', 'w') as outfile:
		for line in log_lines:
			outfile.write(line)

if __name__=='__main__':
	main()
