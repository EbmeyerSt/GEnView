# GEnView
A phylogeny based comparative genomics software to analyze the genetic environment of genes. The user can select one or several taxa and provide one or several reference protein(s). Genomes and plasmids (based on user choice) will be downloaded from the NCBI Assembly/NR database and searched for the respective gene. Alternatively, custom genomes can be provided. User selected stretches (20kbp by default) of the genes genetic environment are extracted, annotated and aligned between all genomes. The sequences are then visualized, enabling comparison of synteny and gene content.

# Current Bug:
We identified a bug in GEnView that throws an error in the visualization command if you had more than one gene in your reference database. A patch is implemented in the newest release (v.0.1.1), until that is updated on bioconda you can circumvent the issue by just searching your target genomes for one reference gene. 

# Installation

GEnView can be downloaded using conda:
`conda install -c bioconda genview`

Alternatively, the code can be downloaded simply through cloning the git repository:

`git clone https://github.com/EbmeyerSt/GEnView.git`

# Dependencies

We recommend to set up a virtual environment and install all dependencies to run GEnView there.

* Python >=3.6
  * pandas >=1.0.3
  * sqlite3 >=2.6
  * Biopython >=1.68
  
The following softwares should be located in your $PATH:
  
* Diamond >=0.9.24
* Prodigal >=2.6.3
* CD-hit >=4.7
* mafft >=7.310
* FastTree >=2.1.9



# Usage

`genview-makedb` (python /path/to/genview/genview_scripts/genview_create_db.py if downloaded manually) automatically downloads genomes (and if specied, plasmids) of the specified species from the NCBI Assembly/NR database and search them for the provided reference protein(s). Up to 20kb are extracted upstream and downstream, annotated and aligned.

An example walk-through can be found in the wiki: https://github.com/EbmeyerSt/GEnView/wiki

```
usage: genview-makedb [-h] -d TARGET_DIRECTORY -db DATABASE
                            [-p PROCESSES] [-id IDENTITY]
                            [-scov SUBJECT_COVERAGE] [--is_db] [--uniprot_db]
                            [--uniprot_cutoff UNIPROT_CUTOFF]
                            [--taxa TAXA [TAXA ...]] [--assemblies]
                            [--plasmids] [--custom] [--save_tmps]
                            [--acc_list ACC_LIST]
                            [--flanking_length FLANKING_LENGTH]

Creates sqlite3 database with genetic environment from genomes containing the provided reference gene(s).

optional arguments:
  -h, --help            show this help message and exit
  -d TARGET_DIRECTORY, --target_directory TARGET_DIRECTORY
                        path to output directory
  -db DATABASE, --database DATABASE
                        fasta/multifasta file containing amino acid sequences of translated genes to be annotated
  -p PROCESSES, --processes PROCESSES
                        number of cores to run the script on
  -id IDENTITY, --identity IDENTITY
                        identity cutoff for hits to be saved to the database (e.g 80 for 80% cutoff)
  -scov SUBJECT_COVERAGE, --subject_coverage SUBJECT_COVERAGE
                        minimum coverage for a hit to be saved to db (e.g 80 for 80% cutoff)
  --is_db               database containing IS, integrons, ISCR sequences
  --uniprot_db          path to database for annotation of surrounding sequences. If unspecified default uniprotKB database will be downloaded to target directory
  --uniprot_cutoff UNIPROT_CUTOFF
                        % identity threshold for annotating orfs aurrounding the target sequence, default 60
  --taxa TAXA [TAXA ...]
                        taxon/taxa names to download genomes for - use "all" do download all available genomes, cannot be specified at the same time as --acc_list
  --assemblies          Search NCBI Assembly database
  --plasmids            Search NCBI Refseq plasmid database
  --custom              Search custom genomes
  --save_tmps           keep temporary files
  --acc_list ACC_LIST   csv file containing one genome accession number per row, cannot be specied at the same time as --taxa
  --flanking_length FLANKING_LENGTH
                        Max length of flanking regions to annotate

  ```

**WARNING**: if you specify --taxa 'all', enview_db.py will attempt to download all available genomes and plasmids - > 4TB of data! Processing of all data at once will take several **days**. Consider downloading smaller subsets instead. If you do not know in which taxa your gene is present, we recommend doing a manual blast at https://blast.ncbi.nlm.nih.gov/Blast.cgi first - This will show you in which taxa your protein is found!


`genview-visualize` (python /path/to/genview/genview_scripts/genview_visualize.py if downloaded manually) takes the previously extracted sequences as input and creates an interactive visualization of the reference gene in its different genetic environments.

```
usage: genview-visualize [-h] -gene GENE -db DB -id ID
                            [-taxa TAXA [TAXA ...]] [--force] [--compressed]

Extract, visualize and annotate genes and genetic environments from genview database

optional arguments:
  -h, --help            show this help message and exit
  -gene GENE            name of gene/orf to extract and visualize
  -db DB                genview database created by genview-create-db
  -id ID                percent identity threshold for genes to extract
  -taxa TAXA [TAXA ...]
                        list of genera and/or species to extract
                        By default all taxa are extracted
  --force               Force new alignment and phylogeny
  --compressed          Compress number of displayed sequences, helpful with large number of identical sequences

```
 specifying `--compressed` will cluster all sequences based on their annotaton profile and remove duplicates. Only the centroids of the resulting clusters will be visualized in that case. Helpful for large numbers of genomes. 

 Note that the output file was tested on Google Chrome (Version 91.0.4472.164) and Firefox (Version 90.0) and may not function properly when opened in older versions other web browsers.


## Custom genomes or contigs

It is possible to use GEnView with custom genomes and/or contigs. To do so, the parameter `--custom` must be specified and the genomes/contigs to be used must be placed in the target directory in a file named **all_assemblies.fna**.
If you have taxonomical information of your sequences this can be provided in a file named **custom_summary.txt** which then also must be placed in the target directory. The format of the file should be tab separated with the columns: contigID, taxonomical ID, taxonomical name, custom information (e.g. country, site, time...).


# Output files

The following output files will be produced in the specified output directory when running `genview-visualize`:

**annotation_meta.csv** - Contains information on annotated genes in the extracted range of the target gene, such as name, position and sequence

**visualization_meta.csv** - Contains information on annotated genes in the extracted range of the target gene, such as name, position and sequence in a format used for creation of the interactive HTML visualization output file.

**yourgenename_contexts._tree_annotated.pdf** - PDF containing a phylogeny based visualization of the target genes genetic environment

**yourgenename_contexts.fna** - FASTA file containing extracted sequence (target gene+genetic environment) for every genome

**yourgenename_contexts.unique.fna** - FASTA file containing extracted sequence (target gene+genetic environment) for every genome after removing duplicates

**yourgenename_contexts.unique.aln** Alignment of unique extracted sequences

**yourgenename_contexts.unique.tree** Tree file created by FastTree2

**interactive_output.html** Interactive visualization of tree and extracted sequence (target gene+genetic environment). Viewed in web browser (Google Chrome 91.0.4472.164 and Firefox 90.0).
