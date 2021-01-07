# EnView
A phylogeny based comparative genomics software to analyze the genetic environment of genes. The user can select one or several taxa and provide one or several reference protein(s). Genomes and plasmids (based on user choice) will be downloaded from the NCBI Assembly/NR database and searched for the respective gene. Up to 20kb of the genes genetic environment are extracted, annotated and aligned between all genomes. The sequences are then visualized, enabling comparison of synteny and gene content.

# Installation

EnView consists of three main scripts - enview_db.py, extract.py and visualize.py

The code can be downloaded simply through cloning the git repository:

`git clone https://github.com/EbmeyerSt/EnView.git`

# Dependencies

We recommend to set up a virtual environment and install all dependencies to run EnView there.

* Python >=3.7.6
  * ete3 >=3.1.1
  * pandas >=1.0.3
  * sqlite3 >=2.6
  
The following softwares should be located in your $PATH:
  
* Diamond >=0.9.24
* Prodigal >=2.6.3
* CD-hit >=4.7
* Usearch >=8.0.1445
* mafft >=7.310
* FastTree >=2.1.9
(optional)
* integron_finder >= 2


# Usage

`enview_db.py` automatically downloads genomes (and if specied, plasmids) of the specified species from the NCBI Assembly/NR database and search them for the provided reference protein(s). Up to 20kb are extracted upstream and downstream, annotated and aligned.

```
usage: enview_db.py [-h] -d TARGET_DIRECTORY -db DATABASE
                                  [-p PROCESSES] [-id IDENTITY] [--erase]
                                  [-scov SUBJECT_COVERAGE] -split
                                  SPLIT [--update] [--is_db IS_DB]
                                  [--taxa TAXA [TAXA ...]] [--assemblies]
                                  [--plasmids] [--integron_finder]


optional arguments:
  -h, --help            show this help message and exit
  -d TARGET_DIRECTORY, --target_directory TARGET_DIRECTORY
                        path to folder containing a folder for each assembly to be processed
  -db DATABASE, --database DATABASE
                        fasta/multifasta file containing genes to be annotated
  -p PROCESSES, --processes PROCESSES
                        of cores to run the script on
  -id IDENTITY, --identity IDENTITY
                        identity cutoff for hits to be saved to the database
  --erase               erase results of previous analysis and create new ones
  -scov SUBJECT_COVERAGE, --subject_coverage SUBJECT_COVERAGE
                        minimum coverage for a hit to be saved to db
  -split SPLIT          number of files to obtain for processing flanking regions (increasing annotation speed)
  --update              downloads new genomes and updates database
  --is_db IS_DB         database containing IS, integrons, ISCR sequences
  --taxa TAXA [TAXA ...]
                        taxon/taxa names to download genomes for - use "all" do download all available genomes
  --assemblies          Search NCBI Assembly database 
  --plasmids            Search NCBI Refseq plasmid database
  --integron_finder     Use integron_finderv2 for integron identification
  ```

**WARNING**: if you specify --taxa 'all', enview_db.py will attempt to download all available genomes and plasmids - > 4TB of data! Processing of all data at once will take several **days**. Consider downloading smaller subsets instead. If you do not know in which taxa your gene is present, we recommend doing a manual blast at https://blast.ncbi.nlm.nih.gov/Blast.cgi first - This will show you in which taxa your protein is found!

`extract.py` extracts the gene coding for the reference protein and its genetic environment as multifasta from the previously created database.

```
usage: extract_genes_sqlite_v7.5.py [-h] [-genes GENES [GENES ...]] -o O -db
                                    DB -id ID
                                    [-taxa TAXA [TAXA ...]]

________________________________________________________________________________

Extract specified genes and flanking regions from db	
________________________________________________________________________________

optional arguments:
  -h, --help            show this help message and exit
  -genes GENES [GENES ...]
                        genes to extract. If list: "qnrA" "qnrB" ...
  -o O                  path to output directory
  -db DB                sqlite db file
  -id ID                percent identity threshold for genes to extract
  -taxa TAXA [TAXA ...]
                        list of genera and/or species to extract.
                        	By default all taxa are extracted
```

`visualize.py` takes the previously extracted sequences as input and creates a pdf visualizing the reference gene in its genetic environment

```
usage: visualize_v8_py3.py [-h] -db DB -o O [--force] [--compressed] [--all]

________________________________________________________________________________

Visualize annotate genes and genetic environments
________________________________________________________________________________

optional arguments:
  -h, --help    show this help message and exit
  -db DB        sqlite3 db containing annotations
  -o O          target directory
  --force       Force new alignment and phylogeny
  --compressed  Compress number of displayed sequences
  --all         Create visualizations for all gene analyses in directory
```
 specifying `--compressed` will cluster all sequences at 95% nucleotide identity and remove duplicates. Only the centroids of the resulting clusters will be visualized. 

# Output files

The following output files will be produced in the specified output directory when running `visualize.py`:

**annotation_meta.csv** - Contains information on annotated genes in the extracted range of the target gene, such as name, position and sequence
**yourgenename_contexts._tree_annotated.pdf** - PDF containing a phylogeny based visualization of the target genes genetic environment
**yourgenename_contexts.fna** - FASTA file containing extracted sequence (target gene+genetic environment) for every genome
**yourgenename_contexts.unique.fna** - FASTA file containing extracted sequence (target gene+genetic environment) for every genome after removing duplicates
**yourgenename_contexts.unique.aln** Alignment of unique extracted sequences
**yourgenename_contexts.unique.tree** Tree file created by FastTree2
