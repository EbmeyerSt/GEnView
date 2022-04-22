# DUE TO A COMING MAJOR UPDATE, THE WIKI IS CURRENTLY UNDER RECONSTRUCTION

# GEnView
A phylogeny based comparative genomics software to analyze the genetic environment of genes. The user can select one or several taxa and provide one or several reference protein(s). Genomes and plasmids (based on user choice) will be downloaded from the NCBI Assembly/NR database and searched for the respective gene. Alternatively, custom genomes can be provided. User selected stretches (20kbp by default) of the genes genetic environment are extracted, annotated and aligned between all genomes. The sequences are then visualized, enabling comparison of synteny and gene content.


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

## genview-makedb

`genview-makedb` (python /path/to/genview/genview_scripts/genview_create_db.py if downloaded manually) automatically downloads genomes (and if specied, plasmids) of the specified species from the NCBI Assembly/NR database and search them for the provided reference protein(s). Up to 20kb are extracted upstream and downstream, annotated and aligned. Detailed parameter explanations can be found below.

An example walk-through can be found in the wiki: https://github.com/EbmeyerSt/GEnView/wiki

```
usage: genview_create_db.py [-h] -d TARGET_DIRECTORY -db DATABASE [-p PROCESSES] [-id IDENTITY]
                            [-scov SUBJECT_COVERAGE] [--update] [--uniprot_db UNIPROT_DB]
                            [--uniprot_cutoff UNIPROT_CUTOFF] [--taxa TAXA [TAXA ...]] [--assemblies] [--plasmids]
                            [--local LOCAL] [--save_tmps] [--accessions ACCESSIONS] [--flanking_length FLANKING_LENGTH]
                            [--kraken2 KRAKEN2] [--log] [--clean]

Creates sqlite3 database with genetic environment from genomes containing the provided reference gene(s).

optional arguments:
  -h, --help            show this help message and exit
  -d TARGET_DIRECTORY, --target_directory TARGET_DIRECTORY
                        path to output directory
  -db DATABASE, --database DATABASE
                        fasta/multifasta file containing amino acid sequences of translated target genes to be annotated
  -p PROCESSES, --processes PROCESSES
                        number of cores to run the script on
  -id IDENTITY, --identity IDENTITY
                        identity cutoff for hits to be saved to the database (e.g 80 for 80% cutoff)
  -scov SUBJECT_COVERAGE, --subject_coverage SUBJECT_COVERAGE
                        minimum coverage for a hit to be saved to db (e.g 80 for 80% cutoff)
  --update              update an existing genview database with new genomes
  --uniprot_db UNIPROT_DB
                        Path to uniprotKB database
  --uniprot_cutoff UNIPROT_CUTOFF
                        % identity threshold for annotating orfs aurrounding the target sequence, default 60
  --taxa TAXA [TAXA ...]
                        taxon/taxa names to download genomes for - use "all" do download all available genomes, cannot be specified at the same time as --accessions
  --assemblies          Search NCBI Assembly database
  --plasmids            Search NCBI Refseq plasmid database
  --local LOCAL         path to local genomes
  --save_tmps           keep temporary files
  --accessions ACCESSIONS
                        csv file containing one genome accession number per row, cannot be specied at the same time as --taxa
  --flanking_length FLANKING_LENGTH
                        Max length of flanking regions to annotate
  --kraken2 KRAKEN2     Path to kraken2 database. Uses kraken2 to classify unclassified sequences.
  --log                 Write log file for debugging
  --clean               Erase files from previous genview runs from target directory

  ```

**WARNING**: if you specify --taxa 'all', genview-makedb will attempt to download all available genomes and plasmids - > 4TB of data! Processing of all data at once will take several **days** on a large server. Consider downloading smaller subsets instead. If you do not know in which taxa your gene is present, we recommend doing a manual blast at https://blast.ncbi.nlm.nih.gov/Blast.cgi first - This will show you in which taxa your reference sequence is found!


**genview-makedb parameter details**

```--update``` | Update a previously created database. When using ```--update```, specify the path to a directory containing the genview database you want to update as the target directory (-d) for this run

```--uniprot_db``` | Path to uniprotKB diamond database for annotation of the target genes flanking sequences. GenView will automatically attempt to download the database if ```--uniprot_db``` is not specified. If the download fails, instructions to manually download the database are provided and GEnView will exit.

```--taxa``` | taxon names to download and search genomes for, e.g ```--taxa 'Aeromonas'``` to download all genomes and/or plasmids of the genus *Aeromonas*, e.g ```--taxa 'Aeromonas caviae' 'Escherichia coli' 'Leclercia adecarboxylata'``` to download all genomes for *A. caviae, E. coli*, and *L. adecarboxylata*.

```--accessions```/ ```--plasmids``` | csv file containing one genome accession number per row, not to be specified at the same time as ```--taxa```. Accessions should either be in Assembly database format, e.g 'GCA_006364295.1' or nucleotide database format for plasmids, e.g 'NZ_01234.1'. If the providing Assembly accessions/plasmid accessions, specify ```--assemblies``` and/or ```--plasmids``` along with ```--accessions```.

```--local``` | Path to local genome files in fasta format. Headers should be as simple as possible, with or without species classification, e.g: '>contig1 Escherichia coli'

```--kraken2``` | Path to kraken2 database, can be used to classify unclassified local sequences. Requires prebuilt kraken2 database at specified path.

```--clean``` | Removes files from previous genview run from target directory.

## genview-visualize

`genview-visualize` (python /path/to/genview/genview_scripts/genview_visualize.py if downloaded manually) takes the previously extracted sequences as input and creates an interactive visualization of the reference gene in its different genetic environments.

```
usage: genview_visualize.py [-h] -gene GENE -db DB -id ID [-nodes NODES] [-taxa TAXA [TAXA ...]] [--force]
                            [--compressed] [--custom_colors CUSTOM_COLORS] [--log]

Extract, visualize and annotate genes and genetic environments from genview database

optional arguments:
  -h, --help            show this help message and exit
  -gene GENE            name of gene/orf to extract and visualize
  -db DB                genview database created by genview-create-db
  -id ID                percent identity threshold for genes to extract
  -nodes NODES          should nodes be connected to genome with solid line (solid), connected by dashed line (dash) or no connection (none)
  -taxa TAXA [TAXA ...]
                        list of genera and/or species to extract
                        By default all taxa are extracted
  --force               Force new alignment and phylogeny
  --compressed          Compress number of displayed sequences, helpful with large number of identical sequences
  --custom_colors CUSTOM_COLORS
                        path to file containing RGB color codes for gene color customization
  --log                 Write log file

```
 specifying `--compressed` will cluster all sequences based on their annotaton profile and remove duplicates. Only the centroids of the resulting clusters will be visualized in that case. Helpful for large numbers of genomes. 

 Note that the output file was tested on Google Chrome (Version 91.0.4472.164) and Firefox (Version 90.0) and may not function properly when opened in older versions other web browsers.
 
 **genview-visualize parameter details**
 
 ```-gene``` | gene names to extract and visualize from previously created database. E.g ```-gene 'PER-1'``` will extract and visualize all found PER-1 genes. ```-gene 'PER-'``` will extract and visualize all target genes containing 'PER-' from the database. ```-gene 'PER-1' 'PER-2'``` will extract and visualize all instances of 'PER-1' and 'PER-2' genes.
 
 ```-taxa``` | taxon names to extract and visualize from previously created database. E.g ```-taxa 'Aeromonas'``` to visualize th eidentified target gene in all sequences of the genus *Aeromonas*, e.g ```-taxa 'Aeromonas caviae' 'Escherichia coli' 'Leclercia adecarboxylata'``` to visualize the target gene for *A. caviae, E. coli*, and *L. adecarboxylata*.
 
 ```-db``` | path to genview database 'genview_database.db', created in the target directory for ```genview-makedb```
 
 ```--force``` | Force new alignment and phylogeny for a previously extracted gene
 
 ```--compressed``` | Cluster all sequences containing the target gene based on their annotation profile and remove dupliocates. Only the centroids of the resulting clusters will be visualized. helpful for large numbers of genomes.
 
 ```--custom_colors``` | Path to file specifying RGB color codes for user selected classes of genes. File should be text file in format:
 
 gene_class1\tkeyword1,keyword2\trgbcolorcode
 gene_class2\tkeyword3, keyword4, keyword3\trgbcolorcode
 gene_class3\tkeyword5\trgbcolorcode
 
 The 'keywords' are words found in the name of the respective genes. E.g if you want to display all efflux pumps in your visualization in blue, the corresponding row in the file specified under ```--custom_colors``` would look like this:
 
 pumps\tefflux, pump\trgb(70,130,180)
 

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
