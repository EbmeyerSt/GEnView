# EnView
A phylogeny based comparative genomics software to analyze the genetic environment of genes. The user can select one or several taxa and provide one or several reference gene(s). Genomes and plasmids (based on user choice) will be downloaded from the NCBI Assembly/NR database and searched for the respective gene. Up to 20kb of the genes genetic environment are extracted, annotated and aligned between all genomes. The sequences are then visualized, enabling comparison of synteny and gene content.

# Installation

EnView consists of three main scripts - enview_db.py, extract.py and visualize.py

The code can be downloaded simply through cloning the git repository:

`git clone https://github.com/EbmeyerSt/EnView.git`

# Dependencies
* Python >=3.7.6
  * ete3 >=3.1.1
  * pandas >=1.0.3
  * sqlite3 >=2.6
  
* Diamond >=0.9.24
* Prodigal >=2.6.3
* CD-hit >=4.7
* Usearch >=8.0.1445
* mafft >=7.310
* FastTree >=2.1.9
(optional)
* integron_finder >= 2


# Usage

# Output files
