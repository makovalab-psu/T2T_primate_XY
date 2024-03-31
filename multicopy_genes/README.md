# Multi-copy gene analysis

This folder contains a jupyter notebook used to identify and cluster homologous genes within and between species. 

## Input Data

 - The input of the method are the following genomic annotations provided by NCBI:
    https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029281585.1/ \
    https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/ \
    https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029289425.1/ \
    https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_028858775.1/ \
    https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_028885655.1/ \
    https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_028885625.1/ \
    https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_028878055.1/

 - Sequence classes annotation: \
    Additional Data file 2 of the publication.

 - Palindrome locations file
    Provided in the `/data` folder alongside the jupyter notebook. 


## Output 
 - fasta files - one per cluster of homology
 - tables of gene counts per species
 - tables of gene counts per sequence class
 - tables of pairwise identities of gene paralogs (within one specie)

Output tables are presented in supplementary tables S36 and S37.

## Requirements

 - All python libraries installed in venv used to run this notebook are listed in `requirements.txt`.

 - Local installation of `blast  2.15.0+`.

<hr/>
Karol Pál, karolpal [dot] jr [at] gmail [dot] com
