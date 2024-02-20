De novo gene analysis
Margaux Aubel, University of Münster, m.aubel@uni-muenster.de

The scripts are for analysis of de novo genes after detection of appropriate candidates. 
TE_denovo.py compares de novo gene candidates with annotated repeats. Inputs needed are 1. a gff file including only the (de novo) genes of interest and 2. the file including annotated repeats in repeatmasker format.
get_all_dn_xy.py extracts aligned regions of de novo genes detected on the Y chromosome in the X chromosome. Input files needed are the .maf files of X and Y chromosome. The exon coordinates of candidate genes are provided as lists in the python script.
