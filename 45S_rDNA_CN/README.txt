
Author
Matthew Borchers
--------------------------
Brief Pipeline Description

This pipeline uses kmer depth in whole-genome sequencing data to estimate the
copy number of 45S rDNA. Regions with G/C content within 1 percent of the
45S reference are identified for normalization purposes. Only kmers which are
unique to these regions are used. Kmers from both regions are counted, and the
median kmer count of the 45S kmers is divided by the median count of the
normalization regoins (assumed to be single copy), yielding copy number.
--------------------------
Running the Pipeline

A config file needs to be created for snakemake with a list of IDs and kmer size
to use. The config for this project was:

=========
# config file for ribosome CN pipeline with GC normzliation method

ID:
  - "bonobo"
  - "gibbon"
  - "gorilla"
  - "s_orangutan"

K:
  - "31"
=========

Within the project directory, users will need to manually create directories
labelled genomes, rDNA_references, jellyfish_files, WGS, and results. Users should
place a fasta file for the reference genome in single line format and labeled
"{ID}_singleline.fa" in the genomes directory. The Snakefile also expects a
symlink of this file in the same directory labeled "{ID}.fa". A genome file should
be created using samtools faidx and placed in the same directory. A file with the
reference sequence for just the 18S portion of the rDNA needs to be placed in the
directory rDNA_references.

To initiate the pipeline, users enter "snakemake count_cn" in the linux console.

A more streamlined, generalized, and feature rich version of this pipeline will
eventually be published, currently under the name "COpy Numbers using Kmers
NOrmalized to Depth of comparable sequence" (CONKORD).

--------------------------
Tool Versions

jellyfish 2.2.10
samtools 1.3.1
htslib 1.3.1
bedtools v2.29.0

--------------------------
