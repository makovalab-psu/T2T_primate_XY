#!/bin/bash
#SBATCH --job-name=make_windows.20230619
#SBATCH --partition=main
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=make_windows.202306169.%j.log
#SBATCH --time=168:00:00

set -e
set -x
pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

references=("chrY.human.WideHardmasked.fasta" "chrY.chimpanzee.WideHardmasked.fasta" "chrY.bonobo.WideHardmasked.fasta" "chrY.gorilla.WideHardmasked.fasta" "chrY.sorang.WideHardmasked.fasta" "chrY.borang.WideHardmasked.fasta" "chrY.gibbon.WideHardmasked.fasta" "chrX.human.WideHardmasked.fasta" "chrX.chimpanzee.WideHardmasked.fasta" "chrX.bonobo.WideHardmasked.fasta" "chrX.gorilla.WideHardmasked.fasta" "chrX.sorang.WideHardmasked.fasta" "chrX.borang.WideHardmasked.fasta" "chrX.gibbon.WideHardmasked.fasta")

threads=1

#CREATE THE 5kb windows
for reference in "${references[@]}"
do
	echo ${reference}

	makeblastdb -in ${reference} -parse_seqids -dbtype nucl
	#prepare bedtools genomes
	samtools faidx ${reference} #index the reference
	awk -v OFS='\t' {'print $1,$2'} ${reference}.fai > ${reference}.genome.txt #use the indexed information to create a genome file
	echo "Bedtools genomes generated."

	#create overlapping windows of the size 5k with 2kb steps in the bed format
	bedtools makewindows -g ${reference}.genome.txt -w 5000 -s 2000 -i srcwinnum >${reference}.windows.5000.bed
	echo "5kb sliding windows with 2kb steps generated."

	#generate windows of the size 5k using a version of the Y chromosome
	bedtools getfasta -fi ${reference} -bed ${reference}.windows.5000.bed >${reference}.windows.5000.fasta
done


echo "Done."
