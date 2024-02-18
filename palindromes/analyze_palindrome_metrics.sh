#!/bin/bash
set -e
set -x

###STEP 1
#load conda environment that contains bedtools installation
#for example the one described in https://github.com/biomonika/HPP/blob/main/palindromes/environment.yml
source /opt/miniconda/etc/profile.d/conda.sh
conda activate alignment #uses bedtools v2.30.0

###STEP 2
#using coordinates of palindromes from palindrover, and the reference Y chromosomes, extract the sequence of palindrome arms
myArray=("mPanTro3" "mPanPan1" "mGorGor1" "mPonAbe1" "mPonPyg2" "mSymSyn1")
for species in ${myArray[@]}; do bedtools getfasta -name -fi chrY.${species}.dip.20221111.fasta -bed ${species}.chrY.palindromes.bed >${species}.palindromeSeq.chrY.fa; done;
for species in ${myArray[@]}; do bedtools getfasta -name -fi chrX.${species}.dip.20221111.fasta -bed ${species}.chrX.palindromes.bed >${species}.palindromeSeq.chrX.fa; done;

###STEP 2
#download perl script to split palindromes into individual files using downloaded tool fasta-splitter
#each palindrome arm will be an individual sequence
#downloaded https://kirill-kryukov.com/study/tools/fasta-splitter/
mkdir -p parts
for a in *palindromeSeq.chr*.fa; do echo $a; perl fasta-splitter.pl --part-size 1 --eol unix --line-length 0 --nopad --out-dir parts --measure count ${a}; done;
#one sequence per file, then they need to be analyzed in pairs to capture palindrome arms correctly
#note that each palindrome will get a new name here, using the numbering of the files

###STEP 3
#align palindrome arms with each other 
#NOTE THAT THE STRETCHER STEP WILL RUN IN PARALLEL AND MIGHT USE LOTS OF RESOURCES
#because the corresponding palindrome arms have corresponding numbers (e.g. part-1 and part-2 are two arms of the same palindrome), we can do this in a loop

cd parts

for species in ${myArray[@]}; do echo ${species}; for i in $(seq 1 2 `ls ${species}.palindromeSeq.chrY*.fa | wc -l`); do echo $i; seq1=${species}.palindromeSeq.chrY.part-${i}.fa; seq2=${species}.palindromeSeq.chrY.part-$((i+1)).fa; echo $seq1; echo $seq2; echo "===="; stretcher -asequence $seq1 -sreverse2 -bsequence $seq2 -outfile ${seq1}.stretcher & done; done; wait;
for species in ${myArray[@]}; do echo ${species}; for i in $(seq 1 2 `ls ${species}.palindromeSeq.chrX*.fa | wc -l`); do echo $i; seq1=${species}.palindromeSeq.chrX.part-${i}.fa; seq2=${species}.palindromeSeq.chrX.part-$((i+1)).fa; echo $seq1; echo $seq2; echo "===="; stretcher -asequence $seq1 -sreverse2 -bsequence $seq2 -outfile ${seq1}.stretcher & done; done; wait;
#create key value pairs, in other words remember how the name from palindrover and new numerical name correspond
for a in *.fa; do echo $a | tr "\n" "\t"; head -n 1 $a | sed s'/>//g'; done >../palindromes.keyvaluepairs.txt

###STEP 4
#parse the alignment results from stretcher, and write down the calculated identity between palindrome arms
#build table with results
#2 columns, filename and identity
for a in *stretcher; do echo $a | tr "\n" "\t" |  sed s'/.stretcher//g'; grep "Identity" $a | tr -s '[:space:]' ' ' | cut -d' ' -f4 | sed s'/)//g' | sed s'/(//g' | sed s'/%//g'; done >../palindroms.identity.txt

###STEP 5
#calculate the GC content of palindrome arms, and write it down
for a in *.fa; do echo $a; geecee -sequence ${a} -outfile ${a}.geecee; done;
for a in *geecee; do echo $a | tr "\n" "\t" |  sed s'/.geecee//g'; grep -v "#" $a | cut -d' ' -f3; done >../palindroms.geecee.txt

cd ..

###STEP 5
#extract the coordinates of the spacers and write down the sequence of the spacers
for a in *palindromes.bed; do echo $a; python extract_spacer.py $a; done;
mkdir -p spacers
for species in ${myArray[@]}; do bedtools getfasta -name -fi chrY.${species}.dip.20221111.fasta -bed ${species}.chrY.palindromes.bedspacers.bed >${species}.spacerSeq.chrY.fa; done;
for species in ${myArray[@]}; do bedtools getfasta -name -fi chrX.${species}.dip.20221111.fasta -bed ${species}.chrX.palindromes.bedspacers.bed >${species}.spacerSeq.chrX.fa; done;


###STEP 6
#split the palindrome spacers into individual files
for a in *spacerSeq.chr*.fa; do echo $a; perl fasta-splitter.pl --part-size 1 --eol unix --line-length 0 --nopad --out-dir spacers --measure count ${a}; done;


###STEP 7
#calculate the GC content of spacers, and write it down
cd spacers
for a in *.fa; do echo $a; geecee -sequence ${a} -outfile ${a}.geecee; done;
for a in *geecee; do echo $a | tr "\n" "\t" |  sed s'/.geecee//g'; grep -v "#" $a | cut -d' ' -f3; done >../spacers.geecee.txt
#create key value pairs
for a in *.fa; do echo $a | tr "\n" "\t"; head -n 1 $a | sed s'/>//g'; done >../spacers.keyvaluepairs.txt

###STEP 8
#calculate spacer lengths
cd ..
for a in *.bedspacers.bed; do echo $a; awk '{print $4"\t"($3-$2)}' $a >>palindromes.spacerlengths.txt; done;
#calculate arm lengths, but only for the first arm
for a in *.palindromes.bed; do echo $a; awk '{print $4"\t"($3-$2)}' $a | egrep "A\b" >>palindromes.armlengths.txt; done;

###STEP 9
#now the data is ready to be analyzed and plotted in R
#use script palindrome_metrics.Rnw



