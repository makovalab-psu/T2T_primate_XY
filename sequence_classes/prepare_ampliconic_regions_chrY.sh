set -e
set -x

#prepare the regions with intrachromosomal similarity in the bed format
#this script reports the identity of the window anywhere on the Y, other than the position it originates from, provided it exists
for a in blastn*WideHardmasked.fasta.txt; do echo $a; Rscript createBedFileFromBlastn.R $a; done;

#merge the neigboring hits, and only keep the consecutive regions larger than 90000
for a in blastn*WideHardmasked.fasta.bed; do cat $a | sort -V -k1,1 -k2,2 | bedtools merge -d 100000 | awk '{ if (($3-$2)>90000) print;}' >${a}.ampliconic.bed; done;

#merge into a single file for circos
cat blastn.long.chrY.human*ampliconic.bed | sed -e 's/chrY/hs1/' >ampliconic.bed.tmp
cat blastn.long.chrY.chimpanzee*ampliconic.bed | sed -e 's/chrY/hs2/' >>ampliconic.bed.tmp
cat blastn.long.chrY.bonobo*ampliconic.bed | sed -e 's/chrY/hs3/' >>ampliconic.bed.tmp
cat blastn.long.chrY.gorilla*ampliconic.bed | sed -e 's/chrY/hs4/' >>ampliconic.bed.tmp
cat blastn.long.chrY.sorang*ampliconic.bed | sed -e 's/chrY/hs5/' >>ampliconic.bed.tmp
cat blastn.long.chrY.borang*ampliconic.bed | sed -e 's/chrY/hs6/' >>ampliconic.bed.tmp
cat blastn.long.chrY.gibbon*ampliconic.bed | sed -e 's/chrY/hs7/' >>ampliconic.bed.tmp

cat ampliconic.bed.tmp | awk '{print $0 "\tfill_color=blue"}' | sed 's/ /\t/g' | sort -V -k1,1 -k2,2 >ampliconic.bed
rm ampliconic.bed.tmp

