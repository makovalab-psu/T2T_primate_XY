set -e
set -x

#prepare the regions with intrachromosomal similarity in the bed format
#this script reports the identity of the window anywhere on the Y, other than the position it originates from, provided it exists
for a in blastn*WideHardmasked.fasta.txt; do echo $a; Rscript createBedFileFromBlastn.R $a; done;

#merge the neigboring hits, and only keep the consecutive regions larger than 90000
for a in blastn*WideHardmasked.fasta.bed; do cat $a | sort -V -k1,1 -k2,2 | bedtools merge -d 100000 | awk '{ if (($3-$2)>90000) print;}' >${a}.ampliconic.bed; done;

#merge into a single file for circos
cat blastn.long.chrX.human*ampliconic.bed | sed -e 's/chrX/hs1/' >ampliconic.bed.tmp
cat blastn.long.chrX.chimpanzee*ampliconic.bed | sed -e 's/chrX/hs2/' >>ampliconic.bed.tmp
cat blastn.long.chrX.bonobo*ampliconic.bed | sed -e 's/chrX/hs3/' >>ampliconic.bed.tmp
cat blastn.long.chrX.gorilla*ampliconic.bed | sed -e 's/chrX/hs4/' >>ampliconic.bed.tmp
cat blastn.long.chrX.sorang*ampliconic.bed | sed -e 's/chrX/hs5/' >>ampliconic.bed.tmp
cat blastn.long.chrX.borang*ampliconic.bed | sed -e 's/chrX/hs6/' >>ampliconic.bed.tmp
cat blastn.long.chrX.gibbon*ampliconic.bed | sed -e 's/chrX/hs7/' >>ampliconic.bed.tmp

cat ampliconic.bed.tmp | awk '{print $0 "\tfill_color=blue"}' | sed 's/ /\t/g' | sort -V -k1,1 -k2,2 >ampliconic.bed
rm ampliconic.bed.tmp

