#!/bin/bash
set -e
set -x

#At the start, we have the following files:

#merged.chrX.repeats.species.bed => merged coordinates of the continous stretches of repeat sequences
#ampliconic.bed => merged coordinates of the regions with high intrachromosomal similarity
#PARs.bed => coordinates of pseudoautozomal regions from Bob Harris
#palindrover.merged.bed => coordinates of merged palindromes from Bob Harris

#Start with the regions annotated as repeats
#Subtract PARs _from the repeats, so that the repeats are shorter than PARs

bedtools subtract -a merged.chrX.repeats.species.bed -b PARs.bed | sort -V -k1,1 -k2,2 >repeats_without_PARs.bed 
#Now we have PARs and repeats

#Subtracts PARs and repeats _from ampliconic sequences, so that ampliconic sequences are shorter now
cat repeats_without_PARs.bed PARs.bed | cut -f1-3 | sort -V -k1,1 -k2,2 >concatenation.repeats_and_PARs.bed

#combine manual ampliconic regions with palindrome annotations
cat ampliconic.bed palindrover.merged.bed | cut -f1-3 | sort -V -k1,1 -k2,2 | bedtools merge >concatenation.ampliconic.bed

bedtools subtract -a concatenation.ampliconic.bed -b concatenation.repeats_and_PARs.bed >ampliconic_subtracted.bed 
#Now we have PARs and repeats and ampliconic sequences


#Find out the complement, call it OTHER
cat repeats_without_PARs.bed PARs.bed ampliconic_subtracted.bed | cut -f1-3 | sort -V -k1,1 -k2,2 >concatenation.repeats_and_PARs_and_ampliconic.bed
bedtools complement -i concatenation.repeats_and_PARs_and_ampliconic.bed -g Xgenomes.txt | awk '{print $0 "\tfill_color=gray"}' >OTHER.bed


#remove unnecessary files
rm -f concatenation.repeats_and_PARs.bed concatenation.repeats_and_PARs_and_ampliconic.bed 

echo "Print potential overlaps to make sure the separation into classes worked:"
bedtools multiinter -header -i PARs.bed ampliconic_subtracted.bed OTHER.bed repeats_without_PARs.bed  | awk '{if ($4>1) print}'

#assign the correct colors for each sequence class
awk 'BEGIN {FS=OFS="\t"} { $4 = "fill_color=#97cb99"; print }' PARs.bed >PARs.txt
awk 'BEGIN {FS=OFS="\t"} { $4 = "fill_color=#89c0e8"; print }' ampliconic_subtracted.bed >ampliconic_subtracted.txt
awk 'BEGIN {FS=OFS="\t"} { $4 = "fill_color=#d9d8d8"; print }' OTHER.bed >OTHER.txt
awk 'BEGIN {FS=OFS="\t"} { $4 = "fill_color=#6a4a1d"; print }' repeats_without_PARs.bed >repeats_subtracted.txt

#Finally, plot the following files
#PARs.txt => coordinates of pseudoautozomal regions from Bob Harris
#ampliconic_subtracted.txt => new coordinates of ampliconic sequences
#repeats_subtracted.txt => new coordinates of repeat sequences
#OTHER.txt => new coordinates of other sequences

cat PARs.txt ampliconic_subtracted.txt OTHER.txt repeats_subtracted.txt | sort -V -k1,1 -k2,2 >circos.all.sequence.classes.bed
python fill_gray_class_if_possible.py #fill the gray
python merge_consecutive_annotations.py #merge the redundant annotation (consecutive rows with the same annotation)
python drop_short_annotations.py #merge short annotations (<=100bp) with their neigborhood regions; output circos.all.sequence.classes.final.bed

#REWRITE ANNOTATION INTO A HUMAN READABLE FORMAT
grep "hs1" circos.all.sequence.classes.final.bed | sed -e 's/hs1/chrX/' -e 's/fill_color=#97cb99/PAR/' -e 's/fill_color=#89c0e8/AMPLICONIC/' -e 's/fill_color=#6a4a1d/SATELLITE/' -e 's/fill_color=#d9d8d8/OTHER/' -e 's/fill_color=#eea9ba/XTR/' >SEQUENCE_CLASSES.chrX.human.bed
grep "hs2" circos.all.sequence.classes.final.bed | sed -e 's/hs2/chrX/' -e 's/fill_color=#97cb99/PAR/' -e 's/fill_color=#89c0e8/AMPLICONIC/' -e 's/fill_color=#6a4a1d/SATELLITE/' -e 's/fill_color=#d9d8d8/OTHER/' -e 's/fill_color=#eea9ba/XTR/' >SEQUENCE_CLASSES.chrX.chimpanzee.bed
grep "hs3" circos.all.sequence.classes.final.bed | sed -e 's/hs3/chrX/' -e 's/fill_color=#97cb99/PAR/' -e 's/fill_color=#89c0e8/AMPLICONIC/' -e 's/fill_color=#6a4a1d/SATELLITE/' -e 's/fill_color=#d9d8d8/OTHER/' -e 's/fill_color=#eea9ba/XTR/' >SEQUENCE_CLASSES.chrX.bonobo.bed
grep "hs4" circos.all.sequence.classes.final.bed | sed -e 's/hs4/chrX/' -e 's/fill_color=#97cb99/PAR/' -e 's/fill_color=#89c0e8/AMPLICONIC/' -e 's/fill_color=#6a4a1d/SATELLITE/' -e 's/fill_color=#d9d8d8/OTHER/' -e 's/fill_color=#eea9ba/XTR/' >SEQUENCE_CLASSES.chrX.gorilla.bed
grep "hs5" circos.all.sequence.classes.final.bed | sed -e 's/hs5/chrX/' -e 's/fill_color=#97cb99/PAR/' -e 's/fill_color=#89c0e8/AMPLICONIC/' -e 's/fill_color=#6a4a1d/SATELLITE/' -e 's/fill_color=#d9d8d8/OTHER/' -e 's/fill_color=#eea9ba/XTR/' >SEQUENCE_CLASSES.chrX.sorang.bed
grep "hs6" circos.all.sequence.classes.final.bed | sed -e 's/hs6/chrX/' -e 's/fill_color=#97cb99/PAR/' -e 's/fill_color=#89c0e8/AMPLICONIC/' -e 's/fill_color=#6a4a1d/SATELLITE/' -e 's/fill_color=#d9d8d8/OTHER/' -e 's/fill_color=#eea9ba/XTR/' >SEQUENCE_CLASSES.chrX.borang.bed
grep "hs7" circos.all.sequence.classes.final.bed | sed -e 's/hs7/chrX/' -e 's/fill_color=#97cb99/PAR/' -e 's/fill_color=#89c0e8/AMPLICONIC/' -e 's/fill_color=#6a4a1d/SATELLITE/' -e 's/fill_color=#d9d8d8/OTHER/' -e 's/fill_color=#eea9ba/XTR/' >SEQUENCE_CLASSES.chrX.gibbon.bed

