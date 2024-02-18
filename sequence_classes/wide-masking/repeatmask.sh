#!/bin/bash
#SBATCH --job-name=repeatmask.20230627
#SBATCH --partition=main
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --output=repeatmask.20230627.%j.log
#SBATCH --time=168:00:00

set -e
set -x
pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

#LINE was identified by greping the original repeatmasker file from Gabby.
#For StSat, I used my coordinated from iteration 1 and NCRF.
#This is because for the chimpanzee, StSat was annotated as PTPCHT, so no StSat could be greped. It only found StSat in gorilla and bonobo. For consistency I decided to go with NCRF.
#I also want to mask GAP and SAT 

#notice that files in the format repeats.species.merged.bed are those with merged repeats revealing long stretches of predominantly satellites

#prepare the repeats track

#MERGE ALL REPEATS TO CREATE REPEAT TRACKS

#first convert .out format to .bed format
python RMOut-to-bed.py CHM13v2.0_ChrXY_FinalRepeatAnnotations_May2023.out RM.human.bed
python RMOut-to-bed.py Chimpanzee_FinalRepeatAnnotations.May2023.out RM.chimpanzee.bed
python RMOut-to-bed.py Bonobo_FinalRepeatAnnotations.May2023.out RM.bonobo.bed
python RMOut-to-bed.py Gorilla_FinalRepeatAnnotations.May2023.out RM.gorilla.bed
python RMOut-to-bed.py Sumatran_FinalRepeatAnnotations.May2023.out RM.sorang.bed
python RMOut-to-bed.py Bornean_FinalRepeatAnnotations.May2023.out RM.borang.bed
python RMOut-to-bed.py Siamang_FinalRepeatAnnotations.May2023.out RM.gibbon.bed

#separate into X and Y chromosome bed files
names=("human" "chimpanzee" "bonobo" "gorilla" "sorang" "borang" "gibbon")
for i in {0..6}; do echo $i; name=${names[i]}; cat RM.${name}.bed | grep "^chrY" >chrY.${name}.bed; done;
for i in {0..6}; do echo $i; name=${names[i]}; cat RM.${name}.bed | grep "^chrX" >chrX.${name}.bed; done;


#only large consecutive regions larger than 250000 will be kept and labeled as repeats
names=("human" "chimpanzee" "bonobo" "gorilla" "sorang" "borang" "gibbon")
for i in {0..6}; do echo $i; name=${names[i]}; cat chrY.${name}.bed | bedtools sort | bedtools merge -d 1000 | awk '{ if (($3-$2)>250000) print;}' | sort -V -k1,1 -k2,2 >repeats.chrY.${name}.merged.bed; done;
for i in {0..6}; do echo $i; name=${names[i]}; cat chrX.${name}.bed | bedtools sort | bedtools merge -d 1000 | awk '{ if (($3-$2)>250000) print;}' | sort -V -k1,1 -k2,2 >repeats.chrX.${name}.merged.bed; done;

#prepare the track for circos as well
#Y CHROMOSOME
cat repeats.chrY.human.merged.bed | sed -e 's/chrY/hs1/' >merged.chrY.repeats.species.bed.tmp
cat repeats.chrY.chimpanzee.merged.bed | sed -e 's/chrY/hs2/' >>merged.chrY.repeats.species.bed.tmp
cat repeats.chrY.bonobo.merged.bed | sed -e 's/chrY/hs3/' >>merged.chrY.repeats.species.bed.tmp
cat repeats.chrY.gorilla.merged.bed | sed -e 's/chrY/hs4/' >>merged.chrY.repeats.species.bed.tmp
cat repeats.chrY.sorang.merged.bed | sed -e 's/chrY/hs5/' >>merged.chrY.repeats.species.bed.tmp
cat repeats.chrY.borang.merged.bed | sed -e 's/chrY/hs6/' >>merged.chrY.repeats.species.bed.tmp
cat repeats.chrY.gibbon.merged.bed | sed -e 's/chrY/hs7/' >>merged.chrY.repeats.species.bed.tmp
cat merged.chrY.repeats.species.bed.tmp | awk '{print $0 "\tfill_color=#6a4a1d"}' | sed 's/ /\t/g' | sort -V -k1,1 -k2,2 >merged.chrY.repeats.species.bed
rm merged.chrY.repeats.species.bed.tmp
#X CHROMOSOME
cat repeats.chrX.human.merged.bed | sed -e 's/chrX/hs1/' >merged.chrX.repeats.species.bed.tmp
cat repeats.chrX.chimpanzee.merged.bed | sed -e 's/chrX/hs2/' >>merged.chrX.repeats.species.bed.tmp
cat repeats.chrX.bonobo.merged.bed | sed -e 's/chrX/hs3/' >>merged.chrX.repeats.species.bed.tmp
cat repeats.chrX.gorilla.merged.bed | sed -e 's/chrX/hs4/' >>merged.chrX.repeats.species.bed.tmp
cat repeats.chrX.sorang.merged.bed | sed -e 's/chrX/hs5/' >>merged.chrX.repeats.species.bed.tmp
cat repeats.chrX.borang.merged.bed | sed -e 's/chrX/hs6/' >>merged.chrX.repeats.species.bed.tmp
cat repeats.chrX.gibbon.merged.bed | sed -e 's/chrX/hs7/' >>merged.chrX.repeats.species.bed.tmp
cat merged.chrX.repeats.species.bed.tmp | awk '{print $0 "\tfill_color=#6a4a1d"}' | sed 's/ /\t/g' | sort -V -k1,1 -k2,2 >merged.chrX.repeats.species.bed
rm merged.chrX.repeats.species.bed.tmp


#GATHER THE HETEROCHROMATIC, REPETITVE PORTIONS OF THE Y, SO THAT THEY CAN BE MASKED AND ONLY EUCHROMATIN LEFT
cat chrY.human.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.human.bed
cat chrY.chimpanzee.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.chimpanzee.bed
cat chrY.bonobo.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.bonobo.bed
cat chrY.gorilla.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.gorilla.bed
cat chrY.sorang.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.sorang.bed
cat chrY.borang.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.borang.bed
cat chrY.gibbon.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.gibbon.bed

#GATHER THE HETEROCHROMATIC, REPETITVE PORTIONS OF THE X, SO THAT THEY CAN BE MASKED AND ONLY EUCHROMATIN LEFT
cat chrX.human.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrX.human.bed
cat chrX.chimpanzee.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrX.chimpanzee.bed
cat chrX.bonobo.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrX.bonobo.bed
cat chrX.gorilla.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrX.gorilla.bed
cat chrX.sorang.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrX.sorang.bed
cat chrX.borang.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrX.borang.bed
cat chrX.gibbon.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrX.gibbon.bed


#Y CHROMOSOME, WIDE MASKING, AFTER WHICH ONLY EUCHROMATIC PROPORTION SHOULD REMAIN
cat repeats.chrY.human.merged.bed sat_gap_interspersed.chrY.human.bed HSAT.chrY.human.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >human.chrY.WideHardmasked.bed
cat repeats.chrY.chimpanzee.merged.bed sat_gap_interspersed.chrY.chimpanzee.bed StSat.chrY.chimpanzee.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >chimpanzee.chrY.WideHardmasked.bed
cat repeats.chrY.bonobo.merged.bed sat_gap_interspersed.chrY.bonobo.bed StSat.chrY.bonobo.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >bonobo.chrY.WideHardmasked.bed
cat repeats.chrY.gorilla.merged.bed sat_gap_interspersed.chrY.gorilla.bed StSat.chrY.gorilla.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >gorilla.chrY.WideHardmasked.bed
cat repeats.chrY.sorang.merged.bed sat_gap_interspersed.chrY.sorang.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >sorang.chrY.WideHardmasked.bed
cat repeats.chrY.borang.merged.bed sat_gap_interspersed.chrY.borang.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >borang.chrY.WideHardmasked.bed
cat repeats.chrY.gibbon.merged.bed sat_gap_interspersed.chrY.gibbon.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >gibbon.chrY.WideHardmasked.bed

#Y CHROMOSOME, hardmask the problematic repeat regions
bedtools maskfasta -fi chrY.human.fasta -bed human.chrY.WideHardmasked.bed -fo chrY.human.WideHardmasked.fasta
bedtools maskfasta -fi chrY.chimpanzee.fasta -bed chimpanzee.chrY.WideHardmasked.bed -fo chrY.chimpanzee..WideHardmasked.fasta
bedtools maskfasta -fi chrY.bonobo.fasta -bed bonobo.chrY.WideHardmasked.bed -fo chrY.bonobo.WideHardmasked.fasta
bedtools maskfasta -fi chrY.gorilla.fasta -bed gorilla.chrY.WideHardmasked.bed -fo chrY.gorilla.WideHardmasked.fasta
bedtools maskfasta -fi chrY.sorang.fasta -bed sorang.chrY.WideHardmasked.bed -fo chrY.sorang.WideHardmasked.fasta
bedtools maskfasta -fi chrY.borang.fasta -bed borang.chrY.WideHardmasked.bed -fo chrY.borang.WideHardmasked.fasta
bedtools maskfasta -fi chrY.gibbon.fasta -bed gibbon.chrY.WideHardmasked.bed -fo chrY.gibbon.WideHardmasked.fasta

#X CHROMOSOME, WIDE MASKING, AFTER WHICH ONLY EUCHROMATIC PROPORTION SHOULD REMAIN
cat repeats.chrX.human.merged.bed sat_gap_interspersed.chrX.human.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >human.chrX.WideHardmasked.bed
cat repeats.chrX.chimpanzee.merged.bed sat_gap_interspersed.chrX.chimpanzee.bed StSat.chrX.chimpanzee.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >chimpanzee.chrX.WideHardmasked.bed
cat repeats.chrX.bonobo.merged.bed sat_gap_interspersed.chrX.bonobo.bed StSat.chrX.bonobo.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >bonobo.chrX.WideHardmasked.bed
cat repeats.chrX.gorilla.merged.bed sat_gap_interspersed.chrX.gorilla.bed StSat.chrX.gorilla.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >gorilla.chrX.WideHardmasked.bed
cat repeats.chrX.sorang.merged.bed sat_gap_interspersed.chrX.sorang.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >sorang.chrX.WideHardmasked.bed
cat repeats.chrX.borang.merged.bed sat_gap_interspersed.chrX.borang.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >borang.chrX.WideHardmasked.bed
cat repeats.chrX.gibbon.merged.bed sat_gap_interspersed.chrX.gibbon.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >gibbon.chrX.WideHardmasked.bed

#X CHROMOSOME, hardmask the problematic repeat regions
bedtools maskfasta -fi chrX.human.fasta -bed human.chrX.WideHardmasked.bed -fo chrX.human.WideHardmasked.fasta
bedtools maskfasta -fi chrX.chimpanzee.fasta -bed chimpanzee.chrX.WideHardmasked.bed -fo chrX.chimpanzee.WideHardmasked.fasta
bedtools maskfasta -fi chrX.bonobo.fasta -bed bonobo.chrX.WideHardmasked.bed -fo chrX.bonobo.WideHardmasked.fasta
bedtools maskfasta -fi chrX.gorilla.fasta -bed gorilla.chrX.WideHardmasked.bed -fo chrX.gorilla.WideHardmasked.fasta
bedtools maskfasta -fi chrX.sorang.fasta -bed sorang.chrX.WideHardmasked.bed -fo chrX.sorang.WideHardmasked.fasta
bedtools maskfasta -fi chrX.borang.fasta -bed borang.chrX.WideHardmasked.bed -fo chrX.borang.WideHardmasked.fasta
bedtools maskfasta -fi chrX.gibbon.fasta -bed gibbon.chrX.WideHardmasked.bed -fo chrX.gibbon.WideHardmasked.fasta


echo "Done."
