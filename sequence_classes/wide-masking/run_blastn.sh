#!/bin/bash
set -e
pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

threads=128

assembly=$1
assembly_name="$(basename ${assembly})"
assembly_name="${assembly%.*}"

reference=$2

#split fasta query for better parallelization
fasta-splitter --n-parts 9 ${assembly}

parts=(1 2 3 4 5 6 7 8 9) #chromosomes to use
for part in "${parts[@]}"; do 
    echo $part;  
    time srun --nodes=1 --ntasks=128 --mem=500gb --time=INFINITE blastn -query ${assembly_name}.part-${part}.fasta -db ${reference} -outfmt 6 -perc_identity 50 -num_threads ${threads} -out merged.blastn.${assembly}.${reference}.part${part}.txt &
done; 
wait

#merge blastn results
cat merged.blastn.${assembly}.${reference}.part*.txt >blastn.${assembly}.${reference}.txt
rm merged.blastn.${assembly}.${reference}.part*.txt
rm ${assembly_name}.part-*.fasta
cat blastn.${assembly}.${reference}.txt | awk '{ if ($4>4750) print;}' >blastn.long.${assembly}.${reference}.txt
#cat blastn.long.${assembly}.${reference}.txt | sort -gk9 | awk '{print $2 "\t" $9 "\t" $10 "\t" $3}' >blastn.long.${assembly}.${reference}.bed

echo "Done."
