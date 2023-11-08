""""Snakefile for Counting rDNA Copy Number"""

# this is where the config file lives
configfile: "config.yml"

import re
import os




rule match_gc:
    input:
        rdna = expand("rDNA_references/rdna_{species}_18S.fa",species=config["ID"])
    params:
        species = config["ID"]
    output:
        match_gcn = expand("matched_windows/matched_windows_subset_{species}.fa", species=config["ID"])
    run:
        for species in params.species:
            shell("python Identify_GC_matched_regions_primates.py -roi rDNA_references/rdna_{species}_18S.fa -assembly genomes/{species}_singleline.fa")
            shell("bedtools getfasta -fi genomes/{species}.fa -bed matched_windows_subset.bed -fo matched_windows_subset_{species}.fa")
            shell("mv matched_windows_subset.bed matched_windows_subset_{species}.bed")
            shell("mv matched_windows.bed matched_windows_{species}.bed")
            shell("mv matched_windows*bed matched_windows/")
            shell("mv matched_windows_subset_{species}.fa matched_windows/")


rule kmerize_rdna:
    input:
        rdna = expand("rDNA_references/rdna_{species}_18S.fa",species=config["ID"])
    params:
        k = config["K"],
        species = config["ID"]
    output:
        gcn = expand("rDNA_references/{species}_18s_k{k}_cn.fa",k=config["K"],species=config["ID"])
    run:
        for species in params.species:
            for k in params.k:
                shell("jellyfish count -m {k} -s 100M -t 5 -C -o rDNA_references/{species}_18s_k{k}_cn.jf rDNA_references/rdna_{species}_18S.fa")
                shell("jellyfish dump rDNA_references/{species}_18s_k{k}_cn.jf > rDNA_references/{species}_18s_k{k}_cn.fa")

rule uniq_kmerize_matched_windows:
    input:
        match_gcn = expand("matched_windows/matched_windows_subset_{species}.fa", species=config["ID"])
    params:
        k = config["K"],
        species = config["ID"]
    output:
        ngcn = expand("matched_windows/matched_windows_subset_{species}_fcn_unique_{k}mers.fa",k=config["K"],species=config["ID"])
    run:
        for species in params.species:
            shell("bedtools sort -i matched_windows/matched_windows_subset_{species}.bed -g genomes/{species}.genome | bedtools complement -i - -g genomes/{species}.genome > matched_windows/all_but_matched_windows_subset_{species}.bed")
            shell("bedtools getfasta -fi genomes/gorilla.fa -bed matched_windows/all_but_matched_windows_subset_{species}.bed > matched_windows/all_but_matched_windows_subset_{species}.fa")
            for k in params.k:
                shell("jellyfish count -m {k} -s 100M -C -t 10 -o matched_windows_subset_{species}_nfcn_{k}mers.jf --if matched_windows/matched_windows_subset_{species}.fa matched_windows/all_but_matched_windows_subset_{species}.fa")
                shell("jellyfish dump -o matched_windows_subset_{species}_nfcn_{k}mers.fa matched_windows_subset_{species}_nfcn_{k}mers.jf")
                shell("grep -A1 -w '>0' matched_windows_subset_{species}_nfcn_{k}mers.fa > matched_windows_subset_{species}_unique_{k}mers.fa")
                shell("grep -v '\-' matched_windows_subset_{species}_unique_{k}mers.fa > matched_windows_subset_{species}_unique_{k}mers.fa1")
                shell("rm matched_windows_subset_{species}_unique_{k}mers.fa")
                shell("mv matched_windows_subset_{species}_unique_{k}mers.fa1 matched_windows_subset_{species}_unique_{k}mers.fa")
                shell("jellyfish count -m {k} -s 100M -C -t 10 -o matched_windows_subset_{species}_fcn_unique_{k}mers.jf --if matched_windows_subset_{species}_unique_{k}mers.fa matched_windows/matched_windows_subset_{species}.fa")
                shell("jellyfish dump -o matched_windows_subset_{species}_fcn_unique_{k}mers.fa matched_windows_subset_{species}_fcn_unique_{k}mers.jf")
                shell("mv matched_windows_subset_{species}*{k}mers.fa matched_windows/")
                shell("mv matched_windows_subset_{species}_* matched_windows/")


rule count_cn:
    input:
        ngcn = expand("matched_windows/matched_windows_subset_{species}_fcn_unique_{k}mers.fa",k=config["K"],species=config["ID"]),
        gcn = expand("rDNA_references/{species}_18s_k{k}_cn.fa",k=config["K"],species=config["ID"])
    params:
        species = config["ID"],
        k = config["K"]
    output:
        CN = expand("results/Copy_Numbers_{species}_k{k}.tsv",species=config["ID"],k=config["K"])
    run:
        for species in params.species:
            for k in params.k:
                shell("jellyfish count -m {k} -s 100M -C -t 15 -o {species}_k{k}_1.jf --if rDNA_references/rdna_{species}_18S.fa WGS/{species}/{species}_1.fastq")
                shell("jellyfish count -m {k} -s 100M -C -t 15 -o {species}_k{k}_2.jf --if rDNA_references/rdna_{species}_18S.fa WGS/{species}/{species}_2.fastq")
                shell("jellyfish dump -o {species}_k{k}_1.fa {species}_k{k}_1.jf")
                shell("jellyfish dump -o {species}_k{k}_2.fa {species}_k{k}_2.jf")
                shell("mv {species}* jellyfish_files/")
                shell("jellyfish count -m {k} -s 100M -C -t 15 -o matched_gc_{species}_k{k}_1.jf --if matched_windows/matched_windows_subset_{species}_unique_{k}mers.fa WGS/{species}/{species}_1.fastq")
                shell("jellyfish dump -o matched_gc_{species}_k{k}_1.fa matched_gc_{species}_k{k}_1.jf")
                shell("jellyfish count -m {k} -s 100M -C -t 15 -o matched_gc_{species}_k{k}_2.jf --if matched_windows/matched_windows_subset_{species}_unique_{k}mers.fa WGS/{species}/{species}_2.fastq")
                shell("jellyfish dump -o matched_gc_{species}_k{k}_2.fa matched_gc_{species}_k{k}_2.jf")
                shell("mv matched_gc_{species}* jellyfish_files/")
                shell("python Call_Copy_Number_GC_Normalization_Version5_primates.py -r1 jellyfish_files/{species}_k{k}_1.fa -r2 jellyfish_files/{species}_k{k}_2.fa -nc1 jellyfish_files/matched_gc_{species}_k{k}_1.fa -nc2 jellyfish_files/matched_gc_{species}_k{k}_2.fa -ngcn matched_windows/matched_windows_subset_{species}_fcn_unique_{k}mers.fa -ID {species} -gcn rDNA_references/{species}_18s_k{k}_cn.fa")
                shell("mv Copy_Numbers.tsv Copy_Numbers_{species}_k{k}.tsv")
                shell("mv Copy_Numbers_{species}_k{k}.tsv results/")
