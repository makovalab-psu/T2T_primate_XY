#!/usr/bin/env Rscript
#this script gets blastn results as an input, and outputs the coordinates of regions with high intrachromosomal similarity per respective basepairs

require(tidyr)
require(dplyr)
options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
print(file)

plotRectangle<-function(species) {
  species<-as.data.frame(read.table(paste0(species,".dotamplicons.bed"),sep="\t",header=FALSE,comment.char = ""))
colnames(species)<-c("chr","start","end")
rect(species$start/coef,50,species$end/coef,100,col=adjustcolor("blue", alpha.f = 0.1),border="NA")
}

parseBlastn <- function(file) { # create a function with the name my_function
    #load the alignment file
  blastn<-as.data.frame(read.table(file,skip=0))
  colnames(blastn)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  
  #add the actual coordinates
  origin<-separate(data = blastn, col = qseqid, into = c("chr", "coordinates"), sep = "\\:")
  origin<-separate(data = origin, col = coordinates, into = c("start", "end"), sep = "\\-")
  origin$start<-as.numeric(as.character(origin$start))
  origin$end<-as.numeric(as.character(origin$end))
  
  #drop self-alignments
  start_within<-((origin$sstart-1)>=origin$start) & ((origin$sstart-1)<=origin$end) #the start coordinate is within the original coordinates from which the sequence is from 
  end_within<-((origin$send)>=origin$start) & ((origin$send)<=origin$end) #the start coordinate is within the original coordinates from which the sequence is from 
  
  origin_without_self<-origin[!(start_within | end_within),] #exclude mappings to the coordinates from which the sequence was subsampled
  
  transformed_identity<-origin_without_self %>% group_by(start) %>% 
  arrange(pident, .by_group = TRUE) %>%
  mutate(max_identity = max(pident))
  
  bed_file<-transformed_identity[,c("chr","start","end","max_identity")]
  bed_file<-bed_file %>% distinct()
  
  return(bed_file)
}

outputfilename<-paste0(file,".bed")
outputfilename<-gsub(".txt.bed", ".bed", outputfilename, ignore.case = TRUE)
bed<-parseBlastn(file)
write.table(bed,file=outputfilename,sep="\t",col.names = FALSE,row.names=FALSE,quote=FALSE)
