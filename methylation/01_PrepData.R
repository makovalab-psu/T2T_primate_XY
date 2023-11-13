
##### 01Start: Obtain 18k points methylation profiles of CATannot genes #####
l=ls()
l=ls()
l=l[which(!ls() %in% c("",""))]
rm(list=l)

# load filepath of methylation and annotation files 
methfiletab<-read.table("Data/Methylation/methylation_filestab.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F)
methfiletab
annotfiletab<-read.table("Data/Annotation/annotation_filestab_CATannot.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F)
annotfiletab

methfiletab$common_name==annotfiletab$common_name

# loop processing for each species
for (i in 1:7) {
  #i=1 i=5
  tempname=annotfiletab$common_name[i]
  tempoutfile=paste("Data/Methylation/gene_CATannot_methprofile_smoothed.",tempname,".txt",sep="")
  
  tempmeth=read.table(methfiletab[methfiletab$common_name==tempname,"methfilename"],header=F,sep="\t",stringsAsFactors=F,check.names=F)
  names(tempmeth)=c("chr","start","end","methylation_percent")
  tempmeth$nucpos=tempmeth$end # Note: GFF/GTF is 1-based. BED is 0-based for the coordinate start and 1-based for the coordinate end, so for BED the "end" is he actual nucleotide position that needs to be checked against GFF/GTF
  head(tempmeth)  
  
  tempannot=read.table(annotfiletab[annotfiletab$common_name==tempname,"annotfilename"],header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="")
  tempannot[tempannot$strand=="+","P_start"]=tempannot[tempannot$strand=="+","start"]-3000
  tempannot[tempannot$strand=="+","P_end"]=tempannot[tempannot$strand=="+","start"]-1
  tempannot[tempannot$strand=="+","GB_start"]=tempannot[tempannot$strand=="+","start"]
  tempannot[tempannot$strand=="+","GB_end"]=tempannot[tempannot$strand=="+","end"]
  tempannot[tempannot$strand=="+","I_start"]=tempannot[tempannot$strand=="+","end"]+1
  tempannot[tempannot$strand=="+","I_end"]=tempannot[tempannot$strand=="+","end"]+3000
  tempannot[tempannot$strand=="+","TSS"]=tempannot[tempannot$strand=="+","start"]
  tempannot[tempannot$strand=="+","TES"]=tempannot[tempannot$strand=="+","end"]
  tempannot[tempannot$strand=="-","P_start"]=tempannot[tempannot$strand=="-","end"]+1
  tempannot[tempannot$strand=="-","P_end"]=tempannot[tempannot$strand=="-","end"]+3000
  tempannot[tempannot$strand=="-","GB_start"]=tempannot[tempannot$strand=="-","start"]
  tempannot[tempannot$strand=="-","GB_end"]=tempannot[tempannot$strand=="-","end"]
  tempannot[tempannot$strand=="-","I_start"]=tempannot[tempannot$strand=="-","start"]-3000
  tempannot[tempannot$strand=="-","I_end"]=tempannot[tempannot$strand=="-","start"]-1
  tempannot[tempannot$strand=="-","TSS"]=tempannot[tempannot$strand=="-","end"]
  tempannot[tempannot$strand=="-","TES"]=tempannot[tempannot$strand=="-","start"]
  
  head(tempannot)
  
  summary(tempannot$end-tempannot$start+1)
  #median human genebody length is ~10kb, so let's have the genebody representation/visualization be 4 times the 3kb promoter.
  tempmedianlength=12000
  
  for (j in 1:length(row.names(tempannot))) {
    #j=1  j=2  j=21
    temptss=tempannot$TSS[j]
    temptes=tempannot$TES[j]
    tempstrand=tempannot$strand[j]
    tempgenelength=(tempannot$end[j]-tempannot$start[j]+1)
    
    tempprom=tempmeth[(tempmeth$chr==tempannot[j,"seqid"])&(tempmeth$nucpos>=tempannot[j,"P_start"])&(tempmeth$nucpos<=tempannot[j,"P_end"]),]
    tempgb=tempmeth[(tempmeth$chr==tempannot[j,"seqid"])&(tempmeth$nucpos>=tempannot[j,"GB_start"])&(tempmeth$nucpos<=tempannot[j,"GB_end"]),]
    tempinter=tempmeth[(tempmeth$chr==tempannot[j,"seqid"])&(tempmeth$nucpos>=tempannot[j,"I_start"])&(tempmeth$nucpos<=tempannot[j,"I_end"]),]
    
    if (tempstrand=="+") {
      tempprom$truedist2TSS=tempprom$nucpos-temptss
      tempprom$relativedist2TSS=tempprom$truedist2TSS
      #head(tempprom)
      #plot(tempprom$relativedist2TSS,tempprom$methylation_percent)
      tempgb$truedist2TSS=tempgb$nucpos-temptss
      tempgb$relativedist2TSS=tempgb$truedist2TSS*tempmedianlength/tempgenelength
      #head(tempgb)
      #tail(tempgb)
      #plot(tempgb$relativedist2TSS,tempgb$methylation_percent)
      tempinter$truedist2TSS=tempinter$nucpos-temptss
      tempinter$relativedist2TSS=tempinter$truedist2TSS-((temptes-temptss)-12000)
      #head(tempinter)
      #tail(tempinter)
      #plot(tempinter$relativedist2TSS,tempinter$methylation_percent)
    }else if (tempstrand=="-") {
      tempprom$truedist2TSS=temptss-tempprom$nucpos
      tempprom$relativedist2TSS=tempprom$truedist2TSS
      #head(tempprom)
      #plot(tempprom$relativedist2TSS,tempprom$methylation_percent)
      tempgb$truedist2TSS=temptss-tempgb$nucpos
      tempgb$relativedist2TSS=tempgb$truedist2TSS*tempmedianlength/tempgenelength
      #head(tempgb)
      #tail(tempgb)
      #plot(tempgb$relativedist2TSS,tempgb$methylation_percent)
      tempinter$truedist2TSS=temptss-tempinter$nucpos
      tempinter$relativedist2TSS=tempinter$truedist2TSS-((temptss-temptes)-12000)
      #head(tempinter)
      #plot(tempinter$relativedist2TSS,tempinter$methylation_percent)
    }else {
      stop("strand error!")
    }
    tempcomb=rbind(tempprom,tempgb,tempinter)
    #head(tempcomb)
    #plot(tempcomb$relativedist2TSS,tempcomb$methylation_percent)
    
    bin=500
    tempcomb2=data.frame(pos=seq(-2999,15000,1))
    for (pos in seq(-2999,15000,1)) {
      #pos=-2999
      tempcomb2[tempcomb2$pos==pos,"methylation_percent"]=mean(tempcomb[(tempcomb$relativedist2TSS>=(pos-bin/2))&(tempcomb$relativedist2TSS<=(pos+bin/2)),"methylation_percent"])
    }
    #head(tempcomb2)
    #plot(tempcomb2$pos,tempcomb2$methylation_percent,pch=20,cex=0.8)
    tempannot[j,paste("pos_",1:18000,sep="")]=tempcomb2$methylation_percent
    #head(tempannot)
    #tempannot[1:5,1:30]
    #tempannot[1:5,18000:18024]
    
    if (j%%10==0) {
      cat(i,"_",j,"_",date(),"\n",sep="")
      write.table(tempannot,tempoutfile,quote=F,sep="\t",row.names=F,col.names=T)
    }
  }  
  write.table(tempannot,tempoutfile,quote=F,sep="\t",row.names=F,col.names=T)
  cat(i,"_",tempname,"_done \n",sep="")
}
##### 01End: Obtain 18k points methylation profiles of CATannot genes #####



##### 02Start: Prep plotdata for repeats and sequence regions #####

### start: Calculate/Plot methylation at repeats and genes ###
l=ls()
l=ls()
l=l[which(!ls() %in% c("",""))]
rm(list=l)

samples=c("HG002","CHM13","Bonobo","Chimpanzee",
          "Gorilla","Bornean","Sumatran","Siamang")
methfiles2=c("Data/Methylation/HG002_hg002XYv2.7_hg002_CpG_ont_guppy6.1.2.Methylation.mergechrXY_withPAR.bedgraph",
            "Data/Methylation/CHM13_chm13v2.0_hg002_CpG_ont_guppy6.1.2.Methylation.mergechrXY_withPAR.bedgraph",
            "Data/Methylation/Bonobo_ONT_Methylation.mergechrXY.filter_withPAR.bedgraph",
            "Data/Methylation/Chimpanzee_ONT_Methylation.mergechrXY.filter_withPAR.bedgraph",
            "Data/Methylation/Gorilla_ONT_Methylation.mergechrXY.filter_withPAR.bedgraph",
            "Data/Methylation/Bornean_ONT_Methylation.mergechrXY.filter_withPAR.bedgraph",
            "Data/Methylation/Sumatran_ONT_Methylation.mergechrXY.filter_withPAR.bedgraph",
            "Data/Methylation/Siamang_ONT_Methylation.mergechrXY.filter_withPAR.bedgraph")
repeatfiles2=c("Data/Repeats/methylation_HG002_XY_eddjoined_FinalAnnotations_June2023_v2.sort.out_withMethCalc.txt",
              "Data/Repeats/methylation_CHM13v2.0_XY_FinalAnnotations_June2023_v2.sort.out_withMethCalc.txt",
              "Data/Repeats/methylation_Bonobo_FinalAnnotations_June2023_v2.sort.out_withMethCalc.txt",
              "Data/Repeats/methylation_Chimpanzee_FinalAnnotations_June2023_v2.sort.out_withMethCalc.txt",
              "Data/Repeats/methylation_Gorilla_FinalAnnotations_June2023_v2.sort.out_withMethCalc.txt",
              "Data/Repeats/methylation_Bornean_FinalAnnotations_June2023_v2.sort.out_withMethCalc.txt",
              "Data/Repeats/methylation_Sumatran_FinalAnnotations_June2023_v2.sort.out_withMethCalc.txt",
              "Data/Repeats/methylation_Siamang_FinalAnnotations_June2023_v2.sort.out_withMethCalc.txt")
geneannotfiles2=c("Data/Annotation/chrXYgenes_CATliftoff_HG002_withMethCalc.txt",
                 "Data/Annotation/chrXYgenes_CATliftoff_CHM13_withMethCalc.txt",
                 "Data/Annotation/chrXYgenes_CATliftoff_Bonobo_withMethCalc.txt",
                 "Data/Annotation/chrXYgenes_CATliftoff_Chimpanzee_withMethCalc.txt",
                 "Data/Annotation/chrXYgenes_CATliftoff_Gorilla_withMethCalc.txt",
                 "Data/Annotation/chrXYgenes_CATliftoff_Bornean_withMethCalc.txt",
                 "Data/Annotation/chrXYgenes_CATliftoff_Sumatran_withMethCalc.txt",
                 "Data/Annotation/chrXYgenes_CATliftoff_Siamang_withMethCalc.txt")
parinfoFile  <- "DATA/PAR/PAR_20230623.txt"
parinfotab<-read.table(parinfoFile,header=F,sep="\t",stringsAsFactors=F,check.names=F,comment.char="#")
parinfotab=parinfotab[,1:4]
names(parinfotab)=c("name","chr","start","end")
parinfotab[grep("Gor",parinfotab$name),"common_name"]="Gorilla"
parinfotab[grep("PanPan",parinfotab$name),"common_name"]="Bonobo"
parinfotab[grep("PanTro",parinfotab$name),"common_name"]="Chimpanzee"
parinfotab[grep("PonAbe",parinfotab$name),"common_name"]="Sumatran"
parinfotab[grep("PonPyg",parinfotab$name),"common_name"]="Bornean"
parinfotab[grep("Sym",parinfotab$name),"common_name"]="Siamang"
parinfotab[grep("CHM13",parinfotab$name),"common_name"]="CHM13"
parinfotab[grep("HG002",parinfotab$name),"common_name"]="HG002"
parinfotab[grep("HG38",parinfotab$name),"common_name"]="HG38"
parinfotab

library(GenomicRanges)

###Prep plotdata for repeats ###
l=ls()
l=ls()
l=l[which(!ls() %in% c("samples","methfiles2","repeatfiles2",
                       "geneannotfiles2","parinfoFile","parinfotab"))]
rm(list=l)

methfiles2
repeatfiles2
geneannotfiles2
parinfotab
for (i in 1:length(samples)) {
  #i=1  i=2  i=4
  tempsample=samples[i]
  tempsample
  
  temprepeat=read.table(repeatfiles2[i],header=T,sep="\t",stringsAsFactors=F,check.names=F)
  tempgene=read.table(geneannotfiles2[i],header=T,sep="\t",stringsAsFactors=F,check.names=F)
  if (i<=2) {
    temprepeat=temprepeat[!((temprepeat$chr=="chrY")&(temprepeat$end>=27449937)),]
    tempgene=tempgene[!((tempgene$seqid=="chrY")&(tempgene$end>=27449937)),]
  }
  head(temprepeat)
  head(tempgene)    
  
  tempparinfo=parinfotab[parinfotab$common_name==tempsample,]  
  tempparinfo$name=sub("^\\S+PAR","PAR",tempparinfo$name)
  tempparinfo$name=sub("[XY]$","",tempparinfo$name)
  tempparinfo$name=sub("PAR$","PAR1",tempparinfo$name)
  tempparinfo
  
  table(temprepeat$in_PAR)
  table(tempgene$in_PAR)
  for (j in 1:length(row.names(tempparinfo))) {
    #j=1
    tempparname=tempparinfo[j,"name"]
    tempparchr=tempparinfo[j,"chr"]
    tempparstart=tempparinfo[j,"start"]
    tempparend=tempparinfo[j,"end"]
    temprepeat[(temprepeat$chr==tempparchr)&(temprepeat$start<=tempparstart)&(temprepeat$end>=tempparend),"in_PAR_b"]=tempparname
    temprepeat[(temprepeat$chr==tempparchr)&(temprepeat$end<=tempparend)&(temprepeat$end>=tempparstart),"in_PAR_b"]=tempparname
    temprepeat[(temprepeat$chr==tempparchr)&(temprepeat$start>=tempparstart)&(temprepeat$start<=tempparend),"in_PAR_b"]=tempparname
    
    tempgene[(tempgene$seqid==tempparchr)&(tempgene$start<=tempparstart)&(tempgene$end>=tempparend),"in_PAR_b"]=tempparname
    tempgene[(tempgene$seqid==tempparchr)&(tempgene$end<=tempparend)&(tempgene$end>=tempparstart),"in_PAR_b"]=tempparname
    tempgene[(tempgene$seqid==tempparchr)&(tempgene$start>=tempparstart)&(tempgene$start<=tempparend),"in_PAR_b"]=tempparname
    
  }
  temprepeat[is.na(temprepeat$in_PAR_b),"in_PAR_b"]="No"
  tempgene[is.na(tempgene$in_PAR_b),"in_PAR_b"]="No"
  table(temprepeat$in_PAR_b)
  table(tempgene$in_PAR_b)
  head(temprepeat)
  head(tempgene)
  
  table(temprepeat$repeat_groups)
  
  keepgroups=c("SINE","LINE","LTR","DNA","Satellite","Simple/LowComplexity")
  temprepeat1=temprepeat[temprepeat$repeat_groups %in% keepgroups,]
  temprepeat1$fordavid_chr=temprepeat1$chr
  temprepeat1[temprepeat1$in_PAR_b!="No","fordavid_chr"] = temprepeat1[temprepeat1$in_PAR_b!="No","in_PAR_b"]
  temprepeat1$fordavid_POS=floor((temprepeat1$start+temprepeat1$end)/2)
  temprepeat1$fordavid_Meth=temprepeat1$mean_methylation
  temprepeat1$fordavid_type=temprepeat1$repeat_groups
  temprepeat1$fordavid_par=sub("PAR","",temprepeat1$in_PAR_b)
  temprepeat1[temprepeat1$fordavid_par=="No","fordavid_par"]="0"
  temprepeat1$fordavid_par=as.numeric(temprepeat1$fordavid_par)
  temprepeat1$fordavid_par
  table(temprepeat1$fordavid_par)
  temprepeat1$fordavid_species=tempsample
  temprepeat1=temprepeat1[,grep("fordavid",names(temprepeat1))]
  head(temprepeat1)
  
  table(tempgene$in_PAR)
  table(tempgene$in_PAR_b)
  tempgene1=tempgene
  tempgene1$fordavid_chr=tempgene1$seqid
  tempgene1[tempgene1$in_PAR_b!="No","fordavid_chr"] = tempgene1[tempgene1$in_PAR_b!="No","in_PAR_b"]
  tempgene1$fordavid_POS=floor((tempgene1$start+tempgene1$end)/2)
  tempgene1$fordavid_Meth=tempgene1$gb_norepeats_meth
  tempgene1$fordavid_type="GenesExclRepeats"
  tempgene1$fordavid_par=sub("PAR","",tempgene1$in_PAR_b)
  tempgene1[tempgene1$fordavid_par=="No","fordavid_par"]="0"
  tempgene1$fordavid_par=as.numeric(tempgene1$fordavid_par)
  tempgene1$fordavid_par
  table(tempgene1$fordavid_par)
  tempgene1$fordavid_species=tempsample
  tempgene1=tempgene1[,grep("fordavid",names(tempgene1))]
  head(tempgene1)
  
  tempcomb=rbind(temprepeat1,tempgene1)
  names(tempcomb)=sub("fordavid_","",names(tempcomb))
  if (mean(tempcomb$Meth,na.rm=T)>1) {
    tempcomb$Meth=tempcomb$Meth/100
  }
  head(tempcomb)
  table(tempcomb$chr)
  hist(tempcomb$POS)
  hist(tempcomb$Meth)
  table(tempcomb$type)
  table(tempcomb$par)
  table(tempcomb$species)
  
  write.table(tempcomb,paste("Data/PlotData/repeats_plotdata_",tempsample,".csv",sep=""),quote=T,sep=",",row.names=F,col.names=T)
  cat(i," ",sep="")  
}

###Prep plotdata for sequence regions ###
l=ls()
l=ls()
l=l[which(!ls() %in% c(""))]
rm(list=l)

samples=c("HG002","CHM13","Bonobo","Chimpanzee",
          "Gorilla","Bornean","Sumatran","Siamang")
uniqclassesX=c("AMPLICONIC_chrX","PAR_chrX","X-ANCESTRAL_chrX")
uniqclassesY=c("AMPLICONIC_chrY","PAR_chrY","XDEG_chrY")
uniqclasses=c(uniqclassesX,uniqclassesY)
for (i in 1:8) {
  #i=1 i=2 i=3
  options(scipen = 999)
  tempsample=samples[i]
  infile=paste("Data/SequenceClasses/seqclassbinmeth_",tempsample,".txt",sep="")
  if (i<=2) {
    infile=paste("Data/SequenceClasses/seqclassbinmeth_20230926_noYq12_",tempsample,".txt",sep="")
  }
  temptab=read.table(infile,header=T,sep="\t",stringsAsFactors=F,check.names=F)
  table(temptab$seqclass)
  temptab=temptab[temptab$seqclass %in% uniqclasses,]
  head(temptab)
  temptab$fordavid_chr=temptab$chr
  temptab$fordavid_POS=floor((temptab$start+temptab$end)/2)
  if (mean(temptab$methylation,na.rm=T)>1) {
    temptab$fordavid_Meth=temptab$methylation/100
  }
  temptab$fordavid_type=sub("\\_chr[XY]$","",temptab$seqclass)
  temptab[(temptab$fordavid_type=="X-ANCESTRAL")|(temptab$fordavid_type=="XDEG"),"fordavid_type"]="X-ANCESTRAL/DEG"
  temptab$fordavid_par=0
  temptab[grep("PAR",temptab$fordavid_type),"fordavid_par"]=1
  temptab$fordavid_species=tempsample
  temptab=temptab[,grep("fordavid",names(temptab))]
  names(temptab)=sub("fordavid_","",names(temptab))
  head(temptab)
  table(temptab$chr)
  hist(temptab$POS)
  hist(temptab$Meth)
  table(temptab$type)
  table(temptab$par)
  temptab$par
  table(temptab$species)
  
  write.table(temptab,paste("Data/PlotData/seqregions_plotdata_",tempsample,".csv",sep=""),quote=T,sep=",",row.names=F,col.names=T)
  cat(i," ",sep="")
}

##### 02End: Prep plotdata for repeats and sequence regions #####




