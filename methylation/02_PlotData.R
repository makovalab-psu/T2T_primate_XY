##### 01Start: Plot repeats and sequence classes #####
l=ls()
l=ls()
l=l[which(!ls() %in% c("",""))]
rm(list=l)

library(ggplot2)
library(ggsignif)
library(scales)

files=dir("Data/PlotData",full.names=T)
files=files[grep("csv$",files)]
files

### Plot repeats ###
repeatfiles=files[grepl("repeats_plotdata",files)]
repeatfiles
for (i in 1:length(repeatfiles)) {
  #i=1
  temptab <- read.csv(repeatfiles[i], stringsAsFactors = FALSE)
  head(temptab)
  if (i==1) {
    plotdata=temptab
  }else {
    plotdata=rbind(temptab,plotdata)
  }
  numrows=length(row.names(plotdata))
  cat(i,"_",numrows," ",sep="")
}
plotdata=plotdata[plotdata$species != "CHM13",]
plotdata[plotdata$species=="HG002","species"]="Human"
plotdata[plotdata$species=="Bornean","species"]="B.orangutan"
plotdata[plotdata$species=="Sumatran","species"]="S.orangutan"
plotdata$chr=sub("chr","",plotdata$chr)
plotdata=plotdata[plotdata$chr != "PAR2",]
plotdata[plotdata$type=="GenesExclRepeats","type"]="Genes\nExcl\nRepeats"
plotdata[plotdata$type=="Simple/LowComplexity","type"]="Simple/\nLow\nComplexity"
plotdata[plotdata$type=="inter","type"]="100kb_bins"
head(plotdata)
table(plotdata$chr)
hist(plotdata$POS)
hist(plotdata$Meth)
table(plotdata$type)
table(plotdata$par)
plotdata$par
table(plotdata$species)
plotdata
head(plotdata)

plotdata=plotdata[plotdata$type!="100kb_bins",]
plotdata[plotdata$chr=="PAR1","chr"]="PAR1   "

plotdata$type=as.factor(plotdata$type)
levels(plotdata$type)
plotdata$type=relevel(plotdata$type,"Genes\nExcl\nRepeats")
plotdata$type=relevel(plotdata$type,"Simple/\nLow\nComplexity")
plotdata$type=relevel(plotdata$type,"Satellite")
plotdata$type=relevel(plotdata$type,"DNA")
plotdata$type=relevel(plotdata$type,"LTR")
plotdata$type=relevel(plotdata$type,"LINE")
plotdata$type=relevel(plotdata$type,"SINE")

plotdata$species=as.factor(plotdata$species)
plotdata$species=relevel(plotdata$species,"Siamang")
plotdata$species=relevel(plotdata$species,"S.orangutan")
plotdata$species=relevel(plotdata$species,"B.orangutan")
plotdata$species=relevel(plotdata$species,"Gorilla")
plotdata$species=relevel(plotdata$species,"Human")
plotdata$species=relevel(plotdata$species,"Chimpanzee")
plotdata$species=relevel(plotdata$species,"Bonobo")

min(plotdata$Meth,na.rm=T)
plotdata
head(plotdata)

windowsFonts()

p1a<-ggplot(plotdata, aes(x=chr, y=Meth)) + geom_boxplot(aes(fill = chr),lwd=0.25,outlier.shape=1,outlier.size=0.5)
p1a<-p1a + facet_grid(cols=vars(type),rows=vars(species), scale="free_x", drop=T, switch="y")
p1a<-p1a +  ylim(0,1.6) 
p1a<-p1a + theme_classic(base_size = 7,base_family="sans") 
p1a<-p1a + geom_signif(comparisons = list(c("PAR1   ", "X"), c("X", "Y"), c("PAR1   ", "Y")), size = 0.3, textsize = 5/.pt, step_increase = 0.2, na.rm = TRUE, test = "wilcox.test", test.args=list(exact = F), map_signif_level=c("***"=0.000001, "**"=0.001, "*"=0.05))
p1a<-p1a + scale_fill_manual(name=plotdata$chr, values=c("#66C2A5", "#FC8D62", "#8DA0CB"))  #c("green","red", "blue"))
p1a<-p1a + labs(title=NULL,x="", y = "Cytosine methylation (mC/C)")
p1a<-p1a + theme(legend.position = "none", panel.spacing = unit(0.1, "in", data = NULL), strip.background = element_blank(), strip.text.y = element_blank())
ggsave("repeats_plot.png",p1a,width=1650,height=1600,units="px",dpi=300)
#saveRDS(p1a, file = "p1a.rds")

p1a_stat<-ggplot(plotdata, aes(x=chr, y=Meth)) + geom_boxplot(aes(fill = chr),lwd=0.25,outlier.shape=1,outlier.size=0.5)
p1a_stat<-p1a_stat + facet_grid(cols=vars(type),rows=vars(species), scale="free_x", drop=T)
p1a_stat<-p1a_stat +  ylim(0,1.6) 
p1a_stat<-p1a_stat + theme_classic(base_size = 7,base_family="sans") 
p1a_stat<-p1a_stat + geom_signif(comparisons = list(c("PAR1   ", "X"), c("X", "Y"), c("PAR1   ", "Y")), size = 0.3, textsize = 5/.pt, step_increase = 0.2, na.rm = TRUE, test = "wilcox.test", test.args=list(exact = F))
p1a_stat<-p1a_stat + scale_fill_manual(name=plotdata$chr, values=c("#66C2A5", "#FC8D62", "#8DA0CB"))  #c("green","red", "blue"))
p1a_stat<-p1a_stat + labs(title=NULL,x="", y = "Cytosine methylation (mC/C)")
p1a_stat<-p1a_stat + theme(legend.position = "none", panel.spacing = unit(0.1, "in", data = NULL))
ggsave("repeats_plot_withstat.png",p1a_stat,width=2000,height=1600,units="px",dpi=300)


### Plot sequence classes ###
seqregionfiles=files[grep("seqregions",files)]
seqregionfiles
for (i in 1:length(seqregionfiles)) {
  #i=1
  temptab <- read.csv(seqregionfiles[i], stringsAsFactors = FALSE)
  head(temptab)
  if (i==1) {
    plotdata2=temptab
  }else {
    plotdata2=rbind(temptab,plotdata2)
  }
  numrows=length(row.names(plotdata2))
  cat(i,"_",numrows," ",sep="")
}
plotdata2=plotdata2[plotdata2$species != "CHM13",]
plotdata2[plotdata2$species=="HG002","species"]="Human"
plotdata2[plotdata2$species=="Bornean","species"]="B.orangutan"
plotdata2[plotdata2$species=="Sumatran","species"]="S.orangutan"
plotdata2=plotdata2[plotdata2$type != "PAR",]
plotdata2[plotdata2$type=="AMPLICONIC","type"]="Ampliconic"
plotdata2[(plotdata2$type=="X-ANCESTRAL/DEG"),"type"]="X-Ancestral"
head(plotdata2)
table(plotdata2$chr)
hist(plotdata2$POS)
hist(plotdata2$Meth)
table(plotdata2$type)
table(plotdata2$par)
plotdata2$par
table(plotdata2$species)
plotdata2=plotdata2[plotdata2$par!=2,]
plotdata2
head(plotdata2)

plotdata2[plotdata2$chr=="chrX","chr"]="\nchrX\n"
plotdata2[plotdata2$chr=="chrY","chr"]="\nchrY\n"
plotdata2[plotdata2$type=="Ampliconic","type"]="Amp."
plotdata2[plotdata2$type=="X-Ancestral","type"]="XAnc."

plotdata2$type=as.factor(plotdata2$type)
levels(plotdata2$type)
plotdata2$type=relevel(plotdata2$type,"XAnc.")
plotdata2$type=relevel(plotdata2$type,"Amp.")

plotdata2$species=as.factor(plotdata2$species)
plotdata2$species=relevel(plotdata2$species,"Siamang")
plotdata2$species=relevel(plotdata2$species,"S.orangutan")
plotdata2$species=relevel(plotdata2$species,"B.orangutan")
plotdata2$species=relevel(plotdata2$species,"Gorilla")
plotdata2$species=relevel(plotdata2$species,"Human")
plotdata2$species=relevel(plotdata2$species,"Chimpanzee")
plotdata2$species=relevel(plotdata2$species,"Bonobo")

min(plotdata2$Meth,na.rm=T)
plotdata2

p1b<-ggplot(plotdata2, aes(x=type, y=Meth)) + geom_boxplot(aes(fill = chr),lwd=0.25,outlier.shape=1,outlier.size=0.5)
p1b<-p1b + facet_grid(cols=vars(chr),rows=vars(species), scale="free_x", drop=T)
p1b<-p1b +  ylim(0.2,1.1)
p1b<-p1b + theme_classic(base_size = 7,base_family="sans")
p1b<-p1b + geom_signif(comparisons = list(c("Amp.", "XAnc.")), size = 0.3, textsize = 5/.pt, step_increase = 0.2, na.rm = TRUE, test = "wilcox.test", test.args=list(exact = F), map_signif_level=c("***"=0.000001, "**"=0.001, "*"=0.05))
p1b<-p1b + scale_fill_manual(name=plotdata2$chr, values=c("#FC8D62", "#8DA0CB"))
p1b<-p1b + labs(title=NULL,x="", y=NULL)
p1b<-p1b + theme(axis.text.x= element_text(size = 5), legend.position = "none", panel.spacing = unit(0.1, "in", data = NULL), strip.background = element_blank())
p1b_ancest<-p1b
ggsave("sequenceclasses_plot.png",p1b_ancest,width=450,height=1600,units="px",dpi=300)

p1b_stat<-ggplot(plotdata2, aes(x=type, y=Meth)) + geom_boxplot(aes(fill = chr),lwd=0.25,outlier.shape=1,outlier.size=0.5)
p1b_stat<-p1b_stat + facet_grid(cols=vars(chr),rows=vars(species), scale="free_x", drop=T)
p1b_stat<-p1b_stat +  ylim(0.2,1.1)
p1b_stat<-p1b_stat + theme_classic(base_size = 7,base_family="sans")
p1b_stat<-p1b_stat + geom_signif(comparisons = list(c("Amp.", "XAnc.")), size = 0.3, textsize = 5/.pt, step_increase = 0.2, na.rm = TRUE, test = "wilcox.test", test.args=list(exact = F))
p1b_stat<-p1b_stat + scale_fill_manual(name=plotdata2$chr, values=c("#FC8D62", "#8DA0CB"))
p1b_stat<-p1b_stat + labs(title=NULL,x="", y = "Cytosine methylation (mC/C)")
p1b_stat<-p1b_stat + theme(axis.text.x= element_text(size = 5), legend.position = "none", panel.spacing = unit(0.1, "in", data = NULL))
ggsave("sequenceclasses_plot_withstat.png",p1b_stat,width=450,height=1600,units="px",dpi=300)

##### 01End: Plot repeats and sequence classes #####


##### 02Start: Plot gene methylation profiles #####
l=ls()
l=ls()
l=l[which(!ls() %in% c("",""))]
rm(list=l)

library(ggplot2)
library(cowplot)

#colors=c("#CCFF00","#1e90ff", "#E31A1C", "#008b00", "#6A3D9A", "#FF7F00","#ff1493")

bonobo=read.table("Data/Methylation/gene_CATannot_methprofile_smoothed.Bonobo.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="")
bornean=read.table("Data/Methylation/gene_CATannot_methprofile_smoothed.Bornean.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="")
chimpanzee=read.table("Data/Methylation/gene_CATannot_methprofile_smoothed.Chimpanzee.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="")
gorilla=read.table("Data/Methylation/gene_CATannot_methprofile_smoothed.Gorilla.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="")
human=read.table("Data/Methylation/gene_CATannot_methprofile_smoothed.Human.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="")
siamang=read.table("Data/Methylation/gene_CATannot_methprofile_smoothed.Siamang.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="")
sumatran=read.table("Data/Methylation/gene_CATannot_methprofile_smoothed.Sumatran.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="")
head(bonobo)
sumatran[1:5,1:30]
sumatran[1:5,18000:18025]
bornean[1:5,1:30]
bornean[1:5,18000:18025]
bonobo[1:5,1:30]
bonobo[1:5,18000:18025]
siamang[1:5,1:30]
siamang[1:5,18000:18025]
chimpanzee[1:5,1:30]
chimpanzee[1:5,18000:18025]
gorilla[1:5,1:30]
gorilla[1:5,18000:18025]
human[1:5,1:30]
human[1:5,18000:18025]
dim(human[!((human$seqid=="chrY")&(human$end>=27449937)),])
human=human[!((human$seqid=="chrY")&(human$end>=27449937)),] #exclude Yq12 region
table(names(sumatran) %in% names(chimpanzee))
table(names(sumatran) %in% names(gorilla))
table(names(sumatran) %in% names(human))
table(names(sumatran) %in% names(bonobo))
table(names(sumatran) %in% names(bornean))
table(names(sumatran) %in% names(siamang))

commongenes=unique(human$gene_name)
commongenes=commongenes[commongenes %in% bonobo$gene_name]
commongenes=commongenes[commongenes %in% bornean$gene_name]
commongenes=commongenes[commongenes %in% chimpanzee$gene_name]
commongenes=commongenes[commongenes %in% gorilla$gene_name]
commongenes=commongenes[commongenes %in% siamang$gene_name]
commongenes=commongenes[commongenes %in% sumatran$gene_name]
commongenes=commongenes[order(commongenes)]
head(commongenes)

human$score="notcommongene"
human[human$gene_name %in% commongenes,"score"]="commongene"
human$score=paste(human$gene_biotype,human$score,human$seqid,sep="_")
bonobo$score="notcommongene"
bonobo[bonobo$gene_name %in% commongenes,"score"]="commongene"
bonobo$score=paste(bonobo$gene_biotype,bonobo$score,bonobo$seqid,sep="_")
bornean$score="notcommongene"
bornean[bornean$gene_name %in% commongenes,"score"]="commongene"
bornean$score=paste(bornean$gene_biotype,bornean$score,bornean$seqid,sep="_")
chimpanzee$score="notcommongene"
chimpanzee[chimpanzee$gene_name %in% commongenes,"score"]="commongene"
chimpanzee$score=paste(chimpanzee$gene_biotype,chimpanzee$score,chimpanzee$seqid,sep="_")
gorilla$score="notcommongene"
gorilla[gorilla$gene_name %in% commongenes,"score"]="commongene"
gorilla$score=paste(gorilla$gene_biotype,gorilla$score,gorilla$seqid,sep="_")
siamang$score="notcommongene"
siamang[siamang$gene_name %in% commongenes,"score"]="commongene"
siamang$score=paste(siamang$gene_biotype,siamang$score,siamang$seqid,sep="_")
sumatran$score="notcommongene"
sumatran[sumatran$gene_name %in% commongenes,"score"]="commongene"
sumatran$score=paste(sumatran$gene_biotype,sumatran$score,sumatran$seqid,sep="_")

table(human$score)
human[1:5,1:30]
human[,grep("pos_",names(human))]=human[,grep("pos_",names(human))]*100
sumatran[1:5,1:30]

#human_orig=human
#bonobo_orig=bonobo
#bornean_orig=bornean
#chimpanzee_orig=chimpanzee
#gorilla_orig=gorilla
#siamang_orig=siamang
#sumatran_orig=sumatran

human=human_orig
bonobo=bonobo_orig
bornean=bornean_orig
chimpanzee=chimpanzee_orig
gorilla=gorilla_orig
siamang=siamang_orig
sumatran=sumatran_orig

human=human[human$gene_biotype=="protein_coding",]
bonobo=bonobo[bonobo$gene_biotype=="protein_coding",]
bornean=bornean[bornean$gene_biotype=="protein_coding",]
chimpanzee=chimpanzee[chimpanzee$gene_biotype=="protein_coding",]
gorilla=gorilla[gorilla$gene_biotype=="protein_coding",]
siamang=siamang[siamang$gene_biotype=="protein_coding",]
sumatran=sumatran[sumatran$gene_biotype=="protein_coding",]

groups=c("protein_coding_commongene_chrX",
         "protein_coding_commongene_chrY")

png("genes_methylation_profile.png",width=700,height=700,units="px",res=300, pointsize=7)
par(mfrow = c(2, 1),oma=c(0,0,0,0),mar=c(1,4,0.5,0.5)+0.1)
colBorang= "#F9CF44"
colSorang= "#FD8A1A"
colBonobo= "#FB2D32"
colChimp= "#C3000E"
colSiamang= "#0E205E"
colGorilla= "#006E81"
colHuman= "#969696"
for (i in 1:2) {
  #i=1
  if (i==1) {
    label="chrX\nmethylation (mC/C)"
  }else if (i==2) {
    label="chrY\nmethylation (mC/C)"
  }
  plot(1,1,xlab="",ylab="",ylim=c(0.25,0.8),xlim=c(-1000,19000),xaxt="n",type="n",main="")
  rect(-4500,0.2,-2500,0.8,col="white",border=NA,xpd=T)
  text(rep(-3300,3),c(0.3,0.5,0.7),c("0.3","0.5","0.7"),xpd=T,cex=6/7,adj=c(0.5,0.5))
  mtext(label,2,1.8)
  abline(v=c(1,3001,15000,18000),lty=2,col="grey")
  text(c(1,3001,15000,18000),rep(0.18,4),c("TSS-3K","TSS","TES","TES+3k"),xpd=T,cex=6/7,family="sans")
  lines(1:18000,apply(human[human$score==groups[i],(length(names(human))-18000+1):length(names(human))],2,mean,na.rm=T)/100,pch=20,cex=0.2,col=colHuman,lwd=0.6)
  lines(1:18000,apply(chimpanzee[chimpanzee$score==groups[i],(length(names(chimpanzee))-18000+1):length(names(chimpanzee))],2,mean,na.rm=T)/100,pch=20,cex=0.2,col=colChimp,lwd=0.6)
  lines(1:18000,apply(bonobo[bonobo$score==groups[i],(length(names(bonobo))-18000+1):length(names(bonobo))],2,mean,na.rm=T)/100,pch=20,cex=0.2,col=colBonobo,lwd=0.6)
  lines(1:18000,apply(gorilla[gorilla$score==groups[i],(length(names(gorilla))-18000+1):length(names(gorilla))],2,mean,na.rm=T)/100,pch=20,cex=0.2,col=colGorilla,lwd=0.6)
  lines(1:18000,apply(bornean[bornean$score==groups[i],(length(names(bornean))-18000+1):length(names(bornean))],2,mean,na.rm=T)/100,pch=20,cex=0.2,col=colBorang,lwd=0.6)
  lines(1:18000,apply(sumatran[sumatran$score==groups[i],(length(names(sumatran))-18000+1):length(names(sumatran))],2,mean,na.rm=T)/100,pch=20,cex=0.2,col=colSorang,lwd=0.6)
  lines(1:18000,apply(siamang[siamang$score==groups[i],(length(names(siamang))-18000+1):length(names(siamang))],2,mean,na.rm=T)/100,pch=20,cex=0.2,col=colSiamang,lwd=0.6)
  cat(i," ",sep="")
}
par(mfrow = c(1, 1),oma=c(4,4,4,4),mar=c(4,4,4,4)+0.1)
graphics.off()

##### 02End: Plot gene methylation profiles #####













