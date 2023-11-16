#! /usr/bin/env Rscript

getParallelVecViaDict <- function(vec, dict) {
	vec2 <- c()
	for(item in vec) {
		vec2 <- c(vec2, dict[[item]])
	}
	return(vec2)
}

library(ggplot2)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 6) {
	stop("Expected 6 args: (1) colors.tsv (2) seq-classes.tsv (3) input.tsv (4) output basename (5) ToLID (6) chrY region of interest or 'null'")
}

colors_tsv_file <- args[1]
#seq.class	hex.color	rgb.color
#PAR	#97CB99	151,203,153
#XDEG	#FFEF57	255,239,87
#XTR	#EEA9BA	238,169,186
#AMPLICONIC	#88C0EA	136,192,234
#aSat(HOR)	#B02026	176,32,38
#SATELLITE	#6A4700	106,71,0
#HET	#777777	119,119,119
#OTHER	#D8D8D8	217,216,216

seqclasses_tsv_file <- args[2]
# BED file with sequence classes. 3 example lines:
#chrY    0   37079   OTHER
#chrY    37079   2526344 PAR
#chrY    2526344 4150212 XDEG

input_tsv_file <- args[3]
# this is lastz's output after being passed through this awk command (TOLID= mPanPan1, mPanPan3, etc.):
#awk 'BEGIN{FS=" +"; OFS="\t"}!/^#/{print $6, $7, $2, $3, gensub(/%$/, "", 1, $9)}' < ${TOLID}.chrY_onto_${TOLID}.chrX.lz > ${TOLID}.chrY_onto_${TOLID}.chrX.forR.tsv
# assuming the lastz output has these columns to begin with:
#name2  zstart2+  end2+     name1  strand2  zstart1   end1      nmatch  id%     cigarx
# thus the resulting TSV has these columns (though without the '%' character at the end of every line):
#zstart1	end1	zstart2+	end2+	id%

output_base <- args[4]
# some output prefix for the filenames. The output of the script will be 3
# plots, each with a PNG and PDF copy. The 3rd plot will be a combination of
# the other two as separate panels, and it will be named as the output_base
# with '.pdf' and '.png' appended to the end. The two sub-plots are named
# similarly, except that they have '.dot' or '.pid' before the final suffix.
# The dot file is the dotplot. The pid file is the percent identity.

tolid <- args[5]
# This is the ToLID (e.g., mPanPan1) for the sample being run. This is just
# used for naming axis labels.

roi <- args[6]
# this allows you to specify a portion of chrY to zoom in on. This was mainly
# for testing/exploration purposes. For the final images for the manuscript,
# "null" was supplied to plot the entire chromosome instead of only a subset.
# If a "roi" (region of interest) were to be supplied, use the form "#-#"
# (i.e., as a range, start through end). Plots may not work well when not
# supploying "null".

output_pdf <- paste(output_base, "pdf", sep='.')
output_png <- paste(output_base, "png", sep='.')
output_dot_pdf <- paste(output_base, "dot", "pdf", sep='.')
output_dot_png <- paste(output_base, "dot", "png", sep='.')
output_pid_pdf <- paste(output_base, "pid", "pdf", sep='.')
output_pid_png <- paste(output_base, "pid", "png", sep='.')

roi.start <- 0
xlimmin <- NA
roi.end   <- Inf
xlimmax <- NA
if (roi != "null") {
	roi.start <- strsplit(roi, '-')[[1]][1]
	roi.end   <- strsplit(roi, '-')[[1]][2]
	if (roi.start == "null") {
		roi.start <- 0
	} else {
		roi.start <- as.numeric(roi.start)/1000000
		xlimmin <- roi.start
	}
	if (roi.end == "null") {
		roi.end <- Inf
	} else {
		roi.end <- as.numeric(roi.end)/1000000
		xlimmax <- roi.end
	}
	if(!is.na(xlimmin) & !is.na(xlimmax)) {
		expansion <- (xlimmax-xlimmin) * 0.05
		if(expansion < 0.001) { expansion <- 0.001 }
		xlimmax <- xlimmax + expansion
		xlimmin <- xlimmin - expansion
		if(xlimmin < 0){xlimmin <- 0}
	}
}

c <- read.csv(colors_tsv_file, sep='\t', header=TRUE, na.strings=c("NA"))
# expecting the following header:
#colnames(c) <- c("seq.class", "hex.color", "rgb.header")

class2hex <- c$hex.color
names(class2hex) <- c$seq.class

s <- read.csv(seqclasses_tsv_file, sep='\t', header=FALSE, na.strings=c("NA"))
# not expecting any header
colnames(s) <- c("chromosome", "start", "end", "seq.class")
s$start <- s$start/1000000
s$end <- s$end/1000000
s$hex.color <- getParallelVecViaDict(s$seq.class, class2hex)
s$seq.class <- as.factor(s$seq.class)

d <- read.csv(input_tsv_file, sep='\t', header=FALSE, na.strings=c("NA"))
# not expecting any header
colnames(d) <- c("qry.start", "qry.end", "ref.start", "ref.end", "pid")

d$qry.start <- d$qry.start/1000000
d$qry.end <- d$qry.end/1000000
d$ref.start <- d$ref.start/1000000
d$ref.end <- d$ref.end/1000000

# subset d based on provided ROI
d <- d[!(d$qry.start > roi.end | d$qry.end < roi.start),]

# subset s based on provided ROI's effect on d
rstart.adj <- min(d$qry.start)
rend.adj <- max(d$qry.end)
s <- s[!(s$start > rend.adj | s$end < rstart.adj),]
s$start[s$start < rstart.adj & s$end > rstart.adj] <- rstart.adj
s$end[s$end > rend.adj & s$start < rend.adj] <- rend.adj
s[1, "start"] <- -Inf
s[length(s$end), "end"] <- Inf

chrom1 <- "Y"
chrom2 <- "X"
tolid1 <- tolid
tolid2 <- tolid

p1 <- ggplot() +
		geom_rect(data=s, mapping=aes(xmin=start, ymin=-Inf, xmax=end, ymax=Inf, fill=seq.class), linewidth=0, alpha=1, show.legend=FALSE) +
		scale_fill_manual(values=class2hex) +
		geom_segment(data=d,  mapping=aes(x=qry.start, xend=qry.end, y=ref.start, yend=ref.end), linewidth=1, na.rm=TRUE, show.legend=FALSE) +
		xlab(paste(tolid1, " chr", chrom1, " (Mb)", sep='')) +
		ylab(paste(tolid2, " chr", chrom2, " (Mb)", sep='')) +
		theme_classic() +
		scale_x_continuous(expand=c(0,0), limits=c(xlimmin, xlimmax)) +
		scale_y_continuous(expand=c(0,0)) +
		theme(axis.title=element_text(size=16, family="serif", face="bold"), axis.text=element_text(size=14, family="serif", color="black"),
				axis.ticks=element_line(color="black")
			)

ggsave(file=output_dot_pdf, plot=p1, device="pdf", units="in", height=8,  width=8, limitsize=FALSE)
ggsave(file=output_dot_png, plot=p1, device="png", units="in", height=8,  width=8, limitsize=FALSE, dpi=300)

p2 <- ggplot() +
		geom_rect(data=s, mapping=aes(xmin=start, ymin=-Inf, xmax=end, ymax=Inf, fill=seq.class), linewidth=0, alpha=1, show.legend=FALSE) +
		scale_fill_manual(values=class2hex) +
		geom_segment(data=d, mapping=aes(x=qry.start, xend=qry.end, y=pid, yend=pid), linewidth=1, na.rm=TRUE, show.legend=FALSE) +
		xlab(paste(tolid1, " chr", chrom1, " (Mb)", sep='')) +
		ylab("% Identity") +
		theme_classic() +
		scale_x_continuous(expand=c(0,0), limits=c(xlimmin, xlimmax)) +
		scale_y_continuous(expand=c(0,0), limits=c(85,100)) +
		theme(axis.title=element_text(size=16, family="serif", face="bold"), axis.text=element_text(size=14, family="serif", color="black"),
				axis.ticks=element_line(color="black")
			)

ggsave(file=output_pid_pdf, plot=p2, device="pdf", units="in", height=2,  width=8, limitsize=FALSE)
ggsave(file=output_pid_png, plot=p2, device="png", units="in", height=2,  width=8, limitsize=FALSE, dpi=300)

p1alt <- p1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())
p <- plot_grid(p1alt, p2, align="v", axis="l", nrow=2, ncol=1, rel_heights=c(4,1))

pdf(NULL)
save_plot(filename=output_pdf, plot=p, ncol=1, nrow=2, base_height=4.5, base_width=8, device="pdf", limitsize=FALSE)
save_plot(filename=output_png, plot=p, ncol=1, nrow=2, base_height=4.5, base_width=8, device="png", limitsize=FALSE, dpi=300)

# quit
q()

