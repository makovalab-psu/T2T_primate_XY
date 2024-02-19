#! /usr/bin/env Rscript

###############################################################################
#
#  Except as indicated otherwise, this file is a 'United States Government
#  Work', and is released in the public domain.
#
#  File 'LICENSE' in this directory (xtr_search) of this distribution
#  contains full conditions and disclaimers.
##

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) {
	stop("Expected 2 args: (1) input.dots (2) output basename. Note, the input.dots file is expected to have the headerline with TOLID1-XorY and TOLID2-XorY, e.g., HG002-chrY	mPanPan1-chrY.")
}
input_dots_file <- args[1]
output_base <- args[2]
output_pdf <- paste(output_base, "pdf", sep='.')
output_png <- paste(output_base, "png", sep='.')

d <- read.csv(input_dots_file, sep='\t', header=TRUE, na.strings=c("NA"))

tolid1 <- colnames(d)[1]
tolid2 <- colnames(d)[2]
chrom1 <- strsplit(tolid1, '\\.')[[1]][2]
chrom2 <- strsplit(tolid2, '\\.')[[1]][2]
tolid1 <- strsplit(tolid1, '\\.')[[1]][1]
tolid2 <- strsplit(tolid2, '\\.')[[1]][1]

colnames(d) <- c("first", "second")

d$first <- d$first/1000000
d$second <- d$second/1000000

p <- ggplot(data=d) +
		annotate("rect", xmin=-Inf,      ymin=-Inf, xmax=2.45832 ,  ymax=Inf, linewidth=0, alpha=1,    fill="#97CB99") + # PAR
		annotate("rect", xmin=2.45832,   ymin=-Inf, xmax=2.727072,  ymax=Inf, linewidth=0, alpha=1,    fill="#FFEF57") + # X-DEG
		annotate("rect", xmin=2.727072,  ymin=-Inf, xmax=5.914561,  ymax=Inf, linewidth=0, alpha=1,    fill="#EEA9BA") + # XTR1
		annotate("rect", xmin=5.914561,  ymin=-Inf, xmax=6.200973,  ymax=Inf, linewidth=0, alpha=1,    fill="#88C0EA") + # AMPL
		annotate("rect", xmin=6.200973,  ymin=-Inf, xmax=6.400875,  ymax=Inf, linewidth=0, alpha=1,    fill="#EEA9BA") + # XTR2
		annotate("rect", xmin=6.400875,  ymin=-Inf, xmax=7.234999,  ymax=Inf, linewidth=0, alpha=1,    fill="#FFEF57") + # X-DEG
		annotate("rect", xmin=7.234999,  ymin=-Inf, xmax=10.389946, ymax=Inf, linewidth=0, alpha=1,    fill="#88C0EA") + # AMPL
		annotate("rect", xmin=10.389946, ymin=-Inf, xmax=10.565750, ymax=Inf, linewidth=0, alpha=1,    fill="#6A4700") + # SAT
		annotate("rect", xmin=10.565750, ymin=-Inf, xmax=10.883085, ymax=Inf, linewidth=0, alpha=1,    fill="#B02026") + # CEN (asat/HOR)
		annotate("rect", xmin=10.883085, ymin=-Inf, xmax=11.711672, ymax=Inf, linewidth=0, alpha=1,    fill="#6A4700") + # SAT
		annotate("rect", xmin=11.711672, ymin=-Inf, xmax=12.656463, ymax=Inf, linewidth=0, alpha=1,    fill="#6A4700") + # SAT
		annotate("rect", xmin=12.656463, ymin=-Inf, xmax=14.891106, ymax=Inf, linewidth=0, alpha=1,    fill="#FFEF57") + # X-DEG
		annotate("rect", xmin=14.891106, ymin=-Inf, xmax=14.964927, ymax=Inf, linewidth=0, alpha=1,    fill="#88C0EA") + # AMPL
		annotate("rect", xmin=14.964927, ymin=-Inf, xmax=16.781335, ymax=Inf, linewidth=0, alpha=1,    fill="#FFEF57") + # X-DEG
		annotate("rect", xmin=16.781335, ymin=-Inf, xmax=16.811393, ymax=Inf, linewidth=0, alpha=1,    fill="#88C0EA") + # AMPL
		annotate("rect", xmin=16.811393, ymin=-Inf, xmax=17.066094, ymax=Inf, linewidth=0, alpha=1,    fill="#FFEF57") + # X-DEG
		annotate("rect", xmin=17.066094, ymin=-Inf, xmax=17.332346, ymax=Inf, linewidth=0, alpha=1,    fill="#88C0EA") + # AMPL
		annotate("rect", xmin=17.332346, ymin=-Inf, xmax=18.362331, ymax=Inf, linewidth=0, alpha=1,    fill="#FFEF57") + # X-DEG
		annotate("rect", xmin=18.362331, ymin=-Inf, xmax=19.776671, ymax=Inf, linewidth=0, alpha=1,    fill="#88C0EA") + # AMPL
		annotate("rect", xmin=19.776671, ymin=-Inf, xmax=20.961203, ymax=Inf, linewidth=0, alpha=1,    fill="#FFEF57") + # X-DEG
		annotate("rect", xmin=20.961203, ymin=-Inf, xmax=21.226263, ymax=Inf, linewidth=0, alpha=1,    fill="#6A4700") + # SAT
		annotate("rect", xmin=21.226263, ymin=-Inf, xmax=22.270252, ymax=Inf, linewidth=0, alpha=1,    fill="#FFEF57") + # X-DEG
		annotate("rect", xmin=22.270252, ymin=-Inf, xmax=27.124000, ymax=Inf, linewidth=0, alpha=1,    fill="#88C0EA") + # AMPL
		annotate("rect", xmin=27.124000, ymin=-Inf, xmax=27.449931, ymax=Inf, linewidth=0, alpha=1,    fill="#D8D8D8") + # OTHER
		annotate("rect", xmin=27.449931, ymin=-Inf, xmax=62.072743, ymax=Inf, linewidth=0, alpha=1,    fill="#777777") + # HET
		annotate("rect", xmin=62.072743, ymin=-Inf, xmax=62.122809, ymax=Inf, linewidth=0, alpha=1,    fill="#D8D8D8") + # OTHER
		annotate("rect", xmin=62.122809, ymin=-Inf, xmax=62.460029, ymax=Inf, linewidth=0, alpha=1,    fill="#97CB99") + # PAR
		geom_path(data=d, mapping=aes(x=first, y=second), na.rm=TRUE, show.legend=FALSE) +
		ggtitle("Dotplot of LASTZ alignments") +
		xlab(paste(tolid1, " chr", chrom1, " (Mb)", sep='')) +
		ylab(paste(tolid2, " chr", chrom2, "(Mb)", sep='')) +
		theme_classic() +
		scale_x_continuous(expand=c(0,0)) +
		scale_y_continuous(expand=c(0,0)) +
		theme(axis.title=element_text(size=16, family="serif", face="bold"), axis.text=element_text(size=14, family="serif", color="black"),
				axis.ticks=element_line(color="black"), plot.title=element_text(size=20, family="serif", face="bold", hjust=0.5)
			)

ggsave(file=output_pdf, plot=p, device="pdf", units="in", height=8, width=8, limitsize=FALSE)
ggsave(file=output_png, plot=p, device="png", units="in", height=8, width=8, limitsize=FALSE, dpi=300)

# quit
q()

