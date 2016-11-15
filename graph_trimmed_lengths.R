#! /usr/bin/Rscript

# given a group of log files from trim_homopolymer.py (STDERR), make nice graphs of the insert length distributions

library(ggplot2)
library(scales)

UMI.LENGTH <- 5

files <- commandArgs(trailingOnly = T)

lengths <- do.call(rbind, lapply(files, function(file) {
	result <- read.table(file, skip = 6 + UMI.LENGTH, col.names = c("length", "reads"))
	result$proportion <- result$reads / sum(result$reads)
	result$file <- file
	return(result)
}))


pdf("insert_lengths.pdf", width = 10, height = 7.5)
ggplot(lengths, aes(length, proportion)) +
		facet_wrap(~ file) + 
		geom_area() +
		xlab("insert length (nt)") +
		ylab("reads") +
		scale_y_continuous(label = percent)
dev.off()

