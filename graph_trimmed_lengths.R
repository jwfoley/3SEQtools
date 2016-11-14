#! /usr/bin/Rscript

# given a group of log files from trim_homopolymer.py (STDERR), make nice graphs of the insert length distributions

library(ggplot2)
library(scales)

UMI.LENGTH <- 5

plots <- lapply(commandArgs(trailingOnly = T), function(log.file) {
	lengths <- read.table(log.file, skip = 6 + UMI.LENGTH, col.names = c("length", "reads"))
	ggplot(lengths, aes(length, reads)) + 
		geom_area() +
		ggtitle(log.file) +
		scale_y_continuous(label = comma)
})

pdf("insert_lengths.pdf", width = 10, height = 7.5)
plots
dev.off()

