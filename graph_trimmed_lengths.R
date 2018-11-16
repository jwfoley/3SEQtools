#! /usr/bin/Rscript

# given a group of log files from trim_homopolymer.py (STDERR), make nice graphs of the insert length distributions

library(ggplot2)
library(scales)

TABLE.HEADER <- "trimmed read length counts"
EXPECTED.FILE.SUFFIX <- "\\.trim\\.log$" # regular expression

files <- commandArgs(trailingOnly = T)

insert.lengths <- do.call(rbind, lapply(files, function(file) {
	raw.file <- readLines(file)
	result <- cbind(
		library = sub(EXPECTED.FILE.SUFFIX, "", basename(file)),
		read.table(textConnection(raw.file), skip = which(raw.file == TABLE.HEADER), col.names = c("length", "reads"))
	)
	result$proportion <- result$reads / sum(result$reads)
	return(result)
}))

insert.length.plot <- ggplot(insert.lengths, aes(length, proportion)) +
		facet_wrap(~ library) + 
		geom_area() +
		xlab("insert length (nt)") +
		ylab("reads") +
		scale_y_continuous(label = percent)

save.image("insert_lengths.RData")
ggsave("insert_lengths.pdf", insert.length.plot, "pdf", width = 10, height = 7.5)

