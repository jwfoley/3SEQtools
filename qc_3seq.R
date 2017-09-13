#! /usr/bin/Rscript

library(reshape2)
library(scales)
library(ggplot2)


# define functions

filename.suffix <- list(
	trim =   ".trim.log",
	align =  ".align.log",
	dedup =  ".dedup.log"
)

category.colors <- c(
	"uniquely aligned" =  "green3",
	"multiply aligned" =  "darkgreen",
	"other" =             "yellow3",
	"too short" =         "yellow2",
	"RT dimer" =          "darkred",
	"PCR dimer" =         "red2",
	"duplicate" =         "orange",
	"non-duplicate" =     "blue"
)

get.filename <- function(library.name, suffix.name) paste0(library.name, filename.suffix[[suffix.name]])


parse.trim.log <- function(library.name) {
	trim.log <- file(get.filename(library.name, "trim"), "rt")
	total.reads <- scan(trim.log, integer(), 1, quiet = T)
	while (TRUE) {
		line <- readLines(trim.log, 1)
		if (length(line) == 0) stop() else if (line == "trimmed read length counts") break
	}
	length.counts <- read.table(trim.log, col.names = c("length", "count"))
	close(trim.log)
	c(
		"total" =      total.reads,
		"PCR dimer" =  sum(subset(length.counts, length < 0)$count),
		"RT dimer" =  subset(length.counts, length == 0)$count
	)
}

parse.align.log <- function(library.name) {
	log.list <- scan(get.filename(library.name, "align"), list(character(), character()), sep = "\t", strip.white = T, fill = T, quiet = T)
	log.vector <- log.list[[2]]
	names(log.vector) <- sub(" \\|$", "", log.list[[1]])
	input.reads <- as.integer(log.vector["Number of input reads"])
	c(
		"input" =             input.reads,
		"uniquely aligned" =  as.integer(log.vector["Uniquely mapped reads number"]),
		"multiply aligned" =  as.integer(log.vector["Number of reads mapped to too many loci"]),
		"too short" =         round(as.numeric(sub("%", "", log.vector["% of reads unmapped: too short"])) / 100 * input.reads),
		"other" =             round(as.numeric(sub("%", "", log.vector["% of reads unmapped: other"])) / 100 * input.reads)
	)
}

parse.dedup.log <- function(library.name) {
	log.list <- scan(get.filename(library.name, "dedup"), list(character(), character()), sep = "\t", strip.white = T, fill = T, quiet = T)
	log.vector <- as.integer(log.list[[1]])
	names(log.vector) <- log.list[[2]]
	result <- c(
		"duplicate" =      as.integer(log.vector["optical duplicates"] + log.vector["PCR duplicates"]),
		"non-duplicate" =  as.integer(log.vector["distinct alignments"] + log.vector["pre-PCR duplicates rescued by UMIs"] + log.vector["pre-PCR duplicates rescued by algorithm"]	)
	)
	stopifnot(sum(result) == log.vector["usable alignments read"])
	result
}


get.read.categories <- function(libraries) t(sapply(libraries, function(library.name) {
	trim.results <- parse.trim.log(library.name)
	align.results <- parse.align.log(library.name)
#	stopifnot(align.results[1] == trim.results[1] - sum(trim.results[-1]))
  c(trim.results, align.results[-1])
}))

get.dedup.counts <- function(libraries) t(sapply(libraries, parse.dedup.log))


plot.read.categories <- function(read.category.counts, normalize = FALSE) {
	result.frame <- melt(read.category.counts[,-1], varnames = c("library", "category"), value.name = "reads", as.is = T)
	result.frame$library <- factor(result.frame$library, levels = rownames(read.category.counts))
	result.frame$category <- factor(result.frame$category, levels = c("PCR dimer", "RT dimer", "too short", "other", "multiply aligned", "uniquely aligned"))
	if (normalize) {
		ggplot(result.frame) +
			geom_col(aes(library, reads, fill = category), position = "fill") +
			scale_y_continuous(label = percent) +
			theme(
				axis.text.x =       element_text(angle = 90, hjust = 1),
				panel.background =  element_blank()
			) +
			scale_fill_manual(values = category.colors)
	} else {
		ggplot(result.frame) +
			geom_col(aes(library, reads, fill = category)) +
			scale_y_continuous(label = comma) +
			theme(
				axis.text.x =         element_text(angle = 90, hjust = 1),
				panel.grid.major.x =  element_blank()
			) +
			scale_fill_manual(values = category.colors)
	}
}

plot.dedup <- function(dedup.counts) {
	result.frame <- melt(dedup.counts, varnames = c("library", "category"), value.name = "reads", as.is = T)
	result.frame$library <- factor(result.frame$library, levels = rownames(dedup.counts))
	result.frame$category <- factor(result.frame$category, levels = c("duplicate", "non-duplicate"))
	ggplot(result.frame) +
		geom_col(aes(library, reads, fill = category)) +
			scale_y_continuous(label = comma) +
			theme(
				axis.text.x =         element_text(angle = 90, hjust = 1),
				panel.grid.major.x =  element_blank()
			) +
			scale_fill_manual(values = category.colors)
}


# run script on libraries provided as command-line arguments
libraries <- commandArgs(trailingOnly = T)
if (length(libraries) > 0) {
	for (suffix in filename.suffix) libraries <- sub(suffix, "", libraries)
	cat("found libraries:\n")
	for (library in libraries) cat(library, "\n")

	read.category.counts <- get.read.categories(libraries)
	rownames(read.category.counts) <- basename(rownames(read.category.counts))
	read.category.count.plot <- plot.read.categories(read.category.counts)
	read.category.percent.plot <- plot.read.categories(read.category.counts, normalize = T)
	
	dedup.counts <- get.dedup.counts(libraries)
	dedup.count.plot <- plot.dedup(dedup.counts)

	save.image("qc_3seq.RData")
	write.table(read.category.counts, "read_category_count.tsv", quote = F, sep = "\t", col.names = NA)
	ggsave("read_category_count.pdf", read.category.count.plot, "pdf", width = 10, height = 7.5)
	ggsave("read_category_percent.pdf", read.category.percent.plot, "pdf", width = 10, height = 7.5)
	ggsave("dedup.pdf", dedup.count.plot, "pdf", width = 10, height = 7.5)
}

