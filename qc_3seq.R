#! /usr/bin/Rscript

library(reshape2)
library(scales)
library(ggplot2)

filename.suffix <- list(
	trim =   "_trim.log",
	align =  "_align.log",
)

category.colors <- c(
	"uniquely aligned" =  "green",
	"multiply aligned" =  "darkgreen",
	"other" =             "yellow3",
	"too short" =         "yellow",
	"no insert" =         "darkred",
	"no header" =         "red"
)

libraries <- sub(filename.suffix$trim, "", Sys.glob(paste0("*", filename.suffix$trim))) # test

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
		"no header" =  sum(subset(length.counts, length < 0)$count),
		"no insert" =  subset(length.counts, length == 0)$count
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

get.read.categories <- function(libraries) t(sapply(libraries, function(library.name) {
	trim.results <- parse.trim.log(library.name)
	align.results <- parse.align.log(library.name)
	stopifnot(align.results[1] == trim.results[1] - sum(trim.results[-1]))
  c(trim.results, align.results[-1])
}))

plot.read.categories <- function(read.category.counts, normalize = FALSE) {
	result.frame <- melt(read.category.counts[,-1], varnames = c("library", "category"), value.name = "reads", as.is = T)
	result.frame$category <- factor(result.frame$category, levels = c("no header", "no insert", "too short", "other", "multiply aligned", "uniquely aligned"))
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

