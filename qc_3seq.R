#! /usr/bin/Rscript

library(reshape2)
library(scales)
library(ggplot2)
library(optparse)

default.dims <- "7.5,10"
graph.theme <- theme(
	axis.text.y =         element_text(hjust = 0),
	axis.ticks.y =        element_blank(),
	panel.grid.major.y =  element_blank()
)

# parse command-line arguments
opt <- parse_args(OptionParser(option_list = list(
	make_option(c("-d", "--dim"), action = "store", default = default.dims, help = paste0("dimensions of output graphs in form W,H (default ", default.dims, ")"))
)), positional_arguments = c(1, Inf))
graph.dims <- as.list(as.numeric(strsplit(opt$options$dim, ",")[[1]]))
names(graph.dims) <- c("width", "height")

# define functions

filename.suffix <- list(
	trim =      ".trim.log",
	align =     ".align.log",
	dedup =     ".dedup.log",
	category =  ".category.log",
	bam =       ".bam",
	bai =       ".bai"
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
	log.file <- get.filename(library.name, "dedup")
	if (! file.exists(log.file)) return(NULL)
	log.list <- scan(log.file, list(character(), character()), sep = "\t", strip.white = T, fill = T, quiet = T)
	log.vector <- log.list[[1]]
	names(log.vector) <- log.list[[2]]
	result <- c(
		"duplicate" =       as.integer(log.vector["optical duplicates"]) + as.integer(log.vector["PCR duplicates"]),
		"non-duplicate" =   as.integer(log.vector["distinct alignments"]) + as.integer(log.vector["pre-PCR duplicates rescued by UMIs"]) + as.integer(log.vector["pre-PCR duplicates rescued by algorithm"]),
		"entropy.before" =  as.numeric(log.vector["library entropy before deduplication"]),
		"entropy.after" =   as.numeric(log.vector["library entropy after deduplication"]),
		"library.size" =    as.integer(log.vector["estimated library size"])
	)
	stopifnot(sum(result["duplicate"], result["non-duplicate"]) == as.integer(log.vector["usable alignments read"]))
	result
}

parse.alignment.categories <- function(library.name) {
	log.file <- get.filename(library.name, "category")
	if (! file.exists(log.file)) return(NULL)
	log.list <- scan(log.file, list(character(), character()), sep = "\t", quote = "", strip.white = T, fill = T, quiet = T)
	log.vector <- as.integer(log.list[[2]])
	names(log.vector) <- log.list[[1]]
	stopifnot(sum(log.vector[-1]) == log.vector[1])
	log.vector
}


get.read.categories <- function(libraries) t(sapply(libraries, function(library.name) {
	trim.results <- parse.trim.log(library.name)
	align.results <- parse.align.log(library.name)
#	stopifnot(align.results[1] == trim.results[1] - sum(trim.results[-1]))
  c(trim.results, align.results[-1])
}))

get.dedup.results <- function(libraries) t(sapply(libraries, parse.dedup.log))

get.alignment.categories <- function(libraries) t(sapply(libraries, parse.alignment.categories))


plot.read.categories <- function(read.category.counts, normalize = FALSE) {
	result.frame <- melt(read.category.counts[,-1], varnames = c("library", "category"), value.name = "reads", as.is = T)
	result.frame$library <- factor(result.frame$library, levels = rev(rownames(read.category.counts))) # reversed for flipped coordinates
	result.frame$category <- factor(result.frame$category, levels = c("PCR dimer", "RT dimer", "too short", "other", "multiply aligned", "uniquely aligned"))
	if (normalize) {
		ggplot(result.frame) +
			geom_col(aes(library, reads, fill = category), position = "fill", width = 1) +
			scale_y_continuous(label = percent, expand = c(0, 0)) +
			coord_flip() +
			scale_fill_manual(values = category.colors) +
			graph.theme +
			theme(panel.background = element_blank())
	} else {
		ggplot(result.frame) +
			geom_col(aes(library, reads, fill = category), width = 1) +
			scale_y_continuous(label = comma, expand = c(0, 0)) +
			coord_flip() +
			scale_fill_manual(values = category.colors) +
			graph.theme
	}
}

plot.dedup <- function(dedup.results) {
	count.frame <- melt(dedup.results[,c("duplicate", "non-duplicate")], varnames = c("library", "category"), value.name = "reads", as.is = T)
	count.frame$library <- factor(count.frame$library, levels = rev(rownames(dedup.results))) # reversed for flipped coordinates
	count.frame$category <- factor(count.frame$category, levels = c("duplicate", "non-duplicate"))
	
	estimate.frame <- data.frame(
		x = 1:nrow(dedup.results) - 1/2,
		xend = 1:nrow(dedup.results) + 1/2,
		y = rev(dedup.results[,"library.size"]) # reversed because plotting coordinate is reversed
	)
	
	ggplot(count.frame) +
		geom_col(aes(library, reads, fill = category), width = 1) +
		scale_y_continuous(label = comma, expand = c(0, 0)) +
		scale_fill_manual(values = category.colors) +
		coord_flip() +
		graph.theme +
		ylab("uniquely aligned reads") + # y because flipped coordinates
		geom_segment(aes(x = x, xend = xend, y = y, yend = y), data = estimate.frame)
}

plot.entropy <- function(dedup.results) {
	entropy.frame <- melt(dedup.results[,c("entropy.before", "entropy.after")], varnames = c("library", "category"), value.name = "entropy", as.is = T)
	entropy.frame$library <- factor(entropy.frame$library, levels = rev(rownames(dedup.results))) # reversed for flipped coordinates
	
	ggplot(entropy.frame) +
		geom_point(aes(library, entropy, color = category)) +
		coord_flip() +
		graph.theme +
		ylab("entropy of reads per position") # y because flipped coordinates
}

plot.alignment.categories <- function(alignment.categories, normalize = TRUE) {
	result.frame <- melt(alignment.categories[,-1], varnames = c("library", "category"), value.name = "alignments", as.is = T)
	result.frame$library <- factor(result.frame$library, levels = rev(rownames(alignment.categories))) # reversed for flipped coordinates
	result.frame$category <- factor(result.frame$category, levels = unique(result.frame$category))
	if (normalize) {
		ggplot(result.frame) +
			geom_col(aes(library, alignments, fill = category), position = "fill", width = 1) +
			scale_y_continuous(label = percent, expand = c(0, 0)) +
			coord_flip() +	
			graph.theme +
			theme(panel.background = element_blank())
	} else {
		ggplot(result.frame) +
			geom_col(aes(library, alignments, fill = category), width = 1) +
			scale_y_continuous(label = comma, expand = c(0, 0)) +
			coord_flip() +	
			graph.theme
	}
}


# run script on libraries provided as command-line arguments
libraries <- opt$arg
if (length(libraries) > 0) {
	libraries <- sub(do.call(paste, c(as.list(gsub("\\.", "\\\\.", paste0(filename.suffix, "$"))), sep = "|")), "", libraries) # convert suffixes to a regex before removing them
	cat("found libraries:\n")
	for (library in libraries) cat(library, "\n")

	read.category.counts <- get.read.categories(libraries)
	rownames(read.category.counts) <- basename(rownames(read.category.counts))
	read.category.count.plot <- plot.read.categories(read.category.counts)
	read.category.percent.plot <- plot.read.categories(read.category.counts, normalize = T)
	
	dedup.results <- get.dedup.results(libraries)
	have.dedup.results <- any(! sapply(dedup.results, is.null))
	if (have.dedup.results) {
		rownames(dedup.results) <- basename(rownames(dedup.results))
		dedup.count.plot <- plot.dedup(dedup.results)
		dedup.entropy.plot <- plot.entropy(dedup.results)
	}
	
	alignment.category.counts <- get.alignment.categories(libraries)
	have.alignment.categories <- any(! sapply(alignment.category.counts, is.null))
	if (have.alignment.categories) {
		rownames(alignment.category.counts) <- basename(rownames(alignment.category.counts))
		alignment.category.plot <- plot.alignment.categories(alignment.category.counts)
	}

	save.image("qc_3seq.RData")
	write.table(read.category.counts, "read_category_count.tsv", quote = F, sep = "\t", col.names = NA)
	ggsave("read_category_count.pdf", read.category.count.plot, "pdf", width = graph.dims$width, height = graph.dims$height)
	ggsave("read_category_percent.pdf", read.category.percent.plot, "pdf", width = graph.dims$width, height = graph.dims$height)
	if (have.dedup.results) {
		ggsave("dedup.pdf", dedup.count.plot, "pdf", width = graph.dims$width, height = graph.dims$height)
		ggsave("entropy.pdf", dedup.entropy.plot, "pdf", width = graph.dims$width, height = graph.dims$height)
	}
	if (have.alignment.categories) ggsave("alignment_categories.pdf", alignment.category.plot, "pdf", width = graph.dims$width, height = graph.dims$height)
}

