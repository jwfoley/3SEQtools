# read gene count files into an integer matrix
read.gene.counts <- function(...) {
	count.vectors <- lapply(..., function(counts.file) {
		raw.data <- scan(counts.file, list(character(), integer()), comment.char = "#", quiet = T)
		count.vector <- raw.data[[2]]
		names(count.vector) <- raw.data[[1]]
		count.vector
	})
	if (length(count.vectors) > 1) for (i in 2:length(count.vectors)) stopifnot(all(names(count.vectors[[i]]) == names(count.vectors[[1]])))
	count.matrix <- do.call(cbind, count.vectors)
	colnames(count.matrix) <- sub("\\.tsv$", "", basename(...))
	count.matrix
}

