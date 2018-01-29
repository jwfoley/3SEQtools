#! /usr/bin/Rscript

# given a group of BAM files from 3SEQ, make an Excel-compatible file containing read counts, transcripts per million (TPM), and DESeq2-normalized (blind rlog) gene expression values
# assumes an Ensembl/GENCODE GTF annotation file

NO.RLOG.FLAG <- "--no-rlog" # use this flag to disable calculating rlog since it can take a long time
TMP.DIR <- "/tmp" # location for featureCounts temporary files


command.args <- commandArgs(trailingOnly = T)
if (length(command.args) < 3) stop(paste0("usage: make_expression_table.R [", NO.RLOG.FLAG, "] annotation.gtf library1.bam library2.bam ..."))
do.rlog <- ! "--no-rlog" %in% command.args
if (! do.rlog) command.args <- command.args[command.args != NO.RLOG.FLAG]

library(parallel)
library(WriteXLS)
library(Rsubread)
library(rtracklayer)

if (do.rlog) {
	library(BiocParallel)
	library(DESeq2)
}


annotation.file <- command.args[1]
bam.files <- command.args[-1]

library.names <- sub("\\.bam$", "", sapply(bam.files, basename))

feature.counts <- featureCounts(
	bam.files,
	annot.ext =            annotation.file,
	isGTFAnnotationFile =  T,
	strandSpecific =       1,
	read2pos =             5,
	nthreads =             detectCores(),
	tmpDir =               TMP.DIR
)

gene.counts <- feature.counts$counts

# read gene IDs + gene symbols from GTF file
annotations <- readGFF(annotation.file, filter = list(type = "gene"), tags = c("gene_id", "gene_name"))
gene.name.lookup <- annotations$gene_name
names(gene.name.lookup) <- annotations$gene_id
gene.id <- rownames(gene.counts)
gene.name <- gene.name.lookup[gene.id]


if (do.rlog) {
	colnames(gene.counts) <- NULL # rlogTransformation doesn't work if there are column names
	dds <- DESeq(DESeqDataSetFromMatrix(gene.counts, colData = data.frame(1:ncol(gene.counts)), ~ 1), fitType = "local", parallel = T, BPPARAM = MulticoreParam(workers = detectCores())) # eats a lot of memory!
	gene.rlogs <- assay(rlog(dds, blind = F)) # takes a long time! but slightly faster than rlog directly on the counts
	colnames(gene.rlogs) <- library.names
	rownames(gene.rlogs) <- rownames(gene.counts)
}

colnames(gene.counts) <- library.names

gene.tpm <- 1E6 * sweep(gene.counts, 2, colSums(gene.counts), "/")


save.image("make_expression_table.RData")

sheets <- list(
	data.frame("gene ID" = gene.id, "gene name" = gene.name, gene.counts, check.names = F),
	data.frame("gene ID" = gene.id, "gene name" = gene.name, round(gene.tpm, 1), check.names = F)
)
sheet.names <- c("count", "TPM")
if (do.rlog) {
	sheets[[3]] <- data.frame("gene ID" = gene.id, "gene name" = gene.name, round(gene.rlogs, 3), check.names = F)
	sheet.names[3] <- "rlog"
}

WriteXLS(
	sheets,
	ExcelFileName =  "gene_expression.xlsx",
	SheetNames =     sheet.names,
	AdjWidth =       T,
	BoldHeaderRow =  T,
	FreezeRow =      1,
	FreezeCol =      2
)

