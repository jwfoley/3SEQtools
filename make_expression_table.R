#! /usr/bin/Rscript

# given a group of BAM files from 3SEQ, make an Excel-compatible file containing read counts, transcripts per million (TPM), and DESeq2-normalized (blind rlog) gene expression values
# assumes an Ensembl/GENCODE GTF annotation file

command.args <- commandArgs(trailingOnly = T)
if (length(command.args) < 3) stop("usage: make_expression_table.R annotation.gtf library1.bam library2.bam ...")

library(parallel)
library(WriteXLS)
library(Rsubread)
library(biomaRt)
library(DESeq2)

ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org")


annotation.file <- command.args[1]
bam.files <- command.args[-1]

library.names <- sub("\\.bam$", "", sapply(bam.files, basename))

feature.counts <- featureCounts(
	bam.files,
	annot.ext =            annotation.file,
	isGTFAnnotationFile =  T,
	strandSpecific =       1,
	nthreads =             detectCores()
)

gene.counts <- feature.counts$counts
colnames(gene.counts) <- NULL # rlogTransformation doesn't work if there are column names

ensembl.gene.id <- sub("\\..*", "", rownames(gene.counts))
gene.name.table <- getBM(
	attributes =  c("ensembl_gene_id", "external_gene_name"),
	filters =     "ensembl_gene_id",
	values =      ensembl.gene.id,
	mart =        ensembl
)
rownames(gene.name.table) <- gene.name.table$ensembl_gene_id
gene.name = gene.name.table[ensembl.gene.id, "external_gene_name"]

gene.rlogs <- rlog(gene.counts, fitType = "local")
colnames(gene.counts) <- library.names
colnames(gene.rlogs) <- library.names
rownames(gene.rlogs) <- rownames(gene.counts)

gene.tpm <- t(apply(gene.counts, 1, function(x) x / colSums(gene.counts) * 1E6))


save.image("make_expression_table.RData")

WriteXLS(
	list(
		data.frame("Ensembl ID" = ensembl.gene.id, "gene name" = gene.name, gene.counts, check.names = F),
		data.frame("Ensembl ID" = ensembl.gene.id, "gene name" = gene.name, round(gene.tpm, 1), check.names = F),
		data.frame("Ensembl ID" = ensembl.gene.id, "gene name" = gene.name, round(gene.rlogs, 3), check.names = F)
	),
	ExcelFileName =  "gene_expression.xlsx",
	SheetNames =     c("count", "TPM", "rlog"),
	AdjWidth =       T,
	BoldHeaderRow =  T,
	FreezeRow =      1,
	FreezeCol =      2
)

