# 3SEQtools

Pipelines for processing and QC of data from Smart-3SEQ or classic 3SEQ. Most scripts will show their usage options if you invoke them with the `--help` flag. Most scripts have hardcoded constants defined at the top so you can customize them to your needs.


## FASTQ generation pipeline

### `make_fastq.sh`

Wrapper to call `bcl2fastq` with preferred settings and simpler usage. Mostly used with NextSeq data so some features might not make sense on other platforms (e.g. `--no-lane-splitting` could be bad on HiSeq).


## FASTQ alteration pipeline

### `umi_homopolymer.py`

Removes the poly(A) tail, moves the UMI to the read metadata, and discards the G-overhang in every read of FASTQ data. Generates a simple log of trimmed read lengths. `align_smart-seq.sh` etc. invoke this but let the aligner do the poly(A) removal instead. Much faster with PyPy.


## Alignment pipelines

### `align_smart-3seq.sh`

Starts from a list of single-end Illumina FASTQ files and aligns them into sorted, indexed, duplicate-marked BAM files. `samtools` and `STAR` are expected to be in your `PATH` and [UMI-dedup](https://github.com/jwfoley/umi-dedup) locally installed. You must provide the path to a prepared STAR genome directory. `pypy` and `pigz` are preferred for performance but not required.

Options
	* `-n [length]`: extract leading UMI of specified length (default 5)
	* `-g [length]`: discard specified number of bases after UMI (default 3)
	* `-t [length]`: truncate reads to specified length
	* `-d`: skip duplicate marking

### `slurm_align_smart-3seq.sh`

Version of `align_smart-3seq.sh` to use on computing clusters powered by SLURM (submits a job array).

### `quasi-align_smart-3seq.sh` and `slurm_quasi-align_smart-3seq.sh`

Experimental pipeline using Salmon to quasi-align the data with special options that make sense for 3SEQ. This doesn't seem to get good results, though.


## Read-counting pipeline

### `make_expression_table.R`

Makes a table of gene expression, including read counts, TPM, and optionally DESeq2 rlog, from a list of BAM files. rlog may take a long time with large numbers of samples so this can be disabled with `--no-rlog`.


## Quality-control scripts

### `qc_3seq.R`

Makes graphs and tables from the output (log files) of `align_smart-3seq.sh` or equivalent.

### `graph_trimmed_lengths.R`

Makes graphs of trimmed insert lengths from the log files of `umi_homopolymer.py`.

### `aggregate_star_logs.py`

Collects the log files from STAR alignments into a simple spreadsheet.

### `trim_fastqc.sh`

Wrapper to run FastQC after trimming the data with `umi_homopolymer.py`.


## Miscellaneous

### `find_reference_homopolymers.py`

Uses a very (too?) simple algorithm to find sites in the reference sequence that might attract off-target priming with the oligo(dT) anchor.

### `blacklist_alignments.py`

Marks or removes reads from a BAM file that are inside unwanted regions in a BED file. Not sure if this is useful.
