#! /bin/bash

# given one or more .fastq.gz files, run FastQC on their contents after processing them with umi_homopolymer.py, but avoid storing the trimmed data more than necessary
# unfortunately FastQC requires seekable input so temporary files must still be created; I am seeking (!) a better workaround

tmp_dir=/tmp
unzip_cmd='pigz -dc'
zip_cmd='pigz --fast' # command to compress the temporary trimmed FASTQ files, e.g. 'pigz --fast', or 'cat' to leave them uncompressed
tmp_suffix='.gz' # suffix for temporary files for FastQC to recognize the format, e.g. '.gz' if zip_cmd involves gzip, or '' if uncompressed
fastqc_cmd='fastqc --noextract --nogroup --quiet --outdir .'
trim_cmd="pypy $(dirname $0)/umi_homopolymer.py" # settings for umi_homopolymer.py go in here

parallel -k "
	name=\$(basename {} .fastq.gz)
	tmp_file=$tmp_dir/\$name$tmp_suffix
	echo Trimming results for \$name: >&2
	$unzip_cmd {} | $trim_cmd | $zip_cmd > \$tmp_file
	$fastqc_cmd \$tmp_file
	rm \$tmp_file
	echo >&2
" ::: "$@"

