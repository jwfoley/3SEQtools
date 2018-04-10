#! /bin/bash

# given the path to an Illumina run folder and optionally a path to a sample sheet, run bcl2fastq
# generates FASTQ files in the current directory and cleans up their names
# the sample sheet argument is necessary unless it's the default $run_folder/SampleSheet.csv
# discards "Undetermined" reads by default unless "-u" is used
# adds an option to ignore read 2 (index read 1) if using Smart-3SEQ 96-plex i5 indexing
# passes any additional arguments to bcl2fastq

bcl2fastq_cmd='bcl2fastq --fastq-compression-level 9 --no-lane-splitting --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0'
undetermined_filename='Undetermined_S0_R1_001.fastq.gz'
default_sample_sheet='SampleSheet.csv'
i5_only_regex='Smart-3SEQ 96-plex i5 indexing'
i5_only_option='--use-bases-mask y*,n*,i*'


sample_sheet=''
discard_undetermined=true
while getopts ":s:u" opt
do
	case $opt in
		s)
			if [ ! -e "$OPTARG" ]; then echo "error: $OPTARG not found" >&2; exit 1; fi
			sample_sheet_arg="--sample-sheet $OPTARG"
			;;
		u)
			discard_undetermined=false
			;;
	esac
done
shift "$((OPTIND-1))"

if [ ! -n "$1" ]
then
	echo "usage: $(basename $0) [-u] [-s sample_sheet] run_folder [...]" >&2
	exit 1
fi
run_folder=$(readlink -f $1)
shift 1

if [ $sample_sheet ]
then
	sample_sheet_arg="--sample-sheet $sample_sheet"
else
	sample_sheet_arg=''
	sample_sheet="$run_folder/$default_sample_sheet" # needed for checking index type
fi


set -euo pipefail

if grep "$i5_only_regex" "$sample_sheet"
then
	base_mask_arg="$i5_only_option"
else
	base_mask_arg=''
fi


#if [ "$discard_undetermined" = true ]; then ln -s /dev/null $undetermined_filename; fi # breaks bcl2fastq's file-size check
$bcl2fastq_cmd -o . --stats-dir $run_folder/Data/Intensities/BaseCalls/Stats/ --reports-dir $run_folder/Data/Intensities/BaseCalls/Reports/ $sample_sheet_arg $base_mask_arg -R $run_folder "$@"

# clean up filenames
if [ "$discard_undetermined" = true ]; then rm $undetermined_filename; fi
for fastq in *.fastq.gz
	do mv $fastq $(basename $fastq | sed -r "s/_S[0-9]+_R1_001\.fastq\.gz$/.fastq.gz/") # this regex might not be foolproof
done

