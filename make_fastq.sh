#! /bin/bash

# given the path to an Illumina run folder and optionally a path to a sample sheet, run bcl2fastq
# generates FASTQ files in the current directory and cleans up their names
# the sample sheet argument is necessary unless it's the default $run_folder/SampleSheet.csv
# discards "Undetermined" reads by default unless "-u" is used

default_undetermined_filename='Undetermined_S0_R1_001.fastq.gz'

bcl2fastq_path=bcl2fastq

sample_sheet_arg=''
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
	echo "usage: $(basename $0) [-u] [-s sample_sheet] run_folder" >&2
	exit 1
fi

run_folder=$(readlink -f $1)


set -euo pipefail

if [ "$discard_undetermined" = true ]; then ln -s /dev/null $default_undetermined_filename; fi
$bcl2fastq_path --fastq-compression-level 9 --no-lane-splitting --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 -o . --stats-dir $run_folder/Data/Intensities/BaseCalls/Stats/ --reports-dir $run_folder/Data/Intensities/BaseCalls/Reports/ $sample_sheet_arg -R $run_folder

# clean up filenames
if [ "$discard_undetermined" = true ]; then rm $default_undetermined_filename; fi
for fastq in *.fastq.gz
	do mv $fastq $(basename $fastq | sed -r "s/_S[0-9]+_R1_001\.fastq\.gz$/.fastq.gz/") # this regex might not be foolproof
done

