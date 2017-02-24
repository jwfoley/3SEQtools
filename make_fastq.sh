#! /bin/bash

# given the path to an Illumina run folder and optionally a path to a sample sheet, run bcl2fastq
# then clean up the names of the FASTQ files and move them to the current directory
# the sample sheet argument is necessary unless it's the default $run_folder/SampleSheet.csv

bcl2fastq_path=bcl2fastq
bcl2fastq_args='--fastq-compression-level 9 --no-lane-splitting --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0'
fastq_subdir='Data/Intensities/BaseCalls'


sample_sheet_arg=''
while getopts ":s:" opt
do
	case $opt in
		s)
			if [ ! -e "$OPTARG" ]; then echo "error: $OPTARG not found" >&2; exit 1; fi
			sample_sheet_arg="--sample-sheet $OPTARG"
			;;
	esac
done
shift "$((OPTIND-1))"

if [ ! -n "$1" ]
then
	echo "usage: $(basename $0) [-s sample_sheet] run_folder" >&2
	exit 1
fi


set -euo pipefail

$bcl2fastq_path $bcl2fastq_args $sample_sheet_arg -R $run_folder

for fastq in $run_folder/$fastq_subdir/*.fastq.gz
do mv $fastq $(basename $fastq | sed -r "s/_S[0-9]+_R1_001\.fastq\.gz$/.fastq.gz/") # this regex might not be foolproof
done

