#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference transcriptome index, perform the special preprocessing for Smart-3SEQ, quasi-align each file to the reference with Salmon
# optionally truncates reads to a specified length before further processing, to simulate results from sequencing fewer cycles
# processes Smart-3SEQ data by trimming first $N_N+$N_G bases: $N_N (the random UMI) is appended to the read name, and $N_G (the G overhang) is discarded
# call from directory where you want results to go

samtools_path=samtools
salmon_path=salmon
umi_trim_path="pypy $(dirname $0)/umi_homopolymer.py -n"
unzip_path='zstd -dc'
N_N=5
N_G=3
N_A=8
N_mismatch=1
salmon_options='-l SF --noLengthCorrection'


truncate_arg=''
make_bam=false
while getopts ":n:d:t:g:" opt
do
	case $opt in
		n)
			if [[ $OPTARG != +([0-9]) ]]; then echo "error: $OPTARG is not a valid UMI length" >&2; exit 1; fi
			N_N="$OPTARG"
			;;
		d)
			if [[ $OPTARG != +([0-9]) ]]; then echo "error: $OPTARG is not a valid discard length" >&2; exit 1; fi
			N_G="$OPTARG"
			;;
		t)
			if [[ $OPTARG != +([0-9]) ]]; then echo "error: $OPTARG is not a valid length to truncate" >&2; exit 1; fi
			truncate_arg="-L $OPTARG"
			;;
		g)
			salmon_options="$salmon_options -g $(readlink -f $OPTARG)"
			;;
	esac
done
shift "$((OPTIND-1))"

if [ ! -n "$2" ]
then
	echo "usage: $(basename $0) [-n umi_length] [-g discard_length] [-t truncate_length] genome_dir file1.fastq.gz file2.fastq.gz file3.fastq.gz ..." >&2
	exit 1
fi


wd=$(pwd)
genome_dir=$(readlink -f $1)
shift 1


set -euo pipefail


for fastq_file in "$@"
do
	cd $wd
	rootname=$(basename $fastq_file .fastq.gz)
	fastq=$(readlink -f $fastq_file)
	
	echo -e "processing $rootname...\n\n" >&2
	
	tmp_dir=$(mktemp -d --suffix .quasi-align_smart-3seq)
	cd $tmp_dir
	$unzip_path $fastq |
		$umi_trim_path -u $N_N -g $N_G -p $N_A -m $N_mismatch $truncate_arg 2> $wd/$rootname.trim.log |
		$salmon_path quant $salmon_options -i $genome_dir -r /dev/stdin -o .
	cp quant.sf $wd/$rootname.quant.tsv
	cp logs/salmon_quant.log $wd/$rootname.quasi-align.log
	
	cd $wd
	rm -rf $tmp_dir
	
	echo -e "$rootname done\n" >&2
done

