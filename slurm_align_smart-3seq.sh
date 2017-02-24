#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference genome, perform the special preprocessing for Smart-3SEQ, align each file to the reference, then index the BAM output
# optionally truncates reads to a specified length before further processing, to simulate results from sequencing fewer cycles
# processes Smart-3SEQ data by trimming first $N_N+$N_G bases: $N_N (the random UMI) is appended to the read name, and $N_G (the G overhang) is discarded
# call from directory where you want results to go

samtools_path=samtools
star_path=STAR
umi_trim_path='pypy ~/3SEQtools/umi_homopolymer.py -n'
dedup_path='python3 ~/umi-bayes/dedup.py -qa naive'
unzip_path='pigz -dc'
tmp_dir='$LOCAL_SCRATCH'
modules='python/3.3.2'
N_thread=8
N_N=5
N_G=3
N_A=8
N_mismatch=1
star_options='--outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterType BySJout --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --clip3pAdapterMMp 0.2' # ENCODE options per manual, except no multimappers reported and poly(A) clipped
time='1:00:00'
mail_type='FAIL'
genome_test='Genome' # check for this file to verify valid genome directory


truncate_arg=''
while getopts ":n:g:t:" opt
do
	case $opt in
		n)
			if [[ $OPTARG != +([0-9]) ]]; then echo "error: $OPTARG is not a valid UMI length" >&2; exit 1; fi
			N_N="$OPTARG"
			;;
		g)
			if [[ $OPTARG != +([0-9]) ]]; then echo "error: $OPTARG is not a valid discard length" >&2; exit 1; fi
			N_G="$OPTARG"
			;;
		t)
			if [[ $OPTARG != +([0-9]) ]]; then echo "error: $OPTARG is not a valid length to truncate" >&2; exit 1; fi
			truncate_arg="-L $OPTARG"
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
if [ ! -e "$genome_dir/$genome_test" ]; then echo "error: $genome_dir is not a valid STAR genome directory" >&2; exit 1; fi
shift 1


for fastq_file in "$@"
do 
	rootname=$(basename $fastq_file .fastq.gz)
	fastq=$(readlink -f $fastq_file)
	echo "#! /bin/bash
		set -euo pipefail
		module load $modules
		cd $tmp_dir
		$unzip_path $fastq |
			$umi_trim_path -u $N_N -g $N_G -p $N_A -m $N_mismatch $truncate_arg 2> $wd/$rootname.trim.log |
			$star_path --genomeDir $genome_dir --readFilesIn /dev/stdin --runThreadN $N_thread --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outBAMcompression 0 $star_options |
			$dedup_path 2> $wd/$rootname.dedup.log |
			tee $wd/$rootname.bam |
			$samtools_path index /dev/stdin $wd/$rootname.bai
		touch $wd/$rootname.bai
		cp Log.final.out $wd/$rootname.align.log
" | sbatch --cpus-per-task=$N_thread --job-name=$rootname\_align --output=$rootname.job.log --time=$time --mail-type=$mail_type
done

