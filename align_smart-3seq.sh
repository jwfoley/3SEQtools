#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference genome, perform the special preprocessing for Smart-3SEQ, align each file to the reference, then index the BAM output
# optionally truncates reads to a specified length before further processing, to simulate results from sequencing fewer cycles
# processes Smart-3SEQ data by trimming first $N_N+$N_G bases: $N_N (the random UMI) is appended to the read name, and $N_G (the G overhang) is discarded
# call from directory where you want results to go
# uses shared memory option in STAR; you may need to increase the kernel's shared memory limits (e.g. "sysctl -w kernel.shmmax=34359738368 kernel.shmall=8388608" in Linux)

samtools_path=samtools
star_path=STAR
umi_trim_path="pypy $(dirname $0)/umi_homopolymer.py -n"
dedup_path='~/umi-dedup/dedup.py -qs'
unzip_path='pigz -dc'
tmp_dir='/tmp/align_smart-3seq_tmp'
N_thread=$(nproc)
bam_mem=2147483648 # maximum bytes of RAM to use for BAM sorting (in addition to the memory usage of the reference index!)
N_N=5
N_G=3
N_A=8
N_mismatch=1
star_options='--outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --clip3pAdapterMMp 0.2' # no multimappers reported, no mismatch filter, and poly(A) clipped


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
shift 1


set -euo pipefail

mkdir -p $tmp_dir
cd $tmp_dir
echo -n 'loading genome into shared memory... ' >&2
$star_path --genomeLoad LoadAndExit --genomeDir $genome_dir > /dev/null
echo 'done' >&2

for fastq_file in "$@"
do
	cd $wd
	rootname=$(basename $fastq_file .fastq.gz)
	fastq=$(readlink -f $fastq_file)
	
	echo -n "processing $rootname... " >&2
	
	cd $tmp_dir
	$unzip_path $fastq |
		$umi_trim_path -u $N_N -g $N_G -p $N_A -m $N_mismatch $truncate_arg 2> $wd/$rootname.trim.log |
		$star_path --genomeLoad LoadAndKeep --genomeDir $genome_dir --readFilesIn /dev/stdin --runThreadN $N_thread --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outBAMcompression 0 --limitBAMsortRAM $bam_mem $star_options |
		$dedup_path 2> $wd/$rootname.dedup.log |
		tee $wd/$rootname.bam |
		$samtools_path index /dev/stdin $wd/$rootname.bai
	touch $wd/$rootname.bai
	cp Log.final.out $wd/$rootname.align.log
	
	echo 'done' >&2
done

$star_path --genomeLoad Remove --genomeDir $genome_dir > /dev/null
cd $wd
rm -rf $tmp_dir

