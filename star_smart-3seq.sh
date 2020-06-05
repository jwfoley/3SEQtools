#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference genome, perform the special preprocessing for Smart-3SEQ, align each file to the reference, then index the BAM output
# optionally truncates reads to a specified length before further processing, to simulate results from sequencing fewer cycles
# processes Smart-3SEQ data by trimming first 8 bases: first 5 (the random UMI) are appended to the read name, and next 3 (the G overhang) are discarded
# writes output files to current working directory
# STAR requires a lot of memory, e.g. about 25 GB for the human genome
# uses shared memory option in STAR; you may need to increase the kernel's shared memory limits (e.g. "sysctl -w kernel.shmmax=34359738368 kernel.shmall=8388608" in Linux)


N_thread=$(nproc)
unzip_path='pigz -dc'
umi_trim_path="pypy $(dirname $0)/umi_homopolymer.py"
umi_trim_options=(
	'-n' # don't trim poly(A); let STAR do that
)
star_path=STAR
star_options=(
	'--readFilesIn /dev/stdin'
	'--genomeLoad LoadAndKeep'
	'--outSAMtype BAM SortedByCoordinate'
	'--outStd BAM_SortedByCoordinate'
	"--runThreadN $N_thread"
	'--limitBAMsortRAM 2147483648' # maximum bytes of RAM to use for BAM sorting (in addition to the memory usage of the reference index!)
	'--outBAMcompression 10'
	'--outFilterMultimapNmax 1'
	'--outFilterMismatchNmax 999' # don't filter alignments with too many mismatches
	'--clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' # trim poly(A)
	'--clip3pAdapterMMp 0.2' # trim poly(A) leniently
)
samtools_path=samtools


truncate_arg=''
while getopts ":n:g:t" opt
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


genome_tmp_dir=$(mktemp -d --suffix .align_smart-3seq.genome)
cd $genome_tmp_dir
echo -n 'loading genome into shared memory... ' >&2
$star_path --genomeLoad LoadAndExit --genomeDir $genome_dir > /dev/null
cd $wd
echo 'done' >&2


for fastq_file in "$@"
do
	rootname=$(basename $fastq_file .fastq.gz)
	fastq=$(readlink -f $fastq_file)
	
	echo -n "aligning $rootname... " >&2
	
	tmp_dir=$(mktemp -d --suffix .align_smart-3seq)
	cd $tmp_dir
	$unzip_path $fastq |
		$umi_trim_path ${umi_trim_options[@]} $truncate_arg 2> $wd/$rootname.trim.log |
		$star_path ${star_options[@]} --genomeDir $genome_dir |
		tee $wd/$rootname.bam |
		$samtools_path index /dev/stdin $wd/$rootname.bai
	touch $wd/$rootname.bai # ensure the index is younger than the BAM to avoid warnings
	cp Log.final.out $wd/$rootname.align.log
	cd $wd
	rm -rf $tmp_dir
	
	echo 'done' >&2
done

cd $genome_tmp_dir
$star_path --genomeLoad Remove --genomeDir $genome_dir > /dev/null
cd $wd
rm -rf $genome_tmp_dir

