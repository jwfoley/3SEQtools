#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference genome, perform the special preprocessing for Smart-3SEQ, align each file to the reference, then index the BAM output
# optionally truncates reads to a specified length before further processing, to simulate results from sequencing fewer cycles
# processes Smart-3SEQ data by trimming first $N_N+$N_G bases: $N_N (the random UMI) is appended to the read name, and $N_G (the G overhang) is discarded
# call from directory where you want results to go
# uses shared memory option in STAR; you may need to increase the kernel's shared memory limits (e.g. "sysctl -w kernel.shmmax=34359738368 kernel.shmall=8388608" in Linux)
# creates temporary BAM files and then deduplicates them all in parallel
# or you can disable deduplication with '-d'

samtools_path=samtools
star_path=STAR
parallel_path=parallel
umi_trim_path="pypy $(dirname $0)/umi_homopolymer.py -n"
dedup_command="$HOME/umi-dedup/dedup.py -qs"
unzip_path='pigz -dc'
N_thread=$(nproc)
bam_mem=2147483648 # maximum bytes of RAM to use for BAM sorting (in addition to the memory usage of the reference index!)
N_N=5
N_G=3
N_A=8
N_mismatch=1
star_options='--outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --clip3pAdapterMMp 0.2' # no multimappers reported, no mismatch filter, and poly(A) clipped


truncate_arg=''
dedup=true
while getopts ":n:g:t:d" opt
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
		d)
			dedup=false
			;;
	esac
done
shift "$((OPTIND-1))"

if [ ! -n "$2" ]
then
	echo "usage: $(basename $0) [-n umi_length] [-g discard_length] [-t truncate_length] [-d] genome_dir file1.fastq.gz file2.fastq.gz file3.fastq.gz ..." >&2
	exit 1
fi


wd=$(pwd)
genome_dir=$(readlink -f $1)
shift 1


set -euo pipefail

if $dedup
then
	loop_suffix='.bam_tmp'
	loop_index_command='cat > /dev/null'
	loop_touch_command=':'
	compression_level=0
else
	loop_suffix='.bam'
	loop_index_command="$samtools_path index /dev/stdin \$wd/\$rootname.bai"
	loop_touch_command='touch $wd/$rootname.bai'
	compression_level=9
fi


genome_tmp_dir=$(mktemp -d --suffix .align_smart-3seq.genome)
cd $genome_tmp_dir
echo -n 'loading genome into shared memory... ' >&2
$star_path --genomeLoad LoadAndExit --genomeDir $genome_dir > /dev/null
echo 'done' >&2


tmp_bam_files=()
for fastq_file in "$@"
do
	cd $wd
	rootname=$(basename $fastq_file .fastq.gz)
	fastq=$(readlink -f $fastq_file)
	outbam="$wd/$rootname$loop_suffix"
	
	echo -n "aligning $rootname... " >&2
	
	tmp_dir=$(mktemp -d --suffix .align_smart-3seq)
	cd $tmp_dir
	$unzip_path $fastq |
		$umi_trim_path -u $N_N -g $N_G -p $N_A -m $N_mismatch $truncate_arg 2> $wd/$rootname.trim.log |
		$star_path --genomeLoad LoadAndKeep --genomeDir $genome_dir --readFilesIn /dev/stdin --runThreadN $N_thread --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outBAMcompression $compression_level --limitBAMsortRAM $bam_mem $star_options |
		tee $outbam |
		eval $loop_index_command
	eval $loop_touch_command
	cp Log.final.out $wd/$rootname.align.log
	cd $wd
	rm -rf $tmp_dir
	tmp_bam_files+=($outbam)
	
	echo 'done' >&2
done

cd $genome_tmp_dir
$star_path --genomeLoad Remove --genomeDir $genome_dir > /dev/null
cd $wd
rm -rf $genome_tmp_dir

if $dedup
then
	echo 'deduplicating alignments... ' >&2
	$parallel_path "
		$dedup_command {} 2> {.}.dedup.log |
			tee {.}.bam |
			$samtools_path index /dev/stdin {.}.bai
		rm {}
		touch {.}.bai
	" ::: ${tmp_bam_files[*]}
	echo 'all done!' >&2
fi

