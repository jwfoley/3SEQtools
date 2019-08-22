#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference genome, perform the special preprocessing for Smart-3SEQ, align each file to the reference with Novoalign 4, then index the BAM output
# optionally truncates reads to a specified length before further processing, to simulate results from sequencing fewer cycles
# processes Smart-3SEQ data by trimming first 8 or more bases: the first 5 are extracted as the UMI (RX and QX tags in SAM format), and the next 3 (or more) G's are discarded
# call from directory where you want results to go
# creates temporary BAM files and then deduplicates them all in parallel
# or you can disable deduplication with '-d'

samtools_path=samtools
novoalign_path=novoalign4
dedup_command="$HOME/umi-dedup/dedup.py -qs"
N_thread=$(nproc)
bam_mem=2147483648 # maximum bytes of RAM to use for BAM sorting (in addition to the memory usage of the reference index!)
fasta_suffix=.fa
index_suffix=.ndx
novoalign_options="-F ILM1.8 -o BAM -c $N_thread --tune NextSeq -5 XXXXXGGGGGGGGGG,8 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
samtools_view_options='-u -@ $N_thread -F 0x4 -F 0x100 -F 0x200 -F 0x800' # remove unaligned or undesirable reads
samtools_sort_options="-l 9 -@ $N_thread -m $(($bam_mem / $N_thread))"
samtools_index_options="-@ $N_thread"


if [ ! -n "$2" ]
then
	echo "usage: $(basename $0) [-d] genome_prefix file1.fastq.gz file2.fastq.gz file3.fastq.gz ..." >&2
	exit 1
fi


wd=$(pwd)
genome_prefix=$(readlink -f $1)
shift 1


set -euo pipefail


for fastq_file in "$@"
do
	cd $wd
	rootname=$(basename $fastq_file .fastq.gz)
	fastq=$(readlink -f $fastq_file)
	
	echo -n "aligning $rootname... " >&2
	
	$novoalign_path $novoalign_options -d $genome_prefix$index_suffix -f $fastq_file 2> $wd/$rootname.novoalign.log |
		$samtools_path view $samtools_view_options |
		$samtools_path sort $samtools_sort_options |
		tee $wd/$rootname.bam |
		$samtools_path index $samtools_index_options /dev/stdin $wd/$rootname.bai
	touch $wd/$rootname.bai	# guarantees the index is newer than the main BAM to prevent warnings later
		
	echo 'done' >&2
done

