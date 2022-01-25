#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference genome, perform the special preprocessing for Smart-3SEQ, align each file to the reference, then either generate a sorted BAM file and index, count gene hits, or both
# preprocesses Smart-3SEQ data by trimming first 8 bases: first 5 (the random UMI) are appended to the read name, and next 3 (the G overhang) are discarded
# trailing P7 adapter sequence and poly(A) from the 1S primer are both removed as well
# writes output files to current working directory
# STAR requires a lot of memory, e.g. about 25 GB for the human genome
# uses shared memory option in STAR; you may need to increase the kernel's shared memory limits (e.g. "sysctl -w kernel.shmmax=34359738368 kernel.shmall=8388608" in Linux)


N_thread=$(nproc)
adapter_trim_path=cutadapt
adapter_trim_options=(
	'-n 2' # trim both sequences from the same read if applicable
	'-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' # TruSeq P7 adapter sequence (if you haven't already trimmed it in bcl2fastq)
	'-a AAAAAAAAACAAAAAAAAACAAAAAAAAAA' # reinforced poly(A) from Smart-3SEQ version 2; for version 1, this is simply a homopolymer, but use a shorter sequence (e.g. AAAAAAAAAA) otherwise sequencing errors may create false negatives
	"-j $N_thread"
)
umi_trim_path="pypy3 $(dirname $0)/extract_umi.py" # if pypy3 is unavailable, use python3 (slower)
umi_trim_options=(
	'-u 5' # UMI is 5 nt
	'-g 3' # G-overhang is 3 nt
	'-l 1' # retain all reads of at least 1 nt (let the aligner decide what to do)
)
star_path=STAR
star_options=(
	'--readFilesIn /dev/stdin'
	'--genomeLoad LoadAndKeep'
	"--runThreadN $N_thread"
	'--scoreGap -999' # don't infer unannotated splice junctions
	'--outFilterMultimapNmax 1'
	'--outFilterMismatchNmax 999' # don't exclude alignments with too many mismatches
)
star_options_bam_true=( # STAR options if we are keeping the BAM file
	'--outSAMtype BAM SortedByCoordinate'
	'--outStd BAM_SortedByCoordinate'
	'--limitBAMsortRAM 2147483648' # maximum bytes of RAM to use for BAM sorting (in addition to the memory usage of the reference index!)
	'--outBAMcompression 10'
)
star_options_bam_false=( # STAR options if we are discarding the BAM file (no need to sort or compress)
	'--outSAMtype BAM Unsorted'
	'--outStd BAM_Unsorted'
	'--outBAMcompression 0'
)
samtools_path=samtools
count_path=featureCounts
count_options=(
	"-T $N_thread"
	'-s 1' # correct strand orientation required
	'--read2pos 5' # reduce each read to its 5'-most position to simplify overlap problems
	'--primary' # use only primary alignments
)


gtf_file=
make_bam=true
while getopts ":f:n" opt
do
	case $opt in
		f)
			gtf_file="$(readlink -f $OPTARG)"
			if [ ! -e $gtf_file ]; then echo "error: $gtf_file not found" >&2; exit 1; fi
			;;
		n)
			make_bam=
			;;
	esac
done
shift "$((OPTIND-1))"

if [ ! -n "$2" ]
then
	echo "
usage: $(basename $0) [-f annotation.gtf] [-n] genome_dir file1.fastq.gz file2.fastq.gz file3.fastq.gz ...

-f annotation.gtf: count hits in genes from this annotation file
-n: don't make BAM files

If '-n' is used without '-f', the only output is log files from the trimmer and aligner.
" >&2
	exit 1
fi


set -euo pipefail


wd=$(pwd)
genome_dir=$(readlink -f $1)
shift 1


if [ $make_bam ]
then
	star_bam_options=${star_options_bam_true[@]}
else
	star_bam_options=${star_options_bam_false[@]}
fi


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
	$adapter_trim_path $fastq 2> $wd/$rootname.adapter.log |
		$umi_trim_path ${umi_trim_options[@]} 2> $wd/$rootname.umi.log |
		$star_path ${star_options[@]} $star_bam_options --genomeDir $genome_dir |
		tee >(
			if [ $make_bam ]
			then
				tee $wd/$rootname.bam |
				$samtools_path index /dev/stdin $wd/$rootname.bai
			else
				cat > /dev/null # no-op to prevent pipe failure
			fi
		) |
		if [ $gtf_file ]
		then
			$count_path ${count_options[@]} -a $gtf_file -o counts 2> $wd/$rootname.count.log
		else
			cat > /dev/null # no-op to prevent pipe failure
		fi
	
	cp Log.final.out $wd/$rootname.align.log
	if [ $make_bam ]; then touch $wd/$rootname.bai; fi # ensure the index is younger than the BAM to avoid warnings
	
	# generate count file with header
	if [ $gtf_file ]
	then
		echo "# $fastq" > $wd/$rootname.counts.tsv
		echo "# $adapter_trim_path ${adapter_trim_options[@]}" >> $wd/$rootname.counts.tsv
		echo "# $umi_trim_path ${umi_trim_options[@]}" >> $wd/$rootname.counts.tsv
		echo "# $star_path ${star_options[@]} $star_bam_options --genomeDir $genome_dir" >> $wd/$rootname.counts.tsv
		echo "# $count_path ${count_options[@]} -a $gtf_file" >> $wd/$rootname.counts.tsv
		tail -n +3 counts | cut -f 1,7 >> $wd/$rootname.counts.tsv
		mv counts.summary $wd/$rootname.counts.summary
	fi
	cd $wd
	rm -rf $tmp_dir
	
	echo 'done' >&2
done

cd $genome_tmp_dir
$star_path --genomeLoad Remove --genomeDir $genome_dir > /dev/null
cd $wd
rm -rf $genome_tmp_dir

echo 'all done' >&2

