#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference genome, perform the special preprocessing for 3SEQ, align each file to the reference, then index the BAM output
# removes A-tail and anything after it, defining A-tail as a run of $N_A A's with $N_mismatch mismatches allowed
# call from directory where you want results to go

samtools_path=samtools
star_path=STAR
homopolymer_trim_path='pypy ~/script/trim_homopolymer.py'
unzip_path='pigz -dc'
tmp_dir='$LOCAL_SCRATCH'
N_thread=8
N_A=8
N_mismatch=1
star_options='--outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterType BySJout --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000' # ENCODE options per manual, except no multimappers reported
time='1:00:00'
mail_type='FAIL'


if [ ! -n "$2" ]
then
	echo "usage: $(basename $0) genome_dir file1.fastq.gz file2.fastq.gz file3.fastq.gz ..."
	exit 1
fi

wd=$(pwd)
genome_dir=$(readlink -f $1)
shift 1


for fastq_file in "$@"
do 
	rootname=$(basename $fastq_file .fastq.gz)
	fastq=$(readlink -f $fastq_file)
	echo "#! /bin/bash
		set -euo pipefail
		cd $tmp_dir
		$unzip_path $fastq |
			$homopolymer_trim_path -p $N_A -m $N_mismatch 2> $wd/$rootname\_trim.log |
			$star_path --genomeDir $genome_dir --readFilesIn /dev/stdin --runThreadN $N_thread --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outBAMcompression 10 $star_options |
			tee $wd/$rootname.bam |
			$samtools_path index /dev/stdin $wd/$rootname.bai
		cp Log.final.out $wd/$rootname\_star.log
" | sbatch --cpus-per-task=$N_thread --job-name=$rootname\_align --output=$rootname\_align_job.log --time=$time --mail-type=$mail_type
done

