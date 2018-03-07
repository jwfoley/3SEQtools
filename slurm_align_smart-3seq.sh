#! /bin/bash

# given a list of single-end Illumina FASTQ files and a reference genome, perform the special preprocessing for Smart-3SEQ, align each file to the reference, then index the BAM output
# optionally truncates reads to a specified length before further processing, to simulate results from sequencing fewer cycles
# processes Smart-3SEQ data by trimming first $N_N+$N_G bases: $N_N (the random UMI) is appended to the read name, and $N_G (the G overhang) is discarded
# call from directory where you want results to go
# submits all tasks to SLURM in an array

job_name=slurm_align_smart-3seq
samtools_path=samtools
star_path=STAR
umi_trim_path="pypy $(dirname $0)/umi_homopolymer.py -n"
dedup_path='~/umi-dedup/dedup.py -qs'
unzip_path='pigz -dc'
tmp_dir='~/tmp'
N_thread=8
N_N=5
N_G=3
N_A=8
N_mismatch=1
star_options='--outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --clip3pAdapterMMp 0.2' # no multimappers reported, no mismatch filter, and poly(A) clipped
time='1:00:00'
mail_type='FAIL'
genome_test='Genome' # check for this file to verify valid genome directory


truncate_arg=''
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
			dedup_path=cat
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

fastq_files=()
for fastq_file in "$@"
do
	if [[ "$fastq_file" != *.fastq.gz ]]
	then
		echo "error: $fastq_file is not a .fastq.gz file" >&2
		exit 1
	fi
	fastq_files+=($(readlink -f $fastq_file))
done

echo "#! /bin/bash
	set -euo pipefail
	fastqs=(${fastq_files[@]})
	fastq=\${fastqs[\$SLURM_ARRAY_TASK_ID]}
	rootname=\$(basename \$fastq .fastq.gz)
	
	this_tmp_dir=$tmp_dir/\$SLURM_ARRAY_JOB_ID\_\$SLURM_ARRAY_TASK_ID
	mkdir \$this_tmp_dir
	cd \$this_tmp_dir
	$unzip_path \$fastq |
		$umi_trim_path -u $N_N -g $N_G -p $N_A -m $N_mismatch $truncate_arg 2> $wd/\$rootname.trim.log |
		$star_path --genomeDir $genome_dir --readFilesIn /dev/stdin --runThreadN $N_thread --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outBAMcompression 0 $star_options |
		$dedup_path 2> $wd/\$rootname.dedup.log |
		tee $wd/\$rootname.bam |
		$samtools_path index /dev/stdin $wd/\$rootname.bai
	touch $wd/\$rootname.bai
	cp Log.final.out $wd/\$rootname.align.log
	
	cd ..
	rm -rf \$this_tmp_dir
" | sbatch --array=0-$((${#fastq_files[@]} - 1)) --cpus-per-task=$N_thread --job-name=$job_name --output=$job_name.log --time=$time --mail-type=$mail_type

