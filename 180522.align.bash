#!/bin/bash
# Contact: Dengfeng Guan, dfguan@hit.edu.cn
# Purpose: GATK Variant calling for mulitple samples using fAstCal as reference
# Date: 180522
set -x
align() {
	fn=$1
	idxbase=$2
	pre_fn=$(basename $fn .bam)
	cd $(dirname $fn)
	samtools fastq -1 "$pre_fn"_1.fq -2 "$pre_fn"_2.fq $fn
	echo "align and sort"
	# align & sort bam & rm duplicates
	#bsub -M4000 -n4 -R"span[hosts=1] select[mem>4000] rusage[mem=4000]" -Jalign_"$pre_fn" -K -oalign_"$pre_fn".o "bwa mem -t4 ${idxbase} ${pre_fn}_1.fq ${pre_fn}_2.fq | samtools sort -O BAM | samtools rmdup - ${pre_fn}.srt.rmdup.bam && rm -f ${pre_fn}_[12].fq"  
	bsub -M4000 -n4 -R"span[hosts=1] select[mem>4000] rusage[mem=4000]" -Jalign_"$pre_fn" -oalign_"$pre_fn".o "bwa mem -t4 ${idxbase} ${pre_fn}_1.fq ${pre_fn}_2.fq | samtools sort -O BAM -o ${pre_fn}.srt.bam && rm -f ${pre_fn}_[12].fq" & 
	# gatk recalibrate or maybe group bams first what I am concerned is the fish genomes can be different	
}

ref=$(readlink -f $1)
seq_dir=$(readlink -f $2)
#bwa index $ref  
export -f align
find $seq_dir -name "*.srt.bam" -o -name "*.srt.rmdup.bam" -exec rm {} \;
find $seq_dir -name "*.bam" | xargs -n 1 -P8 -i bash -c "align {} $ref" 
