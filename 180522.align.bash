#!/bin/bash
# Contact: Dengfeng Guan, dfguan@hit.edu.cn
# Purpose: GATK Variant calling for mulitple samples using fAstCal as reference
# Date: 180522

align() {
	fn=$1
	idxbase=$2
	pre_fn=$(basename $fn .bam)
	samtools fastq -1 "$pre_fn"_1.fq -2 "$pref_fn"_2.fq $fn
	# align and sort bam
	bsub -M4000 -n4 -R"span[hosts=1] select[mem>4000] rusage[mem=4000]" -Jalign_"$pre_fn" -K1 -oalign_"$pre_fn".o "bwa mem -t4 ${idxbase} ${pref_fn}_1.fq ${pref_fn}_2.fq | samtools sort -O BAM -o ${pref_fn}.srt.bam &" 
	# samtools remove duplicates
	samtools rmdup ${pref_fn}.srt.bam > ${pref_fn}.srt.rmdup.bam
	rm -f ${pref_fn}.srt.bam
	# gatk recalibrate or maybe group bams first what I am concerned is the fish genomes can be different	
}
ref=$1
seq_dir=$2
bwa index $ref  
find seq_dir -name "*.bam" | xargs -n 1 -P8 -i sh -c "align {} $ref" 
