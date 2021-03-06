#!/bin/bash
# Contact: Dengfeng Guan, dfguan@hit.edu.cn
# Purpose: GATK Variant calling for mulitple samples using fAstCal as reference
# Date: 180522
set -x
align() {
	fn=$1
	idxbase=$2
	picard_path=$3
	pre_fn=$(basename $fn .bam)
	cd $(dirname $fn)
	samtools fastq -1 "$pre_fn"_1.fq -2 "$pre_fn"_2.fq $fn
	echo "align and sort"
	# align & sort bam & rm duplicates
	#bsub -M4000 -n4 -R"span[hosts=1] select[mem>4000] rusage[mem=4000]" -Jalign_"$pre_fn" -K -oalign_"$pre_fn".o "bwa mem -t4 ${idxbase} ${pre_fn}_1.fq ${pre_fn}_2.fq | samtools sort -O BAM | samtools rmdup - ${pre_fn}.srt.rmdup.bam && rm -f ${pre_fn}_[12].fq"  
	bwa mem -R "@RG\tID:$pre_fn\tSM:$pre_fn" -t8 $idxbase "$pre_fn"_1.fq "$pre_fn"_2.fq | samtools view -h -b - | samtools sort - -O BAM -o "$pre_fn".srt.bam && rm -f "$pre_fn"_[12].fq  
	java -jar -Xmx4G -XX:ParallelGCThreads=4  $picard_path MarkDuplicates I="$pre_fn".srt.bam O="$pre_fn".srt.rmdup.bam M="$pre_fn".rmdup_matrix.txt AS=true  
	# gatk recalibrate or maybe group bams first what I am concerned is the fish genomes can be different	
}
if [ "$#" -lt 3 ]
then
	echo align REF SEQ_DIR PICARD_PATH
else
	ref=$(readlink -f $1)
	seq_dir=$(readlink -f $2)
	picard_path=$(readlink -f $3)
	#bwa index $ref  
	export -f align
	#find $seq_dir \( -name "*.srt.bam" -o -name "*.srt.rmdup.bam" \) -exec rm {} \;
	#find $seq_dir -regextype sed -regex ".*[0-9]\{4\}\.bam" -exec rm {} \;
	#bams=`find $seq_dir -name "*.bam"` 
	#for fl in $bams
	#do
		#bsub  -M6000 -n8 -R"span[hosts=1] select[mem>6000] rusage[mem=6000]" -K -Jalign -oalign.o -ealign.e  align $fl $ref $picard_path & 
	#done
	#wait
	srtbams=`find $seq_dir -name "*.srt.rmdup.bam" | tr '\n' '  '`	
	bsub -n8 -qbasement -M100000 -R"select[mem>100000] rusage[mem=100000] span[hosts=1]" -Jsm_vc -osm_vc.o -esm_vc.e "samtools mpileup -t DP,AD,INFO/AD -C50 -p -m2 -F0.2 -ugf $ref $srtbams | bcftools call -vmOz -f GQ -o smvc.vcf.gz"
fi
