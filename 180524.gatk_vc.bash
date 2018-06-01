#!/bin/bash
# Contact: Dengfeng Guan, dg30@sanger.ac.uk, dfguan@hit.edu.cn
# Purpose: GATK Variant calling for mulitple samples using fAstCal as reference
# Date: 180524
md_vc() {
	cur_dir=$(pwd)
	srt_bam=$1
	pre_fn=$(basename $srt_bam .srt.bam)
	ref=$2
	picard_path=$3
	gatk_path=$4
	cd $(dirname $srt_bam)
	#bsub -M4000 -n4 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -Jmd -omd.o -K -emd.e "java -jar -Xmx4G -XX:ParallelGCThreads=4  $picard_path MarkDuplicates I=${srt_bam} O=${pre_fn}.srt.rmdup.bam M=${pre_fn}.rmdup_matrix.txt AS=true" & 
	#wait
	#bsub -M1000 -R"select[mem>1000] rusage[mem=1000]" -J idx -o idx.o -e idx.e -K "samtools index ${pre_fn}.srt.rmdup.bam" &
	#wait
	bsub -M4000 -n4 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -Jvc_${pre_fn} -ovc_${pre_fn}.o -K -evc_${pre_fn}.e "java -jar -Xmx4G  -XX:ParallelGCThreads=4 ${gatk_path} HaplotypeCaller -R ${ref}  -I ${pre_fn}.srt.rmdup.bam -ERC GVCF -G StandardAnnotation -GAS_StandardAnnotation -O ${pre_fn}.raw.g.vcf" & 
	wait
	cd $cur_dir
}
vc() {
	cur_dir=$(pwd)
	rmdup_srt_bam=$1
	pre_fn=$(basename $rmdup_srt_bam .srt.rmdup.bam)
	ref=$2
	gatk_path=$3
	cd $(dirname $rmdup_srt_bam)
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -Jvc_"$pre_fn" -ovc_"$pre_fn".o -K -evc_"$pre_fn".e "java -jar -Xmx4G  ${gatk_path} HaplotypeCaller -R ${ref}  -I ${rmdup_srt_bam} -ERC GVCF -G StandardAnnotation -GAS_StandardAnnotation -O ${pre_fn}.raw.g.vcf" & 
	cd $cur_dir
	#sleep 100
}

#consolidate variations
cv() {
	gatk_path=$1
	gvcfs=$2
	ref=$3
	out=$4
	#java -jar -Xmx4G $gatk_path GenomicsDBImport $gvcfs --genomicsdb-workspace-path gvcf_db
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -Jcv_gvcf -ocv_gvcf.o -K -ecv_gvcf.e "java -jar -Xmx4G $gatk_path CombineGVCFs $gvcfs -R $ref -O massoko.cmb.g.vcf && java -jar -Xmx4G $gatk_path GenotypeGVCFs -R $ref -V massoko.cmb.g.vcf -O $out" &
}

#filter variations
fltv() {
	gatk_path=$1
	raw_vcf=$2
	ref=$3
	pref=$(basename $raw_vcf .vcf)
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -Jfltv_"$pre_fn" -ofltv_"$pre_fn".o -efltv_"$pre_fn".e "java -jar -Xmx4G ${gatk_path} VariantFiltration -R $ref -V $raw_vcf -filter "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "basic_filtering" -O ${pref}_flted.vcf"
}


export -f vc
#export -f md
export -f cv
export -f fltv
export -f md_vc
if [ "$#" -lt 4 ]
then
	echo gatk_vc REFERENCE BAM_DIR GATK_PATH PICARD_PATH
else
	ref=$(readlink -f $1)
	bam_dir=$(readlink -f $2)
	gatk_path=$(readlink -f $3)
	picard_path=$(readlink -f $4)
	pre_ref=$(basename $ref .fa)
	echo "Create Sequence Dictionary"
	#samtools dict $ref -o "$pre_ref".dict
	echo "Marking Duplicates..."
	find $bam_dir -name "*.srt.bam" | xargs -n 1  -i bash -c "md_vc {} $ref $picard_path $gatk_path &" 
	#wait
	#echo "Run Haplotype Variant Calling"
	#find $bam_dir -name "*.srt.rmdup.bam" | xargs -n 1 -P8 -i bash -c "vc {} $ref $gatk_path" 
	#find $bam_dir -name "*.srt.rmdup.bam" | xargs -n 1 -P8 -i bash -c "vc {} $ref $gatk_path"
	wait
	echo "Consolidating Variations..."
	gvcf_dir=$bam_dir
	gvcfs=""
	find $gvcf_dir -name "*.g.vcf" | xargs -n1 -i gvcfs=$gvcfs"-V {} "
	if [ "$gvcfs" = "" ]
	then
		echo "No GVCF Files found"
	else
		out=massoko.vcf
		cv $gatk_path $gvcfs $ref $out 
		wait
		echo "Filtering Variations..."
		fltv $gatk $out $ref 
	fi
fi
