#!/bin/bash
# Contact: Dengfeng Guan, dfguan@hit.edu.cn
# Purpose: GATK Variant calling for mulitple samples using fAstCal as reference
# Date: 180524

set -x 
vc() {
	rmdup_srt_bam=$1
	pre_fn=$(basename $rmdup_srt_bam .srt.bam)
	ref=$2
	gatk_path=$3
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -Jvc_"$pre_fn" -ovc_"$pre_fn".o -K -evc_"$pre_fn".e "java -jar -Xmx4G  ${gatk_path} HaplotypeCaller -R ${ref}  -I ${rmdup_srt_bam} -ERC GVCF -o ${pre_fn}.raw.g.vcf"  
}
#consolidate variations
cv() {
	gatk_path=$1
	gvcfs=$2
	ref=$3
	out=$4
	java -jar -Xmx4G $gatk_path GenomicsDBImport $gvcfs --genomicsdb-workspace-path gvcf_db
	java -jar -Xmx4G $gatk_path GenotypeGVCFs -R $ref -V gendb://gvcf_db -G StandardAnnotation -newQual -O $out 
}

#filter variations
fltv() {
	gatk_path=$1
	raw_vcf=$2
	ref=$3
	pref=$(basename $raw_vcf .vcf)
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -Jfltv_"$pre_fn" -ofltv_"$pre_fn".o -K -efltv_"$pre_fn".e "java -jar -Xmx4G ${gatk_path} VariantFiltration -R $ref -V $raw_vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "basic_filtering" -o ${pref}_flted.vcf"
 
}



export -f vc
export -f cv
export -f fltv

ref=$(readlink -f $1)
bam_dir=$(readlink -f $2)
gatk_path=$(readlink -f $3)

echo "Run Haplotype Variant Calling"
find $bam_dir -name "*.srt.bam" | xargs -n 1 -P8 -i bash -c "vc {} $ref $gatk_path &" 
wait 

echo "Consolidating Variations..."
gvcf_dir=$bam_dir
gvcfs=""
find $gvcf_dir -name "*.g.vcf" | xargs -n1 -i gvcfs=$gvcfs"-v {}"
if [ $gvcfs = "" ]
then
	echo "No GVCF Files found"
else
	out=massoko.vcf
  	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -Jvc_"$pre_fn" -ovc_"$pre_fn".o -K -evc_"$pre_fn".e -I cv $gatk_path $gvcfs $ref $out &
	wait
	echo "Filtering Variations..."
	fltv $gatk $out $ref 
fi

