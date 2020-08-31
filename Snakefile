#uBAM -> fq -> srt.bam -> markduplicate -> gvcf 
# faidx, dict-> gvcfs -> gt -> filter

import os

configfile: "config.json"

uBAMs_l = config["ubams"]
fpaths = [d[:-4] for d in uBAMs_l]
dirs = [os.path.dirname(d) for d in uBAMs_l]
samples = [os.path.basename(d)[:-4] for d in uBAMs_l]

ref = config["ref"]
picard_pth = config["picard_pth"]
gatk_pth = config["gatk_pth"]

out = config["out"] 

rule all:
	input: out 

rule bwa_idx:
	output: 	
		ref+".amb", ref+".ann", ref+".bwt"
	shell:
		"bwa index {ref}"

rule bam2fq:
	input:
		"{fpath}.bam"
	output:
		o1 = "{fpath}_1.fq", 
		o2 = "{fpath}_2.fq"
	shell:
		"samtools fastq -1 {output.o1} -2 {output.o2} {input}"

rule bwa_map:
	input: 
		ref+".bwt",
		i1 = "{fpath}_1.fq",
		i2 = "{fpath}_2.fq" 
	output:
		srtbam = "{fpath}.srt.bam", 
		lsf_o = "{fpath}_aln.o",
		lsf_e = "{fpath}_aln.e",
	params:
		"@RG\\tID:"+os.path.basename("{fpath}")+"\\tSM:"+os.path.basename("{fpath}")
	shell:
		""" bsub -qlong -M4000 -n4 -R"span[hosts=1] select[mem>4000] rusage[mem=4000]" -K -Jaln -o{output.lsf_o} -e {output.lsf_e} "bwa mem -t4 -R '{params}' {ref} {input.i1} {input.i2} | samtools view -h -b - | samtools sort - -o {output.srtbam}" """

rule bam_idx:
	input:
		srtbam = "{fpath}.srt.rmdup.bam"
	output:	
		bami = "{fpath}.srt.rmdup.bam.bai"
	shell:
		"samtools index {input.srtbam}"
	
rule markdups:
	input: 
		i = "{fpath}.srt.bam",
	output: 
		o_bam = "{fpath}.srt.rmdup.bam",
		o_mat = "{fpath}.rmdup_matrix.txt",
		lsf_o = "{fpath}_md.o",
		lsf_e = "{fpath}_md.e",
		#lsf_j = os.path.dirname("{fpath}")+"/"+"md_"+os.path.basename("{fpath}")[:-4]+".j",
		#lsf_o = os.path.dirname("{fpath}")+"/"+"md_"+os.path.basename("{fpath}")[:-4]+".o",
		#lsf_e = os.path.dirname("{fpath}")+"/"+"md_"+os.path.basename("{fpath}")[:-4]+".e"
	shell:
		""" bsub -n4 -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -K -Jmd -o{output.lsf_o} -e{output.lsf_e} "java -jar -Xmx4G -XX:ParallelGCThreads=4 {picard_pth} MarkDuplicates I={input.i} O={output.o_bam} M={output.o_mat} AS=true" """
rule fa_idx:
	output: 
		dict_fa = ref[:-3]+".dict",
		fai_fa = ref+".fai"
	shell:
		"""
			samtools faidx {ref} && samtools dict {ref} -o {output.dict_fa}
		"""	
rule hapcall: # haplotyecaller
	input:
		#dict_fa = "{ref}[:-3].dict", //error here ref can't be determined
		#fai_fa = "{ref}.fai", 
		dict_fa = ref[:-3]+".dict",
		fai_fa = ref+".fai",
		bami = "{fpath}.srt.rmdup.bam.bai",
		srtbam = "{fpath}.srt.rmdup.bam"
	output:
		ogvcf = "{fpath}.raw.g.vcf.gz",
		#lsf_j = "{fpath}_hc.j",
		lsf_o = "{fpath}_hc.o",
		lsf_e = "{fpath}_hc.e",
		#lsf_j = os.path.dirname("{fpath}")+"/"+"hc_"+os.path.basename("{fpath}")[:-4]+".j",
		#lsf_o = os.path.dirname("{fpath}")+"/"+"hc_"+os.path.basename("{fpath}")[:-4]+".o",
		#lsf_e = os.path.dirname("{fpath}")+"/"+"hc_"+os.path.basename("{fpath}")[:-4]+".e"
	shell:
		""" bsub -n4 -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -Jhc -o{output.lsf_o} -e{output.lsf_e} -K "java -jar -Xmx4G -XX:ParallelGCThreads=4 {gatk_pth} HaplotypeCaller  -I{input.srtbam} -R {ref} -ERC GVCF -G StandardAnnotation -GAS_StandardAnnotation -O {output.ogvcf}" """

rule samtools_vc:
	input: 
		srtbams = expand("{fpath}.srt.rmdup.bam", fpath = fpaths)
	output:
		cmbsvcf = "samtools.cmb.vcf.gz",
		lsf_e = "sm_vc.e",
		lsf_o = "sm_vc.o",
	shell:
		""" bsub -n8 -M20000 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -Jsm_vc -o{output.lsf_o} -e{output.lsf_e} -K "samtools mpileup -t DP,DPR,INFO/DPR -C50 -pm2 -F0.2 –ugf {ref} {input.srtbams} | bcftools call -vmO z -f GQ -o {output.cmbsvcf}" """

rule cmbgvcfs:
	input: 
		vcfs = expand("{fpath}.raw.g.vcf.gz", fpath = fpaths),
	output: 
		o_gvcfs = "cmb.g.vcf",
		#lsf_j = "cmb.job",
		lsf_e = "cmb.e",
		lsf_o = "cmb.o"
	params:
		"-V "+" -V ".join([x+".raw.g.vcf.gz" for x in fpaths])	
	shell:
		""" bsub -n4 -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -Jcmb -o{output.lsf_o} -e{output.lsf_e} -K "java -jar -Xmx4G -XX:ParallelGCThreads=4 {gatk_pth} CombineGVCFs {params} -R {ref} -O {output.o_gvcfs}" """

rule genotypevcf:
	input:
		i_gvcfs = "cmb.g.vcf",
	output:
		o_vcf = "gntp.cmb.g.vcf",
		#lsf_j = "gnv.job",
		lsf_e = "gnv.e",
		lsf_o = "gnv.o"
	shell:	
		""" bsub -M4000 -n4 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -Jgnv -o{output.lsf_o} -K -e{output.lsf_e} "java -jar -Xmx4G -XX:ParallelGCThreads=4 {gatk_pth} GenotypeGVCFs -R {ref} -V {input.i_gvcfs} -O {output.o_vcf}" """	


rule smvc_fltv:
	input: "samtools.cmb.vcf.gz",
	output:
		"samtools.cmb.fltv.vcf.gz"
	shell:
		""" bcftools view -i "INFO/MQB > 0.0001 && QUAL > 30" {input} -O z > {output} """
	
rule gatkv_fltv:
	input:
		"gntp.cmb.g.vcf"
	output:
		"gatk.gntp.cmb.fltv.g.vcf.gz"	
	shell:
		""" bcftools view -i "INFO/InbreedingCoeff > -0.05 && INFO/FS < 20 && QUAL > 300" {input} -O z > {output} """

rule isec_fltv:
	input: "gatk.gntp.cmb.fltv.g.vcf.gz", "samtools.cmb.fltv.vcf.gz",
	output: out
	shell:
		""" bcftools isec -c indels -O z {input[0]} {input[1]} -O v -o {output} """	

#rule fltv: 
	#input:
		#i_vcf = "gntp.cmb.g.vcf"
	#output:
		#out = out,
		#lsf_j = "gnv.job",
		#lsf_e = "fltv.e",
		#lsf_o = "fltv.o"
		
	#shell:
		#""" bsub -K -n4 -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -Jfltv -o{output.lsf_o} -e{output.lsf_e} "java -jar -Xmx4G  -XX:ParallelGCThreads=4 {gatk_pth} VariantFiltration -R {ref} -V {input.i_vcf} -filter 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name "basic_filtering" -O {out}" """ #don't use "QD < 2.0 ......" 

