#uBAM -> fq -> srt.bam -> markduplicate -> gvcf 
# faidx, dict-> gvcfs -> gt -> filter

configfile: "config.json"

uBAMs_l = config["ubams"]
ref = config["ref"]
picard_pth = config["picard_pth"]
gatk_pth = config["gatk_pth"]

f_out = "massoko_fish_flted.vcf"

rule all:
	input: f_out

rule bwa_idx:
	input: ref
	output: 	
		{ref}.amb, {ref}.ann, {ref}.bwt
	shell:
		"bwa index {input}"

rule bam2fq:
	input:
		{dir}/{sample}.bam
	output:
		o1 = {dir}/{sample}_1.fq, 
		o2 = {dir}/{sample}_2.fq
	shell:
		"samtools fastq -1 {o1} -2 {o2} {input}"

rule bwa_map:
	input: 
		i1 = {dir}/{sample}_1.fq,
		i2 = {dir}/{sample}_2.fq
		idx = ref 
	output:
		srtbam = {dir}/{sample}.srt.bam, 
		bami = {dir}/{sample}.srt.bam.bai,
		lsf_j = aln_{sample}.job,
		lsf_o = aln_{sample}.o,
		lsf_e = aln_{sample}.e
	threads: 4
	shell:
	"""	
	bsub -M4000 -n {threads} -R"span[hosts=1] select[mem>4000] rusage[mem=4000]" -J{output.lsf_j} -o{output.lsf_o} -e {output.lsf_e} "bwa mem -t4 {input.idx} {input.i1} {input.i2} | samtools view -h -b - | samtools sort - -o {output.srtbam}"  
	samtools index {output.srtbam}
	"""
rule markdups:
	input: 
		i = {dir}/{sample}.srt.bam,
		picd_pth = {picard_pth}		
	output: 
		o_bam = {dir}/{sample}.srt.rmdup.bam,
		o_mat = {dir}/{sample}.rmdup_matrix.txt,
		lsf_j = md_{sample}.job, #maybe need dir
		lsf_o = md_{sample}.o,
		lsf_e = md_{sample}.e
	shell:
	"""
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -J{output.lsf_j} -o{output.lsf_o} -e{output.lsf_e} "java -jar -Xmx4G {input.picd_pth} MarkDuplicates I={input.i} O={output.o_bam} M={output.o_mat} AS=true"  
	"""
rule fa_idx:
	input: ref
	output: 
		dict_fa = {ref}[:-3].dict,
		fai_fa = {ref}.fai
	shell:
		"""
		samtools faidx {ref}
		samtools dict {ref} -o {output.dict_fa}
		"""	
rule hapcall: # haplotyecaller
	input:
		ref = ref
		dict_fa = {ref}[:-3].dict,
		fai_fa = {ref}.fai,
		bami = {dir}/{sample}.srt.bam.bai,
		srtbam = {dir}/{sample}.srt.bam
	output:
		ogvcf = {dir}/{sample}.raw.g.vcf.gz,
		lsf_j = hc_{sample}.job,
		lsf_o = hc_{sample}.o,
		lsf_e = hc_{sample}.e
	shell:
	"""
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -J{output.lsf_j} -o{output.lsf_o} -e{output.lsf_e} "java -jar -Xmx4G {gatk_pth} HaplotypeCaller  -I{input.srtbam} -R{input.ref} -ERC GVCF -G StandardAnnotation -GAS_StandardAnnotation -O {output.ogvcf}"
	"""
rule cmbgvcfs:
	input: 
		vcfs = expand("{dir}/{sample}.raw.g.vcf.gz", sample = uBAMs_l),
		gatk_pth = gatk_pth,
		ref = ref
	output: 
		o_gvcfs = "massoko.cmb.g.vcf",
		lsf_j = "cmb.job",
		lsf_e = "cmb.e",
		lsf_o = "cmb.o"
	params:
		"-V ".join([x for x in vcfs])	
	shell:
	"""
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -J{output.lsf_j} -o{output.lsf_o} -e{output.lsf_e} "java -jar -Xmx4G {input.gatk_pth} CombineGVCFs {params} -R {input.ref} -O {output.o_gvcfs}" 
	"""
rule genotypevcf:
	input:
		gatk_pth = gatk_pth,
		ref = ref,
		i_gvcfs = "massoko.cmb.g.vcf",
	output:
		o_vcf = "massoko.gntp.cmb.g.vcf",
		lsf_j = "gnv.job",
		lsf_e = "gnv.e",
		lsf_o = "gnv.o"
	shell:	
	"""
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -J{output.lsf_j} -o{output.lsf_o} -e{output.lsf_e} "java -jar -Xmx4G {input.gatk_pth} GenotypeGVCFs -R {input.ref} -V {input.i_gvcfs} -O {output.o_vcf}" 
	"""	
rule fltv: 
	input:
		gatk_pth = gatk_pth,
		ref = ref,
		i_vcf = "massoko.gntp.cmb.g.vcf",
	output:
		f_out = f_out,
		lsf_j = "gnv.job",
		lsf_e = "gnv.e",
		lsf_o = "gnv.o"
		
	shell:
	"""
	bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -J{output.lsf_j} -o{output.lsf_o} -e{output.lsf_e} "java -jar -Xmx4G ${gatk_pth} VariantFiltration -R {ref} -V {i_vcf} -filter "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "basic_filtering" -O {f_out}"
	""" 
