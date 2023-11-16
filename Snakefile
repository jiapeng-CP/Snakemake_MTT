import pandas as pd

samplesDF = pd.read_csv("samplesheet.tsv", sep="\t", index_col="Sample")
SAMPLES = []
sampleFq1 = {}
sampleFq2 = {}


for x in samplesDF.index:
	SAMPLES.append(x)
	sampleFq1[x] = samplesDF.loc[x,'fq1']
	sampleFq2[x] = samplesDF.loc[x,'fq2']

def get_sample(wildcards):
	return(wildcards.sample)

def get_fastq1(wildcards):
	return(samplesDF.loc[wildcards.sample,'fq1'])

def get_fastq2(wildcards):
	return(samplesDF.loc[wildcards.sample,'fq2'])

rule all:
	input:
		expand("map2HOMD/{sample}.json", sample=SAMPLES),
		#expand("map2human/{sample}.json", sample=SAMPLES)
        

rule fastpBfQC:
	input:
		r1 = get_fastq1,
		r2 = get_fastq2

	output:
		r1 = "fastp/{sample}.r1.fq.gz",
		r2 = "fastp/{sample}.r2.fq.gz",
		json = "fastp/{sample}.json",
		html = "fastp/{sample}.html"

	log: "logs/{sample}.fastp.log"
	threads: 8
	params:
		sn = get_sample
	shell:
		"mkdir -p fastp \n"
		"/home/jiapengc/.conda/envs/QC/bin/fastp --in1 {input.r1} "
		"--in2 {input.r2} "
		"--out1 {output.r1} "
		"--out2 {output.r2} "
		"--json {output.BfQC.json} "
		"--html {output.BfQC.html} "
		"--thread 8"
        

rule rRNArm:
	input:
		r1 = "fastp/{sample}.r1.fq.gz",
		r2 = "fastp/{sample}.r2.fq.gz"

	output:
		sam = "map2rRNA/{sample}.sam",
		fq1 = "map2rRNA/{sample}_rRNA_rm.fastq.1.gz",
		fq2 = "map2rRNA/{sample}_rRNA_rm.fastq.2.gz"

	threads: 8

	params: sp = get_sample

	shell:
		"mkdir -p map2rRNA \n"
		#"/home/jiapengc/.conda/envs/biobakery3/bin/bowtie2 -x /home/jiapengc/db/rRNA/rRNA.rfam.silva "
		"/home/jiapengc/.conda/envs/biobakery3/bin/bowtie2 -x /home/jiapengc/db/SILVA_128_LSUParc_SSUParc_ribosomal_RNA/SILVA_128_LSUParc_SSUParc_ribosomal_RNA "
		"-1 {input.r1} -2 {input.r2} "
		"-S {output.sam} "
		"--sensitive --threads 8 "
		"--un-conc-gz map2rRNA/{params.sp}_rRNA_rm.fastq.gz"



rule rRNAstat:
	input:
		sam = "map2rRNA/{sample}.sam"
	output:
		json = "map2rRNA/{sample}.json"
	threads: 4
	params: sp = get_sample
	shell:
		"/home/jiapengc/.conda/envs/QC/bin/samtools view -bhS --threads 4 {input.sam} > map2rRNA/{params.sp}.bam \n"
		"/home/jiapengc/bin/bamstats --cpu 4 --input map2rRNA/{params.sp}.bam > {output.json}"



rule rmHuman:
	input:
		r1 = "map2rRNA/{sample}_rRNA_rm.fastq.1.gz",
		r2 = "map2rRNA/{sample}_rRNA_rm.fastq.2.gz"

	output:
		fq1 = "map2human/{sample}_human_rm.fastq.1.gz",
		fq2 = "map2human/{sample}_human_rm.fastq.2.gz",
		sam = "map2human/{sample}.sam"
	threads: 8
	params: sp = get_sample

	shell:
		"mkdir -p map2human \n"
		"/home/jiapengc/.conda/envs/biobakery3/bin/bowtie2 -x /home/jiapengc/db/chm13.draft_v1.0_plusY/chm13.draft_v1.0_plusY "
		"-1 {input.r1} -2 {input.r2} "
		"-S {output.sam} "
		"--sensitive --threads 8 "
		"--un-conc-gz map2human/{params.sp}_human_rm.fastq.gz"



rule humanStat:
	input:
		sam = "map2human/{sample}.sam"
	output:
		json = "map2human/{sample}.json"
	threads: 4
	params: sp = get_sample
	shell:
		"/home/jiapengc/.conda/envs/QC/bin/samtools view -bhS --threads 4 {input.sam} > map2human/{params.sp}.bam \n"
		"/home/jiapengc/bin/bamstats --cpu 4 --input map2human/{params.sp}.bam > {output.json}"


rule fastpAfQC:
	input:
		r1 = "map2human/{sample}_human_rm.fastq.1.gz",
		r2 = "map2human/{sample}_human_rm.fastq.2.gz"

	output: #https://github.com/OpenGene/fastp/issues/164
		json = "fastp/{sample}.AfQC.json",
		html = "fastp/{sample}.AfQC.html"

	log: "logs/{sample}.fastp.log"
	threads: 8
	params:
		sn = get_sample
	shell:
		"mkdir -p fastp \n"
		"/home/jiapengc/.conda/envs/QC/bin/fastp "
		"--in1 {input.r1} "
		"--in2 {input.r2} "
		"--json {output.json} "
		"--html {output.html} "
		"--thread 8"


rule map2HOMD:
	input:
		fq1 = "map2human/{sample}_human_rm.fastq.1.gz",
		fq2 = "map2human/{sample}_human_rm.fastq.2.gz"

	output:
		sam = "map2HOMD/{sample}.sam"

	threads: 8
	params: sp = get_sample

	shell:
		"mkdir -p map2HOMD \n"
		"/home/jiapengc/.conda/envs/biobakery3/bin/bowtie2 -x /home/jiapengc/db/HOMD/ALL_genomes.fna "
		"-1 {input.fq1} -2 {input.fq2} "
		"-S {output.sam} "
		"--sensitive --threads 8 "



rule HOMDstat:
	input:
		sam = "map2HOMD/{sample}.sam"
	output:
		json = "map2HOMD/{sample}.json"
	threads: 4
	params: sp = get_sample
	shell:
		"/home/jiapengc/.conda/envs/QC/bin/samtools view -bhS --threads 4 {input.sam} > map2HOMD/{params.sp}.bam \n"
		"/home/jiapengc/bin/bamstats --cpu 4 --input map2HOMD/{params.sp}.bam > {output.json} \n"
		#"samtools sort map2HOMD/{params.sp}.bam -o map2HOMD/{params.sp}.s.bam "
		#"/home/artemisl/.conda/envs/biobakery/bin/samtools coverage map2HOMD/{params.sp}.s.bam > {params.sp}.ref.coverage "
