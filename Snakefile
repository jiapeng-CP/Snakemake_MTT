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
		expand("map2rRNA/{sample}.sam", sample=SAMPLES)
        

rule fastp:
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
		"--json {output.json} "
		"--html {output.html} "
		"--thread 8"
        

rule rRNArm:
	input:
		r1 = "fastp/{sample}.r1.fq.gz",
		r2 = "fastp/{sample}.r2.fq.gz"

	output:
		sam = "map2rRNA/{sample}.sam",
		fq = "map2rRNA/{sample}_rRNA_rm.fastq.gz",
		fq1 = "map2rRNA/{sample}_rRNA_rm.fastq.1.gz"

	threads: 8

	shell:
		"mkdir -p map2rRNA \n"
		"/home/jiapengc/.conda/envs/biobakery3/bin/bowtie2 -x /home/jiapengc/db/rRNA/rRNA.rfam.silva "
		"-1 {input.r1} -2 {input.r2} "
		"-S {output.sam} "
		"--sensitive --threads 8 "
		"--un-conc-gz {output.fq}"

#/home/jiapengc/.conda/envs/QC/bin/samtools view -bhS --threads 5 out.sam > out.bam
#/home/jiapengc/bin/bamstats --cpu 5 --input out.bam > bamstat.json