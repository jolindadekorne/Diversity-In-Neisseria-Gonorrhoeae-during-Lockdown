IDS, = glob_wildcards("output_fastp/{id}_out1.fastq.gz")
CONT, = glob_wildcards("contaminated/{id}_out1.fastq.gz")
ALL_SAMPLES = IDS + CONT

import pandas as pd

f = open("contaminated_isolates.txt")
lines = f.readlines()
f.close()

contaminated_isolates = {}

for line in lines:
  isolate = line.rstrip('\n')
  contaminated_isolates[isolate] = isolate

print(contaminated_isolates)

rule all:
	input:
                expand("output_fastp/{sample}_fastp.json", sample = IDS),
                "multiqc_fastp_out",
		expand("coverage_FA1090/{sample}_cov.txt", sample = IDS),
		expand("output_fastp/{sample}output_spades", sample = IDS),
		expand("quast_out/{sample}", sample = IDS),
		"multiqc_quast_out",
                expand("kraken2_out/{sample}_kraken2_report.txt", sample=contaminated_isolates),
		expand("kraken2_out/filtered_reads/{sample}_kraken_extracted_1.fastq", sample = contaminated_isolates),
                expand("kraken2_out/{sample}output_spades", sample = CONT),
		expand("snippy_refFA1090_out/{sample}", sample = IDS),
		expand("snippy_refFA1090_out/{sample}", sample = contaminated_isolates),
		"snippy_refFA1090_out/clean.full.aln",
                "maskrc_gubbins_FA1090_out/masked_gubbins_snippyrefFA1090.aln",
		"snp_dists/snp_dists_masked_gubbins_FA1090.tsv"

rule fastp:
        input:
                fw = "raw_data/{sample}_R1.fastq.gz",
                rv = "raw_data/{sample}_R2.fastq.gz"
        output:
                fw = "output_fastp/{sample}_out1.fastq.gz",
                rv = "output_fastp/{sample}_out2.fastq.gz",
                json = "output_fastp/{sample}_fastp.json",
                html = "output_fastp/{sample}_fastp.html"
	conda: 
		"envs/fastp.yml"
	params:
		compression_level = 9
	log:
		"logs/fastp/fastp_{sample}.log"
	threads: 4
	shell:
		"""
		fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} --html {output.html} --json {output.json} 2>&1>{log}
		"""

rule multiqc_fastp:
	input:
		expand("output_fastp/{sample}_fastp.json", sample=IDS)
	output:
		directory("multiqc_fastp_out")
	conda: 
		"envs/multiqc.yml"
	log:
		"logs/multiqc/multiqc_fastp.log"
	shell:
		"""
		mkdir -p {output}
		multiqc {input} --outdir {output} 2>&1>{log}
		"""

rule coverage:
        input:
                fw = "output_fastp/{sample}_out1.fastq.gz",
                rv = "output_fastp/{sample}_out2.fastq.gz"
        output:
                cov = "coverage_FA1090/{sample}_cov.txt"
        log:
                bam = "logs/bwamem_samtools/{sample}.log",
                cov = "logs/coverage_FA1090/coverage_FA1090_{sample}.log"
        conda:
                "envs/samtools.yml"
        params:
                ref = "ref/NgRefFA1090",
                path = "/home/jdkorne/samtools_1.11/bin",
                bam_out = "temp_bam/{sample}_temp_sorted.bam"
        threads: 6
        shell:
                """
                mkdir -p temp_bam
                mkdir -p coverage_FA1090
		/home/jdkorne/bwa-mem2/bwa-mem2 index -p {params.ref} {params.ref}.fna
		/home/jdkorne/bwa-mem2/bwa-mem2 mem -t {threads} {params.ref} {input.fw} {input.rv} | {params.path}/samtools sort --threads {threads} -o {params.bam_out} 2>&1>{log.bam}
                {params.path}/samtools coverage -o {output.cov} {params.bam_out} 2>&1>{log.cov}
                rm {params.bam_out}
                """

rule spades:
	input:
                fw = "output_fastp/{sample}_out1.fastq.gz",
                rv = "output_fastp/{sample}_out2.fastq.gz"
        output:
                dir = "output_fastp/{sample}output_spades"
        log:
                "logs/spades/spades_{sample}.log"
        conda:
                "envs/spades.yml"
        threads: 6
        shell:
                """
		spades.py --isolate -1 {input.fw} -2 {input.rv} -o {output.dir}
                """
		
rule quast:
        input:
                assembly = "output_fastp/{sample}output_spades/scaffolds.fasta"
        output:
                directory("quast_out/{sample}")
        conda:
                "envs/quast.yml"
        log:
                "logs/quast/quast_{sample}.log"
        shell:
                """
                quast.py -o {output} {input.assembly}
                """

rule multiqc_quast:
        input:
                expand("quast_out/{sample}", sample=IDS)
        output:
                directory("multiqc_quast_out")
        params:
                input = expand("quast_out/{sample}/report.tsv", sample=IDS)
        conda:
                "envs/multiqc.yml"
        log:
                "logs/multiqc_quast/multiqc_quast.log"
        shell:
                """
                multiqc {params.input} -o {output} 2>&1>{log}
                """

rule kraken2:
	input:
		fw = "contaminated/{sample}_out1.fastq.gz",
                rv = "contaminated/{sample}_out2.fastq.gz"
	output:
		report = "kraken2_out/{sample}_kraken2_report.txt",
		output = "kraken2_out/{sample}_kraken2_output.txt"
	conda:
		"envs/kraken2.yml"
	params:
		db = "/home/jdkorne/kraken2db/minikraken2_v1_8GB/"
	log:
		"logs/kraken2/{sample}.log"
	threads: 8
	shell:
		"""
		kraken2 --db {params.db} --gzip-compressed --paired --threads {threads} --report {output.report} --output {output.output} {input.fw} {input.rv}
		"""

rule kraken_extract:
	input:
                fw = "contaminated/{sample}_out1.fastq.gz",
                rv = "contaminated/{sample}_out2.fastq.gz",
                report = "kraken2_out/{sample}_kraken2_report.txt",
                output = "kraken2_out/{sample}_kraken2_output.txt"
	output:
		fw = "kraken2_out/filtered_reads/{sample}_kraken_extracted_1.fastq",
		rv = "kraken2_out/filtered_reads/{sample}_kraken_extracted_2.fastq"
	params:
		script = "kraken2_out/extract_kraken_reads.py",
		tax = 485
	log:
		"logs/kraken2/{sample}_extraction.log"
	shell:
		"""	
		python3 {params.script} -k {input.output} -s1 {input.fw} -s2 {input.rv} --fastq-output -o {output.fw} -o2 {output.rv} -t {params.tax} --include-children -r {input.report}
		"""

rule spades_cont:
	input:
		fw = "kraken2_out/filtered_reads/{sample}_kraken_extracted_1.fastq",
		rv = "kraken2_out/filtered_reads/{sample}_kraken_extracted_2.fastq"
        output:
                dir = "kraken2_out/{sample}output_spades"
        log:
                "logs/spades_cont/spades_{sample}.log"
        conda:
                "envs/spades.yml"
        threads: 6
        shell:
                """
		spades.py --isolate -1 {input.fw} -2 {input.rv} -o {output.dir}
                """

rule snippy_FA1090:
	input:
		fw = "output_fastp/{sample}_out1.fastq.gz",
		rv = "output_fastp/{sample}_out2.fastq.gz"
	output:
		directory("snippy_refFA1090_out/{sample}")
	conda:
		"envs/snippy.yml"
	params:
		outdir = directory("snippy_refFA1090_out"),
		ref = "ref/NgRefFA1090.gbk"
	log:    
		"logs/snippy_refFA1090/snippy_refFA1090_{sample}.log"
	threads: 16
	shell:
		"""
		mkdir -p {params.outdir}
		snippy --cpus {threads} --ref {params.ref} --R1 {input.fw} --R2 {input.rv} --outdir {output} 2>&1>{log}
		"""

rule snippy_FA1090_cont:
        input:
                fw = "kraken2_out/filtered_reads/{sample}_kraken_extracted_1.fastq.gz",
                rv = "kraken2_out/filtered_reads/{sample}_kraken_extracted_2.fastq.gz"
        output:
                directory("snippy_refFA1090_out/{sample}")
        conda:
                "envs/snippy.yml"
        params:
                outdir = directory("snippy_refFA1090_out"),
                ref = "ref/NgRefFA1090.gbk"
        log:
                "logs/snippy_refFA1090/snippy_refFA1090_{sample}.log"
        threads: 16
        shell:
                """
                mkdir -p {params.outdir}
                snippy --cpus {threads} --ref {params.ref} --R1 {input.fw} --R2 {input.rv} --outdir {output} 2>&1>{log}
                """

rule snippy_aln_FA1090:
        input:
                corein = expand("snippy_refFA1090_out/{sample}", sample = ALL_SAMPLES)
        output:
                cleanout = "snippy_refFA1090_out/clean.full.aln"
        conda:
                "envs/snippy.yml"
        params:
                outdir = directory("snippy_refFA1090_out"),
                ref = "ref/NgRefFA1090.gbk",
                prefix = "core",
                cleanin = "core.full.aln"
        shell:
                """
                snippy-core {input.corein} --ref {params.ref} --prefix {params.prefix}
                snippy-clean_full_aln {params.cleanin} > {output.cleanout}
                """

rule gubbins:
	input:
		snippy_FA1090 = "snippy_refFA1090_out/clean.full.aln"
	output:
		FA1090 = "maskrc_gubbins_FA1090_out/masked_gubbins_snippyrefFA1090.aln"
	conda:
		"envs/gubbins.yml"
	params:
		gubbins_out_FA1090 = "gubbins_snippyrefFA1090_out",
                maskrc_out_FA1090 = "maskrc_gubbins_FA1090_out",
                model = "GTRGAMMA",
                prefix_FA1090 = "gubbins_snippyrefFA1090"
	log:    
		gubbins_FA1090 = "logs/gubbins/gubbins_snippyrefFA1090.log",
                maskrc_FA1090 = "logs/maskrc_gubbins_snippyrefFA1090.log"
	shell:
                """
                mkdir -p {params.gubbins_out_FA1090} {params.maskrc_out_FA1090}
		run_gubbins.py {input.snippy_FA1090} --prefix {params.prefix_FA1090} --raxml_model {params.model} 2>&1>{log.gubbins_FA1090}
                python3 script/maskrc-svg.py --gubbins --aln {input.snippy_FA1090} --out {output.FA1090} {params.prefix_FA1090} 2>{log.maskrc_FA1090}
                mv {params.prefix_FA1090}.* {params.gubbins_out_FA1090}
                """

rule snp_dists_snippy:
        input:
                snippy_FA1090 = "snippy_refFA1090_out/clean.full.aln",
                masked_gubbins_FA1090 = "maskrc_gubbins_FA1090_out/masked_gubbins_snippyrefFA1090.aln"
        output:
                dists_snippy_FA1090 = "snp_dists/snp_dists_snippy_FA1090.tsv",
                dists_masked_gubbins_FA1090 = "snp_dists/snp_dists_masked_gubbins_FA1090.tsv"
        params:
                dir = "snp_dists"
        conda:
                "envs/snp_dists.yml"
        shell:
                """
                mkdir -p {params.dir}
                snp-dists -m {input.snippy_FA1090}>{output.dists_snippy_FA1090}
                snp-dists -m {input.masked_gubbins_FA1090}>{output.dists_masked_gubbins_FA1090}
                """

