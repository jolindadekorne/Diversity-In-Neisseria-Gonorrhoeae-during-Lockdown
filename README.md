# Diversity In Neisseria Gonorrhoeae during Lockdown (DINGL)

This repository contains the Snakemake pipeline used for the study on genetic diversity in the _Neisseria gonorrhoeae_ population during the first COVID-19 lockdown in Amsterdam, the Netherlands. 

## Dependencies
This pipeline uses the following dependencies:
- Conda
- Snakemake
- Python3

## Input
This pipeline uses forward and reverse raw Illumina sequencing reads which are located in the folder `raw_data`. The raw data files should be named `{id}_R1.fastq.gz` and `{id}_R2.fastq.gz`. 

## Other files needed
- [maskrc-svg script](https://github.com/kwongj/maskrc-svg); expected path: 'scripts/maskrc-svg.py'
- Reference genome FA1090 is used for calculating coverage and calling variants: NC_002946.2. The reference genome should be located in the same directory as the Snakefile.

## Pipeline 
The pipeline includes the following steps and tools:

| Step     | Tool     |   
| ---------|----------|
| Filter low quality raw reads + trim adapters | fastp |
| Assembly | SPAdes | 
| Assembly quality check | QUAST |
| Read mapping against reference genome FA1090 | BWA-MEM2 |
| Calculation of percentage of bases covered and coverage depth | samtools |
| Call variants using reference FA1090 + create core genome alignment | snippy |
| Remove recombination from variant alignment | Gubbins |
| Mask the recombination sites in the core genome alignment | [maskrc-svg script](https://github.com/kwongj/maskrc-svg) |
| Calculate recombination filtered- and unfiltered SNP distances | snp-dists |
