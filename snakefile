import glob
import os
import pandas as pd
from snakemake.io import expand
from snakemake.utils import R
from snakemake.io import glob_wildcards
import re
from os.path import join

configfile: "config.yaml"

data_dir = "/home/suresh/Desktop/raw_data"
samples = glob_wildcards(data_dir + "/{sample}_r1.fastq.gz").sample

rule all:
    input:
    	expand('results/fastqc/{samples}_r1_fastqc.zip', samples=samples),
    	expand('results/fastqc/{samples}_r2_fastqc.zip', samples=samples),
    	expand('results/cutadapt/{samples}_r1_cutadapt.fastq.gz', samples=samples),
        expand('results/cutadapt/{samples}_r2_cutadapt.fastq.gz', samples=samples),
        expand('results/cutadapt/{samples}_r1_cutadapt_fastqc.zip', samples=samples),
        expand('results/cutadapt/{samples}_r2_cutadapt_fastqc.zip', samples=samples),
        expand('results/hisat2/{samples}.cutadapt.sam', samples=samples),
        expand('results/hisat2/{samples}.cutadapt.bam', samples=samples),
        expand('results/stringtie/{samples}/transcript.gtf', samples=samples),
        expand('results/stringtie/{samples}/gene_abundances.tsv', samples=samples),
        expand('results/stringtie/{samples}/cov_ref.gtf', samples=samples),
        expand('results/stringtie/merge_transcripts.gtf'),
        expand('results/ballgown/SigDE.txt')

include: "rules/fastqc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/stringtie.smk"
include: "rules/gtf_merge.smk"
include: "rules/ballgown.smk"

rule clean:
    shell: """
        rm ../.snakemake/log/[0-9]*.snakemake.log
    """
