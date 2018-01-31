This repository contains snakemake script for analyzing RNAseq samples from Griffith lab [tutorial](http://genomedata.org/rnaseq-tutorial/practical.tar). This script follows Tuxedo2 protocol (HISAT2-STRINGTIE-BALLGOWN). Workflow has following steps:
* Quality check by fastqc
* Adapter trimming by cutadapt
* Quality check post adapter trimming by fastqc
* Alignment with HISAT2
* Transcript quanitification with Stringtie
* Statistical tests and visualization with Ballgown (in R)

Please **note** that **this script is not production ready**. This is for beginners in using snakemake.
