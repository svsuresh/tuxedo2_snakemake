This repository contains snakemake script for analyzing RNAseq samples from Griffith lab [tutorial](https://github.com/griffithlab/rnaseq_tutorial/wiki) and data used in the tutorial is available [here]((http://genomedata.org/rnaseq-tutorial/practical.tar)). This script follows Tuxedo2 protocol (HISAT2-STRINGTIE-BALLGOWN). Workflow has following steps:
* Quality check by fastqc
* Adapter trimming by cutadapt
* Quality check post adapter trimming by fastqc
* Alignment with HISAT2
* Transcript quanitification with Stringtie
* Statistical tests and visualization with Ballgown (in R)

[Tuxedo Workflow](https://drive.google.com/uc?id=1TIwsrxA3w64SoJYh6OasGBo8iVYQMiQf)

Please **note** that **this script is not production ready**. This is for beginners in using snakemake. For generating graph, snakefile should be in sentence case.


