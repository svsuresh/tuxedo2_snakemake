rule hisat2:
    input:
        r1 = 'results/cutadapt/{samples}_r1_cutadapt.fastq.gz',
        r2 = 'results/cutadapt/{samples}_r2_cutadapt.fastq.gz'
    output:
        r1 = 'results/hisat2/{samples}.cutadapt.sam'
    message:
        "------aligning with hisat2....wait.."
    params:
        index= config['ref']['index']
    threads: 4
    shell: """
        hisat2  -p {threads} -x {params.index} --dta  --rna-strandness RF -1 {input.r1}  -2 {input.r2} -S {output.r1}
    """

rule create_bams:
    input:
        r1 = 'results/hisat2/{samples}.cutadapt.sam'
    output:
        r1 = 'results/hisat2/{samples}.cutadapt.bam'
    message:
        "---covnering sam to bam  and indexing the bam files"
    shell: """
    samtools view -bh {input.r1} | samtools sort - -o {output.r1}; samtools index {output.r1}
    """