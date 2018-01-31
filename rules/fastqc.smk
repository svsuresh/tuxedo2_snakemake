rule fastqc:
    input:
        r1 = data_dir+'/{samples}_r1.fastq.gz',
        r2 = data_dir+'/{samples}_r2.fastq.gz'

    output:
        r1 = 'results/fastqc/{samples}_r1_fastqc.zip',
        r2 = 'results/fastqc/{samples}_r2_fastqc.zip'
    priority: 50
    threads: 4
    message:
        "--- running fastqc ---"
    shell: """
        fastqc  {input.r1} {input.r2} -q -f fastq -o results/fastqc/
    """

rule fastqc_after:
    input:
        r1 = 'results/cutadapt/{samples}_r1_cutadapt.fastq.gz',
        r2 = 'results/cutadapt/{samples}_r2_cutadapt.fastq.gz'

    output:
        r1 = 'results/cutadapt/{samples}_r1_cutadapt_fastqc.zip',
        r2 = 'results/cutadapt/{samples}_r2_cutadapt_fastqc.zip'
    threads: 4
    message:
        "--- running fastqc again ---"
    shell: """
        fastqc  {input.r1} {input.r2} -q -f fastq -o results/cutadapt/
    """