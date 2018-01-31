rule cutadapt:
    input:
        r1 = data_dir+'/{samples}1.fastq.gz',
        r2 = data_dir+'/{samples}2.fastq.gz'

    output:
        r1 = 'results/cutadapt/{samples}1_cutadapt.fastq.gz',
        r2 = 'results/cutadapt/{samples}2_cutadapt.fastq.gz'
    message:
        "--- running cutadapt ---"
    threads: 8
    shell: """
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT  -o {output.r1} -p {output.r2} {input.r1} {input.r2}
    """