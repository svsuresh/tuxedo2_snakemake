rule stringtie:
    input:
        r1 = 'results/hisat2/{samples}.cutadapt.bam'
    output:
        r1 = 'results/stringtie/{samples}/transcript.gtf',
        r2 = 'results/stringtie/{samples}/gene_abundances.tsv',
        r3 = 'results/stringtie/{samples}/cov_ref.gtf'
    params:
        gtf= config['ref']['annotation']

    shell: """
    stringtie -G {params.gtf} --rf  -e -B -o {output.r1} -A {output.r2} -C {output.r3} --rf {input.r1}
    """