rule gtf_merge:
    input:
        r1 = expand('results/stringtie/{samples}/transcript.gtf', samples=samples)
    output:
        r1 = 'results/stringtie/merge_transcripts.gtf'
    params:
        gtf= config['ref']['annotation']
    threads: 2
    shell: """
        stringtie -p {threads} --merge  -G {params.gtf} -o {output.r1} {input.r1}
    """