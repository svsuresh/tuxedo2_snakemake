rule ballgown:
    input: 'results/stringtie/merge_transcripts.gtf'
    output: 'results/ballgown/SigDE.txt'
    script:
        "scripts/ballgown.R"