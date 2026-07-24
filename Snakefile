rule fragments_sample:
    input:
        bam="data/results/bam/{sample}.dedup.bam",
    output:
        metrics="data/results/fragments/{sample}_fragment_metrics.csv",
        bed="data/results/fragments/{sample}_fragments.bed",
        plot="data/results/fragments/{sample}_length_distribution.png",
    threads: 2
    params:
        outdir="data/results/fragments"
    shell:
        "python workflow/03_fragmentation_analysis.py --config config/config.yaml --bam {input.bam} --sample {wildcards.sample} --outdir {params.outdir}"
