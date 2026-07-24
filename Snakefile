configfile: "config/config.yaml"

import csv

sample_sheet = config.get("sample_sheet", "config/sample_sheet.csv")
SAMPLES = []
with open(sample_sheet) as fh:
    for row in csv.DictReader(fh):
        SAMPLES.append(row["sample_id"])

rule all:
    input:
        expand("data/results/bam/{sample}.dedup.bam", sample=SAMPLES),
        "data/results/qc/qc_summary.csv"

rule run_qc:
    input:
        sample_sheet
    output:
        "data/results/qc/qc_summary.csv"
    threads: 4
    shell:
        "python workflow/01_qc.py --config config/config.yaml --sample-sheet {input} --workers {threads}"

rule map_sample:
    input:
        r1="data/results/trimmed_fastq/{sample}_R1.trimmed.fastq.gz",
        r2="data/results/trimmed_fastq/{sample}_R2.trimmed.fastq.gz",
    output:
        bam="data/results/bam/{sample}.dedup.bam",
        bai="data/results/bam/{sample}.dedup.bam.bai",
    threads: 8
    params:
        outdir="data/results/bam"
    shell:
        "python workflow/02_mapping.py --config config/config.yaml --sample {wildcards.sample} --r1 {input.r1} --r2 {input.r2} --threads {threads} --outdir {params.outdir}"
