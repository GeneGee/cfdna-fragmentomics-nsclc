#!/usr/bin/env python3
"""
workflow/02_mapping.py

Mapping module for cfDNA fragmentomics pipeline.
Performs: BWA MEM -> samtools sort -> Picard MarkDuplicates -> samtools index
Generates: deduplicated BAM (.dedup.bam), .bai, markduplicates metrics, alignment stats

Usage:
  python workflow/02_mapping.py --config config/config.yaml --sample SAMPLE --r1 path/to/R1.fastq.gz --r2 path/to/R2.fastq.gz --threads 8 --outdir data/results/bam

Supports --dry-run which prints planned commands.
"""

import argparse
import logging
import os
import subprocess
from pathlib import Path

import yaml


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def run_command(cmd, dry_run=False):
    cmd_str = " ".join(map(str, cmd))
    logging.debug(f"Running: {cmd_str}")
    if dry_run:
        print(f"DRY RUN: {cmd_str}")
        return 0, cmd_str
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.debug(res.stdout)
        return res.returncode, res.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd_str}")
        logging.error(e.stderr)
        raise


def main():
    parser = argparse.ArgumentParser(description="Mapping: BWA MEM -> samtools sort -> Picard MarkDuplicates -> index")
    parser.add_argument("--config", required=True, help="Path to config YAML")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--r1", required=True, help="R1 fastq path")
    parser.add_argument("--r2", required=False, default=None, help="R2 fastq path (optional)")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--outdir", required=True, help="Output directory for BAMs and stats")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s: %(message)s")

    cfg = load_config(args.config)
    ref = cfg["reference"]["genome"]

    ensure_dir(args.outdir)
    tmp_prefix = os.path.join(args.outdir, args.sample)

    # Filenames
    sorted_bam = f"{tmp_prefix}.sorted.bam"
    dedup_bam = f"{tmp_prefix}.dedup.bam"
    dedup_bai = f"{dedup_bam}.bai"
    markdup_metrics = f"{tmp_prefix}.markdup.metrics.txt"
    alignment_stats = f"{tmp_prefix}.alignment_stats.txt"

    # Read group
    rg_id = args.sample
    rg_sm = args.sample
    rg_lb = cfg.get("read_group", {}).get("library", "lib1")
    rg_pl = cfg.get("read_group", {}).get("platform", "ILLUMINA")
    rg = f"@RG\\tID:{rg_id}\\tSM:{rg_sm}\\tLB:{rg_lb}\\tPL:{rg_pl}"

    # Build bwa mem command
    bwa_threads = args.threads
    bwa_alg = cfg.get("alignment", {}).get("bwa_algorithm", "mem")
    bwa_opts = cfg.get("alignment", {}).get("bwa_options", "-M")

    bwa_cmd = [
        "bwa",
        bwa_alg,
        bwa_opts,
        "-t",
        str(bwa_threads),
        "-R",
        rg,
        ref,
        args.r1,
    ]
    if args.r2:
        bwa_cmd.append(args.r2)

    # samtools view -> sort
    samtools_sort_cmd = ["samtools", "sort", "-@", str(max(1, args.threads // 2)), "-o", sorted_bam, "-T", f"{tmp_prefix}.tmp"]

    # Picard MarkDuplicates
    picard_cmd = [
        "picard",
        "MarkDuplicates",
        f"I={sorted_bam}",
        f"O={dedup_bam}",
        f"M={markdup_metrics}",
        "CREATE_INDEX=true",
        "VALIDATION_STRINGENCY=LENIENT",
    ]

    # samtools flagstat and idxstats
    flagstat_cmd = ["samtools", "flagstat", dedup_bam]
    idxstats_cmd = ["samtools", "idxstats", dedup_bam]

    # Execute pipeline
    # bwa mem | samtools view -b - | samtools sort -o sorted.bam
    # We'll run bwa mem and pipe to samtools sort via subprocess piping

    logging.info(f"Mapping sample {args.sample}")

    if args.dry_run:
        print("DRY RUN: bwa command:\n", " ".join(bwa_cmd))
        print("DRY RUN: samtools sort command:\n", " ".join(samtools_sort_cmd))
        print("DRY RUN: picard command:\n", " ".join(picard_cmd))
        print("DRY RUN: flagstat:\n", " ".join(flagstat_cmd))
        print("DRY RUN: idxstats:\n", " ".join(idxstats_cmd))
        return

    # Run bwa mem and pipe to samtools sort
    bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
    sort_proc = subprocess.Popen(["samtools", "sort", "-@", str(max(1, args.threads // 2)), "-o", sorted_bam, "-T", f"{tmp_prefix}.tmp"], stdin=bwa_proc.stdout)
    bwa_proc.stdout.close()
    _, bwa_err = bwa_proc.communicate()
    sort_proc.communicate()

    if bwa_proc.returncode not in (0, None):
        logging.error(f"bwa mem failed for sample {args.sample}")
        logging.error(bwa_err.decode() if bwa_err else "")
        raise RuntimeError("bwa mem failed")

    # Check sorted_bam exists
    if not os.path.exists(sorted_bam):
        raise RuntimeError(f"Sorted BAM not found: {sorted_bam}")

    # Run Picard MarkDuplicates
    logging.info("Running Picard MarkDuplicates")
    run_command(picard_cmd)

    # Ensure dedup bam exists
    if not os.path.exists(dedup_bam):
        raise RuntimeError(f"Deduplicated BAM not produced: {dedup_bam}")

    # Index produced by Picard if CREATE_INDEX=true; if not, run samtools index
    if not os.path.exists(dedup_bai):
        logging.info("Index not found after Picard, running samtools index")
        run_command(["samtools", "index", dedup_bam])

    # Run flagstat and idxstats and write to alignment_stats
    with open(alignment_stats, "w") as outf:
        logging.info("Computing samtools flagstat")
        res = subprocess.run(flagstat_cmd, check=True, capture_output=True, text=True)
        outf.write("# samtools flagstat\n")
        outf.write(res.stdout)
        logging.info("Computing samtools idxstats")
        res2 = subprocess.run(idxstats_cmd, check=True, capture_output=True, text=True)
        outf.write("\n# samtools idxstats\n")
        outf.write(res2.stdout)

    logging.info(f"Mapping completed for {args.sample}. Outputs: {dedup_bam}, {markdup_metrics}, {alignment_stats}")


if __name__ == "__main__":
    main()
