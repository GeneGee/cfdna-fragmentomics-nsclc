#!/usr/bin/env python3
"""
workflow/03_fragmentation_analysis.py

Fragmentation analysis module:
- Input: deduplicated BAM file (paired-end)
- Outputs:
  - fragments BED: {sample}_fragments.bed (chrom, start, end, sample:len:gc)
  - fragment_metrics.csv: summary metrics (total fragments, length histogram peaks, GC mean)
  - length_distribution.png

Notes:
- For paired-end BAM, this script iterates read1 records (is_read1) and uses template_length (tlen) to compute fragment coordinates.
- It requires pysam and matplotlib (provided in environment.yml).

Usage:
  python workflow/03_fragmentation_analysis.py --config config/config.yaml --bam data/results/bam/NSCLC_001.dedup.bam --sample NSCLC_001 --outdir data/results/fragments

Supports --dry-run which prints planned actions.
"""

import argparse
import logging
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import yaml
from matplotlib import pyplot as plt


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def compute_gc(seq: str) -> float:
    if not seq:
        return np.nan
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    atgc = sum(seq.count(x) for x in ("A", "T", "G", "C"))
    return (g + c) / atgc if atgc > 0 else np.nan


def process_bam(bam_path, sample_id, cfg, outdir, dry_run=False, max_fragments=None):
    outdir = Path(outdir)
    ensure_dir(outdir)

    fragments_bed = outdir / f"{sample_id}_fragments.bed"
    metrics_csv = outdir / f"{sample_id}_fragment_metrics.csv"
    length_png = outdir / f"{sample_id}_length_distribution.png"

    ref = cfg["reference"]["genome"]

    planned = {
        "bam": str(bam_path),
        "fragments_bed": str(fragments_bed),
        "metrics_csv": str(metrics_csv),
        "length_png": str(length_png),
        "reference": ref,
    }

    if dry_run:
        print("DRY RUN: would generate the following files:")
        for k, v in planned.items():
            print(f"  {k}: {v}")
        return planned

    logging.info(f"Opening BAM: {bam_path}")
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    fasta = pysam.FastaFile(ref)

    counts = 0
    lengths = []
    gc_list = []
    records_out = []

    # Iterate primary reads and take read1 of properly paired reads
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if not read.is_proper_pair:
            continue
        # only process read1 to avoid double counting
        if not read.is_read1:
            continue

        tlen = read.template_length
        if tlen == 0:
            # fallback: use mate positions if available
            try:
                mate = bam.mate(read)
                start = min(read.reference_start, mate.reference_start)
                end = max(read.reference_end, mate.reference_end)
                length = end - start
            except Exception:
                continue
        else:
            # compute start and end using tlen; tlen may be negative depending on orientation
            start = min(read.reference_start, read.reference_start + tlen)
            end = max(read.reference_start, read.reference_start + tlen)
            length = end - start

        if length <= 0:
            continue

        # limit for performance/testing
        if max_fragments is not None and counts >= max_fragments:
            break

        # fetch sequence for GC
        try:
            seq = fasta.fetch(read.reference_name, start, end)
        except Exception:
            seq = ""

        gc = compute_gc(seq)

        lengths.append(length)
        gc_list.append(gc)
        records_out.append((read.reference_name, start, end, f"{sample_id};len={length};gc={gc:.3f}"))
        counts += 1

    bam.close()
    fasta.close()

    if counts == 0:
        logging.warning(f"No fragments extracted from {bam_path}")

    # write fragments bed
    logging.info(f"Writing fragments BED: {fragments_bed}")
    with open(fragments_bed, "w") as fh:
        for chrom, start, end, name in records_out:
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    # compute metrics
    arr = np.array(lengths) if lengths else np.array([])
    metrics = {}
    metrics["sample_id"] = sample_id
    metrics["n_fragments"] = int(len(arr))
    metrics["mean_length"] = float(np.nanmean(arr)) if arr.size else np.nan
    metrics["median_length"] = float(np.nanmedian(arr)) if arr.size else np.nan
    metrics["std_length"] = float(np.nanstd(arr)) if arr.size else np.nan
    metrics["mean_gc"] = float(np.nanmean(gc_list)) if gc_list else np.nan

    # compute fragment ranges from config
    frag_ranges = cfg.get("fragmentation", {}).get("fragment_ranges", {})
    range_counts = {}
    for name, (lo, hi) in frag_ranges.items():
        # hi may be null -> treat as large
        hi_val = hi if hi is not None else 1_000_000
        cnt = int(((arr >= lo) & (arr < hi_val)).sum()) if arr.size else 0
        range_counts[f"range_{name}"] = cnt
    metrics.update(range_counts)

    # nucleosome peak detection (120-180)
    if arr.size:
        bins = np.arange(20, 1000, 1)
        hist, edges = np.histogram(arr, bins=bins)
        # find peak between nucleosome_min and nucleosome_max
        nm = cfg.get("fragmentation", {}).get("nucleosome_min", 120)
        nx = cfg.get("fragmentation", {}).get("nucleosome_max", 180)
        mask = (edges[:-1] >= nm) & (edges[:-1] <= nx)
        if mask.any():
            peak_idx = np.argmax(hist * mask)
            peak_length = int(edges[:-1][peak_idx])
            metrics["nucleosome_peak_length"] = int(peak_length)
            metrics["nucleosome_peak_count"] = int(hist[peak_idx])
        else:
            metrics["nucleosome_peak_length"] = np.nan
            metrics["nucleosome_peak_count"] = 0
    else:
        metrics["nucleosome_peak_length"] = np.nan
        metrics["nucleosome_peak_count"] = 0

    # write metrics
    logging.info(f"Writing metrics CSV: {metrics_csv}")
    dfm = pd.DataFrame([metrics])
    dfm.to_csv(metrics_csv, index=False)

    # plot length distribution
    if arr.size:
        plt.figure(figsize=(8, 4))
        plt.hist(arr, bins=range(20, 500, 1), color="#4C72B0")
        plt.xlabel("Fragment length (bp)")
        plt.ylabel("Count")
        plt.title(f"Fragment length distribution: {sample_id}")
        plt.axvspan(cfg.get("fragmentation", {}).get("nucleosome_min", 120), cfg.get("fragmentation", {}).get("nucleosome_max", 180), color="orange", alpha=0.3, label="nucleosome range")
        plt.legend()
        plt.tight_layout()
        plt.savefig(length_png)
        plt.close()

    logging.info(f"Fragmentation analysis completed for {sample_id}: wrote {counts} fragments")
    return planned


def main():
    parser = argparse.ArgumentParser(description="Fragmentation analysis: extract fragments, length dist, GC, nucleosome peak")
    parser.add_argument("--config", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--max-fragments", type=int, default=None, help="For testing: limit number of fragments to process")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s: %(message)s")
    cfg = load_config(args.config)
    process_bam(args.bam, args.sample, cfg, args.outdir, dry_run=args.dry_run, max_fragments=args.max_fragments)


if __name__ == "__main__":
    main()
