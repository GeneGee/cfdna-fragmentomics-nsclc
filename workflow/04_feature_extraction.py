#!/usr/bin/env python3
"""
workflow/04_feature_extraction.py

Extended feature extraction from fragment BEDs and fragment metrics.
Features include (first version):
 - total_fragments, mean/median/std length, mean_gc
 - proportion/counts in fragment_ranges
 - nucleosome peak features
 - endpoint k-mer preferences (k configurable)
 - promoter / gene body enrichment (requires GTF in config)
 - CpG island overlap counts (requires cpg_bed in config)
 - chromosome-level counts (per-chromosome normalized counts)
 - Normalization options: z-score (default), quantile, RPM

Inputs:
 - fragments_dir: directory with {sample}_fragments.bed and {sample}_fragment_metrics.csv
 - sample_sheet: used to list samples and groups

Outputs:
 - CSV feature matrix at --out

Note: uses pysam, pybedtools. If GTF or CpG bed not provided in config, corresponding features are skipped.
"""

import argparse
import logging
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import yaml
from pybedtools import BedTool
from sklearn.preprocessing import quantile_transform
from scipy import stats


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def read_fragment_metrics(metrics_path):
    if not os.path.exists(metrics_path):
        return None
    try:
        df = pd.read_csv(metrics_path)
        return df.iloc[0].to_dict()
    except Exception:
        return None


def read_fragments_bed(bed_path):
    if not os.path.exists(bed_path):
        return None
    try:
        bt = BedTool(bed_path)
        return bt
    except Exception:
        return None


def extract_kmers_from_endpoints(bedtool, ref_fasta, k=4, max_records=None):
    """Return Counter of k-mers from fragment starts and ends."""
    fasta = pysam.FastaFile(ref_fasta)
    kmers = Counter()
    n = 0
    for interval in bedtool:
        chrom = interval.chrom
        start = int(interval.start)
        end = int(interval.end)
        # start k-mer: from start to start+k
        try:
            seq_start = fasta.fetch(chrom, start, min(start + k, end))
        except Exception:
            seq_start = ""
        # end k-mer: last k bases (end-k to end)
        try:
            seq_end = fasta.fetch(chrom, max(end - k, start), end)
        except Exception:
            seq_end = ""
        for seq in (seq_start, seq_end):
            seq = seq.upper()
            if len(seq) == k and all(c in "ATGC" for c in seq):
                kmers[seq] += 1
        n += 1
        if max_records is not None and n >= max_records:
            break
    fasta.close()
    return kmers


def promoter_genebody_counts(bedtool, gtf_path, promoter_window=(-2000, 500)):
    """Compute counts overlapping promoters and gene bodies from a GTF file.
    Returns counts: promoter_count, gene_body_count
    """
    if not gtf_path or not os.path.exists(gtf_path):
        return None

    # Build promoters BedTool from GTF
    promoters = []
    genes = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = cols
            if feature == "gene":
                genes.append((chrom, int(start) - 1, int(end)))

    if not genes:
        # fallback: parse transcripts
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 9:
                    continue
                chrom, _, feature, start, end, _, strand, _, attrs = cols
                if feature == "transcript":
                    if strand == "+":
                        tss = int(start) - 1
                    else:
                        tss = int(end)
                    pstart = max(0, tss + promoter_window[0])
                    pend = tss + promoter_window[1]
                    promoters.append((chrom, pstart, max(pstart + 1, pend)))
    else:
        for chrom, s, e in genes:
            tss = s
            pstart = max(0, tss + promoter_window[0])
            pend = tss + promoter_window[1]
            promoters.append((chrom, pstart, max(pstart + 1, pend)))

    promoters_bt = BedTool(["\t".join(map(str, p)) for p in promoters])
    genes_bt = BedTool(["\t".join(map(str, g)) for g in genes]) if genes else None

    promoter_count = int(bedtool.intersect(promoters_bt, u=True).count())
    gene_body_count = int(bedtool.intersect(genes_bt, u=True).count()) if genes_bt is not None else None
    return {"promoter_count": promoter_count, "gene_body_count": gene_body_count}


def cpg_counts(bedtool, cpg_bed):
    if not cpg_bed or not os.path.exists(cpg_bed):
        return None
    cpg_bt = BedTool(cpg_bed)
    cpg_overlap = int(bedtool.intersect(cpg_bt, u=True).count())
    return {"cpg_overlap_count": cpg_overlap}


def chromosome_level_counts(bedtool):
    chrom_counts = defaultdict(int)
    for interval in bedtool:
        chrom_counts[interval.chrom] += 1
    return dict(chrom_counts)


def rpm_normalize(df, count_col="total_fragments"):
    # convert count columns to RPM per sample (per million)
    df = df.copy()
    total = df[count_col].replace(0, np.nan)
    for col in df.columns:
        if col.startswith("range_") or col.startswith("prop_") or col.endswith("_count"):
            # convert counts to RPM
            df[col + "_rpm"] = df[col] / (total / 1e6)
    return df


def collect_features_extended(cfg, fragments_dir, sample_sheet_path, out_path, normalize_method=None, kmer_k=4, max_records_kmer=100000):
    sample_df = pd.read_csv(sample_sheet_path)
    features = []

    ref = cfg["reference"].get("genome")
    gtf = cfg["reference"].get("gtf")
    cpg_bed = cfg["reference"].get("cpg_bed") or cfg.get("reference", {}).get("cpg_bed")
    frag_ranges = cfg.get("fragmentation", {}).get("fragment_ranges", {})

    chrom_set = set()

    for _, row in sample_df.iterrows():
        sample_id = row["sample_id"]
        metrics_file = os.path.join(fragments_dir, f"{sample_id}_fragment_metrics.csv")
        bed_file = os.path.join(fragments_dir, f"{sample_id}_fragments.bed")
        metrics = read_fragment_metrics(metrics_file)
        bedtool = read_fragments_bed(bed_file)
        if metrics is None or bedtool is None:
            logging.warning(f"Missing data for {sample_id}, skipping")
            continue

        feat = {"sample_id": sample_id, "group": row.get("group", "unknown")}
        feat["total_fragments"] = int(metrics.get("n_fragments", 0))
        feat["mean_length"] = float(metrics.get("mean_length", np.nan))
        feat["median_length"] = float(metrics.get("median_length", np.nan))
        feat["std_length"] = float(metrics.get("std_length", np.nan))
        feat["mean_gc"] = float(metrics.get("mean_gc", np.nan))
        feat["nucleosome_peak_length"] = metrics.get("nucleosome_peak_length", np.nan)
        feat["nucleosome_peak_count"] = int(metrics.get("nucleosome_peak_count", 0))

        # fragment ranges counts and proportions
        for name, (lo, hi) in frag_ranges.items():
            key = f"range_{name}"
            feat[key] = int(metrics.get(key, 0))
            feat[f"prop_{name}"] = float(feat[key] / feat["total_fragments"]) if feat["total_fragments"] > 0 else np.nan

        # k-mer endpoints
        try:
            kmers = extract_kmers_from_endpoints(bedtool, ref, k=kmer_k, max_records=max_records_kmer)
            # take top N kmers or encode as frequencies of top 20 kmers across samples
            total_k = sum(kmers.values()) if kmers else 0
            for kmer, cnt in kmers.most_common(20):
                feat[f"kmer_{kmer}"] = cnt / total_k if total_k > 0 else 0
        except Exception as e:
            logging.warning(f"Failed to compute k-mers for {sample_id}: {e}")

        # promoter/gene body
        try:
            pg = promoter_genebody_counts(bedtool, gtf)
            if pg:
                feat.update(pg)
        except Exception as e:
            logging.warning(f"Promoter/gene body analysis failed for {sample_id}: {e}")

        # CpG
        try:
            cpg = cpg_counts(bedtool, cpg_bed)
            if cpg:
                feat.update(cpg)
        except Exception as e:
            logging.warning(f"CpG overlap analysis failed for {sample_id}: {e}")

        # chromosome-level counts (normalized later)
        chrom_counts = chromosome_level_counts(bedtool)
        for chrom, cnt in chrom_counts.items():
            feat[f"chr_{chrom}_count"] = cnt
            chrom_set.add(chrom)

        features.append(feat)

    df_feat = pd.DataFrame(features)
    df_feat = df_feat.fillna(0)

    # Ensure consistent kmer columns across samples: fill missing with zero
    # Find all kmer columns
    kmer_cols = sorted([c for c in df_feat.columns if c.startswith("kmer_")])
    # chromosome columns
    chrom_cols = sorted([c for c in df_feat.columns if c.startswith("chr_")])

    # Normalization
    if normalize_method == "zscore":
        num_cols = df_feat.select_dtypes(include=[np.number]).columns.tolist()
        num_cols = [c for c in num_cols if c not in ("total_fragments",)]
        df_feat[num_cols] = df_feat[num_cols].apply(lambda x: stats.zscore(x, nan_policy='omit'))
    elif normalize_method == "quantile":
        num_cols = df_feat.select_dtypes(include=[np.number]).columns.tolist()
        num_cols = [c for c in num_cols if c not in ("total_fragments",)]
        if num_cols:
            arr = df_feat[num_cols].values
            arr_q = quantile_transform(arr, axis=0, copy=True, n_quantiles=min(1000, arr.shape[0]))
            df_feat[num_cols] = arr_q
    elif normalize_method == "rpm":
        df_feat = rpm_normalize(df_feat, count_col="total_fragments")

    ensure_dir(os.path.dirname(out_path))
    df_feat.to_csv(out_path, index=False)
    logging.info(f"Wrote extended feature matrix to {out_path}")
    return df_feat


def rpm_normalize(df, count_col="total_fragments"):
    df = df.copy()
    total = df[count_col].replace(0, np.nan)
    count_cols = [c for c in df.columns if c.endswith("_count") or c.startswith("range_")]
    for col in count_cols:
        df[col + "_rpm"] = df[col] / (total / 1e6)
    return df


def main():
    parser = argparse.ArgumentParser(description="Extended feature extraction from fragment BEDs")
    parser.add_argument("--config", required=True)
    parser.add_argument("--fragments-dir", required=True)
    parser.add_argument("--sample-sheet", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--normalize", choices=["zscore", "quantile", "rpm", "none"], default="zscore")
    parser.add_argument("--kmer-k", type=int, default=4)
    parser.add_argument("--max-kmers", type=int, default=100000)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s: %(message)s")
    cfg = load_config(args.config)
    norm = None if args.normalize == "none" else args.normalize
    collect_features_extended(cfg, args.fragments_dir, args.sample_sheet, args.out, normalize_method=norm, kmer_k=args.kmer_k, max_records_kmer=args.max_kmers)


if __name__ == "__main__":
    main()
