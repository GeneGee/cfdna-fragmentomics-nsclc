#!/usr/bin/env python3
"""
workflow/04_feature_extraction.py

Builds a per-sample feature matrix from fragment BED files produced by workflow/03_fragmentation_analysis.py
Features (first version):
 - total_fragments
 - mean_length, median_length, std_length
 - mean_gc
 - proportion of fragments in configured fragment_ranges (e.g., short, medium, long)
 - nucleosome_peak_length, nucleosome_peak_count

Usage:
  python workflow/04_feature_extraction.py --config config/config.yaml --fragments-dir data/results/fragments --sample-sheet config/sample_sheet.csv --out data/results/features/fragmentomic_features.csv

Produces per-sample CSV and aggregated feature matrix.
"""

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def read_fragment_metrics(metrics_path):
    # metrics CSV written by 03_fragmentation_analysis.py
    if not os.path.exists(metrics_path):
        return None
    try:
        df = pd.read_csv(metrics_path)
        return df.iloc[0].to_dict()
    except Exception:
        return None


def collect_features(cfg, fragments_dir, sample_sheet_path, out_path, normalize=False):
    sample_df = pd.read_csv(sample_sheet_path)
    features = []
    frag_ranges = cfg.get("fragmentation", {}).get("fragment_ranges", {})

    for _, row in sample_df.iterrows():
        sample_id = row["sample_id"]
        metrics_file = os.path.join(fragments_dir, f"{sample_id}_fragment_metrics.csv")
        metrics = read_fragment_metrics(metrics_file)
        if metrics is None:
            logging.warning(f"No metrics for sample {sample_id} (expected at {metrics_file}), skipping")
            continue

        feat = {"sample_id": sample_id, "group": row.get("group", "unknown")}
        # core metrics
        feat["total_fragments"] = int(metrics.get("n_fragments", 0))
        feat["mean_length"] = float(metrics.get("mean_length", np.nan))
        feat["median_length"] = float(metrics.get("median_length", np.nan))
        feat["std_length"] = float(metrics.get("std_length", np.nan))
        feat["mean_gc"] = float(metrics.get("mean_gc", np.nan))
        feat["nucleosome_peak_length"] = metrics.get("nucleosome_peak_length", np.nan)
        feat["nucleosome_peak_count"] = int(metrics.get("nucleosome_peak_count", 0))

        # fragment ranges
        for name in frag_ranges.keys():
            key = f"range_{name}"
            feat[key] = int(metrics.get(key, 0))
            # also proportion
            feat[f"prop_{name}"] = float(feat[key] / feat["total_fragments"]) if feat["total_fragments"] > 0 else np.nan

        features.append(feat)

    df_feat = pd.DataFrame(features)

    # Optional normalization (z-score) for numeric features excluding sample_id and group
    if normalize:
        num_cols = df_feat.select_dtypes(include=[np.number]).columns.tolist()
        num_cols = [c for c in num_cols if c not in ("total_fragments",)]
        df_feat[num_cols] = df_feat[num_cols].apply(lambda x: stats.zscore(x, nan_policy='omit'))

    ensure_dir(os.path.dirname(out_path))
    df_feat.to_csv(out_path, index=False)
    logging.info(f"Wrote feature matrix to {out_path}")
    return df_feat


def main():
    parser = argparse.ArgumentParser(description="Feature extraction from fragment metrics")
    parser.add_argument("--config", required=True)
    parser.add_argument("--fragments-dir", required=True)
    parser.add_argument("--sample-sheet", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s: %(message)s")
    cfg = load_config(args.config)
    collect_features(cfg, args.fragments_dir, args.sample_sheet, args.out, normalize=args.normalize)


if __name__ == "__main__":
    main()
