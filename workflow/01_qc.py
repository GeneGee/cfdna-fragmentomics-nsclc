#!/usr/bin/env python3
"""
workflow/01_qc.py

Quality control module for cfDNA fragmentomics pipeline.
Performs: FastQC (raw), adapter trimming with cutadapt, FastQC (trimmed), and writes a qc_summary.csv.

Supports a --dry-run mode which only prints planned commands (useful for CI/tests).
"""

import argparse
import csv
import logging
import os
import subprocess
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import yaml


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def plan_commands(sample_row, cfg):
    sample_id = sample_row["sample_id"]
    r1 = sample_row["fastq_r1"]
    r2 = sample_row.get("fastq_r2", "")

    results_dir = cfg["paths"]["results_dir"]
    qc_dir = os.path.join(results_dir, "qc")
    fastqc_raw = os.path.join(qc_dir, "fastqc_raw")
    fastqc_trimmed = os.path.join(qc_dir, "fastqc_trimmed")
    trimmed_dir = os.path.join(results_dir, "trimmed_fastq")

    ensure_dir(fastqc_raw)
    ensure_dir(fastqc_trimmed)
    ensure_dir(trimmed_dir)

    # FastQC commands (raw)
    fastqc_cmd_raw = [
        "fastqc",
        "-o",
        fastqc_raw,
        "-t",
        str(cfg["qc"]["fastqc_threads"]),
        r1,
    ]
    if r2:
        fastqc_cmd_raw.append(r2)

    # Cutadapt trim -> trimmed files
    trimmed_r1 = os.path.join(trimmed_dir, f"{sample_id}_R1.trimmed.fastq.gz")
    trimmed_r2 = os.path.join(trimmed_dir, f"{sample_id}_R2.trimmed.fastq.gz") if r2 else ""

    adapter_r1 = cfg["qc"].get("adapter_r1")
    adapter_r2 = cfg["qc"].get("adapter_r2")
    min_length = cfg["qc"].get("min_length", 30)
    quality_threshold = cfg["qc"].get("quality_threshold", 20)

    if r2:
        cutadapt_cmd = [
            "cutadapt",
            "-a",
            adapter_r1,
            "-A",
            adapter_r2,
            "-q",
            str(quality_threshold),
            "-m",
            str(min_length),
            "-o",
            trimmed_r1,
            "-p",
            trimmed_r2,
            r1,
            r2,
        ]
    else:
        cutadapt_cmd = [
            "cutadapt",
            "-a",
            adapter_r1,
            "-q",
            str(quality_threshold),
            "-m",
            str(min_length),
            "-o",
            trimmed_r1,
            r1,
        ]

    # FastQC on trimmed
    fastqc_cmd_trim = [
        "fastqc",
        "-o",
        fastqc_trimmed,
        "-t",
        str(max(1, cfg["qc"].get("fastqc_threads", 1))),
        trimmed_r1,
    ]
    if r2:
        fastqc_cmd_trim.append(trimmed_r2)

    planned = {
        "sample_id": sample_id,
        "r1": r1,
        "r2": r2,
        "trimmed_r1": trimmed_r1,
        "trimmed_r2": trimmed_r2,
        "fastqc_raw_dir": fastqc_raw,
        "fastqc_trimmed_dir": fastqc_trimmed,
        "fastqc_cmd_raw": fastqc_cmd_raw,
        "cutadapt_cmd": cutadapt_cmd,
        "fastqc_cmd_trim": fastqc_cmd_trim,
    }
    return planned


def run_command(cmd, dry_run=False):
    cmd_str = " ".join(map(str, cmd))
    logging.debug(f"Running: {cmd_str}")
    if dry_run:
        return 0, cmd_str
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.debug(res.stdout)
        return res.returncode, res.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd_str}")
        logging.error(e.stderr)
        return e.returncode, e.stderr


def process_sample(args):
    sample_row, cfg, dry_run = args
    planned = plan_commands(sample_row, cfg)
    sample_id = planned["sample_id"]
    # Execute fastqc raw
    rc1, out1 = run_command(planned["fastqc_cmd_raw"], dry_run=dry_run)
    # Execute cutadapt
    rc2, out2 = run_command(planned["cutadapt_cmd"], dry_run=dry_run)
    # Execute fastqc trimmed
    rc3, out3 = run_command(planned["fastqc_cmd_trim"], dry_run=dry_run)

    summary = {
        "sample_id": sample_id,
        "r1": planned["r1"],
        "r2": planned["r2"],
        "trimmed_r1": planned["trimmed_r1"],
        "trimmed_r2": planned["trimmed_r2"],
        "fastqc_raw_dir": planned["fastqc_raw_dir"],
        "fastqc_trimmed_dir": planned["fastqc_trimmed_dir"],
        "rc_fastqc_raw": rc1,
        "rc_cutadapt": rc2,
        "rc_fastqc_trim": rc3,
    }
    return summary


def main():
    parser = argparse.ArgumentParser(description="QC module: FastQC + Cutadapt + FastQC(trimmed)")
    parser.add_argument("--config", required=True, help="Path to config YAML")
    parser.add_argument("--sample-sheet", required=True, help="CSV with sample_id,fastq_r1,fastq_r2,...")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s: %(message)s")

    cfg = load_config(args.config)
    df = pd.read_csv(args.sample_sheet)

    results_dir = cfg["paths"]["results_dir"]
    qc_dir = os.path.join(results_dir, "qc")
    ensure_dir(qc_dir)

    tasks = [(row._asdict() if hasattr(row, "_asdict") else row.to_dict(), cfg, args.dry_run) for _, row in df.iterrows()]

    # multiprocessing pool
    summaries = []
    if args.workers > 1:
        with Pool(processes=args.workers) as pool:
            for s in pool.imap_unordered(process_sample, tasks):
                summaries.append(s)
    else:
        for t in tasks:
            summaries.append(process_sample(t))

    # Write summary CSV
    summary_path = os.path.join(qc_dir, "qc_summary.csv")
    keys = [
        "sample_id",
        "r1",
        "r2",
        "trimmed_r1",
        "trimmed_r2",
        "fastqc_raw_dir",
        "fastqc_trimmed_dir",
        "rc_fastqc_raw",
        "rc_cutadapt",
        "rc_fastqc_trim",
    ]
    with open(summary_path, "w", newline="") as outcsv:
        writer = csv.DictWriter(outcsv, fieldnames=keys)
        writer.writeheader()
        for s in summaries:
            writer.writerow({k: s.get(k, "") for k in keys})

    logging.info(f"QC summary written to {summary_path}")


if __name__ == "__main__":
    main()
