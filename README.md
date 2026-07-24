# cfDNA Fragmentomics Analysis Pipeline for NSCLC

A comprehensive, production-ready workflow for analyzing cell-free DNA (cfDNA) fragmentation patterns from low-pass whole genome sequencing (WGS) data in Non-Small Cell Lung Cancer (NSCLC) and healthy controls.

## Overview

This pipeline performs integrated analysis of cfDNA fragmentomics with the following stages:

1. **Quality Control (QC)** - FastQC, adapter trimming, quality filtering
2. **Alignment** - BWA mapping to reference genome, BAM processing
3. **Fragmentation Analysis** - Fragment length distribution, nucleosome positioning
4. **Feature Extraction** - Comprehensive fragmentomic biomarkers
5. **Comparative Analysis** - Statistical comparison between cancer and control groups

## Key Features

- ✅ Batch processing of multiple samples
- ✅ Comprehensive quality metrics and reporting
- ✅ Fragment-level analysis (length, GC%, end-repair patterns)
- ✅ Nucleosome-scale fragmentation patterns
- ✅ Cancer vs. control statistical comparison
- ✅ Machine learning classification ready
- ✅ Reproducible and containerized analysis

## Quick Start

### Requirements

- Python 3.9+
- Conda/Mamba for environment management
- 100GB+ disk space (for reference genome and results)
- Linux/macOS (or WSL on Windows)

### Installation

```bash
# Clone repository
git clone https://github.com/GeneGee/cfdna-fragmentomics-nsclc.git
cd cfdna-fragmentomics-nsclc

# Create conda environment
conda env create -f environment.yml
conda activate cfdna-frag

# Download reference genome (hg38)
bash scripts/download_reference.sh

# Prepare sample sheet
cp config/sample_sheet_template.csv config/sample_sheet.csv
# Edit config/sample_sheet.csv with your sample information
```

### Running the Pipeline

```bash
# Basic run
python scripts/run_pipeline.py --config config/config.yaml --sample-sheet config/sample_sheet.csv

# With specific modules
python scripts/run_pipeline.py --config config/config.yaml --sample-sheet config/sample_sheet.csv --modules qc mapping fragmentation feature-extraction comparison

# Dry run (show what would be executed)
python scripts/run_pipeline.py --config config/config.yaml --sample-sheet config/sample_sheet.csv --dry-run
```

## Pipeline Modules

### 1. Quality Control (QC)
```bash
python workflow/01_qc.py --sample-sheet config/sample_sheet.csv
```
- FastQC analysis
- Adapter sequence trimming (Cutadapt)
- Low-quality read filtering
- Quality report generation

### 2. Alignment & BAM Processing
```bash
python workflow/02_mapping.py --config config/config.yaml --sample-sheet config/sample_sheet.csv
```
- BWA-MEM alignment to hg38
- BAM sorting and indexing
- Duplicate marking
- Coverage statistics

### 3. Fragmentation Analysis
```bash
python workflow/03_fragmentation_analysis.py --bam-dir data/results/bam --output-dir data/results/fragments
```
- Fragment length distribution
- Nucleosome positioning (120-180bp peaks)
- Monomer/dimer/multimer ratios
- GC content analysis
- End repair patterns

### 4. Feature Extraction
```bash
python workflow/04_feature_extraction.py --fragment-dir data/results/fragments --output-dir data/results/features
```
- 300+ fragmentomic features
- Fragment size ranges (30-50bp, 50-100bp, 100-150bp, etc.)
- Chromosome-level patterns
- Gene body vs. promoter ratios
- Nucleosome positioning scores

### 5. Comparative Analysis
```bash
python workflow/05_comparison_analysis.py --feature-dir data/results/features --sample-sheet config/sample_sheet.csv --output-dir data/results/comparison
```
- Statistical testing (t-test, Mann-Whitney U)
- Multi-group comparison (ANOVA)
- Effect size calculation
- ROC curve analysis
- Feature selection (random forest importance)
- Heatmap visualization

## Input Format

### Sample Sheet (CSV)
```
sample_id,group,fastq_r1,fastq_r2,description
NSCLC_001,cancer,data/raw/NSCLC_001_R1.fastq.gz,data/raw/NSCLC_001_R2.fastq.gz,NSCLC stage III
NSCLC_002,cancer,data/raw/NSCLC_002_R1.fastq.gz,data/raw/NSCLC_002_R2.fastq.gz,NSCLC stage II
Control_001,control,data/raw/Control_001_R1.fastq.gz,data/raw/Control_001_R2.fastq.gz,Healthy volunteer
Control_002,control,data/raw/Control_002_R1.fastq.gz,data/raw/Control_002_R2.fastq.gz,Healthy volunteer
```

## Output Structure

```
data/results/
├── qc/                          # Quality control reports
│   ├── fastqc_raw/
│   ├── fastqc_trimmed/
│   └── qc_summary.csv
├── bam/                         # Aligned BAM files
│   ├── sample_001.bam
│   ├── sample_001.bam.bai
│   └── alignment_stats.txt
├── fragments/                   # Fragment-level data
│   ├── sample_001_fragments.bed
│   ├── sample_001_length_distribution.pdf
│   └── fragment_metrics.csv
├── features/                    # Extracted features
│   ├── fragmentomic_features.csv
│   ├── feature_matrix.h5
│   └── feature_correlation.pdf
└── comparison/                  # Comparative analysis
    ├── statistical_results.csv
    ├── roc_curves.pdf
    ├── feature_importance.pdf
    └── comparison_report.html
```

## Configuration

Edit `config/config.yaml` to customize:
- Reference genome path
- Trimming parameters
- Alignment settings
- Fragment size ranges for analysis
- Statistical test parameters
- Output directories

## Citation

If you use this pipeline, please cite:

```
[Citation to be added after publication]
```

## Support & Troubleshooting

See [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) for common issues.

For questions, please open an issue on GitHub.

## License

MIT License - see LICENSE file for details

## Authors

- GeneGee - Development and maintenance
