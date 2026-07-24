#!/usr/bin/env python3
"""
workflow/05_comparison_analysis.py

Performs group-wise statistical comparison and baseline classification.
Inputs:
 - feature matrix CSV (samples x features)
 - sample_sheet CSV (to get group labels)
Outputs (in out_dir):
 - statistical_results.csv
 - feature_importance.csv
 - roc_curves.png
 - comparison_report.html (basic static report)

Default behavior:
 - two-group t-test per feature, BH FDR correction
 - compute Cohen's d
 - random forest baseline classification with StratifiedCV and AUC
 - output plots and CSVs

Usage:
  python workflow/05_comparison_analysis.py --config config/config.yaml --features data/results/features/fragmentomic_features.csv --sample-sheet config/sample_sheet.csv --out-dir data/results/comparison
"""

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import LabelEncoder
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import yaml


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def cohen_d(x, y):
    # pooled sd
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    if dof <= 0:
        return np.nan
    pooled_sd = np.sqrt(((nx - 1) * np.nanvar(x, ddof=1) + (ny - 1) * np.nanvar(y, ddof=1)) / dof)
    if pooled_sd == 0:
        return np.nan
    return (np.nanmean(x) - np.nanmean(y)) / pooled_sd


def statistical_tests(df_feat, group_col="group", method="ttest"):
    groups = df_feat[group_col].unique()
    if len(groups) != 2 and method in ("ttest", "mannwhitney"):
        raise ValueError("Statistical tests currently support two groups for t-test or Mann-Whitney; use ANOVA for >2 groups")

    grp1, grp2 = groups[0], groups[1]
    a = df_feat[df_feat[group_col] == grp1]
    b = df_feat[df_feat[group_col] == grp2]

    numeric_cols = df_feat.select_dtypes(include=[np.number]).columns.tolist()
    # exclude total_fragments maybe
    if "total_fragments" in numeric_cols:
        numeric_cols.remove("total_fragments")

    results = []
    for col in numeric_cols:
        x = a[col].values
        y = b[col].values
        # handle NaNs by dropping
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]
        if method == "ttest":
            try:
                t_stat, p = stats.ttest_ind(x, y, equal_var=False, nan_policy='omit')
            except Exception:
                p = np.nan
        elif method == "mannwhitney":
            try:
                p = stats.mannwhitneyu(x, y, alternative='two-sided').pvalue
            except Exception:
                p = np.nan
        else:
            p = np.nan

        mean1 = np.nanmean(x) if x.size else np.nan
        mean2 = np.nanmean(y) if y.size else np.nan
        fold = mean1 / mean2 if (mean2 not in (0, np.nan) and mean1 not in (np.nan, 0)) else np.nan
        d = cohen_d(x, y)
        results.append({"feature": col, "group1": grp1, "group2": grp2, "mean1": mean1, "mean2": mean2, "fold_change": fold, "p_value": p, "cohens_d": d})

    df_res = pd.DataFrame(results)
    # FDR correction
    pvals = df_res["p_value"].fillna(1).values
    rejected, qvals, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    df_res["q_value"] = qvals
    df_res["significant"] = rejected
    df_res = df_res.sort_values("q_value")
    return df_res


def random_forest_cv(X, y, n_splits=5, n_jobs=1, n_estimators=100, random_state=0):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    aucs = []
    importances = np.zeros(X.shape[1])
    mean_fpr = np.linspace(0, 1, 100)
    tprs = []

    for train_idx, test_idx in skf.split(X, y):
        Xtr, Xte = X[train_idx], X[test_idx]
        ytr, yte = y[train_idx], y[test_idx]
        clf = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, random_state=random_state)
        clf.fit(Xtr, ytr)
        probs = clf.predict_proba(Xte)[:, 1]
        try:
            auc = roc_auc_score(yte, probs)
        except Exception:
            auc = np.nan
        aucs.append(auc)
        importances += clf.feature_importances_
        fpr, tpr, _ = roc_curve(yte, probs)
        # interpolate tpr at mean_fpr
        tpr_interp = np.interp(mean_fpr, fpr, tpr)
        tprs.append(tpr_interp)

    importances /= n_splits
    mean_auc = np.nanmean(aucs)
    std_auc = np.nanstd(aucs)
    mean_tpr = np.nanmean(tprs, axis=0) if tprs else None
    return {"aucs": aucs, "mean_auc": mean_auc, "std_auc": std_auc, "importances": importances, "mean_fpr": mean_fpr, "mean_tpr": mean_tpr}


def make_plots(df_res, feature_names, rf_res, out_dir):
    ensure_dir(out_dir)
    # Volcano plot (log2 fold vs -log10 p)
    plt.figure(figsize=(6, 5))
    fc = df_res["fold_change"].replace(0, np.nan)
    # avoid divide by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        logfc = np.log2(fc)
    neglogp = -np.log10(df_res["p_value"].replace(0, np.nan))
    plt.scatter(logfc, neglogp, s=10, alpha=0.6)
    plt.xlabel("log2 fold change")
    plt.ylabel("-log10(p)")
    plt.title("Volcano plot")
    plt.axhline(-np.log10(0.05), color='red', linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "volcano.png"))
    plt.close()

    # ROC curve
    if rf_res and rf_res.get("mean_tpr") is not None:
        plt.figure(figsize=(6, 6))
        plt.plot(rf_res["mean_fpr"], rf_res["mean_tpr"], color='b', label=f"Mean ROC (AUC={rf_res['mean_auc']:.2f} +/- {rf_res['std_auc']:.2f})")
        plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC Curve")
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "roc_curve.png"))
        plt.close()

    # Feature importance barplot (top 20)
    if rf_res is not None:
        importances = rf_res["importances"]
        idx = np.argsort(importances)[::-1][:20]
        names = [feature_names[i] for i in idx]
        vals = importances[idx]
        plt.figure(figsize=(8, 6))
        plt.barh(range(len(names))[::-1], vals)
        plt.yticks(range(len(names))[::-1], names)
        plt.xlabel("Feature importance")
        plt.title("Top features (RF)")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "feature_importance.png"))
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="Comparison analysis: stats + RF baseline")
    parser.add_argument("--config", required=True)
    parser.add_argument("--features", required=True, help="Feature matrix CSV")
    parser.add_argument("--sample-sheet", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--test", choices=["ttest", "mannwhitney"], default="ttest")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s: %(message)s")
    cfg = load_config(args.config)
    ensure_dir(args.out_dir)

    df_feat = pd.read_csv(args.features)
    sample_df = pd.read_csv(args.sample_sheet)

    # merge to ensure same order and groups
    df = pd.merge(sample_df[["sample_id", "group"]], df_feat, on="sample_id")

    df_res = statistical_tests(df, group_col="group", method=args.test)
    out_stats = os.path.join(args.out_dir, "statistical_results.csv")
    df_res.to_csv(out_stats, index=False)

    # Prepare data for RF (binary classification required)
    groups = df["group"].unique()
    if len(groups) != 2:
        logging.warning("Random forest classification requires exactly 2 groups; skipping RF")
        rf_res = None
    else:
        X = df.select_dtypes(include=[np.number]).values
        # remove total_fragments column if present
        cols = df.select_dtypes(include=[np.number]).columns.tolist()
        if "total_fragments" in cols:
            tidx = cols.index("total_fragments")
            X = np.delete(X, tidx, axis=1)
            cols.pop(tidx)
        le = LabelEncoder()
        y = le.fit_transform(df["group"].values)
        rf_res = random_forest_cv(X, y, n_splits=args.cv_folds, n_jobs=args.n_jobs, n_estimators=100, random_state=args.random_state)

        # feature importance output
        feat_imp = pd.DataFrame({"feature": cols, "importance": rf_res["importances"]})
        feat_imp = feat_imp.sort_values("importance", ascending=False)
        feat_imp.to_csv(os.path.join(args.out_dir, "feature_importance.csv"), index=False)

    # plots
    make_plots(df_res, cols if 'cols' in locals() else [], rf_res, args.out_dir)

    # basic HTML report (static)
    report_path = os.path.join(args.out_dir, "comparison_report.html")
    with open(report_path, "w") as fh:
        fh.write("<html><head><title>Comparison report</title></head><body>\n")
        fh.write("<h1>Comparison results</h1>\n")
        fh.write(f"<p>Statistical results: <a href=\"statistical_results.csv\">statistical_results.csv</a></p>\n")
        fh.write(f"<p>Feature importance: <a href=\"feature_importance.csv\">feature_importance.csv</a></p>\n")
        fh.write(f"<p>Volcano plot: <img src=\"volcano.png\" width=600></p>\n")
        fh.write(f"<p>ROC curve: <img src=\"roc_curve.png\" width=600></p>\n")
        fh.write(f"<p>Feature importance: <img src=\"feature_importance.png\" width=600></p>\n")
        fh.write("</body></html>\n")

    logging.info(f"Comparison outputs written to {args.out_dir}")


if __name__ == "__main__":
    main()
