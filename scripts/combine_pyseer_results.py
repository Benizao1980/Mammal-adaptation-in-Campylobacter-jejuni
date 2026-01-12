#!/usr/bin/env python3
"""
Combine pyseer GWAS outputs (genes / SNPs / unitigs) into one harmonised table.

- Standardises column names: variant, af, beta, p_lrt, p_filter
- Adds: gwas_type, n_tests, p_bonf, q_fdr (BH within each gwas_type)
- Writes:
  - combined_gwas.all.tsv
  - combined_gwas.sig.tsv (default: q_fdr < 0.05)
  - combined_gwas.top.tsv (top N per type by p_lrt)
"""

from __future__ import annotations
import argparse
import math
import os
import sys
from typing import Optional, Dict, List

import pandas as pd


def bh_fdr(pvals: pd.Series) -> pd.Series:
    """Benjamini-Hochberg FDR; returns q-values aligned to original order."""
    p = pvals.astype(float)
    n = p.notna().sum()
    if n == 0:
        return pd.Series([float("nan")] * len(p), index=p.index)

    # rank only non-NA
    pv = p[p.notna()].sort_values()
    m = len(pv)
    ranks = pd.Series(range(1, m + 1), index=pv.index, dtype=float)
    q = pv * m / ranks
    # enforce monotonicity
    q = q[::-1].cummin()[::-1]
    q = q.clip(upper=1.0)

    out = pd.Series([float("nan")] * len(p), index=p.index)
    out.loc[q.index] = q
    return out


def load_pyseer_tsv(path: str, gwas_type: str) -> pd.DataFrame:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"Missing/empty file: {path}")

    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.shape[0] == 0:
        raise ValueError(f"No rows found in: {path}")

    # Standardise column names we care about
    colmap = {
        "lrt-pvalue": "p_lrt",
        "pvalue": "p_lrt",          # fallback
        "filter-pvalue": "p_filter",
        "variant": "variant",
        "af": "af",
        "beta": "beta",
    }
    # create a mapping only for columns that exist
    rename = {c: colmap[c] for c in df.columns if c in colmap}
    df = df.rename(columns=rename)

    # Ensure required columns exist
    if "variant" not in df.columns:
        # some older outputs might have first col unnamed
        first = df.columns[0]
        df = df.rename(columns={first: "variant"})
    if "p_lrt" not in df.columns:
        raise ValueError(f"{path}: no lrt-pvalue or pvalue column found")

    # coerce numerics
    for c in ["af", "beta", "p_lrt", "p_filter"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # add metadata
    df["gwas_type"] = gwas_type
    df["source_file"] = os.path.basename(path)

    # keep a consistent minimal set + carry through PCs/notes etc.
    # We'll keep all columns but ensure the core ones are front-loaded later.
    return df


def add_multiple_testing(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["n_tests"] = out.groupby("gwas_type")["variant"].transform("count")

    # Bonferroni within each gwas_type
    out["p_bonf"] = out["p_lrt"] * out["n_tests"].astype(float)
    out["p_bonf"] = out["p_bonf"].clip(upper=1.0)

    # BH-FDR within each gwas_type
    out["q_fdr"] = out.groupby("gwas_type", group_keys=False)["p_lrt"].apply(bh_fdr)

    # effect direction convenience
    out["direction"] = out["beta"].apply(
        lambda x: "positive" if pd.notna(x) and x > 0 else ("negative" if pd.notna(x) and x < 0 else "NA")
    )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", help="pyseer genes GWAS TSV (e.g. gwas_genes_*.tsv)")
    ap.add_argument("--snps", help="pyseer SNP GWAS TSV (e.g. gwas_snps_*.tsv)")
    ap.add_argument("--unitigs", help="pyseer unitigs GWAS TSV (optional)")
    ap.add_argument("--out-prefix", default="pyseer/combined_gwas",
                    help="Output prefix (default: pyseer/combined_gwas)")
    ap.add_argument("--sig-q", type=float, default=0.05, help="FDR threshold for sig table (default 0.05)")
    ap.add_argument("--top-n", type=int, default=50, help="Top N per type for top table (default 50)")
    args = ap.parse_args()

    frames: List[pd.DataFrame] = []
    if args.genes:
        frames.append(load_pyseer_tsv(args.genes, "gene"))
    if args.snps:
        frames.append(load_pyseer_tsv(args.snps, "snp"))
    if args.unitigs:
        frames.append(load_pyseer_tsv(args.unitigs, "unitig"))

    if not frames:
        raise SystemExit("Provide at least one of --genes/--snps/--unitigs")

    df = pd.concat(frames, axis=0, ignore_index=True)

    # Add multiple testing stats
    df = add_multiple_testing(df)

    # Front-load key columns
    key_cols = [
        "gwas_type", "variant", "af", "beta", "direction",
        "p_filter", "p_lrt", "p_bonf", "q_fdr", "n_tests", "source_file"
    ]
    cols = key_cols + [c for c in df.columns if c not in key_cols]
    df = df[cols]

    all_out = f"{args.out_prefix}.all.tsv"
    sig_out = f"{args.out_prefix}.sig_q{args.sig_q}.tsv"
    top_out = f"{args.out_prefix}.top{args.top_n}.tsv"

    df.to_csv(all_out, sep="\t", index=False)

    sig = df[df["q_fdr"] < args.sig_q].copy()
    sig = sig.sort_values(["gwas_type", "q_fdr", "p_lrt"], ascending=[True, True, True])
    sig.to_csv(sig_out, sep="\t", index=False)

    top = (
        df.sort_values(["gwas_type", "p_lrt"], ascending=[True, True])
          .groupby("gwas_type", as_index=False)
          .head(args.top_n)
    )
    top.to_csv(top_out, sep="\t", index=False)

    print(f"Wrote: {all_out}")
    print(f"Wrote: {sig_out}  (rows={len(sig)})")
    print(f"Wrote: {top_out}")


if __name__ == "__main__":
    main()
