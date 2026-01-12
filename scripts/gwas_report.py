#!/usr/bin/env python3
"""
GWAS report pack for pyseer outputs with Wes Anderson-style palettes.

Inputs:
- combined_gwas.all.tsv from combine_pyseer_results.py
- optional: genes annotated hits from annotate_gene_hits_pirate.py (for PIRATE-family candidate list)

Outputs (under outdir):
- tables/
  - top_hits_<type>.tsv
  - sig_hits_q<q>_<type>.tsv
  - candidate_families_q<q>.tsv   (if annotated genes provided)
- plots/
  - qq_<type>.svg/.png
  - manhattan_<type>.svg/.png
  - volcano_<type>.svg/.png

Notes:
- For SNP Manhattan: expects variant like "chrom_pos_ref_alt" (e.g., 1_1383_G_A).
- For genes/unitigs Manhattan: uses rank on x-axis.
"""

from __future__ import annotations
import argparse
import os
import re
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---- Wes Anderson-inspired palette (approx) ----
# A set that prints nicely and stays distinct.
WES = {
    "gene":   "#D55E00",  # warm red/orange
    "snp":    "#0072B2",  # deep blue
    "unitig": "#009E73",  # green
    "all":    "#4D4D4D",  # neutral
}
BG = "#FFFFFF"


def ensure_outdirs(outdir: str):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "tables"), exist_ok=True)


def read_combined(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    # numeric coercions
    for c in ["af", "beta", "p_lrt", "p_filter", "p_bonf", "q_fdr"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def save_fig(fig, outbase: str):
    fig.savefig(outbase + ".svg", bbox_inches="tight")
    fig.savefig(outbase + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def qq_plot(df: pd.DataFrame, pcol: str, title: str, color: str):
    p = df[pcol].dropna().astype(float)
    p = p[(p > 0) & (p <= 1)]
    if len(p) < 5:
        return None

    p_sorted = np.sort(p.values)
    n = len(p_sorted)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(p_sorted)

    fig = plt.figure(figsize=(5.2, 5.2))
    ax = fig.add_subplot(111)
    ax.scatter(exp, obs, s=10, alpha=0.8, edgecolor="none", c=color)
    maxv = max(exp.max(), obs.max())
    ax.plot([0, maxv], [0, maxv], linewidth=1)
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title(title)
    ax.grid(True, linewidth=0.4, alpha=0.3)
    return fig


def parse_snp_variant(v: str):
    # common pyseer variant: chrom_pos_ref_alt
    m = re.match(r"^([^_]+)_(\d+)_([ACGTN]+)_([ACGTN]+)$", str(v))
    if not m:
        return None
    chrom = m.group(1)
    pos = int(m.group(2))
    return chrom, pos


def manhattan_like(df: pd.DataFrame, pcol: str, gwas_type: str, title: str, color: str):
    d = df.copy()
    d = d[[c for c in d.columns if c in ["variant", pcol]]].dropna()
    if d.shape[0] < 5:
        return None

    # x-axis:
    if gwas_type == "snp":
        parsed = d["variant"].apply(parse_snp_variant)
        ok = parsed.notna()
        d = d.loc[ok].copy()
        if d.shape[0] < 5:
            return None
        d["chrom"] = parsed[ok].apply(lambda x: x[0])
        d["pos"] = parsed[ok].apply(lambda x: x[1])
        # order by chrom then pos; chrom might be numeric or string
        # make an ordering key that sorts numeric chroms properly
        def chrom_key(x):
            try:
                return (0, int(x))
            except Exception:
                return (1, str(x))
        d["chrom_key"] = d["chrom"].map(chrom_key)
        d = d.sort_values(["chrom_key", "pos"])
        # create cumulative x with simple offsets per chrom
        x = np.arange(len(d))
        xlabel = "SNPs (ordered by position)"
    else:
        # genes/unitigs: just rank order
        d = d.sort_values(pcol)
        x = np.arange(len(d))
        xlabel = f"{gwas_type} variants (ranked)"

    y = -np.log10(d[pcol].astype(float).clip(lower=1e-300))

    fig = plt.figure(figsize=(9.0, 3.8))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, s=8, alpha=0.85, edgecolor="none", c=color)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    ax.grid(True, linewidth=0.4, alpha=0.25)
    return fig


def volcano(df: pd.DataFrame, pcol: str, title: str, color: str):
    d = df.copy()
    d = d.dropna(subset=["beta", pcol])
    if d.shape[0] < 5:
        return None
    x = d["beta"].astype(float)
    y = -np.log10(d[pcol].astype(float).clip(lower=1e-300))

    fig = plt.figure(figsize=(6.6, 4.2))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, s=10, alpha=0.75, edgecolor="none", c=color)
    ax.set_xlabel("Effect size (beta)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    ax.grid(True, linewidth=0.4, alpha=0.25)
    return fig


def write_top_tables(df: pd.DataFrame, outdir: str, gwas_type: str, pcol: str, q: float, top_n: int):
    d = df[df["gwas_type"] == gwas_type].copy()
    d = d.dropna(subset=[pcol])
    d = d.sort_values(pcol, ascending=True)

    top = d.head(top_n)
    top.to_csv(os.path.join(outdir, "tables", f"top_hits_{gwas_type}.tsv"), sep="\t", index=False)

    sig = d[d["q_fdr"] < q]
    sig.to_csv(os.path.join(outdir, "tables", f"sig_hits_q{q}_{gwas_type}.tsv"), sep="\t", index=False)

    return len(d), len(sig)


def candidate_families(annotated_gene_hits: str, outdir: str, pcol: str, q: float):
    """
    Collapse gene allele hits to PIRATE family.
    Requires output from annotate_gene_hits_pirate.py (contains pirate_family + q_fdr).
    """
    df = pd.read_csv(annotated_gene_hits, sep="\t", dtype=str)
    for c in [pcol, "q_fdr", "beta", "af"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    if "pirate_family" not in df.columns:
        raise SystemExit("annotated gene hits file lacks pirate_family column")

    df["is_sig"] = df["q_fdr"] < q
    fam = (
        df.groupby("pirate_family", dropna=False)
          .agg(
              n_alleles=("variant", "count"),
              n_sig=("is_sig", "sum"),
              min_p=(pcol, "min"),
              min_q=("q_fdr", "min"),
              top_variant=("variant", lambda x: x.iloc[df.loc[x.index, pcol].astype(float).argmin()]),
          )
          .reset_index()
          .sort_values(["min_q", "min_p"], ascending=[True, True])
    )
    out = os.path.join(outdir, "tables", f"candidate_families_q{q}.tsv")
    fam.to_csv(out, sep="\t", index=False)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", required=True, help="combined_gwas.all.tsv")
    ap.add_argument("--outdir", default="pyseer/report_pack", help="output directory")
    ap.add_argument("--pcol", default="p_lrt", help="p-value column (default p_lrt)")
    ap.add_argument("--sig-q", type=float, default=0.05, help="FDR threshold (default 0.05)")
    ap.add_argument("--top-n", type=int, default=50, help="top hits per type (default 50)")
    ap.add_argument("--annotated-genes", help="optional: genes.annotated_hits.tsv to produce PIRATE family candidates")
    args = ap.parse_args()

    ensure_outdirs(args.outdir)
    df = read_combined(args.combined)

    if args.pcol not in df.columns:
        raise SystemExit(f"pcol {args.pcol} not found in {args.combined}. Columns: {df.columns.tolist()[:25]}...")

    # per-type outputs
    types = [t for t in ["gene", "snp", "unitig"] if t in set(df["gwas_type"].dropna())]
    summary = []

    for t in types:
        n_tests, n_sig = write_top_tables(df, args.outdir, t, args.pcol, args.sig_q, args.top_n)
        summary.append((t, n_tests, n_sig))

        color = WES.get(t, WES["all"])

        qq = qq_plot(df[df["gwas_type"] == t], args.pcol, f"QQ plot: {t}", color)
        if qq:
            save_fig(qq, os.path.join(args.outdir, "plots", f"qq_{t}"))

        mh = manhattan_like(df[df["gwas_type"] == t], args.pcol, t, f"Manhattan-like: {t}", color)
        if mh:
            save_fig(mh, os.path.join(args.outdir, "plots", f"manhattan_{t}"))

        volc = volcano(df[df["gwas_type"] == t], args.pcol, f"Volcano: {t}", color)
        if volc:
            save_fig(volc, os.path.join(args.outdir, "plots", f"volcano_{t}"))

    # optional: PIRATE family candidate list
    if args.annotated_genes:
        out = candidate_families(args.annotated_genes, args.outdir, args.pcol, args.sig_q)
        print(f"Wrote PIRATE family candidates: {out}")

    # write a tiny text summary
    with open(os.path.join(args.outdir, "tables", "SUMMARY.txt"), "w") as fh:
        fh.write("GWAS report pack summary\n")
        fh.write(f"Input: {args.combined}\n")
        fh.write(f"FDR threshold: {args.sig_q}\n")
        fh.write("\nPer type:\n")
        for t, n_tests, n_sig in summary:
            fh.write(f"  - {t}: tests={n_tests}  sig(q<{args.sig_q})={n_sig}\n")

    print(f"Done. Outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
