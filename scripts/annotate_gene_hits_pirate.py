#!/usr/bin/env python3
"""
Annotate pyseer gene GWAS hits using PIRATE.unique_alleles.tsv

Outputs:
- annotated_hits.tsv: allele-level hits with PIRATE family mapping
- family_summary.tsv: family-level summary (min p, n_sig_alleles, etc.)

Assumptions:
- GWAS 'variant' column matches PIRATE 'allele_name' (often column 0, header 'allele_name')
"""

from __future__ import annotations
import argparse
import os
import sys
import pandas as pd


def detect_allele_and_family_cols(pir: pd.DataFrame) -> tuple[str, str]:
    # allele column
    allele_candidates = ["allele_name", "allele", "variant", "Allele", "allele_id"]
    allele_col = next((c for c in allele_candidates if c in pir.columns), None)
    if allele_col is None:
        allele_col = pir.columns[0]  # best guess

    # family/gene column
    family_candidates = ["gene_name", "gene", "family", "cluster", "gene_family", "g_name"]
    family_col = next((c for c in family_candidates if c in pir.columns), None)

    # if not found, try common PIRATE layout: allele_name, gene_name in col1
    if family_col is None and len(pir.columns) > 1:
        family_col = pir.columns[1]

    return allele_col, family_col


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", required=True, help="pyseer gene GWAS TSV (e.g. gwas_genes_*.tsv)")
    ap.add_argument("--pirate-alleles", required=True, help="PIRATE.unique_alleles.tsv")
    ap.add_argument("--out-prefix", default="pyseer/gene_hits", help="Output prefix")
    ap.add_argument("--pcol", default="lrt-pvalue", help="P-value col in GWAS (default lrt-pvalue)")
    ap.add_argument("--sig-q", type=float, default=0.05, help="FDR threshold to count significant alleles")
    args = ap.parse_args()

    gwas = pd.read_csv(args.gwas, sep="\t", dtype=str)
    pir = pd.read_csv(args.pirate_alleles, sep="\t", dtype=str)

    if "variant" not in gwas.columns:
        gwas = gwas.rename(columns={gwas.columns[0]: "variant"})
    # coerce numeric columns
    for c in ["af", "beta", "lrt-pvalue", "pvalue", "filter-pvalue"]:
        if c in gwas.columns:
            gwas[c] = pd.to_numeric(gwas[c], errors="coerce")

    pcol = args.pcol if args.pcol in gwas.columns else ("lrt-pvalue" if "lrt-pvalue" in gwas.columns else "pvalue")
    if pcol not in gwas.columns:
        raise SystemExit(f"No p-value column found. Columns: {gwas.columns.tolist()}")

    # compute BH-FDR over gene GWAS
    gwas = gwas.copy()
    gwas["rank"] = gwas[pcol].rank(method="min", ascending=True)
    m = gwas[pcol].notna().sum()
    gwas["q_fdr"] = (gwas[pcol] * m / gwas["rank"]).clip(upper=1.0)
    gwas = gwas.sort_values(pcol)
    gwas["q_fdr"] = gwas["q_fdr"][::-1].cummin()[::-1]  # enforce monotonicity
    gwas = gwas.drop(columns=["rank"])

    allele_col, family_col = detect_allele_and_family_cols(pir)

    # shrink PIRATE table to useful cols (keep everything if you prefer)
    keep_cols = [allele_col, family_col]
    for extra in ["threshold", "dosage", "no_loci", "cluster_order"]:
        if extra in pir.columns:
            keep_cols.append(extra)
    pir_small = pir[keep_cols].drop_duplicates()

    ann = gwas.merge(pir_small, how="left", left_on="variant", right_on=allele_col)

    # rename columns to stable names
    ann = ann.rename(columns={allele_col: "pirate_allele_name", family_col: "pirate_family"})

    # allele-level annotated output
    ann_out = f"{args.out_prefix}.annotated_hits.tsv"
    ann.to_csv(ann_out, sep="\t", index=False)

    # family-level summary
    ann["is_sig_fdr"] = ann["q_fdr"] < args.sig_q
    fam = (
        ann.groupby("pirate_family", dropna=False)
           .agg(
               n_alleles=("variant", "count"),
               n_sig_fdr=("is_sig_fdr", "sum"),
               min_p=(pcol, "min"),
               min_q=("q_fdr", "min"),
               top_variant=("variant", lambda x: x.iloc[0]),
               top_beta=("beta", lambda x: x.iloc[0]),
               top_af=("af", lambda x: x.iloc[0]),
           )
           .reset_index()
           .sort_values(["min_q", "min_p"], ascending=[True, True])
    )
    fam_out = f"{args.out_prefix}.family_summary.tsv"
    fam.to_csv(fam_out, sep="\t", index=False)

    print(f"Wrote: {ann_out}")
    print(f"Wrote: {fam_out}")
    missing = ann["pirate_family"].isna().sum()
    if missing:
        print(f"WARNING: {missing} GWAS variants did not map to PIRATE alleles (ID mismatch?)", file=sys.stderr)


if __name__ == "__main__":
    main()
