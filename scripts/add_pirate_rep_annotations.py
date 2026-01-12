#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def parse_pirate_rep_headers(rep_faa: str) -> pd.DataFrame:
    """
    Parse PIRATE representative_sequences.faa headers.
    Example header:
    >g00001;representative_genome=100014_FSIS1605809;locus_tag=...;gene_name=soj;gene_product=...;number_genomes=2327
    """
    rows = []
    with open(rep_faa, "r") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            h = line[1:].strip()

            # Split on ';' into tokens
            toks = h.split(";")
            pirate_family = toks[0].strip()

            kv = {}
            for t in toks[1:]:
                if "=" in t:
                    k, v = t.split("=", 1)
                    kv[k.strip()] = v.strip()

            # Some PIRATE headers may use slightly different keys; be tolerant
            rep_genome = kv.get("representative_genome") or kv.get("rep_genome")
            locus_tag  = kv.get("locus_tag") or kv.get("locus") or kv.get("rep_locus_tag")
            gene_name  = kv.get("gene_name") or kv.get("gene") or kv.get("name")
            product    = kv.get("gene_product") or kv.get("product") or kv.get("annotation")
            n_genomes  = kv.get("number_genomes") or kv.get("n_genomes") or kv.get("genomes")

            # Coerce numeric if possible
            try:
                n_genomes = int(n_genomes) if n_genomes is not None else None
            except ValueError:
                n_genomes = None

            rows.append({
                "pirate_family": pirate_family,
                "pirate_rep_genome": rep_genome,
                "pirate_rep_locus_tag": locus_tag,
                "pirate_gene_name": gene_name,
                "pirate_gene_product": product,
                "pirate_number_genomes": n_genomes,
                "pirate_rep_header": h,
            })

    df = pd.DataFrame(rows)

    # Deduplicate defensively (should be 1 rep per family)
    df = df.drop_duplicates("pirate_family", keep="first").copy()
    return df


def main():
    ap = argparse.ArgumentParser(
        description="Merge PIRATE representative header annotations into family-level GWAS table."
    )
    ap.add_argument("--rep-faa", required=True, help="PIRATE_out/representative_sequences.faa")
    ap.add_argument("--families", required=True, help="Family-level table (e.g. birds_vs_mammal.family_level.annotated.clean.tsv)")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--make-bestnames", action="store_true",
                    help="Add best_gene/best_product using ref_* as primary and pirate_* as fallback")
    args = ap.parse_args()

    fam = pd.read_csv(args.families, sep="\t")
    rep = parse_pirate_rep_headers(args.rep_faa)

    # Avoid annoying _x/_y suffixes if someone already merged partial PIRATE fields earlier
    pirate_cols = [c for c in fam.columns if c.startswith("pirate_")]
    fam2 = fam.drop(columns=pirate_cols, errors="ignore").copy()

    # ensure merge key is present
    if "pirate_family" not in fam2.columns and "pirate_family" in fam.columns:
        fam2 = fam2.copy()
        fam2["pirate_family"] = fam["pirate_family"].values
    out = fam2.merge(rep, on="pirate_family", how="left")

    if args.make_bestnames:
        # best_product
        if "ref_product" in out.columns:
            out["best_product"] = out["ref_product"].where(out["ref_product"].notna(), out["pirate_gene_product"])
        else:
            out["best_product"] = out["pirate_gene_product"]

        # best_gene
        if "ref_gene" in out.columns:
            out["best_gene"] = out["ref_gene"].where(out["ref_gene"].notna(), out["pirate_gene_name"])
        else:
            out["best_gene"] = out["pirate_gene_name"]

    out.to_csv(args.out, sep="\t", index=False)

    # Small sanity summary to stderr
    n_total = len(out)
    n_rep = out["pirate_gene_product"].notna().sum() if "pirate_gene_product" in out.columns else 0
    sys.stderr.write(f"[ok] Wrote {args.out}\n")
    sys.stderr.write(f"[info] Families: {n_total}\n")
    sys.stderr.write(f"[info] PIRATE rep products present: {n_rep}/{n_total}\n")

    if args.make_bestnames:
        n_bestprod = out["best_product"].notna().sum()
        sys.stderr.write(f"[info] best_product present: {n_bestprod}/{n_total}\n")


if __name__ == "__main__":
    main()
