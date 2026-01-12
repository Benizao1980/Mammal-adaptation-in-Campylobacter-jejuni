#!/usr/bin/env python3
import argparse, re
import pandas as pd

def parse_faa_headers(faa_path: str, strain_label: str) -> pd.DataFrame:
    """
    Parse headers from a protein FASTA into a lookup table:
    sseqid -> gene/product/locus_tag (best effort).
    """
    rows = []
    current_id = None
    current_desc = None

    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                hdr = line[1:].strip()
                parts = hdr.split(None, 1)
                sid = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                rows.append({"sseqid": sid, "ref_strain": strain_label, "ref_header": desc})
    df = pd.DataFrame(rows)

    # Best-effort parsing patterns (handles many RefSeq/Prokka-like styles)
    # Try to capture locus_tag, gene, product if present
    def get(pattern, text):
        m = re.search(pattern, text)
        return m.group(1) if m else None

    df["ref_locus_tag"] = df["ref_header"].apply(lambda x: get(r"(?:locus_tag=|locus_tag:)(\S+)", x) or get(r"\b([A-Z]{2,}_[0-9]+)\b", x))
    df["ref_gene"]      = df["ref_header"].apply(lambda x: get(r"(?:gene=|gene:)(\S+)", x))
    df["ref_product"]   = df["ref_header"].apply(lambda x: get(r"(?:product=|product:)(.+)$", x))

    # Fallback: if "product=" not present, use whole header as product-ish text
    df.loc[df["ref_product"].isna(), "ref_product"] = df.loc[df["ref_product"].isna(), "ref_header"].replace("", pd.NA)

    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast", required=True, help="BLAST TSV (qseqid sseqid pident length qlen slen evalue bitscore)")
    ap.add_argument("--ref-faa", required=True, nargs="+", help="Reference protein FASTAs (published annotations)")
    ap.add_argument("--ref-label", required=True, nargs="+", help="Labels for refs (same order as --ref-faa)")
    ap.add_argument("--out", required=True, help="Output annotations TSV")
    ap.add_argument("--min-pident", type=float, default=70.0)
    ap.add_argument("--min-qcov", type=float, default=0.80)
    args = ap.parse_args()

    if len(args.ref_faa) != len(args.ref_label):
        raise SystemExit("ERROR: --ref-faa and --ref-label must be same length")

    # Load blast
    b = pd.read_csv(args.blast, sep="\t", header=None,
                    names=["qseqid","sseqid","pident","length","qlen","slen","evalue","bitscore"])
    b["qcov"] = b["length"] / b["qlen"]

    # Filter to decent hits
    b = b[(b["pident"] >= args.min_pident) & (b["qcov"] >= args.min_qcov)].copy()

    # Choose best hit per qseqid (bitscore desc, then evalue asc)
    b = b.sort_values(["qseqid","bitscore","evalue"], ascending=[True, False, True])
    best = b.drop_duplicates("qseqid", keep="first").copy()

    # Build header lookup across refs
    ref_tables = []
    for faa, lab in zip(args.ref_faa, args.ref_label):
        ref_tables.append(parse_faa_headers(faa, lab))
    ref = pd.concat(ref_tables, ignore_index=True)

    out = best.merge(ref, on="sseqid", how="left")

    # Tidy / select
    out = out[[
        "qseqid","sseqid","ref_strain","ref_locus_tag","ref_gene","ref_product",
        "pident","qcov","evalue","bitscore"
    ]].copy()

    out.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out} (rows={len(out)})")

if __name__ == "__main__":
    main()
