#!/usr/bin/env python3

import pandas as pd

def main(gwas, annotated, out):
    g = pd.read_csv(gwas, sep="\t")
    a = pd.read_csv(annotated, sep="\t")

    df = g.merge(
        a[["variant","pirate_family","product","ref_gene"]],
        on="variant",
        how="left"
    )

    pcol = "lrt-pvalue"
    df["abs_beta"] = df["beta"].abs()

    fam = (
        df.groupby("pirate_family", dropna=False)
          .agg(
              n_alleles_tested=("variant","count"),
              n_alleles_sig=(pcol, lambda x: (x < 0.05).sum()),
              min_p=(pcol,"min"),
              max_abs_beta=("abs_beta","max"),
              mean_beta=("beta","mean"),
              direction=("beta", lambda x:
                  "bird" if x.mean() > 0 else "mammal"
              ),
              product=("product", lambda x: x.dropna().mode().iloc[0] if len(x.dropna()) else "hypothetical"),
              ref_gene=("ref_gene", lambda x: x.dropna().mode().iloc[0] if len(x.dropna()) else None)
          )
          .reset_index()
          .sort_values("min_p")
    )

    fam.to_csv(out, sep="\t", index=False)

if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3])
