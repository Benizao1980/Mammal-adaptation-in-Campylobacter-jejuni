# Pangenome analysis of the Zang et al. dataset

---

## Overview and rationale

To characterise gene content diversity, structural variation, and host-associated genomic features across the *Zang et al.* genome collection, we performed a comprehensive pangenome analysis using **PIRATE** (Pangenome Iterative Refinement and Threshold Evaluation). PIRATE is designed for large bacterial genome collections and explicitly models allelic diversity, gene fragmentation, fission/fusion events, and gene duplication by clustering genes across multiple amino-acid identity thresholds.

This document combines:
- a **step-by-step computational walkthrough**,
- **inline commentary on observed results**, and
- **reproducibility metadata**, including software versions and references.

All analyses were performed on the full dataset and form the basis for downstream comparative, evolutionary, and host-association analyses.

---

## Genome dataset

The analysis comprised **2,327 whole-genome assemblies** derived from the *Zang et al.* dataset. These genomes represent isolates sampled across multiple host species and epidemiological contexts. All assemblies were processed using a uniform annotation and pangenome workflow to minimise technical artefacts.

*Commentary:*  
The large sample size and host diversity provide sufficient power to characterise both rare accessory genes and structural variation within highly prevalent (core) gene families.

---

## Genome annotation with Prokka

Consistent gene annotation is critical for pangenome analysis. Differences in gene calling or annotation standards can artificially inflate accessory gene counts or introduce spurious gene fragmentation. All genomes were therefore annotated de novo using **Prokka**, ensuring consistent gene prediction and functional annotation across the dataset.

### Input preparation

```bash
mkdir -p input
```

```bash
find contigs -maxdepth 1 -type f \
  \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fas" \) \
  -printf "%f\n" | sort > input/InputFiles
```

```bash
wc -l input/InputFiles
head input/InputFiles
```

### High-throughput annotation

Prokka was executed as a SLURM array job, processing one genome per task (up to 500 concurrent jobs):

```slurm
#!/bin/bash
#SBATCH --job-name=prokka_batch
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --array=1-500

source ~/.bashrc
conda activate prokka

WORKDIR="/home/u12/bpascoe/Zang_cdtB"
cd "$WORKDIR"

INPUT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input/InputFiles)
[ -z "$INPUT_FILE" ] && exit 0

PREFIX="${INPUT_FILE%.*}"
OUTDIR="output/${PREFIX}"
mkdir -p "$OUTDIR"

prokka --outdir "$OUTDIR" --prefix "$PREFIX" --cpus 4 --force "contigs/${INPUT_FILE}"
```

*Commentary:*  
This produced one annotated GFF per genome. These GFFs form the direct input for PIRATE and preserve consistent gene boundaries across the dataset.

---

## Pangenome reconstruction with PIRATE

Unlike single-threshold clustering methods, PIRATE iteratively clusters genes across descending amino-acid identity thresholds. This allows separation of true orthologues from divergent alleles and paralogues while retaining information on structural variation.

### Execution

```slurm
#!/bin/bash
#SBATCH --job-name=pirate_zang
#SBATCH --output=slurm_logs/pirate_%A.out
#SBATCH --error=slurm_logs/pirate_%A.err
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=72:00:00
#SBATCH --mem=250G

source ~/.bashrc
conda activate pirate

PIRATE -i pirate_gff -o PIRATE_out -t 16 -a -s "80,85,90,92,94,95,96,97,98"
```

*Commentary:*  
Identity thresholds spanning 80–98% are appropriate for *Campylobacter*, capturing both deep divergence and fine-scale allelic variation.

---

## Global pangenome structure (results)

The inferred pangenome comprised **5,387 gene families** across 2,327 genomes.

Key observations:
- **682 gene families** exhibited allelic diversity (>1 allele)
- **1,701 gene families** showed evidence of gene fission or fusion
- **998 gene families** exhibited gene duplication or loss

*Interpretation:*  
These results indicate a highly dynamic genome architecture. Structural variation is not restricted to rare accessory genes but is also prevalent among widely distributed gene families.

---

## Gene frequency and structural variation

Gene families were stratified by prevalence to assess how structural variation scales with frequency.

| % isolates | # gene clusters | >1 allele | fission/fusion | multicopy |
|-----------:|----------------:|----------:|---------------:|----------:|
| 0–10%      | 3,518           | 249       | 325            | 83        |
| 10–25%     | 228             | 118       | 143            | 105       |
| 25–50%     | 128             | 60        | 92             | 68        |
| 50–75%     | 68              | 29        | 55             | 41        |
| 75–90%     | 50              | 21        | 41             | 35        |
| 90–95%     | 29              | 6         | 25             | 17        |
| 95–100%    | 1,366           | 199       | 1,020          | 649       |

*Interpretation:*  
Rare gene families dominate the accessory genome, consistent with an open pangenome. However, even near-core genes frequently show fission/fusion and copy-number variation, indicating ongoing structural evolution within the core genome.

---

## Visualisation

```bash
python3 pirate_plotting.py --pirate-out PIRATE_out --outdir pirate_plots --all --palette life_aquatic_blues
```

These plots summarise core/accessory structure, allelic diversity, and structural variation across the pangenome.

---

## Host-associated gene analysis

```bash
python3 pirate_host_association.py \
  --pirate-out PIRATE_out \
  --meta meta.tsv \
  --outdir pirate_host_assoc \
  --group-col host \
  --rarefaction-per-group \
  --assoc \
  --min-prevalence 0.01 \
  --max-prevalence 0.99
```

*Interpretation:*  
This analysis identifies accessory gene families whose distribution is significantly associated with host category, while controlling for prevalence and unequal sampling depth.

---

## Core-genome phylogeny with IQ-TREE

To infer a high-quality maximum-likelihood (ML) phylogeny from the PIRATE core alignment, we used IQ-TREE2 on the nucleotide core alignment (PIRATE_out/core_alignment.fasta) under a GTR+F+I+G4 model, with branch support from both ultrafast bootstrap and SH-aLRT. The resulting tree provides the fixed topology required by ClonalFrameML.

```bash
iqtree2 \
  -s PIRATE_out/core_alignment.fasta \
  -m GTR+F+I+G4 \
  -T 40 \
  -B 1000 \
  --alrt 1000 \
  --prefix zang_core \
  --seed 12345 \
  --safe
```

*Commentary:*
This run used IQ-TREE2 v2.3.6, inferred a tree from 2,327 sequences with an alignment length of 1,371,024 bp, and used a fixed seed for reproducibility. 

### Key outputs

- zang_core.treefile — **ML tree** (Newick; used for ClonalFrameML)
- zang_core.iqtree / zang_core.log — **run record and model details**
- zang_core.ckp.gz — **checkpoint file for safe/restart mode**
- zang_core.mldist — **pairwise ML distances** (optional downstream use)

## Recombination-aware inference with ClonalFrameML (CF-ML)

We next used ClonalFrameML to infer recombinant tracts on the fixed IQ-TREE topology, enabling downstream analyses on a recombination-masked alignment (or recombination-aware branch lengths, depending on the application).

### Inputs

- Tree: zang_core.treefile (IQ-TREE output)
- Alignment: PIRATE_out/core_alignment.fasta (core alignment)

Example SLURM job
```slurm
#!/bin/bash
#SBATCH --job-name=cfml_core
#SBATCH --output=slurm_logs/cfml_%A.out
#SBATCH --error=slurm_logs/cfml_%A.err
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
#SBATCH --mem=160G

set -euo pipefail
source ~/.bashrc
conda activate clonalframeml

cd /groups/cooperma/bpascoe/Zang_cdtB
mkdir -p slurm_logs cfml_snps_out

TREE="iqtree/zang_core.treefile"
ALN="PIRATE_out/core_alignment.fasta"
OUT="cfml_core_out/cfml_core"

test -s "$TREE"
test -s "$ALN"

ClonalFrameML "$TREE" "$ALN" "$OUT"
```

### Key CF-ML outputs (by prefix)

- Given OUT="cfml_core_out/cfml_core", ClonalFrameML will produce (among others):
- cfml_core_out/cfml_core.labelled_tree.newick
- cfml_core_out/cfml_core.importation_status.txt
- cfml_core_out/cfml_core.em.txt

These are the required inputs for recombination masking in the next step.

## Mask recombination tracts (CF-ML → masked alignment)

To generate a recombination-masked core alignment for downstream phylogenetic inference and association testing, we applied a masking script that converts CF-ML inferred recombinant segments into per-isolate masked sites (default mask symbol: N).

We used the updated script cfml-maskrc_updated.py, which adds robustness and useful outputs (ID normalization, interval merging, optional ancestral masking, per-isolate masking metrics, optional SVG plotting, and rectangular alignment checks). 


Minimal masking run (extant recombination only)
```bash
python3 cfml-maskrc_updated.py \
  cfml_core_out/cfml_core \
  --aln PIRATE_out/core_alignment.fasta \
  --out cfml_core_out/core_alignment.masked \
  --symbol N \
  --regions cfml_core_out/recomb_regions.tsv \
  --metrics cfml_core_out/recomb_masking_metrics.tsv
```

This produces:

- cfml_core_out/core_alignment.masked.fasta (masked alignment) 
- cfml-maskrc_updated
- cfml_core_out/recomb_regions.tsv (per-isolate recombinant intervals) 
- cfml-maskrc_updated
- cfml_core_out/recomb_masking_metrics.tsv (bp masked, segment counts, fraction masked; can embed CF-ML params if detected) 

*Optional:* also mask ancestral segments

If you want to mask segments inferred on internal branches (not just leaf-specific “extant” imports), add:

```
  --mask-ancestral
```

*Commentary:*
Masking ancestral segments is more aggressive and may be appropriate if your downstream method is sensitive to deep recombination signal; for many applications, extant-only masking is a good default.

Optional: generate an SVG summary plot

```bash
python3 cfml-maskrc_updated.py \
  cfml_core_out/cfml_core \
  --aln PIRATE_out/core_alignment.fasta \
  --out cfml_core_out/core_alignment.masked \
  --svg cfml_core_out/recombination_map.svg \
  --svgsize 1200x800 \
  --svgcolour black \
  --consensus
```

This yields a compact “recombination barcode” view across isolates plus an optional consensus/hotspot track. 

## Downstream: phylogeny on recombination-masked alignment (recommended)

Once masking is complete, re-infer an ML tree on the masked alignment to obtain branch lengths/topology less driven by homologous recombination:

```bash
iqtree2 \
  -s cfml_core_out/core_alignment.masked.fasta \
  -m GTR+F+I+G4 \
  -T 40 \
  -B 1000 \
  --alrt 1000 \
  --prefix zang_core.masked \
  --seed 12345 \
  --safe
```

*Interpretation:*
This “masked-core” tree is typically the best default for host-association, trait mapping, and selection/structure analyses where recombination can inflate apparent homoplasy or distort branch lengths.

---

## Reproducibility, software versions, and resources

All analyses were performed on an HPC cluster using SLURM.

### Software
- Prokka v1.14.x — https://github.com/tseemann/prokka
- PIRATE v1.0.x — https://github.com/SionBayliss/PIRATE
- PIRATE additional analyses scripts <link><link><link>
- IQ-TREE2 v2.3.6 - <link>
- ClonalFrameML v1.0.x - <link>
- ClonalFrameML masking script - <link>
- Python v3.9+ — https://www.python.org
- Conda — https://docs.conda.io

### Reproducibility notes
- All genomes processed with identical parameters
- FigShare link to genomes -
- PubMLST shared proejct - 
- microreact link to tree and isoalte data -
- - Conda environments used throughout

### Please cite: 
TBC <link to manuscript>

---

## Summary

This integrated workflow documents a reproducible pangenome analysis of the *Zang et al.* dataset, combining detailed computational steps with inline biological interpretation. The results demonstrate an open and structurally dynamic pangenome with extensive allelic diversity and host-associated gene content variation.
