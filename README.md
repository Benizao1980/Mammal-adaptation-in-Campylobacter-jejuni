# Panopticon

**Panopticon** is a modular analysis toolbox for exploring **pangenome structure, diversity, and host association** using gene presence/absence data and associated metadata.

It is designed for **downstream analysis** of pangenomes inferred using tools such as **PIRATE**, **Panaroo**, or similar workflows, with a focus on:
- accessory genome structure,
- host and source association,
- lineage-aware statistical testing, and
- interpretable, publication-ready outputs.

Panopticon prioritises **clarity, reproducibility, and statistical rigour**, while remaining flexible across pathogens and datasets.

---

## Key concepts

Panopticon treats the pangenome as a **high-dimensional ecological system**, where genomes are samples, genes are features, and metadata (host, source, lineage) define structured variation.

Rather than relying solely on clustering or visual inspection, Panopticon integrates:
- distance-based ordination,
- rarefaction and accumulation analyses,
- formal hypothesis testing (e.g. PERMANOVA),
- and biologically interpretable summaries of gene content.

This makes it particularly suited for questions around:
- host adaptation,
- zoonotic transmission,
- lineage vs ecology effects,
- and accessory genome diversification.

---

## Core functionality

### Accessory genome structure
- Presence/absence matrix handling (binary gene families)
- Jaccard and related distance metrics
- Ordination methods (PCA, MDS, UMAP)
- Lineage-aware visualisation

### Pangenome composition
- Core / soft-core / shell / cloud summaries
- Gene frequency distributions
- Per-genome gene content statistics

### Rarefaction and accumulation
- Pangenome and core genome accumulation curves
- Stratified rarefaction by host, source, or other metadata
- Confidence intervals via resampling

### Statistical testing
- PERMANOVA for testing host/source effects on gene content
- Support for **restricted permutations** (e.g. within clonal complexes)
- PERMDISP to assess dispersion effects
- Pairwise comparisons with multiple-testing correction

### Annotation and rationalisation (planned)
- Mapping gene families to representative loci
- Collapsing gene families to functional categories (e.g. COGs)
- Functional enrichment analyses across hosts or sources

---

## Inputs

Panopticon is input-agnostic but expects:

### Gene presence/absence matrix
- Rows: samples / genomes  
- Columns: gene families  
- Values: binary (0/1)

This can be derived from:
- PIRATE (`binary_presence_absence.fasta`)
- Panaroo gene presence/absence tables
- Other equivalent formats

### Metadata table
A tab-delimited file containing at minimum:
