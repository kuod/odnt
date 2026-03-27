# 00 — Perturb-seq: Overview

## What problem does Perturb-seq solve?

Traditional pooled CRISPR screens (e.g., genome-wide dropout screens) read out a single phenotype — typically cell survival or a fluorescent reporter. They answer "which gene, when knocked out, changes this one thing?" but tell you nothing about *mechanism*.

Single-cell RNA-seq (scRNA-seq) gives you transcriptome-wide mechanistic signal in individual cells, but is purely observational — you cannot directly assign cause to what you see.

**Perturb-seq combines both:** every cell carries an identifiable CRISPR perturbation *and* a full transcriptome readout. The result is a systematic, causal map from genetic perturbation to transcriptional response, at single-cell resolution.

---

## The core technical insight: two libraries, one cell barcode

A Perturb-seq experiment generates two sequencing libraries from the same cells:

```
Each cell in a droplet
 ├── GEX library      → standard 3' (or 5') scRNA-seq
 │                       reads out ~20,000 genes
 └── Guide capture    → short amplicon from a unique barcode
     library            appended to each guide RNA
                         reads out which guide is in this cell
```

Both libraries are sequenced, demultiplexed by the shared **cell barcode** (10x Genomics 16-nt barcode + UMI), and merged in silico. The result: for every cell, you know (1) its entire transcriptome and (2) which gene was perturbed.

---

## CRISPRa vs. CRISPRi vs. CRISPRko

All three use Cas9 or a Cas9 variant guided to a genomic locus by a single guide RNA (sgRNA). They differ in *what* they do there:

| Modality | Cas9 variant | Effect on target gene | Notes |
|---|---|---|---|
| **CRISPRko** (knockout) | WT Cas9 (nuclease) | Cuts DNA → indels → loss of function | Permanent; editing efficiency varies |
| **CRISPRi** (interference) | dCas9-KRAB | Represses transcription | Reversible; no DNA cutting; targets TSS |
| **CRISPRa** (activation) | dCas9-VPR or dCas9-SAM | Activates transcription | Reversible; targets TSS; overexpresses gene |

**For Perturb-seq, CRISPRi and CRISPRa are generally preferred** because:
- Effects are graded and proportional (vs. binary cut/no-cut in CRISPRko)
- No indel heterogeneity — all cells with the same guide have the same molecular perturbation
- More reproducible across guides targeting the same gene

CRISPRko is appropriate when you specifically need protein loss (e.g., surface proteins, secreted factors).

---

## Perturb-seq vs. related methods

| Method | GEX | Protein | Chromatin | Notes |
|---|---|---|---|---|
| **Perturb-seq** | Yes (3'/5') | No | No | Original; most common |
| **CROP-seq** | Yes | No | No | Nearly identical to Perturb-seq; slightly different guide capture |
| **ECCITE-seq** | Yes | Yes (CITE-seq) | No | Adds surface protein readout |
| **Perturb-ATAC** | No | No | Yes | Chromatin accessibility instead of RNA |
| **Mosaic-seq** | Yes | No | No | Cis-regulatory element perturbations |
| **CRISP-seq** | Yes | No | No | Earlier protocol; less common now |

All of these read out the same information — "which perturbation is in this cell" — using the same feature barcoding principle. The analysis concepts in this course apply to all variants.

---

## Full pipeline schematic

```
Experimental
  ┌──────────────────────────────────────────────────────────────────┐
  │  Design gRNA library  →  Clone into lentiviral vector            │
  │  Transduce Cas9+ cells at low MOI (~0.3)                        │
  │  Capture cells in droplets (10x Genomics Chromium)              │
  └──────────────────────────────────────────────────────────────────┘
           │
           │ Two PCR amplification tracks per cell
           ▼
  ┌───────────────────┐   ┌────────────────────────────┐
  │   GEX library     │   │   Guide capture library    │
  │  (whole mRNA)     │   │  (guide barcode amplicon)  │
  └────────┬──────────┘   └────────────┬───────────────┘
           │ Illumina sequencing       │
           ▼                           ▼
  ┌──────────────────────────────────────────────────────────────────┐
  │  Cell Ranger / STARsolo / kb-python                              │
  │  → GEX count matrix  +  guide UMI count matrix                  │
  │  → Linked by shared cell barcode                                 │
  └──────────────────────────────────────────────────────────────────┘
           │
           ▼
  ┌──────────────────────────────────────────────────────────────────┐
  │  Quality control (05_quality_control.ipynb)                      │
  │  Guide assignment (06_guide_assignment.ipynb)                    │
  │  Normalization + dim reduction (07_normalization_dimred.ipynb)   │
  └──────────────────────────────────────────────────────────────────┘
           │
           ▼
  ┌──────────────────────────────────────────────────────────────────┐
  │  Perturbation effects (08)  →  DE per perturbation vs. NTC      │
  │  Mixscape (09)              →  Filter escape cells               │
  │  Visualization (10)         →  Volcano, heatmap, embedding       │
  │  Genetic interactions (11)  →  Additive model, epistasis         │
  │  GRN + prediction (12)      →  TF activity, GEARS               │
  └──────────────────────────────────────────────────────────────────┘
```

---

## When to use Perturb-seq

**Good fit:**
- You want mechanistic understanding of many perturbations at once
- You care about which pathways / gene programs a perturbation affects
- You are studying context-dependent effects (different cell states)
- You want to map genetic interactions transcriptome-wide

**Not the right tool:**
- You need readout of a non-transcriptional phenotype (protein levels, morphology, metabolomics)
- Budget/scale constraints rule out scRNA-seq costs (~$1–5 per cell sequenced)
- You need high statistical power per perturbation — pseudobulk DE with 50–200 cells/perturbation is underpowered for subtle effects

---

## Key papers

| Paper | What it introduced |
|---|---|
| Dixit et al. 2016, *Cell* | Original Perturb-seq protocol (CRISPRko, dendritic cells) |
| Adamson et al. 2016, *Cell* | CRISPRi Perturb-seq; UPR dissection |
| Norman et al. 2019, *Science* | CRISPRa; genetic interactions; Norman dataset used in this course |
| Replogle et al. 2022, *Cell* | Genome-scale Perturb-seq; Replogle dataset used in this course |
| Heumos et al. 2025, *Nature Methods* | pertpy framework |
| Kernfeld et al. 2025, *Nature Methods* | Benchmarking of perturbation prediction models |

---

**Next:** [01_experimental_design.md](01_experimental_design.md)
