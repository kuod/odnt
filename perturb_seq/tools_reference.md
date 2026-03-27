# Perturb-seq Tools Reference

A comparative overview of tools used across the Perturb-seq analysis pipeline, organized by stage. Each entry covers what the tool does, its strengths and weaknesses, and how it compares to alternatives.

---

## Table of Contents

1. [Raw Data Processing](#1-raw-data-processing)
2. [Guide RNA Assignment](#2-guide-rna-assignment)
3. [Quality Control](#3-quality-control)
4. [Single-Cell Analysis Frameworks](#4-single-cell-analysis-frameworks)
5. [Perturbation Effect Analysis](#5-perturbation-effect-analysis)
6. [Escape Cell Detection](#6-escape-cell-detection)
7. [Differential Expression](#7-differential-expression)
8. [Transcription Factor Activity](#8-transcription-factor-activity)
9. [Perturbation Prediction & GRN](#9-perturbation-prediction--grn)
10. [Databases & Data Resources](#10-databases--data-resources)

---

## 1. Raw Data Processing

Tools that convert FASTQ files into cell × gene count matrices.

### Cell Ranger (10x Genomics)

| | |
|---|---|
| **What it does** | Official 10x pipeline; processes GEX + CRISPR guide capture libraries jointly via `--feature-ref` |
| **Language** | Proprietary (Rust/Python wrapper) |
| **License** | Free for research; no source code |

**Pros**
- Native CRISPR guide capture support — no extra steps to merge GEX and guide matrices
- Best tested with 10x Chromium data; highest compatibility
- Automated cell filtering (EmptyDrops-based)
- Comprehensive HTML QC report; saturation curves, knee plots

**Cons**
- Closed source; cannot inspect or modify the pipeline
- Slow (~4–6 hr for 250M reads) and memory-hungry (32–64 GB RAM)
- Locked to 10x chemistry; poor support for non-standard protocols
- Requires a license agreement for commercial use

**Differentiator:** The only tool with fully integrated, validated CRISPR feature barcoding. If you're using 10x, use this unless you have a specific reason not to.

---

### STARsolo

| | |
|---|---|
| **What it does** | Single-cell mode of the STAR aligner; produces Cell Ranger-compatible output |
| **Language** | C++ |
| **License** | MIT |

**Pros**
- 2–4× faster than Cell Ranger
- Open source; fully auditable
- Supports many chemistries (10x v2/v3, Drop-seq, inDrops, etc.)
- Outputs BAM + splice junction info useful for isoform analysis

**Cons**
- No native CRISPR guide capture support — guide libraries must be processed separately with a custom script
- Requires manual cell barcode whitelist management
- QC reporting is minimal compared to Cell Ranger

**Differentiator:** Best choice when you need open-source flexibility, non-10x chemistries, or speed. Requires extra work to handle guide capture.

---

### kb-python / kallisto|bustools

| | |
|---|---|
| **What it does** | Pseudoalignment-based quantification; `kb count --workflow kite` handles CRISPR feature barcoding |
| **Language** | Python wrapper (C++ core) |
| **License** | BSD |

**Pros**
- Fastest of the three (~15–30 min for 250M reads)
- `kite` workflow natively handles feature barcoding (guides, antibodies)
- Fully open source and scriptable
- Works well for custom/novel chemistries
- Best for ultra-large datasets where Cell Ranger's runtime is prohibitive

**Cons**
- Pseudoalignment: does not produce a BAM file (no splicing analysis)
- Less thoroughly validated for Perturb-seq specifically
- Cell filtering not built in — requires downstream tools (e.g., `barcodeutils`)
- Community support smaller than Cell Ranger

**Differentiator:** Fastest option; ideal for high-throughput or compute-constrained settings. The go-to for non-10x chemistries and for feature barcoding flexibility.

| | Cell Ranger | STARsolo | kb-python |
|---|---|---|---|
| Speed | Slow | Fast | Fastest |
| CRISPR capture | Native | Manual | Native (kite) |
| Open source | No | Yes | Yes |
| Non-10x support | Poor | Good | Excellent |
| BAM output | Yes | Yes | No |

---

## 2. Guide RNA Assignment

Tools that map guide UMI count matrices to per-cell guide identities.

### pertpy — GuideRNAAssignment

| | |
|---|---|
| **What it does** | Poisson mixture model to separate true guide captures from background noise |
| **Part of** | `pertpy` (scverse ecosystem) |

**Pros**
- Probabilistic: provides assignment confidence scores, not just binary calls
- Handles ambient guide contamination explicitly
- Integrates directly with AnnData; output populates `obs` columns
- Actively maintained (scverse consortium)

**Cons**
- Requires the raw guide UMI count matrix (separate from GEX matrix)
- API has changed across pertpy versions — check documentation for your version
- Slower than threshold-based methods for very large datasets

**Differentiator:** Best default choice for assignment when you have the raw guide UMI matrix. Probabilistic model reduces false positives from ambient contamination.

---

### Threshold-based assignment (manual)

| | |
|---|---|
| **What it does** | Assign the top-UMI guide per cell if it exceeds a minimum count; flag cells with two qualifying guides as multiplets |

**Pros**
- Completely transparent; no model assumptions
- Fast; works on any matrix size
- Easy to tune and inspect

**Cons**
- Binary; no confidence score
- Fixed threshold does not adapt to per-cell sequencing depth variation
- Underperforms probabilistic methods when ambient contamination is high

**Differentiator:** Good as a sanity check and for simple datasets with high guide capture depth. Use as a fallback when pertpy is unavailable.

---

### crispat

| | |
|---|---|
| **What it does** | Probabilistic guide assignment using a Bayesian model; published in *Bioinformatics* 2024 |
| **Language** | Python |

**Pros**
- Explicitly models guide-specific capture efficiency (some guides are captured less efficiently)
- Provides posterior probabilities for multiplet detection
- Well-validated on simulated and real data in the publication

**Cons**
- Newer, smaller community than pertpy
- Requires additional installation steps
- Documentation less complete than pertpy

**Differentiator:** Best choice when guides have highly variable capture efficiency or when you need rigorous probabilistic multiplet detection.

| | pertpy | Threshold | crispat |
|---|---|---|---|
| Model | Poisson mixture | None | Bayesian |
| Confidence scores | Yes | No | Yes (posterior) |
| Multiplet handling | Basic | Binary | Rigorous |
| Maturity | High | N/A | Medium |

---

## 3. Quality Control

### scanpy — `sc.pp.calculate_qc_metrics`

Standard scRNA-seq QC. No Perturb-seq-specific logic — the key is to run it **per perturbation**, not globally, to detect perturbation-correlated QC artifacts. See notebook 05.

---

## 4. Single-Cell Analysis Frameworks

### scanpy

| | |
|---|---|
| **What it does** | Core Python single-cell analysis: normalization, HVG selection, PCA, UMAP, clustering, DE |
| **Language** | Python |
| **License** | BSD |

**Pros**
- De facto standard for Python single-cell analysis
- Excellent documentation and large community
- Tight AnnData integration
- Fast sparse matrix operations

**Cons**
- DE methods (Wilcoxon, t-test on single cells) are **not appropriate for Perturb-seq** — use pseudobulk instead
- UMAP implementation is stochastic; results not bit-reproducible without fixed seed

**Differentiator:** Use scanpy for everything except DE (use PyDESeq2 pseudobulk there). The backbone of this entire course.

---

### pertpy

| | |
|---|---|
| **What it does** | Perturb-seq-specific analysis: guide assignment, Mixscape, perturbation distances, Augur, genetic interactions |
| **Language** | Python |
| **License** | BSD |
| **Paper** | Heumos et al. 2025, *Nature Methods* |

**Pros**
- Purpose-built for Perturb-seq; one package covers most analysis steps
- Part of the scverse ecosystem — tight AnnData/scanpy integration
- Actively maintained; published in Nature Methods 2025
- Covers Mixscape, E-distance, Augur, guide assignment, and more

**Cons**
- API changes frequently between versions — code written for v0.6 may break on v0.8+
- Some methods (e.g., Mixscape LPC embedding) require specific input formats that are not well documented
- Heavier dependency footprint than scanpy alone

**Differentiator:** The single most important Perturb-seq-specific Python package. Install this first when setting up a new analysis environment.

---

### muon / MuData

| | |
|---|---|
| **What it does** | Multi-modal AnnData container; holds GEX + guide capture + surface proteins in one object |
| **Language** | Python |
| **License** | BSD |

**Pros**
- Clean API for multi-modal data: `mdata["rna"]`, `mdata["crispr"]`
- `mu.pp.intersect_obs()` handles cell barcode alignment between modalities
- Compatible with all scanpy/pertpy functions (each modality is a standard AnnData)

**Cons**
- Overhead for simple single-modality analyses
- Not all pertpy functions accept MuData directly — may need to extract the RNA modality

**Differentiator:** Use when working with raw Cell Ranger output (which separates GEX and guide matrices). For pre-processed public datasets (Norman, Replogle), MuData is usually not needed.

---

## 5. Perturbation Effect Analysis

### PyDESeq2 (pseudobulk DE)

| | |
|---|---|
| **What it does** | Python reimplementation of DESeq2's negative binomial model for bulk/pseudobulk DE |
| **Language** | Python |
| **License** | MIT |

**Pros**
- Statistically correct for Perturb-seq: aggregating to pseudobulk + NB model avoids pseudoreplication
- Results concordant with R DESeq2 for most use cases
- Pure Python — no R dependency
- Well-maintained by the scverse community

**Cons**
- Slow for genome-scale screens (2,000+ perturbations × 20,000 genes takes hours)
- Requires ≥2 pseudobulk replicates per condition — screens without replicate structure cannot use it
- Less tested than R's DESeq2 for edge cases (very small counts, outliers)

**Differentiator:** The correct DE approach for Perturb-seq. If you use per-cell Wilcoxon instead, your p-values are invalid.

---

### SCEPTRE

| | |
|---|---|
| **What it does** | NB regression on single-cell counts with guide-level replication and permutation calibration |
| **Language** | R |
| **Paper** | Barry et al. 2023, *Genome Biology* |

**Pros**
- Rigorous type I error control: calibrated using NTC guides as empirical null
- Works with single-cell counts directly (no pseudobulking required)
- Handles low cell counts per guide better than pseudobulk
- Specifically designed for regulatory element screens (CRE-to-gene mapping)
- Well-validated on both simulated and real data

**Cons**
- R only; no Python interface (as of 2025)
- Slower than PyDESeq2 for large gene × perturbation matrices
- Resampling calibration can be computationally expensive

**Differentiator:** Gold standard for type I error control, especially for regulatory element screens where false positives are costly. Use when you need conservative, calibrated p-values.

---

### pertpy — E-distance

| | |
|---|---|
| **What it does** | Non-parametric energy distance between the transcriptome distribution of perturbed vs. NTC cells |

**Pros**
- Model-free; no distributional assumptions
- Sensitive to global shifts, not just individual gene changes
- Works directly on PCA/embedding coordinates — fast
- Useful for ranking perturbations before running full DE

**Cons**
- Not interpretable at the gene level — tells you *that* distributions differ, not *which* genes drive it
- Sensitive to embedding choice (PCA dimensions, normalization)
- No p-value in the standard implementation; bootstrap required for significance

**Differentiator:** Use E-distance as a fast first-pass ranking tool before committing to full pseudobulk DE. Complements, rather than replaces, gene-level DE.

---

### pertpy — Augur

| | |
|---|---|
| **What it does** | Trains a classifier to distinguish perturbed from NTC cells; AUC = effect size per perturbation |

**Pros**
- Cell-type-aware: can be run within each cell cluster separately
- No assumption about which genes are important — the classifier finds them
- AUC score is intuitive: 0.5 = no effect, 1.0 = perfect separation from NTC

**Cons**
- Black box: harder to interpret which genes drive the score
- Slower than E-distance (trains a classifier per perturbation)
- Less useful when many perturbations have small, distributed effects

**Differentiator:** Best for multi-cell-type experiments where you want to know which perturbations are most effective *within a specific cell type*.

| | PyDESeq2 | SCEPTRE | E-distance | Augur |
|---|---|---|---|---|
| Output | Gene-level LFC + p-value | Gene-level LFC + p-value | Per-perturbation scalar | Per-perturbation AUC |
| Requires replicates | Yes | No (guide-level) | No | No |
| Type I error control | Good | Excellent | N/A | N/A |
| Speed | Slow | Medium | Fast | Medium |
| Language | Python | R | Python | Python |

---

## 6. Escape Cell Detection

### pertpy — Mixscape

| | |
|---|---|
| **What it does** | Computes perturbation signature per cell, fits Gaussian mixture model to classify KO vs. NP (not perturbed) |
| **Paper** | Papalexi et al. 2021, *Nature Genetics* |

**Pros**
- Directly removes the major source of dilution in Perturb-seq DE analysis
- Unsupervised: no labels required; mixture model adapts to each perturbation
- `perturbation_signature()` function is useful independently of the classification step
- LPC (perturbation UMAP) embedding is a powerful visualization tool

**Cons**
- Assumes a bimodal distribution of perturbation response — fails for weak perturbations where the two modes are not well separated
- Validated primarily on CRISPRi data; less tested for CRISPRa or CRISPRko
- Can over-aggressively exclude cells if the perturbation effect is subtle
- API has evolved; code from the original paper needs updating for current pertpy

**Differentiator:** Essential preprocessing step before DE analysis, especially for CRISPRi screens where escape rates of 20–50% are common. Removing escape cells substantially increases LFC estimates.

---

## 7. Differential Expression

*(See Section 5 for PyDESeq2 and SCEPTRE full entries.)*

### Wilcoxon rank-sum / t-test on single cells

**Do not use for Perturb-seq DE.** These tests treat cells as independent observations, inflating sample size by 100–500× and producing catastrophically low p-values that are statistically invalid. Use pseudobulk (PyDESeq2) or guide-level models (SCEPTRE) instead.

---

## 8. Transcription Factor Activity

### decoupler-py

| | |
|---|---|
| **What it does** | Infers transcription factor and pathway activity from gene expression using prior networks (CollecTRI, PROGENy, DoRothEA) |
| **Language** | Python |
| **License** | GPL-3 |

**Pros**
- Multiple statistical methods: ULM, VIPER, MLM, GSVA — compare and choose
- Access to curated prior networks: CollecTRI (TF), PROGENy (pathway), DoRothEA (TF with confidence levels)
- Works on any expression matrix: single-cell, pseudobulk, or DE LFC matrix
- More sensitive than DE of the TF itself: detects activity changes even when the TF mRNA is unchanged

**Cons**
- Dependent on prior network quality — if a TF's targets are poorly annotated, the score is unreliable
- ULM assumes linear relationship between TF activity and target expression
- CollecTRI human coverage is good; mouse and other organisms are less complete

**Differentiator:** Run decoupler TF activity inference on your DE LFC matrix (perturbation × gene) as a first-pass pathway interpretation step — it costs almost nothing and reveals the regulatory landscape immediately.

---

### gseapy

| | |
|---|---|
| **What it does** | Enrichment analysis (ORA and GSEA) against MSigDB, GO, KEGG, Reactome gene sets |
| **Language** | Python |

**Pros**
- Broad gene set coverage: 50+ databases via Enrichr API
- Both ORA (over-representation) and GSEA (ranked gene list) modes
- Pure Python; no R dependency
- Output is a standard DataFrame — easy to filter and plot

**Cons**
- ORA ignores gene ranking (all significant genes treated equally)
- GSEA requires a ranked gene list — must define a meaningful ranking metric (e.g., LFC)
- Internet connection required to access Enrichr API; use local databases for offline work

**Differentiator:** Use gseapy for quick pathway interpretation of per-perturbation DEG lists. Pair with decoupler for deeper mechanistic insight.

| | decoupler | gseapy |
|---|---|---|
| Input | Expression matrix or LFC matrix | Gene list or ranked gene list |
| Prior knowledge | TF-target / pathway networks | Gene sets (GO, MSigDB, KEGG) |
| Output | Per-sample TF/pathway activity scores | Enriched gene sets with p-values |
| Best for | Regulatory interpretation | Biological pathway annotation |

---

## 9. Perturbation Prediction & GRN

### GEARS

| | |
|---|---|
| **What it does** | Graph neural network using gene-gene interaction graphs to predict transcriptional response to unseen perturbations |
| **Language** | Python (PyTorch + PyTorch Geometric) |
| **Paper** | Roohani et al. 2023, *Nature Biotechnology* |

**Pros**
- Incorporates prior gene-gene interaction knowledge (GO, STRING) as an inductive bias
- Designed for both single and double perturbation prediction
- Tested on the Norman 2019 dataset (canonical benchmark)
- Can generalize to unseen gene combinations by graph message passing

**Cons**
- Requires GPU for practical training times (CPU is prohibitively slow)
- Complex installation: PyTorch + PyG version compatibility is finicky
- **Kernfeld et al. 2025 (Nature Methods)**: does not consistently outperform a simple mean predictor on held-out single perturbations
- Performance advantage is mostly on double perturbations in specific datasets
- Overfits when training data is small

**Differentiator:** Use GEARS when predicting **double perturbation** responses where no experimental data exists — this is where the graph prior provides genuine value. For single perturbation prediction, always compare against the mean predictor first.

---

### scGPT

| | |
|---|---|
| **What it does** | Foundation model for single-cell biology; perturbation prediction is one of several tasks |
| **Language** | Python (PyTorch) |
| **Paper** | Cui et al. 2024, *Nature Methods* |

**Pros**
- Pre-trained on millions of cells; transfer learning requires much less data than training from scratch
- Handles multiple tasks: cell type annotation, perturbation prediction, gene network inference
- Emerging standard for foundation model approaches in single-cell

**Cons**
- Very large model; requires significant GPU memory (≥24 GB for fine-tuning)
- The perturbation prediction task is one of the less-validated applications
- Same benchmarking concern as GEARS: Kernfeld et al. 2025 shows it does not consistently beat baselines
- Pre-training data composition affects transferability; K562 may not be well-represented

**Differentiator:** Most useful when you have downstream tasks beyond perturbation prediction (cell type annotation, trajectory, etc.). Avoid as a perturbation-only tool until benchmarking matures.

---

### scLAMBDA

| | |
|---|---|
| **What it does** | Emerging deep learning model for perturbation response prediction; explicitly handles dosage effects |
| **Language** | Python (PyTorch) |
| **Status** | Pre-print / early release (2024–2025) |

**Pros**
- Models graded (dose-dependent) perturbation effects — useful for CRISPRa titration data
- Incorporates gene ontology and regulatory network structure

**Cons**
- Very new; limited independent validation
- API and model architecture still evolving
- Smaller community than GEARS or scGPT

**Differentiator:** Watch this space, especially for CRISPRa dose-response modeling. Not yet production-ready.

---

### Linear / mean predictor baseline

| | |
|---|---|
| **What it does** | Predict each gene's response as the mean of all training perturbation responses (mean predictor) or via ridge regression on gene features |

**Pros**
- Trivially fast; no GPU required
- Fully interpretable
- Surprisingly competitive: Kernfeld et al. 2025 showed it matches or beats GEARS on most datasets for single perturbation prediction

**Cons**
- No generalization to novel gene biology; extrapolation is naive
- Cannot model non-linear interactions between perturbations

**Differentiator:** **Always run this first.** It is the mandatory baseline against which any deep learning model must be compared. If your model doesn't beat it, the model provides no practical value.

| | GEARS | scGPT | Linear baseline |
|---|---|---|---|
| Architecture | GNN | Transformer | Ridge regression |
| GPU required | Yes | Yes | No |
| Single pert. prediction | Marginal vs. baseline | Marginal vs. baseline | Baseline |
| Double pert. prediction | Better than baseline | Unclear | Weaker |
| Training data needed | Moderate | Pretrained | Minimal |
| Maturity | High | Medium | N/A |

---

## 10. Databases & Data Resources

### scPerturb

| | |
|---|---|
| **What it does** | 44 harmonized, quality-controlled Perturb-seq and related perturbation datasets in a unified format |
| **URL** | https://scperturb.org |
| **Paper** | Peidli et al. 2024, *Nature Methods* |

**Pros**
- Uniform preprocessing: same pipeline applied to all datasets — directly comparable
- Covers diverse cell types, modalities, and perturbation types
- Available via Zenodo and pertpy integration
- Tutorial notebooks provided

**Cons**
- Preprocessing choices are fixed — you cannot adjust them without reprocessing from raw
- Coverage biased toward published datasets in major journals

**Differentiator:** Best resource for **methods development and benchmarking** — when you want to test a new algorithm across many datasets with consistent preprocessing.

---

### PerturBase

| | |
|---|---|
| **What it does** | Database of 122 Perturb-seq/perturbation datasets with built-in DEG analysis (5 DE methods) and visualization |
| **URL** | http://www.perturbase.cn |
| **Paper** | Wang et al. 2025, *Nucleic Acids Research* |

**Pros**
- Largest collection: 122 datasets, ~5 million cells, 24,254 genetic + 230 chemical perturbations
- Built-in DE analysis — no local installation needed for basic exploration
- Covers ATAC-seq perturbation data as well as RNA
- No registration required

**Cons**
- Web-only interface; bulk download requires contacting authors
- DE analysis uses fixed parameters — cannot customize for your specific question
- Primarily a browsing tool, not an analysis platform

**Differentiator:** Best for **quickly finding and exploring** a perturbation dataset for a gene of interest. Use before designing your own screen to check if the perturbation has already been profiled.

---

### Virtual Cells Platform (CZI)

| | |
|---|---|
| **What it does** | CZ Science pre-processed Perturb-seq datasets with browser-based exploration and download |
| **URL** | https://virtualcellmodels.cziscience.com |

**Pros**
- Pre-processed, ready-to-load h5ad files for major datasets (Norman, Replogle, Adamson)
- Browser-based DE analysis without local setup
- Backed by CZ Science; likely stable long-term hosting

**Cons**
- Small number of datasets (focused on flagship studies)
- Fixed preprocessing; cannot change normalization or QC thresholds

**Differentiator:** The easiest way to download the Norman 2019 and Replogle 2022 datasets used in this course. Start here for quick access to canonical benchmark data.

---

### Addgene Perturb-seq libraries

| | |
|---|---|
| **What it does** | Repository of validated, cloned guide RNA libraries for CRISPRi, CRISPRa, and CRISPRko screens |
| **URL** | https://www.addgene.org |

**Key libraries:**

| Library | Modality | Target | Addgene ID |
|---|---|---|---|
| CRISPRi v2 (Horlbeck) | CRISPRi | ~1,700 human TFs | 83969 |
| CRISPRa v2 (Horlbeck) | CRISPRa | ~1,700 human TFs | 83978 |
| Brunello (Doench) | CRISPRko | ~19,000 human genes | 73178 |
| Brie (Doench) | CRISPRko | ~19,000 mouse genes | 73632 |

**Differentiator:** Use CRISPRi v2 / CRISPRa v2 for any Perturb-seq experiment targeting transcription factors in human cells — these libraries have been validated in K562 and provide 5 guides per gene, matching the Weissman lab standard used in the Norman and Replogle datasets.

---

## Quick-Reference: Tool by Pipeline Stage

| Stage | Recommended tool | Alternative |
|---|---|---|
| FASTQ → count matrix (10x) | Cell Ranger | STARsolo |
| FASTQ → count matrix (non-10x) | kb-python | STARsolo |
| Guide assignment | pertpy GuideRNAAssignment | crispat |
| QC | scanpy `calculate_qc_metrics` | — |
| Multi-modal storage | muon MuData | AnnData + concat |
| Normalization / dim. reduction | scanpy | scVI (batch correction) |
| Escape cell removal | pertpy Mixscape | — |
| Pseudobulk DE | PyDESeq2 | edgeR (R) |
| Rigorous DE (low cell count) | SCEPTRE (R) | — |
| Perturbation ranking | pertpy E-distance | Augur |
| TF activity inference | decoupler-py (CollecTRI + ULM) | — |
| Pathway enrichment | gseapy (Enrichr/MSigDB) | decoupler (PROGENy) |
| Genetic interactions | Norman additive model | pertpy GeneticInteraction |
| Perturbation prediction | Linear baseline → GEARS | scGPT |
| Dataset exploration | PerturBase | scPerturb |
| Canonical benchmark data | Virtual Cells Platform | GEO / Figshare |
