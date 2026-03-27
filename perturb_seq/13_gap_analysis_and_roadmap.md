# 13 — Gap Analysis and Course Roadmap

This document is an honest assessment of what the current perturb-101 curriculum covers well, where the gaps are, and what stretch topics would meaningfully extend it. It is written for students who have completed notebooks 00–12 and want to know what is missing or what to learn next.

---

## What Is Well Covered

The core pipeline from experimental concept through GRN inference is complete:

| Stage | Coverage |
|---|---|
| Conceptual foundations (modalities, MOI, guide design) | Deep |
| Raw data processing concepts (Cell Ranger, STARsolo, kb-python) | Moderate |
| Data acquisition and AnnData structure | Deep |
| Quality control (scRNA-seq metrics + per-perturbation visualization) | Deep |
| Guide assignment (pertpy probabilistic model, concordance) | Deep |
| NTC-anchored normalization and dimensionality reduction | Deep |
| Pseudobulk differential expression (PyDESeq2) | Deep |
| Escape cell detection (Mixscape) | Deep |
| Perturbation effect visualization (volcano, heatmap, perturbation UMAP) | Moderate |
| Genetic interaction analysis (additive model, interaction scores) | Deep |
| TF activity inference (decoupler-py, CollecTRI) | Moderate |
| Perturbation prediction with critical benchmarking (Kernfeld 2025) | Moderate |
| Tool comparison reference | Deep |

---

## Implementation Gaps

These are topics that are covered conceptually but not fully demonstrated in runnable code.

### 1. Guide UMI matrix is simulated, not real

**Where:** `05_quality_control.ipynb` (guide UMI histogram), `06_guide_assignment.ipynb` (assignment demo)

**Problem:** The Norman 2019 dataset we use comes pre-assigned — it does not include the raw Cell Ranger guide capture matrix. The UMI histogram in notebook 05 uses hardcoded Poisson parameters to illustrate the bimodal distribution, and notebook 06 runs assignment on a synthetic AnnData.

**What's missing:** A worked example starting from an actual `crispr_filtered_feature_bc_matrix/` directory from Cell Ranger, merging it with the GEX matrix by cell barcode, and running pertpy assignment from scratch.

**Workaround:** The 10x Genomics website hosts publicly downloadable Cell Ranger outputs from their CRISPR Feature Barcoding demo datasets. These can be used as drop-in replacements.

---

### 2. Batch correction not covered

**Where:** After `07_normalization_dimred.ipynb`

**Problem:** Real Perturb-seq screens span multiple transduction events, sequencing runs, or cell passages. Batch effects are the norm, not the exception. The current curriculum assumes a single clean batch, which is appropriate for the Norman 2019 download but unrealistic for new experiments.

**What's missing:** A notebook demonstrating:
- Diagnosing batch effects: UMAP colored by batch, kBET statistic
- Harmony correction (`harmonypy`) for fast batch correction on PCA coordinates
- scVI for joint embedding and batch correction (heavier but more principled)
- Critical question: should you correct batch effects before or after Mixscape?

**Planned notebook:** `13_batch_correction.ipynb`

---

### 3. SCEPTRE not runnable from Python

**Where:** `09_mixscape.ipynb` (conceptual section only)

**Problem:** SCEPTRE is described as the gold standard for rigorous type I error control in CRISPR screen DE analysis, but it is implemented in R with no Python interface. The notebook points to the R package documentation without bridging the gap.

**What's missing:** Either (a) an R notebook demonstrating SCEPTRE on the Norman dataset or (b) guidance on running SCEPTRE from Python via `rpy2` or as a subprocess, then reading the results back into pandas.

**Workaround:** Export the Norman assigned AnnData to an R-readable format (`write_csv` for obs metadata + `write_csv` for a subset of the count matrix) and follow the SCEPTRE vignette at https://katsevich-lab.github.io/sceptre/.

---

### 4. GEARS training deferred to GPU

**Where:** `12_grn_and_prediction.ipynb`

**Problem:** The GEARS section checks for GPU availability and gracefully falls back to the linear baseline. In practice, most students will not have a GPU available and will never see a trained GEARS model.

**What's missing:** A pre-trained GEARS checkpoint for the Norman 2019 dataset that students can load and use for inference without training. The GEARS GitHub repository hosts pre-trained checkpoints — pointing to these would close this gap.

**Workaround:** `model.load_pretrained('norman')` loads the published Norman 2019 checkpoint from the GEARS model zoo. Add this as an alternative to training from scratch.

---

### 5. Interaction significance not tested

**Where:** `11_genetic_interactions.ipynb`

**Problem:** Interaction scores are computed and ranked, but no statistical significance is assigned. The student cannot distinguish a real interaction from one that arises by chance from low cell counts.

**What's missing:** A permutation test: shuffle the NTC assignment labels many times, recompute interaction scores for each shuffle, and use the resulting null distribution to compute empirical p-values. This is standard in the genetic interaction literature.

**Computational note:** Permutation testing across all gene pairs is slow but embarrassingly parallelizable. A bootstrap approximation on the top 20 pairs is tractable on a laptop.

---

### 6. Augur not demonstrated

**Where:** `tools_reference.md` (description only)

**Problem:** Augur (`pertpy`) scores perturbations by how well a classifier can distinguish perturbed from NTC cells within a cell type. It is particularly useful for heterogeneous cell populations where different perturbations affect different subpopulations. The current content describes Augur but never runs it.

**What's missing:** A code block in `08_perturbation_effects.ipynb` or a new visualization notebook section showing:
```python
ag = pt.tl.Augur()
augur_results = ag.load(adata, label_col=PERT_COL, cell_type_col="leiden")
```

---

### 7. Pathway enrichment not scaled

**Where:** `10_visualization.ipynb`

**Problem:** Pathway enrichment (ORA via gseapy) is demonstrated for only the single top perturbation. In a real screen you want enrichment for every perturbation, aggregated into a pathway × perturbation activity matrix.

**What's missing:** A loop over all perturbations, collecting top enriched pathways, building a pathway × perturbation matrix, and clustering it. This is the most interpretable output of a Perturb-seq screen for a biology audience.

---

## Conceptual Gaps

These topics are not addressed anywhere in the current curriculum.

### Statistical power and sample size

The most common question from anyone designing a Perturb-seq experiment is: *How many cells do I need per perturbation?* The current curriculum provides rules of thumb (≥100 cells/perturbation) without derivation or quantitative justification.

A rigorous answer requires:
- Simulating pseudobulk count matrices with known effect sizes
- Running PyDESeq2 on these simulations
- Computing power (fraction of simulated DEGs recovered) as a function of n_cells, effect_size, and dispersion

This is critical knowledge for experimental design and is currently absent.

**Planned notebook:** `15_power_analysis.ipynb`

---

### Multiple testing across perturbations

When testing 2,000 perturbations × 20,000 genes, you are performing 40 million statistical tests. The curriculum covers within-perturbation FDR correction (the standard Benjamini-Hochberg `padj` from DESeq2), but not cross-perturbation multiple testing.

The relevant questions:
- Should you correct for the number of perturbations tested?
- How does this change the threshold for calling a perturbation "significant"?
- What is the right FDR framework when the goal is to rank perturbations (discovery) vs. claim specific interactions (confirmation)?

---

### Negative control calibration

NTC guides are not just for DE reference — they are also a calibration tool. Running DE for each NTC guide vs. the other NTC guides gives you an empirical null distribution of effect sizes and p-values under no perturbation. If your NTC-vs-NTC DE results show many "significant" genes, your DE pipeline is miscalibrated (underdispersed model, batch confounding, etc.).

This is the foundation of SCEPTRE's calibration check and should be demonstrated as a standalone QC step.

---

### Cell-to-cell variability within a perturbation

The curriculum treats all cells with the same guide as a homogeneous group (after Mixscape filtering). But perturbed cells can be heterogeneous: some may be in different cell cycle stages, some may have compensatory adaptations, and some may represent distinct cell states induced by the perturbation.

Tools for exploring this:
- Subcluster cells within each perturbation; find dominant subclusters
- Test whether subcluster identity predicts treatment response
- Compare intra-perturbation variance to inter-perturbation variance

---

### Comparison to bulk CRISPR screens

The curriculum assumes you have already decided to use Perturb-seq. But many biology questions are better addressed by bulk CRISPR dropout screens (e.g., viability screens, FACS-based reporter screens) which are 10–100× cheaper per perturbation.

A decision framework comparing:
- Perturb-seq: transcriptome-wide, mechanistic, expensive, 50–500 cells/perturbation
- MAGeCK dropout screen: one phenotype, cheap, millions of cells/perturbation
- CROPseq vs. arrayed: pooled vs. individual well testing

---

### Cost and scale tradeoffs

Practical considerations not covered:
- **10x Flex** (fixation-based): enables sample multiplexing for multi-condition Perturb-seq; lower per-cell cost
- **Parse Biosciences**: split-pool barcoding alternative to 10x; scales to millions of cells per run; no microfluidics
- **10x CRISPR Feature Barcoding v2**: updated chemistry with higher guide assignment rates
- Rough cost estimates: $1–5 per cell sequenced for GEX; guide capture adds ~$0.30/cell

---

## Stretch Topics: Proposed New Content

These are substantive extensions requiring new notebooks or markdown documents. Listed in recommended implementation order.

---

### `13_batch_correction.ipynb` (Priority 1)

**What:** Batch effect diagnosis and correction for multi-batch Perturb-seq screens.

**Dataset:** Norman 2019, artificially split into two "batches" by randomly assigning cells (no real batch effect — useful for demonstrating that correction should not be applied blindly). Alternatively, the Adamson et al. 2016 dataset has genuine replicate structure.

**Tools:** `harmonypy`, `scvi-tools` (optional, heavier)

**Key questions answered:**
1. How do I know if I have a batch effect? (kBET, UMAP colored by batch, per-batch QC)
2. Harmony: correct on PCA coordinates; re-run UMAP; compare before/after
3. Does batch correction affect perturbation effect estimates? (it should not if done correctly)
4. When should you correct and when should you not? (overcorrection can remove real biology)

---

### `14_dosage_effects.ipynb` (Priority 1)

**What:** Modeling CRISPRi/a guide efficiency as a continuous variable; dose-response analysis.

**Dataset:** Jost et al. 2020 (*Nature Biotechnology*) — GSE132080. 25 target genes, each targeted by 5 guides with 0 to 3 mismatches from the perfect-match guide. The mismatch number controls knockdown efficiency, providing a natural dose series.

**Tools:** `scanpy`, `scipy.stats`, `statsmodels`, `seaborn`

**Key questions answered:**
1. How does transcriptional response scale with guide efficiency?
2. Can you estimate the "full knockdown" effect by extrapolating from partial knockdowns?
3. Which genes are most sensitive to partial perturbation (steep dose-response)?
4. How does this inform guide selection for future experiments?

**Why this matters:** Nearly all CRISPRi experiments have guides with varying efficiency (by design or by chance). Treating them as equivalent (as current notebooks do) dilutes signal. Modeling efficiency as a covariate recovers signal and improves LFC estimates.

---

### `15_power_analysis.ipynb` (Priority 1)

**What:** Simulation-based statistical power analysis for Perturb-seq experimental design.

**Dataset:** No download required. Simulated data from realistic negative binomial parameters estimated from the Norman 2019 dataset.

**Tools:** `numpy`, `scipy.stats`, `pydeseq2`, `matplotlib`

**Simulation approach:**
1. Estimate NB dispersion parameters from NTC cells in Norman 2019
2. Simulate pseudobulk count matrices with known DE genes (effect size = LFC from real perturbations)
3. Run PyDESeq2 on simulated data
4. Repeat 100 times per (n_cells, effect_size) combination
5. Power = fraction of simulated DEGs called significant at padj < 0.05

**Output:** Power curves as a function of cells/perturbation, effect size (LFC), and number of guides (replicates). A practical lookup table for experimental design.

---

### `16_regulatory_element_screens.md` (Priority 2)

**What:** Conceptual overview of CRE-to-gene mapping screens — tiling dCas9-KRAB across enhancers and promoters to identify regulatory elements.

**Key differences from gene-targeting Perturb-seq:**
- Perturbations tile genomic regions (100–500 bp windows), not genes
- Effect sizes are typically smaller (partial repression of one regulatory element)
- Statistical methods must be more conservative (many more tests; smaller effects)
- The question is "which element regulates which gene?" not "what does gene X do?"

**Datasets:** Gasperini et al. 2019 (*Cell*), Fulco et al. 2019 (*Nature Genetics*)

**Tools:** SCEPTRE is the method of choice for CRE screens (rigorous type I error control at small effect sizes). The Activity-by-Contact (ABC) model provides a complementary prior for enhancer-gene assignment.

**Why it matters:** Regulatory element screens are now a major application of Perturb-seq infrastructure, especially for interpreting GWAS hits. The experimental setup is identical to gene-targeting screens but the analysis requires different statistical assumptions.

---

### `17_multimodal_perturbseq.ipynb` (Priority 2)

**What:** ECCITE-seq — simultaneous GEX + surface protein (CITE-seq antibody barcodes) + CRISPR guide capture in the same cells.

**Dataset:** Mimitou et al. 2021 (*Nature Methods*, ECCITE-seq paper) — available on GEO: GSE153056. T cells with 111 CRISPR perturbations + 49 surface proteins measured simultaneously.

**Tools:** `muon`, `scanpy`, `pertpy`

**Key additions over standard Perturb-seq:**
1. Surface protein normalization: CLR (centered log-ratio) normalization for antibody counts
2. MuData: hold GEX + protein + guide in one object
3. Protein-based cell typing: use surface markers to define T cell subsets before analyzing perturbation effects within each subset
4. Cross-modal correlation: which perturbations change surface protein levels vs. transcriptome?

**Why it matters:** Surface protein readout adds an orthogonal measurement that validates transcriptional changes and enables cell type identification independent of GEX clustering. Increasingly common in immune cell screens.

---

### `18_trajectory_and_cellstate.ipynb` (Priority 2)

**What:** Perturb-seq in differentiating systems; analyzing how perturbations alter cell fate trajectories.

**Dataset:** A differentiation-focused Perturb-seq dataset from scPerturb (https://scperturb.org) — e.g., the Morris et al. iPSC-to-cardiomyocyte screen, or the Peidli et al. 2024 collection entry for a differentiation screen.

**Tools:** `cellrank` or `scFates` for trajectory inference; `scanpy` for pseudotime; `pertpy` for perturbation effects within trajectory bins

**Key questions answered:**
1. Does perturbation X accelerate or delay differentiation? (shift in pseudotime distribution)
2. Does perturbation X redirect cells toward an alternative fate? (change in fate probability)
3. Which genes are perturbed most strongly at specific pseudotime points?
4. How to compare perturbation effects in a pseudotime-aware manner (not just global means)

**Why it matters:** The most biologically important Perturb-seq applications — identifying factors that control cell fate — require trajectory-aware analysis. Treating differentiation as static ignores the most interesting dimension.

---

### `19_comparison_to_bulk_screens.md` (Priority 3)

**What:** Decision framework for choosing between Perturb-seq and bulk CRISPR screens.

| Dimension | Perturb-seq | Bulk dropout screen | Arrayed screen |
|---|---|---|---|
| Readout | Transcriptome (20k genes) | One phenotype (growth/reporter) | One phenotype per well |
| Cost per perturbation | High ($5–50) | Very low ($0.001) | Medium ($1–5) |
| Perturbations per run | 50–2,000 | Genome-wide (20,000+) | 96–384 |
| Cells per perturbation | 50–500 | Millions | Thousands |
| Mechanistic insight | High | None | Low |
| Statistical power | Moderate | Very high | High |

**MAGeCK** (for bulk dropout analysis) is briefly introduced as the standard tool, with pointers to its documentation.

---

### `20_emerging_technologies.md` (Priority 3)

**What:** Survey of Perturb-seq extensions not yet mature enough for full tutorials.

**Topics:**
- **In vivo Perturb-seq:** AAV delivery of guide libraries to mouse tissue; brain and liver applications (Swiech et al., Adamson et al. in vivo)
- **snRNA-seq + CRISPR:** Single-nucleus Perturb-seq for post-mitotic cells (neurons) or frozen tissue; higher ambient RNA background
- **Spatial Perturb-seq:** Combining MERFISH or Visium with CRISPR perturbations; very early stage; key challenge is retaining spatial information after perturbation
- **Base and prime editors:** CRISnot-cutting perturbations (A→G, C→T edits); fewer DSBs; better for sensitive cell types; pertpy supports base editor variant effect analysis
- **10x Flex:** Fixation-based scRNA-seq; enables pooling of 4–16 samples with CMO hashtags; cheaper per cell; emerging for multi-condition Perturb-seq
- **Parse Biosciences:** Split-pool barcoding; scales to millions of cells; compatible with CRISPR perturbations via guide capture

---

## Summary: Priority Ranking

| Priority | Content | Type | New data needed? |
|---|---|---|---|
| Fix now | Hash verification in notebook 03 | Code edit | No |
| Fix now | Augur demo in notebook 08 | Code addition | No |
| High | `15_power_analysis.ipynb` | New notebook | No (simulated) |
| High | `13_batch_correction.ipynb` | New notebook | No (split Norman) |
| High | `14_dosage_effects.ipynb` | New notebook | Yes (Jost 2020, GSE132080, small) |
| Medium | `16_regulatory_element_screens.md` | New markdown | No |
| Medium | `17_multimodal_perturbseq.ipynb` | New notebook | Yes (Mimitou 2021, GEO) |
| Medium | `18_trajectory_and_cellstate.ipynb` | New notebook | Yes (scPerturb) |
| Low | `19_comparison_to_bulk_screens.md` | New markdown | No |
| Low | `20_emerging_technologies.md` | New markdown | No |

---

## Recommended Next Steps for Students

If you have finished notebooks 00–12 and want to go deeper:

1. **Immediately actionable:** Read the [pertpy documentation](https://pertpy.readthedocs.io/) — the API evolves quickly and the tutorials there will cover edge cases not shown here.

2. **Most important missing skill:** Statistical power. Run a quick simulation in Python to understand how your planned n_cells/perturbation translates to power for your expected effect size before starting a new experiment.

3. **Most important missing tool:** SCEPTRE for rigorous DE, especially if you are running a regulatory element screen or any screen where false positives are expensive.

4. **Best dataset to explore independently:** The [scPerturb](https://scperturb.org) collection — 44 harmonized datasets, all in the same format as the files used in this course. Pick one with cell types or biology relevant to your work and re-run the pipeline.

5. **Best paper to read next:** Replogle et al. 2022 (*Cell*) supplementary methods — the most detailed description of end-to-end genome-scale Perturb-seq analysis in print.
