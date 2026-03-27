# Glossary

Alphabetical reference for key Perturb-seq terms.

---

**AnnData** (`adata`)
The core data container in the scverse ecosystem (Python package: `anndata`). Stores a cells × genes count matrix (`adata.X`) alongside per-cell metadata (`adata.obs`), per-gene metadata (`adata.var`), embeddings (`adata.obsm`), and unstructured data (`adata.uns`). All notebooks in this course use AnnData as the primary data structure.

**Augur**
A method (implemented in `pertpy`) for scoring how strongly each perturbation changes the transcriptome relative to non-targeting controls, using a machine learning classifier. Useful for ranking perturbations by effect size without running full DE analysis.

**Cell barcode**
A short (16 bp for 10x Genomics v3) oligonucleotide that labels all molecules from one cell during droplet capture. Used to demultiplex sequencing reads to individual cells. Both the GEX and guide capture libraries from the same cell share the same barcode.

**CRISPRa (CRISPR activation)**
CRISPR-based transcriptional activation. Uses a catalytically dead Cas9 (dCas9) fused to transcriptional activators (VPR, SAM system) directed to the TSS of a target gene. Increases transcription without cutting DNA. Reversible.

**CRISPRi (CRISPR interference)**
CRISPR-based transcriptional repression. Uses dCas9 fused to a KRAB repressor domain directed to the TSS. Deposits repressive chromatin marks (H3K9me3) and blocks Pol II. Reversible.

**CRISPRko (CRISPR knockout)**
Canonical CRISPR editing using nuclease-active Cas9 to create double-strand breaks, leading to indels via NHEJ and subsequent loss of protein function. Irreversible; editing efficiency varies by cell.

**CROP-seq**
A variant of Perturb-seq (Datlinger et al. 2017) that captures guide RNA sequences in the 3' UTR of a reporter transcript, enabling co-capture with GEX in the same read. Functionally very similar to standard Perturb-seq; both are now generically called "Perturb-seq."

**dCas9**
Catalytically dead Cas9: point mutations in both nuclease domains (D10A, H840A) abolish DNA cutting while retaining DNA binding. Used for CRISPRi and CRISPRa.

**DEG (differentially expressed gene)**
A gene whose expression level is statistically significantly different between a perturbed condition and a control (NTC) condition. Typically defined by adjusted p-value < 0.05 and |log2 fold change| > some threshold.

**E-distance (energy distance)**
A non-parametric distance metric between two multivariate distributions (Székely & Rizzo 2004). Used in `pertpy` to measure how much the transcriptome distribution of perturbed cells differs from NTC cells. Does not require DE gene selection; sensitive to global distribution shifts.

**Epistasis**
A genetic interaction where the combined effect of two perturbations differs from the sum of their individual effects. In Perturb-seq, epistasis is quantified as: interaction score = observed double-perturbation effect − expected (additive) effect. See notebook 11.

**Escape cells**
Cells that received a guide RNA but show no transcriptional response consistent with the expected perturbation. Caused by failed editing (CRISPRko), incomplete dCas9-KRAB recruitment (CRISPRi), or epigenetic compensation. Identified and removed by Mixscape (notebook 09).

**Feature barcoding**
A 10x Genomics technology that captures non-mRNA molecular labels (guide barcodes, antibody barcodes, CMO hashtags) in the same droplet as the GEX library, linked by the cell barcode. The foundation of Perturb-seq guide capture.

**GEARS**
Graph-based Encoder for Accurate RNA-seq Simulation (Roohani et al. 2023). A deep learning model that uses gene-gene interaction graphs to predict transcriptional responses to unseen perturbations. See notebook 12. Note: Kernfeld et al. 2025 showed that GEARS does not consistently outperform linear baselines.

**Guide assignment**
The step that maps each cell barcode to a guide RNA identity based on guide capture UMI counts. Produces `guide_id`, `gene_target`, and `assignment_confidence` columns in `adata.obs`. See notebook 06.

**Guide concordance**
A QC metric: the correlation between the mean expression profiles of cells assigned to two different guides targeting the same gene. High concordance (r > 0.5) suggests both guides have the expected on-target effect; low concordance suggests one guide has off-target effects.

**Guide purity score**
The fraction of guide capture UMIs in a cell that map to the most-detected guide. A cell with 90% of UMIs from one guide has high purity; a cell with 50/50 split is likely a multiplet. Threshold typically ≥ 0.75.

**HVG (highly variable gene)**
Genes with high variance across cells, selected for dimensionality reduction (PCA) to focus on informative signal. In Perturb-seq, HVGs should be selected from NTC cells only to avoid perturbation-driven genes dominating the embedding (notebook 07).

**Interaction score**
A quantitative measure of genetic interaction: `interaction score = observed_double - expected_double`, where `expected_double = single_A_deviation + single_B_deviation`. Positive scores = synergy; negative scores = buffering/alleviating. See notebook 11.

**lentiGuide / lentiCRISPRv2 / SAM**
Lentiviral vectors for guide RNA delivery:
- **lentiGuide**: guide-only vector; used with cell lines stably expressing Cas9/dCas9
- **lentiCRISPRv2**: all-in-one Cas9 + guide vector; less common for pooled screens (integration of large Cas9 causes packaging issues)
- **SAM (Synergistic Activation Mediator)**: three-component CRISPRa system (dCas9-VP64 + sgRNA2.0 + MPH effector); very high activation efficiency

**LFC (log fold change)**
The log2 ratio of mean expression in perturbed cells vs. NTC cells: `LFC = log2(mean_perturbed / mean_NTC)`. Positive LFC = upregulated; negative LFC = downregulated. Often shrunk toward zero using empirical Bayes shrinkage (apeGLM) to stabilize estimates from low-count genes.

**Mixscape**
A method (Papalexi et al. 2021, implemented in `pertpy`) for classifying cells in a Perturb-seq experiment as: **KO** (successfully perturbed), **NP** (not perturbed / escape), or **NT** (non-targeting control). Uses a Gaussian mixture model on the perturbation signature. See notebook 09.

**MOI (multiplicity of infection)**
The average number of viral particles per cell in a transduction. For Perturb-seq, target MOI ≈ 0.3 to minimize multiplets (cells with two guides). See notebook 01.

**MuData**
A multi-modal AnnData container (Python package: `muon`). Holds multiple modalities (e.g., `rna`, `crispr`, `prot`) with shared cell barcodes and a joint observation space.

**Non-targeting control (NTC)**
Guide RNAs with sequences that have no match in the genome. NTC cells experience guide delivery and Cas9 expression but no gene perturbation. They serve as the reference baseline for all DE comparisons.

**Perturbation signature**
The vector of mean expression differences between cells with a given perturbation and NTC cells. Computed as `mean(perturbed_cells) - mean(NTC_cells)` on the normalized expression matrix. Used by Mixscape as input to the mixture model.

**Pseudobulk**
An analysis approach that aggregates single-cell counts from the same perturbation (and same replicate/sample) into a "pseudo-sample," then applies bulk RNA-seq DE methods (DESeq2, edgeR). Statistically correct for single-cell experiments because it avoids pseudoreplication. See notebook 08.

**PyDESeq2**
A Python reimplementation of the DESeq2 negative binomial model for differential expression. Used in notebook 08 for pseudobulk analysis. Results are concordant with R's DESeq2 for most use cases.

**Safe harbor**
A genomic locus (e.g., *AAVS1* on chr19, *ROSA26* in mouse) that tolerates transgene insertion without disrupting nearby gene expression. Used as a negative control in Perturb-seq: guides targeting the safe harbor induce a cut (CRISPRko) or binding event (CRISPRi/a) with minimal transcriptional consequence.

**SCEPTRE**
A statistical framework (Katsevich et al. 2021; Barry et al. 2023) for CRISPR screen differential expression analysis. Uses a negative binomial model, guide-level replicates as the unit of analysis, and permutation calibration using NTC guides. More conservative than pseudobulk DESeq2; well-calibrated type I error. Available as an R package.

**sgRNA (single guide RNA)**
The guide RNA used with Cas9/dCas9. Consists of the 20-nt spacer fused to the constant scaffold sequence. "Single guide" distinguishes it from the original two-component system (crRNA + tracrRNA). The terms guide RNA and sgRNA are used interchangeably in the Perturb-seq literature.

**TSS (transcription start site)**
The genomic position where RNA polymerase II initiates transcription. CRISPRi and CRISPRa guides should target within ±200 bp of the TSS for maximum effect. TSS coordinates are available from ENCODE, RefSeq, or the FANTOM5 CAGE atlas.

**UMI (unique molecular identifier)**
A short random oligonucleotide (12 bp for 10x v3) ligated to each cDNA molecule during reverse transcription. Used to collapse PCR duplicates: reads with the same cell barcode + UMI + gene assignment represent one original molecule. UMI counts are the raw input to all downstream analysis.
