# 01 — Experimental Design

Before running a Perturb-seq experiment, several design decisions determine whether your analysis will be interpretable. This document covers the key parameters you control.

---

## 1. Perturbation modality

Choose based on your biological question:

| Question | Modality |
|---|---|
| Loss-of-function across a gene set | CRISPRi (preferred) or CRISPRko |
| Overexpression / gain-of-function | CRISPRa |
| Both up and down in one experiment | Mix CRISPRa + CRISPRi libraries (complex) |
| Non-coding elements (enhancers) | CRISPRi or CRISPRko targeting the element |

The Cas9/dCas9 must be stably integrated in your cell line before transducing the guide library. Many labs use inducible dCas9 systems (e.g., tet-on) to avoid constitutive perturbation during the expansion phase.

---

## 2. Multiplicity of infection (MOI)

MOI describes the average number of viral particles per cell. This is the single most important design parameter for data quality.

**Target MOI ≈ 0.3** (for 10x Genomics feature barcoding).

At MOI 0.3, the fraction of cells receiving exactly one guide follows a Poisson distribution:
- P(0 guides) ≈ 74% — these cells are excluded from analysis
- P(1 guide) ≈ 22% — usable cells
- P(≥2 guides) ≈ 4% — multiplets; excluded or treated as double perturbations

At higher MOI, multiplet rates increase dramatically. Multiplets are problematic because you cannot unambiguously assign the observed transcriptome to one perturbation. Even a 10% multiplet rate can introduce substantial noise in downstream DE analysis.

**Practical consequence:** You need to sequence more cells than you think. If only ~22% of transduced cells carry exactly one guide, and you need 100 cells per perturbation for 200 perturbations, you need to capture ~200,000 cells × (1/0.22) ≈ 909,000 cells from the 10x chip — which is not realistic. In practice, you select cells by antibiotic selection (e.g., puromycin resistance on the guide vector) to enrich transduced cells before capture, increasing the fraction of guide-positive cells to 60–80%.

---

## 3. Cells per perturbation

Statistical power for detecting perturbation effects scales with the number of cells per perturbation.

| Cells/perturbation | What you can detect |
|---|---|
| < 20 | Essentially nothing reliable; flag these as underpowered |
| 20–50 | Strong perturbations only (e.g., essential genes in proliferating cells) |
| 50–200 | Most strong and moderate perturbations |
| 200–500 | Subtle effects, cell-state heterogeneity within perturbation |
| > 500 | Rare subpopulations, interaction effects |

For pseudobulk DE (the statistically correct approach; see notebook 08), you also need **replication structure** — multiple independent samples/replicates contributing cells to each perturbation. Without replicates, pseudobulk DE collapses to a one-sample test with inflated type I error.

**Practical recommendation:** Aim for ≥100 cells per perturbation across ≥2 independent transductions.

---

## 4. Guides per gene

Use multiple guides targeting each gene to:
1. Distinguish on-target effects from guide-specific off-targets
2. Increase confidence: if two independent guides produce the same transcriptional signature, it is almost certainly on-target
3. Allow guide concordance as an internal QC metric (notebook 06)

**Standard:** 3–5 guides per gene for CRISPRi/a; 2–3 for CRISPRko (given lower off-target rates with modern guide scoring).

Guides targeting the same gene but producing discordant transcriptional signatures indicate guide-specific off-targets or highly variable editing efficiency — those guides should be excluded from pooled analysis.

---

## 5. Controls

Every Perturb-seq library must include control guides. There are several types:

| Control type | Description | Use |
|---|---|---|
| Non-targeting control (NTC) | Sequences with no match in the genome | Baseline transcriptome; reference for DE analysis; QC for guide assignment |
| Safe-harbor control | Guides targeting a known neutral locus (e.g., *AAVS1*, *ROSA26*) | Controls for guide delivery without perturbation |
| Scrambled guide | Random sequence; may have off-targets | Less preferred; NTC is cleaner |
| Positive control | Guides targeting well-characterized genes (e.g., *MYC*, essential ribosomal genes) | Verify the screen is working; benchmark DE pipeline |

**Rules of thumb:**
- Include ≥10 distinct NTC guides so their cells serve as a high-powered baseline
- NTCs should comprise 10–20% of cells in your experiment
- Positive controls help you calibrate the signal-to-noise before drawing conclusions about unknown perturbations

---

## 6. Sequencing depth

Two libraries require two sequencing runs (or a mixed run), each at different depths:

| Library | Recommended depth | Reason |
|---|---|---|
| GEX | ~20,000–25,000 reads/cell | Saturate the expressed transcriptome |
| Guide capture | ~3,000–5,000 reads/cell | UMI counts are low; guide barcode is short |

Under-sequencing the guide capture library leads to failed guide assignments (notebook 06). Over-sequencing GEX costs money without proportional gain in detected genes beyond ~25k reads/cell (check saturation curves in Cell Ranger QC report).

For 10x 3' v3 chemistry, 10,000 cells at 25k reads/cell = 250M GEX reads + 50M guide capture reads.

---

## 7. Single vs. dual perturbation design

**Single perturbation:** One guide per cell; straightforward DE analysis; used in most genome-scale screens.

**Dual perturbation:** Two guides per cell; used to map **genetic interactions** (epistasis). Examples:
- Norman et al. 2019 (this course): CRISPRa, K562, systematic pairwise perturbations of 105 genes
- Horlbeck et al. 2018: CRISPRi dual-guide library for interaction mapping

For dual-guide screens:
- MOI must be calibrated so most cells get exactly one guide from library A and one from library B
- Guide capture library must distinguish guide identity from two separate pools
- Analysis requires an additive null model (notebook 11): the expected double-perturbation effect is the sum of the two single-perturbation effects; deviation is the interaction

**Practical note:** The Norman 2019 dataset used in this course is a mixed single + dual guide experiment. Single-guide cells are used to learn single-perturbation effects; dual-guide cells are used to measure interactions.

---

## 8. Cell line selection

Perturb-seq works in any Cas9/dCas9-expressing cell line or primary cells. Common choices:

| Cell type | Cas9 delivery | Typical use |
|---|---|---|
| K562 (CML, suspension) | Lentiviral integration | Canonical benchmark; Norman, Replogle datasets |
| RPE1 (retinal, adherent) | Lentiviral integration | Non-cancerous diploid; Replogle dataset |
| HEK293T | Lentiviral | Easy to transduce; many protocols |
| Primary T cells | Electroporation (RNP) | Immune screen; harder protocol |
| iPSC-derived neurons | Lentiviral + differentiation | Disease modeling; long timelines |

Cell line choice affects interpretation: K562 is a cancer cell line with aneuploidy; effects of tumor-suppressor perturbation may differ from diploid cells.

---

**Next:** [02_guide_design.md](02_guide_design.md)
