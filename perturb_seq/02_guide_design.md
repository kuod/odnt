# 02 — Guide RNA Design

Guide RNA design determines perturbation efficiency and specificity. Poor guides produce high escape rates (notebook 09) and confound downstream analysis. This document covers the design principles for CRISPRi/a and CRISPRko, and how guides acquire the barcodes that make Perturb-seq possible.

---

## 1. Guide RNA structure

A single guide RNA (sgRNA) for SpCas9 consists of:

```
5'—[20 nt spacer]—[scaffold (≈80 nt)]—3'
```

- **Spacer**: the 20-nucleotide sequence complementary to your target; this is what you design
- **Scaffold**: the constant sequence that folds into the Cas9-binding hairpin; do not modify
- **PAM**: not part of the guide, but Cas9 requires `NGG` immediately 3' of the target in the genome

For Cas12a (used in some screens): PAM is `TTTV` 5' of the target; spacer is 20–24 nt.

---

## 2. Choosing target positions

Target position requirements differ by modality:

### CRISPRko

- Target any constitutive exon that, when disrupted by indels, causes frameshift/truncation
- Prefer early exons (exon 2–4) and coding sequence, not UTRs
- Avoid the last exon (truncated C-terminal protein may retain partial function)
- Avoid alternatively spliced regions unless targeting a specific isoform

### CRISPRi

- Target the **proximal promoter**: within ±200 bp of the transcription start site (TSS)
- The dCas9-KRAB complex sterically blocks Pol II and deposits repressive H3K9me3
- Multiple guides tiling the TSS window are additive; 3–5 guides/gene is standard
- Avoid targeting the antisense strand of a nearby gene

### CRISPRa

- Target the same TSS window as CRISPRi (−200 to +50 bp relative to TSS)
- For the SAM system (synergistic activation mediator), place 2–3 guides within −200 to −50 bp upstream
- Multiple guides act synergistically for activation; using 2–3 guides per gene is strongly recommended

---

## 3. On-target scoring

After identifying all candidate 20-mers with the right PAM, score them for predicted efficiency.

### Doench 2016 Rule Set 2 (RuleSet2)

The most widely used guide scoring model for SpCas9. Trained on dropout screen data from human and mouse cell lines. Scores guides 0–1 (higher = more active).

- Python: `rs2score` function from the `azimuth` package, or use the CRISPOR web tool
- Recommended threshold: RS2 score ≥ 0.5 for CRISPRko

### CRISPRscan

A complementary scoring algorithm (Moreno-Mateos et al. 2015) based on in vitro cleavage data. Useful as a second scoring layer; pick guides that score well on both.

### CRISPRi/a-specific scoring

Standard cutting efficiency scores are less predictive for CRISPRi/a. Instead:
- Prefer guides within the optimal TSS window (position matters more than sequence score)
- Use experimentally validated guide libraries from the Weissman lab (published CRISPRi v2 and CRISPRa v2 libraries; available on Addgene)
- The Weissman CRISPRi/a libraries have been validated in K562 cells — the same cell line used in this course

---

## 4. Off-target prediction

Every 20-mer spacer can potentially direct Cas9 to similar genomic sequences. Off-target editing introduces confounders in Perturb-seq: a transcriptional response you attribute to the guide's intended target gene may actually be caused by an unintended edit elsewhere.

### CRISPOR (recommended)

Web tool and command-line tool: http://crispor.gi.ucsf.edu

- Input: target gene name or coordinates
- Output: all candidate guides with RS2 scores + off-target prediction (MIT score, CFD score)
- Select guides with CFD off-target score < 0.2 and no predicted off-targets with < 3 mismatches in coding regions

### Cas-OFFinder

Command-line tool for exhaustive off-target site enumeration allowing mismatches and bulges. Use for final validation of your top guide candidates.

### Practical mitigation

The best mitigation is using **multiple independent guides per gene** and requiring concordance. If two guides with different predicted off-targets produce the same transcriptional signature, the common signal is almost certainly on-target. This is why guide concordance analysis (notebook 06) is a core QC step.

---

## 5. Library design considerations

### Number of guides per gene

- CRISPRi/a: 5 guides/gene (the Weissman lab standard)
- CRISPRko: 2–3 guides/gene
- NTC: ≥10 distinct non-targeting sequences (use random sequences with no 4-bp seed match to any human/mouse gene)

### Safe harbors and positive controls

Add guides for:
- Safe harbor (AAVS1 site): 3–5 guides
- Positive control genes with known strong transcriptional effects: 3–5 guides per gene
  - CRISPRi positive control: *MYC* knockdown (strong growth phenotype)
  - CRISPRa positive control: *CDKN1A* (p21) activation

### Library size vs. screen power

More perturbations → fewer cells per perturbation at fixed total cell count. A practical ceiling for a standard 10x Genomics run capturing ~10,000 cells:

```
Total usable cells ≈ 10,000 × (fraction transduced) × (fraction with 1 guide)
                  ≈ 10,000 × 0.7 × 0.85 ≈ 6,000

At 100 cells/perturbation → 60 perturbations maximum
At 50 cells/perturbation  → 120 perturbations maximum
```

Genome-scale screens (like Replogle 2022) require millions of cells and specialized protocols (e.g., 10x Genomics CRISPR Feature Barcoding with very high cell capture).

---

## 6. Feature barcoding: how guides become detectable

Standard scRNA-seq cannot distinguish sgRNA sequences from mRNA (sgRNAs are not polyadenylated and are short). Perturb-seq uses **feature barcoding** to solve this:

```
Guide RNA architecture for feature barcoding:

5'—U6 promoter—[20 nt spacer]—[scaffold]—[unique 8-20 nt guide barcode]—PolyA signal—3'

                                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                            This sequence is captured
                                            by the GEX poly-T primer
                                            in the 10x droplet
```

The unique barcode appended 3' of the scaffold is:
- Short enough to amplify efficiently
- Unique per guide (not per gene — one barcode per guide sequence)
- Distinct from any mRNA sequence

During 10x library preparation, a separate PCR amplification step enriches for guide barcode fragments. This is sequenced as the **CRISPR guide capture library**.

### Alternatives

- **Direct capture (lentiGuide):** Some protocols amplify the U6-sgRNA region directly
- **Feature barcoding v2 (10x):** Updated chemistry for higher guide assignment rates
- **ECCITE-seq / MULTI-seq:** Extend feature barcoding to include antibody barcodes (surface proteins)

---

## 7. Pre-experiment quality check: amplicon sequencing of the guide library

Before transducing cells, verify your synthesized guide library has:
1. **Even representation**: every guide is present at approximately equal frequency (< 5× variation across guides)
2. **No dropped guides**: all guides in the design are present
3. **No unexpected sequences**: synthesis errors or contamination

Protocol: PCR amplify the guide spacer region from the plasmid pool, sequence by Illumina MiSeq or NGS, count reads per guide. A Lorenz curve of guide representation should show minimal inequality; a Gini coefficient < 0.2 is acceptable.

---

## 8. Validated library resources

You do not have to design guides from scratch. Widely used, validated Perturb-seq libraries:

| Library | Modality | Species | Genes | Source |
|---|---|---|---|---|
| CRISPRi v2 (Horlbeck et al. 2016) | CRISPRi | Human | ~1,700 TFs | Addgene #83969 |
| CRISPRa v2 (Horlbeck et al. 2016) | CRISPRa | Human | ~1,700 TFs | Addgene #83978 |
| Genome-wide CRISPRi (Weissman lab) | CRISPRi | Human | ~18,000 | Addgene pooled libraries |
| Brunello (Doench et al. 2016) | CRISPRko | Human | ~19,000 | Addgene #73178 |
| Brie (Doench et al. 2016) | CRISPRko | Mouse | ~19,000 | Addgene #73632 |

---

**Next:** [03_data_acquisition.ipynb](03_data_acquisition.ipynb)
