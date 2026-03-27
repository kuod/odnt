# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**perturb-101** is a sequential, end-to-end tutorial curriculum for Perturb-seq analysis — combining CRISPR perturbations with single-cell RNA-seq. It is a pure documentation/tutorial project with no software package or build system.

## Environment Setup

```bash
conda env create -f environment.yml
conda activate perturb-101
jupyter lab
```

There are no build, test, or lint commands — each notebook is self-contained and run interactively in JupyterLab.

## Curriculum Structure

13 numbered units (00–12 plus gap analysis at 13), mixing Markdown guides and Jupyter notebooks:

- **00–02, 04** — Markdown conceptual guides (overview, experimental design, guide RNA design, raw processing)
- **03** — Notebook: data download (Norman 2019 CRISPRa, Replogle 2022 CRISPRi datasets via `pooch`)
- **05–12** — Notebooks: QC → guide assignment → normalization → perturbation effects → escape detection → visualization → genetic interactions → GRN inference
- **13** — Gap analysis and roadmap

Notebooks must be run in order; 03 downloads raw data to `data/` (git-ignored). Notebooks 08 and 09 can run in parallel after 07.

## Key Architecture Decisions

**Data containers:** AnnData (`.h5ad`) for single-modality, MuData for multi-modal (GEX + guide matrices). Intermediate objects are written to `data/` between notebooks and not committed.

**Two primary datasets:**
- Norman et al. 2019 (*Science*): CRISPRa, K562, 236 perturbations, ~109k cells
- Replogle et al. 2022 (*Cell*): CRISPRi, K562, ~2,057 genes, ~650k cells

**Core analysis stack:** `scanpy` + `pertpy` (guide assignment, Mixscape, distances) + `pydeseq2` (pseudobulk DE) + `decoupler-py` (TF activity).

**Perturbation modalities covered:** CRISPRko (Cas9 indels), CRISPRi (dCas9-KRAB repression), CRISPRa (dCas9-VP64/SAM activation).

## Known Gaps (from 13_gap_analysis_and_roadmap.md)

- Guide UMI matrix in notebook 06 is simulated, not real Cell Ranger output
- Batch correction not covered (harmonypy/scvi-tools are commented out in environment.yml)
- SCEPTRE (rigorous permutation-based DE) not runnable from Python
- GEARS deep learning deferred (requires GPU); only a linear baseline is demonstrated
- Statistical power analysis absent
- Multiple testing correction across perturbations not addressed

## Supporting Reference Files

- **glossary.md** — Definitions for all domain-specific terms
- **tools_reference.md** — Pros/cons/differentiators for every tool across the pipeline
