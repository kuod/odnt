# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**perturb-101** is a tutorial curriculum repository with three series:
1. **`perturb_seq/`** — sequential end-to-end Perturb-seq analysis (CRISPR + scRNA-seq)
2. **`theis_ecosystem/`** — hands-on refresher for the full scverse / Theis lab tool stack
3. **`spatial_transcriptomics/`** — spatial omics from Visium to sub-cellular resolution

Pure documentation/tutorial project; no software package or build system.

## Environment Setup

```bash
conda env create -f environment.yml
conda activate perturb-101
jupyter lab
```

There are no build, test, or lint commands — each notebook is self-contained and run interactively in JupyterLab.

## Repository Structure

```
perturb-101/
├── perturb_seq/          # Perturb-seq curriculum (00–13)
│   ├── 00–02, 04         # Markdown conceptual guides
│   ├── 03                # Data download (Norman 2019, Replogle 2022 via pooch)
│   ├── 05–12             # QC → guide assignment → normalization → DE → Mixscape
│   │                     #   → visualization → genetic interactions → GRN
│   ├── 13                # Differential abundance (Milo)
│   ├── 14                # Batch correction (Harmony + scVI)
│   └── 13_gap            # Gap analysis and roadmap
│   ├── glossary.md
│   └── tools_reference.md
├── theis_ecosystem/      # scverse ecosystem refresher (T00–T08)
│   ├── T00               # AnnData + Scanpy + doublet detection + ambient RNA
│   ├── T01               # scvi-tools (deep generative models, integration)
│   ├── T02               # Squidpy (spatial transcriptomics)
│   ├── T03               # CellRank 2 (trajectory inference)
│   ├── T04               # Muon + MOFA+ (multi-modal)
│   ├── T05               # LIANA (cell-cell communication, scverse core 2025)
│   ├── T06               # cell2location (spatial deconvolution)
│   ├── T07               # Moscot (optimal transport trajectories + spatial)
│   └── T08               # scCODA (compositional cell type analysis)
├── spatial_transcriptomics/  # Spatial omics series (SP00–SP04)
│   ├── SP00              # SpatialData container (Zarr, coordinate systems)
│   ├── SP01              # Visium in depth (spatial QC, regional DE, H&E features)
│   ├── SP02              # Spatial domain detection (BANKSY, GraphST)
│   ├── SP03              # Sub-cellular resolution (seqFISH+, Xenium data model)
│   └── SP04              # Multi-sample integration (Moscot alignment, Harmony)
├── data/                 # git-ignored; created at runtime
├── environment.yml
└── README.md
```

Perturb-seq notebooks must run in order; 03 downloads raw data to `../data/` (git-ignored). Notebooks 08 and 09 can run in parallel after 07. All intermediate `.h5ad` files land in `data/` at repo root.

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

## Supporting Reference Files (in `perturb_seq/`)

- **glossary.md** — Definitions for all domain-specific terms
- **tools_reference.md** — Pros/cons/differentiators for every tool across the pipeline
