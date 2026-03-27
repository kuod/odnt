# Theis Lab / scverse Ecosystem Tutorial Series

A hands-on refresher tour of the [scverse](https://scverse.org/) ecosystem — the suite of interoperable Python packages for single-cell omics anchored by the Theis Lab at Helmholtz Munich.

## Why this series?

All tools share a common data contract: **AnnData / MuData** objects that pass cleanly between packages. Once you know the ecosystem's conventions, moving from QC → integration → spatial → trajectories → multi-modal is mostly learning new function signatures, not new mental models.

## Notebooks

| # | Notebook | Tool(s) | Dataset | Key concepts |
|---|----------|---------|---------|--------------|
| T00 | AnnData + Scanpy | `anndata`, `scanpy` | PBMC 3k (10x) | Data structure anatomy, full preprocessing → clustering → DE pipeline |
| T01 | Deep Generative Models | `scvi-tools` | PBMC 3k (simulated batches) + Heart Cell Atlas | scVI latent space, batch integration, scANVI label transfer |
| T02 | Spatial Transcriptomics | `squidpy` | Mouse brain Visium H&E | Spatial graphs, nhood enrichment, Moran's I, ligand-receptor |
| T03 | Trajectory Inference | `cellrank` (v2), `scvelo` | Mouse pancreas (E15.5) | RNA velocity, transition kernels, macrostates, fate probabilities |
| T04 | Multi-Modal Analysis | `muon`, `mofapy2` | PBMC CITE-seq (RNA + protein) | MuData structure, per-modality QC, MOFA+ factors, WNN embedding |

## Setup

Install the full ecosystem on top of the base `perturb-101` environment:

```bash
conda activate perturb-101
pip install scvi-tools>=1.1 squidpy>=1.6 cellrank>=2.0 scvelo>=0.3 \
            scirpy>=0.22 mofapy2>=0.7 spatialdata>=0.2
```

Or recreate from the updated `environment.yml` in the repo root.

## How the packages relate

```
            ┌─────────────────────────────────────────────┐
            │             AnnData / MuData                │
            │        (the shared data contract)           │
            └───────┬──────────────────────┬──────────────┘
                    │                      │
          ┌─────────▼──────┐    ┌──────────▼──────────┐
          │     scanpy     │    │      muon            │
          │ (core pp/viz)  │    │  (multi-modal)       │
          └─────┬──────────┘    └──────────┬───────────┘
                │                          │
      ┌─────────┼──────────┐               │ MOFA+
      │         │          │               │
 ┌────▼───┐ ┌───▼────┐ ┌───▼──────┐  ┌────▼──────┐
 │scvi-   │ │squidpy │ │cellrank  │  │mofapy2    │
 │tools   │ │(space) │ │(traject.)│  │(factors)  │
 └────────┘ └────────┘ └──────────┘  └───────────┘
```

## Notes

- Run notebooks in order T00 → T04; each builds on the previous AnnData conventions.
- T01 requires ~4 GB RAM for scVI training; reduce `max_epochs` if needed.
- T03 requires `scvelo` which has a heavy PyTorch dependency.
- All datasets are downloaded automatically on first run.
