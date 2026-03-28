# Spatial Transcriptomics Deep Dive

A focused series on spatial transcriptomics — from the modern data container to sub-cellular resolution platforms. Builds on the basics in `theis_ecosystem/T02` (squidpy) and `theis_ecosystem/T06` (cell2location); those are prerequisites, not repeated here.

## Series map

| # | Notebook | Tool(s) | Dataset | Key concepts |
|---|----------|---------|---------|--------------|
| SP00 | SpatialData — The Modern Container | `spatialdata`, `spatialdata-io` | Synthetic blobs + mouse brain Visium | Zarr store, coordinate systems, images/labels/shapes/points/tables, squidpy compat |
| SP01 | Visium Analysis In Depth | `scanpy`, `squidpy` | Mouse brain Visium H&E | Spatial QC, spatial DE between tissue regions, H&E image feature extraction |
| SP02 | Spatial Domain Detection | `banksy-py` / `GraphST` | Mouse brain Visium H&E | Tissue niche discovery using spatial context; expression-only vs. spatially-aware clustering |
| SP03 | Sub-cellular Resolution | `squidpy`, `spatialdata` | Mouse brain seqFISH+ | Single-molecule transcript model, cell polygons, transcript density, single-cell spatial analysis |
| SP04 | Multi-sample Spatial Integration | `moscot`, `harmonypy` | Mouse brain Visium (two simulated slides) | Slide alignment, batch effects, cross-sample spatial DE, consensus annotations |

## Prerequisites

- `theis_ecosystem/T00` — AnnData fundamentals
- `theis_ecosystem/T02` — squidpy basics (spatial neighbors, Moran's I, nhood enrichment, ligrec)
- `theis_ecosystem/T06` — cell2location spatial deconvolution

## Installation

```bash
conda activate perturb-101
pip install spatialdata spatialdata-io spatialdata-plot \
            banksy-py graphst \
            squidpy>=1.6 moscot
```

## What's NOT covered here (see theis_ecosystem instead)

| Topic | Where |
|-------|-------|
| Spatial neighbors, Moran's I, nhood enrichment, ligrec | `theis_ecosystem/T02` |
| cell2location two-stage deconvolution | `theis_ecosystem/T06` |
| Moscot SpatialAlignmentProblem concept | `theis_ecosystem/T07` |
| LIANA cell-cell communication | `theis_ecosystem/T05` |

## Technology landscape

```
Resolution:     ←────────────────────────────────────────────→
                Sub-cellular          Spot-level     Tissue
                (single-molecule)     (bulk per spot)

Platforms:      Xenium  MERFISH       Visium          Slide-seq
                Vizgen  seqFISH+      Visium HD       ST (old)
                CosMx   FISH-seq      Stereo-seq

Python tools:   spatialdata           squidpy         squidpy
                squidpy               scanpy          scanpy
```
