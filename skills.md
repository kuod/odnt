# skills.md — Tutorial Series Generation Spec

This file defines how to generate a new tutorial deep-dive series in this repo. Reference it when adding a new topic (e.g., "bulk RNA-seq", "scATAC-seq", "proteomics", "spatial transcriptomics extended") so the output is consistent with existing series.

---

## Audience Profile

- **Background:** Traditional bioinformatics (NGS pipelines, bulk RNA-seq, genome assembly, etc.)
- **Goal:** Get up to speed on a specific sc-omics topic or tool ecosystem
- **Calibration:** Assume familiarity with sequencing, count matrices, R/Bioconductor concepts — but introduce scRNA-seq-specific ideas from first principles
- **Style:** Practitioner, not student. Skip hand-holding; lead with working code and explain design decisions, not syntax

---

## Repo Conventions

### Directory layout
Each series gets its own folder at repo root:
```
perturb-101/
├── perturb_seq/        # existing: Perturb-seq curriculum
├── theis_ecosystem/    # existing: scverse tool refresher
├── <new_series>/       # new series goes here
│   ├── README.md
│   └── <prefix><NN>_<slug>.ipynb
└── data/               # git-ignored; all series write here
```

### Naming conventions
| Item | Convention | Example |
|------|-----------|---------|
| Series folder | `snake_case` | `spatial_extended/`, `atac_seq/` |
| Notebook prefix | Short uppercase abbreviation | `SP` for spatial, `AT` for ATAC |
| Notebook filename | `<PREFIX><NN>_<slug>.ipynb` | `SP00_spatialdata_intro.ipynb` |
| Data paths | Always `../data/<filename>` | `../data/visium_processed.h5ad` |
| Intermediate files | `.h5ad` for AnnData, `.h5mu` for MuData | `../data/pbmc_cite.h5mu` |

### Notebook cell conventions
- First cell: **markdown header** (title, tools, dataset, paper DOI link)
- Second cell: **imports + version print**
- Each major section: **markdown header** (`## N. Section Title`) before the code
- Last cell: **markdown summary table** of key functions → result locations
- Code style: prefer concise, commented cells (~15–25 lines max per code cell)
- Never use bare `data/` paths — always `../data/`

---

## Series Structure Template

A complete series = `README.md` + 4–8 notebooks. Typical arc:

```
00 — Foundation / Data Container        (always first; establishes the AnnData/object pattern)
01 — Core Analysis Pipeline             (the "hello world" workflow end-to-end)
02 — Advanced Method A                  (tool-specific, builds on 01)
03 — Advanced Method B
...
NN — Integration / Synthesis            (connect back to other series or broader context)
```

Minimum viable series: **3 notebooks** (foundation, core pipeline, one advanced topic).
Maximum recommended: **8–10 notebooks** (beyond that, split into sub-series).

---

## Notebook Anatomy

Every notebook follows this skeleton:

```markdown
# <PREFIX><NN> — <Title>

**Tool(s):** `package-name`
**Dataset:** <Name> (<species>, <assay>, <cell count>)
**Paper:** [Author et al. YYYY, Journal](https://doi.org/<DOI>)

---

## What problem does this solve?
(2–4 sentences: the biological or technical motivation)

## How does it work?
(Conceptual diagram or table comparing approaches)
```

```python
# imports + version print cell
```

```markdown
## 1. Load Data
```

```python
# dataset loading — auto-download with pooch/built-in dataset function
# never hardcode a URL the user must manually visit
```

```markdown
## 2. [Core step]
```

...

```markdown
---
## Summary

| Step | Function | Result stored in |
|------|----------|-----------------|
...

**Next:** <NN+1> — <next notebook title>
```

---

## Dataset Selection Criteria

Pick datasets that satisfy all three:

1. **Auto-downloadable** — via a built-in function (`sc.datasets.*`, `sq.datasets.*`, `scvi.data.*`, `pooch.retrieve`) with no manual browser steps
2. **Small** — runs on a laptop without GPU in <10 min. Target: <10k cells for notebooks that train models; any size for read-only exploration
3. **Standard** — a dataset the community recognizes (PBMC 3k, mouse brain Visium, pancreas E15.5, etc.) so learners can cross-reference other tutorials

Preferred datasets by topic area:

| Topic | Dataset | Access |
|-------|---------|--------|
| General scRNA-seq | PBMC 3k | `sc.datasets.pbmc3k()` |
| Integration / batches | Heart Cell Atlas subsampled | `scvi.data.heart_cell_atlas_subsampled()` |
| Spatial (Visium) | Mouse brain H&E | `sq.datasets.visium_hne_adata()` |
| Trajectory | Mouse pancreas E15.5 | `scv.datasets.pancreas()` |
| CITE-seq | PBMC 5k multiome | `mu.datasets.pbmc5k()` |
| Perturb-seq | Norman 2019 CRISPRa | downloaded in `perturb_seq/03` |
| scATAC-seq | 10k PBMC ATAC | `snapatac2.datasets.pbmc5k()` |
| Spatial HD / in situ | Mouse brain Xenium | `spatialdata_io` loaders |

---

## Topic Specification Format

When requesting a new series, provide this spec block. Claude will generate the full series from it.

```
## New Series Request

**Topic:** <one-line description>
**Series folder:** <snake_case folder name>
**Notebook prefix:** <2-3 uppercase letters>
**Audience goal:** <what the learner can do after this series>

**Notebooks:**
| # | Title | Tool(s) | Dataset | Core concept |
|---|-------|---------|---------|--------------|
| 00 | ... | ... | ... | ... |
| 01 | ... | ... | ... | ... |

**Key tools to cover:**
- `tool-name` — one-line description, link to paper

**Datasets:**
- Primary: <name> (<access function>)
- Reference (if needed): <name>

**Connections to existing series:**
- Prereq: <series/notebook> (e.g., "theis_ecosystem/T00 for AnnData basics")
- Downstream: <what this enables>

**Known gaps / explicit non-goals:**
- <thing to skip and why>
```

---

## Quality Checklist

Before committing a new series, verify:

- [ ] All notebooks parse as valid JSON: `python3 -c "import json; json.load(open('NB.ipynb'))"`
- [ ] No bare `data/` paths — all use `../data/`
- [ ] Each notebook has a title markdown cell with tool, dataset, paper DOI
- [ ] Each notebook ends with a summary table of key functions
- [ ] Datasets load without hardcoded URLs (use built-in loaders or `pooch`)
- [ ] `README.md` in the series folder has the full notebook table
- [ ] `environment.yml` at repo root has commented install instructions for new deps
- [ ] `CLAUDE.md` at repo root has the new folder in the structure diagram
- [ ] `clean.sh` exists in the series folder listing all `../data/` files the series writes

## clean.sh Convention

Every series folder must include a `clean.sh` that:
- **Defaults to dry run** — prints files and sizes that would be deleted, exits without touching anything
- **Deletes with `--force`** — `./clean.sh --force` removes all listed files

Template (copy and fill in `FILES=()`):

```bash
#!/usr/bin/env bash
# clean.sh — remove data files generated by the <series> series
# Usage: ./clean.sh [--force]
set -euo pipefail
FORCE=false; [[ "${1:-}" == "--force" ]] && FORCE=true
DATA_DIR="$(cd "$(dirname "$0")/../data" 2>/dev/null && pwd)" || { echo "data/ not found"; exit 0; }

FILES=(
  "$DATA_DIR/<file1>.h5ad"   # created by notebook NN
  "$DATA_DIR/<file2>.h5ad"
)

FOUND=(); for f in "${FILES[@]}"; do [[ -e "$f" ]] && FOUND+=("$f"); done
[[ ${#FOUND[@]} -eq 0 ]] && { echo "Nothing to clean."; exit 0; }

echo "<series> data files to remove:"
for f in "${FOUND[@]}"; do echo "  [$(du -sh "$f" | cut -f1)]  $f"; done
[[ "$FORCE" == false ]] && { echo -e "\nDry run — run with --force to delete."; exit 0; }
for f in "${FOUND[@]}"; do rm -rf "$f" && echo "Deleted: $f"; done
```

List **only files that `../data/`** — do not list library caches (squidpy, scvi-tools, etc. manage their own cache dirs). Include both raw downloads and derived intermediates since all can be regenerated.

---

## Existing Series Reference

| Series | Prefix | Notebooks | Primary tools |
|--------|--------|-----------|---------------|
| `perturb_seq/` | `NN_` | 00–14 | scanpy, pertpy, pydeseq2, mixscape, milopy |
| `theis_ecosystem/` | `T0N_` | T00–T08 | anndata, scanpy, scvi-tools, squidpy, cellrank, muon, liana, cell2location, moscot, sccoda |
| `spatial_transcriptomics/` | `SP0N_` | SP00–SP04 | spatialdata, squidpy, banksy-py, GraphST, moscot |

---

## Example: Adding a New Series

**Request:**
```
## New Series Request

**Topic:** Single-cell ATAC-seq analysis
**Series folder:** atac_seq/
**Notebook prefix:** AT
**Audience goal:** Process raw fragment files through peak calling,
  dimensionality reduction, motif enrichment, and linking peaks to genes

**Notebooks:**
| # | Title | Tool(s) | Dataset |
|---|-------|---------|---------|
| AT00 | ATAC Data Structures | snapatac2, muon | PBMC 5k ATAC |
| AT01 | QC and Preprocessing | snapatac2 | PBMC 5k ATAC |
| AT02 | Dimensionality Reduction and Clustering | snapatac2 | PBMC 5k ATAC |
| AT03 | Peak Calling and Motif Enrichment | snapatac2, pyfaidx | PBMC 5k ATAC |
| AT04 | Gene Activity and RNA-ATAC Integration | snapatac2, muon | PBMC 10k Multiome |

**Connections:** prereq theis_ecosystem/T04 for MuData; downstream enables theis_ecosystem/T06 cell2location with ATAC
```

Claude will generate all notebooks, README, and environment updates from this spec.
