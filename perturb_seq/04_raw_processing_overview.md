# 04 — Raw Data Processing Overview

This document covers the computational pipeline that converts raw FASTQ files from a Perturb-seq run into the count matrices you load into Python. You will not run this pipeline in this course — it requires HPC access, a genome index, and tens of GB of FASTQ data — but understanding it is essential for interpreting what is (and is not) in your AnnData object.

---

## What comes off the sequencer

A Perturb-seq run produces **two sets of FASTQs** per sample:

```
GEX library:
  R1: 28 bp  — cell barcode (16 bp) + UMI (12 bp, 10x v3)
  R2: 90+ bp — cDNA read (mapped to transcriptome)

Guide capture library:
  R1: 28 bp  — same cell barcode + UMI structure
  R2: 28 bp  — guide barcode sequence (short amplicon)
```

Both libraries share the same cell barcode whitelist, which is how they are later merged by cell.

---

## Tool 1: Cell Ranger (10x Genomics)

Cell Ranger is the official 10x Genomics pipeline. It handles both GEX and guide capture.

### GEX processing

```bash
cellranger count \
  --id=sample_gex \
  --transcriptome=/path/to/refdata-gex-GRCh38 \
  --fastqs=/path/to/gex_fastqs \
  --sample=SampleName \
  --localcores=16 \
  --localmem=64
```

Output: `filtered_feature_bc_matrix/` containing `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`

### CRISPR guide capture

For guide capture, you need a **feature reference CSV** listing every guide barcode:

```csv
id,name,read,pattern,sequence,feature_type
guide_CBFB_1,CBFB_guide1,R2,5PNNNNNNNNNN(BC)NNNNNNN,ACGTACGT,CRISPR Guide Capture
...
```

```bash
cellranger count \
  --id=sample_crispr \
  --transcriptome=/path/to/refdata-gex-GRCh38 \
  --fastqs=/path/to/gex_fastqs,/path/to/crispr_fastqs \
  --sample=GEX_SampleName,CRISPR_SampleName \
  --feature-ref=/path/to/feature_reference.csv \
  --localcores=16
```

This produces a combined output where the feature matrix includes both gene counts and guide UMI counts. In Python:

```python
import scanpy as sc
adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5", gex_only=False)
# adata now contains both gene features and CRISPR guide features
# Split them:
gex = adata[:, adata.var["feature_types"] == "Gene Expression"]
guides = adata[:, adata.var["feature_types"] == "CRISPR Guide Capture"]
```

---

## Tool 2: STARsolo

STARsolo is the single-cell mode of the STAR aligner. More configurable than Cell Ranger, free from licensing constraints, and produces compatible output.

### Build genome index

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/star_index \
  --genomeFastaFiles /path/to/GRCh38.fa \
  --sjdbGTFfile /path/to/gencode.v44.gtf \
  --runThreadN 16
```

### GEX mapping

```bash
STAR --soloType CB_UMI_Simple \
  --soloCBwhitelist 3M-february-2018.txt \  # 10x v3 whitelist
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloFeatures Gene \
  --readFilesIn R2.fastq.gz R1.fastq.gz \  # Note: R2 first (cDNA), R1 second (barcode)
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --genomeDir /path/to/star_index \
  --outFileNamePrefix sample_ \
  --runThreadN 16
```

For guide capture: STARsolo does not have native feature barcoding support. Process guide FASTQs separately with a custom script that matches R2 to your feature reference using fuzzy matching (≤1 mismatch), then join by cell barcode.

**Load STARsolo output into scanpy:**

```python
adata = sc.read_10x_mtx("Solo.out/Gene/filtered/")
```

---

## Tool 3: kallisto | bustools via kb-python

The `kb-python` wrapper provides the most streamlined command-line interface for quantification. For CRISPR feature barcoding, use the `kite` workflow.

### Install

```bash
pip install kb-python
```

### Build KITE index for guide barcodes

KITE (Kallisto Indexing and Tag Extraction) treats guide barcodes like small "transcripts":

```bash
# Create a FASTA of guide barcode sequences
kb ref \
  -i index.idx \
  -g t2g.txt \
  -f1 feature_barcodes.fa \
  --workflow kite \
  feature_barcodes.fa
```

### Quantify guide capture library

```bash
kb count \
  --workflow kite \
  -i index.idx \
  -g t2g.txt \
  -x 10xv3 \
  -o guide_output/ \
  R1.fastq.gz R2.fastq.gz
```

This produces `cells_x_genes.bus`, processed to `output.mtx` — a cell × guide UMI count matrix.

**Load into AnnData:**

```python
import anndata
import scipy.io
import pandas as pd

guide_matrix = scipy.io.mmread("guide_output/counts_unfiltered/cells_x_genes.mtx").T.tocsr()
barcodes = pd.read_csv("guide_output/counts_unfiltered/cells_x_genes.barcodes.txt", header=None)[0]
features = pd.read_csv("guide_output/counts_unfiltered/cells_x_genes.genes.txt", header=None)[0]

guide_adata = anndata.AnnData(
    X=guide_matrix,
    obs=pd.DataFrame(index=barcodes),
    var=pd.DataFrame(index=features)
)
```

---

## Multi-modal AnnData with MuData

After processing both libraries separately, merge by cell barcode using `muon`:

```python
import muon as mu

mdata = mu.MuData({
    "rna": gex_adata,
    "crispr": guide_adata
})

# Intersect cell barcodes
mu.pp.intersect_obs(mdata)

# Access each modality
mdata["rna"]    # AnnData with gene counts
mdata["crispr"] # AnnData with guide UMI counts
```

For this course, the downloaded Norman 2019 and Replogle 2022 datasets already have guide metadata integrated in `adata.obs`. You will not need to merge modalities manually.

---

## Output format comparison

| Tool | Output format | Filtered/unfiltered cells | Notes |
|---|---|---|---|
| Cell Ranger | `filtered_feature_bc_matrix/` (MTX), `.h5` | Both | Most compatible with 10x data |
| STARsolo | `Solo.out/Gene/filtered/` (MTX) | Both | No native CRISPR support |
| kb-python | `counts_unfiltered/` (MTX + BUS) | Unfiltered; use knee plot | Best for non-10x chemistries |

All three can be loaded into AnnData via `sc.read_10x_mtx()` or `anndata.read_mtx()`.

---

## Cell filtering (knee plot / EmptyDrops)

All three tools report an "unfiltered" matrix (all barcodes detected) and a "filtered" matrix (barcodes above the expected cell threshold). The filtering uses:

- **Cell Ranger**: their internal `cellranger_cells` algorithm
- **STARsolo**: knee plot heuristic or EmptyDrops
- **kb-python**: no filtering by default; use `kneepoint` or `emptydrops` from the `barcodeutils` package

For Perturb-seq, use the **filtered** (cell-containing) matrix. Empty droplets will have near-zero guide UMIs, which creates noise in guide assignment.

---

## Approximate runtimes (10,000 cells, 250M GEX reads)

| Tool | Wall time | Memory |
|---|---|---|
| Cell Ranger | 4–6 hours | 32–64 GB |
| STARsolo | 1–2 hours | 32 GB |
| kb-python | 15–30 min | 16 GB |

---

**Next:** [05_quality_control.ipynb](05_quality_control.ipynb)
