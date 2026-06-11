# P. falciparum Alba Family scRNA-seq Analysis Pipeline

**Author:** Prabhav   
**Date:** June 2026

---

## What This Project Does

This pipeline analyzes **single-cell RNA sequencing (scRNA-seq)** data from the malaria parasite *Plasmodium falciparum* (NF54 strain) to visualize the expression of the **PfAlba protein family (Alba1–4)** across the parasite's intraerythrocytic life cycle stages: ring, trophozoite, schizont, and gametocyte.

The dataset is extracted from the **Malaria Cell Atlas** (Dogga et al., *Science* 2024) — a landmark atlas of ~37,000 *P. falciparum* cells covering the complete intraerythrocytic cycle including sexual development. The biological interpretation of results is grounded in **Acharya et al., *Microbiology Spectrum* 2025** — a study from Dr. Shruti Sridhar Vembar's lab at IBAB that directly characterizes PfAlba2, PfAlba3, and PfAlba4 function during asexual blood stage development.

---

## Biological Question

> *How are PfAlba1, PfAlba2, PfAlba3, and PfAlba4 expressed across the ring, trophozoite, schizont, and gametocyte stages of the P. falciparum NF54 intraerythrocytic cycle?*

---

## About the Dataset

| Property | Details |
|----------|---------|
| **Source paper** | Dogga et al., *Science* 384, eadj4088 (2024) |
| **Dataset ID** | pf-ch10x-set4 |
| **Sequencing platform** | 10x Genomics Single Cell 3′ v3.1 |
| **Lab strains profiled** | NF54 and 7G8 |
| **This analysis uses** | NF54 strain only |
| **Data available at** | [malariacellatlas.org](https://www.malariacellatlas.org) |
| **ENA accession** | PRJEB58790 |

The raw counts and cell metadata were downloaded from the Malaria Cell Atlas resource. The metadata already contains pre-computed stage annotations (`STAGE_LR`, `STAGE_HR`), cluster assignments (`CLUSTER`), and dimensionality reduction coordinates (`PC_1/2/3`, `UMAP_1/2/3`) assigned by Dogga et al. using a consensus method combining Seurat clustering, scmap reference mapping, and marker gene identification. **No re-clustering or UMAP computation is needed in this pipeline** — we use the existing annotations directly.

---

## Key Findings

| Gene | Total UMIs | Cells Expressing | % of Cells | Pattern |
|------|-----------|-----------------|-----------|---------|
| PfAlba1 | 15,055 | 9,442 / 35,529 | 26.6% | High, trophozoite peak |
| PfAlba2 | 28,279 | 15,234 / 35,529 | 42.9% | Highest, ring + trophozoite |
| PfAlba3 | 510 | 484 / 35,529 | 1.4% | Near-absent — scRNA-seq dropout |
| PfAlba4 | 138 | 52 / 35,529 | 0.15% | Near-absent — scRNA-seq dropout |

Alba3 and Alba4 low detection is confirmed as **scRNA-seq dropout**, not a pipeline error. Two lines of evidence support this:

1. **Raw count verification:** Alba3 max = 3 UMIs per cell; Alba4 detected in only 52 cells with no stage preference — random scatter across all 4 stages.
2. **Orthogonal bulk RNA-seq (TPM):** TPM analysis of the same parasite shows Alba3 peaks at 40–44 hpi and Alba4 at late schizont/merozoite. The discrepancy is explained by dropout — each droplet captures only ~5–15% of a cell's mRNA, so low-abundance transcripts frequently get zero counts.

This is further supported by Dogga et al. 2024, which showed that **female gametocytes accumulate translationally repressed mRNAs** held in mRNP complexes (module E1, Fig. 1D). Alba3 and Alba4 transcripts may be held in such complexes — present but not actively translating and therefore at very low steady-state levels detectable by scRNA-seq.

---

## Pipeline Flowchart

```
[ Raw Counts CSV ]                               [ Metadata CSV ]
  (40k Cells x 5k Genes)                    (Stage, Strain, Cluster — from
  from Malaria Cell Atlas                    Dogga et al. 2024 annotation)
            |                                               |
            | 1. data.table::fread()                        | 1. base::read.csv()
            |    (Fast memory read)                         |    (row.names = 1)
            v                                               |
  [ Dense Data Frame ]                                      |
            |                                               |
            | 2. Extract gene column                        |
            | 3. as.matrix()                                |
            v                                               |
   [ Standard Matrix ]                                      |
            |                                               |
            | 4. Matrix::Matrix(sparse = TRUE)              |
            |    (Drops all zeros from RAM)                 |
            v                                               |
    [ Sparse Matrix ]                                       |
            |                                               |
            | 5. CreateSeuratObject()                       |
            |    (Applies base QC filters)                  |
            v                                               |
  [ Base Seurat Object ] <----------------------------------+
            |                       6. AddMetaData()
            |                          (Attaches Dogga et al. stage labels
            |                           by matching barcodes)
            v
[ Annotated Seurat Object ]
   (Contains all 40k cells)
            |
            | 7. subset(subset = STRAIN == "NF54")
            |    (Isolates NF54 — one of two lab strains
            |     profiled in Dogga et al. 2024)
            v
 [ Final NF54 Seurat Object ]
   (35,529 cells · 5,274 genes)
            |
            | 8. rm() & gc()
            |    (Flush old data from RAM)
            v
[ NormalizeData() ]
   LogNormalize · scale.factor = 10,000
            |
            | 9. Verify gene names in rownames()
            |    (Dataset uses hyphens: PF3D7-0813600)
            v
[ Confirm PfAlba IDs ]
   PF3D7-0813600 / -1006800 / -1418300 / -0600200
            |
            | 10. Visualize expression across STAGE_LR
            |     (Pre-annotated by Dogga et al. — no re-clustering needed)
            v
  ┌─────────────────┐     ┌─────────────────┐
  │   DotPlot()     │     │   VlnPlot()     │
  │ % expressing +  │     │  distribution   │
  │ avg expression  │     │   per stage     │
  │ group.by=STAGE  │     │  group.by=STAGE │
  └─────────────────┘     └─────────────────┘
            |
            v
  [ Bar Plot — Mean Expression by Stage ]
  (All 4 Albas + focused Alba3/Alba4 plot)
```

---

## Project Structure

```
pf-alba-scrna/
│
├── README.md                    # This file
│
├── main_pipeline.R              # Complete pipeline in one file (run this)
│
├── scripts/                     # Modular sub-scripts (one per section)
│   ├── 01_load_data.R           # Read CSV files, build sparse matrix
│   ├── 02_build_seurat.R        # Create Seurat object, add metadata, subset NF54
│   ├── 03_normalize.R           # Normalize data, verify gene names
│   ├── 04_alba_analysis.R       # Expression statistics and cross-validation
│   └── 05_visualize.R           # All plots (VlnPlot, DotPlot, BarPlot)
│
├── data/                        # Input files go here
│   ├── pf-ch10x-set4-ch10x-raw.csv     # Raw count matrix (from MCA)
│   └── pf-ch10x-set4-ch10x-data.csv    # Cell metadata (from MCA)
│
├── plots/                       # All output PDFs saved here
│   ├── alba_vlnplot.pdf
│   ├── alba_dotplot_unscaled.pdf
│   ├── alba_barplot.pdf
│   └── alba_barplot_low.pdf
│
└── results/
    └── RESULTS.md               # Documented findings and interpretation
```

---

## Requirements

### R Version
R 4.6.0 or higher

### Required Packages

```r
install.packages("BiocManager")
BiocManager::install(c("Seurat", "data.table", "Matrix"))
install.packages(c("ggplot2", "dplyr", "tidyr"))
```



---

## Input Files

| File | Description | Format |
|------|-------------|--------|
| `pf-ch10x-set4-ch10x-raw.csv` | Raw UMI count matrix from Malaria Cell Atlas | Rows = genes (PF3D7- IDs), Columns = cell barcodes |
| `pf-ch10x-set4-ch10x-data.csv` | Cell metadata with Dogga et al. stage annotations | Rows = barcodes, Columns = STRAIN, STAGE_LR, CLUSTER, UMAP etc. |

> **Important:** Gene names in this dataset use **hyphens**, not underscores.  
> Correct: `PF3D7-0813600` | Wrong: `PF3D7_0813600`  
> This is specific to how the Malaria Cell Atlas formatted the gene IDs in this export.

---

## Output Files

| File | Description |
|------|-------------|
| `alba_vlnplot.pdf` | Violin plot — expression distribution per stage |
| `alba_dotplot_unscaled.pdf` | Dot plot — % cells expressing + avg expression (unscaled) |
| `alba_barplot.pdf` | Bar plot — mean expression of all 4 Albas by stage |
| `alba_barplot_low.pdf` | Bar plot — Alba3 + Alba4 only, rescaled y-axis |

---

## How to Run

### Option 1: Run the full pipeline at once

```r
source("main_pipeline.R")
```

### Option 2: Run step by step (recommended for first time)

```r
source("scripts/01_load_data.R")
source("scripts/02_build_seurat.R")
source("scripts/03_normalize.R")
source("scripts/04_alba_analysis.R")
source("scripts/05_visualize.R")
```

> Set working directory to the project folder first.  
> Check: `getwd()` | Set: `setwd("path/to/pf-alba-scrna")`

---

## Key Technical Notes

1. **Gene name format:** This dataset uses hyphens (`PF3D7-0813600`), not underscores. Always verify with `head(rownames(obj), 5)` before searching.
2. **Stage annotations are pre-computed:** The `STAGE_LR` column comes from Dogga et al. 2024's consensus annotation pipeline. No PCA, UMAP, or clustering is needed.
3. **DotPlot scaling:** Always use `scale = FALSE`. With only 4 stage groups, default z-score scaling produces misleading -1 to +1 color ranges.
4. **ScaleData not required:** Skipped intentionally. Only normalization is needed for expression visualization. ScaleData + PCA would be needed only for discovering unknown cell types — which is not our goal here.
5. **Alba3/4 dropout:** Confirmed by raw count verification AND bulk TPM data. Not a pipeline error.
6. **Why no re-clustering:** Dogga et al. used a rigorous consensus method (Seurat + scmap + marker genes) to assign stages. Repeating this would require the same reference datasets and would not improve on their annotation.

---

## References

### Data Source
Dogga SK, Rop JC, Cudini J, Farr E, Dara A, Ouologuem D, Djimdé AA, Talman AM, Lawniczak MKN.  
*A single-cell atlas of sexual development in Plasmodium falciparum.*  
**Science** 384, eadj4088 (2024). DOI: [10.1126/science.adj4088](https://doi.org/10.1126/science.adj4088)  
Dataset: `pf-ch10x-set4` | Available at: [malariacellatlas.org](https://www.malariacellatlas.org) | ENA: PRJEB58790

### Biological Context
Acharya D, Bavikatte AN, Ashok VV, Hegde SR, Macpherson CR, Scherf A, Vembar SS.  
*Ectopic overexpression of Plasmodium falciparum DNA-/RNA-binding Alba proteins misregulates virulence gene homeostasis during asexual blood development.*  
**Microbiology Spectrum** 13(3) (2025). DOI: [10.1128/spectrum.00885-24](https://doi.org/10.1128/spectrum.00885-24)
