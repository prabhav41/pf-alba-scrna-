# =============================================================================
# FILE:    main_pipeline.R
# PROJECT: P. falciparum Alba Family scRNA-seq Analysis
# AUTHOR:  Prabhav
# DATE:    June 2026
#
# PURPOSE:
#   This script analyzes single-cell RNA-seq data from Plasmodium falciparum
#   NF54 strain to visualize expression of the PfAlba protein family (Alba1-4)
#   across life cycle stages: ring, trophozoite, schizont, gametocyte.
#
# INPUT FILES REQUIRED (place in data/ folder):
#   - pf-ch10x-set4-ch10x-raw.csv   (raw count matrix: genes x cells)
#   - pf-ch10x-set4-ch10x-data.csv  (cell metadata: barcodes x annotations)
#
# OUTPUT FILES GENERATED (saved to plots/ folder):
#   - alba_vlnplot.pdf
#   - alba_dotplot_unscaled.pdf
#   - alba_barplot.pdf
#   - alba_barplot_low.pdf
#
# HOW TO RUN:
#   Set working directory to project root, then:
#   source("main_pipeline.R")
# =============================================================================


# =============================================================================
# SECTION 1: LOAD LIBRARIES
# All packages must be installed before running. See README.md for install
# instructions. Each package serves a specific purpose:
#   - Seurat:      single-cell analysis framework
#   - data.table:  fast CSV reading (fread is faster than read.csv for large files)
#   - Matrix:      sparse matrix format to save RAM
#   - ggplot2:     plotting
#   - dplyr:       data manipulation (group_by, summarise, filter)
#   - tidyr:       reshaping data (pivot_longer)
# =============================================================================

library(Seurat)
library(data.table)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

message("=== Libraries loaded ===")


# =============================================================================
# SECTION 2: DEFINE INPUT FILE PATHS
# Change these paths if your files are in a different location.
# Both files must be in the same working directory, or provide full paths.
# =============================================================================

raw_counts_file <- "data/pf-ch10x-set4-ch10x-raw.csv"
meta_file       <- "data/pf-ch10x-set4-ch10x-data.csv"
plots_dir       <- "plots/"

# Create plots directory if it doesn't exist
if (!dir.exists(plots_dir)) dir.create(plots_dir)


# =============================================================================
# SECTION 3: READ RAW COUNT MATRIX
#
# The raw CSV has:
#   - Rows    = genes (PlasmoDB IDs, e.g. PF3D7-0813600)
#   - Columns = cell barcodes (e.g. D0_AAACCCAGTCTGTCCT-1)
#   - Values  = UMI counts (how many times each gene was detected per cell)
#
# fread() from data.table is used because this file is ~500MB and read.csv
# would be much slower and use more RAM.
#
# We then extract the gene names from column 1, convert to a standard matrix,
# and assign gene names as row names.
# =============================================================================

message("Reading raw count matrix... (this may take 1-2 minutes)")

raw_df <- fread(raw_counts_file)

# Column 1 contains gene IDs — extract before converting to numeric matrix
genes <- raw_df[[1]]

# Drop the gene ID column and convert remaining data to a numeric matrix
raw_mat <- as.matrix(
  raw_df[, -1, with = FALSE]
)

# Assign gene names as row names so Seurat can look up genes by name
rownames(raw_mat) <- genes

message("Raw matrix dimensions: ", nrow(raw_mat), " genes x ", ncol(raw_mat), " cells")


# =============================================================================
# SECTION 4: CONVERT TO SPARSE MATRIX
#
# scRNA-seq data is very sparse — most gene-cell combinations are zero
# (the gene was not detected in that cell). Storing all these zeros in a
# standard matrix wastes RAM. A sparse matrix stores only the non-zero values.
#
# For this dataset (~5200 genes x 40000 cells = 208 million values), this
# reduces RAM usage from ~4GB to ~600MB.
#
# After creating the sparse matrix, we delete the dense matrix and run
# garbage collection to free up that RAM immediately.
# =============================================================================

message("Converting to sparse matrix to save RAM...")

sparse_mat <- Matrix(
  raw_mat,
  sparse = TRUE
)

# Free the dense matrix from RAM — we no longer need it
rm(raw_df, raw_mat)
gc()

message("Sparse matrix created. RAM freed.")


# =============================================================================
# SECTION 5: CREATE SEURAT OBJECT
#
# Seurat is the standard framework for single-cell analysis in R.
# CreateSeuratObject() does two things:
#   1. Stores the count matrix in a structured object
#   2. Applies basic quality control filters:
#      - min.cells = 3:    only keep genes detected in at least 3 cells
#      - min.features = 200: only keep cells that express at least 200 genes
#
# Note: For P. falciparum (~5200 genes total), min.features=200 means we keep
# cells expressing at least 3.8% of the genome. This is appropriate for parasite
# data but more strict proportionally than for human data.
# =============================================================================

message("Creating Seurat object...")

seurat_obj <- CreateSeuratObject(
  counts    = sparse_mat,
  project   = "Plasmodium_scRNA",
  min.cells = 3,
  min.features = 200
)

message("Seurat object created: ", ncol(seurat_obj), " cells, ", nrow(seurat_obj), " genes")


# =============================================================================
# SECTION 6: LOAD AND ATTACH METADATA
#
# The metadata CSV contains biological annotations for each cell barcode:
#   - STRAIN:    which parasite strain (NF54, 3D7, etc.)
#   - STAGE_LR:  life cycle stage (ring, trophozoite, schizont, gametocyte)
#   - STAGE_HR:  higher resolution stage annotation
#   - CLUSTER:   pre-computed cluster assignment
#   - DAY, HOST, SOURCE: experimental metadata
#   - PC_1/2/3, UMAP_1/2/3: pre-computed dimensionality reduction coordinates
#
# row.names = 1 ensures the barcode column becomes the row name,
# which Seurat uses to match cells between the count matrix and metadata.
# =============================================================================

message("Loading cell metadata...")

meta_df <- read.csv(
  meta_file,
  row.names = 1
)

seurat_obj <- AddMetaData(
  object   = seurat_obj,
  metadata = meta_df
)

message("Metadata attached. Columns available: ", paste(colnames(meta_df), collapse=", "))


# =============================================================================
# SECTION 7: SUBSET NF54 STRAIN
#
# The dataset contains cells from multiple P. falciparum strains.
# We isolate only NF54 strain cells for this analysis because:
#   1. NF54 is the strain used in Acharya et al. 2025 (the reference paper)
#   2. Mixing strains would confound expression analysis
#
# The STRAIN column in metadata was added in Section 6.
# =============================================================================

message("Subsetting to NF54 strain only...")

nf54_obj <- subset(
  seurat_obj,
  subset = STRAIN == "NF54"
)

# Free the full object from RAM — we only need NF54 cells
rm(sparse_mat, meta_df, seurat_obj)
gc()

message("NF54 subset: ", ncol(nf54_obj), " cells, ", nrow(nf54_obj), " genes")

# Check stage distribution
message("Cells per life cycle stage:")
print(table(nf54_obj$STAGE_LR))


# =============================================================================
# SECTION 8: NORMALIZE DATA
#
# Raw UMI counts cannot be directly compared between cells because:
#   - Some cells have more total UMIs (they were sequenced more deeply)
#   - A gene with 10 counts in a cell with 1000 total UMIs is MORE expressed
#     than 10 counts in a cell with 10000 total UMIs
#
# LogNormalize does three steps per cell:
#   Step 1: Divide each gene's count by the cell's total UMI count
#   Step 2: Multiply by scale.factor (10,000) — this puts values on a
#           comparable scale regardless of sequencing depth
#   Step 3: Take log1p (natural log + 1) — compresses the range and makes
#           the distribution more normal (log(0+1) = 0, so zeros stay zero)
#
# This creates the "data" layer in the Seurat object (separate from raw "counts").
# We only need this step — no ScaleData, PCA, or UMAP required because
# stages are already annotated in the metadata.
# =============================================================================

message("Normalizing data...")

nf54_obj <- NormalizeData(
  nf54_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

message("Normalization complete.")


# =============================================================================
# SECTION 9: DEFINE ALBA GENE IDs AND VERIFY THEY EXIST IN DATASET
#
# IMPORTANT: This dataset uses hyphens in gene IDs, not underscores.
#   Correct:   PF3D7-0813600
#   Incorrect: PF3D7_0813600  <- this will return "not found"
#
# We verify all 4 genes exist in rownames(nf54_obj) before plotting.
# If any are missing, valid_alba will only contain the ones found.
#
# PlasmoDB identifiers for PfAlba family:
#   PF3D7-0813600 = PfAlba1 (master translational regulator)
#   PF3D7-1006800 = PfAlba2 (forms complex with Alba1, represses var genes)
#   PF3D7-1418300 = PfAlba3 (endonuclease activity, strong var repressor)
#   PF3D7-0600200 = PfAlba4 (MAF1 domain, tolerates overexpression)
# =============================================================================

message("Defining PfAlba gene IDs...")

alba_genes <- c(
  "PF3D7-0813600",   # PfAlba1
  "PF3D7-1006800",   # PfAlba2
  "PF3D7-1418300",   # PfAlba3
  "PF3D7-0600200"    # PfAlba4
)

# Keep only genes that actually exist in our dataset
valid_alba <- alba_genes[
  alba_genes %in% rownames(nf54_obj)
]

message("Genes found in dataset: ", paste(valid_alba, collapse = ", "))

# Human-readable name mapping for plot labels
alba_labels <- c(
  "PF3D7-0813600" = "PfAlba1",
  "PF3D7-1006800" = "PfAlba2",
  "PF3D7-1418300" = "PfAlba3",
  "PF3D7-0600200" = "PfAlba4"
)


# =============================================================================
# SECTION 10: EXPRESSION STATISTICS — RAW COUNT CROSS-VALIDATION
#
# Before trusting the plots, we verify from the raw counts directly:
#   - How many total UMIs does each Alba gene accumulate?
#   - How many cells detect each gene at all?
#   - What is the maximum UMI count in any single cell?
#
# This is important because normalized plots can look misleading.
# If a gene shows "low" expression in plots, we verify it's genuinely low
# and not an artifact of normalization or scaling.
# =============================================================================

message("=== Raw Count Verification ===")

for (gene in valid_alba) {

  # Pull raw counts for this gene across all NF54 cells
  raw <- GetAssayData(
    nf54_obj,
    layer = "counts"
  )[gene, ]

  message(alba_labels[gene], " (", gene, ")")
  message("  Total UMIs across all cells : ", sum(raw))
  message("  Cells with >0 counts        : ", sum(raw > 0), " / ", ncol(nf54_obj))
  message("  Maximum UMI in single cell  : ", max(raw))
  message("------------------------------------------")
}


# =============================================================================
# SECTION 11: ALBA4 STAGE DISTRIBUTION CHECK
#
# PfAlba4 was detected in only 52 cells. We check if those cells cluster
# in a specific life cycle stage or are randomly distributed.
# Random distribution across stages = not stage-enriched = consistent with
# genuine low expression rather than stage-specific biology.
# =============================================================================

message("Checking which stages express PfAlba4...")

alba4_counts <- GetAssayData(
  nf54_obj,
  layer = "counts"
)["PF3D7-0600200", ]

stage_of_expressing <- nf54_obj$STAGE_LR[alba4_counts > 0]

message("PfAlba4-expressing cells per stage:")
print(table(stage_of_expressing))


# =============================================================================
# SECTION 12: VIOLIN PLOT
#
# VlnPlot shows the full distribution of normalized expression values
# for each gene, split by life cycle stage.
#
# pt.size = 0: hides individual cell dots. With 35k cells, showing
#   individual dots creates a completely black plot.
# ncol = 2: arranges the 4 gene plots in a 2x2 grid
#
# We rename each subplot title from the PlasmoDB ID to the gene name
# using a for loop that modifies each subplot individually.
# =============================================================================

message("Generating Violin Plot...")

vln <- VlnPlot(
  nf54_obj,
  features  = valid_alba,
  group.by  = "STAGE_LR",
  pt.size   = 0,
  ncol      = 2
)

# Replace PlasmoDB IDs with human-readable gene names as plot titles
for (i in seq_along(valid_alba)) {
  vln[[i]] <- vln[[i]] + ggtitle(alba_labels[valid_alba[i]])
}

ggsave(
  file.path(plots_dir, "alba_vlnplot.pdf"),
  vln,
  width  = 10,
  height = 8
)

message("Violin plot saved.")


# =============================================================================
# SECTION 13: DOT PLOT
#
# DotPlot is the preferred plot for stage-wise expression in publications.
# Each dot encodes TWO pieces of information simultaneously:
#   - Dot SIZE:  percentage of cells in that stage that express the gene
#   - Dot COLOR: average normalized expression level in expressing cells
#
# scale = FALSE: shows actual expression values (not z-scores).
#   With only 4 stage groups, z-score scaling (the default) produces
#   misleading results and triggers a warning message.
#
# scale_x_discrete(labels = alba_labels): replaces PlasmoDB IDs on the
#   x-axis with readable gene names (PfAlba1, PfAlba2, etc.)
# =============================================================================

message("Generating Dot Plot...")

dot <- DotPlot(
  nf54_obj,
  features = valid_alba,
  group.by = "STAGE_LR",
  scale    = FALSE           # Show real expression, not z-scores
) +
  scale_x_discrete(labels = alba_labels) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Alba Family Expression by Stage")

ggsave(
  file.path(plots_dir, "alba_dotplot_unscaled.pdf"),
  dot,
  width  = 8,
  height = 5
)

message("Dot plot saved.")


# =============================================================================
# SECTION 14: PREPARE DATA FOR BAR PLOTS
#
# VlnPlot and DotPlot work directly from the Seurat object.
# For a bar plot we need to build the data manually:
#
#   Step 1: Pull the normalized expression matrix for just the 4 Alba genes
#   Step 2: Transpose it (cells become rows, genes become columns)
#   Step 3: Add the STAGE_LR metadata as a column
#   Step 4: Convert from wide format to long format (pivot_longer)
#           Wide:  one column per gene
#           Long:  one row per gene-cell combination
#   Step 5: Calculate mean expression per gene per stage (group_by + summarise)
# =============================================================================

message("Preparing data for bar plots...")

# Pull normalized expression for just the 4 Alba genes
norm_data <- GetAssayData(
  nf54_obj,
  layer = "data"            # "data" = normalized, "counts" = raw
)[valid_alba, ]

# Build tidy data frame: rows = cells, columns = gene expression + stage
alba_df <- data.frame(
  t(as.matrix(norm_data)),
  Stage = nf54_obj$STAGE_LR
)

# Rename columns from PlasmoDB IDs to gene names
colnames(alba_df)[1:4] <- c("PfAlba1", "PfAlba2", "PfAlba3", "PfAlba4")

# Convert from wide to long format for ggplot2
alba_long <- pivot_longer(
  alba_df,
  cols      = c("PfAlba1", "PfAlba2", "PfAlba3", "PfAlba4"),
  names_to  = "Gene",
  values_to = "Expression"
)

# Calculate mean expression per gene per stage
alba_summary <- alba_long %>%
  group_by(Stage, Gene) %>%
  summarise(
    Mean_Expression = mean(Expression),
    .groups = "drop"
  )


# =============================================================================
# SECTION 15: BAR PLOT — ALL FOUR ALBA GENES
#
# Shows mean normalized expression of all 4 Alba genes side by side,
# grouped by life cycle stage.
#
# This plot makes it immediately obvious that Alba1 and Alba2 dominate,
# while Alba3 and Alba4 are nearly invisible on this scale.
# That is why we also generate a separate focused plot for Alba3/4.
# =============================================================================

message("Generating bar plot (all 4 Albas)...")

bar <- ggplot(
  alba_summary,
  aes(x = Stage, y = Mean_Expression, fill = Gene)
) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(
    title = "PfAlba Family Mean Expression by Life Cycle Stage",
    x     = "Life Cycle Stage",
    y     = "Mean Normalized Expression",
    fill  = "Gene"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(plots_dir, "alba_barplot.pdf"),
  bar,
  width  = 10,
  height = 6
)

message("Full bar plot saved.")


# =============================================================================
# SECTION 16: BAR PLOT — ALBA3 AND ALBA4 FOCUSED
#
# Because Alba3 and Alba4 expression is so low compared to Alba1/2,
# they are essentially invisible in the combined bar plot.
# This separate plot rescales the y-axis to show only Alba3 and Alba4,
# making their (low but non-zero) expression pattern visible.
# =============================================================================

message("Generating focused bar plot (Alba3 + Alba4 only)...")

alba_low <- alba_summary %>%
  filter(Gene %in% c("PfAlba3", "PfAlba4"))

bar_low <- ggplot(
  alba_low,
  aes(x = Stage, y = Mean_Expression, fill = Gene)
) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(
    title = "PfAlba3 and PfAlba4 Expression by Life Cycle Stage",
    x     = "Life Cycle Stage",
    y     = "Mean Normalized Expression",
    fill  = "Gene"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(plots_dir, "alba_barplot_low.pdf"),
  bar_low,
  width  = 8,
  height = 5
)

message("Focused bar plot saved.")


# =============================================================================
# PIPELINE COMPLETE
# =============================================================================

message("")
message("=== PIPELINE COMPLETE ===")
message("All plots saved to: ", normalizePath(plots_dir))
message("")
message("Output files:")
message("  - alba_vlnplot.pdf")
message("  - alba_dotplot_unscaled.pdf")
message("  - alba_barplot.pdf")
message("  - alba_barplot_low.pdf")
