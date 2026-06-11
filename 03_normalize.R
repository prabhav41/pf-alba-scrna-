# =============================================================================
# FILE:    scripts/01_load_data.R
# PURPOSE: Read raw CSV count matrix and convert to sparse matrix
#
# INPUT:   data/pf-ch10x-set4-ch10x-raw.csv
# OUTPUT:  sparse_mat (R object in memory, passed to 02_build_seurat.R)
#
# NEXT:    Run scripts/02_build_seurat.R
# =============================================================================

library(data.table)
library(Matrix)

raw_counts_file <- "data/pf-ch10x-set4-ch10x-raw.csv"

# --- Step 1: Read CSV using fread (fast, memory-efficient) ---
message("Reading raw counts... Do not interrupt.")
raw_df <- fread(raw_counts_file)

# --- Step 2: Extract gene names from column 1 ---
genes   <- raw_df[[1]]

# --- Step 3: Convert remaining columns to numeric matrix ---
raw_mat <- as.matrix(raw_df[, -1, with = FALSE])
rownames(raw_mat) <- genes

message("Matrix loaded: ", nrow(raw_mat), " genes x ", ncol(raw_mat), " cells")

# --- Step 4: Convert to sparse format (drops all zero values from RAM) ---
message("Converting to sparse matrix...")
sparse_mat <- Matrix(raw_mat, sparse = TRUE)

# --- Step 5: Free dense matrix from RAM ---
rm(raw_df, raw_mat)
gc()

message("Done. Object 'sparse_mat' ready for 02_build_seurat.R")
