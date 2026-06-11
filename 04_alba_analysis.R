# =============================================================================
# FILE:    scripts/03_normalize.R
# PURPOSE: Normalize count data and verify Alba gene names exist in dataset
#
# REQUIRES: nf54_obj (from 02_build_seurat.R)
# OUTPUT:   nf54_obj with normalized "data" layer added
#           valid_alba: vector of confirmed gene IDs
#           alba_labels: named vector mapping IDs to gene names
#
# NEXT:     Run scripts/04_alba_analysis.R
# =============================================================================

library(Seurat)

# --- Step 1: LogNormalize ---
# Divides by total UMIs per cell, multiplies by 10000, takes log1p
# This makes expression values comparable across cells of different depths
message("Normalizing data...")
nf54_obj <- NormalizeData(
  nf54_obj,
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
message("Normalization complete.")

# --- Step 2: Define PfAlba gene IDs ---
# IMPORTANT: This dataset uses HYPHENS not underscores
# Wrong: PF3D7_0813600  |  Correct: PF3D7-0813600
alba_genes <- c(
  "PF3D7-0813600",   # PfAlba1
  "PF3D7-1006800",   # PfAlba2
  "PF3D7-1418300",   # PfAlba3
  "PF3D7-0600200"    # PfAlba4
)

# --- Step 3: Verify genes exist in dataset (prevents silent errors) ---
valid_alba <- alba_genes[alba_genes %in% rownames(nf54_obj)]
message("Verified Alba genes in dataset: ", paste(valid_alba, collapse = ", "))

if (length(valid_alba) == 0) {
  stop("No Alba genes found. Check gene ID format with: head(rownames(nf54_obj), 5)")
}

# --- Step 4: Create human-readable name mapping ---
alba_labels <- c(
  "PF3D7-0813600" = "PfAlba1",
  "PF3D7-1006800" = "PfAlba2",
  "PF3D7-1418300" = "PfAlba3",
  "PF3D7-0600200" = "PfAlba4"
)

message("Done. Objects 'valid_alba' and 'alba_labels' ready for 04_alba_analysis.R")
