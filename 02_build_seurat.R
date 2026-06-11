# =============================================================================
# FILE:    scripts/04_alba_analysis.R
# PURPOSE: Cross-validate expression from raw counts. Compute per-stage
#          mean expression for bar plots.
#
# REQUIRES: nf54_obj, valid_alba, alba_labels (from 03_normalize.R)
# OUTPUT:   alba_summary (data frame for bar plots)
#           Console output: UMI stats and stage distributions
#
# NEXT:     Run scripts/05_visualize.R
# =============================================================================

library(Seurat)
library(dplyr)
library(tidyr)

# --- Step 1: Raw count verification for each Alba gene ---
# Pulls directly from raw counts layer (before normalization)
# This is the ground truth check — if counts are near-zero here,
# the gene is genuinely lowly expressed, not a normalization artifact.
message("=== Raw Count Verification ===")

for (gene in valid_alba) {
  raw <- GetAssayData(nf54_obj, layer = "counts")[gene, ]

  message(alba_labels[gene], " (", gene, ")")
  message("  Total UMIs across all cells : ", sum(raw))
  message("  Cells with >0 counts        : ", sum(raw > 0), " / ", ncol(nf54_obj))
  message("  Max UMI in a single cell    : ", max(raw))
  message("---")
}

# --- Step 2: Stage distribution of cells ---
message("=== Cells per life cycle stage ===")
print(table(nf54_obj$STAGE_LR))

# --- Step 3: Check which stages express PfAlba4 ---
# Alba4 is detected in only 52 cells. Is it enriched in any stage?
# If scattered randomly, there is no stage preference.
alba4_counts <- GetAssayData(nf54_obj, layer = "counts")["PF3D7-0600200", ]
stage_of_expressing <- nf54_obj$STAGE_LR[alba4_counts > 0]
message("=== PfAlba4 expressing cells per stage ===")
print(table(stage_of_expressing))

# --- Step 4: Build tidy data frame from normalized expression ---
norm_data <- GetAssayData(nf54_obj, layer = "data")[valid_alba, ]

alba_df <- data.frame(
  t(as.matrix(norm_data)),
  Stage = nf54_obj$STAGE_LR
)
colnames(alba_df)[1:4] <- c("PfAlba1", "PfAlba2", "PfAlba3", "PfAlba4")

# --- Step 5: Convert to long format for ggplot2 ---
alba_long <- pivot_longer(
  alba_df,
  cols      = c("PfAlba1", "PfAlba2", "PfAlba3", "PfAlba4"),
  names_to  = "Gene",
  values_to = "Expression"
)

# --- Step 6: Calculate mean expression per gene per stage ---
alba_summary <- alba_long %>%
  group_by(Stage, Gene) %>%
  summarise(Mean_Expression = mean(Expression), .groups = "drop")

message("Done. Object 'alba_summary' ready for 05_visualize.R")
