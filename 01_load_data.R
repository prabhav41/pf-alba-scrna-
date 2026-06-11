# =============================================================================
# FILE:    scripts/05_visualize.R
# PURPOSE: Generate all expression plots and save to plots/ directory
#
# REQUIRES: nf54_obj, valid_alba, alba_labels, alba_summary
#           (from 03_normalize.R and 04_alba_analysis.R)
# OUTPUT:   4 PDF files in plots/ directory
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)

plots_dir <- "plots/"
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# =============================================================================
# PLOT 1: VIOLIN PLOT
# Shows distribution of normalized expression per stage per gene.
# Good for seeing the spread (are most cells zero? bimodal?)
# =============================================================================
message("Generating Violin Plot...")

vln <- VlnPlot(
  nf54_obj,
  features = valid_alba,
  group.by = "STAGE_LR",
  pt.size  = 0,              # No dots — 35k cells would black out the plot
  ncol     = 2
)

# Replace PlasmoDB ID titles with gene names
for (i in seq_along(valid_alba)) {
  vln[[i]] <- vln[[i]] + ggtitle(alba_labels[valid_alba[i]])
}

ggsave(file.path(plots_dir, "alba_vlnplot.pdf"), vln, width = 10, height = 8)
message("Saved: alba_vlnplot.pdf")

# =============================================================================
# PLOT 2: DOT PLOT (primary publication figure)
# Dot size = % cells expressing the gene
# Dot color = average normalized expression level
# scale = FALSE is critical — with 4 groups, z-score scaling is misleading
# =============================================================================
message("Generating Dot Plot...")

dot <- DotPlot(
  nf54_obj,
  features = valid_alba,
  group.by = "STAGE_LR",
  scale    = FALSE
) +
  scale_x_discrete(labels = alba_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Alba Family Expression by Stage")

ggsave(file.path(plots_dir, "alba_dotplot_unscaled.pdf"), dot, width = 8, height = 5)
message("Saved: alba_dotplot_unscaled.pdf")

# =============================================================================
# PLOT 3: BAR PLOT — ALL 4 ALBA GENES
# Mean normalized expression per stage, all genes side by side
# Alba1/2 dominate; Alba3/4 nearly invisible at this scale
# =============================================================================
message("Generating Bar Plot (all 4 Albas)...")

bar <- ggplot(alba_summary, aes(x = Stage, y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(
    title = "PfAlba Family Mean Expression by Life Cycle Stage",
    x     = "Life Cycle Stage",
    y     = "Mean Normalized Expression",
    fill  = "Gene"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "alba_barplot.pdf"), bar, width = 10, height = 6)
message("Saved: alba_barplot.pdf")

# =============================================================================
# PLOT 4: BAR PLOT — ALBA3 AND ALBA4 ONLY
# Focused plot with rescaled y-axis so Alba3/4 patterns are visible
# =============================================================================
message("Generating Bar Plot (Alba3 + Alba4 focused)...")

alba_low <- alba_summary %>% filter(Gene %in% c("PfAlba3", "PfAlba4"))

bar_low <- ggplot(alba_low, aes(x = Stage, y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(
    title = "PfAlba3 and PfAlba4 Expression by Life Cycle Stage",
    x     = "Life Cycle Stage",
    y     = "Mean Normalized Expression",
    fill  = "Gene"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "alba_barplot_low.pdf"), bar_low, width = 8, height = 5)
message("Saved: alba_barplot_low.pdf")

message("=== All plots saved to plots/ directory ===")
