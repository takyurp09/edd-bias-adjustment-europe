# =============================================================================
# 06_generate_conference_figures.R
# Additional figures for conference presentation and paper appendix
# =============================================================================

source("00_config.R")

library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(patchwork)
library(scales)

cat("\n========================================\n")
cat("GENERATING CONFERENCE FIGURES\n")
cat("========================================\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")
dskill <- read_parquet(file.path(OUT_TABLES, "dskill_paired.parquet"))
metrics_pooled <- read_parquet(file.path(OUT_TABLES, "metrics_pooled_paired.parquet"))

# Load shapefile for maps
gdf <- st_read(SHAPEFILE_PATH, quiet = TRUE)

cat("  Data loaded successfully\n\n")

# =============================================================================
# FIGURE 1: SIMPLIFIED ΔSKILL BY THRESHOLD
# Cleaner version than Fig_dskill_mae_by_threshold_crop.png
# =============================================================================

cat("FIGURE 1: Simplified ΔSkill by threshold...\n")

fig1_data <- dskill %>%
  group_by(source_id, edd_bin) %>%
  summarise(
    mean_dskill = mean(dskill_mae, na.rm = TRUE),
    se_dskill = sd(dskill_mae, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    threshold = as.numeric(gsub("edd_", "", edd_bin)),
    improves = mean_dskill > 0
  )

p1 <- ggplot(fig1_data, aes(x = threshold, y = mean_dskill, color = source_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_dskill - se_dskill, ymax = mean_dskill + se_dskill),
                width = 0.5, linewidth = 0.8) +
  scale_x_continuous(breaks = EDD_THRESHOLDS) +
  scale_color_manual(values = c("GFDL-ESM4" = "#E69F00",
                                 "MPI-ESM1-2-HR" = "#56B4E9",
                                 "EC-Earth3" = "#009E73")) +
  labs(
    title = "Bias-adjustment skill across temperature thresholds",
    subtitle = "Mean ΔSkill (MAE improvement) averaged across all regions, crops, and periods",
    x = "EDD threshold (°C)",
    y = "ΔSkill MAE (°C·days)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_FIGURES, "conference_fig1_dskill_simplified.png"),
       p1, width = 10, height = 6, dpi = FIG_DPI)

cat("  Saved: conference_fig1_dskill_simplified.png\n\n")


# =============================================================================
# FIGURE 2: THRESHOLD PERFORMANCE MATRIX
# Heatmap showing ΔSkill by model × threshold
# =============================================================================

cat("FIGURE 2: Threshold performance matrix...\n")

fig2_data <- fig1_data %>%
  mutate(
    model_label = factor(source_id, levels = PAIRED_MODELS),
    threshold_label = paste0("EDD-", threshold)
  )

p2 <- ggplot(fig2_data, aes(x = threshold_label, y = model_label, fill = mean_dskill)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(mean_dskill, 0)), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = COL_WORSENS, mid = "white", high = COL_IMPROVES,
    midpoint = 0,
    limits = c(-20, 150),
    name = "ΔSkill MAE\n(°C·days)"
  ) +
  labs(
    title = "Model performance matrix: Bias-adjustment skill by threshold",
    subtitle = "Positive (blue) = improvement; Negative (red) = degradation",
    x = "Temperature threshold",
    y = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(OUT_FIGURES, "conference_fig2_performance_matrix.png"),
       p2, width = 10, height = 5, dpi = FIG_DPI)

cat("  Saved: conference_fig2_performance_matrix.png\n\n")


# =============================================================================
# FIGURE 3: CROP SENSITIVITY RANKING
# Which crops benefit most from bias adjustment?
# =============================================================================

cat("FIGURE 3: Crop sensitivity ranking...\n")

fig3_data <- dskill %>%
  filter(edd_bin %in% c("edd_0", "edd_8")) %>%
  group_by(crop, edd_bin) %>%
  summarise(
    mean_dskill = mean(dskill_mae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    threshold = paste0("EDD-", gsub("edd_", "", edd_bin)),
    crop_label = tools::toTitleCase(gsub("_", " ", crop))
  )

p3 <- ggplot(fig3_data, aes(x = reorder(crop_label, mean_dskill), y = mean_dskill, fill = threshold)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("EDD-0" = "#4575b4", "EDD-8" = "#91bfdb")) +
  coord_flip() +
  labs(
    title = "Crop-specific benefits of bias adjustment",
    subtitle = "Mean ΔSkill (MAE improvement) at moderate thresholds",
    x = "Crop",
    y = "ΔSkill MAE (°C·days)",
    fill = "Threshold"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(OUT_FIGURES, "conference_fig3_crop_ranking.png"),
       p3, width = 10, height = 6, dpi = FIG_DPI)

cat("  Saved: conference_fig3_crop_ranking.png\n\n")


# =============================================================================
# FIGURE 4: ERROR DECOMPOSITION
# Bias vs MAE vs RMSE across thresholds
# =============================================================================

cat("FIGURE 4: Error decomposition...\n")

fig4_data <- metrics_pooled %>%
  filter(source_id == "GFDL-ESM4") %>%
  group_by(source_type, edd_bin) %>%
  summarise(
    mean_bias = mean(abs(bias), na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(threshold = as.numeric(gsub("edd_", "", edd_bin))) %>%
  pivot_longer(cols = c(mean_bias, mean_mae, mean_rmse),
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric_label = case_when(
      metric == "mean_bias" ~ "|Bias|",
      metric == "mean_mae" ~ "MAE",
      metric == "mean_rmse" ~ "RMSE"
    ),
    source_label = ifelse(source_type == "raw", "Raw CMIP6", "Bias-adjusted")
  )

p4 <- ggplot(fig4_data, aes(x = threshold, y = value, color = metric_label, linetype = source_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = EDD_THRESHOLDS) +
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("Raw CMIP6" = "dashed", "Bias-adjusted" = "solid")) +
  labs(
    title = "Error metrics across thresholds (GFDL-ESM4)",
    subtitle = "Comparison of raw vs bias-adjusted performance",
    x = "EDD threshold (°C)",
    y = "Error magnitude (°C·days)",
    color = "Metric",
    linetype = "Dataset"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_FIGURES, "conference_fig4_error_decomposition.png"),
       p4, width = 10, height = 6, dpi = FIG_DPI)

cat("  Saved: conference_fig4_error_decomposition.png\n\n")


# =============================================================================
# FIGURE 5: DISTRIBUTION COMPARISON
# Violin plots showing EDD distributions for raw vs bias-adj vs ERA5
# =============================================================================

cat("FIGURE 5: Distribution comparison (sample)...\n")

# Sample one region/crop/period for clarity
sample_region <- dskill %>%
  filter(edd_bin == "edd_8") %>%
  slice(1) %>%
  select(adm_code, crop, period, irrigation) %>%
  as.list()

# This would require going back to the spine data
# For now, create a placeholder note
cat("  NOTE: Distribution comparison requires spine-level data\n")
cat("  Skipping for now - can add if needed\n\n")


# =============================================================================
# FIGURE 6: SPATIAL COVERAGE MAP
# Show which ADM1 regions have data for each crop
# =============================================================================

cat("FIGURE 6: Spatial coverage map...\n")

coverage_data <- dskill %>%
  filter(edd_bin == "edd_8") %>%
  group_by(adm_code, crop) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  mutate(has_data = n_obs > 0) %>%
  group_by(adm_code) %>%
  summarise(n_crops = sum(has_data), .groups = "drop")

gdf_coverage <- gdf %>%
  left_join(coverage_data, by = c("GID_1" = "adm_code")) %>%
  mutate(n_crops = ifelse(is.na(n_crops), 0, n_crops))

p6 <- ggplot(gdf_coverage) +
  geom_sf(aes(fill = n_crops), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "plasma", name = "Number of\ncrops") +
  labs(
    title = "Spatial coverage: Crops analyzed by ADM1 region",
    subtitle = "Number of crops with harvested area data per region"
  ) +
  theme_void(base_size = 14) +
  theme(legend.position = "right")

ggsave(file.path(OUT_FIGURES, "conference_fig6_spatial_coverage.png"),
       p6, width = 12, height = 8, dpi = FIG_DPI)

cat("  Saved: conference_fig6_spatial_coverage.png\n\n")


# =============================================================================
# FIGURE 7: MODEL COMPARISON (ALL 7 MODELS)
# Line plot showing raw MAE for all available models
# =============================================================================

cat("FIGURE 7: Multi-model comparison...\n")

# This requires metrics_pooled_raw_all
metrics_raw_all <- read_parquet(file.path(OUT_TABLES, "metrics_pooled_raw_all.parquet"))

fig7_data <- metrics_raw_all %>%
  group_by(source_id, edd_bin) %>%
  summarise(
    mean_mae = mean(mae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(threshold = as.numeric(gsub("edd_", "", edd_bin)))

p7 <- ggplot(fig7_data, aes(x = threshold, y = mean_mae, color = source_id)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = EDD_THRESHOLDS) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Raw CMIP6 model comparison",
    subtitle = "Mean absolute error vs ERA5 across all models",
    x = "EDD threshold (°C)",
    y = "MAE (°C·days)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_FIGURES, "conference_fig7_model_comparison.png"),
       p7, width = 12, height = 6, dpi = FIG_DPI)

cat("  Saved: conference_fig7_model_comparison.png\n\n")


# =============================================================================
# FIGURE 8: IMPROVEMENT FREQUENCY
# Bar chart showing % of regions where bias adjustment improves performance
# =============================================================================

cat("FIGURE 8: Improvement frequency...\n")

fig8_data <- dskill %>%
  group_by(source_id, edd_bin) %>%
  summarise(
    pct_improved = mean(improved_mae, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    threshold = as.numeric(gsub("edd_", "", edd_bin)),
    threshold_label = paste0("EDD-", threshold)
  )

p8 <- ggplot(fig8_data, aes(x = threshold_label, y = pct_improved, fill = source_id)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("GFDL-ESM4" = "#E69F00",
                                "MPI-ESM1-2-HR" = "#56B4E9",
                                "EC-Earth3" = "#009E73")) +
  labs(
    title = "Frequency of bias-adjustment improvement",
    subtitle = "Percentage of ADM1-crop-period combinations where bias adjustment reduces MAE",
    x = "Temperature threshold",
    y = "% Improved",
    fill = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(OUT_FIGURES, "conference_fig8_improvement_frequency.png"),
       p8, width = 10, height = 6, dpi = FIG_DPI)

cat("  Saved: conference_fig8_improvement_frequency.png\n\n")


# =============================================================================
# DONE
# =============================================================================

cat("========================================\n")
cat("CONFERENCE FIGURES GENERATION COMPLETE\n")
cat("========================================\n\n")

cat("Figures saved in:", OUT_FIGURES, "\n")
cat("\nGenerated figures:\n")
cat("  1. conference_fig1_dskill_simplified.png\n")
cat("  2. conference_fig2_performance_matrix.png\n")
cat("  3. conference_fig3_crop_ranking.png\n")
cat("  4. conference_fig4_error_decomposition.png\n")
cat("  5. conference_fig6_spatial_coverage.png\n")
cat("  6. conference_fig7_model_comparison.png\n")
cat("  7. conference_fig8_improvement_frequency.png\n\n")

cat("These figures can be used in:\n")
cat("  - Conference presentation slides\n")
cat("  - Paper appendix\n")
cat("  - Supplementary materials\n\n")
