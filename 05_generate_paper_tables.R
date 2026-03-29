# =============================================================================
# 05_generate_paper_tables.R
# Extract exact summary statistics for paper tables and LaTeX
# =============================================================================

source("00_config.R")

library(arrow)
library(dplyr)
library(tidyr)
library(knitr)
library(xtable)

cat("\n========================================\n")
cat("GENERATING PAPER TABLES\n")
cat("========================================\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")
dskill <- read_parquet(file.path(OUT_TABLES, "dskill_paired.parquet"))
metrics_pooled <- read_parquet(file.path(OUT_TABLES, "metrics_pooled_paired.parquet"))
metrics_raw_all <- read_parquet(file.path(OUT_TABLES, "metrics_pooled_raw_all.parquet"))

cat("  dskill rows:", nrow(dskill), "\n")
cat("  metrics_pooled rows:", nrow(metrics_pooled), "\n")
cat("  metrics_raw_all rows:", nrow(metrics_raw_all), "\n\n")

# =============================================================================
# TABLE 1: SUMMARY ΔSKILL BY THRESHOLD × MODEL
# Mean across all ADM1 regions, crops, irrigation types, and periods
# =============================================================================

cat("TABLE 1: Mean ΔSkill by threshold × model...\n")

table1 <- dskill %>%
  group_by(source_id, edd_bin) %>%
  summarise(
    mean_dskill_mae = mean(dskill_mae, na.rm = TRUE),
    median_dskill_mae = median(dskill_mae, na.rm = TRUE),
    sd_dskill_mae = sd(dskill_mae, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(source_id, edd_bin)

write.csv(table1, file.path(OUT_TABLES, "paper_table1_dskill_summary.csv"), row.names = FALSE)
cat("  Saved: paper_table1_dskill_summary.csv\n")

# Pivot for LaTeX (wide format: models as rows, thresholds as columns)
table1_wide <- table1 %>%
  select(source_id, edd_bin, mean_dskill_mae) %>%
  pivot_wider(names_from = edd_bin, values_from = mean_dskill_mae) %>%
  arrange(source_id)

write.csv(table1_wide, file.path(OUT_TABLES, "paper_table1_wide.csv"), row.names = FALSE)
cat("  Saved: paper_table1_wide.csv\n\n")


# =============================================================================
# TABLE 2: RAW VS BIAS-ADJUSTED MAE COMPARISON
# =============================================================================

cat("TABLE 2: Raw vs Bias-adjusted MAE...\n")

table2 <- metrics_pooled %>%
  filter(source_id %in% PAIRED_MODELS) %>%
  group_by(source_id, source_type, edd_bin) %>%
  summarise(
    mean_mae = mean(mae, na.rm = TRUE),
    median_mae = median(mae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = source_type, values_from = c(mean_mae, median_mae)) %>%
  arrange(source_id, edd_bin)

write.csv(table2, file.path(OUT_TABLES, "paper_table2_mae_comparison.csv"), row.names = FALSE)
cat("  Saved: paper_table2_mae_comparison.csv\n\n")


# =============================================================================
# TABLE 3: QUANTILE ERRORS (Q90, Q95, Q99)
# =============================================================================

cat("TABLE 3: Quantile errors...\n")

table3 <- metrics_pooled %>%
  filter(source_id %in% PAIRED_MODELS) %>%
  group_by(source_id, source_type, edd_bin) %>%
  summarise(
    mean_q90_err = mean(q90_err, na.rm = TRUE),
    mean_q95_err = mean(q95_err, na.rm = TRUE),
    mean_q99_err = mean(q99_err, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source_id, source_type, edd_bin)

write.csv(table3, file.path(OUT_TABLES, "paper_table3_quantile_errors.csv"), row.names = FALSE)
cat("  Saved: paper_table3_quantile_errors.csv\n\n")


# =============================================================================
# TABLE 4: TAIL MEAN ERRORS
# =============================================================================

cat("TABLE 4: Tail mean errors...\n")

table4 <- metrics_pooled %>%
  filter(source_id %in% PAIRED_MODELS) %>%
  group_by(source_id, source_type, edd_bin) %>%
  summarise(
    mean_tail_mean_err = mean(tail_mean_err, na.rm = TRUE),
    median_tail_mean_err = median(tail_mean_err, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source_id, source_type, edd_bin)

write.csv(table4, file.path(OUT_TABLES, "paper_table4_tail_errors.csv"), row.names = FALSE)
cat("  Saved: paper_table4_tail_errors.csv\n\n")


# =============================================================================
# TABLE 5: KS STATISTICS
# =============================================================================

cat("TABLE 5: KS statistics...\n")

table5 <- metrics_pooled %>%
  filter(source_id %in% PAIRED_MODELS) %>%
  group_by(source_id, source_type, edd_bin) %>%
  summarise(
    mean_ks_stat = mean(ks_stat, na.rm = TRUE),
    median_ks_stat = median(ks_stat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source_id, source_type, edd_bin)

write.csv(table5, file.path(OUT_TABLES, "paper_table5_ks_stats.csv"), row.names = FALSE)
cat("  Saved: paper_table5_ks_stats.csv\n\n")


# =============================================================================
# TABLE 6: CROP-SPECIFIC ΔSKILL AT KEY THRESHOLDS
# =============================================================================

cat("TABLE 6: Crop-specific ΔSkill at key thresholds...\n")

table6 <- dskill %>%
  filter(edd_bin %in% c("edd_0", "edd_8", "edd_28")) %>%
  group_by(source_id, crop, edd_bin) %>%
  summarise(
    mean_dskill_mae = mean(dskill_mae, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = edd_bin, values_from = mean_dskill_mae) %>%
  arrange(source_id, crop)

write.csv(table6, file.path(OUT_TABLES, "paper_table6_crop_dskill.csv"), row.names = FALSE)
cat("  Saved: paper_table6_crop_dskill.csv\n\n")


# =============================================================================
# TABLE 7: PERIOD-SPECIFIC ΔSKILL
# =============================================================================

cat("TABLE 7: Period-specific ΔSkill...\n")

table7 <- dskill %>%
  group_by(source_id, period, edd_bin) %>%
  summarise(
    mean_dskill_mae = mean(dskill_mae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source_id, period, edd_bin)

write.csv(table7, file.path(OUT_TABLES, "paper_table7_period_dskill.csv"), row.names = FALSE)
cat("  Saved: paper_table7_period_dskill.csv\n\n")


# =============================================================================
# TABLE 8: IRRIGATION-SPECIFIC ΔSKILL
# =============================================================================

cat("TABLE 8: Irrigation-specific ΔSkill...\n")

table8 <- dskill %>%
  group_by(source_id, irrigation, edd_bin) %>%
  summarise(
    mean_dskill_mae = mean(dskill_mae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source_id, irrigation, edd_bin)

write.csv(table8, file.path(OUT_TABLES, "paper_table8_irrigation_dskill.csv"), row.names = FALSE)
cat("  Saved: paper_table8_irrigation_dskill.csv\n\n")


# =============================================================================
# TABLE 9: RAW MODEL COMPARISON (ALL 7 MODELS)
# =============================================================================

cat("TABLE 9: Raw model MAE comparison (all models)...\n")

table9 <- metrics_raw_all %>%
  group_by(source_id, edd_bin) %>%
  summarise(
    mean_mae = mean(mae, na.rm = TRUE),
    median_mae = median(mae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source_id, edd_bin)

write.csv(table9, file.path(OUT_TABLES, "paper_table9_raw_all_mae.csv"), row.names = FALSE)
cat("  Saved: paper_table9_raw_all_mae.csv\n\n")


# =============================================================================
# TABLE 10: CORRELATION COMPARISON
# =============================================================================

cat("TABLE 10: Correlation comparison...\n")

# Extract correlation data from dskill table
table10 <- dskill %>%
  group_by(source_id, edd_bin) %>%
  summarise(
    mean_corr_raw = mean(corr_raw, na.rm = TRUE),
    mean_corr_ba = mean(corr_ba, na.rm = TRUE),
    mean_dskill_corr = mean(dskill_corr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source_id, edd_bin)

write.csv(table10, file.path(OUT_TABLES, "paper_table10_correlation.csv"), row.names = FALSE)
cat("  Saved: paper_table10_correlation.csv\n\n")


# =============================================================================
# SUMMARY STATS FOR ABSTRACT/TEXT
# =============================================================================

cat("SUMMARY STATS: Key numbers for abstract and text...\n\n")

summary_stats <- list()

# Overall ΔSkill range at moderate thresholds
summary_stats$dskill_edd0_range <- dskill %>%
  filter(edd_bin == "edd_0") %>%
  summarise(
    min = round(min(dskill_mae, na.rm = TRUE), 0),
    max = round(max(dskill_mae, na.rm = TRUE), 0),
    mean = round(mean(dskill_mae, na.rm = TRUE), 0)
  )

summary_stats$dskill_edd8_range <- dskill %>%
  filter(edd_bin == "edd_8") %>%
  summarise(
    min = round(min(dskill_mae, na.rm = TRUE), 0),
    max = round(max(dskill_mae, na.rm = TRUE), 0),
    mean = round(mean(dskill_mae, na.rm = TRUE), 0)
  )

summary_stats$dskill_edd28_range <- dskill %>%
  filter(edd_bin == "edd_28") %>%
  summarise(
    min = round(min(dskill_mae, na.rm = TRUE), 1),
    max = round(max(dskill_mae, na.rm = TRUE), 1),
    mean = round(mean(dskill_mae, na.rm = TRUE), 1)
  )

# Raw MAE baseline
summary_stats$raw_mae_baseline <- metrics_pooled %>%
  filter(source_type == "raw", edd_bin == "edd_0") %>%
  summarise(
    mean = round(mean(mae, na.rm = TRUE), 0),
    min = round(min(mae, na.rm = TRUE), 0),
    max = round(max(mae, na.rm = TRUE), 0)
  )

# Model-specific ΔSkill at key thresholds
summary_stats$dskill_by_model <- dskill %>%
  filter(edd_bin %in% c("edd_0", "edd_8", "edd_28")) %>%
  group_by(source_id, edd_bin) %>%
  summarise(mean_dskill = round(mean(dskill_mae, na.rm = TRUE), 0), .groups = "drop")

# Quantile error improvements
summary_stats$quantile_improvements <- metrics_pooled %>%
  filter(source_id == "GFDL-ESM4", edd_bin == "edd_0") %>%
  group_by(source_type) %>%
  summarise(
    mean_q90 = round(mean(q90_err, na.rm = TRUE), 0),
    mean_q95 = round(mean(q95_err, na.rm = TRUE), 0),
    .groups = "drop"
  )

# Tail error magnitudes
summary_stats$tail_error_magnitudes <- metrics_pooled %>%
  filter(edd_bin == "edd_0") %>%
  group_by(source_id, source_type) %>%
  summarise(
    mean_tail_err = round(mean(tail_mean_err, na.rm = TRUE), 0),
    .groups = "drop"
  )

# KS statistic improvements
summary_stats$ks_improvements <- metrics_pooled %>%
  filter(edd_bin %in% c("edd_0", "edd_28")) %>%
  group_by(source_id, source_type, edd_bin) %>%
  summarise(
    mean_ks = round(mean(ks_stat, na.rm = TRUE), 2),
    .groups = "drop"
  )

# Save summary stats as R object
saveRDS(summary_stats, file.path(OUT_TABLES, "summary_stats_for_text.rds"))
cat("  Saved: summary_stats_for_text.rds\n")

# Also save as readable text file
sink(file.path(OUT_TABLES, "summary_stats_for_text.txt"))
cat("SUMMARY STATISTICS FOR PAPER TEXT\n")
cat("==================================\n\n")

cat("ΔSkill at EDD-0 (moderate threshold):\n")
print(summary_stats$dskill_edd0_range)
cat("\n")

cat("ΔSkill at EDD-8 (moderate threshold):\n")
print(summary_stats$dskill_edd8_range)
cat("\n")

cat("ΔSkill at EDD-28 (extreme threshold):\n")
print(summary_stats$dskill_edd28_range)
cat("\n")

cat("Raw MAE baseline at EDD-0:\n")
print(summary_stats$raw_mae_baseline)
cat("\n")

cat("ΔSkill by model and threshold:\n")
print(summary_stats$dskill_by_model)
cat("\n")

cat("Quantile errors (GFDL-ESM4, EDD-0):\n")
print(summary_stats$quantile_improvements)
cat("\n")

cat("Tail mean errors:\n")
print(summary_stats$tail_error_magnitudes)
cat("\n")

cat("KS statistics:\n")
print(summary_stats$ks_improvements)
cat("\n")

sink()
cat("  Saved: summary_stats_for_text.txt\n\n")


# =============================================================================
# LATEX TABLE GENERATION
# =============================================================================

cat("Generating LaTeX tables...\n")

# LaTeX Table 1: Summary ΔSkill (wide format)
latex_table1 <- table1_wide %>%
  mutate(across(where(is.numeric), ~ round(.x, 0)))

xtab1 <- xtable(latex_table1,
                caption = "Mean $\\Delta$Skill (MAE improvement, °C$\\cdot$days) by model and EDD threshold, averaged across all ADM1 regions, crops, irrigation types, and periods.",
                label = "tab:summary_dskill")

print(xtab1,
      file = file.path(OUT_TABLES, "latex_table1_dskill.tex"),
      include.rownames = FALSE,
      sanitize.text.function = identity,
      booktabs = TRUE)

cat("  Saved: latex_table1_dskill.tex\n")

# LaTeX Table 2: MAE comparison (selected thresholds)
latex_table2 <- table2 %>%
  filter(edd_bin %in% c("edd_0", "edd_8", "edd_12", "edd_28")) %>%
  select(source_id, edd_bin, mean_mae_raw, mean_mae_biasadj) %>%
  mutate(across(where(is.numeric), ~ round(.x, 0))) %>%
  arrange(source_id, edd_bin)

xtab2 <- xtable(latex_table2,
                caption = "Raw vs.\ bias-adjusted mean absolute error (MAE, °C$\\cdot$days) at selected EDD thresholds.",
                label = "tab:mae_comparison")

print(xtab2,
      file = file.path(OUT_TABLES, "latex_table2_mae.tex"),
      include.rownames = FALSE,
      sanitize.text.function = identity,
      booktabs = TRUE)

cat("  Saved: latex_table2_mae.tex\n")

# LaTeX Table 3: Crop-specific ΔSkill (GFDL-ESM4 only, for main text)
latex_table3 <- table6 %>%
  filter(source_id == "GFDL-ESM4") %>%
  select(crop, edd_0, edd_8, edd_28) %>%
  mutate(across(where(is.numeric), ~ round(.x, 0)))

xtab3 <- xtable(latex_table3,
                caption = "Crop-specific $\\Delta$Skill (MAE improvement, °C$\\cdot$days) for GFDL-ESM4 at key thresholds.",
                label = "tab:crop_dskill")

print(xtab3,
      file = file.path(OUT_TABLES, "latex_table3_crop.tex"),
      include.rownames = FALSE,
      sanitize.text.function = identity,
      booktabs = TRUE)

cat("  Saved: latex_table3_crop.tex\n\n")


# =============================================================================
# DONE
# =============================================================================

cat("========================================\n")
cat("PAPER TABLES GENERATION COMPLETE\n")
cat("========================================\n\n")

cat("Tables saved in:", OUT_TABLES, "\n")
cat("\nNext steps:\n")
cat("  1. Review CSV files for accuracy\n")
cat("  2. Use LaTeX .tex files in paper\n")
cat("  3. Reference summary_stats_for_text.txt for exact numbers\n\n")
