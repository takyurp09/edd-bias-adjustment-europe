# =============================================================================
# 02_metrics.R
# Compute all evaluation metrics
# NOTE: adm_name is now in GROUP_COLS_POOLED for true ADM1-level metrics
# =============================================================================

source("00_config.R")

library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# =============================================================================
# 0. HELPERS
# =============================================================================

compute_pooled_metrics <- function(df, model_val_col, era5_val_col = "era5_value") {

  df2 <- df %>%
    mutate(
      err     = .data[[model_val_col]] - .data[[era5_val_col]],
      abs_err = abs(err),
      sq_err  = err^2
    )

  metrics <- df2 %>%
    group_by(across(all_of(GROUP_COLS_POOLED))) %>%
    summarise(
      n    = n(),
      bias = mean(err,     na.rm = TRUE),
      mae  = mean(abs_err, na.rm = TRUE),
      rmse = sqrt(mean(sq_err, na.rm = TRUE)),

      q90_err = quantile(.data[[model_val_col]], 0.90, na.rm=TRUE) -
                quantile(.data[[era5_val_col]],  0.90, na.rm=TRUE),
      q95_err = quantile(.data[[model_val_col]], 0.95, na.rm=TRUE) -
                quantile(.data[[era5_val_col]],  0.95, na.rm=TRUE),
      q99_err = quantile(.data[[model_val_col]], 0.99, na.rm=TRUE) -
                quantile(.data[[era5_val_col]],  0.99, na.rm=TRUE),

      era5_p95          = quantile(.data[[era5_val_col]], 0.95, na.rm=TRUE),
      exceed_freq_era5  = mean(.data[[era5_val_col]]  > era5_p95, na.rm=TRUE),
      exceed_freq_model = mean(.data[[model_val_col]] > era5_p95, na.rm=TRUE),
      exceed_freq_err   = exceed_freq_model - exceed_freq_era5,

      tail_thresh     = quantile(.data[[era5_val_col]], TAIL_QUANTILE, na.rm=TRUE),
      tail_mean_era5  = mean(.data[[era5_val_col]][
                          .data[[era5_val_col]] > tail_thresh], na.rm=TRUE),
      tail_mean_model = mean(.data[[model_val_col]][
                          .data[[era5_val_col]] > tail_thresh], na.rm=TRUE),
      tail_mean_err   = tail_mean_model - tail_mean_era5,

      .groups = "drop"
    ) %>%
    select(-era5_p95, -exceed_freq_era5, -exceed_freq_model,
           -tail_thresh, -tail_mean_era5, -tail_mean_model)

  # KS statistic via group_map
  ks_df <- df2 %>%
    group_by(across(all_of(GROUP_COLS_POOLED))) %>%
    group_map(function(grp, keys) {
      ks_val <- tryCatch(
        ks.test(grp[[model_val_col]], grp[[era5_val_col]])$statistic,
        error = function(e) NA_real_
      )
      bind_cols(keys, tibble(ks_stat = as.numeric(ks_val)))
    }, .keep = FALSE) %>%
    bind_rows()

  metrics %>% left_join(ks_df, by = GROUP_COLS_POOLED)
}


compute_corr_metrics <- function(df, model_val_col, era5_val_col = "era5_value") {

  annual <- df %>%
    group_by(across(all_of(c(GROUP_COLS_POOLED, "year")))) %>%
    summarise(
      model_annual = mean(.data[[model_val_col]], na.rm = TRUE),
      era5_annual  = mean(.data[[era5_val_col]],  na.rm = TRUE),
      .groups = "drop"
    )

  annual %>%
    group_by(across(all_of(GROUP_COLS_POOLED))) %>%
    summarise(
      n_years = n(),
      corr    = tryCatch(
        cor(model_annual, era5_annual, use = "complete.obs"),
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
}


# =============================================================================
# 1. PAIRED ANALYSIS
# =============================================================================

cat("\n=== PAIRED ANALYSIS ===\n")

spine_paired <- read_parquet(file.path(OUT_SPINE, "spine_paired.parquet"))
cat("Loaded spine_paired:", nrow(spine_paired), "rows\n")

era5_ref <- spine_paired %>%
  filter(source_type == "era5") %>%
  select(source_id, adm_code, adm_name, poly_idx, year,
         period, crop, irrigation, calendar, edd_bin,
         era5_value = edd_value)

# JOIN_COLS now includes adm_name for ADM1-level join
JOIN_COLS <- c("source_id", "adm_code", "adm_name", "poly_idx", "year",
               "period", "crop", "irrigation", "calendar", "edd_bin")

# --- Raw ---
cat("\nComputing pooled metrics: raw (paired)...\n")
raw_paired <- spine_paired %>%
  filter(source_type == "raw") %>%
  inner_join(era5_ref, by = JOIN_COLS)

metrics_pooled_raw <- compute_pooled_metrics(raw_paired, "edd_value") %>%
  mutate(source_type = "raw")
gc()

cat("Computing correlation metrics: raw (paired)...\n")
metrics_corr_raw <- compute_corr_metrics(raw_paired, "edd_value") %>%
  mutate(source_type = "raw")
rm(raw_paired); gc()

# --- Biasadj ---
cat("\nComputing pooled metrics: biasadj (paired)...\n")
biasadj_paired <- spine_paired %>%
  filter(source_type == "biasadj") %>%
  inner_join(era5_ref, by = JOIN_COLS)

metrics_pooled_biasadj <- compute_pooled_metrics(biasadj_paired, "edd_value") %>%
  mutate(source_type = "biasadj")
gc()

cat("Computing correlation metrics: biasadj (paired)...\n")
metrics_corr_biasadj <- compute_corr_metrics(biasadj_paired, "edd_value") %>%
  mutate(source_type = "biasadj")
rm(biasadj_paired, era5_ref, spine_paired); gc()

# --- Stack ---
metrics_pooled_paired <- bind_rows(metrics_pooled_raw, metrics_pooled_biasadj)
metrics_corr_paired   <- bind_rows(metrics_corr_raw,   metrics_corr_biasadj)

cat("Pooled metrics (paired) rows:", nrow(metrics_pooled_paired), "\n")
cat("Corr   metrics (paired) rows:", nrow(metrics_corr_paired),   "\n")

# --- ΔSkill ---
cat("\nComputing ΔSkill...\n")

METRIC_COLS <- c("bias","mae","rmse",
                 "q90_err","q95_err","q99_err",
                 "exceed_freq_err","tail_mean_err","ks_stat")

raw_wide <- metrics_pooled_raw %>%
  select(all_of(c(GROUP_COLS_POOLED, METRIC_COLS))) %>%
  rename_with(~ paste0("raw_", .x), all_of(METRIC_COLS))

ba_wide <- metrics_pooled_biasadj %>%
  select(all_of(c(GROUP_COLS_POOLED, METRIC_COLS))) %>%
  rename_with(~ paste0("ba_", .x), all_of(METRIC_COLS))

dskill <- raw_wide %>%
  inner_join(ba_wide, by = GROUP_COLS_POOLED) %>%
  mutate(
    dskill_bias        = raw_bias        - ba_bias,
    dskill_mae         = raw_mae         - ba_mae,
    dskill_rmse        = raw_rmse        - ba_rmse,
    dskill_q90         = abs(raw_q90_err)         - abs(ba_q90_err),
    dskill_q95         = abs(raw_q95_err)         - abs(ba_q95_err),
    dskill_q99         = abs(raw_q99_err)         - abs(ba_q99_err),
    dskill_exceed_freq = abs(raw_exceed_freq_err)  - abs(ba_exceed_freq_err),
    dskill_tail_mean   = abs(raw_tail_mean_err)    - abs(ba_tail_mean_err),
    dskill_ks          = raw_ks_stat      - ba_ks_stat,
    improved_mae   = dskill_mae  > 0,
    improved_rmse  = dskill_rmse > 0,
    improved_q95   = dskill_q95  > 0,
    improved_tail  = dskill_tail_mean > 0
  )

corr_raw <- metrics_corr_raw %>%
  select(all_of(c(GROUP_COLS_POOLED, "corr"))) %>%
  rename(corr_raw = corr)

corr_ba <- metrics_corr_biasadj %>%
  select(all_of(c(GROUP_COLS_POOLED, "corr"))) %>%
  rename(corr_ba = corr)

dskill <- dskill %>%
  left_join(corr_raw, by = GROUP_COLS_POOLED) %>%
  left_join(corr_ba,  by = GROUP_COLS_POOLED) %>%
  mutate(
    dskill_corr   = corr_ba - corr_raw,
    improved_corr = dskill_corr > 0
  )

cat("ΔSkill table rows:", nrow(dskill), "\n")

cat("\nSaving paired metric tables...\n")
write_parquet(metrics_pooled_paired, file.path(OUT_TABLES, "metrics_pooled_paired.parquet"))
write_parquet(metrics_corr_paired,   file.path(OUT_TABLES, "metrics_corr_paired.parquet"))
write_parquet(dskill,                file.path(OUT_TABLES, "dskill_paired.parquet"))
cat("Saved metrics_pooled_paired.parquet ✓\n")
cat("Saved metrics_corr_paired.parquet   ✓\n")
cat("Saved dskill_paired.parquet         ✓\n")

rm(metrics_pooled_raw, metrics_pooled_biasadj,
   metrics_corr_raw, metrics_corr_biasadj,
   raw_wide, ba_wide, corr_raw, corr_ba); gc()


# =============================================================================
# 2. RAW-ALL ANALYSIS
# =============================================================================

cat("\n=== RAW-ALL ANALYSIS ===\n")

spine_raw_all <- read_parquet(file.path(OUT_SPINE, "spine_raw_all.parquet"))
cat("Loaded spine_raw_all:", nrow(spine_raw_all), "rows\n")

era5_ref_all <- spine_raw_all %>%
  filter(source_type == "era5") %>%
  select(source_id, adm_code, adm_name, poly_idx, year,
         period, crop, irrigation, calendar, edd_bin,
         era5_value = edd_value)

metrics_pooled_list <- list()
metrics_corr_list   <- list()

for (model in ALL_ESGF_MODELS) {

  cat("\nProcessing raw model:", model, "\n")

  model_df <- spine_raw_all %>%
    filter(source_type == "raw", source_id == model) %>%
    inner_join(
      era5_ref_all %>% filter(source_id == model),
      by = JOIN_COLS
    )

  if (nrow(model_df) == 0) {
    cat("  No data found for", model, "— skipping\n")
    next
  }

  cat("  Rows after join:", nrow(model_df), "\n")

  metrics_pooled_list[[model]] <- compute_pooled_metrics(model_df, "edd_value") %>%
    mutate(source_type = "raw")

  metrics_corr_list[[model]] <- compute_corr_metrics(model_df, "edd_value") %>%
    mutate(source_type = "raw")

  rm(model_df); gc()
}

metrics_pooled_raw_all <- bind_rows(metrics_pooled_list)
metrics_corr_raw_all   <- bind_rows(metrics_corr_list)

cat("\nRaw-all pooled metrics rows:", nrow(metrics_pooled_raw_all), "\n")
cat("Raw-all corr   metrics rows:", nrow(metrics_corr_raw_all),   "\n")

cat("\nSaving raw-all metric tables...\n")
write_parquet(metrics_pooled_raw_all, file.path(OUT_TABLES, "metrics_pooled_raw_all.parquet"))
write_parquet(metrics_corr_raw_all,   file.path(OUT_TABLES, "metrics_corr_raw_all.parquet"))
cat("Saved metrics_pooled_raw_all.parquet ✓\n")
cat("Saved metrics_corr_raw_all.parquet   ✓\n")

rm(spine_raw_all, era5_ref_all,
   metrics_pooled_list, metrics_corr_list); gc()

cat("\n--- 02_metrics.R complete ---\n")
