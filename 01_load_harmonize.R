# =============================================================================
# 01_load_harmonize.R
# Load ERA5, ESGF, and bias-adjusted parquets
# Parse source_id, tag source_type, output clean spine parquets
#
# Outputs:
#   spine/spine_paired.parquet   — 3 matched models, all 3 source types
#   spine/spine_raw_all.parquet  — all 7 esgf models vs ERA5 (model by model)
# =============================================================================

source("00_config.R")

library(arrow)
library(dplyr)
library(tidyr)
library(stringr)

# =============================================================================
# 1. LOAD RAW PARQUETS
# =============================================================================

cat("\n--- Loading parquets ---\n")

df_era5    <- read_parquet(DATA_PATHS$era5)
df_esgf    <- read_parquet(DATA_PATHS$esgf)
df_biasadj <- read_parquet(DATA_PATHS$biasadj)

cat("ERA5    :", nrow(df_era5),    "rows\n")
cat("ESGF    :", nrow(df_esgf),    "rows\n")
cat("Biasadj :", nrow(df_biasadj), "rows\n")

# =============================================================================
# 2. PARSE / ASSIGN source_id
# =============================================================================

cat("\n--- Parsing source_id ---\n")

df_era5 <- df_era5 %>%
  mutate(source_id   = "ERA5",
         source_type = "era5")

df_esgf <- df_esgf %>%
  mutate(source_id   = str_remove(run_id, ESGF_RUNID_SUFFIX),
         source_type = "raw")

df_biasadj <- df_biasadj %>%
  mutate(source_type = "biasadj")

cat("ESGF source_ids    :", paste(sort(unique(df_esgf$source_id)),    collapse=", "), "\n")
cat("Biasadj source_ids :", paste(sort(unique(df_biasadj$source_id)), collapse=", "), "\n")

# =============================================================================
# 3. VERIFY EXPECTED MODELS
# =============================================================================

cat("\n--- Verifying model coverage ---\n")

missing_paired_esgf <- setdiff(PAIRED_MODELS, unique(df_esgf$source_id))
missing_paired_ba   <- setdiff(PAIRED_MODELS, unique(df_biasadj$source_id))
missing_all_esgf    <- setdiff(ALL_ESGF_MODELS, unique(df_esgf$source_id))

if (length(missing_paired_esgf) > 0) warning("Missing paired models in ESGF    : ", paste(missing_paired_esgf, collapse=", ")) else cat("All paired models in ESGF ✓\n")
if (length(missing_paired_ba)   > 0) warning("Missing paired models in biasadj : ", paste(missing_paired_ba,   collapse=", ")) else cat("All paired models in biasadj ✓\n")
if (length(missing_all_esgf)    > 0) warning("Missing from ALL_ESGF_MODELS     : ", paste(missing_all_esgf,    collapse=", ")) else cat("All esgf models present ✓\n")

# =============================================================================
# 4. STANDARDIZE COLUMNS
# =============================================================================

cat("\n--- Standardizing columns ---\n")

KEEP_COLS <- c("source_id", "source_type",
               "adm_code", "adm_name", "poly_idx",
               "year", "period", "crop", "irrigation", "calendar",
               "days", EDD_BINS)

df_era5    <- df_era5    %>% select(all_of(KEEP_COLS))
df_esgf    <- df_esgf    %>% select(all_of(KEEP_COLS))
df_biasadj <- df_biasadj %>% select(all_of(KEEP_COLS))

# =============================================================================
# 5. PIVOT TO LONG FORMAT
# =============================================================================

cat("\n--- Pivoting to long format ---\n")

to_long <- function(df) {
  df %>%
    pivot_longer(cols = all_of(EDD_BINS),
                 names_to  = "edd_bin",
                 values_to = "edd_value") %>%
    mutate(edd_bin = factor(edd_bin, levels = EDD_BINS))
}

df_era5_long    <- to_long(df_era5);    rm(df_era5);    gc()
df_esgf_long    <- to_long(df_esgf);    rm(df_esgf);    gc()
df_biasadj_long <- to_long(df_biasadj); rm(df_biasadj); gc()

cat("ERA5 long    :", nrow(df_era5_long),    "rows\n")
cat("ESGF long    :", nrow(df_esgf_long),    "rows\n")
cat("Biasadj long :", nrow(df_biasadj_long), "rows\n")

# =============================================================================
# 6. BUILD SPINE: PAIRED
# Only 3 matched models — manageable size
# =============================================================================

cat("\n--- Building paired spine ---\n")

LONG_KEY_COLS <- c("adm_code", "poly_idx", "year", "period",
                   "crop", "irrigation", "calendar", "edd_bin")

esgf_paired    <- df_esgf_long    %>% filter(source_id %in% PAIRED_MODELS)
biasadj_paired <- df_biasadj_long %>% filter(source_id %in% PAIRED_MODELS)

# Find common keys
common_keys <- esgf_paired %>%
  distinct(across(all_of(c("source_id", LONG_KEY_COLS)))) %>%
  rename(paired_model = source_id) %>%
  inner_join(
    biasadj_paired %>%
      distinct(across(all_of(c("source_id", LONG_KEY_COLS)))) %>%
      rename(paired_model = source_id),
    by = c("paired_model", LONG_KEY_COLS)
  )

cat("Common key combinations (paired):", nrow(common_keys), "\n")

# Filter each source to common keys
esgf_paired <- esgf_paired %>%
  rename(paired_model = source_id) %>%
  semi_join(common_keys, by = c("paired_model", LONG_KEY_COLS)) %>%
  rename(source_id = paired_model)

biasadj_paired <- biasadj_paired %>%
  rename(paired_model = source_id) %>%
  semi_join(common_keys, by = c("paired_model", LONG_KEY_COLS)) %>%
  rename(source_id = paired_model)

# ERA5: replicate only for the 3 paired models (not 7) — much smaller
era5_paired <- bind_rows(
  lapply(PAIRED_MODELS, function(m) {
    df_era5_long %>%
      semi_join(
        common_keys %>% filter(paired_model == m),
        by = LONG_KEY_COLS
      ) %>%
      mutate(source_id = m)
  })
)

spine_paired <- bind_rows(esgf_paired, biasadj_paired, era5_paired)
rm(esgf_paired, biasadj_paired, era5_paired, common_keys); gc()

cat("Paired spine rows:", nrow(spine_paired), "\n")
cat("source_type breakdown:\n")
print(table(spine_paired$source_type))

write_parquet(spine_paired, file.path(OUT_SPINE, "spine_paired.parquet"))
cat("Saved spine_paired.parquet ✓\n")
rm(spine_paired); gc()

# =============================================================================
# 7. BUILD SPINE: RAW ALL
# Process ONE model at a time — never load all 7 × ERA5 into memory at once
# =============================================================================

cat("\n--- Building raw-all spine (model by model) ---\n")

spine_raw_all_list <- list()

for (model in ALL_ESGF_MODELS) {

  cat("  Processing:", model, "\n")

  model_esgf <- df_esgf_long %>% filter(source_id == model)

  # Find common keys with ERA5 for this model
  common_keys_m <- model_esgf %>%
    distinct(across(all_of(LONG_KEY_COLS))) %>%
    inner_join(
      df_era5_long %>% distinct(across(all_of(LONG_KEY_COLS))),
      by = LONG_KEY_COLS
    )

  model_esgf <- model_esgf %>%
    semi_join(common_keys_m, by = LONG_KEY_COLS)

  era5_m <- df_era5_long %>%
    semi_join(common_keys_m, by = LONG_KEY_COLS) %>%
    mutate(source_id = model)

  spine_raw_all_list[[model]] <- bind_rows(model_esgf, era5_m)
  rm(model_esgf, era5_m, common_keys_m); gc()

  cat("    Rows:", nrow(spine_raw_all_list[[model]]), "\n")
}

spine_raw_all <- bind_rows(spine_raw_all_list)
rm(spine_raw_all_list, df_esgf_long, df_biasadj_long, df_era5_long); gc()

cat("Raw-all spine rows:", nrow(spine_raw_all), "\n")
cat("source_id breakdown:\n")
print(table(spine_raw_all$source_id[spine_raw_all$source_type == "raw"]))

write_parquet(spine_raw_all, file.path(OUT_SPINE, "spine_raw_all.parquet"))
cat("Saved spine_raw_all.parquet ✓\n")

cat("\n--- 01_load_harmonize.R complete ---\n")
