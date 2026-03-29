# =============================================================================
# 00_config.R
# Central configuration for EDD bias-adjustment evaluation pipeline
# =============================================================================

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------

BASE_PATH <- "/Users/tahmid/Library/CloudStorage/GoogleDrive-tahmid@udel.edu/Other computers/My Laptop/UDel/Taky_research/ci26_biasadj_cmip6/code"

DATA_PATHS <- list(
  era5    = file.path(BASE_PATH, "era5/src/output/europe/europe_edd_adm1_seasonal_panel.parquet"),
  esgf    = file.path(BASE_PATH, "esgf/src/output/europe/europe_edd_adm1_seasonal_panel.parquet"),
  biasadj = file.path(BASE_PATH, "getedd_cil-gdpcir/src/output/europe/europe_edd_adm1_seasonal_panel.parquet")
)

SHAPEFILE_PATH <- file.path(BASE_PATH, "../data/shapefiles/gadm41_EUR_shp/gadm41_EUR_1.shp")

OUT_DIR    <- file.path(BASE_PATH, "../outputs")
OUT_TABLES <- file.path(OUT_DIR, "tables")
OUT_FIGURES<- file.path(OUT_DIR, "figures")
OUT_SPINE  <- file.path(OUT_DIR, "spine")

for (d in c(OUT_TABLES, OUT_FIGURES, OUT_SPINE)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------------------------------------------------------
# MODELS
# -----------------------------------------------------------------------------

PAIRED_MODELS <- c("GFDL-ESM4", "MPI-ESM1-2-HR", "EC-Earth3")

ALL_ESGF_MODELS <- c(
  "GFDL-ESM4", "MPI-ESM1-2-HR", "EC-Earth3",
  "CNRM-CM6-1", "CNRM-CM6-1-HR", "IPSL-CM6A-LR", "MRI-ESM2-0"
)

ESGF_RUNID_SUFFIX <- "_historical"

# -----------------------------------------------------------------------------
# ANALYSIS DIMENSIONS
# -----------------------------------------------------------------------------

YEAR_START <- 1994
YEAR_END   <- 2014

EDD_BINS       <- c("edd_0","edd_4","edd_8","edd_12","edd_28","edd_30","edd_32")
EDD_THRESHOLDS <- c(0, 4, 8, 12, 28, 30, 32)

PERIODS    <- c("plant", "between", "harvest")
CROPS      <- c("wheat","maize","rice","barley","potato","sugar_beet","rapeseed")
IRRIGATION <- c("rainfed", "irrigated")

KEY_THRESHOLDS <- c("edd_8", "edd_30")

TAIL_QUANTILE  <- 0.90
EVAL_QUANTILES <- c(0.90, 0.95, 0.99)

# -----------------------------------------------------------------------------
# FIGURE SETTINGS
# -----------------------------------------------------------------------------

FIG_WIDTH_MAIN  <- 10
FIG_HEIGHT_MAIN <- 5
FIG_WIDTH_MAP   <- 12
FIG_HEIGHT_MAP  <- 8
FIG_DPI         <- 300

COL_IMPROVES <- "steelblue"
COL_WORSENS  <- "firebrick"
COL_NEUTRAL  <- "grey70"

MODEL_LABELS <- c(
  "GFDL-ESM4"     = "GFDL-ESM4",
  "MPI-ESM1-2-HR" = "MPI-ESM1-2-HR",
  "EC-Earth3"     = "EC-Earth3"
)

# -----------------------------------------------------------------------------
# KEY COLUMNS
# -----------------------------------------------------------------------------

KEY_COLS <- c("adm_code", "poly_idx", "year", "period", "crop", "irrigation", "calendar")

# adm_name added so metrics are at true ADM1 level (needed for spatial map)
GROUP_COLS_POOLED <- c("source_id", "adm_code", "adm_name",
                       "period", "crop", "irrigation", "edd_bin")

GROUP_COLS_CORR <- c("source_id", "adm_code", "adm_name",
                     "period", "crop", "irrigation", "edd_bin", "year")
