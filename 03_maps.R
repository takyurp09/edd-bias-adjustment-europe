# =============================================================================
# 03_maps.R
# One spatial map: ΔSkill MAE at edd_8, true ADM1 level
# 3 panels side by side (one per paired model)
# =============================================================================

source("00_config.R")

library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(sf)
library(ggplot2)

sf_use_s2(FALSE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("\n--- Loading data ---\n")

shp <- st_read(SHAPEFILE_PATH, quiet = TRUE) %>%
  st_make_valid() %>%
  st_simplify(dTolerance = 0.05, preserveTopology = TRUE) %>%
  st_transform(crs = 4326) %>%
  select(GID_1, adm_name = NAME_1, geometry)

cat("Shapefile polygons:", nrow(shp), "\n")

dskill <- read_parquet(file.path(OUT_TABLES, "dskill_paired.parquet"))
cat("dskill rows:", nrow(dskill), "\n")
cat("dskill columns:", paste(names(dskill), collapse=", "), "\n")

# =============================================================================
# 2. PREPARE MAP DATA
# adm_name is now in GROUP_COLS_POOLED so dskill has true ADM1-level data
# =============================================================================

cat("\n--- Preparing map data ---\n")

# Aggregate calendars first
aggregate_to_crop <- function(df) {
  non_num    <- names(df)[!sapply(df, is.numeric)]
  group_vars <- setdiff(non_num, c("calendar", "n"))
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
}

dskill_map <- dskill %>%
  aggregate_to_crop() %>%
  filter(as.character(edd_bin) == "edd_8") %>%
  group_by(source_id, adm_name) %>%
  summarise(
    dskill_mae = mean(dskill_mae, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(source_id = factor(source_id, levels = PAIRED_MODELS))

cat("Map data rows:", nrow(dskill_map), "\n")
cat("Unique adm_name in dskill_map:", length(unique(dskill_map$adm_name)), "\n")

# Join to shapefile using adm_name
map_sf <- shp %>%
  left_join(dskill_map, by = "adm_name")

cat("map_sf rows:", nrow(map_sf), "\n")
cat("Non-NA dskill_mae:", sum(!is.na(map_sf$dskill_mae)), "\n")

# =============================================================================
# 3. DRAW MAP
# =============================================================================

cat("\n--- Drawing map ---\n")

# Full data range — no cap
max_val <- max(abs(dskill_map$dskill_mae), na.rm = TRUE)
lim     <- ceiling(max_val / 5) * 5

p_map <- ggplot() +
  geom_sf(data = map_sf,
          aes(fill = dskill_mae),
          colour    = "white",
          linewidth = 0.1) +
  scale_fill_gradient2(
    low      = COL_WORSENS,
    mid      = "white",
    high     = COL_IMPROVES,
    midpoint = 0,
    limits   = c(-lim, lim),
    name     = "ΔSkill MAE\n(Raw − Biasadj)",
    na.value = "grey85"
  ) +
  coord_sf(
    xlim        = c(-25, 60),
    ylim        = c(27, 72),
    expand      = FALSE,
    default_crs = st_crs(4326)
  ) +
  facet_wrap(~ source_id, nrow = 1) +
  labs(
    title    = "Spatial distribution of bias-adjustment skill — EDD threshold 8°C",
    subtitle = "Positive (blue) = bias adjustment reduces MAE vs ERA5; Negative (red) = bias adjustment increases MAE",
    caption  = "Averaged across crops, irrigation types, and growing season periods (1994–2014). Grey = regions with no matching data."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid       = element_line(color = "grey90", linewidth = 0.2),
    axis.text        = element_text(size = 7),
    axis.title       = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

ggsave(file.path(OUT_FIGURES, "Fig_map_dskill_edd8.png"),
       p_map, width = 14, height = 6, dpi = FIG_DPI, bg = "white")
ggsave(file.path(OUT_FIGURES, "Fig_map_dskill_edd8.pdf"),
       p_map, width = 14, height = 6)
cat("Saved Fig_map_dskill_edd8 ✓\n")

cat("\n--- 03_maps.R complete ---\n")
