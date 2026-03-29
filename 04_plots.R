# =============================================================================
# 04_plots.R — All non-spatial figures
# =============================================================================

source("00_config.R")

library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(scales)

# =============================================================================
# 0. HELPERS
# =============================================================================

aggregate_to_crop <- function(df) {
  non_num    <- names(df)[!sapply(df, is.numeric)]
  group_vars <- setdiff(non_num, c("calendar", "n"))
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
}

THR_TO_POS <- setNames(seq_along(EDD_THRESHOLDS) - 1, sort(EDD_THRESHOLDS))

add_thr_label <- function(df) {
  df %>%
    mutate(
      thr      = as.numeric(str_extract(as.character(edd_bin), "[0-9]+")),
      axis_pos = THR_TO_POS[as.character(thr)]
    )
}

scale_x_edd_bar <- function() {
  scale_x_continuous(
    breaks = unname(THR_TO_POS),
    labels = names(THR_TO_POS),
    expand = expansion(mult = 0.05)
  )
}

geom_axis_break <- function() {
  annotate("rect",
           xmin = 3.35, xmax = 3.65,
           ymin = -Inf, ymax = Inf,
           fill = "white", colour = NA, alpha = 0.95)
}

theme_edd <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      legend.position  = "top",
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey40", size = base_size - 2),
      plot.caption     = element_text(colour = "grey50", size = base_size - 3)
    )
}

save_fig <- function(p, name, width = FIG_WIDTH_MAIN, height = FIG_HEIGHT_MAIN) {
  ggsave(file.path(OUT_FIGURES, paste0(name, ".png")),
         p, width = width, height = height, dpi = FIG_DPI)
  ggsave(file.path(OUT_FIGURES, paste0(name, ".pdf")),
         p, width = width, height = height)
  cat("Saved", name, "✓\n")
}

# Simple net bar: one bar per threshold, colored by sign of net ΔSkill
# No split blue/red — just the mean across all ADM1 regions
make_net_bars <- function(df, x_col, y_col, facet_formula,
                          title, subtitle = NULL, caption = NULL) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 fill = .data[[y_col]] >= 0)) +
    geom_col(width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.6) +
    scale_fill_manual(
      values = c("TRUE" = COL_IMPROVES, "FALSE" = COL_WORSENS),
      labels = c("TRUE" = "Improves", "FALSE" = "Worsens"),
      name   = NULL
    ) +
    scale_x_edd_bar() +
    geom_axis_break() +
    facet_grid(facet_formula) +
    labs(title = title, subtitle = subtitle,
         x = "EDD threshold (°C)",
         y = "ΔSkill MAE (Raw − Biasadj)",
         caption = caption) +
    theme_edd()
}

# Prepare net bar data
prep_net_data <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(c("source_id", group_vars, "thr", "axis_pos")))) %>%
    summarise(
      net = mean(dskill_mae, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(source_id = factor(source_id, levels = PAIRED_MODELS))
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("\n--- Loading metric tables ---\n")

dskill       <- read_parquet(file.path(OUT_TABLES, "dskill_paired.parquet"))
metrics_corr <- read_parquet(file.path(OUT_TABLES, "metrics_corr_paired.parquet"))
metrics_raw  <- read_parquet(file.path(OUT_TABLES, "metrics_pooled_raw_all.parquet"))
metrics_pair <- read_parquet(file.path(OUT_TABLES, "metrics_pooled_paired.parquet"))

cat("dskill rows      :", nrow(dskill), "\n")
cat("metrics_corr rows:", nrow(metrics_corr), "\n")
cat("metrics_raw rows :", nrow(metrics_raw), "\n")
cat("metrics_pair rows:", nrow(metrics_pair), "\n")

# country lookup: adm_code is already country-level in our data
country_lookup <- dskill %>%
  distinct(adm_code) %>%
  mutate(country = adm_code)

# =============================================================================
# 2. PREPARE
# =============================================================================

cat("\n--- Preparing data ---\n")

dskill_crop <- dskill %>%
  aggregate_to_crop() %>%
  add_thr_label()

metrics_pair_crop <- metrics_pair %>%
  aggregate_to_crop() %>%
  add_thr_label() %>%
  mutate(source_id = factor(source_id, levels = PAIRED_MODELS))

metrics_raw_crop <- metrics_raw %>%
  aggregate_to_crop() %>%
  add_thr_label() %>%
  mutate(source_id = factor(source_id, levels = ALL_ESGF_MODELS))

# =============================================================================
# 3. MAIN FIGURE: ΔSkill MAE by threshold × crop × model
# =============================================================================

cat("\n--- Fig: ΔSkill MAE by threshold, crop, model ---\n")

net_crop <- prep_net_data(dskill_crop, "crop") %>%
  mutate(crop = factor(crop, levels = CROPS))

p_main <- make_net_bars(
  net_crop, x_col = "axis_pos", y_col = "net",
  facet_formula = source_id ~ crop,
  title    = "Bias-adjustment skill: MAE improvement vs ERA5",
  subtitle = "Positive (blue) = bias adjustment reduces MAE; Negative (red) = bias adjustment increases MAE",
  caption  = "Net mean ΔSkill across all ADM1 regions, irrigation types, and growing season periods (1994–2014)"
)

save_fig(p_main, "Fig_dskill_mae_by_threshold_crop", width = 16, height = 8)

# =============================================================================
# 4. ΔSkill MAE BY PERIOD
# =============================================================================

cat("\n--- Fig: ΔSkill MAE by period ---\n")

net_period <- prep_net_data(dskill_crop, "period") %>%
  mutate(period = factor(period, levels = PERIODS))

p_period <- make_net_bars(
  net_period, x_col = "axis_pos", y_col = "net",
  facet_formula = source_id ~ period,
  title    = "Bias-adjustment skill by growing season period",
  subtitle = "Positive (blue) = bias adjustment reduces MAE; Negative (red) = bias adjustment increases MAE",
  caption  = "High thresholds (≥28°C) near zero: EDD accumulation negligible across most of Europe.\nNet mean ΔSkill across all ADM1 regions, crops, and irrigation types (1994–2014)"
)

save_fig(p_period, "Fig_dskill_mae_by_period", width = 12, height = 7)

# =============================================================================
# 5. ΔSkill MAE BY IRRIGATION
# =============================================================================

cat("\n--- Fig: ΔSkill MAE by irrigation ---\n")

net_irr <- prep_net_data(dskill_crop, "irrigation") %>%
  mutate(irrigation = factor(irrigation, levels = IRRIGATION))

p_irr <- make_net_bars(
  net_irr, x_col = "axis_pos", y_col = "net",
  facet_formula = source_id ~ irrigation,
  title    = "Bias-adjustment skill by irrigation type",
  subtitle = "Positive (blue) = bias adjustment reduces MAE; Negative (red) = bias adjustment increases MAE",
  caption  = "High thresholds (≥28°C) near zero: EDD accumulation negligible across most of Europe.\nNet mean ΔSkill across all ADM1 regions, crops, and growing season periods (1994–2014)"
)

save_fig(p_irr, "Fig_dskill_mae_by_irrigation", width = 10, height = 7)

# =============================================================================
# 6. QUANTILE ERRORS
# =============================================================================

cat("\n--- Fig: Quantile errors ---\n")

qerr <- metrics_pair_crop %>%
  group_by(source_id, source_type, thr, axis_pos) %>%
  summarise(
    q90_err = mean(q90_err, na.rm = TRUE),
    q95_err = mean(q95_err, na.rm = TRUE),
    q99_err = mean(q99_err, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(q90_err, q95_err, q99_err),
               names_to = "quantile", values_to = "error") %>%
  mutate(
    quantile    = factor(toupper(str_replace(str_replace(quantile,"_err",""),"q","Q")),
                         levels = c("Q90","Q95","Q99")),
    source_type = factor(source_type,
                         levels = c("raw","biasadj"),
                         labels = c("Raw CMIP6","Bias-adjusted"))
  )

p_qerr <- ggplot(qerr, aes(x = axis_pos, y = error,
                            colour = source_type, group = source_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_colour_manual(values = c("Raw CMIP6" = "grey40",
                                 "Bias-adjusted" = COL_IMPROVES)) +
  scale_x_edd_bar() +
  facet_grid(source_id ~ quantile) +
  labs(
    title    = "Quantile errors vs ERA5 (model quantile − ERA5 quantile)",
    subtitle = "Closer to zero = better reproduction of EDD distribution tails",
    x = "EDD threshold (°C)", y = "Quantile error (°C·days)", colour = NULL,
    caption = "Averaged across ADM1 regions, crops, irrigation types, and periods (1994–2014)"
  ) +
  theme_edd()

save_fig(p_qerr, "Fig_quantile_error", width = 12, height = 7)

# =============================================================================
# 7. TAIL MEAN ERROR
# =============================================================================

cat("\n--- Fig: Tail mean error ---\n")

tail_err <- metrics_pair_crop %>%
  group_by(source_id, source_type, thr, axis_pos) %>%
  summarise(tail_mean_err = mean(tail_mean_err, na.rm = TRUE), .groups = "drop") %>%
  mutate(source_type = factor(source_type,
                              levels = c("raw","biasadj"),
                              labels = c("Raw CMIP6","Bias-adjusted")))

p_tail <- ggplot(tail_err, aes(x = axis_pos, y = tail_mean_err,
                               colour = source_type, group = source_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_colour_manual(values = c("Raw CMIP6" = "grey40",
                                 "Bias-adjusted" = COL_IMPROVES)) +
  scale_x_edd_bar() +
  facet_wrap(~ source_id) +
  labs(
    title    = "Tail mean error vs ERA5 (top 10% EDD values)",
    subtitle = "Closer to zero = better reproduction of extreme heat exposure",
    x = "EDD threshold (°C)", y = "Tail mean error (°C·days)", colour = NULL,
    caption = "Averaged across ADM1 regions, crops, irrigation types, and periods (1994–2014)"
  ) +
  theme_edd()

save_fig(p_tail, "Fig_tail_mean_error", width = 12, height = 5)

# =============================================================================
# 8. KS STATISTIC
# =============================================================================

cat("\n--- Fig: KS statistic ---\n")

ks_df <- metrics_pair_crop %>%
  group_by(source_id, source_type, thr, axis_pos) %>%
  summarise(ks_stat = mean(ks_stat, na.rm = TRUE), .groups = "drop") %>%
  mutate(source_type = factor(source_type,
                              levels = c("raw","biasadj"),
                              labels = c("Raw CMIP6","Bias-adjusted")))

p_ks <- ggplot(ks_df, aes(x = axis_pos, y = ks_stat,
                           colour = source_type, group = source_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Raw CMIP6" = "grey40",
                                 "Bias-adjusted" = COL_IMPROVES)) +
  scale_x_edd_bar() +
  facet_wrap(~ source_id) +
  labs(
    title    = "Distributional distance vs ERA5 (KS statistic)",
    subtitle = "Lower = better reproduction of full EDD distribution",
    x = "EDD threshold (°C)", y = "KS statistic", colour = NULL,
    caption = "Averaged across ADM1 regions, crops, irrigation types, and periods (1994–2014)"
  ) +
  theme_edd()

save_fig(p_ks, "Fig_ks_stat", width = 12, height = 5)

# =============================================================================
# 9. INTERANNUAL CORRELATION ΔSkill
# =============================================================================

cat("\n--- Fig: Correlation ΔSkill ---\n")

corr_df <- dskill_crop %>%
  group_by(source_id, thr, axis_pos) %>%
  summarise(dskill_corr = mean(dskill_corr, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    sign_corr = ifelse(dskill_corr >= 0, "Improves", "Worsens"),
    source_id = factor(source_id, levels = PAIRED_MODELS)
  )

p_corr <- ggplot(corr_df, aes(x = axis_pos, y = dskill_corr, fill = sign_corr)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  scale_fill_manual(values = c("Improves" = COL_IMPROVES, "Worsens" = COL_WORSENS)) +
  scale_x_edd_bar() +
  facet_wrap(~ source_id) +
  labs(
    title    = "Bias-adjustment effect on interannual correlation vs ERA5",
    subtitle = "Positive = bias adjustment improves year-to-year variability agreement",
    x = "EDD threshold (°C)", y = "ΔCorr (Biasadj − Raw)", fill = NULL,
    caption = "Note: ΔCorr magnitudes are small (< 0.05) — bias adjustment has negligible effect on interannual variability.\nAveraged across ADM1 regions, crops, irrigation types, and periods (1994–2014)"
  ) +
  theme_edd()

save_fig(p_corr, "Fig_corr_dskill", width = 12, height = 5)

# =============================================================================
# 10. RAW-ALL MODELS: MAE vs ERA5
# =============================================================================

cat("\n--- Fig: Raw MAE all models ---\n")

raw_global <- metrics_raw_crop %>%
  group_by(source_id, thr, axis_pos) %>%
  summarise(mae = mean(mae, na.rm = TRUE), .groups = "drop")

p_raw_all <- ggplot(raw_global, aes(x = axis_pos, y = mae,
                                    colour = source_id, group = source_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_edd_bar() +
  labs(
    title    = "Raw CMIP6 MAE vs ERA5 — all models",
    subtitle = "Lower = better agreement with ERA5",
    x = "EDD threshold (°C)", y = "MAE vs ERA5 (°C·days)", colour = "Model",
    caption = "Averaged across ADM1 regions, crops, irrigation types, and periods (1994–2014)"
  ) +
  theme_edd()

save_fig(p_raw_all, "Fig_raw_mae_all_models", width = 10, height = 5)

# =============================================================================
# 11. COUNTRY HEATMAP — one figure per model
# adm_code is country-level in our data so no aggregation needed
# =============================================================================

cat("\n--- Fig: Country heatmap (one per model) ---\n")

dskill_country <- dskill_crop %>%
  group_by(source_id, adm_code, crop, thr) %>%
  summarise(dskill_mae = mean(dskill_mae, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    thr_label = factor(thr, levels = sort(EDD_THRESHOLDS)),
    crop      = factor(crop, levels = CROPS),
    source_id = factor(source_id, levels = PAIRED_MODELS)
  )

# Sort countries by mean ΔSkill across all thresholds and crops
country_order <- dskill_country %>%
  group_by(adm_code) %>%
  summarise(mean_dskill = mean(dskill_mae, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_dskill)) %>%
  pull(adm_code)

dskill_country <- dskill_country %>%
  mutate(adm_code = factor(adm_code, levels = country_order))

# Use full data range — no cap, no squish
max_val <- max(abs(dskill_country$dskill_mae), na.rm = TRUE)
lim     <- ceiling(max_val / 10) * 10

for (model in PAIRED_MODELS) {

  df_model <- dskill_country %>% filter(source_id == model)

  p_heat <- ggplot(df_model,
                   aes(x = thr_label, y = adm_code, fill = dskill_mae)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low      = COL_WORSENS,
      mid      = "white",
      high     = COL_IMPROVES,
      midpoint = 0,
      limits   = c(-lim, lim),
      na.value = "grey88",
      name     = "ΔSkill MAE\n(Raw − Biasadj)"
    ) +
    facet_wrap(~ crop, nrow = 1) +
    labs(
      title    = paste0("Country-level ΔSkill MAE — ", model),
      subtitle = "Positive (blue) = bias adjustment reduces error; Negative (red) = bias adjustment increases error",
      x        = "EDD threshold (°C)",
      y        = NULL,
      caption  = "Countries sorted by mean ΔSkill across all thresholds and crops (1994–2014).\nGrey = no harvested area data for that crop-country combination."
    ) +
    theme_edd(base_size = 10) +
    theme(
      axis.text.y   = element_text(size = 7),
      panel.spacing = unit(0.3, "lines")
    )

  safe_model <- str_replace_all(model, "[^A-Za-z0-9]", "_")
  save_fig(p_heat, paste0("Fig_country_heatmap_", safe_model),
           width = 16, height = 10)
}

# =============================================================================
# 12. APPENDIX: ΔSkill by calendar
# =============================================================================

cat("\n--- Fig (appendix): ΔSkill by calendar ---\n")

metrics_cal <- read_parquet(file.path(OUT_TABLES, "metrics_pooled_paired.parquet"))

dskill_cal <- metrics_cal %>%
  filter(source_type %in% c("raw","biasadj")) %>%
  select(source_id, source_type, adm_code, period,
         crop, irrigation, calendar, edd_bin, mae) %>%
  pivot_wider(names_from = source_type, values_from = mae,
              names_prefix = "mae_") %>%
  mutate(dskill_mae = mae_raw - mae_biasadj) %>%
  add_thr_label() %>%
  mutate(source_id = factor(source_id, levels = PAIRED_MODELS)) %>%
  group_by(source_id, calendar, thr, axis_pos) %>%
  summarise(mean_dskill_mae = mean(dskill_mae, na.rm = TRUE), .groups = "drop") %>%
  mutate(sign_mae = ifelse(mean_dskill_mae >= 0, "Improves", "Worsens"))

p_cal <- ggplot(dskill_cal,
                aes(x = axis_pos, y = mean_dskill_mae, fill = sign_mae)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("Improves" = COL_IMPROVES, "Worsens" = COL_WORSENS)) +
  scale_x_edd_bar() +
  geom_axis_break() +
  facet_grid(source_id ~ calendar) +
  labs(
    title    = "Bias-adjustment skill by crop calendar (appendix)",
    x        = "EDD threshold (°C)",
    y        = "ΔSkill MAE (Raw − Biasadj)",
    fill     = NULL,
    caption  = "Averaged across ADM1 regions, irrigation types, and periods (1994–2014)"
  ) +
  theme_edd(base_size = 10)

save_fig(p_cal, "Fig_dskill_mae_by_calendar_appendix", width = 18, height = 9)

cat("\n--- 04_plots.R complete ---\n")
cat("All figures saved to:", OUT_FIGURES, "\n")
