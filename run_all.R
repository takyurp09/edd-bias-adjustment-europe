# =============================================================================
# run_all.R
# Master script — runs the full EDD bias-adjustment evaluation pipeline
#
# Usage:
#   Rscript run_all.R
#
# Order:
#   01_load_harmonize.R  — load, parse, harmonize, save spines
#   02_metrics.R         — compute all evaluation metrics
#   03_maps.R            — spatial figures
#   04_plots.R           — non-spatial figures
# =============================================================================

start_time <- proc.time()

steps <- list(
  list(script = "01_load_harmonize.R", label = "Load & Harmonize"),
  list(script = "02_metrics.R",        label = "Compute Metrics"),
  list(script = "03_maps.R",           label = "Map Figures"),
  list(script = "04_plots.R",          label = "Non-spatial Figures")
)

cat("\n", strrep("=", 60), "\n")
cat("  EDD Bias-Adjustment Evaluation Pipeline\n")
cat(strrep("=", 60), "\n\n")

for (i in seq_along(steps)) {

  step   <- steps[[i]]
  label  <- step$label
  script <- step$script

  cat(strrep("-", 60), "\n")
  cat(sprintf("Step %d / %d : %s\n", i, length(steps), label))
  cat(sprintf("Script     : %s\n", script))
  cat(strrep("-", 60), "\n")

  step_start <- proc.time()

  tryCatch(
    {
      source(script)
      elapsed <- (proc.time() - step_start)[["elapsed"]]
      cat(sprintf("\n✓ Step %d complete in %.1f seconds\n\n", i, elapsed))
    },
    error = function(e) {
      cat(sprintf("\n✗ Step %d FAILED: %s\n", i, conditionMessage(e)))
      cat("Pipeline stopped. Fix the error above and re-run.\n\n")
      stop(e)
    }
  )
}

total <- (proc.time() - start_time)[["elapsed"]]

cat("\n", strrep("=", 60), "\n")
cat(sprintf("  Pipeline complete in %.1f seconds (%.1f minutes)\n",
            total, total / 60))
cat(strrep("=", 60), "\n\n")
