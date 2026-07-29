library(changepoint)
library(ggplot2)

detect_oscillation_onset <- function(hc_df, burn_in = 20) {
  hc_df <- hc_df[order(hc_df$t), ]
  conv_data <- subset(hc_df, dist_to_limit > 1e-8)
  
  fit <- lm(log(dist_to_limit) ~ t, data = conv_data)
  conv_data$resid <- residuals(fit)
  
  # drop the initial transient — early epochs are naturally high-variance
  # because dist_to_limit is large before the trajectory has approached
  # (H*, C*) at all, which isn't "oscillation"
  conv_data_trim <- conv_data[-(1:min(burn_in, nrow(conv_data) - 1)), ]
  
  cpt <- cpt.var(conv_data_trim$resid, method = "AMOC")
  onset_idx <- cpts(cpt)
  
  onset_epoch <- if (length(onset_idx) == 0) NA else conv_data_trim$t[onset_idx]
  
  list(onset_epoch = onset_epoch, conv_data = conv_data_trim)
}

plot_oscillation_diagnostic <- function(conv_data, onset_epoch, run_tag) {
  p1 <- ggplot(conv_data, aes(x = t, y = log(dist_to_limit))) +
    geom_line(color = "grey30", linewidth = 0.4) +
    geom_point(size = 0.8, alpha = 0.6) +
    { if (!is.na(onset_epoch))
      geom_vline(xintercept = onset_epoch, color = "firebrick",
                 linewidth = 0.8, linetype = "dashed") } +
    labs(title = run_tag, x = "Epoch", y = expression(log(d[to~limit]))) +
    theme_minimal(base_size = 11)
  
  p2 <- ggplot(conv_data, aes(x = t, y = resid)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_line(color = "steelblue", linewidth = 0.4) +
    { if (!is.na(onset_epoch))
      geom_vline(xintercept = onset_epoch, color = "firebrick",
                 linewidth = 0.8, linetype = "dashed") } +
    labs(x = "Epoch", y = "Residual (log-fit)") +
    theme_minimal(base_size = 11)
  
  list(dist_plot = p1, resid_plot = p2)
}

# ------------------------------------------------------------
# Spot-check a handful of files before running the full batch
# ------------------------------------------------------------
spot_check_files <- csv_files[sample(seq_along(csv_files), min(6, length(csv_files)))]

for (f in spot_check_files) {
  hc_df <- read.csv(f)
  run_tag <- sub("\\.csv$", "", sub("^results_", "", f))
  
  onset <- detect_oscillation_onset(hc_df, min_seg = 20)
  plots <- plot_oscillation_diagnostic(onset$conv_data, onset$onset_epoch, run_tag)
  
  ggsave(sprintf("diag_dist_%s.pdf", run_tag), plots$dist_plot, width = 7, height = 4)
  ggsave(sprintf("diag_resid_%s.pdf", run_tag), plots$resid_plot, width = 7, height = 4)
  
  cat(run_tag, "-> onset:", onset$onset_epoch,
      "| n_changepoints:", onset$n_changepoints, "\n")
}

# ------------------------------------------------------------
# Full batch (run only after spot-checks look right)
# ------------------------------------------------------------
osc_results <- lapply(csv_files, function(f) {
  hc_df <- read.csv(f)
  run_tag <- sub("\\.csv$", "", sub("^results_", "", f))
  
  onset <- detect_oscillation_onset(hc_df, min_seg = 30)
  
  parts <- regmatches(f, regexec("results_(\\d+)k_([0-9.]+)_emb(\\d+)_(\\d+)bits\\.csv", f))[[1]]
  
  data.frame(
    epochs_k    = as.numeric(parts[2]),
    live_pct    = as.numeric(parts[3]),
    emb         = as.numeric(parts[4]),
    bits        = as.numeric(parts[5]),
    onset_epoch = onset$onset_epoch,
    n_cpts      = onset$n_changepoints
  )
})

osc_summary <- do.call(rbind, osc_results)
write.csv(osc_summary, "oscillation_onset_summary.csv", row.names = FALSE)

