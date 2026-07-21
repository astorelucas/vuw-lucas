source("Code/R/Programs/Functions.R")
library(ordinalpatterns)

emb_dims <- c(4, 5, 6)
bits_vals <- c(16, 32, 64)
epochs_k <- c(3, 5, 10)
live_pct <- c(0.1, 0.3, 0.5, 0.8)

seed <- 12345678
N_tail <- 300 # epochs used to estimate (H*, C*)

for (k in epochs_k) {
  for (p in live_pct) {
    rds_file <- sprintf("Data/ca_history%sk_%s.rds", k, p)
    if (!file.exists(rds_file)) {
      cat("Missing file, skipping:", rds_file, "\n")
      next
    }
    results_CAsimulation <- readRDS(rds_file)

    for (emb_used in emb_dims) {
      for (b in bits_vals) {
        run_tag <- sprintf("%sk_%s_emb%s_%sbits", k, p, emb_used, b)
        cat("-----> Running:", run_tag, "\n")

        hxc_out <- run_gol_hxc_analysis_from_simulation_pht(
          results_CAsimulation,
          emb = emb_used,
          group_size = b,
          seed = seed
        )

        hc_df <- data.frame(do.call(rbind, hxc_out$hc_results))
        hc_df$epoch <- 1:nrow(hc_df)
        hc_df$H <- as.numeric(hc_df$H)
        hc_df$C <- as.numeric(hc_df$C)
        hc_df <- hc_df[!is.na(hc_df$H) & !is.na(hc_df$C), ]
        hc_df$t <- hc_df$epoch
        hc_df <- hc_df[order(hc_df$t), ]

        hc_df$velocity <- c(NA, sqrt(diff(hc_df$H)^2 + diff(hc_df$C)^2))

        H_star <- mean(tail(hc_df$H, N_tail))
        C_star <- mean(tail(hc_df$C, N_tail))
        hc_df$dist_to_limit <- sqrt((hc_df$H - H_star)^2 + (hc_df$C - C_star)^2)

        conv_data <- subset(hc_df, dist_to_limit > 1e-8)
        write.csv(
          conv_data,
          file = sprintf("results_%s.csv", run_tag),
          row.names = FALSE
        )
      }
    }
  }
}
