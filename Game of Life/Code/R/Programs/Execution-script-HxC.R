source("Game of Life/Code/R/Programs/Functions.R")
library(progress)
#library(remotes)
#remotes::install_github("labepi/ordinalpatterns")
library(ordinalpatterns)

# Definição de Espaço Paramétrico
print("Starting sweep...")
emb_dims <- c(5,6)
bits_vals <- c(8,16,32)
epochs_k <- c(10)
live_pct <- c(0.1,0.3,0.5,0.8)

n_repetitions<- 100

seed <- 12345678
N_tail <- 1000 # epochs used to estimate (H*, C*)

total_runs <- length(epochs_k) * length(live_pct) * length(emb_dims) * length(bits_vals)

pb <- progress_bar$new(
  format = "  [:bar] :percent | :current/:total | run: :run_tag | elapsed: :elapsed | eta: :eta",
  total = total_runs,
  width = 100
)

timing_log <- tibble::tibble(
  epochs_k = integer(),
  live_pct = numeric(),
  emb = integer(),
  bits = integer(),
  start_time = character(),
  end_time = character(),
  elapsed_sec = numeric()
)

log_file <- "timing_log.csv"
readr::write_csv(timing_log, log_file)  # cria o arquivo com cabeçalho

for (k in epochs_k) {
  for (p in live_pct) {
    for (i in 1:n_repetitions){
      rds_file <- sprintf("Data/ca-simulations/%sK/ca_history_%sk_%s_%s.rds", k, k, p, i)
      if (!file.exists(rds_file)) {
        cat("Missing file, skipping:", rds_file, "\n")
        next
      }
      results_CAsimulation <- readRDS(rds_file)
      
      for (emb_used in emb_dims) {
        for (b in bits_vals) {
          run_tag <- sprintf("%sk_%s_emb%s_%sbits_%s", k, p, emb_used, b,i)
          
          t_start <- Sys.time()
          cat(sprintf("[%s] -----> Running: %s\n",
                      format(t_start, "%H:%M:%S"), run_tag))
          
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
          saveRDS(
            conv_data,
            file = sprintf("Data/results/%sK/new_experiments/results_%s.rds", k, run_tag)
          )
          
          t_end <- Sys.time()
          elapsed <- as.numeric(difftime(t_end, t_start, units = "secs"))
          
          readr::write_csv(
            tibble::tibble(
              epochs_k = k, live_pct = p, emb = emb_used, bits = b,
              start_time = format(t_start, "%Y-%m-%d %H:%M:%S"),
              end_time = format(t_end, "%Y-%m-%d %H:%M:%S"),
              elapsed_sec = elapsed
            ),
            log_file, append = TRUE
          )
          
          cat(sprintf("      finalizado em %.2fs\n", elapsed))
        }
      }
    }
  }
}