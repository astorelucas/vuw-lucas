# ============================================================
# 2D Cellular Automaton Simulator
# ============================================================
n_rows <- 256
n_cols <- 256
n_steps <- c(2999, 4999, 9999, 14999)
p_alive <- c(0.1, 0.3, 0.5, 0.8) # percentage_alive in initial_state

for (n_steps in n_steps_vals) {
  for (p in p_alive_vals) {
    
    k <- round((n_steps + 1) / 1000)
    
    init_mat <- make_random_init(n_rows, n_cols, p)
    
    results_CAsimulation <- simulate_ca(
      initial_mat = init_mat,
      n_steps     = n_steps
    )
    
    cat("Simulation complete:", length(results_CAsimulation), "time steps generated for",
        sprintf("%sk_%s", k, p), "\n")
    
    out_file <- sprintf("ca_history_%sk_%s.rds", k, p)
    saveRDS(results_CAsimulation, file = out_file, compress = "xz")
    cat("Saved full history to", out_file, "\n\n")
    
    rm(results_CAsimulation, init_mat); gc()
  }
}