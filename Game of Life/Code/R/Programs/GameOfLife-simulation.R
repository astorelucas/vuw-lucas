source("Game of Life/Code/R/Programs/Functions.R")

# ============================================================
# 2D Cellular Automaton Simulator
# ============================================================
n_rows <- 1024
n_cols <- 1024
n_steps <- 9999
p_alive <- c(0.1,0.3,0.5,0.8) # percentage_alive in initial_state
n_repeticoes <- 100
epoch <- 10

tempo_total_inicio <- Sys.time()

for (p in p_alive) {
  for (i in 1:n_repeticoes){
    seed <- as.integer(p * 1000 + i)
    cat("--> Running:","rep=",i," for p-alive ", p ,"\n")
    
    tempo_rep_inicio <- Sys.time()
    
    init_mat <- make_random_init(n_rows, n_cols, p, seed=seed)
    
    results_CAsimulation <- simulate_ca(
      initial_mat = init_mat,
      n_steps     = n_steps
    )
    
    out_file <- sprintf("Data/ca-simulations/%sK/new_simulations/ca_history_%sk_%s_%s.rds", epoch, epoch, p,i)
    saveRDS(results_CAsimulation, file = out_file, compress = "xz")
    cat("Saved full history to", out_file, "\n")
    
    tempo_rep_fim <- Sys.time()
    duracao_rep <- as.numeric(difftime(tempo_rep_fim, tempo_rep_inicio, units = "secs"))
    cat(sprintf("Tempo de execução (p=%s, rep=%s): %.2f segundos (%.2f min)\n\n", p, i, duracao_rep, duracao_rep / 60))
    
    rm(results_CAsimulation, init_mat); gc()
  }
}

tempo_total_fim <- Sys.time()
duracao_total <- as.numeric(difftime(tempo_total_fim, tempo_total_inicio, units = "secs"))
cat(sprintf("\n=== Tempo total de execução: %.2f segundos (%.2f min / %.2f h) ===\n",
            duracao_total, duracao_total / 60, duracao_total / 3600))