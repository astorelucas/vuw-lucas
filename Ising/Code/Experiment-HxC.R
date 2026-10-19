library(ggplot2)
library(ggthemes)
library(here)
library(Rcpp)

sourceCpp(here("Code", "ising.cpp"))  
source(here("Code", "Functions.R"))

# ------------------------------------------------------------
# Variáveis referentes à simulação
# ------------------------------------------------------------

beta_vals   <- c(log(1 + sqrt(2)) / 2) * 1.5
generations <- c(1000)
methods     <- c("GB", "MH", "SW")
L           <- 512L

# Variáveis referentes aos padrões ordinais
embdim <- c(6)
bits   <- c(16)
seed  <- "1000-99-01-1.5"

step_fun <- list(MH = metropolis_step,
                 GB = gibbs_step,
                 SW = swendsen_wang_step)

order_coords <- hilbert_curve_order(L)
coords  <- 0:(L * L - 1L)                 # distâncias d ao longo da curva

# Same initial condition
initial <- matrix(
  sample(c(-1L, 1L), L * L, replace = TRUE, prob = c(0.999999, 1-0.999999)),
  nrow = L, ncol = L
)

out_dir <- here("Results", sprintf("seed%s", seed))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(initial, file.path(out_dir, sprintf("initial_seed%d.rds", seed)))

record_every <- 1
max_gen      <- max(generations)
combos       <- expand.grid(embdim = embdim, bits = bits)
gens_rec     <- seq(0, max_gen, by = record_every)

combo_file <- function(run_dir, e, b) {
  file.path(run_dir, sprintf("emb%d_bits%d.csv", e, b))
}
seed <- 901
for (im in seq_along(methods)) {
  m    <- methods[im]
  step <- step_fun[[m]]
  
  for (ib in seq_along(beta_vals)) {
    beta    <- beta_vals[ib]
    run_dir <- file.path(out_dir, m, sprintf("beta%.4f", beta))
    
    last <- combos[nrow(combos), ]
    if (file.exists(combo_file(run_dir, last$embdim, last$bits))) {
      message("Pulando (já existe): ", m, " beta = ", sprintf("%.4f", beta))
      next
    }
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    
    message(sprintf("Rodando %s, beta = %.4f (%d/%d métodos, %d/%d betas)",
                    m, beta, im, length(methods), ib, length(beta_vals)))
    t0 <- Sys.time()
    set.seed(seed + 1000 * im + ib)
    spins <- initial + 0L
    
    H_mat <- matrix(NA_real_, length(gens_rec), nrow(combos))
    C_mat <- matrix(NA_real_, length(gens_rec), nrow(combos))
    r <- 0L
    
    pb <- txtProgressBar(min = 0, max = max_gen, style = 3)
    
    for (g in 0:max_gen) {
      if (g > 0) spins <- step(spins, beta)
      
      if (g %% record_every == 0) {
        r <- r + 1L
        for (ic in seq_len(nrow(combos))) {
          hxc_out <- run_hxc_from_simulation_ph(
            spins,
            emb          = combos$embdim[ic],
            group_size   = combos$bits[ic],
            seed         = seed,
            order_coords = order_coords
          )
          hc <- hxc_out$hc_results[[1]]
          H_mat[r, ic] <- hc[["H"]]
          C_mat[r, ic] <- hc[["C"]]
        }
      }
      setTxtProgressBar(pb, g)
    }
    close(pb)
    
    # um CSV por combinação, só com generation, H, C
    message(sprintf("Salvando %s, beta = %.2f, embdim = %d, bits = %d", m, beta, combos$embdim[ic], combos$bits[ic]))
    for (ic in seq_len(nrow(combos))) {
      write.csv(
        data.frame(generation = gens_rec, H = H_mat[, ic], C = C_mat[, ic]),
        combo_file(run_dir, combos$embdim[ic], combos$bits[ic]),
        row.names = FALSE
      )
    }
    
    # uma linha por execução com os metadados
    run_info <- data.frame(
      seed = seed, L = L, method = m, beta = beta,
      max_gen = max_gen, record_every = record_every,
      run_time_s = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    info_file <- file.path(out_dir, "run_info.csv")
    write.table(run_info, info_file, sep = ",", row.names = FALSE,
                col.names = !file.exists(info_file),
                append = file.exists(info_file))
  }
}
