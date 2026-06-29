# ============================================================
# Game of Life - Sliding Window Binary Encoding -> HxC Analysis
# ============================================================
# For each epoch t:
#   1. Evolve the GoL map one step (map_{t-1} -> map_t)
#   2. Slide a non-overlapping 4x4 window across the map
#   3. Flatten each window to a 16-bit binary vector
#   4. Convert that binary vector to a real (decimal) number
#   5. Collect these numbers into a time series for epoch t
#   6. Compute permutation entropy (H) and statistical
#      complexity (C) on that time series via OPprob/HShannon/StatComplexity
#
# Requires: StatOrdPattHxC package (for OPprob, HShannon, StatComplexity)
# install.packages("StatOrdPattHxC")  # if not already installed
# ============================================================

library(StatOrdPattHxC)

# ---------------------------------------------------------------
# Helper: convert a binary vector to its decimal (real) value.
# Convention: LAST element is the least significant bit.
#   e.g. c(0, 0, 1) -> 1 ,  c(1, 0, 0) -> 4 ,  c(1, 1, 1) -> 7
# ---------------------------------------------------------------
binary_to_real <- function(window_binary) {
  n <- length(window_binary)
  weights <- 2^((n - 1):0)   # MSB first, LSB last
  sum(window_binary * weights)
}

# ---------------------------------------------------------------
# Helper: count live neighbors for every cell (toroidal boundary)
# Reused from the Game of Life script.
# ---------------------------------------------------------------
count_neighbors <- function(mat) {
  nrow_m <- nrow(mat)
  ncol_m <- ncol(mat)
  neighbor_count <- matrix(0, nrow = nrow_m, ncol = ncol_m)
  
  offsets <- list(
    c(-1, -1), c(-1, 0), c(-1, 1),
    c(0, -1),            c(0, 1),
    c(1, -1),  c(1, 0),  c(1, 1)
  )
  
  for (off in offsets) {
    shifted_rows <- ((1:nrow_m - 1 + off[1]) %% nrow_m) + 1
    shifted_cols <- ((1:ncol_m - 1 + off[2]) %% ncol_m) + 1
    neighbor_count <- neighbor_count + mat[shifted_rows, shifted_cols]
  }
  
  neighbor_count
}

# ---------------------------------------------------------------
# Helper: one Game of Life update step (map_{t-1} -> map_t)
# ---------------------------------------------------------------
next_generation <- function(mat) {
  neighbors <- count_neighbors(mat)
  survive <- (mat == 1) & (neighbors %in% c(2, 3))
  born    <- (mat == 0) & (neighbors == 3)
  
  new_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  new_mat[survive | born] <- 1
  new_mat
}

# ---------------------------------------------------------------
# Helper: compute permutation entropy (H) and statistical
# complexity (C) for a time series, following the same recipe
# used for the sea-level station HxC analysis:
#   OPprob() -> normalize -> HShannon() / StatComplexity()
# Returns a named vector c(H = ..., C = ...)
# ---------------------------------------------------------------
compute_hxc <- function(x, emb) {
  op_probs <- OPprob(x, emb = emb)
  op_probs <- op_probs / sum(op_probs)
  
  h <- HShannon(op_probs)
  c <- StatComplexity(op_probs)
  
  c(H = h, C = c)
}

# ---------------------------------------------------------------
# Main routine
#
# Inputs:
#   seed         - RNG seed for reproducible initial map
#   k            - image dimension is 2^k (e.g. k = 8 -> 256x256)
#   T            - number of epochs to simulate
#   emb          - embedding dimension passed to OPprob() (ordinal pattern length)
#   window_size  - side length of the sliding window (default 4)
#   stride       - step size between windows (default 4, i.e. non-overlapping)
#
# Returns a list with:
#   maps          - list of GoL matrices, one per epoch (t = 1..T)
#   time_series   - list of numeric time series, one per epoch
#   hc_results    - list of c(H, C) named vectors, one per epoch
# ---------------------------------------------------------------
run_gol_hxc_analysis <- function(seed, k, T, emb, window_size = 4, stride = 4) {
  
  n <- 2^k  # grid dimension
  
  set.seed(seed)
  map <- matrix(sample(c(0, 1), n * n, replace = TRUE), nrow = n, ncol = n)
  
  # window start positions: 1, 1+stride, ... such that window fits within n
  starts <- seq(1, n - window_size + 1, by = stride)
  
  maps        <- vector("list", T)
  time_series_list <- vector("list", T)
  hc_results  <- vector("list", T)
  
  for (t in 1:T) {
    
    # map_t from map_{t-1}
    map <- next_generation(map)
    
    time.series <- numeric(0)
    
    for (i in starts) {
      for (j in starts) {
        window.binary <- as.vector(map[i:(i + window_size - 1),
                                       j:(j + window_size - 1)])
        window.real <- binary_to_real(window.binary)
        time.series <- c(time.series, window.real)
      }
    }
    
    hc <- compute_hxc(time.series, emb = emb)
    
    maps[[t]]             <- map
    time_series_list[[t]] <- time.series
    hc_results[[t]]       <- hc
  }
  
  list(
    maps        = maps,
    time_series = time_series_list,
    hc_results  = hc_results
  )
}

# ============================================================
# Example usage
# ============================================================

results <- run_gol_hxc_analysis(
  seed = 123456789,
  k    = 8,    # 2^8 = 256x256 grid
  T    = 2999,
  emb  = 3     # embedding dimension passed to OPprob()
)

cat("Number of epochs simulated:", length(results$maps), "\n")
cat("Length of time series per epoch:", length(results$time_series[[1]]), "\n")
cat("H, C for epoch 1:\n")
print(results$hc_results[[1]])

# Build a data frame of H and C across all epochs, ready for plotting
# on the complexity-entropy plane (same style as the station-level plot):
hc_df <- data.frame(do.call(rbind, results$hc_results))
hc_df$epoch <- 1:nrow(hc_df)
print(hc_df)

library(ggplot2)
library(ggrepel)
library(ggthemes)

hc_df$epoch <- factor(hc_df$epoch)  # discrete labels/colors, one per epoch
hc_df$H <- as.numeric(hc_df$H)
hc_df$C <- as.numeric(hc_df$C)
# ACF Acrescentar: remover linhas em que H ou C são NA

emb_used <- 3  # set to whatever 'emb' you used above, for the LinfLsup filter

# ACF No gráfico, identificar um a cada K pontos,
#     K é variável de entrada, e.g., 500
ggplot(hc_df, aes(x = H, y = C, color = epoch)) +
  geom_line(data = subset(LinfLsup, Dimension == emb_used & Side == "Lower"),
            aes(x = H, y = C), color = "grey40", linewidth = 0.6,
            inherit.aes = FALSE) +
  geom_line(data = subset(LinfLsup, Dimension == emb_used & Side == "Upper"),
            aes(x = H, y = C), color = "grey40", linewidth = 0.6,
            inherit.aes = FALSE) +
  geom_point(size = 3, alpha = 0.85) +
  # geom_text_repel(aes(label = epoch),
  #                 size = 3,
  #                 fontface = "bold",
  #                 max.overlaps = Inf,
  #                 box.padding = 0.4,
  #                 point.padding = 0.3,
  #                 segment.color = "grey60",
  #                 segment.size = 0.4,
  #                 seed = 42) +
  scale_x_continuous(limits = c(min(hc_df$H), max(hc_df$H))) +
  scale_y_continuous(limits = c(min(hc_df$C), max(hc_df$C))) +
  labs(title = "HxC Complexity-Entropy Plane",
       subtitle = paste0("Game of Life — Embedding dimension ", emb_used),
       x = "Permutation Entropy (H)",
       y = "Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey40", size = 10),
        axis.title = element_text(size = 11))

# ACF Organizar o códig: funções e programas nas respectivas pastas
#     Salvar os dados na pasta "Data"