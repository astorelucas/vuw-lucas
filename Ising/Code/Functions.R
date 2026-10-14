library(ggplot2)
library(ggrepel)
library(ggthemes)
library(viridis)
library(StatOrdPattHxC)
library(patchwork)
library(dplyr)
library(png)
library(ordinalpatterns)

read_hxc <- function(method, beta, emb, bits, root = out_dir) {
  read.csv(file.path(root, method, sprintf("beta%.4f", beta),
                     sprintf("emb%d_bits%d.csv", emb, bits)))
}

binary_to_real <- function(window_binary) {
  n <- length(window_binary)
  weights <- 2^((n - 1):0) # MSB first, LSB last
  sum(window_binary * weights)
}

compute_hxc <- function(x, emb) {
  op_probs <- op_pd(x, D = emb, tau=1)
  
  probs <- op_probs$probabilities
  h_raw <- permutation_entropy(probs) # raw entropy, up to ln(D!)
  h <- h_raw / log(factorial(emb))    # normalize to [0,1]
  c <- statistical_complexity(probs)
  
  c(H = h, C = c)
}

# Converte distâncias d da curva de Hilbert em (x, y) numa grade n x n.
# n deve ser potência de 2; d é um vetor de inteiros em 0..n^2-1.
hilbert_d2xy <- function(n, d) {
  t <- as.integer(d)
  x <- integer(length(t))
  y <- integer(length(t))
  s <- 1L
  
  while (s < n) {
    rx <- bitwAnd(bitwShiftR(t, 1L), 1L)
    ry <- bitwAnd(bitwXor(t, rx), 1L)
    
    # onde ry == 0: reflete se rx == 1 e depois troca x <-> y
    rot  <- ry == 0L
    flip <- rot & (rx == 1L)
    x[flip] <- s - 1L - x[flip]
    y[flip] <- s - 1L - y[flip]
    tmp    <- x[rot]
    x[rot] <- y[rot]
    y[rot] <- tmp
    
    x <- x + s * rx
    y <- y + s * ry
    t <- bitwShiftR(t, 2L)
    s <- s * 2L
  }
  
  if (length(t) == 1L) c(x = x, y = y) else cbind(x = x, y = y)
}

# Full traversal order for an n x n grid as a matrix of 1-indexed
# (row, col) pairs, ordered along the Hilbert curve (d = 0 .. n^2-1).
hilbert_curve_order <- function(n) {
  n <- as.integer(n)
  stopifnot(bitwAnd(n, n - 1L) == 0L) # n must be a power of 2
  coords <- t(vapply(0:(n^2 - 1), function(d) hilbert_d2xy(n, d), numeric(2)))
  cbind(row = coords[, 2] + 1L, col = coords[, 1] + 1L)
}


run_hxc_from_simulation_ph <- function(
    maps,                  # lista de matrizes (uma por geração) ou uma única matriz
    emb,
    group_size = 16,
    seed,
    order_coords = NULL    # opcional: calcule uma vez fora e reutilize
) {
  # aceita uma única matriz
  if (is.matrix(maps)) maps <- list(maps)
  
  # preserva o gerador de números aleatórios global
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  
  n_maps <- length(maps)
  n      <- nrow(maps[[1]])
  
  if (is.null(order_coords)) order_coords <- hilbert_curve_order(n)
  n_total <- nrow(order_coords)
  stopifnot(n_total %% group_size == 0)
  
  time_series_list <- vector("list", n_maps)
  hc_results       <- vector("list", n_maps)
  
  for (t in seq_len(n_maps)) {
    map  <- maps[[t]]
    vals <- map[order_coords]            # valores na ordem da curva de Hilbert
    
    # spins -1/+1 -> bits 0/1 (só converte se houver -1)
    if (any(vals < 0L)) vals <- (vals + 1L) %/% 2L
    
    # cada coluna da matriz = um grupo de group_size valores consecutivos
    time.series <- apply(matrix(vals, nrow = group_size), 2, binary_to_real)
    
    time_series_list[[t]] <- time.series
    hc_results[[t]]       <- compute_hxc(time.series, emb = emb)
  }
  
  # não devolve `maps`, para não carregar as matrizes na saída
  list(
    time_series = time_series_list,
    hc_results  = hc_results
  )
}

read_gens_list <- function(file, from, to, L = 512) {
  nbytes <- L * L / 8
  n <- to - from + 1
  con <- file(file, "rb"); on.exit(close(con))
  seek(con, where = from * nbytes)
  lapply(seq_len(n), function(i) {
    bits <- rawToBits(readBin(con, "raw", n = nbytes))
    matrix(as.integer(bits) * 2L - 1L, L, L)
  })
}
