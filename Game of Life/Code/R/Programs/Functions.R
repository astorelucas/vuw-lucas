library(ggplot2)
library(ggrepel)
library(ggthemes)
library(viridis)
library(StatOrdPattHxC)
library(patchwork)
library(dplyr)
library(png)

# ---------------------------------------------------------------
# Helper: convert a binary vector to its decimal (real) value.
# Convention: LAST element is the least significant bit.
#   e.g. c(0, 0, 1) -> 1 ,  c(1, 0, 0) -> 4 ,  c(1, 1, 1) -> 7
# ---------------------------------------------------------------
binary_to_real <- function(window_binary) {
  n <- length(window_binary)
  weights <- 2^((n - 1):0) # MSB first, LSB last
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
    c(-1, -1),
    c(-1, 0),
    c(-1, 1),
    c(0, -1),
    c(0, 1),
    c(1, -1),
    c(1, 0),
    c(1, 1)
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
next_generation_generic <- function(mat, born, survive) {
  neighbors <- count_neighbors(mat)

  survive_mask <- (mat == 1) & (neighbors %in% survive)
  born_mask <- (mat == 0) & (neighbors %in% born)

  new_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  new_mat[survive_mask | born_mask] <- 1
  new_mat
}

game_of_life_rule <- function() {
  list(
    born = 3,
    survive = c(2, 3)
  )
}


simulate_ca <- function(initial_mat, n_steps) {
  history <- vector("list", n_steps + 1)
  history[[1]] <- initial_mat
  current <- initial_mat

  parsed <- game_of_life_rule()
  born <- parsed$born
  survive <- parsed$survive

  for (t in 1:n_steps) {
    current <- next_generation_generic(current, born, survive)

    history[[t + 1]] <- current
  }

  history
}

make_random_init <- function(n_rows, n_cols, p_alive) {
  matrix(
    sample(
      c(0, 1),
      n_rows * n_cols,
      replace = TRUE,
      prob = c(1 - p_alive, p_alive)
    ),
    nrow = n_rows,
    ncol = n_cols
  )
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

run_gol_hxc_analysis_selfcontained <- function(
  seed,
  k,
  T,
  emb,
  window_size = 4,
  stride = 4
) {
  n <- 2^k # grid dimension

  set.seed(seed)
  map <- matrix(sample(c(0, 1), n * n, replace = TRUE), nrow = n, ncol = n)

  # window start positions: 1, 1+stride, ... such that window fits within n
  starts <- seq(1, n - window_size + 1, by = stride)

  maps <- vector("list", T)
  time_series_list <- vector("list", T)
  hc_results <- vector("list", T)

  for (t in 1:T) {
    # map_t from map_{t-1}
    map <- next_generation(map)

    time.series <- numeric(0)

    for (i in starts) {
      for (j in starts) {
        window.binary <- as.vector(map[
          i:(i + window_size - 1),
          j:(j + window_size - 1)
        ])
        window.real <- binary_to_real(window.binary)
        time.series <- c(time.series, window.real)
      }
    }

    hc <- compute_hxc(time.series, emb = emb)

    maps[[t]] <- map
    time_series_list[[t]] <- time.series
    hc_results[[t]] <- hc
  }

  list(
    maps = maps,
    time_series = time_series_list,
    hc_results = hc_results
  )
}

run_gol_hxc_analysis_from_simulation <- function(
  maps,
  emb,
  window_size,
  stride,
  seed
) {
  set.seed(seed)

  T <- length(maps)
  n <- nrow(maps[[1]])

  starts <- seq(1, n - window_size + 1, by = stride)

  time_series_list <- vector("list", T)
  hc_results <- vector("list", T)

  for (t in 1:T) {
    map <- maps[[t]]

    time.series <- numeric(0)
    for (i in starts) {
      for (j in starts) {
        window.binary <- as.vector(map[
          i:(i + window_size - 1),
          j:(j + window_size - 1)
        ])
        window.real <- binary_to_real(window.binary)
        time.series <- c(time.series, window.real)
      }
    }

    hc <- compute_hxc(time.series, emb = emb)

    time_series_list[[t]] <- time.series
    hc_results[[t]] <- hc
  }

  list(
    maps = maps,
    time_series = time_series_list,
    hc_results = hc_results
  )
}


# Peano-Hilbert trajectory

# Convert Hilbert curve distance d to (x, y) for an n x n grid.
# n must be a power of 2.
hilbert_d2xy <- function(n, d) {
  rx <- 0L
  ry <- 0L
  t <- d
  x <- 0L
  y <- 0L
  s <- 1L
  while (s < n) {
    rx <- bitwAnd(bitwShiftR(t, 1L), 1L)
    ry <- bitwAnd(bitwXor(t, rx), 1L)
    if (ry == 0L) {
      if (rx == 1L) {
        x <- s - 1L - x
        y <- s - 1L - y
      }
      tmp <- x
      x <- y
      y <- tmp
    }
    x <- x + s * rx
    y <- y + s * ry
    t <- bitwShiftR(t, 2L)
    s <- s * 2L
  }
  c(x = x, y = y)
}

# Full traversal order for an n x n grid as a matrix of 1-indexed
# (row, col) pairs, ordered along the Hilbert curve (d = 0 .. n^2-1).
hilbert_curve_order <- function(n) {
  stopifnot(bitwAnd(n, n - 1L) == 0L) # n must be a power of 2
  coords <- t(vapply(0:(n^2 - 1), function(d) hilbert_d2xy(n, d), numeric(2)))
  cbind(row = coords[, 2] + 1L, col = coords[, 1] + 1L)
}

run_gol_hxc_analysis_from_simulation_pht <- function(
  maps,
  emb,
  group_size = 16,
  seed
) {
  set.seed(seed)
  T <- length(maps)
  n <- nrow(maps[[1]])

  order_coords <- hilbert_curve_order(n) # computed once, same for every t
  n_groups <- nrow(order_coords) %/% group_size
  stopifnot(nrow(order_coords) %% group_size == 0)

  time_series_list <- vector("list", T)
  hc_results <- vector("list", T)

  for (t in 1:T) {
    map <- maps[[t]]

    # read cell values in Hilbert-curve order
    vals <- map[order_coords]

    # every consecutive 16 values (along the curve) -> one real number
    time.series <- numeric(n_groups)
    for (g in seq_len(n_groups)) {
      idx <- ((g - 1) * group_size + 1):(g * group_size)
      time.series[g] <- binary_to_real(vals[idx])
    }

    hc <- compute_hxc(time.series, emb = emb)

    time_series_list[[t]] <- time.series
    hc_results[[t]] <- hc
  }

  list(
    maps = maps,
    time_series = time_series_list,
    hc_results = hc_results
  )
}

plot_time_series_epoch <- function(result, hc_df, epoch, title = NULL) {
  ts <- result$time_series[[epoch]]

  df <- data.frame(
    index = seq_along(ts),
    value = ts
  )

  # pull the matching row from hc_df
  row <- hc_df[hc_df$epoch == epoch, ]

  if (nrow(row) == 0) {
    stop(sprintf("No entry found in hc_df for epoch %d", epoch))
  }

  if (is.null(title)) {
    title <- sprintf("Time series - Epoch %d", epoch)
  }

  subtitle <- sprintf(
    "H = %.3f | C = %.3f | velocity = %s | dist_to_limit = %.2f",
    row$H,
    row$C,
    ifelse(is.na(row$velocity), "NA", formatC(row$velocity, digits = 2)),
    row$dist_to_limit
  )

  ggplot(df, aes(x = index, y = value)) +
    geom_line(color = "steelblue", linewidth = 0.4) +
    geom_point(color = "steelblue", size = 0.8, alpha = 0.6) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Window index (Hilbert order)",
      y = "Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      panel.grid.minor = element_blank()
    )
}

plot_time_series_epoch_faceted <- function(result, hc_df, epoch, n_panels = 4) {
  ts <- result$time_series[[epoch]]
  n <- length(ts)
  panel <- cut(seq_along(ts), breaks = n_panels, labels = FALSE)

  df <- data.frame(index = seq_along(ts), value = ts, panel = factor(panel))

  ggplot(df, aes(x = index, y = value)) +
    geom_line(color = "steelblue", linewidth = 0.25) +
    facet_wrap(~panel, scales = "free_x", ncol = 1) +
    labs(
      title = sprintf("Time series - Epoch %d", epoch),
      x = "Window index (Hilbert order)",
      y = "binary_to_real value"
    ) +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_blank(), panel.grid.minor = element_blank())
}

analyze_convergence <- function(hc_df) {
  df <- hc_df %>%
    filter(!is.na(velocity)) %>%
    mutate(epoch = as.numeric(as.character(epoch)))

  # sanity check antes de seguir
  stopifnot(is.numeric(df$epoch))
  if (any(is.na(df$epoch))) {
    warning(
      "Some epoch values could not be coerced to numeric — check hc_df$epoch"
    )
  }

  # --- Plot 1: velocity vs epoch ---
  p_velocity <- ggplot(df, aes(x = epoch, y = velocity)) +
    geom_line(color = "steelblue", linewidth = 0.4) +
    geom_point(color = "steelblue", size = 0.8, alpha = 0.6) +
    scale_y_log10() +
    labs(
      title = "Convergence velocity over epochs",
      x = "Epoch",
      y = "Velocity (log scale)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  # --- Plot 2: dist_to_limit vs epoch ---
  p_dist <- ggplot(df, aes(x = epoch, y = dist_to_limit)) +
    geom_line(color = "darkorange", linewidth = 0.4) +
    geom_point(color = "darkorange", size = 0.8, alpha = 0.6) +
    scale_y_log10() +
    labs(
      title = "Distance to H×C boundary over epochs",
      x = "Epoch",
      y = "dist_to_limit (log scale)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  # --- Regressões log-lineares ---
  fit_dist <- lm(log(dist_to_limit) ~ epoch, data = df)
  fit_vel <- lm(log(velocity) ~ epoch, data = df)

  cat("=== dist_to_limit decay fit ===\n")
  print(summary(fit_dist))
  cat("\nDecay rate k (dist_to_limit):", -coef(fit_dist)[["epoch"]], "\n")
  cat("Half-life (epochs):", log(2) / -coef(fit_dist)[["epoch"]], "\n\n")

  cat("=== velocity decay fit ===\n")
  print(summary(fit_vel))
  cat("\nDecay rate k (velocity):", -coef(fit_vel)[["epoch"]], "\n")

  df$dist_fit <- exp(predict(fit_dist))
  p_dist <- p_dist +
    geom_line(
      data = df,
      aes(y = dist_fit),
      color = "black",
      linetype = "dashed",
      linewidth = 0.5
    )

  df$vel_fit <- exp(predict(fit_vel))
  p_velocity <- p_velocity +
    geom_line(
      data = df,
      aes(y = vel_fit),
      color = "black",
      linetype = "dashed",
      linewidth = 0.5
    )

  combined <- p_velocity /
    p_dist +
    plot_annotation(title = "GoL H×C convergence dynamics")

  list(
    plot = combined,
    fit_dist = fit_dist,
    fit_velocity = fit_vel
  )
}

diagnose_convergence_residuals <- function(
  hc_df,
  fit_dist,
  fit_velocity,
  outlier_z = 2
) {
  df <- hc_df %>%
    filter(!is.na(velocity)) %>%
    mutate(epoch = as.numeric(as.character(epoch)))

  # resíduos (alinhados por linha, já que ambos os fits usam o mesmo df filtrado)
  df$resid_dist <- residuals(fit_dist)
  df$resid_vel <- residuals(fit_velocity)

  # padroniza os resíduos (z-score) pra definir outliers de forma comparável
  df$z_dist <- scale(df$resid_dist)[, 1]
  df$z_vel <- scale(df$resid_vel)[, 1]

  df$outlier_dist <- abs(df$z_dist) > outlier_z
  df$outlier_vel <- abs(df$z_vel) > outlier_z
  df$outlier_both <- df$outlier_dist & df$outlier_vel

  # --- Plot 1: resíduos de dist_to_limit ao longo das épocas ---
  p1 <- ggplot(df, aes(x = epoch, y = resid_dist)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
    geom_point(aes(color = outlier_dist), size = 1, alpha = 0.6) +
    scale_color_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "firebrick")
    ) +
    labs(
      title = "Residuals: log(dist_to_limit) ~ epoch",
      x = "Epoch",
      y = "Residual",
      color = "Outlier"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")

  # --- Plot 2: resíduos de velocity ao longo das épocas ---
  p2 <- ggplot(df, aes(x = epoch, y = resid_vel)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
    geom_point(aes(color = outlier_vel), size = 1, alpha = 0.6) +
    scale_color_manual(
      values = c("FALSE" = "darkorange", "TRUE" = "firebrick")
    ) +
    labs(
      title = "Residuals: log(velocity) ~ epoch",
      x = "Epoch",
      y = "Residual",
      color = "Outlier"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")

  # --- Plot 3: sobreposição - épocas outlier em ambos os modelos ---
  p3 <- ggplot(df, aes(x = epoch, y = dist_to_limit)) +
    geom_line(color = "grey70", linewidth = 0.3) +
    geom_point(
      data = filter(df, outlier_dist & !outlier_vel),
      aes(y = dist_to_limit),
      color = "steelblue",
      size = 1.5
    ) +
    geom_point(
      data = filter(df, outlier_vel & !outlier_dist),
      aes(y = dist_to_limit),
      color = "darkorange",
      size = 1.5
    ) +
    geom_point(
      data = filter(df, outlier_both),
      aes(y = dist_to_limit),
      color = "firebrick",
      size = 2
    ) +
    labs(
      title = "dist_to_limit trajectory with flagged reorganization epochs",
      subtitle = "blue = dist outlier only | orange = velocity outlier only | red = both",
      x = "Epoch",
      y = "dist_to_limit"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )

  combined <- (p1 / p2) | p3

  # tabela resumo
  n_dist_outliers <- sum(df$outlier_dist)
  n_vel_outliers <- sum(df$outlier_vel)
  n_both <- sum(df$outlier_both)

  cat("=== Outlier summary (|z| >", outlier_z, ") ===\n")
  cat("dist_to_limit outliers:", n_dist_outliers, "\n")
  cat("velocity outliers:     ", n_vel_outliers, "\n")
  cat("overlap (both):        ", n_both, "\n")
  cat(
    "overlap rate (of dist outliers):",
    round(100 * n_both / max(n_dist_outliers, 1), 1),
    "%\n"
  )
  cat(
    "overlap rate (of velocity outliers):",
    round(100 * n_both / max(n_vel_outliers, 1), 1),
    "%\n\n"
  )

  # teste simples de associação: outliers coincidem mais do que o esperado ao acaso?
  tab <- table(df$outlier_dist, df$outlier_vel)
  print(tab)
  cat("\nFisher's exact test (association between outlier flags):\n")
  print(fisher.test(tab))

  list(
    df = df,
    plot = combined,
    outliers_both = filter(df, outlier_both) %>%
      select(epoch, H, C, velocity, dist_to_limit)
  )
}

analyze_ccf_residuals <- function(diag_df, max_lag = 20) {
  df <- diag_df %>% arrange(epoch)

  # checa se a série de epochs é regular (sem buracos), pré-requisito pro ccf ter sentido
  gaps <- diff(df$epoch)
  if (any(gaps != 1)) {
    warning(
      "epoch não é sequencial/regular — CCF pode não ser interpretável diretamente. ",
      "Considere reindexar ou interpolar."
    )
  }

  resid_dist_ts <- ts(df$resid_dist)
  resid_vel_ts <- ts(df$resid_vel)

  # --- CCF entre os resíduos ---
  ccf_result <- ccf(
    resid_vel_ts,
    resid_dist_ts,
    lag.max = max_lag,
    plot = FALSE
  )

  ccf_df <- data.frame(
    lag = ccf_result$lag[, 1, 1],
    acf = ccf_result$acf[, 1, 1]
  )

  # limite de significância aproximado (95%, ~ ±1.96/sqrt(n))
  n <- length(resid_dist_ts)
  ci <- 1.96 / sqrt(n)

  p_ccf <- ggplot(ccf_df, aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_hline(
      yintercept = c(-ci, ci),
      color = "steelblue",
      linetype = "dashed"
    ) +
    geom_segment(
      aes(xend = lag, yend = 0),
      color = "darkorange",
      linewidth = 0.6
    ) +
    geom_point(color = "darkorange", size = 1.5) +
    labs(
      title = "Cross-correlation: velocity residuals vs dist_to_limit residuals",
      subtitle = "Positive lag: velocity leads dist_to_limit | dashed lines = 95% significance",
      x = "Lag (epochs)",
      y = "Cross-correlation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )

  # lags significativos
  sig_lags <- ccf_df %>% filter(abs(acf) > ci)

  cat("=== CCF summary ===\n")
  cat("n observations:", n, "\n")
  cat("95% significance threshold: ±", round(ci, 4), "\n\n")

  if (nrow(sig_lags) > 0) {
    cat("Significant lags found:\n")
    print(sig_lags[order(-abs(sig_lags$acf)), ])
    cat("\nInterpretation guide:\n")
    cat(
      "  lag > 0: velocity residual at t predicts dist_to_limit residual at t+lag\n"
    )
    cat(
      "  lag < 0: dist_to_limit residual at t predicts velocity residual at t+|lag|\n"
    )
  } else {
    cat(
      "Nenhum lag ultrapassou o limite de significância — não há evidência de\n"
    )
    cat("relação defasada linear entre os resíduos das duas séries.\n")
  }

  list(
    ccf_df = ccf_df,
    plot = p_ccf,
    sig_lags = sig_lags,
    n = n,
    ci = ci
  )
}

plot_dist_to_limit <- function(hc_df, fit_dist = NULL, log_scale = FALSE) {
  df <- hc_df %>%
    mutate(epoch = as.numeric(as.character(epoch)))

  p <- ggplot(df, aes(x = epoch, y = dist_to_limit)) +
    geom_point(size = 1.6, alpha = 0.5, color = "grey20")

  if (!is.null(fit_dist)) {
    pred <- predict(fit_dist, newdata = df, se.fit = TRUE)
    df$fit <- exp(pred$fit)
    df$fit_lo <- exp(pred$fit - 1.96 * pred$se.fit)
    df$fit_hi <- exp(pred$fit + 1.96 * pred$se.fit)

    p <- p +
      geom_ribbon(
        data = df,
        aes(y = fit, ymin = fit_lo, ymax = fit_hi),
        fill = "firebrick",
        alpha = 0.15,
        color = NA
      ) +
      geom_line(data = df, aes(y = fit), color = "firebrick", linewidth = 0.8)
  } else {
    p <- p +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "firebrick",
        linewidth = 0.8,
        fill = "firebrick",
        alpha = 0.15
      )
  }

  p <- p +
    labs(x = "Epoch", y = expression(d[to ~ limit])) +
    theme_tufte() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 10),
      axis.title = element_text(size = 15),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.2)
    ) +
    geom_vline(xintercept = 688, linetype = "dotted", color = "grey40") +
    annotate(
      "text",
      x = 688,
      y = max(df$dist_to_limit),
      label = "half-life",
      angle = 90,
      vjust = -0.5,
      size = 3,
      color = "grey40"
    )

  if (log_scale) {
    p <- p + scale_y_log10()
  }

  p
}

# ---------------------------------------------------------------
# 1. Find the smallest tau for which the final tau epochs are
#    stationary (ADF test), searching over a grid of candidate taus.
# ---------------------------------------------------------------
find_tau_stationary <- function(
  hc_df,
  tau_grid = seq(300, 1500, by = 100),
  alpha = 0.05,
  vars = c("H", "C")
) {
  hc_df <- hc_df %>% mutate(epoch = as.numeric(as.character(epoch)))
  T_max <- max(hc_df$epoch)

  results <- lapply(tau_grid, function(tau) {
    window_df <- hc_df %>% filter(epoch > T_max - tau)

    # ADF test needs a reasonable number of observations;
    # skip windows too short to test reliably
    if (nrow(window_df) < 15) {
      return(data.frame(tau = tau, p_H = NA, p_C = NA, stationary = NA))
    }

    p_vals <- sapply(vars, function(v) {
      tryCatch(
        adf.test(window_df[[v]], alternative = "stationary")$p.value,
        error = function(e) NA
      )
    })

    data.frame(
      tau = tau,
      p_H = p_vals[1],
      p_C = p_vals[2],
      stationary = all(p_vals < alpha, na.rm = TRUE)
    )
  }) %>%
    bind_rows()

  # smallest tau for which BOTH H and C pass the stationarity test
  tau_star <- results %>%
    filter(stationary) %>%
    slice_min(tau, n = 1) %>%
    pull(tau)

  list(tau_star = tau_star, table = results)
}
# ---------------------------------------------------------------
# 2. Fit the exponential-decay model for a given tau and return
#    r, half-life, and the fit object.
# ---------------------------------------------------------------
fit_decay_for_tau <- function(hc_df, tau) {
  hc_df <- hc_df %>% mutate(epoch = as.numeric(as.character(epoch)))
  T_max <- max(hc_df$epoch)

  asym <- hc_df %>%
    filter(epoch > T_max - tau) %>%
    summarise(H_star = mean(H), C_star = mean(C))

  df <- hc_df %>%
    mutate(dist_to_limit = sqrt((H - asym$H_star)^2 + (C - asym$C_star)^2))

  fit <- lm(log(dist_to_limit) ~ epoch, data = df)

  r <- -coef(fit)[["epoch"]]
  half_life <- log(2) / r

  list(
    tau = tau,
    r = r,
    half_life = half_life,
    fit = fit,
    df = df,
    H_star = asym$H_star,
    C_star = asym$C_star
  )
}
# ---------------------------------------------------------------
# 3. Sensitivity sweep: refit across a range of tau and report
#    how much r / half-life move.
# ---------------------------------------------------------------
sensitivity_sweep <- function(hc_df, tau_values = c(50, 100, 200, 400)) {
  results <- lapply(tau_values, function(tau) {
    res <- fit_decay_for_tau(hc_df, tau)
    data.frame(tau = tau, r = res$r, half_life = res$half_life)
  }) %>%
    bind_rows()

  results <- results %>%
    mutate(
      pct_dev_r = 100 * (r - mean(r)) / mean(r),
      pct_dev_hl = 100 * (half_life - mean(half_life)) / mean(half_life)
    )

  results
}


# Building video

results_CAsimulation <- ca_history3k_0.8

dir.create("frames_3k_0-8", showWarnings = FALSE)

for (t in seq_along(results_CAsimulation)) {
  mat <- results_CAsimulation[[t]]
  img <- array(1 - mat, dim = c(nrow(mat), ncol(mat), 1))
  writePNG(img, target = sprintf("frames_3k_0-8/frame_%04d.png", t - 1))
}

# Debuging

mat <- ca_history3k_0.8[[1]]

pct_alive <- 100 * mean(mat == 1)
cat(sprintf(
  "Alive cells: %.4f%% (%d out of %d)\n",
  pct_alive,
  sum(mat == 1),
  length(mat)
))
