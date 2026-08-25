epoch <- 10
p_alive <- 0.5
embdim <- 6
bits <- 32
i<-2
N_tail <- 1000

file <- sprintf("Data/results/%sK/new_results/emb6/results_%sk_%sbits_%s_%s.csv", epoch, epoch, bits, p_alive, i)

hc_df <- read_csv(file, show_col_types = FALSE)

# ---------------------------------------------------------------
# Velocity of convergence in the HxC plane
# ---------------------------------------------------------------

hc_df <- hc_df[order(as.numeric(as.character(hc_df$step))), ]
hc_df$H <- as.numeric(hc_df$H)
hc_df$C <- as.numeric(hc_df$C)
hc_df$t <- as.numeric(as.character(hc_df$step))

# 1) Step-wise velocity: Euclidean displacement per epoch in (H, C) space
hc_df$velocity <- c(NA, sqrt(diff(hc_df$H)^2 + diff(hc_df$C)^2))

# 2) Estimate the limit point (H*, C*) from the tail of the run
#    (only valid if the system has actually settled to a fixed point by then)

H_star <- mean(tail(hc_df$H, N_tail))
C_star <- mean(tail(hc_df$C, N_tail))
cat("(H*, C*) =", H_star, C_star, "\n")

hc_df$dist_to_limit <- sqrt((hc_df$H - H_star)^2 + (hc_df$C - C_star)^2)

# 3) Fit exponential decay: d(t) ~ d0 * exp(-r*t)  ->  log(d) = log(d0) - r*t
conv_data <- subset(hc_df, dist_to_limit > 1e-8)
#write.csv(conv_data, file = "conv_data.csv")
fit <- lm(log(dist_to_limit) ~ t, data = conv_data)
r_rate <- -coef(fit)[["t"]]

cat("Convergence rate r =", r_rate, "\n")
cat("Half-life (epochs) =", log(2) / r_rate, "\n")
summary(fit)$r.squared

#title = "Convergence to Asymptotic State", subtitle = "Log-distance to estimated limit point (H*, C*) over time",
ggplot(conv_data, aes(x = t, y = log(dist_to_limit))) +
  geom_point(size = 1.6, alpha = 0.5, color = "grey20") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", 
              linewidth = 0.8, fill = "firebrick", alpha = 0.15) +
  labs(
      title = sprintf("%dk generations | p(alive) = %d%%  |  dim = %d  |  bits = %d",epoch, p_alive*100, embdim, bits),
      x = "Epoch",
       y = expression(log(d[to~limit]))) +
  theme_tufte() +
  theme(plot.title = element_text(size = 14),
        plot.subtitle = element_text(color = "grey40", size = 10),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey94", linewidth = 0.2))



# Convergence analysis
result_analysis <- analyze_convergence(hc_df)
result_analysis$plot          # mostra os 2 painéis empilhados
summary(result_analysis$fit_dist)
summary(result_analysis$fit_velocity)

diag <- diagnose_convergence_residuals(hc_df, result_analysis$fit_dist, result_analysis$fit_velocity)
diag$plot
diag$outliers_both  

ccf_analysis <- analyze_ccf_residuals(diag$df, max_lag = 20)
ccf_analysis$plot
ccf_analysis$sig_lags

# Plot distance
plot_dist_to_limit <- function(hc_df, fit_dist = NULL, log_scale = FALSE) {
  df <- hc_df %>%
    mutate(step = as.numeric(as.character(step)))
  
  p <- ggplot(df, aes(x = step, y = dist_to_limit)) +
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


plot_dist_to_limit(hc_df)
