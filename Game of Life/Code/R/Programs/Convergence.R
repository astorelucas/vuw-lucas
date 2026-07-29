# ---------------------------------------------------------------
# Velocity of convergence in the HxC plane
# ---------------------------------------------------------------

hc_df <- hc_df[order(as.numeric(as.character(hc_df$epoch))), ]
hc_df$H <- as.numeric(hc_df$H)
hc_df$C <- as.numeric(hc_df$C)
hc_df$t <- as.numeric(as.character(hc_df$epoch))

# 1) Step-wise velocity: Euclidean displacement per epoch in (H, C) space
hc_df$velocity <- c(NA, sqrt(diff(hc_df$H)^2 + diff(hc_df$C)^2))

# 2) Estimate the limit point (H*, C*) from the tail of the run
#    (only valid if the system has actually settled to a fixed point by then)

H_star <- mean(tail(hc_df$H, N_tail))
C_star <- mean(tail(hc_df$C, N_tail))

hc_df$dist_to_limit <- sqrt((hc_df$H - H_star)^2 + (hc_df$C - C_star)^2)

# 3) Fit exponential decay: d(t) ~ d0 * exp(-r*t)  ->  log(d) = log(d0) - r*t
conv_data <- subset(hc_df, dist_to_limit > 1e-8)
write.csv(conv_data, file = "conv_data.csv")
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
  labs(x = "Epoch",
       y = expression(log(d[to~limit]))) +
  theme_tufte() +
  theme(plot.title = element_text(face = "bold", size = 14),
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

plot_dist_to_limit(hc_df)