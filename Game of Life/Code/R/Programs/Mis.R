
hxc_out <- run_gol_hxc_analysis(
  seed = 123456789,
  k    = 8,    # 2^8 = 256x256 grid
  T    = 3000,
  emb  = emb_used     # embedding dimension passed to OPprob()
)

# --- 1. Fit only the upper envelope (local maxima) for a cleaner rate ---
library(zoo)
win <- 15  # rolling window for local-max detection; tune to taste
is_local_max <- conv_data$dist_to_limit == 
  rollapply(conv_data$dist_to_limit, width = win, FUN = max, 
            align = "center", fill = NA)
envelope_data <- conv_data[which(is_local_max), ]

fit_envelope <- lm(log(dist_to_limit) ~ t, data = envelope_data)
r_rate_envelope <- -coef(fit_envelope)[["t"]]
cat("Envelope-based convergence rate:", r_rate_envelope, "\n")
cat("Envelope half-life:", log(2) / r_rate_envelope, "\n")

# --- 2. Quantify the oscillation period via autocorrelation of the residual ---
resid_series <- resid(fit)  # detrended log-distance
acf_result <- acf(resid_series, lag.max = 500, plot = TRUE)
# look for the first strong peak after lag 0 -> that's your oscillator period

# --- 3. Visual: overlay envelope fit on the original plot ---
plot(conv_data$t, log(conv_data$dist_to_limit), pch = 1, cex = 0.5,
     xlab = "Epoch", ylab = "log(dist_to_limit)")
abline(fit, col = "red", lwd = 1)                  # original overall fit
abline(fit_envelope, col = "blue", lwd = 2)        # envelope-only fit
legend("topright", legend = c("all points", "upper envelope"),
       col = c("red", "blue"), lwd = c(1, 2))

# Tau sensivity

# # Step 1: find tau via stationarity criterion
stat_result <- find_tau_stationary(hc_df)
tau_star <- stat_result$tau_star

# Step 2: fit at tau_star
main_fit <- fit_decay_for_tau(hc_df, tau_star)
main_fit$r
main_fit$half_life

# Step 3: sensitivity sweep to report in the text
sens <- sensitivity_sweep(hc_df, tau_values = c(300, 400, 500,600, 700, 800, 900, tau_star))
sens
max(abs(sens$pct_dev_hl))   

ggplot(stat_result$table, aes(x = tau)) +
  geom_line(aes(y = p_H, color = "p (H)"), linewidth = 0.7) +
  geom_line(aes(y = p_C, color = "p (C)"), linewidth = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 400, linetype = "dotted", color = "grey40") +
  annotate("text", x = 400, y = 0.75, label = expression(tau*"* = 400"),
           angle = 90, vjust = -0.5, size = 3, color = "grey40") +
  scale_color_manual(values = c("p (H)" = "#2a78d6", "p (C)" = "#e34948")) +
  labs(x = expression(tau~"(epochs)"), y = "ADF p-value", color = NULL) +
  theme_tufte(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.major = element_line(color = "grey88", linewidth = 0.3))

