library(ggplot2)
library(ggthemes)

# ============================================================
# Game of Life -  HxC Analysis
# ============================================================

emb_used <- 3  # Embedding dimention
bits <- 32
seed <- 12345678
K <- 500  # K = Passo para plot (se K=500 a cada 500 épocas, coloca label no plot)
N_tail <- 300 # número de épocas para cáculo de (H*,C*)

# Using Peano-Hilbert trajectory
hxc_out <-run_gol_hxc_analysis_from_simulation_pht(
  results_CAsimulation, 
  emb = emb_used, 
  group_size = bits, 
  seed = seed)

cat("Number of epochs simulated:", length(hxc_out$maps), "\n")
cat("Length of time series per epoch:", length(hxc_out$time_series[[1]]), "\n")
cat("H, C for epoch 1:\n")
print(hxc_out$hc_results[[1]])

# Build a data frame of H and C across all epochs, ready for plotting
# on the complexity-entropy plane (same style as the station-level plot):
hc_df <- data.frame(do.call(rbind, hxc_out$hc_results))

hc_df$epoch <- 1:nrow(hc_df)
hc_df$epoch <- factor(hc_df$epoch)  # discrete labels/colors, one per epoch
hc_df$H <- as.numeric(hc_df$H)
hc_df$C <- as.numeric(hc_df$C)
hc_df <- hc_df[!is.na(hc_df$H) & !is.na(hc_df$C), ] # remover linhas em que H ou C são NA

# ---------------------------------------------------------------
# Plotting HxC
# ---------------------------------------------------------------

hc_df$t <- as.numeric(as.character(hc_df$epoch))

# Always include first and last epoch, plus every K-th one, for labeling
label_epochs <- unique(c(min(hc_df$t), 
                         hc_df$t[hc_df$t %% K == 0], 
                         max(hc_df$t)))
hc_df_labels <- subset(hc_df, t %in% label_epochs)

ggplot(hc_df, aes(x = H, y = C, color = epoch)) +
  geom_line(data = subset(LinfLsup, Dimension == emb_used & Side == "Lower"),
            aes(x = H, y = C), color = "grey40", linewidth = 0.6,
            inherit.aes = FALSE) +
  geom_line(data = subset(LinfLsup, Dimension == emb_used & Side == "Upper"),
            aes(x = H, y = C), color = "grey40", linewidth = 0.6,
            inherit.aes = FALSE) +
  geom_point(size = 1, alpha = 0.85) +
  geom_text_repel(data = hc_df_labels,
                  aes(label = t),
                  size = 4,
                  fontface = "bold",
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  point.padding = 0.4,
                  segment.color = "grey60",
                  segment.size = 0.4,
                  seed = 42) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(
       x = "Permutation Entropy (H)",
       y = "Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey92", linewidth = 0.2))


# ---------------------------------------------------------------
# Plotting Time Series
# ---------------------------------------------------------------
# pode alterar epoch para ver a série temporal gerada na época
plot_time_series_epoch(hxc_out, hc_df, epoch = 3000)
plot_time_series_epoch_faceted(hxc_out, hc_df, epoch = 150, n_panels = 4) 

