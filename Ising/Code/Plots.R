
source(here("Code", "Functions.R"))

# um método
d <- read_hxc("MH", 1, emb = 6, bits = 16, )

# os três métodos juntos; a coluna method só existe em memória, para o gráfico
d <- do.call(rbind, lapply(methods, function(m)
  cbind(method = m, read_hxc(m, 0.5, 4, 16))))

ggplot(d, aes(H, C, colour = generation)) +
  geom_path() +
  facet_wrap(~ method) +
  scale_colour_viridis_c() +
  theme_few()




library(readr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(here)

method_names <- c(MH = "Metropolis", GB = "Gibbs", SW = "Swendsen-Wang")

plot_hxc_methods <- function(beta, emb, bits,
                             methods = c("MH", "GB", "SW"),
                             label_gens = c(0, 5, 50, 100),
                             root = out_dir) {
  
  # read the same configuration for each method and stack them
  hc_df <- do.call(rbind, lapply(methods, function(m) {
    file <- file.path(root, m, sprintf("beta%.4f", beta),
                      sprintf("emb%d_bits%d.csv", emb, bits))
    d <- as.data.frame(read_csv(file, show_col_types = FALSE))
    d$method <- m
    d
  }))
  
  hc_df$H <- as.numeric(hc_df$H)
  hc_df$C <- as.numeric(hc_df$C)
  hc_df <- hc_df[!is.na(hc_df$H) & !is.na(hc_df$C), ]
  
  hc_df$t      <- hc_df$generation
  hc_df$epoch  <- factor(hc_df$generation)
  hc_df$method <- factor(hc_df$method, levels = methods,
                         labels = method_names[methods])
  
  hc_df_labels <- subset(hc_df, t %in% label_gens)
  max_gen      <- max(hc_df$generation)
  
  ggplot(hc_df, aes(x = H, y = C, color = epoch)) +
    geom_line(data = subset(LinfLsup, Dimension == emb & Side == "Lower"),
              aes(x = H, y = C), color = "grey40", linewidth = 0.6,
              inherit.aes = FALSE) +
    geom_line(data = subset(LinfLsup, Dimension == emb & Side == "Upper"),
              aes(x = H, y = C), color = "grey40", linewidth = 0.6,
              inherit.aes = FALSE) +
    geom_point(size = 1, alpha = 0.85) +
    geom_label_repel(
      data = hc_df_labels,
      aes(label = t),
      size = 4,
      fontface = "bold",
      box.padding = 0.25,
      point.padding = 0.1,
      segment.color = "grey60",
      segment.size = 0.4,
      seed = 42,
      nudge_y = 0.05
    ) +
    facet_wrap(~ method, nrow = 1) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 0.5)) +
    labs(
      title = sprintf("%d generations  |  beta = %.4f  |  dim = %d  |  bits = %d",
                      max_gen, beta, emb, bits),
      x = "Permutation Entropy (H)",
      y = "Statistical Complexity (C)") +
    theme_tufte() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          strip.text   = element_text(size = 14, face = "bold"),
          axis.title   = element_text(size = 15),
          axis.text.y  = element_text(size = 12),
          axis.text.x  = element_text(size = 12),
          panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
          panel.grid.minor = element_line(color = "grey92", linewidth = 0.2))
}

p <- plot_hxc("MH", beta = 0.5, emb = 4, bits = 16)


# save one figure per method and beta
dir.create(here("Figures"), showWarnings = FALSE)
for (b in beta_vals) {
  p <- plot_hxc_methods(b, emb = 4, bits = 16)
  ggsave(here("ising", "figures",
              sprintf("hxc_methods_beta%.4f_emb4_bits16.png", b)),
         p, width = 15, height = 5, dpi = 300, bg = "white")
}
