library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(patchwork)

## ----------------------------------------
## upsilon estimation - results
## ----------------------------------------

epoch <- 10
p_alive <- 0.5 #0.3 and 0.5 was evaluated
embdim <- 6
bits <- 8
i <- 2

N_tail_values <- c(1000, 2000, 4000, 6000, 8000)

file <- sprintf("Data/results/%sK/thimming results/results_%sk_%s_emb%s_%sbits.csv",
                epoch, epoch, p_alive, embdim, bits)

hc_df <- read_csv(file, show_col_types = FALSE)

hc_df <- hc_df[order(as.numeric(as.character(hc_df$step))), ]
hc_df$H <- as.numeric(hc_df$H)
hc_df$C <- as.numeric(hc_df$C)
hc_df$t <- as.numeric(as.character(hc_df$step))

hc_df$velocity <- c(NA, sqrt(diff(hc_df$H)^2 + diff(hc_df$C)^2))

estimate_convergence <- function(hc_df, N_tail) {
  
  n_rows <- nrow(hc_df)
  if (N_tail > n_rows) {
    warning(sprintf(
      "N_tail = %d > nrow(hc_df) = %d; using all %d rows instead.",
      N_tail, n_rows, n_rows
    ))
    N_tail <- n_rows
  }
  
  H_star <- mean(tail(hc_df$H, N_tail))
  C_star <- mean(tail(hc_df$C, N_tail))
  
  dist_to_limit <- sqrt((hc_df$H - H_star)^2 + (hc_df$C - C_star)^2)
  df <- hc_df %>% mutate(dist_to_limit = dist_to_limit)
  
  conv_data <- subset(df, dist_to_limit > 1e-8)
  fit <- lm(log(dist_to_limit) ~ t, data = conv_data)
  
  r_rate <- -coef(fit)[["t"]]
  half_life <- log(2) / r_rate
  r_squared <- summary(fit)$r.squared
  
  tibble(
    N_tail = N_tail,
    H_star = H_star,
    C_star = C_star,
    r_rate = r_rate,
    half_life = half_life,
    r_squared = r_squared
  )
}

results_table <- map_dfr(N_tail_values, ~ estimate_convergence(hc_df, .x))

print(results_table)

# Uncomment to save:
# write_csv(results_table, "ntail_sensitivity_results.csv")

## ----------------------------------------
## Graph Final(H) vs Final(C)
## ----------------------------------------

embeddings <- c(5)
bits <- 16
epochs <- 3
p_alive <- 0.8

## Gerar 3 imagens
summary_results <- tibble()
epochs <- c(3, 5, 10)
embeddings <- c(3, 4, 5, 6)
bits <- c(2, 4, 8, 16, 32)
p_alive_vals <- c(0.1, 0.3, 0.5, 0.8)

for (ep in epochs) {
  for (pa in p_alive_vals) {
    for (emb in embeddings) {
      for (b in bits) {
        
        if (ep == 10) {
          if (pa == 0.8) {
            file <- sprintf(
              "Data/results/%sK/thimming results/results_%sk_%s_emb%s_%sbits_1.csv",
              ep, ep, pa, emb, b
            )
          } else {
            file <- sprintf(
              "Data/results/%sK/thimming results/results_%sk_%s_emb%s_%sbits.csv",
              ep, ep, pa, emb, b
            )
          }
        } else {
          file <- sprintf(
            "Data/results/%sK/results_%sk_%s_emb%s_%sbits.csv",
            ep, ep, pa, emb, b
          )
        }
        
        df <- read_csv(file, show_col_types = FALSE)
        
        tail_epochs <- if (ep == 3) 1800 else if (ep == 5) 3000 else 6000
        
        final_df <- df %>% slice_tail(n = tail_epochs)
        n <- nrow(final_df)
        
        summary_results <- bind_rows(
          summary_results,
          tibble(
            epoch = factor(ep),
            p_alive = factor(pa),
            embedding = factor(emb),
            bits = factor(b),
            
            mean_H = mean(final_df$H),
            ci_H = qt(0.975, n - 1) * sd(final_df$H) / sqrt(n),
            
            mean_C = mean(final_df$C),
            ci_C = qt(0.975, n - 1) * sd(final_df$C) / sqrt(n)
          )
        )
      }
    }
  }
}

ggplot(summary_results,
       aes(mean_H, mean_C, color = embedding)) +
  
  geom_path(aes(group = embedding),
            linewidth = 0.6,
            alpha = 0.7,
            linetype = "dashed") +
  
  geom_errorbar(aes(ymin = mean_C - ci_C,
                    ymax = mean_C + ci_C),
                width = 0.005) +
  
  geom_errorbarh(aes(xmin = mean_H - ci_H,
                     xmax = mean_H + ci_H),
                 height = 0.005) +
  
  geom_point(aes(shape = bits), size = 4) +
  
  facet_grid(p_alive ~ epoch) +
  theme_bw() +
  labs(
    x = "Final H",
    y = "Final C",
    color = "Embedding",
    shape = "Hilbert bits"
  )

#write.csv(summary_results, "hxc_mean_results.csv", row.names = FALSE)


# -------------------------------------------------------
# Graph Generations to reach H min
# -------------------------------------------------------


summary_results <- tibble()
epochs <- c(3, 5, 10)
embeddings <- c(3, 4, 5, 6)
bits <- c(2, 4, 8, 16, 32)
p_alive_vals <- c(0.1, 0.3, 0.5, 0.8)

for (ep in epochs) {
  for (pa in p_alive_vals) {
    for (emb in embeddings) {
      for (b in bits) {
        
        file <- sprintf(
          "Data/results/%sK/thimming results/results_%sk_%s_emb%s_%sbits.csv",
          ep, ep, pa, emb, b
        )
        
        df <- read_csv(file, show_col_types = FALSE)
        
        if ("step" %in% names(df)) {
          df <- rename(df, "epoch" = "step")
        }
        
        
        if (pa != 0.8){
          step_min <- which.min(df$H)
          minH <- min(df$H)
          
          step_max <- max(df$epoch)
          step_last <- df[df$epoch == step_max, "H"] %>% pull(H)  # valor do H no último passo
        }else{
          step_min <- which.max(df$H)
          minH <- max(df$H)
          
          step_max <- max(df$epoch)
          step_last <- df[df$epoch == step_max, "H"] %>% pull(H)  # valor do H no último passo
        }
        
        summary_results <- bind_rows(
          summary_results,
          tibble(
            epoch = factor(ep),
            p_alive = factor(pa),
            embedding = factor(emb),
            bits = factor(b),
            step_last = step_last,
            step_max = step_max,
            min_H = minH,
            step_min_H = step_min
          )
        )
      }
    }
  }
}

summary_results <- summary_results %>%
  mutate(
    ratio_min_H = step_min_H / step_max,          # quão perto do fim (1 = no último step)
    dist_to_end = step_max - step_min_H            # quantos steps antes do fim
  )

summary_results %>%
  filter(p_alive != 0.8) %>%
  ggplot(aes(x = epoch, y = step_min_H, fill = epoch)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  #geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  theme_bw() +
  labs(
    x = "Epoch (k)",
    y = "Generations until reaches minimum value",
    #title = "Distance (in steps) between H minimum and end of training"
  ) +
  theme(legend.position = "none")



# ---------------------------------------------------------------
# Build the combined data frame from your two sensitivity tables - from results
# ---------------------------------------------------------------

sens_df <- tibble::tribble(
  ~LC,   ~upsilon, ~H_star, ~C_star, ~r_rate,   ~half_life, ~r_squared,
  "30%", 1000,     0.284,   0.233,   0.000722,  961,        0.943,
  "30%", 2000,     0.284,   0.233,   0.000722,  960,        0.943,
  "30%", 4000,     0.285,   0.234,   0.000532,  1304,       0.906,
  "30%", 6000,     0.287,   0.234,   0.000457,  1516,       0.664,
  "30%", 8000,     0.291,   0.236,   0.000274,  2526,       0.304,
  "50%", 1000,     0.276,   0.228,   0.000818,  848,        0.947,
  "50%", 2000,     0.276,   0.228,   0.000818,  848,        0.947,
  "50%", 4000,     0.276,   0.228,   0.000813,  853,        0.944,
  "50%", 6000,     0.277,   0.228,   0.000583,  1189,       0.769,
  "50%", 8000,     0.281,   0.230,   0.000294,  2356,       0.328
)


upsilon_chosen <- 2000
max_upsilon <- max(sens_df$upsilon)

p_r2 <- ggplot(sens_df, aes(x = upsilon, y = r_squared, color = LC, group = LC)) +
  annotate("rect", xmin = -Inf, xmax = upsilon_chosen, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.4) +
  geom_hline(yintercept = 0.9, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = upsilon_chosen, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_point(data = subset(sens_df, upsilon == upsilon_chosen),
             size = 3.5, shape = 21, fill = "white", stroke = 1.1) +
  annotate("label", x = upsilon_chosen, y = min(sens_df$r_squared) + 0.02,
           label = paste0("", "\U03C5", " = ", upsilon_chosen),
           size = 3.2, hjust = -0.05, label.size = 0, fill = "white") +
  labs(x = expression(upsilon), y = expression(R^2), color = "LC") +
  theme_tufte() +
  theme(
    axis.title = element_text(size = 14),
    legend.position = "top"
  )


p_halflife <- ggplot(sens_df, aes(x = upsilon, y = half_life, color = LC, group = LC)) +
  annotate("rect", xmin = -Inf, xmax = upsilon_chosen, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.4) +
  geom_vline(xintercept = upsilon_chosen, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_point(data = subset(sens_df, upsilon == upsilon_chosen),
             size = 3.5, shape = 21, fill = "white", stroke = 1.1) +
  labs(x = expression(upsilon), y = expression(t[1/2])) +
  theme_tufte() +
  theme(
    axis.title = element_text(size = 14),
    legend.position = "none"  # avoid duplicating the legend
  )


combined_plot <- p_r2 + p_halflife +
  plot_annotation(
    title = NULL,
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

combined_plot



# ---------
# Violin plot
# ----------

summary_results <- tibble()
epochs <- c(3, 5, 10)
embeddings <- c(3, 4, 5, 6)
bits <- c(2, 4, 8, 16, 32)
p_alive_vals <- c(0.1, 0.3, 0.5, 0.8)

for (ep in epochs) {
  for (pa in p_alive_vals) {
    for (emb in embeddings) {
      for (b in bits) {
        
        file <- sprintf(
          "Data/results/%sK/thimming results/results_%sk_%s_emb%s_%sbits.csv",
          ep, ep, pa, emb, b
        )
        
        df <- read_csv(file, show_col_types = FALSE)
        
        if ("step" %in% names(df)) {
          df <- rename(df, "epoch" = "step")
        }
        
        step_min <- which.min(df$H)
        minH <- min(df$H)
        
        step_max <- max(df$epoch)
        step_last <- df[df$epoch == step_max, "H"] %>% pull(H)
        
        summary_results <- bind_rows(
          summary_results,
          tibble(
            epoch = factor(ep),
            p_alive = factor(pa),
            embedding = factor(emb),
            bits = factor(b),
            step_last = step_last,
            step_max = step_max,
            min_H = minH,
            step_min_H = step_min
          )
        )
      }
    }
  }
}

summary_results <- summary_results %>%
  mutate(
    ratio_min_H = step_min_H / step_max,
    dist_to_end = step_max - step_min_H
  )

plot_data <- summary_results %>% filter(!p_alive %in% c(0.1, 0.8))

n_per_group <- plot_data %>% count(epoch)

ggplot(plot_data, aes(x = epoch, y = step_min_H, fill = epoch)) +
  geom_violin(alpha = 0.5, color = NA, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA,
               color = "grey20", linewidth = 0.4) +
  geom_jitter(width = 0.05, alpha = 0.35, size = 1.2, color = "grey20") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.5,
               fill = "white", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Epoch (k)",
    y = "Generations until reaching minimum H"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

plot_data %>%
  filter(epoch %in% c("5", "10"), step_min_H > 4000) %>%
  select(epoch, p_alive, embedding, bits, step_min_H) %>%
  arrange(epoch, desc(step_min_H))

# ggsave("step_min_H_violin.pdf", width = 6, height = 4.5)