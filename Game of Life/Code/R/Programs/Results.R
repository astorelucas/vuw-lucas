library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
install.packages("sp")
library(sp)
library(StatOrdPattHxC)
library(ggrepel)
library(ggthemes)

## ----------------------------------------
## Graph HxC
## ----------------------------------------

epoch <- 10
p_alive <- 0.1
embdim <- 6
bits <- 16
i<-1

file <- sprintf("Data/results/%sK/new_results/emb6/results_%sk_%s_emb%s_%sbits_%s.csv", epoch, epoch, p_alive, embdim, bits,i)
#file <- sprintf("Data/results/%sK/results_%sk_%sbits_%s_%s.csv",epoch, epoch, bits, p_alive, i)
#file <- sprintf("Data/results/%dK/thimming results/results_%dk_%s_emb%d_%dbits.csv", epoch, epoch, p_alive, embdim, bits)
#file <- sprintf("Data/results/%dK/thimming results/results_3k_0.3_emb3_2bits.csv",epoch )

K <- 5000 # Show label for each K steps

hc_df <- read_csv(file, show_col_types = FALSE)
#hc_df <- readRDS(file)
hc_df <- as.data.frame(hc_df)
hc_df$epoch <- 1:nrow(hc_df)
hc_df$epoch <- factor(hc_df$epoch)  # discrete labels/colors, one per epoch
hc_df$H <- as.numeric(hc_df$H)
hc_df$C <- as.numeric(hc_df$C)
#hc_df <- hc_df[!is.na(hc_df$H) & !is.na(hc_df$C), ] # remover linhas em que H ou C são NA
hc_df$t <- as.numeric(as.character(hc_df$epoch))
#label_epochs <- c(1,50,100,200, 500, 1000, 3000, 8000, 10000)  # whatever epochs you want
label_epochs <- c(1,5000, 10000)  # whatever epochs you want
#label_epochs <- unique(c(min(hc_df$t), 
#                         hc_df$t[hc_df$t %% K == 0], 
#                         max(hc_df$t)))
hc_df_labels <- subset(hc_df, t %in% label_epochs)
ggplot(hc_df, aes(x = H, y = C, color = epoch)) +
  geom_line(data = subset(LinfLsup, Dimension == embdim & Side == "Lower"),
            aes(x = H, y = C), color = "grey40", linewidth = 0.6,
            inherit.aes = FALSE) +
  geom_line(data = subset(LinfLsup, Dimension == embdim & Side == "Upper"),
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
    nudge_y = 0.05,    
  ) +
  coord_cartesian(clip = "off") +  # permite rótulos saírem levemente da área do plot
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(
    title = sprintf("%dk generations | p(alive) = %d%%  |  dim = %d  |  bits = %d",epoch, p_alive*100, embdim, bits),
    x = "Permutation Entropy (H)",
    y = "Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey92", linewidth = 0.2))


## ----------------------------------------
## Convergence
## ----------------------------------------

## ============================================================
## Convergence to the Regime -- H x C plane (simple version)
##
##   H* = mean(H_t) over last v generations
##   C* = mean(C_t) over last v generations
##   d_t = sqrt((H_t - H*)^2 + (C_t - C*)^2)
##   log(d_t) = log(d_0) - r*t + eps_t   (OLS)
##   r = -slope
##   t_half = ln(2) / r
##
## Usage:
##   df <- read.csv("your_file.csv")   # needs columns: step, H, C
##   res <- convergence_to_regime(df, v = 200)
##   res$r
##   res$half_life
## ============================================================

convergence_to_regime <- function(df, v, t0 = 1) {
  
  df <- df[order(df$step), ]
  T <- nrow(df)
  
  if (v <= 1 || v >= T) stop("v must satisfy 1 < v < T")
  
  # asymptotic point: mean of last v generations
  tail_idx <- (T - v + 1):T
  H_star <- mean(df$H[tail_idx])
  C_star <- mean(df$C[tail_idx])
  
  # distance to asymptotic point
  d_t <- sqrt((df$H - H_star)^2 + (df$C - C_star)^2)
  
  # log-linear fit (drop d_t == 0, undefined in log space)
  keep <- df$step >= t0 & d_t > 0
  fit <- lm(log(d_t[keep]) ~ df$step[keep])
  
  slope <- coef(fit)[2]
  r <- -as.numeric(slope)
  half_life <- log(2) / r
  
  list(
    H_star = H_star,
    C_star = C_star,
    d_t = d_t,
    fit = fit,
    r = r,
    half_life = half_life,
    r_squared = summary(fit)$r.squared
  )
}

p_alive <- c(0.1,0.3,0.5,0.8)
embdim <- c(5,6)
bits <- c(8,16,32)
i<-c(1,2,3)
v <- 3000

results <- data.frame()

for (ed in embdim) {
  for (b in bits) {
    for (p in p_alive) {
      for (rep in i) {
        
        if(ed ==5){
          file <- sprintf("Data/results/10K/new_results/emb%s/results_10k_%sbits_%s_%s.csv",
                          ed, b, p, rep)
        }else{
          file <- sprintf("Data/results/10K/new_results/emb%s/results_10k_%s_emb%s_%sbits_%s.csv",
                          ed, p, ed, b, rep)
        }

        
        row <- data.frame(
          embdim = ed, bits = b, p_alive = p, rep = rep, file = file,
          H_star = NA, C_star = NA, r = NA, half_life = NA, r_squared = NA,
          status = "ok", error = NA
        )
        
        if (!file.exists(file)) {
          row$status <- "file_not_found"
          row$error <- "file does not exist"
          results <- rbind(results, row)
          next
        }
        
        res <- tryCatch({
          df <- read.csv(file)
          convergence_to_regime(df, v = v)
        }, error = function(e) e)
        
        if (inherits(res, "error")) {
          row$status <- "error"
          row$error <- conditionMessage(res)
        } else {
          row$H_star    <- res$H_star
          row$C_star    <- res$C_star
          row$r         <- res$r
          row$half_life <- res$half_life
          row$r_squared <- res$r_squared
        }
        
        results <- rbind(results, row)
        
        cat(sprintf("[%s] emb=%s bits=%s p=%s rep=%s\n", row$status, ed, b, p, rep))
      }
    }
  }
}

write.csv(results, "convergence_results_all.csv", row.names = FALSE)
cat("\nSaved:", nrow(results), "rows -> convergence_results_all.csv\n")
cat("Failed:", sum(results$status != "ok"), "\n")



# ---------------------------------------------------------------
# Velocity of convergence in the HxC plane
# ---------------------------------------------------------------
N_tail <- 3000
epoch <- 10
p_alive <- 0.1
embdim <- 6
bits <- 16
i<-1

file <- sprintf("Data/results/%sK/new_results/emb6/results_%sk_%s_emb%s_%sbits_%s.csv", epoch, epoch, p_alive, embdim, bits,i)

hc_df <- read_csv(file, show_col_types = FALSE)


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
    x = "Generation",
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
