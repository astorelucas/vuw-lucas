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
p_alive <- 0.8
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



# -------------------------------------------------------
# Graph 100 execucoes to reach H min
# -------------------------------------------------------

emb <- 5
b <- 16
epochs <- 10
p_alive <- 0.8
tail_epochs <- 1400
n_repeticoes <- 101

summary_results <- data.frame()

for (i in 2:n_repeticoes){
  
  file <- sprintf("Data/results/%sK/emb5_0.8_16bits/results_%sk_%s_emb%s_%sbits_%s.csv", epochs, epochs,p_alive, emb, b, i)
  
  df <- read_csv(file, show_col_types = FALSE)
  
  final_df <- df %>% slice_tail(n = tail_epochs)
  
  n <- nrow(final_df)
  
  summary_results <- bind_rows(
    summary_results,
    tibble(
      embedding = emb,
      bits = b,
      
      mean_H = mean(final_df$H),
      ci_H = qt(0.975, n - 1) * sd(final_df$H) / sqrt(n),
      
      mean_C = mean(final_df$C),
      ci_C = qt(0.975, n - 1) * sd(final_df$C) / sqrt(n)
    )
  )
}

ggplot(summary_results,
       aes(mean_H, mean_C,
           color = factor(embedding))) +
  geom_errorbar(aes(ymin = mean_C - ci_C,
                    ymax = mean_C + ci_C),
                width = 0) +
  geom_errorbarh(aes(xmin = mean_H - ci_H,
                     xmax = mean_H + ci_H),
                 height = 0) +
  geom_point(aes(shape = factor(bits)), size = 4) +
  theme_bw() +
  labs(
    x = "Final H",
    y = "Final C",
    color = "Embedding",
    shape = "Hilbert bits"
  )
summary_results
