library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
install.packages("sp")
library(sp)

#palive80dim5bits8
## Graph HxC

epoch <- 10
p_alive <- 0.8
embdim <- 6
bits <- 32
#file <- sprintf("Data/results/%sK/emb5_0.8_16bits/results_%sk_%s_emb%s_%sbits_1.csv", epochs, epochs,p_alive, emb, b, i)

file <- sprintf("Data/results/%dK/results_%dk_%s_emb%d_%dbits.csv", epoch, epoch, p_alive, embdim, bits)
K <- 2000

hc_df <- read_csv(file, show_col_types = FALSE)
hc_df$epoch <- 1:nrow(hc_df)
hc_df$epoch <- factor(hc_df$epoch)  # discrete labels/colors, one per epoch
hc_df$H <- as.numeric(hc_df$H)
hc_df$C <- as.numeric(hc_df$C)
hc_df <- hc_df[!is.na(hc_df$H) & !is.na(hc_df$C), ] # remover linhas em que H ou C são NA
hc_df$t <- as.numeric(as.character(hc_df$epoch))
label_epochs <- unique(c(min(hc_df$t), 
                         hc_df$t[hc_df$t %% K == 0], 
                         max(hc_df$t)))
hc_df_labels <- subset(hc_df, t %in% label_epochs)
#hc_df_labels <- hc_df_labels %>%
 # mutate(
  #  nudge_x = case_when(
   #   t == 2000 ~ -0.05,
    #  TRUE ~ 0
    #),
    #nudge_y = case_when(
    #  t == 2000 ~ -0.04,   # empurra o rótulo pra baixo, longe da curva/pontos
    #  TRUE ~ 0
    #),
    #nudge_x = case_when(
    #  t == 4000 ~ -0.05,
    #  TRUE ~ 0
    #),
    #nudge_y = case_when(
    #  t == 4000 ~ -0.04,   # empurra o rótulo pra baixo, longe da curva/pontos
    #  TRUE ~ 0
    #)
  #)

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
                  box.padding = 0.6,
                  point.padding = 0.4,
                  #nudge_x = hc_df_labels$nudge_x,
                  #nudge_y = hc_df_labels$nudge_y,
                  segment.color = "grey60",
                  segment.size = 0.8,
                  seed = 42) +
  coord_cartesian(clip = "off") +  # permite rótulos saírem levemente da área do plot
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(
    title = sprintf("p(alive) = %d%%  |  dim = %d  |  bits = %d", p_alive*100, embdim, bits),
    x = "Permutation Entropy (H)",
    y = "Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey92", linewidth = 0.2))

# Graph

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


plot_dist_to_limit(hc_df)

## Graph Final(H) vs Final(C)

embeddings <- c(5)
bits <- 16
epochs <- 3
p_alive <- 0.8
tail_epochs <- 500
n_repeticoes <- 101

summary_results <- data.frame()

for (i in 2:n_repeticoes){
  for (emb in embeddings) {
    for (b in bits) {
      
      #file <- sprintf("Data/results/%sK/results_%sk_%s_emb%s_%sbits.csv", epochs, epochs,p_alive, emb, b)
      file <- sprintf("Data/results/%sK/results_%sk_%s_emb%s_%sbits_%s.csv", epochs, epochs,p_alive, emb, b, i)
      
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
  }
}

ggplot(summary_results,
       aes(mean_H, mean_C,
           color = factor(embedding))) +
  geom_path(aes(group = factor(embedding)),
            linewidth = 0.6, alpha = 0.7,linetype = "dashed",) +
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

summary_results <- tibble()

for (emb in embeddings) {
  for (b in bits) {
    
    file <- sprintf(
      "Data/results/%sK/results_%sk_0.1_emb%s_%sbits.csv",
      epochs, epochs, emb, b
    )
    
    df <- read_csv(file, show_col_types = FALSE)
    
    final_df <- df %>% slice_tail(n = tail_epochs)
    
    H <- final_df$H
    C <- final_df$C
    
    meanH <- mean(H)
    meanC <- mean(C)
    
    ## Radius
    radius <- mean(sqrt((H - meanH)^2 +
                          (C - meanC)^2))
    
    ## Trajectory length
    traj_length <- sum(
      sqrt(diff(H)^2 + diff(C)^2)
    )
    
    ## Mean velocity in HxC
    mean_velocity <- mean(
      sqrt(diff(H)^2 + diff(C)^2)
    )
    
    ## Convex hull area
    hull <- chull(H, C)
    hull_pts <- cbind(H[hull], C[hull])
    
    x <- H[hull]
    y <- C[hull]
    
    # close polygon
    x <- c(x, x[1])
    y <- c(y, y[1])
    
    hull_area <- 0.5 * abs(
      sum(x[-1] * y[-length(y)] -
            x[-length(x)] * y[-1])
    )
    
    ## Covariance
    covHC <- cov(H, C)
    
    summary_results <- bind_rows(
      summary_results,
      tibble(
        embedding = emb,
        bits = b,
        
        mean_H = meanH,
        mean_C = meanC,
        
        sd_H = sd(H),
        sd_C = sd(C),
        
        cov_HC = covHC,
        
        radius = radius,
        
        trajectory_length = traj_length,
        
        mean_velocity = mean_velocity,
        
        hull_area = hull_area
      )
    )
  }
}

ggplot(summary_results,
       aes(mean_H,
           mean_C,
           color=factor(embedding))) +
  
  geom_path(aes(group=embedding),
            linewidth=.7) +
  
  geom_errorbar(
    aes(ymin=mean_C-sd_C,
        ymax=mean_C+sd_C),
    width=0
  ) +
  
  geom_errorbarh(
    aes(xmin=mean_H-sd_H,
        xmax=mean_H+sd_H),
    height=0
  ) +
  
  geom_point(aes(shape=factor(bits)),
             size=4) +
  
  theme_bw()

ggplot(summary_results,
       aes(bits,
           radius,
           fill=factor(embedding)))+
  geom_col(position="dodge")+
  theme_bw()

ggplot(summary_results,
       aes(bits,
           trajectory_length,
           color=factor(embedding),
           group=embedding))+
  geom_line()+
  geom_point(size=3)+
  theme_bw()

ggplot(summary_results,
       aes(bits,
           hull_area,
           color=factor(embedding),
           group=embedding))+
  geom_line()+
  geom_point(size=3)+
  theme_bw()

plot_df <- tibble()

for (emb in embeddings) {
  for (b in bits) {
    
    file <- sprintf(
      "Data/results/%sK/results_%sk_0.4_emb%s_%sbits.csv",
      epochs, epochs, emb, b
    )
    
    df <- read_csv(file, show_col_types = FALSE)
    
    final_df <- df %>%
      slice_tail(n = tail_epochs) %>%
      mutate(
        embedding = factor(emb),
        bits = factor(b)
      )
    
    plot_df <- bind_rows(plot_df, final_df)
  }
}

ggplot(plot_df,
       aes(H, C,
           colour = embedding)) +
  
  stat_ellipse(
    aes(group = interaction(embedding, bits)),
    level = 0.95,
    linewidth = 1
  ) +
  
  stat_summary(
    aes(group = interaction(embedding,bits),
        shape = bits),
    fun = mean,
    geom = "point",
    size = 4
  ) +
  
  theme_bw()

stat_ellipse(
  aes(fill=embedding,
      group=interaction(embedding,bits)),
  geom="polygon",
  alpha=.20,
  colour=NA
)

## Gerar 3 imagens

summary_results <- tibble()
epochs <- c(3,5,10)
tail_epochs <- 1400

embeddings <- c(3,4,5,6)
bits <- c(2, 4, 8, 16, 32)
p_alive <- 0.8
tail_epochs <- 1400


for (ep in epochs) {
  for (emb in embeddings) {
    for (b in bits) {
      
      if (ep == 10 && p_alive == 0.8){
        file <- sprintf(
          "Data/results/%sK/results_%sk_%s_emb%s_%sbits_1.csv",
          ep, ep, p_alive, emb, b
        )
      }else{
        file <- sprintf(
          "Data/results/%sK/results_%sk_%s_emb%s_%sbits.csv",
          ep, ep, p_alive, emb, b
        )
      }

      df <- read_csv(file, show_col_types = FALSE)
      
      if (ep == 3){
        tail_epochs <- 400
      } 
      if (ep == 5){
        tail_epochs <- 700
      }
      final_df <- df %>%
        slice_tail(n = tail_epochs)
      
      n <- nrow(final_df)
      
      summary_results <- bind_rows(
        summary_results,
        tibble(
          epoch = factor(ep),
          embedding = factor(emb),
          bits = factor(b),
          
          mean_H = mean(final_df$H),
          ci_H = qt(0.975, n-1) * sd(final_df$H)/sqrt(n),
          
          mean_C = mean(final_df$C),
          ci_C = qt(0.975, n-1) * sd(final_df$C)/sqrt(n)
        )
      )
    }
  }
}
facet_wrap(~epoch, nrow = 1)


ggplot(summary_results,
       aes(mean_H,
           mean_C,
           color = embedding)) +
  
  geom_path(aes(group = embedding),
            linewidth = 0.6,
            alpha = 0.7,
            linetype = "dashed") +
  
  geom_errorbar(aes(ymin = mean_C-ci_C,
                    ymax = mean_C+ci_C),
                width = 0) +
  
  geom_errorbarh(aes(xmin = mean_H-ci_H,
                     xmax = mean_H+ci_H),
                 height = 0) +
  
  geom_point(aes(shape = bits),
             size = 4) +
  
  facet_wrap(~epoch, nrow = 1) +
  
  theme_bw()


## 100 execucoes

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
