library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)

# --- monta o data frame a partir da sua tabela ---
df <- tribble(
  ~bits, ~L,    ~Tamanho, ~H_mean, ~C_mean, ~H_sd,   ~C_sd,   ~n_rep,
  8,     128,   2048,     0.2783,  0.2368,  0.0199,  0.0137,  3,
  8,     256,   8192,     0.2804,  0.2341,  0.0122,  0.0076,  3,
  8,     512,   32768,    0.2979,  0.2397,  0.0068,  0.0032,  3,
  8,     1024,  131072,   0.2989,  0.2386,  0.0041,  0.0023,  3,
  8,     2048,  524288,   0.2872,  0.2454,  0.0045,  0.0045,  3,
  
  16,    128,   1024,     0.3956,  0.3193,  0.0245,  0.0142,  3,
  16,    256,   4096,     0.4047,  0.3117,  0.01514, 0.00756, 3,
  16,    512,   16384,    0.4074,  0.3064,  0.0095,  0.0040,  3,
  16,    1024,  65536,    0.4110,  0.3045,  0.0056,  0.0017,  3,
  16,    2048,  262144,   0.4152,  0.3055,  0.0047,  0.0015,  3,
  
  32,    128,   512,      0.5408,  0.4012,  0.0321,  0.0121,  3,
  32,    256,   2048,     0.5655,  0.3764,  0.0193,  0.0048,  3,
  32,    512,   8192,     0.5913,  0.3549,  0.0163,  0.0044,  3,
  32,    1024,  32768,    0.5815,  0.3503,  0.0043,  0.0015,  3,
  32,    2048,  131072,   0.5981,  0.3438,  0.0097,  0.0027,  3
) %>%
  mutate(
    bits   = factor(bits, levels = c(8, 16, 32)),
    H_se   = H_sd / sqrt(n_rep),   # erro padrão da média, não o dp bruto
    C_se   = C_sd / sqrt(n_rep)
  )
# --- formato longo, para facetar H* e C* juntos ---
df_long <- df %>%
  pivot_longer(
    cols = c(H_mean, C_mean),
    names_to = "metric",
    values_to = "mean"
  ) %>%
  mutate(
    se = ifelse(metric == "H_mean", H_se, C_se),
    metric = recode(metric, H_mean = "H*", C_mean = "C*"),
    metric = factor(metric, levels = c("H*", "C*"))  # <- garante H* à esquerda
  )

ggplot(df_long, aes(x = L, y = mean, color = bits, group = bits)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = bits),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.05) +
  scale_x_continuous(trans = "log2", breaks = c(256, 512, 1024, 2048)) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "L",
    y = "Average",
    color = "Bits",
    fill = "Bits",
    #title = "Convergência de H* e C* com o tamanho da grade (L)",
    #subtitle = "Barras = erro padrão da média (dp / \u221an)"
  ) +
  theme_tufte(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 13)
  )
