# ============================================================
# Análise de distribuição de H* e C* (R)
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # para combinar gráficos lado a lado

# --- 1. Dados ---
df <- data.frame(
  rep = 1:30,
  H_star = c(
    0.2988388420242297, 0.2575637282735801, 0.29175957164136224, 0.29675734019977695,
    0.2811459660530557, 0.26348950637881535, 0.28317087331744845, 0.279708295251985,
    0.3075146892104318, 0.27477403477598683, 0.26026971179509756, 0.2522875039044994,
    0.2736103984276758, 0.25806036377534475, 0.24516365551762745, 0.2951249153457096,
    0.2691718073372512, 0.27448439205702607, 0.2677319104265286, 0.2772540182048604,
    0.29207996794573887, 0.2981924003391143, 0.30540077974661667, 0.23490949267584024,
    0.2883375110910137, 0.2638976748150433, 0.2554843210166479, 0.30390272693364945,
    0.2801033977250085, 0.31915127725232106
  ),
  C_star = c(
    0.25202037789547993, 0.22359451552352447, 0.24361100577090603, 0.24633445332481105,
    0.23868531738738863, 0.22552656324805026, 0.23877814898355762, 0.23855621765717008,
    0.258563689915975, 0.2344912943346862, 0.2250345545463719, 0.22332518438021237,
    0.23256191228980735, 0.22136303328364432, 0.2135841317012521, 0.24624659576431987,
    0.23042269835224913, 0.23402388415449707, 0.23173173830580487, 0.2354749391402692,
    0.24762419720861317, 0.2518656640886432, 0.2601889067670046, 0.2065821343062338,
    0.24536886466402752, 0.22745474586098705, 0.22095260460676794, 0.24708826703266057,
    0.23682391127580615, 0.26626993624501855
  )
)

# --- 2. Estatísticas descritivas ---
resumo <- function(x, nome) {
  cat("\n---", nome, "---\n")
  cat(sprintf("n        : %d\n", length(x)))
  cat(sprintf("média    : %.4f\n", mean(x)))
  cat(sprintf("mediana  : %.4f\n", median(x)))
  cat(sprintf("desvio-p : %.4f\n", sd(x)))
  cat(sprintf("erro-pad : %.4f\n", sd(x) / sqrt(length(x))))
  cat(sprintf("mín/máx  : %.4f / %.4f\n", min(x), max(x)))
  
  # Shapiro-Wilk (normalidade)
  sw <- shapiro.test(x)
  veredito <- if (sw$p.value > 0.05) "parece normal (p > 0.05)" else "foge da normalidade (p <= 0.05)"
  cat(sprintf("Shapiro-Wilk: W=%.4f, p=%.4f -> %s\n", sw$statistic, sw$p.value, veredito))
}

resumo(df$H_star, "H_star")
resumo(df$C_star, "C_star")

# Correlação entre H* e C*
ct <- cor.test(df$H_star, df$C_star)
cat(sprintf("\nCorrelação de Pearson entre H* e C*: r=%.4f, p=%.6f\n", ct$estimate, ct$p.value))

# --- 3. Histogramas + densidade ---
p_hist_H <- ggplot(df, aes(x = H_star)) +
  geom_histogram(aes(y = after_stat(density)), bins = 8,
                 fill = "#4C72B0", color = "white", alpha = 0.85) +
  geom_density(color = "#2C3E60", linewidth = 1) +
  geom_vline(xintercept = mean(df$H_star), linetype = "dashed", color = "black") +
  labs(title = "Histogram of H*", x = "H*", y = "Density") +
  theme_minimal()

p_hist_C <- ggplot(df, aes(x = C_star)) +
  geom_histogram(aes(y = after_stat(density)), bins = 8,
                 fill = "#DD8452", color = "white", alpha = 0.85) +
  geom_density(color = "#A85A2E", linewidth = 1) +
  geom_vline(xintercept = mean(df$C_star), linetype = "dashed", color = "black") +
  labs(title = "Histogram of C*", x = "C*", y = "Density") +
  theme_minimal()

p_hist_H + p_hist_C
ggsave("histogramas_H_C.png", width = 10, height = 4, dpi = 200)

# --- 4. Boxplots ---
df_long <- df %>% pivot_longer(cols = c(H_star, C_star), names_to = "variavel", values_to = "valor")

p_box <- ggplot(df_long, aes(x = variavel, y = valor, fill = variavel)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.08, color = "black", alpha = 0.5, size = 1.8) +
  scale_fill_manual(values = c("H_star" = "#4C72B0", "C_star" = "#DD8452")) +
  labs(title = "Boxplot de H* e C*", x = NULL, y = "Valor") +
  theme_minimal() +
  theme(legend.position = "none")

p_box
ggsave("boxplots_H_C.png", width = 6, height = 4.5, dpi = 200)

# --- 5. QQ-plots (normalidade) ---
p_qq_H <- ggplot(df, aes(sample = H_star)) +
  stat_qq(color = "#4C72B0") +