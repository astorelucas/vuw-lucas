# Multivariate regression

library(ggplot2)
library(ggthemes)

data <- read.csv("Dataset/sensor_pairs_model.csv")

# valor absoluto do Pearson e Spearman
data$pearson_abs  <- abs(data$pearson_r)
data$spearman_abs <- abs(data$spearman_r)

# EXPLORAÇÃO: Pearson vs Spearman
ggplot(data, aes(x = pearson_abs, y = spearman_abs)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), col = "red") +
  geom_point() +
  labs(title = "|Pearson| vs |Spearman|",
       x = "|Pearson r|", y = "|Spearman r|") +
  theme_clean() +
  theme(aspect.ratio = 1)

ggplot(data, aes(x = dist_km, y = pearson_abs)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red", se = TRUE) +
  labs(title = "|Pearson r| vs distância geográfica",
       x = "Distância (km)", y = "|Pearson r|") +
  theme_clean()

ggplot(data, aes(x = dist_km, y = pearson_abs)) +
  geom_point() +
  geom_smooth(method = "loess", col = "blue", se = TRUE) +
  geom_smooth(method = "lm",    col = "red",  se = FALSE, linetype = "dashed") +
  labs(title = "|Pearson r| vs distância — linear vs loess",
       x = "Distância (km)", y = "|Pearson r|") +
  theme_clean()

## Simple linear regression on all regressors

LinearRegression <- lm(
  pearson_abs ~ dist_km + spearman_abs,
  data = data
)
summary(LinearRegression)

LinearRegression_p <- lm(
  spearman_abs ~ dist_km + pearson_abs,
  data = data
)
summary(LinearRegression_p)

# STEPWISE LM — AIC, ambas as direções

full_model <- lm(
  pearson_abs ~ dist_km + spearman_abs,
  data = data
)
null_model <- lm(pearson_abs ~ 1, data = data)

stepwise_model <- step(
  null_model,
  scope     = list(lower = null_model, upper = full_model),
  direction = "both"
)
summary(stepwise_model)
plot(stepwise_model)

# --------------------------------------------------
# STEPWISE GLM BETA — |pearson_r| em (0, 1)
# mais apropriado que Gamma para proporções
# --------------------------------------------------
install.packages("betareg")

library(betareg)

# betareg exige y estritamente em (0,1) — afastar de 0 e 1 se necessário
epsilon <- 1e-4
data$pearson_abs_adj <- pmin(pmax(data$pearson_abs, epsilon), 1 - epsilon)

full_beta_model <- betareg(
  pearson_abs_adj ~ dist_km + spearman_abs,
  data = data
)
summary(full_beta_model)

# stepwise manual para betareg (não suporta step() nativo)
# comparar modelos por AIC
m_full  <- betareg(pearson_abs_adj ~ dist_km + spearman_abs, data = data)
m_dist  <- betareg(pearson_abs_adj ~ dist_km,                data = data)
m_spear <- betareg(pearson_abs_adj ~ spearman_abs,           data = data)
m_null  <- betareg(pearson_abs_adj ~ 1,                      data = data)

cat("\n=== AIC comparação Beta regression ===\n")
cat("Full (dist + spearman) :", AIC(m_full),  "\n")
cat("Só dist_km             :", AIC(m_dist),  "\n")
cat("Só spearman_abs        :", AIC(m_spear), "\n")
cat("Null                   :", AIC(m_null),  "\n")

best_beta <- list(m_full, m_dist, m_spear, m_null)[
  which.min(c(AIC(m_full), AIC(m_dist), AIC(m_spear), AIC(m_null)))
][[1]]

cat("\nModelo Beta selecionado:\n")
summary(best_beta)
plot(best_beta)


# --------------------------------------------------
# COMPARAÇÃO LM vs Beta
# --------------------------------------------------
cat("\n=== AIC: LM stepwise vs Beta regression ===\n")
cat("LM stepwise :", AIC(stepwise_model), "\n")
cat("Beta melhor :", AIC(best_beta),      "\n")

