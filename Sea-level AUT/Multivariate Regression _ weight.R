library(ggplot2)
library(ggthemes)
library(betareg)

data <- read.csv("sensor_pairs_model.csv")

# remover pares sem edge no grafo KNN
data <- data[!is.na(data$weight), ]

data$pearson_abs  <- abs(data$pearson_r)
data$spearman_abs <- abs(data$spearman_r)

# EXPLORAÇÃO
ggplot(data, aes(x = dist_km, y = weight)) +
  geom_point() +
  geom_smooth(method = "lm",    col = "red",  se = TRUE) +
  geom_smooth(method = "loess", col = "blue", se = FALSE, linetype = "dashed") +
  labs(title = "weight vs distância", x = "dist_km", y = "weight") +
  theme_clean()

ggplot(data, aes(x = pearson_abs, y = weight)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red", se = TRUE) +
  labs(title = "weight vs |Pearson|", x = "|pearson_r|", y = "weight") +
  theme_clean()

ggplot(data, aes(x = w_dist, y = w_corr)) +
  geom_point() +
  labs(title = "componentes do weight: w_dist vs w_corr") +
  theme_clean()

# REGRESSÃO LINEAR — todos os preditores

lm_full <- lm(
  weight ~ dist_km + pearson_abs + spearman_abs + w_dist + w_corr,
  data = data
)
summary(lm_full)

# --------------------------------------------------
# STEPWISE LM
# --------------------------------------------------
null_model <- lm(weight ~ 1, data = data)
full_model <- lm(
  weight ~ dist_km + pearson_abs + spearman_abs + w_dist + w_corr,
  data = data
)

stepwise_lm <- step(
  null_model,
  scope     = list(lower = null_model, upper = full_model),
  direction = "both"
)
summary(stepwise_lm)
plot(stepwise_lm)


# BETA REGRESSION — weight em (0, 1)

epsilon <- 1e-4
data$weight_adj <- pmin(pmax(data$weight, epsilon), 1 - epsilon)

full_beta <- betareg(
  weight_adj ~ dist_km + pearson_abs + spearman_abs + w_dist + w_corr,
  data = data
)

null_beta <- betareg(weight_adj ~ 1, data = data)

# stepwise manual por AIC
candidates <- c("dist_km", "pearson_abs", "spearman_abs", "w_dist", "w_corr")

stepwise_beta_aic <- function(data, candidates) {
  selected <- c()
  best_aic <- AIC(betareg(weight_adj ~ 1, data = data))
  improved <- TRUE
  
  while (improved) {
    improved <- FALSE
    
    # forward
    remaining <- setdiff(candidates, selected)
    for (feat in remaining) {
      trial   <- c(selected, feat)
      formula <- as.formula(paste("weight_adj ~", paste(trial, collapse = " + ")))
      aic     <- AIC(betareg(formula, data = data))
      if (aic < best_aic) {
        best_aic <- aic
        best_feat <- feat
        improved <- TRUE
      }
    }
    if (improved) {
      selected <- c(selected, best_feat)
      cat("  +", best_feat, " AIC=", round(best_aic, 3), "\n")
    }
    
    # backward
    if (length(selected) > 1) {
      for (feat in selected) {
        trial   <- setdiff(selected, feat)
        formula <- as.formula(paste("weight_adj ~", paste(trial, collapse = " + ")))
        aic     <- AIC(betareg(formula, data = data))
        if (aic < best_aic) {
          best_aic <- aic
          selected <- trial
          improved <- TRUE
          cat("  -", feat, " AIC=", round(best_aic, 3), "\n")
        }
      }
    }
  }
  return(selected)
}

cat("=== Stepwise Beta (AIC) ===\n")
selected_vars <- stepwise_beta_aic(data, candidates)
cat("Variáveis selecionadas:", selected_vars, "\n")

formula_beta <- as.formula(
  paste("weight_adj ~", paste(selected_vars, collapse = " + "))
)
stepwise_beta <- betareg(formula_beta, data = data)
summary(stepwise_beta)
plot(stepwise_beta)


# COMPARAÇÃO FINAL

cat("\n=== AIC: LM stepwise vs Beta stepwise ===\n")
cat("LM stepwise  :", AIC(stepwise_lm),   "\n")
cat("Beta stepwise:", AIC(stepwise_beta),  "\n")

cat("\n=== R² ===\n")
cat("LM   R²       :", summary(stepwise_lm)$r.squared,        "\n")
cat("Beta pseudo-R²:", summary(stepwise_beta)$pseudo.r.squared, "\n")