# Multivariate regression

library(ggplot2)
library(ggthemes)

data <- read.csv("edges_with_corrs_all.csv")

data <- data[data$weight<1000, ]

ggplot(data) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), col="red") +
  geom_point(aes(x=pearson_corr_flow, y=spearman_corr_flow)) +
  theme_clean() +
  theme(aspect.ratio = 1)

ggplot(data, aes(x=pearson_corr_speed, y=spearman_corr_speed)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), col="red") +
  geom_point() +
  theme_clean() +
  theme(aspect.ratio = 1)

ggplot(data, aes(x=pearson_corr_occupy, y=spearman_corr_occupy)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), col="red") +
  geom_point() +
  theme_clean() +
  theme(aspect.ratio = 1)

## Simple linear regression on all regressors

LinearRegression <- lm(data=data,
                       weight~pearson_corr_flow+spearman_corr_flow+
                         pearson_corr_speed+spearman_corr_speed+
                         pearson_corr_occupy+spearman_corr_occupy
                         )
summary(LinearRegression)

## Stepwise regression linear model

# 1. Fit a full model (all predictors)
full_model <- lm(data=data,
                 weight~pearson_corr_flow+spearman_corr_flow+
                   pearson_corr_speed+spearman_corr_speed+
                   pearson_corr_occupy+spearman_corr_occupy)

# 2. Fit an intercept-only model (no predictors)
null_model <- lm(data=data,
                 weight~1)

# 3. Perform stepwise selection (both directions)
stepwise_model <- step(null_model, 
                       scope = list(lower = null_model, upper = full_model), 
                       direction = "both")

# 4. View results
summary(stepwise_model)
plot(stepwise_model)

## Stepwise regression GLM gamma with reciprocal link function <<<REFAZER>>>

# 1. Fit a full model (all predictors)
full_gamma_model <- glm(weight ~ pearson_corr_flow + spearman_corr_flow +
                                pearson_corr_speed + spearman_corr_speed +
                                pearson_corr_occupy + spearman_corr_occupy,
                              data = data, 
                              family = Gamma(link = "inverse"))


# 2. Fit an intercept-only model (no predictors)
unitary_gamma_model <- glm(weight ~ 1,
                        data = data, 
                        family = Gamma(link = "inverse"))

# 3. Perform stepwise selection (both directions)
stepwise_gamma_model <- step(full_gamma_model, 
                       scope = list(lower = unitary_gamma_model, upper = full_gamma_model), 
                       direction = "both", trace = TRUE)

# 4. View results
summary(stepwise_gamma_model)
plot(stepwise_gamma_model)

