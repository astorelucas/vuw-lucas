library(tidyverse)


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