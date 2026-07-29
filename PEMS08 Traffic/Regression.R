library(ggplot2)
library(ggthemes)

## Data loading and cleansing

data <- read.csv("edges_with_corrs.csv")
names(data) <- c("s1", "s2", "impedance", "Pearson", "Spearman")

# Removed largest impedance
data <- data[data$impedance<1000, ]

# Removed saturated (=1) correlations
data <- data[data$Pearson<1, ]
data <- data[data$Spearman<1, ]

# EDA

ggplot(data) +
  geom_point(aes(
    x=Pearson, 
    y=impedance)
             ) +
  theme_clean()


## Regression analysis

ggplot(data, aes(x=Pearson, y=impedance)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = inverse.gaussian(link = "1/mu^2")), # or "poisson", etc.
              se = TRUE) + # Set to FALSE to remove the confidence interval
  theme_clean()
