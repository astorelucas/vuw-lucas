
library(ggplot2)
library(reshape2)  # or use tidyr::pivot_longer


mat <- ca_history10k_0.5[[10000]]


df <- melt(mat)  # columns: Var1 (row), Var2 (col), value
names(df) <- c("row", "col", "state")

ggplot(df, aes(x = col, y = row, fill = factor(state))) +
  geom_raster() +
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +
  scale_y_reverse() +  # so row 1 is at top
  coord_fixed() +      # square cells
  theme_void() +
  theme(legend.position = "none")
