library(plotly)

nbins <- 100

hr <- range(emb6_bits16$H)
cr <- range(emb6_bits16$C)

# log-spaced bin edges, since the data itself is log-spread
bx <- 10^seq(log10(hr[1]), log10(hr[2]), length.out = nbins + 1)
by <- 10^seq(log10(cr[1]), log10(cr[2]), length.out = nbins + 1)

counts <- table(
  cut(emb6_bits16$H, breaks = bx, include.lowest = TRUE),
  cut(emb6_bits16$C, breaks = by, include.lowest = TRUE)
)
counts <- matrix(as.numeric(counts), nrow = nbins, ncol = nbins)

# geometric mean as bin center — the right "middle" for a log-spaced bin
xmid <- sqrt(bx[-1] * bx[-(nbins + 1)])
ymid <- sqrt(by[-1] * by[-(nbins + 1)])

plot_ly(
  x = ~xmid, y = ~ymid, z = ~t(counts),
  type = "surface",
  colorscale = list(c(0, "white"), c(0.5, "steelblue"), c(1, "darkblue"))
) |>
  layout(
    title = "2D histogram of H and C (log-log axes)",
    scene = list(
      xaxis = list(title = "H", type = "log"),
      yaxis = list(title = "C", type = "log"),
      zaxis = list(title = "count")
    )
  )

