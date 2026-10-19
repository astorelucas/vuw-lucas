source(here("Code", "Functions.R"))
library(readr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(here)

method_names <- c(MH = "Metropolis", GB = "Gibbs", SW = "Swendsen-Wang")

emb <- 6
bits <- 16
beta_vals <- (0.6610)

for (b in beta_vals) {
  p <- plot_hxc_methods(b, emb = emb, bits = bits, root = here("Results/seed1000-99-01-15"))
  ggsave(here("Figures",
              sprintf("hxc_methods_beta%.2f_emb%d_%dbits.png", b, emb, bits)),
         p, width = 15, height = 5, dpi = 400, bg = "white")
}

# ---- alejandro 

library(plotly)

nbins <- 30

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

# ---- alejandro 

library(data.table)
library(plotly)
library(htmlwidgets)
library(here)

# ---- parametros --------------------------------------------------------------
base <- "Results/seed1000_90-10"
beta <- "beta1.0000"
arq  <- "emb6_bits16.csv"

arquivos <- c("Metropolis"    = file.path(base, "MH",   beta, arq),
              "Gibbs"         = file.path(base, "GB", beta, arq),
              "Swendsen-Wang" = file.path(base, "SW",   beta, arq))
tempo_log <- TRUE   # TRUE: eixo z = log10(geracao + 1)  (abre o comeco da dinamica)
# FALSE: eixo z linear
marcos    <- c(0, 5, 50, 100, 500, 1000)  # rotulos do eixo z (como no seu grafico 2D)

# ---- leitura -------------------------------------------------------------------
read_traj <- function(f) {
  d <- fread(f)[, 1:3]
  setnames(d, c("generation", "H", "C"))
  d[order(generation)]
}
dados <- lapply(arquivos, read_traj)

zf   <- function(g) if (tempo_log) log10(g + 1) else g
gmin <- min(sapply(dados, function(d) min(d$generation)))
gmax <- max(sapply(dados, function(d) max(d$generation)))
zmin <- zf(gmin)

zaxis <- list(title = if (tempo_log) "Geração (escala log)" else "Geração")
if (tempo_log) {
  zaxis$tickvals <- zf(marcos)
  zaxis$ticktext <- as.character(marcos)
}
cena <- function() list(
  xaxis = list(title = "H", range = c(0, 1)),
  yaxis = list(title = "C", range = c(0, 0.5)),
  zaxis = zaxis,
  aspectmode = "manual", aspectratio = list(x = 1, y = 0.7, z = 0.9),
  camera = list(eye = list(x = 1.6, y = -1.6, z = 0.9))
)

# ---- trajetoria de um metodo em uma cena -------------------------------------
add_traj <- function(p, d, scene, show_bar = FALSE) {
  p %>%
    # trajetoria colorida pelo tempo
    add_trace(type = "scatter3d", mode = "lines+markers", scene = scene,
              x = d$H, y = d$C, z = zf(d$generation),
              line   = list(color = d$generation, colorscale = "Viridis",
                            width = 4, cmin = gmin, cmax = gmax),
              marker = list(size = 2, color = d$generation, colorscale = "Viridis",
                            cmin = gmin, cmax = gmax, showscale = show_bar,
                            colorbar = list(title = "Geração")),
              text = paste0("geração ", d$generation),
              hoverinfo = "text+x+y", showlegend = FALSE) %>%
    # sombra no piso = o plano H x C original
    add_trace(type = "scatter3d", mode = "lines", scene = scene,
              x = d$H, y = d$C, z = rep(zmin, nrow(d)),
              line = list(color = "gray", width = 2),
              hoverinfo = "skip", showlegend = FALSE)
}

# ---- figura 1: tres cenas lado a lado ------------------------------------------
nomes <- names(dados)
p1 <- plot_ly()
for (i in seq_along(dados)) {
  p1 <- add_traj(p1, dados[[i]], scene = c("scene", "scene2", "scene3")[i],
                 show_bar = (i == length(dados)))
}
dom <- list(c(0, 0.33), c(0.335, 0.665), c(0.67, 1))
sc  <- lapply(dom, function(x) c(cena(), list(domain = list(x = x, y = c(0, 1)))))
p1 <- p1 %>% layout(
  scene = sc[[1]], scene2 = sc[[2]], scene3 = sc[[3]],
  annotations = lapply(seq_along(nomes), function(i)
    list(text = paste0("<b>", nomes[i], "</b>"),
         x = mean(dom[[i]]), y = 1, xref = "paper", yref = "paper",
         showarrow = FALSE, font = list(size = 16)))
)
p1
saveWidget(p1, "hxc_3d_lado_a_lado.html", selfcontained = TRUE)

# ---- figura 2: metodos sobrepostos na mesma cena (cor = metodo) --------------------
cores <- c("#E4572E", "#17BEBB", "#2E294E")
p2 <- plot_ly()
for (i in seq_along(dados)) {
  d <- dados[[i]]
  p2 <- p2 %>% add_trace(type = "scatter3d", mode = "lines+markers",
                         x = d$H, y = d$C, z = zf(d$generation), name = nomes[i],
                         line = list(color = cores[i], width = 4),
                         marker = list(size = 2, color = cores[i]),
                         text = paste0("geração ", d$generation),
                         hoverinfo = "text+x+y")
}
p2 <- p2 %>% layout(scene = cena())
p2
saveWidget(p2, "hxc_3d_sobreposto.html", selfcontained = TRUE)









library(tidyverse)
library(boot)

files <- c(A = "sim_A.csv", B = "sim_B.csv", C = "sim_C.csv")  # ajuste os nomes
upsilon <- 0.8   # janela final = últimos 20%
eps     <- 0.05

# ---------- bootstrap em blocos para (H*, C*) ----------
block_boot_mean <- function(x, R = 2000) {
  l <- max(2, round(length(x)^(1/3) * 2))      # tamanho de bloco heurístico
  b <- tsboot(x, mean, R = R, l = l, sim = "fixed")
  c(est = mean(x), lo = quantile(b$t, .025), hi = quantile(b$t, .975))
}

analyze <- function(path, name) {
  d <- read_csv(path, show_col_types = FALSE) |> arrange(generation)
  T_ <- nrow(d)
  tail_idx <- ceiling(upsilon * T_):T_
  Hs <- mean(d$H[tail_idx]); Cs <- mean(d$C[tail_idx])
  
  d <- d |> mutate(dist = sqrt((H - Hs)^2 + (C - Cs)^2))
  
  # --- fase de decaimento: de t=1 até onde dist cai a ~ 2x o piso ---
  floor_d <- median(d$dist[tail_idx])
  t_end   <- which(d$dist < 2 * floor_d)[1]
  if (is.na(t_end) || t_end < 10) t_end <- floor(T_ * 0.1)
  fit_lm  <- lm(log(dist) ~ generation, data = d[1:t_end, ] |> filter(dist > 0))
  r_lm    <- -coef(fit_lm)[["generation"]]
  
  # --- ajuste não linear com piso ---
  fit_nls <- try(nls(dist ~ A * exp(-r * generation) + dinf, data = d,
                     start = list(A = d$dist[1], r = max(r_lm, 1e-4), dinf = floor_d),
                     control = nls.control(maxiter = 200, warnOnly = TRUE)), silent = TRUE)
  r_nls <- if (inherits(fit_nls, "try-error")) NA else coef(fit_nls)[["r"]]
  dinf  <- if (inherits(fit_nls, "try-error")) NA else coef(fit_nls)[["dinf"]]
  
  # --- tempo de assentamento: última vez que dist > eps * dist[1] ---
  thr <- eps * d$dist[1]
  t_settle <- (max(which(d$dist > thr), 0)) + 1
  
  # --- bootstrap em blocos de H* e C* ---
  bH <- block_boot_mean(d$H[tail_idx]); bC <- block_boot_mean(d$C[tail_idx])
  
  # --- variância e autocorrelação na janela final ---
  acf_H <- acf(d$H[tail_idx], plot = FALSE, lag.max = 100)$acf[-1]
  tau_H <- which(acf_H < exp(-1))[1]    # lag onde ACF cai abaixo de 1/e
  
  # --- tendência residual no platô ---
  slope_H <- coef(lm(H ~ generation, data = d[tail_idx, ]))[2]
  slope_C <- coef(lm(C ~ generation, data = d[tail_idx, ]))[2]
  
  # --- comprimento do caminho no plano ---
  path_len <- sum(sqrt(diff(d$H)^2 + diff(d$C)^2))
  
  list(
    data = d |> mutate(sim = name),
    summary = tibble(sim = name, Hstar = Hs, Cstar = Cs,
                     H_lo = bH["lo.2.5%"], H_hi = bH["hi.97.5%"],
                     C_lo = bC["lo.2.5%"], C_hi = bC["hi.97.5%"],
                     r_loglin = r_lm, thalf_loglin = log(2) / r_lm,
                     r_nls = r_nls, thalf_nls = log(2) / r_nls, d_inf = dinf,
                     t_settle = t_settle, sd_H = sd(d$H[tail_idx]),
                     sd_C = sd(d$C[tail_idx]), tau_H = tau_H,
                     slope_H = slope_H, slope_C = slope_C, path_len = path_len)
  )
}

res <- imap(files, ~ analyze(.x, .y))
summ <- map_dfr(res, "summary")
traj <- map_dfr(res, "data")
print(summ, width = Inf)

# ---------- distância entre regimes (par a par) vs incerteza ----------
pairs <- combn(names(files), 2, simplify = FALSE)
pair_tbl <- map_dfr(pairs, function(p) {
  a <- summ |> filter(sim == p[1]); b <- summ |> filter(sim == p[2])
  tibble(par = paste(p, collapse = " vs "),
         dist_regime = sqrt((a$Hstar - b$Hstar)^2 + (a$Cstar - b$Cstar)^2),
         # os ICs se sobrepõem nas duas dimensões?
         overlap_H = !(a$H_hi < b$H_lo | b$H_hi < a$H_lo),
         overlap_C = !(a$C_hi < b$C_lo | b$C_hi < a$C_lo))
})
print(pair_tbl)

# ---------- distância de Fréchet / DTW entre trajetórias ----------
# install.packages(c("dtw","princurve")) se necessário
library(dtw)
traj_mat <- map(res, ~ as.matrix(.x$data[, c("H", "C")]))
dtw_tbl <- map_dfr(pairs, function(p) {
  tibble(par = paste(p, collapse = " vs "),
         dtw = dtw(traj_mat[[p[1]]], traj_mat[[p[2]]], distance.only = TRUE)$normalizedDistance)
})
print(dtw_tbl)

# ---------- gráficos ----------
ggplot(traj, aes(generation, dist, colour = sim)) +
  geom_line(alpha = .7) + scale_y_log10() +
  labs(y = "d_t (log)", title = "Convergência ao regime") + theme_minimal()

ggplot(traj, aes(H, C, colour = sim)) +
  geom_path(alpha = .6) +
  geom_point(data = summ, aes(Hstar, Cstar), size = 3, shape = 4, stroke = 1.5) +
  theme_minimal() + labs(title = "Trajetórias no plano H×C")
