L <- 1024
beta <- log(1+sqrt(2))/2

# Same initial condition
initial <- matrix(
  sample(c(-1L, 1L), L * L, replace = TRUE),
  L, L
)

spins_MH <- initial
spins_GB <- initial
spins_SW <- initial

plot_spin <- function(spins, main = "") {
  image(
    t(spins[nrow(spins):1, ]),
    col = c("black", "white"),
    axes = FALSE,
    asp = 1,
    main = main
  )
}

# Open one graphics device
quartz(width = 12, height = 4)

for (t in 1:10000) {
  
  spins_MH <- metropolis_step(spins_MH, beta)
  spins_GB <- gibbs_step(spins_GB, beta)
  spins_SW <- swendsen_wang_step(spins_SW, beta)
  
  par(
    mfrow = c(1, 3),
    mar = c(1, 1, 3, 1)
  )
  
  plot_spin(
    spins_MH,
    sprintf("Metropolis\nSweep %d", t)
  )
  
  plot_spin(
    spins_GB,
    sprintf("Gibbs\nSweep %d", t)
  )
  
  plot_spin(
    spins_SW,
    sprintf("Swendsen-Wang\nSweep %d", t)
  )
  
  Sys.sleep(1)
}
