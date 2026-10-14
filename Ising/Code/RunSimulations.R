library(here)
library(Rcpp)
sourceCpp(here("Code", "ising.cpp"))

set.seed(2343443)   # any integer; record it with your results

L <- 512
beta <- log(1 + sqrt(2)) / 2
initial <- matrix(sample(c(-1L, 1L), L * L, replace = TRUE), L, L)

out_dir <- here("Data", "Simulations")

files <- c(MH = file.path(out_dir, "spins_MH.bin"),
           GB = file.path(out_dir, "spins_GB.bin"),
           SW = file.path(out_dir, "spins_SW.bin"))
cons <- lapply(files, file, open = "wb")

# +1 -> TRUE, -1 -> FALSE, packed 8 spins per byte
save_gen <- function(con, spins) {
  writeBin(packBits(as.vector(spins) == 1L, type = "raw"), con)
  flush(con)   # so nothing is lost if you interrupt the run
}

# generation 0 = initial condition
for (m in names(files)) save_gen(cons[[m]], initial)

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

tryCatch({
  for (t in 1:1000) {
    spins_MH <- metropolis_step(spins_MH, beta)
    spins_GB <- gibbs_step(spins_GB, beta)
    spins_SW <- swendsen_wang_step(spins_SW, beta)
    
    save_gen(cons$MH, spins_MH)
    save_gen(cons$GB, spins_GB)
    save_gen(cons$SW, spins_SW)
    
    
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
}, finally = {
  lapply(cons, close)
})
