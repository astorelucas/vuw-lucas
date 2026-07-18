# ============================================================
# Generalized 2D Cellular Automaton Simulator in R
# ============================================================
# Supports any outer-totalistic B/S ("Birth/Survival") rule on a
# toroidal (periodic) Moore neighborhood grid, plus Brian's Brain,
# which is a 3-state automaton and needs its own transition logic.
#
# Rules implemented:
#   Conway's Game of Life   B3/S23
# ============================================================

# Grid dimensions (kept modest since several of these rules grow
# explosively - e.g. Seeds - and 3000 steps at 256x256 for seven
# rules is heavy; adjust n_rows/n_cols/n_steps to taste)
n_rows <- 256
n_cols <- 256
n_steps <- 2999
perturbation_at <- NULL
p_alive <- 0.1 # percentage_alive in initial_state

# ------------------------------------------------------------
# 1. Neighbor counting (toroidal / wrap-around boundary)
# ------------------------------------------------------------
count_neighbors <- function(mat) {
  nrow_m <- nrow(mat)
  ncol_m <- ncol(mat)
  neighbor_count <- matrix(0, nrow = nrow_m, ncol = ncol_m)
  
  offsets <- list(
    c(-1, -1), c(-1, 0), c(-1, 1),
    c(0, -1),            c(0, 1),
    c(1, -1),  c(1, 0),  c(1, 1)
  )
  
  for (off in offsets) {
    shifted_rows <- ((1:nrow_m - 1 + off[1]) %% nrow_m) + 1
    shifted_cols <- ((1:ncol_m - 1 + off[2]) %% ncol_m) + 1
    neighbor_count <- neighbor_count + mat[shifted_rows, shifted_cols]
  }
  
  neighbor_count
}

# ------------------------------------------------------------
# 2. Parse a "B.../S..." rule string into birth/survive vectors
# ------------------------------------------------------------
# e.g. "B3/S23"   -> list(born = 3,       survive = c(2, 3))
#      "B2/S"     -> list(born = 2,       survive = integer(0))
#      "B3678/S34678" -> list(born = c(3,6,7,8), survive = c(3,4,6,7,8))
parse_rule <- function(rule_string) {
  parts <- strsplit(rule_string, "/")[[1]]
  b_part <- sub("^B", "", parts[1])
  s_part <- sub("^S", "", parts[2])
  
  born    <- if (nchar(b_part) == 0) integer(0) else as.integer(strsplit(b_part, "")[[1]])
  survive <- if (nchar(s_part) == 0) integer(0) else as.integer(strsplit(s_part, "")[[1]])
  
  list(born = born, survive = survive)
}

# ------------------------------------------------------------
# 3. Generic outer-totalistic transition (covers GoL, HighLife,
#    Day & Night, Seeds, Maze, Replicator)
# ------------------------------------------------------------
next_generation_generic <- function(mat, born, survive) {
  neighbors <- count_neighbors(mat)
  
  survive_mask <- (mat == 1) & (neighbors %in% survive)
  born_mask    <- (mat == 0) & (neighbors %in% born)
  
  new_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  new_mat[survive_mask | born_mask] <- 1
  new_mat
}

# ------------------------------------------------------------
# 4. Brian's Brain (3-state: 0 = off/ready, 1 = firing/on,
#    2 = dying/refractory)
#    - A ready cell (0) fires (1) if exactly 2 neighbors are firing.
#      Only "firing" (state 1) cells count as neighbors.
#    - A firing cell (1) always moves to dying (2).
#    - A dying cell (2) always moves to off (0).
# ------------------------------------------------------------
next_generation_brian <- function(mat) {
  firing <- matrix(as.integer(mat == 1), nrow = nrow(mat), ncol = ncol(mat))
  firing_neighbors <- count_neighbors(firing)
  
  new_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  
  becomes_firing <- (mat == 0) & (firing_neighbors == 2)
  becomes_dying  <- (mat == 1)                     # every firing cell decays
  becomes_off    <- (mat == 2)                     # every dying cell resets
  
  new_mat[becomes_firing] <- 1
  new_mat[becomes_dying]  <- 2
  new_mat[becomes_off]    <- 0
  
  new_mat
}

# ------------------------------------------------------------
# 5. Unified simulation driver
# ------------------------------------------------------------
# rule_name: one of "gol", "highlife", "daynight", "seeds",
#            "maze", "brian", "replicator"
# perturb_step: time step at which one random cell is flipped
#               (for Brian's Brain, "flip" toggles between 0 and 1)
RULES <- list(
  gol        = "B3/S23"
)

simulate_ca <- function(initial_mat, n_steps, rule_name = "gol", perturb_step = NULL) {
  history <- vector("list", n_steps + 1)
  history[[1]] <- initial_mat
  current <- initial_mat
  
  is_brian <- identical(rule_name, "brian")
  
  if (!is_brian) {
    parsed <- parse_rule(RULES[[rule_name]])
    born <- parsed$born
    survive <- parsed$survive
  }
  
  for (t in 1:n_steps) {
    current <- if (is_brian) {
      next_generation_brian(current)
    } else {
      next_generation_generic(current, born, survive)
    }
    
    if (!is.null(perturb_step) && t == perturb_step) {
      rand_row <- sample(1:nrow(current), 1)
      rand_col <- sample(1:ncol(current), 1)
      old_value <- current[rand_row, rand_col]
      
      if (is_brian) {
        # toggle only between off(0) and firing(1); leave dying(2) alone
        new_value <- if (old_value == 0) 1 else if (old_value == 1) 0 else old_value
      } else {
        new_value <- 1 - old_value
      }
      
      current[rand_row, rand_col] <- new_value
      cat(sprintf("Perturbation applied at t = %d: cell (%d, %d) flipped from %d to %d\n",
                  t, rand_row, rand_col, old_value, new_value))
    }
    
    history[[t + 1]] <- current
  }
  
  history
}

# ------------------------------------------------------------
# 6. Grid printing (supports 2-state and 3-state grids)
# ------------------------------------------------------------
print_grid <- function(mat) {
  chars <- matrix(".", nrow = nrow(mat), ncol = ncol(mat))
  chars[mat == 1] <- "#"
  chars[mat == 2] <- "~"   # dying state, only used by Brian's Brain
  
  for (i in 1:nrow(chars)) {
    cat(paste(chars[i, ], collapse = " "), "\n")
  }
  cat("\n")
}

# ============================================================
# Example usage: run all seven automata and save each history
# ============================================================

set.seed(12345678, kind = "Mersenne-Twister")

make_random_init <- function(n_rows, n_cols, p_alive = 0.5) {
  matrix(
    sample(c(0, 1), n_rows * n_cols, replace = TRUE, prob = c(1 - p_alive, p_alive)),
    nrow = n_rows, ncol = n_cols
  )
}

# Suggested initial densities per rule (dense random grids work for
# GoL/HighLife/Day&Night/Maze; explosive rules like Seeds and
# Replicator are usually seeded sparse to stay readable)
initial_density <- list(
  gol        = p_alive
)

rule_names <- c("gol")

all_results_CA_simu <- list()

for (rn in rule_names) {
  cat("========================================\n")
  cat("Running rule:", rn, "\n")
  cat("========================================\n")
  
  init_mat <- make_random_init(n_rows, n_cols, p_alive = initial_density[[rn]])
  
  results_CAsimulation <- simulate_ca(
    initial_mat  = init_mat,
    n_steps      = n_steps,
    rule_name    = rn,
    perturb_step = perturbation_at
  )
  
  all_results_CA_simu[[rn]] <- results_CAsimulation
  
  cat("Simulation complete:", length(results_CAsimulation), "time steps generated.\n")
  
  out_file <- sprintf("ca_history3k_%s.rds", rn)
  saveRDS(results_CAsimulation, file = out_file)
  cat("Saved full history to", out_file, "\n\n")
}

cat("All rules simulated. Access any run via all_results[[\"gol\"]], all_results[[\"brian\"]], etc.\n")
