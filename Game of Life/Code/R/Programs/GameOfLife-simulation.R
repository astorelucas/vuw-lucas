# ============================================================
# Conway's Game of Life - Cellular Automaton Simulation in R
# ============================================================
# Simulates the Game of Life and stores each time step as a
# binary matrix (0 = dead, 1 = alive) in a list.
# ============================================================

# --- Function: count live neighbors for every cell (with periodic/toroidal boundary) ---
count_neighbors <- function(mat) {
  nrow_m <- nrow(mat)
  ncol_m <- ncol(mat)
  neighbor_count <- matrix(0, nrow = nrow_m, ncol = ncol_m)
  
  # 8 relative neighbor offsets (Moore neighborhood)
  offsets <- list(
    c(-1, -1), c(-1, 0), c(-1, 1),
    c(0, -1),            c(0, 1),
    c(1, -1),  c(1, 0),  c(1, 1)
  )
  
  for (off in offsets) {
    # shift matrix using modulo arithmetic -> toroidal (wrap-around) boundary
    shifted_rows <- ((1:nrow_m - 1 + off[1]) %% nrow_m) + 1
    shifted_cols <- ((1:ncol_m - 1 + off[2]) %% ncol_m) + 1
    neighbor_count <- neighbor_count + mat[shifted_rows, shifted_cols]
  }
  
  neighbor_count
}

# --- Function: apply Game of Life rules to get the next generation ---
next_generation <- function(mat) {
  neighbors <- count_neighbors(mat)
  
  # Rules:
  # - A live cell with 2 or 3 live neighbors survives
  # - A dead cell with exactly 3 live neighbors becomes alive
  # - All other cells die or stay dead
  survive <- (mat == 1) & (neighbors %in% c(2, 3))
  born    <- (mat == 0) & (neighbors == 3)
  
  new_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  new_mat[survive | born] <- 1
  new_mat
}

# --- Function: run the simulation for n_steps and return all matrices ---
simulate_game_of_life <- function(initial_mat, n_steps) {
  history <- vector("list", n_steps + 1)
  history[[1]] <- initial_mat
  
  current <- initial_mat
  for (t in 1:n_steps) {
    current <- next_generation(current)
    history[[t + 1]] <- current
  }
  
  history
}

# --- Function: print a matrix using # for alive and . for dead (optional visualization) ---
print_grid <- function(mat) {
  chars <- ifelse(mat == 1, "#", ".")
  for (i in 1:nrow(chars)) {
    cat(paste(chars[i, ], collapse = " "), "\n")
  }
  cat("\n")
}

# ============================================================
# Example usage
# ============================================================

set.seed(123456789, kind="Mersenne-Twister")

# Grid dimensions
n_rows <- 256
n_cols <- 256

# Random initial state (about 30% alive)
initial_state <- matrix(
  sample(c(0, 1), n_rows * n_cols, replace = TRUE, prob = c(0.7, 0.3)),
  nrow = n_rows, ncol = n_cols
)

# Run simulation for 20 time steps
n_steps <- 3000
results <- simulate_game_of_life(initial_state, n_steps)

# 'results' is a list of (n_steps + 1) matrices, one per time point (0..n_steps)
# Each matrix is n_rows x n_cols of 0s and 1s

cat("Simulation complete:", length(results), "time steps generated.\n")
cat("Each time step is a", n_rows, "x", n_cols, "matrix of 0/1 values.\n\n")

# Print first and last grid as a quick visual sanity check
cat("Initial state (t = 0):\n")
print_grid(results[[1]])

cat("Final state (t =", n_steps, "):\n")
print_grid(results[[n_steps + 1]])

# Example: access the matrix at time t = 5
# results[[6]]  # remember: index 1 = t = 0, so t = 5 is index 6

# Optional: save all matrices to a single .rds file for later use
saveRDS(results, file = "game_of_life_history.rds")
cat("Full history saved to game_of_life_history.rds\n")

# Optional: save each matrix as a separate CSV (uncomment to use)
# for (t in seq_along(results)) {
#   write.csv(results[[t]], file = sprintf("gol_t%03d.csv", t - 1), row.names = FALSE)
# }

