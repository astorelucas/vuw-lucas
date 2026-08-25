## ============================================================
## Convergence to the Regime -- H x C plane (simple version)
##
##   H* = mean(H_t) over last v generations
##   C* = mean(C_t) over last v generations
##   d_t = sqrt((H_t - H*)^2 + (C_t - C*)^2)
##   log(d_t) = log(d_0) - r*t + eps_t   (OLS)
##   r = -slope
##   t_half = ln(2) / r
##
## Usage:
##   df <- read.csv("your_file.csv")   # needs columns: step, H, C
##   res <- convergence_to_regime(df, v = 200)
##   res$r
##   res$half_life
## ============================================================

convergence_to_regime <- function(df, v, t0 = 1) {
  
  df <- df[order(df$step), ]
  T <- nrow(df)
  
  if (v <= 1 || v >= T) stop("v must satisfy 1 < v < T")
  
  # asymptotic point: mean of last v generations
  tail_idx <- (T - v + 1):T
  H_star <- mean(df$H[tail_idx])
  C_star <- mean(df$C[tail_idx])
  
  # distance to asymptotic point
  d_t <- sqrt((df$H - H_star)^2 + (df$C - C_star)^2)
  
  # log-linear fit (drop d_t == 0, undefined in log space)
  keep <- df$step >= t0 & d_t > 0
  fit <- lm(log(d_t[keep]) ~ df$step[keep])
  
  slope <- coef(fit)[2]
  r <- -as.numeric(slope)
  half_life <- log(2) / r
  
  list(
    H_star = H_star,
    C_star = C_star,
    d_t = d_t,
    fit = fit,
    r = r,
    half_life = half_life,
    r_squared = summary(fit)$r.squared
  )
}

## Example:
p_alive <- c(0.1,0.3,0.5,0.8)
embdim <- c(5,6)
bits <- c(8,16,32)
i<-c(1,2,3)

 
source("convergence_simple.R")   # defines convergence_to_regime()

p_alive <- c(0.1, 0.3, 0.5, 0.8)
embdim  <- c(5, 6)
bits    <- c(8, 16, 32)
i       <- c(1, 2, 3)

v <- 6000

results <- data.frame()

for (ed in embdim) {
  for (b in bits) {
    for (p in p_alive) {
      for (rep in i) {
        
        file <- sprintf("Data/results/10K/new_results/emb%s/results_10k_%sbits_%s_%s.csv",
                        ed, b, p, rep)
        
        row <- data.frame(
          embdim = ed, bits = b, p_alive = p, rep = rep, file = file,
          H_star = NA, C_star = NA, r = NA, half_life = NA, r_squared = NA,
          status = "ok", error = NA
        )
        
        if (!file.exists(file)) {
          row$status <- "file_not_found"
          row$error <- "file does not exist"
          results <- rbind(results, row)
          next
        }
        
        res <- tryCatch({
          df <- read.csv(file)
          convergence_to_regime(df, v = v)
        }, error = function(e) e)
        
        if (inherits(res, "error")) {
          row$status <- "error"
          row$error <- conditionMessage(res)
        } else {
          row$H_star    <- res$H_star
          row$C_star    <- res$C_star
          row$r         <- res$r
          row$half_life <- res$half_life
          row$r_squared <- res$r_squared
        }
        
        results <- rbind(results, row)
        
        cat(sprintf("[%s] emb=%s bits=%s p=%s rep=%s\n", row$status, ed, b, p, rep))
      }
    }
  }
}

write.csv(results, "convergence_results_all.csv", row.names = FALSE)
cat("\nSaved:", nrow(results), "rows -> convergence_results_all.csv\n")
cat("Failed:", sum(results$status != "ok"), "\n")
