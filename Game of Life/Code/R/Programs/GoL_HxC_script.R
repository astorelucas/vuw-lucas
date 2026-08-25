library(here)
## ------------------------------------------------------------------
## Simulação CA (Game of Life) + série temporal (curva de Hilbert) +
## HxC a cada passo, salvos em um .rds por (p, repetição)
## ------------------------------------------------------------------

n_rows <- 256
n_cols <- 256
n_steps <- 10000
p_alive <- c(0.5) # percentage_alive no estado inicial
n_repeticoes <- 3
epoch <- 10

# --- AJUSTE ESTE VALOR CONFORME O SEU PIPELINE ------------------------
emb <- 6          # dimensão de embedding usada em compute_hxc()
# -----------------------------------------------------------------

group_bits <- c(8) # bit-widths testados para agrupar a série

# A regra do GoL não muda a cada passo/repetição -> calcular uma única vez
parsed  <- game_of_life_rule()
born    <- parsed$born
survive <- parsed$survive

# A ordem da curva de Hilbert depende só do tamanho da grade -> calcular
# uma única vez fora de todos os loops (antes era recalculada, com "n"
# de step no lugar do tamanho da grade -- bug corrigido aqui).
order_coords <- hilbert_curve_order(n_rows)
n_total <- nrow(order_coords)

# `vals` (a sequência de bits ao longo da curva) é a MESMA para qualquer
# bit-width -- só o agrupamento em time.series muda. Por isso é
# calculada uma única vez por passo, e reaproveitada pelos 3 bit-widths,
# em vez de recalcular a simulação inteira 3x (uma por bit-width).
vals <- integer(n_total)

# estruturas que dependem do bit-width: uma entrada por valor em
# group_bits, pré-alocadas uma única vez
n_groups_list   <- setNames(vector("list", length(group_bits)), as.character(group_bits))
timeseries_list <- setNames(vector("list", length(group_bits)), as.character(group_bits))
for (gb in group_bits) {
  stopifnot(n_total %% gb == 0)
  n_groups_list[[as.character(gb)]] <- n_total %/% gb
  timeseries_list[[as.character(gb)]] <- numeric(n_total %/% gb)
}

# a cada quantos passos forcar o garbage collector a devolver memoria
# ao SO (o R normalmente NAO faz isso sozinho em loops longos e
# apertados, mesmo descartando objetos -- costuma ser a causa real de
# uso de memoria crescente nesse tipo de simulacao)
gc_every <- 200

tempo_total_inicio <- Sys.time()

for (p in p_alive) {
  for (i in 1:n_repeticoes) {
    
    cat("--> Running:", " rep=", i, " for p-alive ", p, "\n")
    tempo_rep_inicio <- Sys.time()
    
    seed <- as.integer(p * 1000 + i)
    
    # uma matriz de resultados (step, H, C) POR bit-width -- vals é
    # calculado uma vez por passo, so a serie/HxC repete por bit-width
    results_list <- setNames(
      lapply(group_bits, function(gb) {
        matrix(NA_real_, nrow = n_steps, ncol = 3,
               dimnames = list(NULL, c("step", "H", "C")))
      }),
      as.character(group_bits)
    )
    
    for (n in 1:n_steps) {
      
      if (n == 1) {
        # primeiro passo da repetição: estado inicial aleatório
        # (antes o teste era `if (i != 1)`, o que misturava a lógica de
        # repetição com a de passo -- corrigido para checar `n`)
        current <- make_random_init(n_rows, n_cols, p, seed = seed)
      } else {
        # evolui a matriz do passo anterior (t-1) um passo no GoL (t)
        current <- next_generation_generic(current, born, survive)
      }
      
      # integer (0/1) ocupa 1/2 da memoria de double para a mesma grade
      # -- se next_generation_generic ja devolver integer/logical, essa
      # linha nao custa nada (storage.mode so converte se precisar)
      storage.mode(current) <- "integer"
      
      # série de bits binária da matriz via curva de Hilbert
      # (antes usava `map`, variável que não existia -- corrigido para `current`)
      # calculada uma unica vez por passo, reaproveitada pelos 3 bit-widths
      vals[] <- current[order_coords]
      
      # para cada bit-width: agrupa `vals` e calcula HxC
      for (gb in group_bits) {
        gb_chr <- as.character(gb)
        n_groups <- n_groups_list[[gb_chr]]
        time.series <- timeseries_list[[gb_chr]]
        
        for (g in seq_len(n_groups)) {
          idx <- ((g - 1) * gb + 1):(g * gb)
          time.series[g] <- binary_to_real(vals[idx])
        }
        
        # HxC do passo n para este bit-width
        hc <- compute_hxc(time.series, emb = emb)
        
        # guarda o resultado desse passo
        # (ajuste hc$H / hc$C se compute_hxc() devolver outros nomes de campo)
        results_list[[gb_chr]][n, ] <- c(n, hc["H"], hc["C"])
      }
      
      # libera memoria periodicamente -- em loops de milhoes de
      # iteracoes o R tende a acumular memoria "solta" sem devolver
      # ao SO; forcar isso de tempos em tempos evita que o uso cresca
      # sem parar ao longo das 100 repeticoes
      if (n %% gc_every == 0) gc(verbose = FALSE)
    }
    
    # limpa a grade grande antes de comecar a proxima repeticao
    rm(current)
    gc(verbose = FALSE)
    
    tempo_rep_fim <- Sys.time()
    cat("   tempo repeticao:",
        round(difftime(tempo_rep_fim, tempo_rep_inicio, units = "secs"), 2),
        "s\n")
    
    # salva um arquivo por bit-width ao final de cada repetição
    for (gb in group_bits) {
      out_file <- sprintf(
        "Data/results/%sK/new_results/emb%s/results_%sk_%sbits_%s_%s.rds",
        epoch, emb, epoch, gb, p, i
      )
      saveRDS(results_list[[as.character(gb)]], file = out_file, compress = "xz")
      cat("Saved full history to", out_file, "\n")
    }
  }
}

tempo_total_fim <- Sys.time()
cat("Tempo total:",
    round(difftime(tempo_total_fim, tempo_total_inicio, units = "mins"), 2),
    "min\n")
