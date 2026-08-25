import os
import time
import numpy as np
import torch
import pandas as pd
from tqdm.auto import tqdm
from functions import hilbert_idx, compute_hxc, make_random_init, next_generation, hilbert_curve_order


if torch.cuda.is_available():
    device = torch.device("cuda")
elif torch.backends.mps.is_available():
    device = torch.device("mps")
else:
    device = torch.device("cpu")

print("Device:", device)

# ----------------------------------------------------------------------
# Parameters — set these for the ONE table line you're running
# ----------------------------------------------------------------------
n_rows = 2048
n_cols = 2048
n_steps = 5000
p_alive = 0.5
n_repeticoes = 30
emb = 6
bits = 32
tail_frac = 0.6  # fraction of trajectory used to compute (H*, C*)

epoch_label = n_steps // 1000  # e.g. 10000 -> "10k", avoids manual mismatch with n_steps

assert (n_rows * n_cols) % bits == 0, "grid size not divisible by bits"

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------
coords = hilbert_curve_order(n_rows)
hi = hilbert_idx(n_cols, coords)
n_total = len(hi)

weights = (2 ** np.arange(bits - 1, -1, -1)).astype(np.uint32)

output_dir = f"Data/results/ts-analysis/emb{emb}"
os.makedirs(output_dir, exist_ok=True)

tail_start = int((1 - tail_frac) * n_steps)  # step index where the tail window begins

tempo_total_inicio = time.time()

# ----------------------------------------------------------------------
# Main loop: one repetition = one full CA trajectory.
# Only H,C values from the final tail_frac of the trajectory are kept
# (used to compute this repetition's converged (H*,C*)); nothing per-rep
# is written to disk.
# ----------------------------------------------------------------------
final_H = []
final_C = []

for rep in range(1, n_repeticoes + 1):
    print(f"\nRunning rep={rep}")
    inicio = time.time()
    seed = int(rep * 1000)

    current = None

    for step in tqdm(range(n_steps), desc=f"rep={rep}", unit="step", leave=True):

        if step == 0:
            current = make_random_init(n_rows, n_cols, p_alive, seed)
        else:
            current = next_generation(current)

        grid = current.squeeze().cpu().numpy().astype(np.uint8)

        vals = grid.flatten()[hi]
        groups = vals.reshape(-1, bits)
        ts = groups @ weights

        H, C = compute_hxc(ts, emb)

        if step >= tail_start:
            final_H.append(H)
            final_C.append(C)

    print("time:", round(time.time() - inicio, 2), "seconds")

# ----------------------------------------------------------------------
# Aggregate: (H*, C*) as mean over the final tail_frac of each trajectory,
# pooled across all repetitions
# ----------------------------------------------------------------------
mean_H = np.mean(final_H)
mean_C = np.mean(final_C)
std_H = np.std(final_H)
std_C = np.std(final_C)

print(f"\n(H*, C*) = ({mean_H:.4f}, {mean_C:.4f})  +- ({std_H:.4f}, {std_C:.4f})")

summary_path = f"{output_dir}/summary_{n_rows}x{n_cols}_{bits}bits.csv"
pd.DataFrame([{
    "bits": bits,
    "L": n_rows,
    "seq_len": n_total // bits,
    "n_reps": n_repeticoes,
    "H_star": mean_H,
    "C_star": mean_C,
    "H_star_std": std_H,
    "C_star_std": std_C,
}]).to_csv(summary_path, index=False)
print("Saved summary:", summary_path)

print("\nTotal:", round((time.time() - tempo_total_inicio) / 60, 2), "minutes")