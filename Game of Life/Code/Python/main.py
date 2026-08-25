import os
import time
import gzip
import pickle
import numpy as np
import torch
import pandas as pd
from tqdm.auto import tqdm
from functions import hilbert_idx, binary_to_real, compute_hxc, make_random_init, next_generation, hilbert_curve_order


if torch.cuda.is_available():
    device = torch.device("cuda")
elif torch.backends.mps.is_available():
    device = torch.device("mps")
else:
    device = torch.device("cpu")


print("Device:", device)

# Parameters
n_rows = 1024
n_cols = 1024
n_steps = 10000
p_alive = [0.1,0.3,0.5,0.8]
n_repeticoes = 10
epoch = 10
emb = 6
group_bits = [8, 16, 32]

coords = hilbert_curve_order(n_rows)
hi = hilbert_idx(n_cols, coords)
n_total = len(hi)
weights = {
    gb: (2 ** np.arange(gb - 1, -1, -1)).astype(np.uint32)
    for gb in group_bits
}

tempo_total_inicio = time.time()

# Main loop
for p in p_alive:
    for rep in range(1,n_repeticoes+1):

        print(f"\nRunning p={p}, rep={rep}")
        inicio = time.time()
        seed = int(p*1000 + rep)
        results = {gb: np.zeros((n_steps, 3)) for gb in group_bits}
        current = None
        for step in tqdm(range(n_steps), desc=f"p={p}, rep={rep}", unit="step", leave=True):

            # Initial state
            if step == 0:
                current = make_random_init(n_rows, n_cols, p, seed)
            else:
                current = next_generation(current)

            grid = (current.squeeze().cpu().numpy().astype(np.uint8))
            vals = grid.flatten()[hi]

            for gb in group_bits:
                groups = vals.reshape(-1, gb)

                ts = groups @ weights[gb]
                
                H,C = compute_hxc(ts, emb)
                results[gb][step] = [
                    step,
                    H,
                    C
                ]

        output_dir = f"Data/results/new_results/emb{emb}"

        os.makedirs(output_dir, exist_ok=True)

        for gb in group_bits:

            filename = (f"{output_dir}/results_{epoch}k_{gb}bits_{p}_{rep}.csv")

            df = pd.DataFrame(results[gb], columns=["step", "H", "C"])

            df.to_csv(filename, index=False)

            print("Saved:", filename)

        print("time:", round(time.time()-inicio, 2), "seconds")


print("\nTotal:", round((time.time()-tempo_total_inicio)/60, 2), "minutes")