import numpy as np
import torch
import torch.nn.functional as F
import ordpy
from hilbertcurve.hilbertcurve import HilbertCurve

if torch.cuda.is_available():
    device = torch.device("cuda")
elif torch.backends.mps.is_available():
    device = torch.device("mps")
else:
    device = torch.device("cpu")

# Game of Life
kernel = torch.tensor(
    [[1,1,1],
     [1,0,1],
     [1,1,1]],
    dtype=torch.float32,
    device=device
).view(1,1,3,3)


def make_random_init(rows, cols, p_alive, seed):
    torch.manual_seed(seed)
    return (
        torch.rand(
            (1,1,rows,cols),
            device=device
        ) < p_alive
    )

@torch.no_grad()
def next_generation(mat):
    x = mat.float()
    x = F.pad(
        x,
        (1,1,1,1),
        mode="circular"
    )
    neighbors = F.conv2d(x, kernel)

    alive = mat
    born = (
        (~alive) &
        (neighbors == 3)
    )
    survive = (
        alive &
        ((neighbors == 2) |
         (neighbors == 3))
    )
    return born | survive

# Hilbert curve
def hilbert_curve_order(size):

    order = int(np.log2(size))
    hilbert = HilbertCurve(order, 2)

    coords = np.empty((size * size, 2), dtype=np.int64)

    for i in range(size * size):
        coords[i] = hilbert.point_from_distance(i)

    return coords

# transforma coordenadas em índices numpy
def hilbert_idx(n_cols, coords):
    return (
        coords[:,0] * n_cols +
        coords[:,1]
    ).astype(np.int64)

def binary_to_real(bits):
    value = 0
    for b in bits:
        value = value*2 + int(b)
    return value

def compute_hxc(series, emb):
    H, C = ordpy.complexity_entropy(
        series,
        dx=emb,
        dy=1,
        taux=1,
        tauy=1
    )

    return H, C