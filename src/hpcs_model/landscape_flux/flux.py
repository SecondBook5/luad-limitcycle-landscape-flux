from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from sklearn.neighbors import NearestNeighbors


def flux_residual(
    Z2: np.ndarray, U_grid: np.ndarray, grid_x: np.ndarray, grid_y: np.ndarray, drift: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    grad_x, grad_y = np.gradient(U_grid, grid_x, grid_y, edge_order=1)
    interp_gx = RegularGridInterpolator((grid_x, grid_y), grad_x)
    interp_gy = RegularGridInterpolator((grid_x, grid_y), grad_y)
    gradU = np.column_stack([interp_gx(Z2), interp_gy(Z2)])
    residual = drift + gradU
    mag = np.linalg.norm(residual, axis=1)

    grid = np.stack(np.meshgrid(grid_x, grid_y, indexing="ij"), axis=-1).reshape(-1, 2)
    nn = NearestNeighbors(n_neighbors=1).fit(Z2)
    _, idx = nn.kneighbors(grid)
    res_grid = residual[idx[:, 0]].reshape(len(grid_x), len(grid_y), 2)
    curl = np.gradient(res_grid[:, :, 1], grid_x, axis=0) - np.gradient(res_grid[:, :, 0], grid_y, axis=1)
    return residual, mag, curl
