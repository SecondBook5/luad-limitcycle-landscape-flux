from __future__ import annotations

from typing import Any, Dict, Tuple

import numpy as np
from sklearn.neighbors import NearestNeighbors


def estimate_density(Z: np.ndarray, cfg: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray]:
    k = int(cfg["landscape_flux"]["density_k"])
    nbrs = NearestNeighbors(n_neighbors=k).fit(Z)
    distances, _ = nbrs.kneighbors(Z)
    radius = np.mean(distances, axis=1) + 1e-8
    density = 1.0 / (radius**Z.shape[1])
    return density, radius
