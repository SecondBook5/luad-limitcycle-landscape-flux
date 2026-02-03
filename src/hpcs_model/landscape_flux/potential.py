from __future__ import annotations

from typing import Any, Dict

import numpy as np


def compute_potential(density: np.ndarray, cfg: Dict[str, Any]) -> np.ndarray:
    eps = float(cfg["landscape_flux"]["eps"])
    return -np.log(density + eps)
