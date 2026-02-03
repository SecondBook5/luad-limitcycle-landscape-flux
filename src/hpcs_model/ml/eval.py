from __future__ import annotations

from typing import Any, Dict

import numpy as np
import torch

from hpcs_model.ml.models import DriftMLP


def evaluate_drift_model(
    Z: np.ndarray,
    drift: np.ndarray,
    condition: np.ndarray,
    cfg: Dict[str, Any],
    model_state: dict,
) -> float:
    model = DriftMLP(
        input_dim=Z.shape[1],
        condition_dim=int(np.max(condition)) + 1,
        hidden_dim=int(cfg["ml"]["hidden_dim"]),
    )
    model.load_state_dict(model_state)
    model.eval()
    mask = ~np.isnan(drift).any(axis=1)
    with torch.no_grad():
        z = torch.tensor(Z[mask], dtype=torch.float32)
        cond = torch.tensor(condition[mask], dtype=torch.long)
        target = torch.tensor(drift[mask], dtype=torch.float32)
        pred = model(z, cond)
        mse = torch.mean((pred - target) ** 2).item()
    return mse
