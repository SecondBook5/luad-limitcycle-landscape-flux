from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader

from hpcs_model.ml.datasets import DriftDataset
from hpcs_model.ml.models import DriftMLP


def train_drift_model(
    Z: np.ndarray,
    drift: np.ndarray,
    condition: np.ndarray,
    cfg: Dict[str, Any],
    output_dir: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    dataset = DriftDataset(Z, drift, condition)
    loader = DataLoader(
        dataset, batch_size=int(cfg["ml"]["batch_size"]), shuffle=True, drop_last=False
    )

    model = DriftMLP(
        input_dim=Z.shape[1],
        condition_dim=int(np.max(condition)) + 1,
        hidden_dim=int(cfg["ml"]["hidden_dim"]),
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=float(cfg["ml"]["lr"]))
    loss_fn = nn.MSELoss()

    model.train()
    for _ in range(int(cfg["ml"]["epochs"])):
        for z, cond, target in loader:
            pred = model(z, cond)
            loss = loss_fn(pred, target)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

    model_path = output_dir / "drift_model.pt"
    torch.save(model.state_dict(), model_path)
    return model_path
