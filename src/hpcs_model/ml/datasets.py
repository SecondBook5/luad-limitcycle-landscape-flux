from __future__ import annotations

from typing import Tuple

import numpy as np
import torch
from torch.utils.data import Dataset


class DriftDataset(Dataset):
    def __init__(self, Z: np.ndarray, drift: np.ndarray, condition: np.ndarray):
        mask = ~np.isnan(drift).any(axis=1)
        self.Z = torch.tensor(Z[mask], dtype=torch.float32)
        self.drift = torch.tensor(drift[mask], dtype=torch.float32)
        self.condition = torch.tensor(condition[mask], dtype=torch.long)

    def __len__(self) -> int:
        return self.Z.shape[0]

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        return self.Z[idx], self.condition[idx], self.drift[idx]
