from __future__ import annotations

import torch
from torch import nn


class DriftMLP(nn.Module):
    def __init__(self, input_dim: int, condition_dim: int, hidden_dim: int = 128):
        super().__init__()
        self.embed = nn.Embedding(condition_dim, hidden_dim)
        self.net = nn.Sequential(
            nn.Linear(input_dim + hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim),
        )

    def forward(self, z: torch.Tensor, condition: torch.Tensor) -> torch.Tensor:
        emb = self.embed(condition)
        x = torch.cat([z, emb], dim=1)
        return self.net(x)
