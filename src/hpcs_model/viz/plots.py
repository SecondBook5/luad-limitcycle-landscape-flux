from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import numpy as np


def plot_embedding(Z2: np.ndarray, labels: np.ndarray, output: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(Z2[:, 0], Z2[:, 1], c=labels, s=6, cmap="tab20")
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(sc, ax=ax, pad=0.02)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def plot_potential(grid_x: np.ndarray, grid_y: np.ndarray, U_grid: np.ndarray, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    c = ax.contourf(grid_x, grid_y, U_grid.T, levels=30, cmap="viridis")
    fig.colorbar(c, ax=ax)
    ax.set_title("Potential landscape")
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def plot_drift(Z2: np.ndarray, drift: np.ndarray, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.quiver(Z2[:, 0], Z2[:, 1], drift[:, 0], drift[:, 1], angles="xy", scale_units="xy", scale=1)
    ax.set_title("Drift field")
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def plot_flux_magnitude(Z2: np.ndarray, mag: np.ndarray, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(Z2[:, 0], Z2[:, 1], c=mag, s=6, cmap="magma")
    fig.colorbar(sc, ax=ax)
    ax.set_title("Flux residual magnitude")
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
