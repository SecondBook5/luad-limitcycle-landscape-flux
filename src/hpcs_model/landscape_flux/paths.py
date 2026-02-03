from __future__ import annotations

from typing import List, Tuple

import networkx as nx
import numpy as np


def dominant_paths(Z2: np.ndarray, U: np.ndarray, k: int = 3) -> List[Tuple[int, int, List[int]]]:
    basins = np.argsort(U)[:k]
    G = nx.Graph()
    for i in range(len(Z2)):
        for j in range(i + 1, len(Z2)):
            if np.linalg.norm(Z2[i] - Z2[j]) < 0.5:
                G.add_edge(i, j, weight=float(U[j]))
    paths = []
    for i in basins:
        for j in basins:
            if i >= j:
                continue
            if nx.has_path(G, i, j):
                path = nx.shortest_path(G, i, j, weight="weight")
                paths.append((int(i), int(j), path))
    return paths
