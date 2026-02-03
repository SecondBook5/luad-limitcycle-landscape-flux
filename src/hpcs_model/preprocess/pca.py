from __future__ import annotations

from typing import Any, Dict

import scanpy as sc


def run_pca(adata, cfg: Dict[str, Any]) -> None:
    sc.tl.pca(adata, n_comps=int(cfg["preprocess"]["pca_n"]))
    sc.pp.neighbors(adata, n_neighbors=int(cfg["preprocess"]["neighbors_k"]))
    embedding = cfg["preprocess"].get("embedding", "umap")
    if embedding == "umap":
        sc.tl.umap(adata)
