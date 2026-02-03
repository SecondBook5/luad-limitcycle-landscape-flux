from __future__ import annotations

from typing import Any, Dict

import scanpy as sc


def select_hvg(adata, cfg: Dict[str, Any]) -> None:
    sc.pp.highly_variable_genes(adata, n_top_genes=int(cfg["preprocess"]["hvg_n"]))
    adata._inplace_subset_var(adata.var["highly_variable"])
