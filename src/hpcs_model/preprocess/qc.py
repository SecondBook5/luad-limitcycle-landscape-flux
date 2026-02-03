from __future__ import annotations

from typing import Any, Dict

import scanpy as sc


def qc_filter(adata, cfg: Dict[str, Any]) -> None:
    sc.pp.filter_cells(adata, min_genes=int(cfg["preprocess"]["min_genes"]))
    sc.pp.filter_genes(adata, min_cells=int(cfg["preprocess"]["min_cells"]))
