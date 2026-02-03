from __future__ import annotations

import scanpy as sc


def normalize_log1p(adata) -> None:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
