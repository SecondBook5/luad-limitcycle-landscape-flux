from __future__ import annotations

from typing import Any, Dict

import scanpy as sc


def run_leiden(adata, cfg: Dict[str, Any]) -> None:
    cluster_key = cfg["clustering"]["cluster_key"]
    sc.tl.leiden(adata, resolution=float(cfg["clustering"]["leiden_resolution"]), key_added=cluster_key)
