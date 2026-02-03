from __future__ import annotations

from typing import Any, Dict

import pandas as pd
import scanpy as sc


def rank_markers(adata, cfg: Dict[str, Any]) -> pd.DataFrame:
    cluster_key = cfg["clustering"]["cluster_key"]
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")
    df = sc.get.rank_genes_groups_df(adata, None)
    return df
