from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import scipy.sparse as sp


def aggregate_to_cluster_graph(
    adata, cfg: Dict[str, Any], coupling_summaries: List[Dict[str, Any]], output_dir: Path
) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)
    cluster_key = cfg["clustering"]["cluster_key"]
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
    edges = []
    for summary in coupling_summaries:
        path = output_dir / summary["file"]
        coupling = sp.load_npz(path).tocoo()
        time_t = summary["time_t"]
        time_t1 = summary["time_t1"]
        src_clusters = clusters[adata.obs["timepoint"] == time_t]
        dst_clusters = clusters[adata.obs["timepoint"] == time_t1]

        src_idx = src_clusters[coupling.row]
        dst_idx = dst_clusters[coupling.col]
        df = pd.DataFrame(
            {
                "src_cluster": src_idx,
                "dst_cluster": dst_idx,
                "mass": coupling.data,
            }
        )
        agg = df.groupby(["src_cluster", "dst_cluster"], as_index=False)["mass"].sum()
        agg["time_t"] = time_t
        agg["time_t1"] = time_t1
        edges.append(agg)

    if not edges:
        return pd.DataFrame()
    edge_df = pd.concat(edges, ignore_index=True)
    edge_df["mass_norm"] = edge_df.groupby(["time_t", "src_cluster"])["mass"].transform(
        lambda x: x / x.sum()
    )
    edge_df = edge_df[edge_df["mass"] >= float(cfg["graph"]["edge_threshold"])]
    edge_df.to_csv(output_dir / "cluster_transition_edges.csv", index=False)
    return edge_df
