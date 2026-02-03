from __future__ import annotations

from typing import Any, Dict

import numpy as np
import pandas as pd
import networkx as nx


def compute_plasticity(edge_df: pd.DataFrame) -> pd.DataFrame:
    if edge_df.empty:
        return pd.DataFrame()
    out_mass = edge_df.groupby("src_cluster")["mass_norm"].sum().rename("out_mass")
    in_mass = edge_df.groupby("dst_cluster")["mass_norm"].sum().rename("in_mass")
    clusters = sorted(set(edge_df["src_cluster"]).union(edge_df["dst_cluster"]))
    base = pd.DataFrame({"cluster": clusters})
    base = base.merge(out_mass, left_on="cluster", right_index=True, how="left")
    base = base.merge(in_mass, left_on="cluster", right_index=True, how="left")
    base = base.fillna(0.0)
    base["plasticity"] = base["out_mass"] + base["in_mass"]

    G = nx.DiGraph()
    for _, row in edge_df.iterrows():
        G.add_edge(row["src_cluster"], row["dst_cluster"], weight=row["mass_norm"])
    pagerank = nx.pagerank(G, weight="weight") if G.number_of_nodes() else {}
    base["pagerank"] = base["cluster"].map(pagerank).fillna(0.0)
    return base


def select_hpcs(plasticity_df: pd.DataFrame) -> str:
    if plasticity_df.empty:
        raise ValueError("Plasticity table is empty.")
    best = plasticity_df.sort_values(["plasticity", "pagerank"], ascending=False).iloc[0]
    return str(best["cluster"])


def bootstrap_hpcs(edge_df: pd.DataFrame, cfg: Dict[str, Any]) -> pd.DataFrame:
    n_iters = int(cfg["bootstrap"]["n_iters"])
    sample_fraction = float(cfg["bootstrap"]["sample_fraction"])
    results = []
    rng = np.random.default_rng(42)
    for _ in range(n_iters):
        sampled = edge_df.sample(frac=sample_fraction, replace=True, random_state=int(rng.integers(0, 1e9)))
        plasticity_df = compute_plasticity(sampled)
        hpcs = select_hpcs(plasticity_df)
        results.append(hpcs)
    counts = pd.Series(results).value_counts().rename_axis("cluster").reset_index(name="count")
    counts["frequency"] = counts["count"] / counts["count"].sum()
    return counts
