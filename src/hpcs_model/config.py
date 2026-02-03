from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import yaml


def load_config(path: str | Path) -> Dict[str, Any]:
    cfg_path = Path(path)
    with cfg_path.open("r") as fh:
        cfg = yaml.safe_load(fh) or {}

    cfg.setdefault("paths", {})
    paths = cfg["paths"]
    paths.setdefault("data_root", "data")
    paths.setdefault("output_root", "outputs")
    paths.setdefault("input_h5ad", "data/raw/adata_raw.h5ad")

    cfg.setdefault("preprocess", {})
    preprocess = cfg["preprocess"]
    preprocess.setdefault("min_cells", 200)
    preprocess.setdefault("min_genes", 200)
    preprocess.setdefault("hvg_n", 2000)
    preprocess.setdefault("pca_n", 50)
    preprocess.setdefault("neighbors_k", 15)
    preprocess.setdefault("embedding", "umap")

    cfg.setdefault("clustering", {})
    clustering = cfg["clustering"]
    clustering.setdefault("leiden_resolution", 1.0)
    clustering.setdefault("cluster_key", "cluster")

    cfg.setdefault("ot", {})
    ot_cfg = cfg["ot"]
    ot_cfg.setdefault("epsilon", 0.05)
    ot_cfg.setdefault("metric", "sqeuclidean")
    ot_cfg.setdefault("seed", 42)

    cfg.setdefault("graph", {})
    graph_cfg = cfg["graph"]
    graph_cfg.setdefault("edge_threshold", 0.0)

    cfg.setdefault("landscape_flux", {})
    lf_cfg = cfg["landscape_flux"]
    lf_cfg.setdefault("density_k", 30)
    lf_cfg.setdefault("grid_size", 50)
    lf_cfg.setdefault("eps", 1e-8)

    cfg.setdefault("ml", {})
    ml_cfg = cfg["ml"]
    ml_cfg.setdefault("hidden_dim", 128)
    ml_cfg.setdefault("epochs", 50)
    ml_cfg.setdefault("lr", 1e-3)
    ml_cfg.setdefault("batch_size", 256)
    ml_cfg.setdefault("seed", 42)

    cfg.setdefault("bootstrap", {})
    bootstrap = cfg["bootstrap"]
    bootstrap.setdefault("n_iters", 50)
    bootstrap.setdefault("sample_fraction", 0.8)

    return cfg


def output_paths(cfg: Dict[str, Any]) -> Dict[str, Path]:
    output_root = Path(cfg["paths"]["output_root"])
    return {
        "root": output_root,
        "figures": output_root / "figures",
        "tables": output_root / "tables",
        "models": output_root / "models",
        "logs": output_root / "logs",
    }
