from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import yaml


def repo_root(start: Path | None = None) -> Path:
    current = (start or Path.cwd()).resolve()
    for _ in range(8):
        if (current / "configs" / "params.yaml").exists():
            return current
        current = current.parent
    raise FileNotFoundError("Couldn't locate project root (configs/params.yaml).")


def load_config(root: Path | None = None) -> Dict[str, Any]:
    base = repo_root(root)
    with open(base / "configs" / "params.yaml", "r") as fh:
        params = yaml.safe_load(fh) or {}

    params.setdefault("paths", {})
    paths = params["paths"]
    paths.setdefault("anndata", "data/interim/mm_timecourse.h5ad")
    paths.setdefault("figures", "outputs/figures")
    paths.setdefault("tables", "outputs/tables")
    paths.setdefault("logs", "outputs/logs")
    paths.setdefault("raw_timecourse", "data/raw/scrna/timecourse_gse154989")
    paths.setdefault("raw_scatac", "data/raw/scatac/gse154965")
    paths.setdefault("raw_bulk_atac", "data/raw/bulk_atac/gse154966")

    params.setdefault("markers", {})
    params["markers"].setdefault("s_file", "configs/markers/S_mouse.txt")
    params["markers"].setdefault("g2m_file", "configs/markers/G2M_mouse.txt")

    params.setdefault("model", {})
    params["model"].setdefault("seed", 42)
    params["model"].setdefault("phase_bins", params.get("cells", {}).get("phase_bins", 60))

    return params
