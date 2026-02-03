from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import anndata as ad

from hpcs_model.config import output_paths
from hpcs_model.io.validate_inputs import validate_adata


def load_anndata(cfg: Dict[str, Any]) -> ad.AnnData:
    input_path = Path(cfg["paths"]["input_h5ad"])
    if not input_path.exists():
        raise FileNotFoundError(f"Input AnnData not found: {input_path}")
    adata = ad.read_h5ad(input_path)
    validate_adata(adata)

    outputs = output_paths(cfg)
    outputs["root"].mkdir(parents=True, exist_ok=True)
    adata.write(outputs["root"] / "adata_raw.h5ad")
    return adata
