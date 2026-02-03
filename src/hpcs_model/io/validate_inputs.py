from __future__ import annotations

from typing import Iterable

import anndata as ad
import pandas as pd


REQUIRED_COLUMNS = ("timepoint", "genotype", "sample_id")


def _require_columns(obs: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [col for col in columns if col not in obs.columns]
    if missing:
        raise ValueError(f"Missing required obs columns: {', '.join(missing)}")


def validate_adata(adata: ad.AnnData) -> None:
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("AnnData must contain cells and genes.")
    _require_columns(adata.obs, REQUIRED_COLUMNS)

    timepoint = adata.obs["timepoint"]
    if not pd.api.types.is_categorical_dtype(timepoint):
        ordered = pd.Categorical(timepoint, ordered=True)
        adata.obs["timepoint"] = ordered
    if adata.obs["timepoint"].nunique() < 2:
        raise ValueError("At least two timepoints are required for OT couplings.")
