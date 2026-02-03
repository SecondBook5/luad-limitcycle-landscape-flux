import anndata as ad
import numpy as np
import pandas as pd

from hpcs_model.io.validate_inputs import validate_adata


def test_validate_missing_columns():
    X = np.random.rand(5, 3)
    adata = ad.AnnData(X)
    adata.obs["timepoint"] = ["t0"] * 5
    adata.obs["genotype"] = ["T"] * 5
    try:
        validate_adata(adata)
    except ValueError as exc:
        assert "sample_id" in str(exc)
    else:
        raise AssertionError("Expected missing sample_id to raise.")


def test_validate_timepoints_ordered():
    X = np.random.rand(5, 3)
    adata = ad.AnnData(X)
    adata.obs["timepoint"] = ["t0"] * 5
    adata.obs["genotype"] = ["T"] * 5
    adata.obs["sample_id"] = ["s1"] * 5
    try:
        validate_adata(adata)
    except ValueError as exc:
        assert "At least two timepoints" in str(exc)
