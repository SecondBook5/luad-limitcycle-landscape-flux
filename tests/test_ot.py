import numpy as np
import anndata as ad

from hpcs_model.dynamics.ot_coupling import compute_couplings


def test_ot_coupling_marginals(tmp_path):
    X = np.random.rand(10, 5)
    adata = ad.AnnData(X)
    adata.obsm["X_pca"] = np.random.rand(10, 3)
    adata.obs["timepoint"] = ["t0"] * 5 + ["t1"] * 5
    cfg = {"ot": {"epsilon": 0.1, "metric": "sqeuclidean"}}
    summaries = compute_couplings(adata, cfg, tmp_path)
    assert summaries
