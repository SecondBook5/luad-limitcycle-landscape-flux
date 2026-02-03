from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
import ot


def _time_pairs(timepoints: pd.Categorical) -> List[Tuple[str, str]]:
    categories = list(timepoints.categories)
    return list(zip(categories[:-1], categories[1:]))


def compute_couplings(adata, cfg: Dict[str, Any], output_dir: Path) -> List[Dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    timepoints = adata.obs["timepoint"]
    pairs = _time_pairs(timepoints)
    Z = adata.obsm["X_pca"]
    summaries = []
    for t0, t1 in pairs:
        idx0 = np.where(timepoints == t0)[0]
        idx1 = np.where(timepoints == t1)[0]
        if len(idx0) == 0 or len(idx1) == 0:
            continue
        Z0 = Z[idx0]
        Z1 = Z[idx1]
        cost = ot.dist(Z0, Z1, metric=cfg["ot"]["metric"]) ** 2
        a = np.ones(len(idx0)) / len(idx0)
        b = np.ones(len(idx1)) / len(idx1)
        coupling = ot.sinkhorn(a, b, cost, reg=float(cfg["ot"]["epsilon"]))
        coupling_sparse = sp.csr_matrix(coupling)

        fname = f"coupling_{t0}_to_{t1}.npz"
        sp.save_npz(output_dir / fname, coupling_sparse)
        summary = {
            "time_t": str(t0),
            "time_t1": str(t1),
            "n_source": int(len(idx0)),
            "n_target": int(len(idx1)),
            "row_sum_mean": float(np.mean(coupling.sum(axis=1))),
            "col_sum_mean": float(np.mean(coupling.sum(axis=0))),
            "file": fname,
        }
        summaries.append(summary)
    return summaries
