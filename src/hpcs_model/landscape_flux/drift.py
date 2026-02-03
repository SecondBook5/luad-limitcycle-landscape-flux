from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import scipy.sparse as sp


def drift_from_couplings(
    Z: np.ndarray,
    timepoints: np.ndarray,
    time_pairs: List[Tuple[str, str]],
    coupling_files: List[Path],
) -> np.ndarray:
    drift = np.full_like(Z, np.nan, dtype=float)
    for (t0, t1), file in zip(time_pairs, coupling_files):
        idx0 = np.where(timepoints == t0)[0]
        idx1 = np.where(timepoints == t1)[0]
        if len(idx0) == 0 or len(idx1) == 0:
            continue
        coupling = sp.load_npz(file).tocsr()
        row_sums = np.array(coupling.sum(axis=1)).ravel()
        row_sums[row_sums == 0] = 1.0
        weights = coupling.multiply(1.0 / row_sums[:, None])
        pred_next = weights @ Z[idx1]
        drift[idx0] = pred_next - Z[idx0]
    return drift
