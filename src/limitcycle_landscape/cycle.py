from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import FastICA, PCA

from limitcycle_landscape import cc_limitcycle
from limitcycle_landscape.config import load_config, repo_root
from limitcycle_landscape.markers import filter_present, load_markers


@dataclass
class CycleFitResult:
    basis: str
    center: np.ndarray
    radius: float
    gap_rad: float
    residual_median: float


def normalize_log1p(adata: ad.AnnData) -> None:
    import scanpy as sc

    X_counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    adata.X = X_counts.copy()
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)


def _zscore_submatrix(X: np.ndarray) -> np.ndarray:
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0) + 1e-8
    return (X - mu) / sd


def _dense_matrix(X) -> np.ndarray:
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def compute_module_scores(
    adata: ad.AnnData,
    s_markers: Iterable[str],
    g2m_markers: Iterable[str],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.Index]:
    genes = pd.Index(adata.var_names.astype(str))
    s_present = pd.Index(filter_present(s_markers, genes))
    g2m_present = pd.Index(filter_present(g2m_markers, genes))
    if len(s_present) < 5 or len(g2m_present) < 5:
        raise ValueError("Too few cell-cycle markers found in gene list.")

    use_genes = s_present.union(g2m_present)
    X_sub = _dense_matrix(adata[:, use_genes].X)
    Z = _zscore_submatrix(X_sub)

    col_index = pd.Index(list(use_genes))

    def avg_score(cols: pd.Index) -> np.ndarray:
        pos = col_index.get_indexer(cols)
        return Z[:, pos].mean(axis=1) if len(cols) else np.zeros(Z.shape[0])

    zS = avg_score(s_present)
    zG2M = avg_score(g2m_present)
    return zS, zG2M, Z, use_genes, col_index


def _ring_weights(zS: np.ndarray, zG2M: np.ndarray) -> np.ndarray:
    cycle_strength = np.sqrt(zS**2 + zG2M**2)
    p60, p95 = np.percentile(cycle_strength, [60, 95])
    w_raw = np.clip((cycle_strength - p60) / (p95 - p60 + 1e-8), 0, 1) ** 2
    w_fit = 0.15 + 0.85 * w_raw
    return (w_fit / (w_fit.mean() + 1e-12)).astype(float)


def fit_cycle_from_markers(
    root: Path | None = None,
    *,
    adata: ad.AnnData | None = None,
    config: Dict[str, object] | None = None,
    markers: Dict[str, Iterable[str]] | None = None,
    output_path: Path | None = None,
    top_frac: float = 0.05,
) -> Tuple[ad.AnnData, CycleFitResult]:
    base = repo_root(root)
    params = config or load_config(base)
    if adata is None:
        adata = ad.read_h5ad(base / params["paths"]["anndata"])

    normalize_log1p(adata)

    if markers is None:
        s_markers = load_markers(base / params["markers"]["s_file"])
        g2m_markers = load_markers(base / params["markers"]["g2m_file"])
    else:
        s_markers = markers.get("s", [])
        g2m_markers = markers.get("g2m", [])

    zS, zG2M, Z, use_genes, col_index = compute_module_scores(adata, s_markers, g2m_markers)
    adata.obs["zS"] = pd.Series(zS, index=adata.obs_names, dtype=float)
    adata.obs["zG2M"] = pd.Series(zG2M, index=adata.obs_names, dtype=float)
    adata.obsm["cycle_xy_modules"] = np.vstack([zS, zG2M]).T

    w_fit = _ring_weights(zS, zG2M)
    k = max(1, int(top_frac * len(zS)))
    idxS = np.argpartition(zS, -k)[-k:]
    idxM = np.argpartition(zG2M, -k)[-k:]

    def _cmean(a):
        return float(np.angle(np.mean(np.exp(1j * a))) % (2 * np.pi))

    def circle_score(gap, residuals):
        return float(gap) - 0.8 * float(np.median(residuals))

    cands = []

    Y = adata.obsm["cycle_xy_modules"]
    c, r = cc_limitcycle.fit_circle(Y, weights=w_fit)
    th = cc_limitcycle.orient_phases(
        cc_limitcycle.angles_from_center(Y, c), zS, zG2M, top_frac=top_frac
    )
    res = np.abs(np.linalg.norm(Y - c[None, :], axis=1) - r)
    gap = (_cmean(th[idxM]) - _cmean(th[idxS]) + 2 * np.pi) % (2 * np.pi)
    cands.append(("modules", circle_score(gap, res), Y, c, r, th, res, gap))

    XS = Z[:, col_index.get_indexer(pd.Index(filter_present(s_markers, use_genes)))]
    XM = Z[:, col_index.get_indexer(pd.Index(filter_present(g2m_markers, use_genes)))]

    pS = PCA(n_components=1, random_state=int(params["model"]["seed"]))
    pM = PCA(n_components=1, random_state=int(params["model"]["seed"]))
    Y_bp = np.vstack([pS.fit_transform(XS).ravel(), pM.fit_transform(XM).ravel()]).T
    c, r = cc_limitcycle.fit_circle(Y_bp, weights=w_fit)
    th = cc_limitcycle.orient_phases(
        cc_limitcycle.angles_from_center(Y_bp, c), zS, zG2M, top_frac=top_frac
    )
    res = np.abs(np.linalg.norm(Y_bp - c[None, :], axis=1) - r)
    gap = (_cmean(th[idxM]) - _cmean(th[idxS]) + 2 * np.pi) % (2 * np.pi)
    cands.append(("pc1S_pc1M", circle_score(gap, res), Y_bp, c, r, th, res, gap))

    pc1S = pS.fit_transform(XS).ravel()
    beta = (pc1S[:, None].T @ XM) / (pc1S[:, None].T @ pc1S[:, None] + 1e-12)
    XM_resid = XM - np.outer(pc1S, beta.ravel())
    pc1M_resid = pM.fit_transform(XM_resid).ravel()
    Y_br = np.vstack([pc1S, pc1M_resid]).T
    c, r = cc_limitcycle.fit_circle(Y_br, weights=w_fit)
    th = cc_limitcycle.orient_phases(
        cc_limitcycle.angles_from_center(Y_br, c), zS, zG2M, top_frac=top_frac
    )
    res = np.abs(np.linalg.norm(Y_br - c[None, :], axis=1) - r)
    gap = (_cmean(th[idxM]) - _cmean(th[idxS]) + 2 * np.pi) % (2 * np.pi)
    cands.append(("pc1S_vs_M|S", circle_score(gap, res), Y_br, c, r, th, res, gap))

    pca = PCA(n_components=2, random_state=int(params["model"]["seed"]))
    Y_p = pca.fit_transform(Z)
    c, r = cc_limitcycle.fit_circle(Y_p, weights=w_fit)
    th = cc_limitcycle.orient_phases(
        cc_limitcycle.angles_from_center(Y_p, c), zS, zG2M, top_frac=top_frac
    )
    res = np.abs(np.linalg.norm(Y_p - c[None, :], axis=1) - r)
    gap = (_cmean(th[idxM]) - _cmean(th[idxS]) + 2 * np.pi) % (2 * np.pi)
    cands.append(("pca_union", circle_score(gap, res), Y_p, c, r, th, res, gap))

    ica = FastICA(n_components=2, random_state=int(params["model"]["seed"]))
    Y_i = ica.fit_transform(Z)
    c, r = cc_limitcycle.fit_circle(Y_i, weights=w_fit)
    th = cc_limitcycle.orient_phases(
        cc_limitcycle.angles_from_center(Y_i, c), zS, zG2M, top_frac=top_frac
    )
    res = np.abs(np.linalg.norm(Y_i - c[None, :], axis=1) - r)
    gap = (_cmean(th[idxM]) - _cmean(th[idxS]) + 2 * np.pi) % (2 * np.pi)
    cands.append(("ica_union", circle_score(gap, res), Y_i, c, r, th, res, gap))

    best_name, best_score, Ybest, cbest, rbest, thbest, resbest, gapbest = max(
        cands, key=lambda x: x[1]
    )

    adata.obsm["cycle_xy"] = Ybest
    adata.obs["theta"] = pd.Series(thbest, index=adata.obs_names, dtype=float)
    adata.obs["radius"] = pd.Series(
        np.linalg.norm(Ybest - cbest[None, :], axis=1),
        index=adata.obs_names,
        dtype=float,
    )
    adata.uns["cycle_fit"] = {
        "basis": best_name,
        "center": [float(cbest[0]), float(cbest[1])],
        "radius": float(rbest),
        "gap_S_to_M_rad": float(gapbest),
        "residual_median": float(np.median(resbest)),
    }

    result = CycleFitResult(
        basis=best_name,
        center=np.array(cbest, dtype=float),
        radius=float(rbest),
        gap_rad=float(gapbest),
        residual_median=float(np.median(resbest)),
    )

    if output_path is None:
        output_path = base / params["paths"]["anndata"]
    adata.write(output_path)
    return adata, result


def bin_atac_by_theta(
    adata: ad.AnnData,
    atac_layer: str = "X_atac",
    bins: int = 60,
) -> Dict[str, object]:
    if "theta" not in adata.obs:
        raise ValueError("Missing theta in adata.obs. Run fit_cycle_from_markers first.")
    if atac_layer not in adata.layers:
        raise ValueError(f"ATAC layer '{atac_layer}' not found in adata.layers.")

    theta = adata.obs["theta"].to_numpy()
    bin_edges = np.linspace(0.0, 2 * np.pi, bins + 1)
    bin_ids = np.digitize(theta, bin_edges) - 1
    bin_ids = np.clip(bin_ids, 0, bins - 1)

    X = adata.layers[atac_layer]
    if sp.issparse(X):
        X = X.tocsr()

    binned = []
    counts = []
    for b in range(bins):
        idx = np.where(bin_ids == b)[0]
        counts.append(len(idx))
        if len(idx) == 0:
            if sp.issparse(X):
                binned.append(sp.csr_matrix((1, X.shape[1]), dtype=np.float32))
            else:
                binned.append(np.zeros((1, X.shape[1]), dtype=np.float32))
            continue
        if sp.issparse(X):
            sub = X[idx]
            mean = sub.mean(axis=0)
            binned.append(sp.csr_matrix(mean))
        else:
            binned.append(np.mean(X[idx], axis=0, keepdims=True))

    if sp.issparse(X):
        binned_matrix = sp.vstack(binned).tocsr()
    else:
        binned_matrix = np.vstack(binned)

    payload = {
        "bin_edges": bin_edges,
        "counts": np.array(counts, dtype=int),
        "matrix": binned_matrix,
    }
    adata.uns["atac_pseudobulk"] = payload
    return payload
