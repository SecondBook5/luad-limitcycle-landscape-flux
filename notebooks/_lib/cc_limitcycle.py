# cc_limitcycle.py

from __future__ import annotations
import math
from typing import Tuple, Iterable
import numpy as np
from numpy.typing import ArrayLike
from sklearn.mixture import GaussianMixture

def zscore_rows(X: np.ndarray) -> np.ndarray:
    """Return per-gene z-scored matrix (cells × genes)."""
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0) + 1e-8
    return (X - mu) / sd

def fit_circle(Y: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Least-squares circle fit (Taubin-style).
    Y: (n×2) array of points [S_score, G2M_score].
    Returns center c (2,) and radius r.
    """
    x = Y[:, 0]; y = Y[:, 1]
    A = np.c_[x, y, np.ones_like(x)]
    b = x**2 + y**2
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    a, bcoef, d = sol
    c = np.array([a/2.0, bcoef/2.0], dtype=float)
    r = float(np.sqrt((c**2).sum() + d))
    return c, r

def angles_from_center(Y: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Angles θ in [0, 2π) from points Y relative to center c."""
    V = Y - c[None, :]
    th = np.arctan2(V[:, 1], V[:, 0]) % (2.0 * np.pi)
    return th

def orient_phases(theta: np.ndarray, S: np.ndarray, G2M: np.ndarray) -> np.ndarray:
    """
    Fix rotation/reflection: place high-S near 0 and ensure S→G2/M direction is positive.
    """
    # rotation: center the mean θ of top-10% S at 0
    k = max(1, int(0.10 * len(S)))
    idxS = np.argsort(S)[-k:]
    rot = np.angle(np.mean(np.exp(1j * theta[idxS])))
    th = (theta - rot) % (2*np.pi)

    # reflection: ensure mean(G2M) > mean(S) along positive direction
    idxM = np.argsort(G2M)[-k:]
    muS = np.angle(np.mean(np.exp(1j * th[idxS])))
    muM = np.angle(np.mean(np.exp(1j * th[idxM])))
    # unwrap small arc distance S→M
    d = (muM - muS + 2*np.pi) % (2*np.pi)
    if d > np.pi:
        th = (-th) % (2*np.pi)
    return th

def fit_periodic(theta: np.ndarray, y: ArrayLike, K: int = 1, ridge: float = 0.0):
    """
    Low-harmonic Fourier regression y ~ [1, cosθ, sinθ, (cos2θ, sin2θ)].
    Returns dict with coeffs, R2, amplitude (first harmonic), and peak phase.
    """
    theta = np.asarray(theta)
    y = np.asarray(y).reshape(-1, 1)
    Phi = [np.ones_like(theta)]
    Phi += [np.cos(theta), np.sin(theta)]
    if K == 2:
        Phi += [np.cos(2*theta), np.sin(2*theta)]
    Phi = np.vstack(Phi).T  # n × p
    A = Phi.T @ Phi
    if ridge > 0:
        A = A + ridge * np.eye(A.shape[0])
    beta = np.linalg.solve(A, Phi.T @ y)
    yhat = Phi @ beta
    ssr = float(((yhat - y)**2).sum())
    sst = float(((y - y.mean())**2).sum() + 1e-12)
    R2 = 1.0 - ssr/sst
    A1 = float(beta[1]); C1 = float(beta[2])
    amp = float(np.sqrt(A1**2 + C1**2))
    peak = float(np.arctan2(-C1, A1) % (2*np.pi))
    return {"beta": beta.ravel(), "R2": R2, "amplitude": amp, "peak": peak}

def annulus_disk_gmm(radius: np.ndarray, seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
    """
    1D GMM on radii: component with higher mean = ring (cycling).
    Returns (proba_ring, label_ring_binary).
    """
    r = radius.reshape(-1, 1)
    gm = GaussianMixture(n_components=2, random_state=seed).fit(r)
    means = gm.means_.ravel()
    ring_comp = int(np.argmax(means))
    proba = gm.predict_proba(r)[:, ring_comp]
    call = (proba >= 0.5).astype(int)
    return proba, call

def circular_w1_hist(p: np.ndarray, q: np.ndarray) -> float:
    """
    W1 on the circle for histogram densities p, q (length m, sum to 1).
    Uses the 'cut' trick: min over all rotations of line W1 of their cumulative diffs.
    """
    p = np.asarray(p, float); q = np.asarray(q, float)
    assert p.shape == q.shape
    m = len(p)
    # cumulative difference on line
    def w1_linear(a, b):
        diff_cum = np.cumsum(a - b)
        return np.sum(np.abs(diff_cum)) / m
    # try all circular cuts (via rotations)
    best = np.inf
    for k in range(m):
        best = min(best, w1_linear(np.roll(p, k), q))
    return float(best)

def piecewise_warp(theta: np.ndarray, a0: float, a1: float, stretch: float = 1.25) -> np.ndarray:
    """
    Stretch arc [a0, a1] by factor 'stretch', compress remainder to keep total 2π.
    Returns warped angles in [0, 2π).
    """
    L = 2*np.pi
    a0 = a0 % L; a1 = a1 % L
    if a1 <= a0:
        a1 += L
    g1 = a1 - a0
    g2 = L - g1
    g1p = min(L - 1e-6, g1 * stretch)
    g2p = L - g1p
    s1 = g1p / g1
    s2 = g2p / g2
    th = theta.copy().astype(float)
    # map to unwrapped, then rewrap after piecewise linear scaling
    t = th.copy()
    t[th < a0] += L
    # three regions: before a0, within [a0,a1], after a1
    out = np.empty_like(t)
    # region 1: [a0,a1]
    m = (t >= a0) & (t <= a1)
    out[m] = (t[m] - a0) * s1 + a0
    # region 2: (a1, a0+2π]
    m2 = (t > a1)
    out[m2] = (t[m2] - a1) * s2 + a0 + g1p
    # region 3: [a0-2π, a0)
    m3 = (t < a0)
    out[m3] = (t[m3] - (a0 - L)) * s2 + a0 - g2p
    return (out % L)
