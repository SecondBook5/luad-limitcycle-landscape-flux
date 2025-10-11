# cc_limitcycle.py
# Utilities for the limit-cycle (cell-cycle) geometry, densities, barriers, and descriptive flux.
# Dependencies: numpy, sklearn (GaussianMixture). No SciPy required.

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, List

import numpy as np
from numpy.typing import ArrayLike
from sklearn.mixture import GaussianMixture


# ---------------------------------------------------------------------
# Basic helpers (angles, zscores)
# ---------------------------------------------------------------------

def zscore_rows(X: np.ndarray) -> np.ndarray:
    """
    Per-gene z-score of a cell-by-gene matrix.
    Args:
        X: (n_cells × n_genes) array of normalized expression values.
    Returns:
        Z: (n_cells × n_genes) array where each gene column is standardized:
           Z[:, g] = (X[:, g] - mean_g) / (std_g + 1e-8)
    Notes:
        Adds 1e-8 to std to prevent division-by-zero if a gene is constant.
    """
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0) + 1e-8  # numerical safety
    return (X - mu) / sd


def wrap_angle(theta: np.ndarray) -> np.ndarray:
    """
    Wrap angles into [0, 2π).
    Args:
        theta: array-like angles in radians (any real values).
    Returns:
        Wrapped angles in [0, 2π).
    """
    return np.mod(theta, 2.0 * np.pi)


def circ_mean(theta: np.ndarray, weights: Optional[np.ndarray] = None) -> float:
    """
    Circular mean of angles in [0, 2π).
    Args:
        theta: 1D array of angles (radians).
        weights: optional non-negative weights per angle (same length as theta).
    Returns:
        Mean angle (radians) in [0, 2π).
    Method:
        Uses complex representation: mean_angle = arg(Σ w_k e^{i θ_k}).
    """
    theta = np.asarray(theta)
    if weights is None:
        weights = np.ones_like(theta)
    wsum = float(np.sum(weights))
    if not np.isfinite(wsum) or wsum <= 0:
        weights = np.ones_like(theta)
        wsum = float(len(theta))
    w = weights / (wsum + 1e-12)
    v = np.exp(1j * theta)
    mu = np.angle(np.sum(w * v))
    return float(mu % (2.0 * np.pi))


def theta_grid(nbins: int = 60) -> np.ndarray:
    """
    Build a uniform θ grid in [0, 2π) with nbins points.
    """
    return np.linspace(0.0, 2.0 * np.pi, nbins, endpoint=False)


# ---------------------------------------------------------------------
# Circle fit & phase mapping
# ---------------------------------------------------------------------

def fit_circle(Y: np.ndarray, weights: Optional[np.ndarray] = None) -> Tuple[np.ndarray, float]:
    """
    Algebraic least-squares circle fit to 2D points.
    Args:
        Y: (n × 2) array; columns are [S_score, G2M_score] (or any 2D coordinates).
        weights: optional non-negative weights per point (length n).
    Returns:
        c: (2,) center [cx, cy]
        r: scalar radius
    Derivation:
        Solve linear system [2x, 2y, 1] @ [cx, cy, t] = x^2 + y^2 in least squares;
        then r = sqrt(t + cx^2 + cy^2).
    """
    Y = np.asarray(Y, float)
    x = Y[:, 0]
    y = Y[:, 1]
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]
    b = x**2 + y**2
    if weights is not None:
        # Weighted least squares should scale rows by sqrt(weights)
        w = np.asarray(weights).reshape(-1, 1)
        s = np.sqrt(np.maximum(w, 0.0))
        Aw = A * s
        bw = b * s.ravel()
        sol, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    else:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy, t = sol
    c = np.array([cx, cy], dtype=float)
    r = float(np.sqrt(max(t + cx**2 + cy**2, 0.0)))
    return c, r


def angles_from_center(Y: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Convert 2D points to angles relative to a center.
    Args:
        Y: (n × 2) array of points.
        c: (2,) center (cx, cy).
    Returns:
        theta: (n,) angles in [0, 2π), computed by atan2(y - cy, x - cx).
    """
    V = Y - c[None, :]
    th = np.arctan2(V[:, 1], V[:, 0])
    return wrap_angle(th)


def orient_phases(theta: np.ndarray, S: np.ndarray, G2M: np.ndarray, top_frac: float = 0.10) -> np.ndarray:
    """
    Fix rotation and reflection of phases to match biology.
    Goal:
        - Put S-high cells near 0 radians.
        - Ensure G2/M-high cells lie ahead of S along the positive direction (S → G2/M).
    Args:
        theta: (n,) raw angles from angles_from_center.
        S: (n,) S-phase scores (used to find “high S” cells).
        G2M: (n,) G2/M scores (used to find “high G2/M” cells).
        top_frac: fraction (0–1) of top-scoring cells to define S-high and G2M-high sets.
    Returns:
        th: (n,) reoriented angles in [0, 2π).
    Steps:
        1) Rotate so the circular mean of S-high is at 0.
        2) If G2/M mean falls “behind” S (distance > π), reflect: th ← -th.
    """
    k = max(1, int(top_frac * len(S)))
    idxS = np.argpartition(S, -k)[-k:]
    idxM = np.argpartition(G2M, -k)[-k:]

    # 1) rotation
    muS = circ_mean(theta[idxS])
    th = wrap_angle(theta - muS)

    # 2) reflection if needed (ensure G2/M lies ahead of S along positive arc)
    muM = circ_mean(th[idxM])
    d = (muM - 0.0) % (2.0 * np.pi)  # arc distance from S(=0) to G2/M
    if d > np.pi:
        th = wrap_angle(-th)
    return th


# ---------------------------------------------------------------------
# Ring vs disk (cycling vs quiescent) classification
# ---------------------------------------------------------------------

def annulus_disk_gmm(
    radius: np.ndarray,
    seed: int = 42,
    covariance_type: str = "full",
    n_init: int = 3
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    """
    Two-component 1D GMM on radii; the higher-mean component is the ring (cycling).
    Args:
        radius: (n,) distances from circle center.
        seed: random seed for reproducibility.
        covariance_type: sklearn GMM covariance_type (default "full").
        n_init: number of initializations for GMM (stability).
    Returns:
        proba_ring: (n,) posterior probability of being on the ring.
        call: (n,) binary call (1 if proba_ring ≥ 0.5 else 0).
        diagnostics: dict with component means, stds, and mixture weights.
    Notes:
        This is intentionally 1D to keep the decision interpretable and robust.
    """
    r = radius.reshape(-1, 1)
    gm = GaussianMixture(n_components=2, covariance_type=covariance_type, random_state=seed, n_init=n_init)
    gm.fit(r)
    means = gm.means_.ravel()
    ring_comp = int(np.argmax(means))  # larger mean radius ⇒ ring
    proba = gm.predict_proba(r)[:, ring_comp]
    call = (proba >= 0.5).astype(int)

    # Extract per-component std (√variance) and weights for logging/reporting.
    def _std_of_comp(k: int) -> float:
        cov = gm.covariances_[k]
        return float(np.sqrt(cov).ravel()[0])

    diags = {
        "mu_ring": float(means[ring_comp]),
        "mu_disk": float(means[1 - ring_comp]),
        "sigma_ring": _std_of_comp(ring_comp),
        "sigma_disk": _std_of_comp(1 - ring_comp),
        "weight_ring": float(gm.weights_[ring_comp]),
        "weight_disk": float(gm.weights_[1 - ring_comp]),
    }
    return proba, call, diags


# ---------------------------------------------------------------------
# Circular KDE (von Mises kernel) and bandwidth selection
# ---------------------------------------------------------------------

def vm_kernel(delta: np.ndarray, kappa: float) -> np.ndarray:
    """
    von Mises kernel: K_kappa(Δ) = exp(kappa * cos(Δ)) / (2π I0(kappa)).
    Args:
        delta: array of angle differences (radians).
        kappa: concentration (>0). Larger kappa ⇒ narrower kernel.
    Returns:
        Kernel values with same shape as delta.
    """
    # numpy.i0 is the modified Bessel function of the first kind, order 0.
    return np.exp(kappa * np.cos(delta)) / (2.0 * np.pi * np.i0(kappa))


def vm_kde(theta_grid: np.ndarray, theta: np.ndarray, weights: Optional[np.ndarray], kappa: float) -> np.ndarray:
    """
    Circular KDE on a grid using von Mises kernels.
    Args:
        theta_grid: (G,) grid of angles where density will be evaluated.
        theta: (N,) sample angles.
        weights: optional (N,) non-negative weights. If None, uses uniform weights.
        kappa: concentration parameter for the kernel.
    Returns:
        rho: (G,) estimated density on theta_grid; normalized so that ∑ rho * Δθ = 1.
    """
    thg = theta_grid.reshape(-1, 1)  # G × 1
    th = theta.reshape(1, -1)        # 1 × N
    delta = thg - th                 # G × N
    K = vm_kernel(delta, kappa)      # G × N

    if weights is None:
        w = np.ones((1, theta.shape[0]))
    else:
        w = weights.reshape(1, -1)

    rho = (K * w).sum(axis=1) / (w.sum() + 1e-12)
    # normalize so sum(rho) * Δθ = 1 for cross-stage comparability
    dtheta = 2.0 * np.pi / len(theta_grid)
    rho = rho / (np.sum(rho) * dtheta + 1e-12)
    return rho.ravel()


def select_kappa_cv(theta: np.ndarray, weights: Optional[np.ndarray] = None, grid: Optional[np.ndarray] = None) -> float:
    """
    Select von Mises concentration kappa by leave-one-out (LOO) log-likelihood.
    Args:
        theta: (N,) angles.
        weights: optional (N,) weights (non-negative). Normalized internally.
        grid: optional array of kappa candidates (e.g., np.linspace(1, 60, 30)).
    Returns:
        kappa*: scalar that maximizes LOO log-likelihood on the provided grid.
    Complexity:
        O(N^2 * |grid|). Use on 1D phase data; typical N is manageable.
    """
    th = theta.reshape(1, -1)                 # 1 × N
    if weights is None:
        w = np.ones_like(theta)
    else:
        w = weights
    w = w / (w.sum() + 1e-12)

    if grid is None:
        grid = np.linspace(1.0, 60.0, 30)     # default search range

    best_ll = -np.inf
    best_k = float(grid[0])

    for kappa in grid:
        # LOO density for each point i:
        # p_i ∝ Σ_{j≠i} w_j * K(θ_i - θ_j)
        delta = th.T - th                      # N × N (i - j)
        K = vm_kernel(delta, kappa)            # N × N
        np.fill_diagonal(K, 0.0)               # remove self-kernel
        p = (K * w).sum(axis=1) / (1.0 - w + 1e-12)
        p = np.clip(p, 1e-300, None)           # avoid log(0)
        ll = float(np.sum(np.log(p)))
        if ll > best_ll:
            best_ll = ll
            best_k = float(kappa)
    return best_k


# ---------------------------------------------------------------------
# Cylinder KDE & ring-restricted potential
# ---------------------------------------------------------------------

def silverman_bandwidth(x: np.ndarray) -> float:
    """
    Silverman’s rule-of-thumb bandwidth for 1D Gaussian KDE.
    Args:
        x: 1D data array.
    Returns:
        h: scalar bandwidth ≈ 1.06 * std(x) * n^{-1/5}
    """
    x = np.asarray(x, float)
    sd = np.std(x)
    n = len(x)
    return 1.06 * sd * (n ** (-1.0 / 5.0) + 1e-12)


def ring_density_and_potential(
    theta_grid: np.ndarray,
    theta: np.ndarray,
    radius: np.ndarray,
    weights: Optional[np.ndarray],
    kappa: float,
    r_star: float,
    tau: float,
    h_r: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Ring-restricted steady density P_ring(θ) and potential U(θ) = -ln P_ring(θ).
    Args:
        theta_grid: (G,) grid of angles.
        theta: (N,) sample angles for cycling cells.
        radius: (N,) radii of the same cells.
        weights: optional (N,) weights (e.g., cycling probabilities^gamma).
        kappa: concentration for von Mises kernel (phase bandwidth).
        r_star: single representative ring radius (e.g., from circle fit).
        tau: ring thickness (std) to define annulus; use σ_ring from GMM.
        h_r: optional radial KDE bandwidth; if None, uses Silverman’s rule on radii.
    Returns:
        (P_ring, U):
            P_ring: (G,) density restricted to annulus around r_star, normalized to integrate to 1.
            U:      (G,) potential = -ln P_ring (clipped for numerical safety).
    Method:
        Product kernel: von Mises in θ and Gaussian in r.
        The ring restriction uses ψ(r) = N(r_star | r, τ^2), which integrates out analytically:
            ∫ N(r | r_i, h_r^2) ψ(r) dr = N(r_star | r_i, h_r^2 + τ^2).
        Thus:
            P_ring(θ) ∝ Σ_i w_i * K_kappa(θ - θ_i) * N(r_star | r_i, h_r^2 + τ^2).
    """
    if h_r is None:
        h_r = silverman_bandwidth(radius)
    var_eff = h_r**2 + tau**2

    # Precompute the radial factor per cell (constant across θ once r_star is fixed).
    rad_fac = (1.0 / np.sqrt(2.0 * np.pi * var_eff)) * np.exp(-0.5 * (r_star - radius) ** 2 / var_eff)

    # Combine optional weights with radial factor
    if weights is None:
        w = rad_fac
    else:
        w = rad_fac * weights

    # Phase KDE with effective per-point weights (w)
    p_ring = vm_kde(theta_grid, theta, w, kappa)

    # Normalize explicitly so ∫ P_ring dθ = 1
    dtheta = 2.0 * np.pi / len(theta_grid)
    p_ring = p_ring / (np.sum(p_ring) * dtheta + 1e-12)
    p_ring = np.clip(p_ring, 1e-300, None)  # avoid -inf in log

    U = -np.log(p_ring)
    return p_ring, U


def center_potential(
    radius: np.ndarray,
    weights: Optional[np.ndarray] = None,
    h_r: Optional[float] = None
) -> float:
    """
    U_center from P_ss(θ, r=0) averaged over θ.
    For the product kernel, ∫ K_kappa dθ = 1, so:
      P_center ∝ Σ_i w_i * N(0 | r_i, h_r^2)
    Returns:
      U_center = −ln P_center.
    """
    r = np.asarray(radius, float)
    if h_r is None:
        h_r = silverman_bandwidth(r)
    w = np.ones_like(r) if weights is None else np.asarray(weights)
    var = h_r ** 2
    gauss0 = (1.0 / np.sqrt(2.0 * np.pi * var)) * np.exp(-0.5 * (r ** 2) / var)
    p_center = float(np.sum(w * gauss0) / (np.sum(w) + 1e-12))
    p_center = max(p_center, 1e-300)
    return float(-np.log(p_center))


def center_barrier(U_theta: np.ndarray, U_center: float) -> float:
    """
    Barrier_center = U_center − min_θ U(θ).
    """
    return float(U_center - float(np.min(U_theta)))


# ---------------------------------------------------------------------
# Barrier finding on the ring
# ---------------------------------------------------------------------

@dataclass
class Barrier:
    """Container for a ring-wise barrier (saddle) and its adjacent basin minima."""
    theta_saddle: float
    theta_basin_left: float
    theta_basin_right: float
    height_left: float
    height_right: float


def find_barriers(
    theta_grid: np.ndarray,
    U: np.ndarray,
    smooth: bool = True,
) -> List[Barrier]:
    """
    Identify ring-wise saddles (local maxima of U) and adjacent basin minima.
    Args:
        theta_grid: (G,) uniform grid of angles (assumed equally spaced).
        U: (G,) potential on the ring.
        smooth: if True, apply a small circular moving average to reduce spurious extrema.
    Returns:
        List of Barrier objects with saddle angle, neighboring basins, and heights.
    Method:
        - Optionally smooth U with a periodic window (radius w=3).
        - Use discrete periodic first/second derivatives to flag maxima.
        - For each maximum, scan left and right to nearest local minima.
    Caveats:
        This is a simple, robust detector; for noisy profiles, consider spline
        smoothing or persistence-based peak picking as refinements.
    """
    th = theta_grid
    u = U.copy()

    if smooth:
        # Periodic moving average with window radius w (total window size = 2w+1)
        w = 3
        kern = np.ones(2 * w + 1, float)
        u = np.convolve(np.r_[u[-w:], u, u[:w]], kern, mode="same")[w:-w] / kern.sum()

    # Uniform-grid assertion
    dtheta = th[1] - th[0]
    if not np.allclose(np.diff(th), dtheta):
        raise ValueError("theta_grid must be a uniform grid in [0, 2π).")

    # Discrete periodic derivatives (central differences)
    up = np.roll(u, -1)
    um = np.roll(u, 1)
    du = (up - um) / (2.0 * dtheta)
    d2u = (up - 2 * u + um) / (dtheta ** 2)

    eps = 1e-12
    saddles: List[Barrier] = []
    for i in range(len(u)):
        i_prev = (i - 1) % len(u)
        # Sign-change test for a local maximum (saddle along the ring)
        if (du[i_prev] > eps) and (du[i] <= eps) and (d2u[i] < -eps):
            # Find nearest minima to the left
            j = i
            while True:
                j = (j - 1) % len(u)
                if (u[(j - 1) % len(u)] >= u[j]) and (u[(j + 1) % len(u)] >= u[j]):
                    basin_left = j
                    break
                if j == i:  # fallback if flat
                    basin_left = j
                    break

            # Find nearest minima to the right
            k = i
            while True:
                k = (k + 1) % len(u)
                if (u[(k - 1) % len(u)] >= u[k]) and (u[(k + 1) % len(u)] >= u[k]):
                    basin_right = k
                    break
                if k == i:
                    basin_right = k
                    break

            theta_saddle = th[i]
            theta_left = th[basin_left]
            theta_right = th[basin_right]
            height_left = float(u[i] - u[basin_left])
            height_right = float(u[i] - u[basin_right])

            saddles.append(Barrier(theta_saddle, theta_left, theta_right, height_left, height_right))
    return saddles


def barriers_to_array(barriers: List[Barrier]) -> np.ndarray:
    """
    Convert a list of Barrier objects to a numeric array:
    columns = [theta_saddle, theta_basin_left, theta_basin_right, height_left, height_right]
    """
    if len(barriers) == 0:
        return np.zeros((0, 5), dtype=float)
    return np.array(
        [[b.theta_saddle, b.theta_basin_left, b.theta_basin_right, b.height_left, b.height_right] for b in barriers],
        dtype=float
    )


# ---------------------------------------------------------------------
# Relative speed proxy from density
# ---------------------------------------------------------------------

def relative_speed_from_density(rho_theta: np.ndarray) -> np.ndarray:
    """
    Relative speed along the ring (dimensionless).
    Args:
        rho_theta: (G,) phase density over θ (non-negative).
    Returns:
        v_rel: (G,) relative speed, proportional to 1 / rho(θ), normalized to mean 1.
    Rationale:
        On a 1D loop at stationarity, continuity implies ρ(θ) v(θ) = const ⇒ v ∝ 1/ρ.
    """
    v = 1.0 / (rho_theta + 1e-12)
    return v / (v.mean() + 1e-12)


# ---------------------------------------------------------------------
# Weighted periodic regression (low-harmonic Fourier)
# ---------------------------------------------------------------------

def fit_periodic(
    theta: np.ndarray,
    y: ArrayLike,
    K: int = 1,
    ridge: float = 0.0,
    weights: Optional[np.ndarray] = None,
) -> Dict[str, float]:
    """
    Fit y ≈ b + a1 cosθ + c1 sinθ (+ a2 cos2θ + c2 sin2θ if K≥2) by weighted least squares.
    Args:
        theta: (n,) angles.
        y: (n,) or (n,1) response values (e.g., gene or pathway scores).
        K: number of harmonics (1 or 2). Default 1.
        ridge: L2 penalty (adds ridge * I to ΦᵀΦ). Default 0.
        weights: optional (n,) non-negative weights for each sample (e.g., cycling probabilities^γ).
    Returns:
        dict with:
          - "beta": coefficients [b, a1, c1, (a2, c2 if K≥2)]
          - "R2": weighted coefficient of determination
          - "amplitude": sqrt(a1^2 + c1^2) (first harmonic)
          - "peak": atan2(-c1, a1) in [0, 2π) (phase of first harmonic maximum)
    Notes:
        Weighted mean used in R² if weights are provided.
    """
    theta = np.asarray(theta)
    y = np.asarray(y).reshape(-1, 1)

    # Design matrix Φ
    Phi_cols = [np.ones_like(theta), np.cos(theta), np.sin(theta)]
    if K >= 2:
        Phi_cols += [np.cos(2 * theta), np.sin(2 * theta)]
    Phi = np.vstack(Phi_cols).T  # n × p

    # Normal equations with optional weights
    if weights is None:
        A = Phi.T @ Phi
        b = Phi.T @ y
        ybar = y.mean()
        denom = ((y - ybar) ** 2).sum() + 1e-12
    else:
        w = weights.reshape(-1, 1)
        A = Phi.T @ (w * Phi)
        b = Phi.T @ (w * y)
        ybar = float((w.T @ y) / (w.sum() + 1e-12))
        denom = (w * (y - ybar) ** 2).sum() + 1e-12

    if ridge > 0:
        A = A + ridge * np.eye(A.shape[0])

    try:
        beta = np.linalg.solve(A, b)  # p × 1
    except np.linalg.LinAlgError:
        beta, *_ = np.linalg.lstsq(A, b, rcond=None)

    yhat = Phi @ beta
    num = (yhat - y) ** 2
    if weights is None:
        ssr = float(num.sum())
    else:
        ssr = float((w * num).sum())

    R2 = 1.0 - ssr / float(denom)
    A1 = float(beta[1])                     # coefficient of cosθ
    C1 = float(beta[2])                     # coefficient of sinθ
    amp = float(np.sqrt(A1**2 + C1**2))
    peak = float(np.arctan2(-C1, A1) % (2.0 * np.pi))

    return {"beta": beta.ravel(), "R2": R2, "amplitude": amp, "peak": peak}


# ---------------------------------------------------------------------
# Circular W1 distance for histograms
# ---------------------------------------------------------------------

def circular_w1_hist(p: np.ndarray, q: np.ndarray) -> float:
    """
    Wasserstein-1 (Earth-Mover’s) distance between two circular histograms.
    Args:
        p, q: (m,) non-negative arrays that sum to 1 (probability mass on m bins).
    Returns:
        Scalar W1 distance (rotation-invariant via “cut” trick).
    Method:
        For each circular rotation (cut), compute line W1 on unwrapped CDFs,
        and take the minimum over cuts.
    Complexity:
        O(m^2) due to trying all m rotations; m=60 is fast in practice.
    """
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    assert p.shape == q.shape
    m = len(p)

    # Defensive normalization
    ps = p / (p.sum() + 1e-12)
    qs = q / (q.sum() + 1e-12)

    def w1_linear(a, b):
        diff_cum = np.cumsum(a - b)
        return np.sum(np.abs(diff_cum)) / m

    best = np.inf
    for k in range(m):
        # Rotate p by k bins relative to q and compute line W1
        best = min(best, w1_linear(np.roll(ps, k), qs))
    return float(best)


# ---------------------------------------------------------------------
# Counterfactual: piecewise angular warp
# ---------------------------------------------------------------------

def piecewise_warp(theta: np.ndarray, a0: float, a1: float, stretch: float = 1.25) -> np.ndarray:
    """
    Stretch a specific arc [a0, a1] by factor 'stretch', compress the complement so total length is 2π.
    Args:
        theta: (n,) original angles in [0, 2π).
        a0, a1: arc endpoints (radians). If a1 < a0, the arc wraps through 2π.
        stretch: factor > 0. Values >1 lengthen [a0, a1]; others compress it.
    Returns:
        (n,) warped angles in [0, 2π).
    Notes:
        - Piecewise linear, monotone mapping on S^1.
        - Useful to emulate interventions (e.g., CDK inhibition lengthening G1).
    """
    L = 2 * np.pi
    a0 = a0 % L
    a1 = a1 % L
    if a1 <= a0:
        a1 += L  # unwrap to ensure a1 > a0

    g1 = a1 - a0  # original arc length
    g2 = L - g1   # remaining arc
    g1p = min(L - 1e-6, g1 * stretch)  # new arc length (cap below 2π)
    g2p = L - g1p
    s1 = g1p / g1  # scale for [a0, a1]
    s2 = g2p / g2  # scale for complement

    th = theta.copy().astype(float)
    t = th.copy()
    # Unwrap angles below a0 forward by 2π so we can treat three linear regions
    t[th < a0] += L

    out = np.empty_like(t)
    # Region inside [a0, a1]
    m = (t >= a0) & (t <= a1)
    out[m] = (t[m] - a0) * s1 + a0
    # Region after a1
    m2 = (t > a1)
    out[m2] = (t[m2] - a1) * s2 + a0 + g1p
    # Region before a0 (wrapped)
    m3 = (t < a0)
    out[m3] = (t[m3] - (a0 - L)) * s2 + a0 - g2p

    return wrap_angle(out)
