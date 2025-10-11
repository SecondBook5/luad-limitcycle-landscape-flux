# Model Notes (Steps 1–4 of 10; discovery-based)

This document specifies **what we will measure and how**, without presupposing any hypothesis. It implements **Steps 1–4** of the 10-step recipe and sets up later analyses (e.g., OT-based routing). All mathematical symbols are defined where they first appear.

---

## 0) Scope, conventions, and data objects

* **System.** KP mouse model of lung adenocarcinoma (LUAD) sampled across tumor evolution (Marjanović SuperSeries).
* **Primary data.** Single-cell RNA-seq (**scRNA-seq**) for all computations in Steps 1–4. Single-cell ATAC-seq (**scATAC-seq**) may be referenced later but is not used here.
* **Coordinate system.** We form a 2-D **cell-cycle plane** per cell from standardized S-phase and G2/M module scores (defined in §4.1). A cell’s coordinate is
  $$
  y_i \equiv (z_{i,S},, z_{i,M}) \in \mathbb{R}^2,
  $$
  where $i$ indexes cells, $z_{i,S}$ and $z_{i,M}$ are S and G2/M score z-averages.
* **Angle (phase).** The **cell-cycle phase** is the polar angle $\theta\in[0,2\pi)$ computed by $\mathrm{atan2}(\text{y},\text{x})$ and wrapped to $[0,2\pi)$.
* **Radius.** The **cycle radius** $r$ is the Euclidean distance from the circle center (defined in §4.2); it quantifies cycle program strength.
* **Notation.** Overbars denote means (e.g., $\bar\theta$). Hats denote estimates (e.g., $\hat\rho$). Weights $w_i\ge 0$ are per-cell and sum to $\sum_i w_i$. “Bootstrap” means resampling cells with replacement.

**AnnData fields we write:**
`.obsm["cycle_xy"]=y_i`, `.obs["theta"]=\theta_i`, `.obs["radius"]=r_i`, `.obs["cycling_prob"]=p_i`, `.obs["cycling_call"]∈{0,1}`.

---

## 1) Problem statement (exploratory, no hypothesis)

We **describe** cell-cycle organization and circulation across tumor stages.

**Stage-wise measurements we will report** (discovery-oriented):

1. **Cycling prevalence.** Fraction of cells on the cycle (“ring”) versus quiescent (“disk”).
2. **Phase occupancy.** The **phase density** $\rho(\theta)$ around the ring, indicating dwell structure along the cycle.
3. **Landscape and barriers.** A **ring potential** $U(\theta)$ derived from steady density; **barrier heights** near checkpoint-like regions (G1/S and G2/M) and a **center barrier** (ring vs center stability).
4. **Flux components (descriptive).**

   * **Within-stage speed shape:** relative speed $\tilde v(\theta)$ inferred from $\rho(\theta)$ (no absolute time).
   * **Inter-stage transport on the ring:** signed mass transport between consecutive stages using **optimal transport (OT)** on the circle (§4.6).
   * **Ring↔center leakage:** change in cycling mass between stages (§4.7).
5. **Rotation-invariant contrasts.** Circular Wasserstein-1 distance between stage-wise phase distributions (§4.9).

All quantities are **descriptive estimands**, not tests of a prior hypothesis.

---

## 2) Available information (datasets, features, methods, limits)

* **Datasets (GEO SuperSeries).**

  * **GSE154989**: scRNA-seq **time course** (core analysis here).
  * **GSE154978**: scRNA-seq **TIGIT-sorted** (later: to calibrate/validate a high-plasticity cell-state signature).
  * **GSE154965**: **scATAC-seq** (optional context later).
* **Gene programs (for the cycle plane).** S-phase set $\mathcal S$ (e.g., *Mcm2–7, Pcna*); G2/M set $\mathcal M$ (e.g., *Top2a, Ccnb1, Aurkb, Mki67*).
* **Method anchors used now.**

  * **Limit-cycle geometry** in a 2-D plane (S vs G2/M).
  * **Landscape on the ring**: potential $U(\theta)=-\ln P_{\text{ring}}(\theta)$ and **barriers** along $\theta$.
  * **Flux framing** compatible with Li–Wang: within-stage *relative* speed from density; between stages, **OT-based transport** around the ring (no parametric kinetics).
* **Limits.**

  * Snapshots only: speeds are **relative**; no absolute minutes.
  * OT transport depends on latent space and cost; we will disclose parameters and run sensitivity checks later.

---

## 3) Mental model (intuitive picture)

Cells form a **ring** (cycling) around an inner **disk** (quiescent). The **angle** $\theta$ indicates position along the cell cycle; the **radius** $r$ indicates strength of cycling programs. High **phase density** $\rho(\theta)$ marks **slow/dwell** segments (checkpoint-like); low $\rho(\theta)$ marks **fast** segments. A **landscape** $U(\theta) = -\ln P_{\text{ring}}(\theta)$ summarizes dwell structure; **saddles** on the ring correspond to **barriers**. Flux is handled two ways: (i) within a stage, relative speed from $\rho(\theta)$; (ii) between stages, **empirical transport** via OT on the circle.

---

## 4) Modeling concept (algorithms and equations)

### 4.1 Preprocessing and module scores

Let $X\in\mathbb{R}^{n\times G}$ be normalized log-expression (cells $\times$ genes).

**Per-gene standardization (z-scores).**
$$
Z_{ig} ;=; \frac{X_{ig}-\mu_g}{\sigma_g},
$$
with $\mu_g,\sigma_g$ the mean and standard deviation of gene $g$ across cells (robust median/MAD also acceptable).

**Module scores (equal-weight means).**
$$
z_{i,S} ;=; \frac{1}{|\mathcal S|}\sum_{g\in\mathcal S} Z_{ig},
\qquad
z_{i,M} ;=; \frac{1}{|\mathcal M|}\sum_{g\in\mathcal M} Z_{ig}.
$$

Define the **cell-cycle plane** coordinate:
$$
y_i \equiv (z_{i,S},, z_{i,M}) \in \mathbb{R}^2.
$$

**If the ring is diffuse:** compute the top **PCA** or **ICA** components on $\mathcal S\cup\mathcal M$ and use those two axes instead; downstream steps remain identical.

---

### 4.2 Principal circle (phase and radius)

Fit a **least-squares circle** to the points ${y_i}$:
$$
(c^{*},,r^{*}) ;=; \arg\min_{c\in\mathbb{R}^{2},, r>0}
\sum_{i}\big(,|y_i - c| - r,\big)^{2},
$$
where $c^{*}$ is the center and $r^{*}$ the radius.

**Per-cell phase and radius**
$$
\theta_i ;=; \mathrm{atan2}!\big((y_i-c^{*})_{2},,(y_i-c^{*})_{1}\big)\ \bmod 2\pi,
\qquad
r_i ;=; |y_i - c^{*}|.
$$

**Orientation and anchoring.** Compute circular means for S-high (top decile in $z_{i,S}$) and G2/M-high (top decile in $z_{i,M}$) cells: $\bar\theta_S,\bar\theta_M$. Rotate by $\delta=-\bar\theta_S$ so S is near $0$: $\theta_i\leftarrow(\theta_i+\delta)\bmod 2\pi$. If $(\bar\theta_M-\bar\theta_S)\bmod 2\pi<0$, flip orientation via $\theta_i\leftarrow(2\pi-\theta_i)\bmod 2\pi$.

**Diagnostics.**

* **S→G2/M gap** $(\bar\theta_M-\bar\theta_S)\bmod 2\pi\approx \pi/2$ (acceptable if $\ge 0.8$ rad).
* **Circle residuals** $\big||y_i-c^{*}|-r^{*}\big|$ (report median ± MAD).

---

### 4.3 Cycling vs. quiescent (ring–disk mixture on radius)

Model radii ${r_i}$ with a **two-component Gaussian mixture** (ring vs. disk):
$$
r_i \sim w,\mathcal{N}(\mu_{\text{ring}},\sigma_{\text{ring}}^2);+;(1-w),\mathcal{N}(\mu_{\text{disk}},\sigma_{\text{disk}}^2), \quad 0<w<1.
$$

Fit by **Expectation–Maximization** (EM), initialized with $k$-means on ${r_i}$; stop at $\Delta\log\mathcal L<10^{-6}$ or 500 iterations.

**Posterior (soft weight) and call.**
$$
p_i ;\equiv; \Pr(\text{ring}\mid r_i)
;=; \frac{w,\mathcal{N}(r_i\mid\mu_{\text{ring}},\sigma_{\text{ring}}^2)}{w,\mathcal{N}(r_i\mid\mu_{\text{ring}},\sigma_{\text{ring}}^2);+;(1-w),\mathcal{N}(r_i\mid\mu_{\text{disk}},\sigma_{\text{disk}}^2)}.
$$
Use $p_i$ as a **weight**; define a binary **cycling call** by $p_i\ge 0.5$.

**Optional calibration with proliferation markers.** Let $m_i$ be a standardized marker score (e.g., MKI67/TOP2A). Fit **ridge-penalized** logistic regression
$$
\mathrm{logit},p_i ;=; \beta_0 + \beta_r,r_i + \beta_m,m_i,
$$
with penalty chosen by cross-validation; replace $p_i$ by the calibrated $\hat p_i$ if it improves marker concordance.

**Diagnostics.** Positive correlation of $p_i$ with MKI67/TOP2A/AURKB; stage-wise cycling fractions with binomial CIs.

---

### 4.4 Phase density, ring potential, and barriers

We estimate a **steady density** $P_{ss}(\theta,r)$ on the **cylinder** (phase $\theta$ and radius $r$), then **restrict** to the ring to obtain $U(\theta)$.

#### 4.4.1 Circular phase density (von Mises KDE)

For cycling-weighted cells (weights $w_i=p_i^\gamma$, default $\gamma=1$), use a **von Mises** kernel (circular analogue of a Gaussian). The kernel at angular difference $\Delta$ is
$$
K_\kappa(\Delta) ;=; \frac{e^{,\kappa\cos\Delta}}{2\pi,I_0(\kappa)},
$$
where $\kappa>0$ is the **concentration** (inverse bandwidth) and $I_0$ is the **modified Bessel function of the first kind** (order 0).

The **phase density** is
$$
\hat\rho(\theta) ;=; \frac{1}{\sum_i w_i}\sum_i w_i,K_\kappa(\theta-\theta_i).
$$

**Choosing $\kappa$.** Prefer leave-one-out log-likelihood maximization; a standard resultant-length rule is available as a fallback.

#### 4.4.2 Cylinder density and ring restriction

Estimate the **joint steady density** with a **product kernel**:
$$
\hat P_{ss}(\theta,r) ;=; \frac{1}{\sum_i w_i}\sum_i w_i; K_\kappa(\theta-\theta_i),\cdot, \mathcal{N}!\big(r\mid r_i,,h_r^2\big),
$$
where $h_r$ is the **radial bandwidth** (Silverman’s rule or CV).

Restrict to the ring with an **annulus weight** centered at $r^{*}$:
$$
\hat P_{\text{ring}}(\theta) ;=; \int \hat P_{ss}(\theta,r);\psi(r),dr,
\qquad
\psi(r);\propto; \mathcal{N}!\big(r\mid r^{*},,\tau^2\big),
$$
taking $\tau=\sigma_{\text{ring}}$ from the mixture’s ring component.

Define the **ring potential**
$$
U(\theta) ;=; -\ln \hat P_{\text{ring}}(\theta),
$$
after flooring $\hat P_{\text{ring}}$ at a small $\epsilon$ (e.g., $10^{-9}$) to avoid numerical issues.

#### 4.4.3 Barrier locations and heights

1. Smooth $U(\theta)$ with a **periodic cubic spline** on a uniform grid in $\theta$.
2. Find zeros of $U'(\theta)$; classify with $U''(\theta)$: minima are **basins**, maxima along the ring are **saddles** (barriers).
3. Report **barrier heights** as the saddle minus the neighboring basin minimum:
   $$
   \text{Barrier}*{\text{G1/S}} ;=; U(\theta*{\text{saddle,G1/S}}) - \min_{\theta\in\text{G1 arc}} U(\theta),
   $$
   $$
   \text{Barrier}*{\text{G2/M}} ;=; U(\theta*{\text{saddle,G2/M}}) - \min_{\theta\in\text{S/G2 arc}} U(\theta).
   $$
   For **ring vs center**,
   $$
   \text{Barrier}*{\text{Center}} ;=; U*{\text{center}} - \min_\theta U(\theta),
   \qquad
   U_{\text{center}} ;\equiv; -\ln \hat P_{ss}(\theta^{*},,0)\ \text{(averaged over $\theta^{*}$)}.
   $$
4. **Uncertainty.** Bootstrap cells within stage; recompute $\hat\rho$, $U$, and barrier heights; report percentile CIs.

---

### 4.5 Within-stage relative speed (continuity-based proxy)

On a 1-D ring at stationarity, the **probability flux** is $J(\theta)=\rho(\theta),v(\theta)$, and the continuity equation implies $\frac{dJ}{d\theta}=0$; thus $J$ is constant (no absolute time scale). Consequently,
$$
\tilde v(\theta) ;\propto; \frac{1}{\hat\rho(\theta)},
\qquad
\langle \tilde v \rangle_\theta ;=; 1,
$$
which serves as an **intra-stage relative speed** profile. (We do **not** interpret its magnitude across stages absolutely; we use the *shape* within a stage and compare normalized profiles across stages.)

---

### 4.6 Inter-stage transport on the ring (OT-based flux)

To capture **empirical circulation** between consecutive time points without assuming kinetics, we couple stage $t$ to $t{+}1$ via **optimal transport** on the circle.

**Setup.**

* Discretize the circle into $K$ bins with centers ${\theta_k}$.
* Let $p^{(t)}\in\Delta^{K-1}$ and $p^{(t+1)}\in\Delta^{K-1}$ be the **phase histograms** at stages $t$ and $t{+}1$ (weighted by $w_i=p_i^\gamma$).
* Define a **circular cost** $C_{kl}=\min!\big(|\theta_k-\theta_l|,,2\pi-|\theta_k-\theta_l|\big)$.

**Coupling.** Solve **entropically regularized** (Sinkhorn) transport to obtain a plan $\Pi^{(t\to t+1)}\in\mathbb{R}_{\ge 0}^{K\times K}$ with row sums $p^{(t)}$ and column sums $p^{(t+1)}$ (or **unbalanced OT** if total mass differs).

**Net directional transport (clockwise vs counter-clockwise).**
Map each pair $(k,l)$ to a signed angular displacement $\Delta\theta_{kl}\in(-\pi,\pi]$ (positive = clockwise). Define the **signed transport**
$$
\Phi^{(t)} ;=; \sum_{k,l}\Pi^{(t\to t+1)}*{kl};\Delta\theta*{kl}.
$$
Report:

* **Normalized directionality** $\Phi^{(t)}/\sum_{k,l}\Pi^{(t\to t+1)}_{kl}$ (in radians per stage; sign indicates net direction).
* **Local transport field** on the ring by aggregating $\Pi$ along the $\theta$ bins to visualize where gains/losses occur.

This **inter-stage flux component** is purely empirical and discovery-oriented: it reveals whether mass preferentially advances around the ring, stalls, or even shows reverse flow between specific stages.

---

### 4.7 Ring↔center leakage (cycling mass budget)

Let $M_{\text{ring}}^{(t)}=\sum_{i\in\text{stage }t} p_i$ and $M_{\text{tot}}^{(t)}$ the total (weighted) cell mass. The **cycling fraction** at stage $t$ is
$$
f_{\text{ring}}^{(t)} ;=; \frac{M_{\text{ring}}^{(t)}}{M_{\text{tot}}^{(t)}}.
$$
Between $t$ and $t{+}1$, define the **leakage (net)** as
$$
\Delta f_{\text{ring}}^{(t)} ;=; f_{\text{ring}}^{(t+1)} - f_{\text{ring}}^{(t)}.
$$
A positive $\Delta f_{\text{ring}}^{(t)}$ indicates **influx** from the center to the ring (or net growth of cycling mass); negative indicates **outflux** to the center. Report with bootstrap CIs.

---

### 4.8 Periodic genes / pathways (optional but recommended)

For any gene (or pathway score) $x_{ig}$, fit a **weighted Fourier** model against phase:
$$
x_{ig};\approx; b_g ;+; a_{1g}\cos\theta_i ;+; a_{2g}\sin\theta_i
\quad \big(+\ a_{3g}\cos 2\theta_i ;+; a_{4g}\sin 2\theta_i\big),
$$
by **ridge-penalized** least squares with weights $w_i=p_i^\gamma$ (default $\gamma=1$). Summaries:
$$
A_g ;=; \sqrt{a_{1g}^2+a_{2g}^2},
\qquad
\phi_g ;=; \mathrm{atan2}(-a_{2g},,a_{1g}).
$$
Prefer one harmonic; allow two if BIC improves by $\ge 10$.

---

### 4.9 Circular distance between stages (Wasserstein-1)

To compare phase distributions for stages $A$ and $B$, compute a **rotation-invariant** circular Earth-Mover’s distance:

1. Bin $\theta$ into $K$ bins to obtain histograms $p^A, p^B$.
2. For each circular cut index $c$, unwrap at $c$, form cumulative distributions $F_c^A,F_c^B$.
3. Compute
   $$
   W_1^{(c)} ;=; \frac{1}{K}\sum_{k=1}^{K} \big|,F_c^A(k)-F_c^B(k),\big|,
   $$
   and report $W_1=\min_c W_1^{(c)}$.

---

## Data-quality checks (non-inferential)

* **Geometry.** S→G2/M gap near $\pi/2$ (≥ 0.8 rad acceptable). Circle residuals small (median ± MAD).
* **Markers.** $p_i$ positively correlates with MKI67/TOP2A; stage-wise cycling fractions look plausible.
* **Landscape.** $U(\theta)$ exhibits two ring saddles consistent with G1/S and G2/M; fast arcs (low $\hat\rho$, high $\tilde v$) align with mitotic marker peaks.
* **Flux components.**

  * **Within-stage:** report $\tilde v(\theta)$ (normalized mean 1); interpret shape only.
  * **Inter-stage:** report signed OT transport $\Phi^{(t)}$ and local gain/loss maps around $\theta$.
  * **Leakage:** report $\Delta f_{\text{ring}}^{(t)}$ with bootstrap CIs.
* **Between-stage distance.** Circular $W_1$ as a rotation-invariant descriptor.

**Implementation notes.** Wrap $\theta$ to $[0,2\pi)$. Use vector sums for circular means. Floor densities before logs. Use **periodic** cubic splines ($U(0)=U(2\pi)$ and matched derivatives). Fix random seeds for GMM, bootstraps, and OT.

---

## Inputs, outputs, and notebook mapping

**Inputs.**

* `data/interim/scrna/timecourse.h5ad` (ingested in `01_ingest_build_anndata.ipynb`).
* Marker lists in `configs/markers/` (S, G2/M).
* (Optional) scPrisma-enhanced cycle coordinates.

**Outputs (from Steps 2–4 notebooks).**

* Per cell: $(\theta_i,\ r_i,\ p_i)$ and the binary cycling call.
* Landscape: $U(\theta)$; $\text{Barrier}*{\text{G1/S}}$, $\text{Barrier}*{\text{G2/M}}$, $\text{Barrier}_{\text{Center}}$.
* Densities/speeds: $\hat\rho(\theta)$, $\tilde v(\theta)$.
* Flux components: inter-stage OT transport $\Pi^{(t\to t+1)}$, signed directionality $\Phi^{(t)}$, leakage $\Delta f_{\text{ring}}^{(t)}$.
* Comparisons: circular $W_1$ between stages.
* (Optional) Periodic features: $(A_g,\ \phi_g)$ per gene/pathway.

**Notebook mapping.**

* `02_modules_circle_phase.ipynb` → §4.1–4.2
* `03_mixture_cycling.ipynb` → §4.3
* `04_landscape_flux_liwang.ipynb` → §4.4–4.7, §4.9
* `04b_periodic_genes.ipynb` (optional) → §4.8

**Default parameters** (see `configs/params.yaml`).

* Phase bins: $K=60$
* Weight exponent: $\gamma=1.0$
* Circular KDE: $\kappa$ by cross-validated likelihood
* Radial bandwidth: $h_r$ by Silverman’s rule (or CV)
* GMM EM tolerance: $\Delta\log\mathcal L<10^{-6}$ or 500 iters
* OT: Sinkhorn $\varepsilon$ default $0.05\times$ median cost; unbalanced $\tau$ if needed (document values)
* Fourier ridge: CV; max harmonics $K=2$

---

## Why this is discovery-aligned

* We **do not** assert hypotheses; we **expose** structure: cycle prevalence, dwell/barriers, within-stage speed shape, empirical transport between stages, and ring↔center leakage.
* The **flux components** are descriptive: within-stage speed from density (no time scale), inter-stage directionality from OT, and a simple mass budget for cycling vs quiescent cells.

---

*End of model.md.*
