# Model Notes (Steps 1–4 of 10)

This is the authoritative description of our model and how to reproduce it. It covers **Steps 1–4** of the 10-step recipe. Later steps (tests, refinement, OT routing, predictions) build directly on this.

---

## 0) Scope, conventions, and data objects

* **System.** KP mouse model of **lung adenocarcinoma** (LUAD) sampled across tumor evolution (Marjanović *et al.* SuperSeries GSE152607).
* **Primary data.** Single-cell RNA-seq (**scRNA-seq**) is used for all computations in Steps 1–4. Single-cell ATAC-seq (**scATAC-seq**) is optional context later (not used here).
* **Coordinates.** We build a 2-D **cell-cycle plane** with axes given by standardized **S-phase** and **G2/M** module scores. Each cell is a point
  $$y_i = (z_{i,S},, z_{i,M}) \in \mathbb{R}^2.$$
* **Angles.** The **cell-cycle phase** is an angle $\theta \in [0,2\pi)$. We compute angles with $\mathrm{atan2}(y,x)$ and wrap them to $[0,2\pi)$.
* **Notation.** A bar denotes a mean (e.g., $\bar\theta$); a hat denotes an estimate (e.g., $\hat\rho$). “Bootstrap” means resampling cells **with** replacement.

**AnnData fields we create:**
`.obsm["cycle_xy"]=y_i`, `.obs["theta"]=\theta_i`, `.obs["radius"]=r_i`, `.obs["cycling_prob"]=p_i`, `.obs["cycling_call"]∈{0,1}`.

---

## 1) Define the question (and the estimands)

**Biological goal.** Quantify how proliferation changes with tumor stage. Specifically, for each stage we want:

1. **Cycling fraction:** what fraction of cells are proliferating.
2. **Phase occupancy:** the **phase density** $\rho(\theta)$ (where on the cycle cells accumulate).
3. **Checkpoints/barriers:** a **landscape** $U(\theta)$ along the ring and the **barrier heights** near **G1/S** and **G2/M**.
4. **Relative speed:** a phase-dependent speed proxy $\tilde v(\theta)$ around the ring.

We also hypothesize a **High-Plasticity Cell State (HPCS)** functions as a routing **hub**. Its **in-/out-flux** across **adjacent time points** will be quantified later via **optimal transport (OT)**; Steps 1–4 provide the ring and landscape we need first.

**Operational hypotheses.** Later stages have: stronger **loop drive** (more circulation), a **lower G2/M barrier**, a **higher cycling fraction**, and—later—greater **HPCS** in/out-flux.

**Estimands (what we will actually report now).**

* **Per cell:** phase $\theta_i$, radius $r_i$, cycling probability $p_i$, cycling call.
* **Per stage:** phase density $\rho(\theta)$, relative speed $\tilde v(\theta)\propto 1/\rho(\theta)$, and ring landscape $U(\theta)=-\ln P_{\text{ring}}(\theta)$ with barrier heights
  $\text{Barrier}*{\text{G1/S}}$, $\text{Barrier}*{\text{G2/M}}$, and a ring-vs-center stability measure $\text{Barrier}_{\text{Center}}$.

---

## 2) Seek available information (data, features, methods, limits)

**Provenance (SuperSeries GSE152607).**

* **GSE154989** — scRNA-seq **time course** (core dataset for Steps 1–4).
* **GSE154978** — scRNA-seq **TIGIT-sorted** (TIGIT$^+$/$^-$); later used to calibrate/validate HPCS thresholds.
* **GSE154965** — **scATAC-seq**; optional for motif/gene-activity context (not used in Steps 1–4).

**Gene programs.** S-phase set $\mathcal S$ (e.g., *Mcm2–7, Pcna*), G2/M set $\mathcal M$ (e.g., *Top2a, Ccnb1, Aurkb, Mki67*). An **HPCS signature** (includes *Slc4a11*) will be used later for hub calling.

**Method anchors we implement now.**

* **Limit cycle**: the cell cycle is a closed orbit. We embed it explicitly in 2-D (S vs G2/M), fit a **principal circle**, and map cells to **phase** $\theta$ and **radius** $r$.
* **Landscape + barriers**: from **steady density** along the ring we define a **potential** $U=-\ln P_{\text{ring}}$; **local maxima** along the ring are **barriers** (checkpoints).
* **Speed proxy**: a steady-state **continuity** argument on the loop gives $\tilde v(\theta)\propto 1/\rho(\theta)$.

**Limits.** We use **snapshots**, so speeds are **relative** (not minutes). Later OT results depend on the chosen **latent space**; we will report sensitivity then.

---

## 3) Mental model (why this captures the biology)

Think of a **roundabout**. The **angle** $\theta$ is where a cell sits; the **radius** $r$ shows how strongly it expresses cycling programs (outer **ring** = cycling; inner **disk** = quiescent). Where the **phase density** $\rho(\theta)$ is **high**, cells **dwell** (slow segments, “checkpoints”); where $\rho(\theta)$ is **low**, cells pass quickly (fast segments). The **landscape** $U(\theta)$ formalizes this as $U=-\ln P_{\text{ring}}$: **barriers** are ring-wise maxima of $U$. With tumor age we expect more cells on the ring, faster segments to expand, and the **G2/M barrier** to drop.

---

## 4) Modeling concept (algorithms and equations)

### 4.1 Preprocessing and module scores

Start from normalized log-expression $X\in\mathbb{R}^{n\times G}$ (cells $\times$ genes). Standardize each gene $g$ across cells so markers have comparable scale:

$$
Z_{ig} = \frac{X_{ig}-\mu_g}{\sigma_g},
$$

where $\mu_g$ and $\sigma_g$ are mean and standard deviation of gene $g$ across cells (robust median/MAD also acceptable).

Compute **module scores** by equal-weight averaging over the S and G2/M sets:

$$
z_{i,S}=\frac{1}{|\mathcal S|}\sum_{g\in\mathcal S} Z_{ig},
\qquad
z_{i,M}=\frac{1}{|\mathcal M|}\sum_{g\in\mathcal M} Z_{ig}.
$$

This yields each cell’s coordinate in the **cell-cycle plane**:

$$
y_i = (z_{i,S},, z_{i,M}) \in \mathbb{R}^2.
$$

**If the ring is weak:** compute the top **PCA**/**ICA** components on $\mathcal S\cup\mathcal M$ and use those two axes instead; downstream steps are unchanged.

---

### 4.2 Principal circle (map cells to phase and radius)

Fit a **least-squares circle** to the points ${y_i}$:

$$
(c^*, r^*) ;=; \arg\min_{c\in\mathbb{R}^2,; r>0}
\sum_i \big(|y_i - c| - r\big)^2,
$$

where $c^*$ is the center and $r^*$ the radius. For each cell:

* **Phase (angle)**
  $$\theta_i ;=; \mathrm{atan2}\big((y_i-c^*)_2,\ (y_i-c^*)_1\big)\ \bmod 2\pi,$$
* **Radius**
  $$r_i ;=; |y_i - c^*|.$$

**Orientation/anchoring.** We want **S near $0$**, **G2/M later**:

1. Identify S-high (top decile of $z_{i,S}$) and M-high (top decile of $z_{i,M}$) cells.
2. Compute circular means $\bar\theta_S$ and $\bar\theta_M$.
3. Rotate all phases by $\delta=-\bar\theta_S$ (so S is near $0$).
4. If $(\bar\theta_M-\bar\theta_S)\bmod 2\pi<0$, flip orientation via $\theta_i\leftarrow (2\pi-\theta_i)\bmod 2\pi$.

**Diagnostics.** S→G2/M phase gap $\approx \pi/2$ radians (acceptable if $\ge 0.8$ rad). Circle residuals $\big||y_i-c^*|-r^*\big|$ should be small (report median ± MAD).

---

### 4.3 Cycling vs quiescent (ring–disk mixture on radius)

We classify cells with a **Gaussian Mixture Model (GMM)** on radii ${r_i}$:

$$
r_i \sim
w,\mathcal N(\mu_{\text{ring}},\sigma_{\text{ring}}^2)
;+;
(1-w),\mathcal N(\mu_{\text{disk}},\sigma_{\text{disk}}^2),
\quad 0<w<1.
$$

Fit by **Expectation–Maximization** (init $k$-means on ${r_i}$; stop when $\Delta\log\mathcal L<10^{-6}$ or after 500 iters). The **posterior** probability of being on the ring is

$$
p_i ;=;
\frac{w,\mathcal N(r_i\mid \mu_{\text{ring}},\sigma_{\text{ring}}^2)}
{w,\mathcal N(r_i\mid \mu_{\text{ring}},\sigma_{\text{ring}}^2)
+ (1-w),\mathcal N(r_i\mid \mu_{\text{disk}},\sigma_{\text{disk}}^2)}.
$$

Keep $p_i$ as a **soft weight**; define a **binary call** by $p_i\ge 0.5$.

**Optional calibration with markers.** Regress $p_i$ on $r_i$ and a standardized proliferation marker score $m_i$ (e.g., MKI67/TOP2A) via ridge-penalized logistic regression:

$$
\mathrm{logit}(p_i) = \beta_0 + \beta_r, r_i + \beta_m, m_i,
$$

tuning the penalty by cross-validation. Replace $p_i$ by $\hat p_i$ if this improves correlation with markers.

**Diagnostics.** $p_i$ should correlate positively with MKI67/TOP2A/AURKB. Cycling fractions by stage should be plausible (export binomial CIs).

---

### 4.4 Phase density, ring landscape, and barrier finding

#### 4.4.1 Phase density on the circle (circular KDE)

Among **cycling** cells (use weights $w_i=p_i^\gamma$, default $\gamma=1$), estimate the **phase density** with a **von Mises** kernel:

$$
\hat\rho(\theta) = \frac{1}{\sum_i w_i}\sum_i w_i,
\underbrace{\frac{e^{\kappa\cos(\theta-\theta_i)}}{2\pi,I_0(\kappa)}}*{K*\kappa(\theta-\theta_i)},
$$

where $\kappa>0$ is the **concentration** (inverse bandwidth) and $I_0$ is the modified Bessel function of order 0.

**Choosing $\kappa$.** Prefer leave-one-out log-likelihood maximization over a grid; alternatively map the resultant length to $\kappa$ via the standard circular-statistics approximation.

#### 4.4.2 Joint steady density on the cylinder and ring restriction

Estimate a **joint steady density** over phase and radius with a **product kernel**:

$$
\hat P_{ss}(\theta,r)=\frac{1}{\sum_i w_i}\sum_i w_i,
K_\kappa(\theta-\theta_i)\cdot \mathcal N(r\mid r_i, h_r^2),
$$

where $h_r$ is the Gaussian **radial bandwidth** (Silverman’s rule or CV). To restrict to the ring, weight by an annulus centered at $r^*$ (use the ring component’s $\sigma_{\text{ring}}$ from the GMM):

$$
\hat P_{\text{ring}}(\theta)
= \int \hat P_{ss}(\theta,r),\underbrace{\mathcal N(r\mid r^*, \tau^2)}*{\psi(r)},dr,
\quad \tau:=\sigma*{\text{ring}}.
$$

Define the **ring potential** as

$$
U(\theta) = -\ln \hat P_{\text{ring}}(\theta).
$$

Lower $U$ means more probable phases; **local maxima** along the ring are **barriers**.

#### 4.4.3 Barrier heights (G1/S, G2/M, center)

1. Smooth $U(\theta)$ with a **periodic cubic spline** on a uniform grid in $\theta$.
2. Find zeros of $U'(\theta)$ and classify by $U''(\theta)$: minima are **basins**, maxima are **saddles** (barriers along the ring).
3. Define barrier heights as saddle minus neighboring basin minimum; e.g.,

* **G1/S barrier**
  $$\text{Barrier}*{\text{G1/S}} = U(\theta*{\text{saddle,G1/S}}) - \min_{\theta\in\text{G1 arc}} U(\theta).$$

* **G2/M barrier**
  $$\text{Barrier}*{\text{G2/M}} = U(\theta*{\text{saddle,G2/M}}) - \min_{\theta\in\text{S/G2 arc}} U(\theta).$$

* **Center barrier** (ring vs center stability)
  $$\text{Barrier}*{\text{Center}} = U*{\text{center}} - \min_\theta U(\theta),\quad
  U_{\text{center}} := -\ln \hat P_{ss}(\theta^*, r{=}0)\ \text{(averaged over $\theta^*$)}.$$

4. **Uncertainty.** Bootstrap cells within stage; recompute $\hat\rho$, $U$, and barriers; report percentile CIs.

#### 4.4.4 Relative speed from continuity (why $\tilde v \propto 1/\rho$)

On a 1-D loop parameterized by $\theta$, stationarity implies the **continuity equation**:

$$
\frac{d}{d\theta}\big(\rho(\theta),v(\theta)\big) = 0
\ \Rightarrow\
\rho(\theta),v(\theta)=\text{const}
\ \Rightarrow
\tilde v(\theta)\propto \frac{1}{\rho(\theta)}.
$$

We report $\tilde v(\theta)$ normalized to mean 1 over $\theta$ (a **relative** speed, not minutes).

---

### 4.5 Periodic gene/pathway fits (optional but recommended)

For gene (or pathway score) $x_{ig}$, fit a **weighted Fourier** model:

$$
x_{ig}\approx b_g + a_{1g}\cos\theta_i + a_{2g}\sin\theta_i
\quad (+\ a_{3g}\cos 2\theta_i + a_{4g}\sin 2\theta_i),
$$

using **ridge-penalized** least squares with weights $w_i=p_i^\gamma$ (default $\gamma=1$). Summaries:

$$
A_g=\sqrt{a_{1g}^2+a_{2g}^2},
\qquad
\phi_g=\mathrm{atan2}(-a_{2g},,a_{1g}).
$$

Default to 1 harmonic; allow 2 if BIC improves by ≥ 10.

---

### 4.6 Circular distance between stages (Wasserstein-1)

To compare phase distributions for stages $A$ and $B$, compute a **circular Earth-Mover’s distance**:

1. Bin $\theta$ into $K$ equal bins to get histograms $p^A, p^B$.
2. For each circular **cut** $c\in{1,\dots,K}$, unwrap at $c$, form cumulative sums $F_c^A,F_c^B$.
3. Compute
   $$W_1^{(c)}=\frac{1}{K}\sum_{k=1}^{K}\big|F_c^A(k)-F_c^B(k)\big|,$$
   and report $W_1=\min_c W_1^{(c)}$ (rotation-invariant).

---

## Practical checks and diagnostics

* **Geometry.** S→G2/M gap $\approx \pi/2$ (≥ 0.8 rad ok). Circle residuals small (report median ± MAD).
* **Markers.** $p_i$ correlates positively with MKI67/TOP2A; cycling fractions by stage are sensible.
* **Landscape.** $U(\theta)$ shows two ring saddles (≈ G1/S and G2/M). Fast arcs (low $\rho$, high $\tilde v$) align with mitotic marker peaks.
* **Between-stage.** Circular $W_1$ distinguishes early vs late.

Always wrap phases to $[0,2\pi)$; compute circular means with vector sums; floor densities before logs ($\hat P\leftarrow \max(\hat P,10^{-9})$); fix seeds for GMM and bootstraps; use **periodic** cubic splines (enforce $U(0)=U(2\pi)$ and matching derivatives).

---

## Summary of key concepts

* The **ring** captures the limit-cycle geometry; the **landscape** $U$ and **barriers** convert intuitive dwell/flow into measurable quantities.
* The continuity-based **speed proxy** ties density to dynamics without time-lapse or RNA velocity.
* These same phase coordinates and ring restriction support later **OT** (for HPCS routing) and **counterfactuals** (e.g., CDK-inhibition warp of the G1 arc).

---

