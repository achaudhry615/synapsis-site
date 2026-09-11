# Γ-geometry measurement specification

Implemented in `code/nsref/gamma.py`; validated by `code/scripts/validate_pipeline.py`.
Every number in the calibration tables below was produced by that validation run,
not assumed.

## 1. Fields

$\omega=\nabla\times u$ and the velocity gradients are computed spectrally from
the dealiased $\hat u$. The stretching density is evaluated as

$$\Gamma_{\rm net}=\omega\cdot S\omega=\sum_{i,j}\omega_i\omega_j\,\partial_j u_i,$$

which is exact (the antisymmetric part of $\nabla u$ is annihilated by the
symmetric $\omega_i\omega_j$) and avoids materialising the $(3,3,N^3)$ strain
tensor. Record three fields, always: $\Gamma^+=\max(\Gamma_{\rm net},0)$,
$\Gamma^-=\max(-\Gamma_{\rm net},0)$, and the net. **Reporting $\int\Gamma^+$
alone overstates production**; the enstrophy balance (A3) contains the net.

Do **not** normalise by $|\omega|$: $|\omega|^2\,\xi^{\mathsf T}S\xi\equiv\omega\cdot S\omega$,
and the normalised form injects $0/0$ noise at vorticity nulls.

## 2. Thresholding

**Primary (capture-fraction).** $E^{(\varepsilon)}$ = the smallest-volume
super-level set carrying $(1-\varepsilon)$ of $\int\Gamma^+$, for
$\varepsilon\in\{0.1,0.2\}$. Scale-free; no absolute level to guess.

**Secondary (theorem-aligned, absolute).** $\Omega_\lambda=\{|\omega|>\lambda M\}$
for $\lambda\in\{0.1,0.25,0.5\}$, and $E_\lambda=\{\Gamma>\lambda M^3\}$ — **but
not at those values of $\lambda$.** Theorems C3/C4/C5 are stated for these sets,
so they must be reported; the thresholds need care.

**Measured volume fraction of $E_\lambda$ at $64^3$:**

| family | $\lambda=0.5$ | $0.25$ | $0.1$ | $10^{-2}$ | $10^{-3}$ | $10^{-4}$ | capture 0.9 |
|---|---|---|---|---|---|---|---|
| Taylor–Green | 0 | 0 | 0 | 5.4e-2 | 3.0e-1 | 3.9e-1 | 2.2e-1 |
| Kida–Pelz | 0 | 0 | 0 | 1.1e-2 | 3.1e-1 | 3.7e-1 | 2.4e-1 |
| antiparallel tubes | 0 | 0 | 0 | 0 | 0 | 1.8e-4 | 3.6e-2 |
| multiscale random | 0 | 0 | 1.2e-4 | 1.2e-1 | 4.1e-1 | 4.9e-1 | 2.6e-1 |

**$\lambda\in\{0.1,0.25,0.5\}$ gives an empty set for every family.** The reason
is structural: $\Gamma\le\|S\|_\infty M^2=C_{\log}M^3$ with $C_{\log}$ measured at
$0.25$–$0.71$ on these fields, and attaining that bound requires $|\omega|=M$ and
$\alpha$ maximal *at the same point*. The usable absolute range is
$\lambda\sim10^{-3}$–$10^{-2}$, and it is strongly family-dependent: the
antiparallel tubes need $\lambda\le10^{-4}$, because their stretching is
genuinely depleted ($\mathcal R_E\approx0.01$).

**Consequence:** an absolute $\Gamma$ threshold cannot be fixed once for all
families without either emptying the set or, worse, emptying it for exactly the
depleted cases the program is trying to characterise — which would bias the
campaign toward the flows with the *least* depletion. The capture threshold is
well-posed at every instant for every family and is the primary rule; absolute
thresholds are reported as a secondary, per-family-calibrated diagnostic.

Quantile thresholds are recorded but never used for dimension estimates.

## 3. Scale decomposition

$\ell_j=L\,2^{-j}$ for $j=1,\dots,\log_2 N-3$; scales below $8\Delta x$ are
excluded (dealiasing and discretisation pollution). Filter the **velocity**
with a spectral Gaussian, then recompute $\bar\omega_\ell,\bar S_\ell,\bar\Gamma_\ell$
(Germano convention) — do not filter $\Gamma$ directly. Report both if in doubt;
they disagree at small $\ell$, and C3/C4/C5 use the unfiltered sets while D2 uses
the filtered ones.

**$\sum_j\int_{E_j}\Gamma_j$ is not a decomposition of $\int\Gamma$.** $\Gamma$ is
cubic in the fields, so scale-local pieces do not sum to the total. Use $E_j$ for
*geometry per scale* and the unfiltered net for the *budget*.

## 4. The binding resolution constraint (derived, not assumed)

The hypothesis of Theorem C3 lives at $r\sim\ell_\nu=(\nu/M)^{1/2}$. With
$k_{\max}=\pi N/L$,

$$\frac{\eta}{\Delta x}=\frac{k_{\max}\eta}{\pi},\qquad
\frac{\ell_\nu}{\Delta x}=\frac{\ell_\nu}{\eta}\cdot\frac{k_{\max}\eta}{\pi}.$$

**$\ell_\nu/\Delta x$ depends only on $k_{\max}\eta$, not on $N$.** Refining the
grid at fixed $k_{\max}\eta$ buys no measurability whatsoever. Measured on the
$64^3$ Taylor–Green run, $\ell_\nu/\eta\approx0.63$, giving
(`code/scripts/resolution_budget.py`):

| target $k_{\max}\eta$ | $\eta/\Delta x$ | $\ell_\nu/\Delta x$ | verdict |
|---|---|---|---|
| 2 (standard DNS) | 0.64 | **0.40** | $r\sim\ell_\nu$ is sub-grid; $\delta\equiv1$, no information |
| 4 | 1.27 | 0.80 | still sub-grid |
| 10 | 3.18 | 2.01 | marginal: a 5-point segment |
| ~20 | 6.4 | ~4 | measurable |

This was confirmed empirically: on the $64^3$ run at $\nu=0.02$,
$\ell_\nu/\Delta x\in[1.00,1.27]$ and $\delta(r\!\sim\!\ell_\nu)=1.000$ at every
snapshot — exactly as predicted, and carrying no information.

**Consequence for the campaign.** Measuring the geometric hypotheses requires
$k_{\max}\eta\sim10$–$20$, i.e. roughly $5$–$10\times$ the linear resolution of a
"well-resolved" DNS at the same $\nu$ ($10^2$–$10^3\times$ the cost), or a
correspondingly larger $\nu$ at fixed $N$. The geometry measurement and the
blow-up search therefore **compete for the same grid** and cannot be satisfied by
one run. See `campaign_4gpu.md` §2.

## 5. Estimators and their calibration

### 5.1 Sparseness $\delta(r)$ — the theorem-relevant primary (C3)

For each sampled point of the set and each of the 13 canonical cubic-lattice
directions, compute the occupied fraction of the segment of half-length $r$;
take the **minimum over directions** (Lemma C1: thin in *some* direction). Report
the distribution over points — mean, median, p95, max — never a single number,
and **always as a curve in $r$**: $\delta(r)$ is identically 1 below the
transverse feature size and saturates at the global volume fraction above the
object.

Calibration at $64^3$, capture 0.9, $r=16$ cells:

| field | true $D$ | $\delta_{p95}$ |
|---|---|---|
| sheet | 2 | 0.091 |
| tube | 1 | 0.152 |
| uniform random | 3 | 0.606 |

Thin structures separate from space-filling ones by a factor $\sim4$-$7$, and the
ordering sheet < tube < uniform is the expected one (a sheet is thin along one
direction, a tube along two but with a larger transverse extent at this $r$).
**$\delta$ is also the most $N$-robust estimator in the suite** — which is
convenient, since Proposition B4 says it is the only one the theorems use.

### 5.2 Direction regularity $\beta$ (C4)

$\big\langle|\xi(x+h)-\xi(x)|^2\big\rangle_{x\in E}^{1/2}\sim h^\beta$, fitted over
$h\in[4\Delta x,\ell_\nu]$ where that range exists, else over the resolvable
range with the range reported. Theorem C4 needs $\beta\ge1/2$.

Calibration: a smooth (ABC) direction field returns $\beta=0.90$ ($r^2=0.94$)
against a true value of 1. The estimator is biased **low** by ~10%, so a measured
$\beta$ near 0.5 must not be read as a failure of C4 without accounting for it.

### 5.3 Mean-oscillation modulus $\mu_\xi(r)$ (C5) — the sharpest target

$$\mu_\xi(r)=\Big\langle \tfrac{1}{|B_r|}\int_{B_r}\big|\xi-\langle\xi\rangle_{B_r}\big|\Big\rangle_{B_r\cap E\neq\emptyset}.$$

A16 asks whether $\xi$ lies locally in the log-weighted bmo space, i.e. whether
$\mu_\xi(r)\,\log(L/r)$ stays **bounded** as $r\to0$. Report that product, not
$\mu_\xi$ alone. This hypothesis is strictly weaker than Hölder-$1/2$ and is the
one worth targeting: a negative $\beta$ measurement is nearly guaranteed by
intermittency and says little, whereas a failure of the log-modulus bound would
be genuinely informative.

**The log weight must be referenced to the outer scale $L$.** Weighting by
$|\log r|$ instead is a trap that we walked into and caught in calibration:
$|\log r|$ vanishes at $r=1$ and is non-monotone across it, so in a box of size
$2\pi$ the product **inverts the verdict**, reporting a smooth field as growing
and a white-noise field as decaying. Use $\log(L/r)$, which is positive and
monotone for $r<L$.

Calibration (`validate_pipeline.py` §4b), $64^3$:

| direction field | true regularity | $\mu(r)$ log-log slope | $\mu\log(L/r)$ as $r\to0$ |
|---|---|---|---|
| ABC (Beltrami, smooth) | Lipschitz | $+1.13$ | $0.61\to0.15$, **bounded** |
| band-limited random | smooth | $+0.71$ | shrinks, bounded |
| white noise | none | $+0.04$ | $1.39\to3.18$, **unbounded** |

The log-log slope of $\mu(r)$ ($\approx1$ for Lipschitz, $\approx0$ for no
regularity) is the robust discriminator and should be reported alongside the
product.

### 5.4 Dimensions — secondary diagnostics only (Prop. B4)

**Box scales must be fixed in physical units, never in grid cells.** This is the
single largest source of spurious resolution dependence. Measured:

| estimator | blob ($D$=0) | tube ($D$=1) | sheet ($D$=2) | uniform ($D$=3) |
|---|---|---|---|---|
| $D_\infty$ at $64^3$ | 0.53 | 1.17 | 2.09 | 2.74 |
| $D_\infty$ at $128^3$, **matched physical scales** | 0.47 | 1.14 | 2.07 | 2.89 |
| $D_\infty$ at $128^3$, grid-relative scales | 1.36 | 1.62 | 2.31 | 2.74 |

With matched physical scales the maximum drift across a factor-2 refinement is
**0.147**; with grid-relative scales the same quantities drift by up to **0.83**
and the ordering degrades. An apparent "dimension that changes with resolution"
is therefore an artefact until the scaling range is shown to be physically
anchored.

$D_\infty$ is biased low by $0.3$–$0.5$ and is the best-behaved of the
dimension estimators. $D_0$ must be computed on a **thresholded mask**: on a
smooth positive field it returns $3$ identically, since such a field has full
support.

### 5.5 The ratio that decides the balance (Prop. B2)

$$\mathcal R_E=\frac{\int_E\omega\cdot S\omega}{\nu\int_E|\nabla\omega|^2}.$$

Report $\mathcal R_E$, not the global $\mathcal R$: by Corollary B2.1 the global
ratio exceeds 1 **iff** $\Omega'>0$ and is therefore uninformative.

### 5.6 Alignment

Joint pdf of $|\cos\angle(\xi,e_i)|$ for the three strain eigenvectors,
conditioned on the set. The classical expectation is alignment with the
intermediate eigenvector $e_2$; this is the geometric depletion mechanism of A7
made visible, and `NS-CH01` is redefined as *quantifying it on $E$* rather than
as an FTLE statement.

## 6. Real law vs. numerical artefact (preregistered discriminators)

A measured geometric law counts as real only if **all** of these hold:

1. At fixed **physical** $\ell$ in the common resolved range, the estimator is
   $N$-independent within bootstrap CI across a factor-2 refinement.
2. The scaling range's endpoints track $\ell_\nu$ and $L$, **not** $8\Delta x$.
3. $\delta$ at fixed $r/\ell_\nu$ is $N$-independent.
4. $\beta$ and $\mu_\xi|\log r|$ are $N$-independent.

It is an artefact if: the estimator drifts with $N$ at fixed physical $\ell$; the
lower end of the scaling range tracks $\Delta x$; $D_\infty\to0$ because the
maximum sits in a single cell (check against the spectrally-refined maximum); or
$\delta$ changes when $r$ crosses a few $\Delta x$.

## 7. Uncertainty

Spatial CI by bootstrap over the 8 (or 64) sub-cubes; temporal CI by bootstrap
over snapshots. **Never quote a dimension, $\delta$, $\beta$ or $\mu_\xi$ without
a CI.** Bootstrap over sub-cubes underestimates uncertainty when the structure
spans the box; note this where it applies.
