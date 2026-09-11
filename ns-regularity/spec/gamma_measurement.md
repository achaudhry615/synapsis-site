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

## 4. The measurement scale — **revised after citation verification**

An earlier version of this section concluded that the Theorem C3 hypothesis is
effectively unmeasurable, on the grounds that it lives at $r\sim\ell_\nu$ and
$\ell_\nu/\Delta x\approx0.4$ at the standard DNS criterion $k_{\max}\eta\ge2$.
**That conclusion was wrong, and wrong because the scale was wrong.**

A13 states the criterion at a scale comparable to the **radius of spatial
analyticity** $R_{\rm an}$, not at $\ell_\nu$. The two are of the same order in
the sense that both are viscous scales, but they differ by a substantial
constant, and that constant decides measurability.

### 4.1 Measuring $R_{\rm an}$

A function analytic in a strip of width $R$ has Fourier coefficients decaying
like $e^{-R|k|}$, so $E(k)\sim e^{-2Rk}$. Fit $\log E(k)=a-2Rk$ over a window
above the energy-containing range and below the $2/3$ dealiasing cutoff
(`diagnostics.analyticity_radius`). **Report $r^2$ and refuse the estimate when
the fit is poor** — for a band-limited initial condition (Taylor–Green,
Kida–Pelz at $t=0$) there is no exponential range at all and the fit returns
$r^2\approx0.04$ with a nonsensical negative $R$. The estimator correctly flags
this rather than returning a number.

### 4.2 What was measured

At $64^3$, $\nu=0.02$, once the flow has developed ($t\gtrsim1$, $r^2\ge0.96$):

| family | $R_{\rm an}/\Delta x$ | $\ell_\nu/\Delta x$ | $R_{\rm an}/\ell_\nu$ |
|---|---|---|---|
| Taylor–Green | 4.8–8.0 | 1.17–1.24 | 3.8–6.9 |
| antiparallel tubes | 10.2–12.0 | 1.44–2.20 | 5.5–7.1 |
| Kida–Pelz | 2.8 | 0.65 | 4.3 |

**$R_{\rm an}\approx4$–$7\,\ell_\nu$, and $R_{\rm an}/\Delta x\approx3$–$12$.** The
theorem's scale is therefore comfortably resolved at $64^3$ while $\ell_\nu$ is
not. Evaluating sparseness at $\ell_\nu$ — as the earlier draft did — used a
scale $4$–$7\times$ too small, which is exactly why $\delta$ came back pinned at
$1.000$ (a segment of one or two cells is fully occupied by construction, for any
set).

### 4.3 Consequence

With both corrections applied — the six component super-level sets (§2) and the
analyticity-radius scale — the C3 hypothesis becomes **measurable at standard
resolution**. First informative measurement, antiparallel tubes at $64^3$:

$$\delta_{\rm worst,\,p95}\big(r\sim R_{\rm an}\big)\approx0.28\text{–}0.32,
\qquad R_{\rm an}/\Delta x\approx11,\qquad r^2=1.00 .$$

Whether $0.3$ clears the theorem's threshold is **unknown**, because $\delta_0$ is
still unverified (`prereg/FREEZE_BLOCKERS.md` item 2). So H2 remains untestable —
but the quantity it is about now exists, which it did not before.

**The earlier claim that the geometry measurement and the blow-up search compete
for the same grid is therefore withdrawn** in the form it was stated. They still
prefer different viscosities, but not by the order of magnitude previously
asserted, and no special over-resolution is needed to see $\delta$ at the
theorem's scale. $k_{\max}\eta\ge2$ remains required for the run to be trustworthy
at all.

### 4.4 What still binds

Two honest caveats:

1. $R_{\rm an}$ from spectral decay is a **global** quantity; the criterion wants a
   uniform lower bound on a local analyticity radius. A global fit can only
   overestimate the local worst case, so $\delta(r\sim R_{\rm an}^{\rm global})$
   is evaluated at a scale at least as large as the theorem wants — the
   conservative direction for a *sparseness* (thinness) claim is the larger $r$,
   so this is not obviously safe and should be stated as a caveat, not waved away.
2. The fit window is a choice. Report it, and report $r^2$ with every $R_{\rm an}$.

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

**Under-resolution drives $\beta\to0$, the white-noise value.** Measured on the
broadband IC at $32^3$ with $\nu=0.02$ ($k_{\max}\eta=0.77$, a factor 6.7 short of
resolved): $\beta=0.06$–$0.09$. Taken at face value that is a catastrophic
violation of the Theorem C4 hypothesis $\beta\ge1/2$; it is entirely a numerical
artefact, because an under-resolved direction field *is* noise at the grid scale.
**Any measured $\beta$ below $0.5$ must be checked against $k_{\max}\eta$ before it
is interpreted**, and the preregistration's H4 falsification condition is
conditional on the resolution gate passing. The same caution applies to
$\mu_\xi$, whose log-log slope has the same white-noise limit.

`diagnostics.ic_resolution_check` reports $k_{\max}\eta$ at $t=0$ and the
viscosity that would be required, and `run_campaign.py` prints a warning before
integrating. The check belongs *before* the run: for a broadband initial spectrum
the minimum of $k_{\max}\eta$ over the whole run is at $t=0$, and the flow becomes
better resolved as it decays — so a run can be corrupt from the first step and
look fine at the end.

## 7. Uncertainty

Spatial CI by bootstrap over the 8 (or 64) sub-cubes; temporal CI by bootstrap
over snapshots. **Never quote a dimension, $\delta$, $\beta$ or $\mu_\xi$ without
a CI.** Bootstrap over sub-cubes underestimates uncertainty when the structure
spans the box; note this where it applies.
