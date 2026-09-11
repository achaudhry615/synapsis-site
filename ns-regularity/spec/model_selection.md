# Model selection for $\|\omega\|_\infty(t)$

Implemented in `code/nsref/fit.py`; validated in `validate_pipeline.py` §5.

## Why four models

A power law fits almost any convex increasing curve over a short window. Fitting
$A(T^*-t)^{-\alpha}$ alone and reporting the resulting $T^*$ is the single most
common way to manufacture a blow-up that is not there. The comparison against
alternatives is what makes the fit informative, so **all four models are always
fitted and all four are always reported, including when the power law loses.**

| id | model | params |
|---|---|---|
| M1 | $A(T^*-t)^{-\alpha}$ | 3 |
| M2 | $A e^{bt}$ | 2 |
| M3 | $A\exp(bt^{c})$ (stretched exponential) | 3 |
| M4 | $A\exp(Be^{ct})$ (double exponential) | 3 |

## Procedure

1. **Data.** $y=\log\|\omega\|_\infty(t)$, with $\|\omega\|_\infty$ from
   *spectrally interpolated* maxima (`diagnostics.omega_max_refined`, 2× zero
   padding). The grid maximum hops between cells as a vortex translates,
   injecting ~1% noise that destroys nonlinear fits. Never fit grid maxima.
2. **Window.** $[t_1,t_{\rm end}]$ with $t_{\rm end}$ the last time
   `HIGH_K_TAIL_RESOLVED` holds. Data past loss of resolution is not data.
3. **Fit.** Levenberg–Marquardt in log space (bounded: $\alpha\in[0,20]$,
   $T^*>t_{\rm end}$). A derivative-free Nelder–Mead fallback is provided so the
   pipeline runs without SciPy.
4. **Selection.** AICc and BIC. M1 is *preferred* only if it beats **every**
   alternative by $\Delta\mathrm{BIC}\ge10$.
5. **Uncertainty.** Wild (Rademacher) residual bootstrap, $B=2000$ in production,
   percentile CIs on $T^*$ and $\alpha$.
6. **Stability.** Refit on windows truncated to $\{1.0,0.9,0.8\}\times t_{\rm end}$
   (Hou–Li drift test) and across $N\in\{128,256,512\}$.

## Preregistered acceptance

| quantity | threshold | rationale |
|---|---|---|
| power-law BIC margin | $\ge10$ over every alternative | decisive evidence |
| $\alpha$ CI lower bound | $\ge1$ | Leray lower bound (A12) forbids $\alpha<1$ |
| $T^*$ drift across windows | $\le2\%$ | a real $T^*$ does not move |
| $T^*$ drift across $N$ | $\le2\%$ | a real $T^*$ is not a resolution artefact |
| $\alpha$ CI width | $\le0.2$ | otherwise the exponent is not determined |
| extrapolation ratio $(T^*-t_{\rm end})/(t_{\rm end}-t_1)$ | $\le0.5$ | do not extrapolate past the data |

## Validation

On synthetic data with $T^*=2.0$, $\alpha=1.6$ and 1% multiplicative noise:

- preferred model = `power_law`, BIC margin **200.6**
- $T^*=2.0008$, CI$_{95}$ $[1.996,2.006]$ (true value inside)
- $\alpha=1.6037$, CI$_{95}$ $[1.591,1.609]$
- window drift $0.25\%$

On synthetic **exponential** data (no singularity):

- preferred model = `exponential`, power-law BIC margin $-58.7$ → correctly rejected
- the spurious $T^*$ drifts by only $1.2\%$ across windows

**That last number matters.** Window drift alone does *not* discriminate: an
exponential can yield a stable-looking spurious $T^*$. The BIC margin is the
load-bearing criterion, and drift is a secondary check. A campaign that used
drift alone would have accepted the exponential.
