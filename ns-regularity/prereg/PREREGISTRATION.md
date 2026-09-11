# Preregistration — NS-Γ geometric measurement campaign

**Status: DRAFT. Not yet frozen.** Freezing requires (i) the citation
verification pass of `paper/references.md` to complete, so that H2's threshold
can be set from the published value rather than from our own data, and (ii) a git
tag `prereg-v1` recording the SHA-256 below. Until then no run executed against
this document may be described as preregistered.

## 0. Why this document exists

The two prior gauntlet rounds reported "zero surviving blow-up candidates" from
ensembles that had no power to find one. That is not a null result; it is an
absence of measurement. Fixing the hypotheses, thresholds and rejection rules
*before* the run is what separates a falsifiable campaign from a fishing
expedition, and it is the only reason a null outcome here would mean anything.

## 1. Hypotheses

Each is stated with a numeric falsification condition.

### H1 — null search (secondary endpoint)
**Prediction:** zero cases satisfy `STRICT_SURVIVING_BLOWUP_CANDIDATE`.
**Falsified by:** any survivor that replicates under a second solver and a
different seed.
**Note:** H1 is *expected to hold*. Confirming it is uninformative and will be
reported as such. It is preregistered so that a surprise cannot be dismissed.

### H2 — sparseness (Theorem C3)
**Prediction:** for every IC and every $t$ in the resolved window, the 95th
percentile of the per-point 1D sparseness of $\Omega_{0.5}(t)$ at $r=\ell_\nu$
satisfies $\delta_{p95}\le\delta_0$.
**$\delta_0$: UNSET.** It must be taken from the published threshold of the
sparseness criterion (A13/A14). **Setting $\delta_0$ from our own measurements
would make H2 unfalsifiable and is forbidden.** If the published threshold cannot
be established, H2 is reported as *descriptive only* and not as a test.
**Falsified by:** any resolved run exceeding $\delta_0$ with $N$-independence
demonstrated per `gamma_measurement.md` §6.

### H3 — concentration geometry
**Prediction:** $D_\infty$ of the $\Gamma$ measure lies in $[0.8,1.5]$
(tube-like) and is $N$-independent at matched physical scales.
**Falsified by:** $D_\infty\ge2$ with CI excluding 1.5, or drift with $N$ at
matched physical scales exceeding the calibrated 0.15.
**Caveat fixed in advance:** $D_\infty$ is biased low by 0.3–0.5
(`gamma_measurement.md` §5.4). The interval above is stated for the *estimator*,
not the true dimension.

### H4 — direction regularity (Theorem C4)
**Prediction:** $\beta\ge0.5$ on a neighbourhood of $E_\lambda$.
**Falsified by:** $\beta<0.5$ with CI excluding 0.5 at two consecutive $N$.
**Caveat fixed in advance:** the estimator is biased low by ~10% (it returns 0.90
on a field whose true $\beta$ is 1). A measured 0.45 is therefore *not* a
falsification; the bias-corrected comparison is preregistered here.

### H5 — depletion trend (Conjecture D2)
**Prediction:** $\rho(\ell_\nu,t)$ is non-increasing as $\nu$ decreases along
Ladder G, i.e. $\theta\ge0$.
**Falsified by:** monotone increase across the three viscosities with
non-overlapping CIs.

### H6 — log-modulus (Theorem C5) — the primary scientific question
**Prediction:** $\mu_\xi(r)\,|\log r|$ remains bounded as $r$ decreases over the
resolved range, and does not grow as $\nu$ decreases.
**Falsified by:** monotone growth of $\mu_\xi|\log r|$ in $r\to0$ that persists
under refinement.
**This is the hypothesis whose falsification would be most informative**, because
it is the weakest of the three sufficient conditions: if the log-modulus bound
fails on the stretching-intense set, the A16 route to depletion closes and the
geometry that defeats it is worth characterising.

## 2. Fixed in advance

- **Solver:** integrating-factor RK4, rotational form, 2/3 dealiasing, FP64,
  CFL ≤ 0.4. No hyperviscosity.
- **Seeds:** fixed per IC family; recorded in the run log meta record.
- **IC parameter tables:** frozen in `code/nsref/ic.py` at the tagged commit.
- **Analysis code version:** the git SHA recorded in every run log's meta record.
- **Gates:** `prereg/gates.json`, 14 gates / 28 criteria, fail-closed.
- **Thresholds:** as in `gates.json` and `spec/model_selection.md`.
- **Estimator calibrations and their known biases:** `spec/gamma_measurement.md`
  §5, measured before the campaign on fields of known geometry.

## 3. Stopping rules

- A case stops at $t_{\rm end}$ = the last time `HIGH_K_TAIL_RESOLVED` holds.
  Data past loss of resolution is discarded, not fitted.
- A case escalates to the next $N$ only if all applicable gates passed at the
  current $N$.
- The campaign stops when Ladder G completes at $512^3$ for all families, or on
  exhaustion of the compute budget — whichever first. Partial completion is
  reported as partial.

## 4. What will be reported regardless of outcome

All six hypotheses with measured values and CIs; the full gate table per case
including failures; all four fitted models per case including when the power law
loses; every case that was run, including those that crashed or failed gates;
and the deviations log.

**In particular:** a null result for H1 will be reported as *uninformative*, not
as evidence for regularity. The inference "no candidate found, therefore no
singularity exists" is invalid and will not be made.

## 5. Deviations

Any change to this document or to `gates.json` after freezing requires a dated
entry in `DEVIATIONS.md` stating what changed and why. Changes to thresholds
invalidate any candidate claim made under the old values.

## 6. Integrity

On freezing, record here: `git tag prereg-v1`, the commit SHA, and
`sha256sum PREREGISTRATION.md gates.json`.

```
SHA-256 (PREREGISTRATION.md): <set on freeze>
SHA-256 (gates.json):         <set on freeze>
git tag:                      <set on freeze>
```
