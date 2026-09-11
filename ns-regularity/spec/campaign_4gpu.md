# 4-GPU campaign specification

## 1. What this campaign is for

**It is not a blow-up hunt.** Paper §6.3: with $k_{\max}\eta\ge2$ enforced,
$512^3$ admits $\nu\sim10^{-4}$–$10^{-3}$ and $\mathrm{Re}_\lambda\sim10^2$.
Viscous Navier–Stokes at that viscosity will not develop a singularity, so a
campaign whose primary endpoint is "find a surviving candidate" is designed to
return null whatever the truth, and its null result is uninformative — exactly
the defect of the two prior gauntlet rounds.

| endpoint | content | expected result |
|---|---|---|
| **primary** | measure $\delta(r/\ell_\nu)$, $\beta$, $\mu_\xi(r)|\log r|$, $\mathcal R_E$, $\rho(\ell_\nu)$ and their trends as $\nu\downarrow$ and under refinement | informative either way |
| **secondary** | blow-up candidate search through the full gate ladder | **preregistered expected count: 0** |

The secondary endpoint is retained only so the gates are exercised and a surprise
would be caught.

## 2. The resolution conflict, and how to resolve it

From `gamma_measurement.md` §4: $\ell_\nu/\Delta x$ depends only on
$k_{\max}\eta$, **not on $N$**. Measuring $\delta$ at $r\sim\ell_\nu$ needs
$k_{\max}\eta\sim10$–$20$; a blow-up search wants the smallest $\nu$ the grid
allows, i.e. $k_{\max}\eta\approx2$. **These cannot be satisfied by one run.**

Run two ladders, labelled, and never mix their outputs:

**Ladder G (geometry — primary).** $k_{\max}\eta\ge10$ enforced. $N=128,256,512$
with $\nu$ from the table in `gamma_measurement.md` §4 ($3.5\times10^{-2}$,
$8.6\times10^{-3}$, $2.2\times10^{-3}$ at $\Omega_v\approx1$; recompute from the
measured $\varepsilon$ at runtime). Purpose: are the C3/C4/C5 hypotheses
satisfied, and do they strengthen or weaken as $\nu$ falls?

**Ladder C (convergence — gates).** Fixed $\nu$ per IC, $N=128\to256\to512$.
Purpose: `SPATIAL_REFINEMENT_CONVERGED`, `TIME_STEP_CONVERGED`,
`PRECISION_ESCALATION_CONSISTENT`, `DOMAIN_SIZE_CONSISTENT` and the $T^*$
stability comparison. Escalate a case to the next $N$ only if it passed every
applicable gate at the current $N$.

## 3. Initial conditions

Seeded and deterministic. Random multiscale spectra are *not* adversarial; the
families below have either a theorem or a published benchmark behind them.

| family | why | status in `nsref/ic.py` |
|---|---|---|
| Taylor–Green | benchmark; analytic $t=0$ values verify the solver | implemented |
| ABC / near-Beltrami | zero nonlinearity limit; $\beta$ calibration | implemented |
| anti-parallel tubes | the Kerr geometry; the classical near-singular candidate | implemented (Biot–Savart from imposed $\omega$) |
| smoothed vortex sheet | $D=2$ concentration | implemented |
| multiscale random (seeded) | regression case for the prior run's ADV023 | implemented |
| **Lu–Doering / Ayala–Protas extreme states** | *variationally optimal* enstrophy growth — saturate the A4 bound | **TODO: gradient ascent, or start from the published analytic candidates** |
| Kida–Pelz high symmetry | classical near-singular; symmetry gate applies | implemented (exact analytic $E=3/8$, $\Omega=33/8$, $\|\omega\|_\infty=8$) |
| **colliding vortex rings** | Kerr's $\sqrt{\nu}$ circulation scaling is a solver-validation target | **TODO** |

The two remaining TODO families are the ones most likely to matter and are the
first extension; the implemented six suffice to exercise the pipeline.

**Note on $\nu$:** at $64^3$ the smooth families (antiparallel tubes) run at
$k_{\max}\eta\approx14$–$22$ with $\nu=0.02$ — far over-resolved, and hence the
only cases in which $\delta(r\!\sim\!\ell_\nu)$ was measurable at all
($\ell_\nu/\Delta x\approx2$–$3$). That is the resolution trade-off of §2 showing
up in practice: over-resolution is what buys the geometry measurement.

## 4. Solver and hardware

- Port `nsref/solver.py` to CuPy (`code/gpu/ns_cupy.py`): identical algorithm —
  rotational form $P[u\times\omega]$, 2/3 truncation, integrating-factor RK4,
  FP64, CFL $\le0.4$. **No hyperviscosity, ever.**
- RK2 (used in the prior rounds) is not acceptable: its $O(\Delta t^2)$ phase
  error corrupts $\|\omega\|_\infty(t)$ and therefore every $T^*$ fit.
- $128^3$/$256^3$: one case per GPU, embarrassingly parallel across ICs.
- $512^3$ FP64: ~20 real fields × 1.07 GB ≈ 22–30 GB with FFT workspace. Single
  GPU if ≥40 GB VRAM; otherwise FP32 at $512^3$ with
  `PRECISION_ESCALATION_CONSISTENT` mandatory at $256^3$. A slab decomposition
  with NCCL is possible but is not recommended — the ladder is cheap enough that
  the engineering risk is not worth it.
- Cost: ~1 s/step at $512^3$, 1–3k steps per case → ~1 h/case. **I/O is the
  bottleneck**, not compute: a $512^3$ velocity snapshot is ~3 GB. Store spectra,
  $\Gamma$ statistics and the geometry tables at every log step; full snapshots
  at ~10 preselected times only.
- `MULTIPLE_SOLVERS_AGREE`: an independent second solver (e.g. spectralDNS with
  3/2-padding instead of 2/3-truncation) at $128^3$/$256^3$. Agreement:
  $\|\omega\|_\infty$ within $10^{-3}$, $\int\Gamma$ within $10^{-2}$.

## 5. Gates

Defined numerically in `prereg/gates.json` (14 gates, 28 criteria), evaluated by
`nsref/gates.py`. Gates **fail closed**: a quantity that was not measured is not
a pass. A `STRICT_SURVIVING_BLOWUP_CANDIDATE` requires every gate at $512^3$ plus
the model-selection criteria in `model_selection.md`.

On a nonzero count: independent replication with a second solver and a different
seed is required before any public statement stronger than "unexplained numerical
behaviour". **A surviving candidate is not a proof and must never be reported as
one.**

## 6. What to do with the result

- **If $\delta(r/\ell_\nu)$ stays well below the A13 threshold and $\mu_\xi|\log r|$
  stays bounded, both strengthening as $\nu\downarrow$:** evidence for the
  depletion mechanism of A16. The analytic target becomes proving the log-modulus
  bound propagates — i.e. attacking `NS-G01` through C5 rather than C3.
- **If $\mu_\xi|\log r|$ grows under refinement on the $\Gamma$-set:** the more
  interesting outcome. It localises where a singularity would have to hide, and
  the next campaign targets that geometry specifically.
- **Either way**, report the trend, the CIs, and the gates. "No candidate found"
  is not a result.
