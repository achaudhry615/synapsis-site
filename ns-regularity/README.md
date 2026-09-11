# NS-Γ: a proof-search architecture for scale-local vortex-stretching depletion

## Status statement (read before anything else)

**Global regularity for the 3D incompressible Navier–Stokes equations is an open
problem. Nothing in this repository is a proof of it, and nothing here should be
cited as one.**

What this repository contains is a *proof-search architecture*: a precisely
stated set of conditional theorems reducing regularity to measurable geometric
hypotheses, an explicit ledger of which inequalities fail to close and why, a
measurement pipeline for those hypotheses, and a preregistered numerical
campaign that can falsify them.

Every claim in this repository carries one of five labels:

| Label | Meaning |
|---|---|
| **A** | Standard published mathematics, cited. Not our contribution. |
| **B** | No-go result: a route we prove *cannot* close, with the obstruction identified. |
| **C** | Conditional theorem: correct, proved here by reduction to A-items, with an explicitly unproved hypothesis. |
| **D** | Conjecture: stated so that it can be false. Not proved, not assumed. |
| **N** | Numerical measurement or preregistered prediction. Evidence, never proof. |

There is no label for "proved". That is deliberate.

## What this program actually targets

The originating question was: *does concentration geometry enforce depletion of
the globally integrated vortex stretching?* The honest answer developed in
`paper/ns_gamma_codimension.md` is:

1. **Not through codimension alone.** We prove (Part B) that a hypothesis of the
   form "the stretching-intense set has box dimension `D`" cannot improve the
   textbook enstrophy bound `Ω' ≤ CΩ³/ν³` for any value of `D`. The closing
   inequality is scale-invariant at the dissipation scale: the volume of the
   concentration set cancels between the stretching and dissipation terms. This
   is a negative result about the method, and it is the reason the program was
   re-anchored.

2. **Plausibly through transverse thinness and direction regularity.** Dimension
   is neither necessary nor sufficient; *transverse thinness* (1D sparseness in
   the sense of Grujić) and *regularity of the vorticity direction field* `ξ` on
   the intense set are the two hypotheses that genuinely appear in published
   sufficient conditions. Those are what the campaign measures.

## The measurable quantities

| Symbol | Definition | Connected theorem |
|---|---|---|
| `Γ(x,t)` | `max(ω·Sω, 0)` — pointwise positive vortex-stretching density | defines the sets below |
| `δ(r)` | 1D sparseness ratio of super-level sets at scale `r ~ ℓ_ν` | Grujić sparseness criterion (A9) |
| `β` | Hölder exponent of `ξ` on a neighbourhood of the intense set | Constantin–Fefferman / Beirão da Veiga–Berselli, critical `β = 1/2` (A6) |
| `μ_ξ(r)` | local mean oscillation of `ξ`, and its modulus as `r → 0` | log-weighted bmo depletion (A13) |
| `R_E(t)` | stretching/dissipation ratio *restricted to the intense set* | the object the enstrophy balance actually turns on (B2) |

`R_E` replaces the global ratio `R` used in the originating runs. A global
`R > 1` is generic in any enstrophy-increasing flow and carries no information.

## Layout

```
paper/     analytic skeleton: standard results, no-go ledger, conditional theorems, conjectures
spec/      measurement specification, 4-GPU campaign design, model-selection procedure
prereg/    preregistered hypotheses and numeric gates, fixed before any production run
graph/     machine-readable proof-obligation graph (nodes, edges, status)
code/      runnable reference implementation (numpy); GPU path specified separately
results/   run outputs (JSONL logs, gate tables)
```

## What this repository actually establishes

Nothing about regularity. What it does establish, and what took the work:

**One analytic result (proved, and machine-checked).** Proposition B1: a
box-dimension hypothesis on the stretching-intense set cannot improve the
enstrophy bound for *any* codimension. The $D=0$ case reproduces Doering–Gibbon
exactly, every $D>0$ is strictly worse, and $D\ge1$ — tubes and sheets, i.e.
every scenario of concern — yields no bound at all. The exponents are verified
numerically in `validate_pipeline.py` §6. This answers the originating question
in the negative and is the reason the program was re-anchored on transverse
thinness and direction regularity.

**Three conditional theorems** (C3, C4, C5), each a correct reduction of
regularity to a numerically measurable hypothesis, with the gap to a theorem
isolated in a single named node (`NS-G01`).

**Five measurement findings**, each derived rather than assumed, and each of
which would have silently corrupted the campaign:

1. **$\ell_\nu/\Delta x$ depends only on $k_{\max}\eta$, not on $N$.** The standard
   DNS criterion $k_{\max}\eta\ge2$ leaves the Theorem-C3 hypothesis *sub-grid*,
   where $\delta\equiv1$ and carries no information. Confirmed on the $64^3$ runs.
   Measuring it needs $k_{\max}\eta\sim10$–$20$ — so the geometry measurement and
   the blow-up search pull $\nu$ in opposite directions by an order of magnitude
   and **cannot be done by one run**.
2. **The C5 modulus must be weighted by $\log(L/r)$, not $|\log r|$.** The latter
   vanishes at $r=1$ and inverted the verdict in calibration — reporting a smooth
   field as unbounded and white noise as bounded, in the campaign's primary endpoint.
3. **Dimension box scales must be fixed in physical units.** Grid-relative ranges
   make estimates drift by up to 0.83 across a factor-2 refinement; matched
   physical scales hold them to 0.147.
4. **The proposed absolute $\Gamma$ thresholds emptied the set for every family**,
   and would have emptied it soonest for the *most depleted* flows — biasing the
   campaign against the very effect it exists to measure.
5. **$D_0$ on a smooth positive field is 3 identically** (full support); it must be
   computed on a thresholded mask.

**A validated implementation.** 33/33 pipeline checks, 13/13 unit tests. The
sharpest check: ABC flow is Beltrami ($\nabla\times u=u$), so $u\times\omega\equiv0$
and $u(t)=u_0e^{-\nu t}$ is an *exact* solution of the full equations — the solver
reproduces it to $8.3\times10^{-16}$ over 30 RK4 steps.

## Evidential status of the prior runs

The two Monte-Carlo rounds that motivated this repository (16 cases at 16³–32³;
24 adversarial cases at 24³–64³, peak vorticity growth 1.537×, max FTLE 9.93,
zero surviving candidates) are recorded in the proof graph as
`NS-X01: SMOKE_TEST_ONLY`. They validate that a pipeline runs. They carry **no
evidential weight** for or against regularity, for the reasons given in
`paper/ns_gamma_codimension.md` §6:

- 64³ admits no inertial range, so the structure-function exponents `ζ_p` are noise;
- 1.537× peak-vorticity growth is less than undriven Taylor–Green produces;
- a transient `R > 1` is implied by `Ω' > 0` and is not a singularity signal;
- FTLE measures Lagrangian chaos, which is known to coexist with global regularity.

Reporting them as "no blow-up candidate found" is true but vacuous: the ensemble
had no power to find one.

## Build order

1. `graph/` — fixes the vocabulary
2. `paper/` §0–§3 — the no-go ledger, written *before* the conditional theorems
3. `paper/` §4–§5 — conditional theorems and conjectures, which must follow the ledger
4. `code/nsref/gamma.py` + validation on synthetic fields of known geometry
5. `code/nsref/solver.py` + dynamic validation
6. `spec/gamma_measurement.md` — written after the code, to match what is computed
7. `prereg/` — numeric gates, fixed before any production run
8. `spec/campaign_4gpu.md`

## Reproducing

```
cd code
python3 scripts/validate_pipeline.py      # 33/33 required before trusting anything else
python3 -m unittest discover -s tests     # 13/13
python3 scripts/resolution_budget.py      # what resolution the geometry actually needs
python3 scripts/run_campaign.py --ic taylor_green --N 64 --nu 0.02 --tend 6
python3 scripts/summarize_results.py      # cross-case table
python3 -m gpu.ns_cupy --selftest         # CPU/GPU equivalence (needs CUDA for the real test)
```

`run_campaign.py --ic` accepts `taylor_green`, `abc`, `kida_pelz`,
`antiparallel_tubes`, `vortex_sheet`, `multiscale_random`.

Requires numpy; SciPy is used for the Levenberg–Marquardt fits when present and
a Nelder–Mead fallback is used when it is not.
