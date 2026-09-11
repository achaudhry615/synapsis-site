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
python3 scripts/validate_pipeline.py     # all PASS required before trusting anything else
python3 scripts/run_tg64.py              # 64³ Taylor–Green, full Γ pipeline
python3 scripts/evaluate_gates.py results/<run>.jsonl
```
