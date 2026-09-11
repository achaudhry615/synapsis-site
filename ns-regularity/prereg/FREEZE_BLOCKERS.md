# What blocks freezing the preregistration

The preregistration cannot be frozen from this runtime. Every publisher host is
egress-blocked (arXiv, nature.com, IOPscience, link.springer.com, orbit.dtu.dk,
semanticscholar.org, doi.org, pmc.ncbi.nlm.nih.gov, personal university pages —
all refuse the CONNECT tunnel), and `WebFetch` is allowlisted to essentially
nothing. Only a web *search* tool is available, whose summariser can read pages
we cannot fetch; that was enough for bibliographic data and for some stated
hypotheses, and it already caught one substantive error (below).

Verified entries: **18 of 39** (`paper/references.md`).

## Blockers, in priority order

Each is a specific question, answerable by reading one paper.

### 1. Albritton–Bradshaw, *Nonlinearity* **35** (2022) 2858 — arXiv:2110.02187
**The single most important unread source.** Two questions:
- (a) What exactly does their "simple proof that sufficiently sparse
  Navier–Stokes solutions do not develop singularities" assume? Theorem C3 is a
  reduction to it, so its hypothesis *is* C3's hypothesis.
- (b) What is their verdict on the claims that a priori sparseness estimates
  reduce the scaling gap? They examine Bradshaw–Farhat–Grujić (2019) and
  Grujić–Xu; Bradshaw co-authored the 2019 paper. **Until this is read, A14 stays
  CONTESTED and is relied upon nowhere.**

### 2. Grujić, *Nonlinearity* **26**(1) (2013) 289 — arXiv:1111.0217
- (a) The numerical threshold $\delta_0$ below which sparseness implies
  regularity. **`H2`'s threshold is UNSET pending this**, and setting it from our
  own measurements would make H2 unfalsifiable — so H2 is currently descriptive,
  not a test.
- (b) The exact constant relating the sparseness scale to the radius of spatial
  analyticity, and the a priori lower bound on that radius in terms of
  $\|\omega\|_\infty$. Confirmed so far: the scale is *comparable to the
  analyticity radius*, not equal to $\ell_\nu$ by construction.

### 3. Constantin, *SIAM Review* **36** (1994) 73–98
The exact kernel bound $|D(\hat y,\xi',\xi)|\le c|\sin\angle(\xi',\xi)|$ in the
singular-integral representation of $\alpha$. This is the mechanism every
surviving route depends on and the structure Tao's averaging destroys. Used only
qualitatively here; no constant from it enters any argument.

### 4. Beirão da Veiga–Berselli, *Diff. Integral Eq.* **15** (2002) 345–356
Confirmed: Lipschitz can be relaxed to $1/2$-Hölder, giving
$\omega\in L^\infty_tL^2_x$ and hence smoothness — this is $\theta_c=1/2$ in
Theorem C4. **Unconfirmed:** the precise region and uniformity over which the
Hölder bound is required. C4 states it on an $r$-neighbourhood of $E_\lambda$;
that form needs checking.

### 5. Grujić, arXiv:2607.08866 (July 2026)
Theorem C5 reduces to it. Need: the exact definition of the log-weighted space,
and the precise $L^{3/2,\infty}$ concentration hypothesis. Also worth reading:
arXiv:2609.05720 (Sept 2026), on decay of local mean oscillations of the
vorticity direction — directly relevant to the $\mu_\xi$ measurement.

### 6. Rafner et al., *Sci. Rep.* **11** (2021) 8824 — arXiv:1909.10408
The only prior DNS measurement of the sparseness scale. Needed for a prior on
$\delta_0$ and as a calibration target for our estimator. Reported qualitatively
that in the Kida-vortex case the sparseness scale was actualised "well beyond the
guaranteed a priori bound and just beyond the critical bound sufficient for
diffusion to fully engage"; the numbers behind that were not obtainable.

## Why this list is worth keeping

The verification pass already paid for itself. It found that Grujić's criterion
is stated on the super-level sets of **the positive and negative parts of the
vorticity components** — six sets — and not on $\{|\omega|>\lambda M\}$, which is
what Theorem C3 and the code had used. Measured at $64^3$ on the smoothed vortex
sheet, the wrong definition gave $\delta_{p95}=0.061$ against $0.152$ for the
worst component set: it reported the configuration as **thinner than it is**, and
so was biased toward falsely confirming the hypothesis. One paper, read
properly, would have prevented it. Five more are unread.
