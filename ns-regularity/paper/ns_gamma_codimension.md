# Scale-local vortex stretching, concentration geometry, and the 3D Navier–Stokes regularity problem

**A conditional-theorem skeleton with an explicit no-go ledger.**

---

## §0. Status

Global regularity for the 3D incompressible Navier–Stokes equations is an open
problem. **This document does not prove it, does not claim to prove it, and
contains no unconditional new theorem.**

What it does contain:

- **Part A** (§2): standard published results, stated with citations, used as given.
- **Part B** (§3): a *no-go ledger* — four propositions establishing that the
  codimension route as originally posed cannot close, with the obstruction
  located precisely in each case. These are proved.
- **Part C** (§4): three *conditional theorems* — correct reductions of regularity
  to geometric hypotheses that are numerically measurable. Each is proved modulo
  an explicitly stated, explicitly unproved hypothesis.
- **Part D** (§5): two *conjectures*, stated so that they can be false, together
  with the identification of the missing lemma in each.
- **Part N** (§6): what the numerics can and cannot establish.

The gap between Part C and a theorem is node `NS-G01` in `graph/nodes.csv`:
*derive the geometric hypothesis from the equations*. No path to it is proposed
here. It is the Millennium problem in different coordinates.

---

## §1. Setting and notation

On the periodic box $\mathbb{T}^3_L=[0,L)^3$ with kinematic viscosity $\nu>0$,

$$\partial_t u + (u\cdot\nabla)u = -\nabla p + \nu\Delta u,\qquad \nabla\cdot u = 0,$$

with $u(\cdot,0)=u_0$ smooth, mean zero, and $u$ smooth on a maximal interval
$[0,T^*)$.

| Symbol | Definition |
|---|---|
| $\omega$ | $\nabla\times u$ |
| $S$ | $\tfrac12(\nabla u + \nabla u^{\mathsf T})$, the strain-rate tensor |
| $\xi$ | $\omega/|\omega|$, defined where $\omega\neq 0$ |
| $\alpha$ | $\xi\cdot S\xi$, the stretching eigenvalue |
| $M(t)$ | $\|\omega(t)\|_{L^\infty}$ |
| $\Omega(t)$ | $\tfrac12\|\omega(t)\|_{L^2}^2$, enstrophy |
| $Y(t)$ | $\|\nabla\omega(t)\|_{L^2}^2$, palinstrophy |
| $\ell_\nu(t)$ | $(\nu/M(t))^{1/2}$, the vorticity-based dissipation scale |
| $\Gamma$ | $\max(\omega\cdot S\omega,\,0) = |\omega|^2\max(\alpha,0)$ |
| $\Omega_\lambda(t)$ | $\{x : |\omega(x,t)| > \lambda M(t)\}$ |
| $E_\lambda(t)$ | $\{x : \Gamma(x,t) > \lambda M(t)^3\}$ |

**Remark 1.1 (why $\Gamma$ and not $\xi^{\mathsf T}S\xi$).** The originating
formulation used $\Gamma=|\omega|^2\max(\xi^{\mathsf T}S\xi,0)$ with $\xi$
normalised. The normalisation is redundant — $|\omega|^2\,\xi^{\mathsf T}S\xi
=\omega\cdot S\omega$ identically — and it introduces a spurious $0/0$ at
vorticity nulls, which are generic and which dominate the measured statistics if
not excluded by hand. Use $\Gamma=\max(\omega\cdot S\omega,0)$. Record
$\Gamma^-=\max(-\omega\cdot S\omega,0)$ separately: $\int\Gamma^+$ alone
overestimates production, and the *net* $\int\omega\cdot S\omega$ is the quantity
appearing in the enstrophy balance.

---

## §2. Part A — standard results

Each item is used as given. Items tagged **[VERIFY]** are ones whose exact
hypotheses or constants could not be checked against the primary source from the
runtime in which this document was drafted (arXiv and nature.com were
unreachable); they are stated from secondary sources and must be confirmed
against the published statements before any of Part C is relied upon.

**A1 (Leray energy identity).** $\frac{d}{dt}\tfrac12\|u\|_2^2 = -\nu\|\nabla u\|_2^2 = -2\nu\Omega$.

**A2 (vorticity equation).** $\partial_t\omega + (u\cdot\nabla)\omega = S\omega + \nu\Delta\omega$, and along particle paths
$$\frac{D|\omega|}{Dt} = \alpha|\omega| + \nu\,\xi\cdot\Delta\omega .$$

**A3 (enstrophy balance).** Exactly,
$$\Omega'(t) = \int \omega\cdot S\omega\,dx \;-\; \nu Y(t). \tag{2.1}$$
This identity, not any inequality, is what the whole program turns on.

**A4 (supercritical enstrophy bound).** $\Omega' \le -\tfrac{\nu}{2}Y + C\nu^{-3}\Omega^3$ (Doering–Gibbon). The instantaneous rate is saturated by explicit vortex states (Lu–Doering 2008), so the inequality is sharp as an instantaneous statement; whether it is saturable *along a trajectory* is open (Ayala–Protas 2017; Kang–Yun–Protas 2020).

**A5 (Calderón–Zygmund and the logarithmic loss).** $\|S\|_{L^p}\le C_p\|\omega\|_{L^p}$ for $1<p<\infty$, with $C_p\to\infty$ as $p\to\infty$. At $p=\infty$ the Beale–Kato–Majda-type inequality
$$\|S\|_{L^\infty}\;\le\; C\,\|\omega\|_{L^\infty}\Big(1+\log^+\tfrac{\|\omega\|_{H^s}}{\|\omega\|_{L^\infty}}\Big)\;=:\;C_{\log}\,\|\omega\|_{L^\infty},\quad s>5/2, \tag{2.2}$$
replaces it. $C_{\log}$ is not bounded a priori; it is measurable.

**A6 (BKM criterion).** $T^*<\infty \iff \int_0^{T^*}\|\omega(t)\|_{L^\infty}dt=\infty$ (Beale–Kato–Majda 1984).

**A7 (Constantin's representation).** $\alpha(x)=\frac{3}{4\pi}\,\mathrm{P.V.}\!\int D\big(\hat y,\xi(x+y),\xi(x)\big)\,|\omega(x+y)|\,\frac{dy}{|y|^3}$, with the kernel satisfying $|D|\le c\,|\sin\angle(\xi(x+y),\xi(x))|$ (Constantin 1994). **This is the mechanism.** Stretching is not merely large-where-vorticity-is-large; it is suppressed wherever the direction field is locally coherent. Every viable depletion route in the literature goes through this bound.

**A8 (Constantin–Fefferman).** If $\xi$ is Lipschitz uniformly on the region where $|\omega|$ is large, the solution stays regular (Constantin–Fefferman 1993).

**A9 (Beirão da Veiga–Berselli).** Lipschitz can be relaxed to Hölder-$\beta$ with $\beta\ge 1/2$ [VERIFY: *Diff. Int. Eq.* **15** (2002) 345–356]. **$\beta=1/2$ is scale-critical and is the only place in this entire program where a genuine exponent threshold exists.** Any statement of the form "$\theta>\theta_c$" must ultimately be a statement about $\beta$, not about a dimension.

**A10 (Caffarelli–Kohn–Nirenberg).** The space-time singular set $\Sigma$ satisfies $\mathcal P^1(\Sigma)=0$ (CKN 1982; Lin 1998). **Caution:** $\mathcal P^1$ is *parabolic* Hausdorff measure of the singular set *at* $T^*$. It is not the box dimension of a pre-singular super-level set, and the two must not be conflated. In particular CKN does not exclude tube-like ($D=1$) or sheet-like ($D=2$) concentration prior to $T^*$, which is exactly the geometry of every scenario of concern.

**A11 (LPS / ESS).** Ladyzhenskaya–Prodi–Serrin conditions; $u\in L^\infty_t L^3_x$ implies regularity (Escauriaza–Seregin–Šverák 2003), and $\|u(t)\|_{L^3}\to\infty$ at a blow-up time (Seregin 2012).

**A12 (Leray lower bound).** $\|\omega(t)\|_{L^\infty}\ge c\,(T^*-t)^{-1}$. Hence any power-law fit $\|\omega\|_\infty\sim A(T^*-t)^{-a}$ to credible blow-up data must return $a\ge 1$.

**A13 (Grujić sparseness criterion).** [VERIFY: Grujić 2013, arXiv:1111.0217] A local anisotropic geometric-measure condition: if the suitably cut vorticity super-level sets are *1D $\delta$-sparse* at a scale comparable to the radius of spatial analyticity (itself bounded below in terms of $\|\omega\|_\infty$), no finite-time blow-up occurs. The proof is via analyticity plus the harmonic-measure majorisation principle: sparse sets carry small harmonic measure, which forbids growth of the $L^\infty$ norm. **This is the rigorous form of "concentration geometry implies regularity", and it is the correct anchor for this program.** [VERIFY the exact exponent relating the sparseness scale to $\|\omega\|_\infty$, and the threshold value of $\delta$.]

**A14 (algebraic reduction of the scaling gap).** Bradshaw–Farhat–Grujić, *ARMA* (2019): within the sparseness framework the scaling gap is reduced by an *algebraic* factor — the first non-logarithmic reduction since the 1960s. [VERIFY exact exponents.] Related: Farhat–Grujić–Leitmeyer (volumetric sparseness, $B^{-1}_{\infty,\infty}$); Grujić–Xu (asymptotic criticality).

**A15 (Cheskidov–Shvydkoy).** Regularity criteria in $B^{-1}_{\infty,\infty}$ [VERIFY]. This is the **only** rigorous bridge from spectral/band-limited data to a regularity theorem, and it is the reason the multifractal branch of the original program is downgraded (§6.2).

**A16 (log-weighted bmo depletion).** [VERIFY: Grujić, arXiv:2607.08866, July 2026] For critical point singularities exhibiting $L^{3/2,\infty}$ concentration of vorticity: if the vorticity *direction* lies locally in a logarithmically weighted BMO space $\mathrm{bmo}_{1/|\log r|}$ — a space that fails the Dini condition and therefore *permits wild oscillatory defects* — then vortex stretching is depleted. The stretching eigenvalue is recast as a singular-integral commutator and controlled by a localised Coifman–Rochberg–Weiss estimate with dyadic BMO tail bounds; the stretching potential vanishes as a logarithmic envelope on shrinking super-level sets, forcing $|\omega|$ into a subcritical Lorentz–Zygmund space. **This hypothesis is strictly weaker than Hölder continuity of $\xi$ and is numerically measurable**, which makes it the sharpest available target for the campaign (§4.3).

**A17 (obstructions).**
- *Tao 2016 (averaged Navier–Stokes)*: an averaged equation with the same energy identity and scaling blows up in finite time. Consequence: **no proof using only the energy identity, scaling, and abstract harmonic analysis can succeed.** Note this is *consistent* with — indeed points toward — a kernel-specific alignment route (A7), since Tao's averaging destroys exactly the kernel cancellation structure that A7 exploits.
- *Elgindi 2021*: finite-time blow-up for 3D Euler with $C^{1,\alpha}$ velocity.
- *Chen–Hou*: computer-assisted proof of finite-time blow-up for 2D Boussinesq and 3D axisymmetric Euler with smooth data **and boundary** [VERIFY: PNAS 2025]. Does not settle the periodic viscous problem, but it removes "blow-up is impossible for smooth data" as an intuition.

---

## §3. Part B — the no-go ledger

This section is the honest core of the document. Each proposition is a negative
result about the *method* originally proposed.

### Proposition B1 (a box-dimension hypothesis cannot improve the enstrophy bound)

Take $L=1$. Let $E(t)$ be any measurable set carrying the stretching, i.e.
$\int\omega\cdot S\omega \le \int_E |\omega|^2|S| + (\text{negligible})$.
Hölder with exponents $(3,6,2)$ — precisely $\tfrac{2}{6}+\tfrac16+\tfrac12=1$ — gives

$$\int_E |\omega|^2|S|\;\le\;\|\omega\|_{L^6}^2\,\|S\|_{L^6}\,|E|^{1/2}.$$

By Sobolev $\|\omega\|_{L^6}\le C\|\nabla\omega\|_{L^2}$, and by A5 with $p=6$,
$\|S\|_{L^6}\le C\|\omega\|_{L^6}\le C\|\nabla\omega\|_{L^2}$. Hence

$$\int_E|\omega|^2|S| \;\le\; C\,Y^{3/2}\,|E|^{1/2}. \tag{3.1}$$

Now impose the codimension hypothesis. The natural length built from the
available norms is $\ell=(\Omega/Y)^{1/2}$, and a set of box dimension $D$
resolved at scale $\ell$ inside a unit envelope has volume
$|E|\sim \ell^{\,3-D}$:

$$\textbf{(H}_D\textbf{)}\qquad |E(t)| \;\le\; C\,(\Omega/Y)^{(3-D)/2}. \tag{3.2}$$

Substituting (3.2) into (3.1) and then into the balance (2.1):

$$\Omega' \;\le\; C\,\Omega^{(3-D)/4}\,Y^{(3+D)/4} \;-\;\nu Y.$$

Young's inequality absorbs the first term into $\nu Y$ **iff** $(3+D)/4<1$, i.e.
$D<1$, and then with conjugate exponents $\tfrac{4}{3+D},\tfrac{4}{1-D}$:

$$\boxed{\;\Omega' \;\le\; -\tfrac{\nu}{2}Y \;+\; C\,\nu^{-\frac{3+D}{1-D}}\;\Omega^{\frac{3-D}{1-D}}\;}\tag{3.3}$$

| $D$ | exponent $\frac{3-D}{1-D}$ | $\nu$ power | verdict |
|---|---|---|---|
| $0$ | $3$ | $\nu^{-3}$ | **exactly Doering–Gibbon (A4)** |
| $1/2$ | $5$ | $\nu^{-7}$ | strictly worse |
| $3/4$ | $9$ | $\nu^{-15}$ | much worse |
| $D\to1^-$ | $\to\infty$ | $\to\infty$ | degenerate |
| $D\ge1$ | — | — | **Young fails; no bound at all** |

**Consequences.**

1. The $D=0$ case *reproduces the textbook bound exactly*, constants and
   $\nu$-power included. This is a correctness check on (3.3), and it is also the
   whole point: **the standard supercritical bound already behaves as though the
   stretching were maximally concentrated.** Codimension has nothing left to buy.
2. For every $D>0$ the hypothesis makes the bound *worse*, monotonically.
3. For $D\ge1$ — tubes, sheets, i.e. every scenario anyone actually worries about
   (Kerr anti-parallel tubes, Luo–Hou) — the route yields nothing whatsoever.

**Therefore the question "does concentration geometry enforce depletion of the
globally integrated vortex stretching?" has the answer "not as posed, and not for
any value of the codimension."** $\square$

### Proposition B2 (the concentration volume cancels at the dissipation scale)

Let $M=\|\omega\|_\infty$ and $\ell_\nu=(\nu/M)^{1/2}$, the scale at which
stretching and diffusion balance. On a concentration set $E$ resolved at
$\ell_\nu$:

- stretching: $|\omega|^2\alpha \le \rho_E M^3$ pointwise on $E$, total $\rho_E M^3|E|$;
- dissipation: $|\nabla\omega|\sim M/\ell_\nu$ on $E$, so $\nu|\nabla\omega|^2 \sim \nu M^2/\ell_\nu^2 = \nu M^2\cdot M/\nu = M^3$, total $\kappa_E M^3|E|$.

Both terms carry **the same power of $M$ and the same factor $|E|$.** The volume
cancels identically. What remains is the dimensionless ratio

$$\mathcal R_E(t)\;:=\;\frac{\int_E \omega\cdot S\omega\,dx}{\nu\int_E|\nabla\omega|^2dx}\;=\;\frac{\rho_E}{\kappa_E}\;=\;O(1).$$

**No geometric quantity with dimensions of a volume, a dimension, or a
codimension can decide this balance.** Only an $O(1)$ depletion constant can, and
the only published mechanism producing one is the directional cancellation in the
kernel of A7. $\square$

**Corollary B2.1 (on the originating diagnostic).** The globally integrated ratio
$\mathcal R(t)=\int|\omega|^2\xi^{\mathsf T}S\xi\big/\nu\int|\nabla\omega|^2$
reported in the motivating runs satisfies $\mathcal R>1 \iff \Omega'>0$, by (2.1).
It is therefore a restatement of "enstrophy is currently increasing" and carries
**no** information about singularity formation: it exceeds 1 in every transiently
amplifying flow, including flows known to be globally regular. The measured value
$\mathcal R\approx2.857$ is not evidence of anything. The set-restricted
$\mathcal R_E$ is the informative object, and it is not a restatement of the balance.

### Proposition B3 (the scaling gap in Sobolev language)

Enstrophy dissipation supplies $\int_0^T\!\|\nabla\omega\|_2^2 <\infty$, i.e.
control of $\omega$ in $L^2_tH^1_x$. Continuation via BKM (A6) requires
$\|\omega\|_{L^\infty}$, and in 3D $H^s\hookrightarrow L^\infty$ needs $s>3/2$.
The deficit is exactly one half derivative, uniformly in $\nu$. Every known
partial result closes a logarithmic fraction of this gap; A14 closes an algebraic
fraction. None closes it. $\square$

### Proposition B4 (dimension is neither necessary nor sufficient for thinness)

The hypothesis that actually appears in A13 is *transverse thinness* — along some
direction through each point, the set occupies a small fraction of a segment of
length $r$ — not dimension.

- **Not sufficient:** a fat Cantor set has positive Lebesgue measure and box
  dimension 3 in the product construction, yet is nowhere dense; conversely sets
  of small box-counting volume can be equidistributed at scale $r$ and fail 1D
  sparseness at that scale.
- **Not necessary:** a straight vortex tube of radius $a$ has $D=1$ — the *worst*
  case for B1 — but is perfectly 1D-sparse transversally at any $r\gg a$.

So the geometric object to measure is the sparseness ratio $\delta(r)$, **not**
$D$. This is why the measurement spec (`spec/gamma_measurement.md`) treats $D_q$
as a secondary diagnostic and $\delta(r)$ as the theorem-relevant primary. $\square$

---

## §4. Part C — conditional theorems

Each theorem below is a *correct reduction*. None is unconditional. The
hypotheses are stated in the form the numerics actually measure.

### Lemma C1 (transverse thinness $\Rightarrow$ 1D sparseness)

Let $A\subset\mathbb T^3$, $x_0\in A$, $r>0$. Say $A$ is *$\delta$-thin at
$(x_0,r)$* if there is a unit vector $n$ with
$\big|\{s\in(-r,r): x_0+sn\in A\}\big|_1 \le 2\delta r$.
Then $A$ is 1D $\delta$-sparse at $(x_0,r)$ in the sense of A13. *Proof:*
immediate from the definitions. $\square$

Trivial, but it fixes the object the pipeline estimates: for each sampled point of
the set, minimise the occupied fraction over a fixed set of directions.

### Lemma C2 ($\Gamma$-sets embed in $\omega$-sets)

$$E_\lambda(t)\;\subset\;\Omega_{\lambda'}(t),\qquad \lambda'=\big(\lambda/C_{\log}\big)^{1/2},$$
with $C_{\log}$ from (2.2).

*Proof.* On $E_\lambda$, $\lambda M^3 < \Gamma \le |\omega|^2\|S\|_\infty \le C_{\log}M|\omega|^2$,
so $|\omega|^2 > (\lambda/C_{\log})M^2$. $\square$

**Consequence.** Sparseness of $\omega$-super-level sets implies sparseness of
$\Gamma$-super-level sets (at a shifted threshold). Hence the $\Gamma$-version of
the sparseness criterion is a **strictly stronger hypothesis** than A13's, so
"$\Gamma$-sparseness $\Rightarrow$ regularity" does **not** follow from A13 and is
recorded as Conjecture D1, not as a theorem. This is the single most important
bookkeeping point in the document: the natural strengthening runs the wrong way.

### Theorem C3 (thinness route) — node `NS-G02`

*Suppose* there exist $\lambda\in(0,1)$, $\delta<\delta_0$ and $c>0$ such that for
all $t$ in a left-neighbourhood of $T^*$, every point of $\Omega_\lambda(t)$ is
$\delta$-thin at scale $r=c\,\ell_\nu(t)$, with $\delta_0$ and the admissible
range of $c$ as in A13. *Then* $T^*$ is not a singular time.

*Proof.* Lemma C1 converts thinness to 1D $\delta$-sparseness at the stated
scale; A13 applies verbatim. $\square$

**This is a reduction, not a new theorem.** Its content is entirely that
$\delta(r)$ at $r\sim\ell_\nu$ is the right thing to measure. The hypothesis is
node `NS-G01` and is unproved.

### Theorem C4 (direction-regularity route) — node `NS-G04`

*Suppose* that for $t$ near $T^*$: (i) $\xi\in C^\beta$ uniformly on an
$r$-neighbourhood of $E_\lambda(t)$ with $\beta\ge1/2$ and $r\gtrsim\ell_\nu$;
and (ii) $E_\lambda(t)$ carries all but an $\varepsilon$-fraction of the positive
stretching, $\int_{E_\lambda}\Gamma \ge (1-\varepsilon)\int\Gamma$ with
$\varepsilon$ small enough that the complement's contribution is absorbed by
$\tfrac{\nu}{2}Y$ via (3.1). *Then* $T^*$ is not a singular time.

*Proof.* (ii) reduces the balance (2.1) to the contribution of $E_\lambda$ up to a
term absorbed by dissipation; on $E_\lambda$, (i) is the hypothesis of A9. $\square$

**This is where "$\theta>\theta_c$" lives honestly: $\theta=\beta$ is the Hölder
exponent of the direction field, and $\theta_c=1/2$.** It is not a codimension and
there is no threshold in any codimension, by B1.

### Theorem C5 (log-modulus route) — node `NS-G06`

*Suppose* $\xi$ restricted to a neighbourhood of $E_\lambda(t)$ lies in
$\mathrm{bmo}_{1/|\log r|}$ uniformly on $[t_0,T^*)$, and the vorticity
concentration is of the $L^{3/2,\infty}$ type assumed in A16. *Then* stretching
is depleted in the sense of A16 and $T^*$ is not singular.

*Proof.* Reduction to A16. $\square$

**C5 is the weakest hypothesis of the three** — bmo with a logarithmic weight
admits oscillatory defects that Hölder-$1/2$ forbids — and it is the one the
campaign should target, because a *negative* measurement (the modulus failing to
decay logarithmically on the intense set) would be informative, whereas a
negative Hölder measurement is nearly guaranteed by intermittency and would tell
us little.

### Where these close NS-C01

C3, C4, C5 each close `NS-C01` **iff** their hypothesis is shown to follow from
the equations for all $t<T^*$. That derivation is `NS-G01`. Nothing in this
document advances it. Per A17, any such derivation must use the kernel structure
of A7; a proof relying only on energy, scaling, and abstract harmonic analysis is
excluded by Tao's counterexample.

---

## §5. Part D — conjectures

### Conjecture D1 ($\Gamma$-sparseness) — node `NS-G03`

Sparseness of the *stretching-intense* sets $E_\lambda(t)$ at scale
$c\,\ell_\nu(t)$ suffices for regularity.

**Missing lemma (L).** By Lemma C2 the implication runs the wrong way, so D1 does
not reduce to A13. What is needed: points of $\Omega_\lambda\setminus E_\lambda$
(high vorticity, currently low stretching) cannot become high-$\Gamma$ within a
time $\sim 1/M$. That is a bound on the *alignment rate* $D\alpha/Dt$ along
Lagrangian trajectories — restricted-Euler / Vieillefosse-type control of the
strain along particle paths, with the pressure Hessian retained rather than
dropped. **This is the extra ingredient.** It is not supplied here.

### Conjecture D2 (depletion exponent) — node `NS-G05`

With $\rho(\ell,t)=\int\Gamma[\bar u_\ell]\big/\int|\bar\omega_\ell|^2|\bar S_\ell|$
the coarse-grained depletion ratio, there is $\theta>0$ with
$\rho(\ell_\nu(t),t)\le C(\ell_\nu/L)^\theta$ uniformly in $t$.

**D2 does not close the enstrophy route even if true.** By B2 the balance is
decided by the $O(1)$ ratio $\mathcal R_E$, and a bound on $\rho$ at a single
scale does not bound $\mathcal R_E$. D2 is a *diagnostic prediction* — a
falsifiable statement about the flow's geometry that the campaign can test — and
it is recorded as such, not as a step in a proof.

---

## §6. Part N — what the numerics can and cannot do

### 6.1 The prior runs carry no evidential weight

The two motivating rounds (16 cases at $16^3$–$32^3$; 24 adversarial cases at
$24^3$–$64^3$) are recorded as `NS-X01: SMOKE_TEST_ONLY`.

- **Peak vorticity growth $1.537\times$ is not large.** Undriven Taylor–Green
  amplifies more than this before its enstrophy peak. A growth factor below ~10
  at $64^3$ distinguishes nothing.
- **$\mathcal R>1$ transiently is vacuous** (Corollary B2.1).
- **$\mathrm{FTLE}\approx9.93$ is not a singularity diagnostic.** Lagrangian chaos
  is known to coexist with global regularity; FTLE is a *search coordinate*, not
  a criterion. This was correctly identified in the originating write-up and is
  restated here only to fix its status.
- **$\zeta_p$ at $64^3$ is noise.** There is no inertial range at $64^3$; the
  fitted exponents measure the dissipation range and the forcing of the initial
  spectrum, not intermittency.
- **Zero surviving candidates is the expected outcome of a search with no power.**
  Reporting it as a result invites the reader to treat it as evidence for
  regularity. It is not.

### 6.2 The multifractal branch is downgraded, not abandoned

`NS-F01` (structure-function exponents $\to$ regularity) has **no theorem behind
it** and is marked `DIAGNOSTIC_ONLY`. It is replaced by `NS-F02`: the
Cheskidov–Shvydkoy number $\sup_q 2^{-q}\|\Delta_q u\|_\infty/\nu$ over the
dissipation range, which is connected to an actual regularity criterion (A15).
Measure that instead of $\zeta_p$.

### 6.3 What a resolution ladder can and cannot establish

At the resolutions reachable on four GPUs, with $k_{\max}\eta\ge2$ enforced,
$512^3$ forces $\nu\sim2\text{–}3\times10^{-3}$ and $\mathrm{Re}_\lambda\sim
100\text{–}200$. **Viscous Navier–Stokes at that viscosity will not blow up.** A
campaign whose primary endpoint is "find a surviving blow-up candidate" is
therefore designed to return null regardless of the truth of the conjecture, and
its null result will again be uninformative.

**The campaign must therefore be re-scoped**, and is, in
`spec/campaign_4gpu.md`:

- **Primary endpoint:** measurement of $\delta(r/\ell_\nu)$, $\beta$, the
  mean-oscillation modulus $\mu_\xi(r)$, $\mathcal R_E$, and $\rho(\ell_\nu)$,
  together with their trends as $\nu$ decreases and their stability under
  refinement. These are the hypotheses of C3, C4, C5.
- **Secondary endpoint:** the blow-up candidate search, with a preregistered
  expected count of **zero**, retained so that the gates are exercised and a
  surprise would be caught.

A measured logarithmic decay of $\mu_\xi(r)$ on the $\Gamma$-set that strengthens
as $\nu\to0$ is evidence *for* the depletion mechanism of A16. A measured failure
of that decay, stable under refinement, would localise where a singularity would
have to hide. Either outcome is publishable; "no candidate found" is not.

---

## §7. Honest conclusion

The question as originally posed — *does concentration geometry, in the sense of
a codimension, enforce depletion of integrated vortex stretching?* — is answered
**no** by Proposition B1, for every value of the codimension, with the $D=0$ case
recovering the textbook bound exactly. That is a real if unwelcome result, and it
is the main new content of this document.

The question survives in a changed form. Replace *codimension* by *transverse
thinness* and *direction regularity*, and there are three published sufficient
conditions (A13, A9, A16) to which regularity genuinely reduces — Theorems C3, C4,
C5. All three hypotheses are measurable. None is proved, and the derivation of any
of them from the equations is `NS-G01`, which is the Millennium problem.

What can honestly be built, and is built here, is: the reduction, the ledger of
what cannot work and why, a measurement pipeline for the surviving hypotheses,
and a preregistered campaign that can falsify them. That is a proof *search*
architecture. It is not a proof, and the distinction is the point.
