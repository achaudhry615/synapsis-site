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

## §0.1 Which Clay problem this addresses

Fefferman's official statement poses four assertions. With $\nu>0$, $n=3$ and
$f\equiv0$:

| | statement | domain |
|---|---|---|
| (A) | existence and smoothness | $\mathbb R^3$ |
| **(B)** | **existence and smoothness** | **$\mathbb R^3/\mathbb Z^3$ (periodic)** |
| (C) | breakdown | $\mathbb R^3$ |
| **(D)** | **breakdown** | **$\mathbb R^3/\mathbb Z^3$ (periodic)** |

This document works on the **periodic, unforced** case, so the conditional
theorems of Part C are statements aimed at **(B)**, and the numerical search
described in Part N is aimed at **(D)**. A proof of any one of the four suffices
for the prize; nothing here proves any of them.

Two points from the official statement are worth carrying explicitly, because
both bear on claims made elsewhere in this document.

**On CKN and dimension (A10).** The partial-regularity theorem is stated for the
*parabolic* analogue of Hausdorff dimension, built from cylinders
$Q_r=B_r\times I_r$ with $B_r\subset\mathbb R^3$ a ball of radius $r$ and
$I_r\subset\mathbb R$ an interval of length $r^2$; the conclusion is
$\mathcal P^1(E)=0$ for the singular set $E$, and in particular $E$ cannot
contain a space-time curve $\{x=\phi(t)\}$. Fefferman calls it "the best partial
regularity theorem known so far" and notes "it appears to be very hard to go
further." This confirms the caveat in A10: $\mathcal P^1$ is a parabolic measure
of the space-time singular set, **not** a box dimension of a pre-singular
super-level set, and the two must not be conflated — which is precisely the
conflation Proposition B4 warns against.

**On method.** The statement closes: *"Standard methods from PDE appear
inadequate to settle the problem. Instead, we probably need some deep, new
ideas."* Read alongside Tao's averaged-Navier–Stokes counterexample (A17), which
makes the same point in a sharper and more specific form — no argument using only
the energy identity, scaling and abstract harmonic analysis can work — this is
the reason Part B spends its effort *ruling routes out*. Anything that survives
must use structure the averaged equation destroys.

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

**A13 (Grujić sparseness criterion).** [VERIFIED bibliographically: Grujić, *Nonlinearity* **26** (2013) 289, arXiv:1111.0217.] A local anisotropic geometric-measure condition preventing finite-time singularity formation, via *1D $\delta$-sparseness* of the regions of intense fluid activity.

**Two corrections to how this was stated in an earlier draft of this document, both found during citation verification:**

1. **The super-level sets are those of the positive and negative parts of the vorticity COMPONENTS $\omega_i^{\pm}$ — six sets — not of $|\omega|$.** The framework is built on "the scale of sparseness of the super-level sets of the positive and negative parts of the vorticity components". A sparseness statement about $\{|\omega|>\lambda M\}$ is a *different* hypothesis and is not what the theorem assumes.
2. **The scale is comparable to the radius of spatial analyticity** (measured in $L^\infty$) — the natural dissipation scale produced by analytic smoothing in the viscous case — *not* $\ell_\nu$ by definition, as written before. The two are of the same order, but the identification needs an argument.

Mechanism: the harmonic-measure majorisation principle (log-convexity of the modulus of an analytic function), applied on a disc of radius comparable to the analyticity radius. If the intense-vorticity regions are sparse enough the associated harmonic measure is small enough to forbid further growth of the $L^\infty$ norm — a "self-improving" bound on the sup-norm.

**Independent second proof.** Albritton–Bradshaw, *Nonlinearity* **35** (2022) 2858 (arXiv:2110.02187) give "a simple proof that sufficiently sparse Navier–Stokes solutions do not develop singularities", explicitly as an alternative to Grujić's analyticity/harmonic-measure argument. The *sufficiency* direction — sparse $\Rightarrow$ regular — is therefore proved twice by different methods. That is what Theorem C3 reduces to.

[STILL UNVERIFIED: the threshold value of $\delta_0$ and the constant relating the sparseness scale to the analyticity radius; primary sources were egress-blocked. `H2`'s threshold remains UNSET for this reason.]

**Numerical consequence of correction (1).** Measured at $64^3$, $\lambda=1/2$, $r=16$ cells: for the smoothed vortex sheet, $\delta_{p95}=0.061$ on the $|\omega|$-set but $0.152$ on the worst of the six component sets — a factor **2.5 larger**. The wrong definition reported the configuration as *thinner than it is*, i.e. it was biased **toward falsely confirming** the C3 hypothesis. `gamma.sparseness_grujic` now computes all six component sets and reports the worst, since the criterion requires every one of them to be thin.

**A14 (algebraic reduction of the scaling gap) — CONTESTED, not relied upon.** Bradshaw–Farhat–Grujić, *ARMA* **231** (2019) 1983–2005 (arXiv:1704.05546): within the sparseness framework, an a priori bound algebraically better than the energy-level bound while keeping the regularity criterion at the classical level — the first non-logarithmic reduction of the scaling gap since the 1960s.

**However:** Albritton–Bradshaw, *Nonlinearity* **35** (2022) 2858, give as their second goal that they "analyze the claims that a priori estimates on the sparseness of the vorticity and higher velocity derivatives reduce the 'scaling gap' in the regularity problem" — a critical examination of Bradshaw–Farhat–Grujić (2019) and Grujić–Xu, **co-authored by one of the authors of the 2019 paper**. We could not obtain it (egress-blocked) and do not know the verdict. **A14 is therefore recorded as contested and is used nowhere in Parts B, C or D.** Nothing in this document depends on the scaling gap having been reduced.

This sharpens the framing. The *sufficiency* direction (sparse $\Rightarrow$ regular, A13) is solid. What is contested is the *a priori* direction — whether Navier–Stokes actually produces sparseness at the required scale. **That a priori direction is precisely `NS-G01`.** The literature's open question and this document's open node are the same question.

Related: Farhat–Grujić–Leitmeyer (volumetric sparseness, $B^{-1}_{\infty,\infty}$); Grujić–Xu (asymptotic criticality). Both inherit the caveat.

**A15 (Cheskidov–Shvydkoy).** Regularity criteria in $B^{-1}_{\infty,\infty}$ [VERIFY]. This is the **only** rigorous bridge from spectral/band-limited data to a regularity theorem, and it is the reason the multifractal branch of the original program is downgraded (§6.2).

**A16 (log-weighted bmo depletion).** [VERIFIED bibliographically; hypotheses from the abstract: Grujić, arXiv:2607.08866, 13 July 2026.] For critical point singularities exhibiting $L^{3/2,\infty}$ spatial concentration of vorticity: if the vorticity *direction* lies locally in a logarithmically weighted BMO space $\mathrm{bmo}_{1/|\log r|}$ — a space failing the Dini condition and therefore permitting wild oscillatory defects — the nonlinear vortex stretching is fundamentally depleted. The chain is:

1. isolate a *unidirectional geometric cancellation* and recast the stretching eigenvalue as a **singular-integral commutator**;
2. control it by a localised **Coifman–Rochberg–Weiss** estimate with dyadic BMO tail bounds, so the stretching potential vanishes as a logarithmic envelope on shrinking super-level sets;
3. an interpolated **De Giorgi** energy method forces $|\omega|$ into a subcritical **Lorentz–Zygmund** space;
4. the logarithmic gain is **transferred to the velocity field, forcing the geometric scale of local 1D sparseness below the uniform radius of spatial analyticity**, which averts blow-up.

**Step 4 is structurally important and corrects how this document previously organised Part C.** A16 does not terminate in its own independent regularity mechanism: it terminates by *establishing the hypothesis of A13* — sparseness below the analyticity radius — and then invokes it. So the three routes are **not parallel**:

$$\text{C5 hypothesis} \;\Longrightarrow\; \text{C3 hypothesis} \;\Longrightarrow\; \text{regularity.}$$

C5 is an *upstream* result whose value is that it derives C3's geometric hypothesis from a hypothesis on the direction field. That is exactly the shape `NS-G01` asks for — except that it replaces one unproved hypothesis with another, rather than deriving either from the equations. It is nonetheless the most substantive progress in the literature toward NS-G01 that we located, and it independently justifies measuring $\delta$ at $r\sim R_{\rm an}$: that is precisely the quantity step 4 drives below threshold.

This also means the C5 hypothesis is the weakest of the three *and* implies the others' geometric conclusion, which makes $\mu_\xi(r)\log(L/r)$ the single most informative thing the campaign can measure.

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

### Proposition B3 (the scaling gap) — **corrected**

**An earlier version of this proposition was wrong** and is corrected here. It
asserted that "enstrophy dissipation supplies $\int_0^T\|\nabla\omega\|_2^2<\infty$,
i.e. control of $\omega$ in $L^2_tH^1_x$." That is not an a priori statement:
$\int Y\,dt$ is exactly what is *not* controlled — an enstrophy bound is the
thing being sought. What the energy identity actually gives is

$$\int_0^T\|\nabla u\|_2^2\,dt=\int_0^T\|\omega\|_2^2\,dt\le \frac{E_0}{\nu},
\qquad\text{i.e. }\omega\in L^2_tL^2_x\text{ only.}$$

Continuation via BKM (A6) requires $\|\omega\|_{L^\infty}$, and in 3D
$H^s\hookrightarrow L^\infty$ needs $s>3/2$. **The half-derivative conclusion
stands; the route to it stated before did not.** $\square$

The cleanest form of the same deficit is a time integral (verified in
`code/scripts/verify_ledger.py`). Analyticity gives $|\omega|\ge M/2$ on a ball of
radius $\sim\ell_\nu$, so $\|\omega\|_2^2\gtrsim M^2\ell_\nu^3=\nu^{3/2}M^{1/2}$,
and therefore

$$\boxed{\;\int_0^{T^*}\!\!M(t)^{1/2}dt\;\lesssim\;E_0\,\nu^{-5/2}\quad\text{a priori},
\qquad\text{while BKM needs}\quad\int_0^{T^*}\!\!M(t)\,dt<\infty.\;}$$

Under $u_\lambda(x,t)=\lambda u(\lambda x,\lambda^2t)$: $M^{1/2}dt\sim\lambda^{-1}$
(**supercritical**) while $M\,dt\sim\lambda^{0}$ (**critical**). **The entire
problem is one factor of $M^{1/2}$ inside a time integral.** Every route in
Part C and Part E is an attempt to buy that factor.

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
all $t$ in a left-neighbourhood of $T^*$, **every one of the six sets
$\{\omega_i^{\pm}>\lambda\|\omega_i^{\pm}\|_\infty\}$, $i=1,2,3$,** has all of its
points $\delta$-thin at scale $r=c\,R_{\rm an}(t)$, where $R_{\rm an}$ is the
uniform lower bound on the radius of spatial analyticity and $\delta_0,c$ are as
in A13. *Then* $T^*$ is not a singular time.

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

**C5 is the weakest hypothesis of the three, and it implies C3's geometric
conclusion** (A16 step 4). bmo with a logarithmic weight admits oscillatory
defects that Hölder-$1/2$ forbids, so C5's hypothesis is strictly weaker than
C4's; and because A16 terminates by forcing the sparseness scale below the
analyticity radius, C5 does not merely run parallel to C3 but *supplies its
hypothesis*. The implication structure is

$$\text{C5 hyp.}\;\Rightarrow\;\text{C3 hyp.}\;\Rightarrow\;\text{regularity},
\qquad \text{C4 hyp.}\;\Rightarrow\;\text{regularity (independently).}$$

The campaign should therefore target $\mu_\xi(r)\log(L/r)$ above all: a *negative*
measurement (the modulus failing to decay) is informative, whereas a negative
Hölder measurement is nearly guaranteed by intermittency and would tell us little.

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

## §5.5 Part E — retrosynthetic ledger for `NS-G01`

Working backwards from regularity through candidate disconnections. **Nothing
here is progress toward a proof.** The product is *death certificates*: for each
sub-lemma, the exact step that fails and the obstruction that kills it. Claims
marked ✔ are machine-checked in `code/scripts/verify_ledger.py` (11/11).

### E1. Scaling audit — applied first, because it is the cheapest filter

Under $u_\lambda(x,t)=\lambda u(\lambda x,\lambda^2 t)$, a hypothesis $Q$ with
$Q(u_\lambda)=\lambda^aQ(u)$ is *supercritical* if $a<0$, *critical* if $a=0$,
*subcritical* if $a>0$.

| hypothesis | $a$ | class | verdict |
|---|---|---|---|
| energy $\|u\|_2^2$ | $-1$ | supercritical | the only a priori input |
| $\int M^{1/2}dt$ | $-1$ | supercritical | a priori ✔ |
| $\int M\,dt$ (BKM) | $0$ | critical | the target ✔ |
| **C4**: $[\xi]_{C^\beta}\le C$ | $+\beta$ | **subcritical for every $\beta>0$** | **not a derivable target** |
| **C3**: $\delta$-sparse at $c\,R_{\rm an}$ | $0$ | critical | admissible |
| **C5**: $\xi\in\mathrm{bmo}_{1/|\log r|}$ | $0$, with a **log gain** | critical | admissible |

Two consequences the table forces:

1. **C4 is dead as a *derivation* target on scaling grounds alone.** Its Hölder
   constant carries a fixed length $L_c=C^{-2}$; at the dissipation scale the
   hypothesis demands $|\Delta\xi|\lesssim(\ell_\nu/L_c)^{1/2}\to0$ — a bound that
   *strengthens under zoom-in*, which no supercritical input can produce. C4
   remains perfectly valid as a **reduction**; it is simply not something one
   derives.
2. **Only C3 and C5 are scale-admissible**, and C5 is admissible precisely
   because it asks for a *logarithmic* rather than algebraic gain — the smallest
   currency that pays the $M^{1/2}$ deficit at the Leray rate
   $M\gtrsim(T^*-t)^{-1}$.

### E2. Why $\beta=1/2$, settled ✔

With $|D|\le C|y|^\beta$ in Constantin's kernel (A7), interpolation using
$3/2-\beta$ derivatives and Young give

$$\Omega'\le-\tfrac{\nu}{2}Y+C(\beta,\nu)\,\Omega^{\frac{3+2\beta}{1+2\beta}},
\qquad \beta=0\to3,\;\; \tfrac14\to\tfrac73,\;\; \tfrac12\to2,\;\; 1\to\tfrac53 .$$

Gronwall closes against the a priori $\int\Omega\,dt\le E_0/2\nu$ **iff the
exponent is $\le2$, i.e. iff $\beta\ge1/2$.** That is the reason for the
threshold, and it is a statement about the *consequence*, not the hypothesis —
consistent with E1(1).

### E3. The alignment ODE, and why "Lemma L" dies pointwise ✔

Exact Lagrangian identity (Euler part; $W\xi=\tfrac12\omega\times\xi=0$ kills the
rotation term):

$$\frac{D\alpha}{Dt}=-\alpha^2+|P^{\perp}_\xi S\xi|^2-\xi\cdot H\xi+\nu(\cdots),
\qquad H=\nabla^2p .$$

Write $F=|P^\perp S\xi|^2-\xi\cdot H\xi$. Pure Riccati $D\alpha/Dt\le-\alpha^2$
would give regularity, so everything is the sign and size of $F-\alpha^2$.

**Death certificate (pointwise version).** For a **Burgers vortex** —
$S=\mathrm{diag}(-\gamma/2,-\gamma/2,\gamma)$, $\xi=e_z$, $p=-\tfrac12x^\top S^2x$ —
one has $P^\perp S\xi=0$ and $\xi\cdot H\xi=-\gamma^2=-\alpha^2$, hence

$$F-\alpha^2=0\quad\text{identically.}$$

**The pressure Hessian exactly cancels the Riccati damping in the canonical
stretched vortex** ✔. Separately, at $t=0$ one may choose smooth data with
$\xi_0$ not an eigenvector of $S_0$ at the maximum, giving $F-\alpha^2\ge c'M^2$.
So any scale-invariant pointwise "$F-\alpha^2\le\varepsilon M^2$" is **false for
generic smooth data**. Only an asymptotic ($t\to T^*$) version survives, and it
amounts to asserting that the blow-up profile self-organises to be
Burgers-trivial at the singular point — a statement *about* the singularity, not
derivable from a priori bounds. **OPEN-HARD.**

**Corollary (why restricted Euler is not evidence).** Replacing $H$ by its
isotropic part gives the Vieillefosse–Cantwell ODE, which blows up for almost all
data. Since true NS is regular for the same data on short times, that blow-up is
an artefact of discarding the deviatoric Hessian $R\otimes R(|\omega|^2/2-|S|^2)$
— a sign-indefinite Riesz transform of a quadratic, with no Constantin-type
diagonal cancellation.

### E4. The commutator identity, and the sharp form of "coherence ⇒ regularity"

Since Constantin's kernel vanishes on the diagonal ($\det(\hat y,\xi,\xi)=0$),
splitting $\xi_k(x+y)|\omega(x+y)|$ about $\xi_k(x)$ annihilates the leading piece
and leaves a **pure commutator**

$$\alpha=\xi_i\xi_j\,[T_{ijk},\xi_k](|\omega|),\qquad
\|\alpha\|_{L^p}\le C_p\|\xi\|_{\rm BMO}\|\omega\|_{L^p}
\;\;\text{(Coifman–Rochberg–Weiss)} .$$

Because $|\xi|=1$ forces $\|\xi\|_{\rm BMO}\le2$ always, the content is entirely
in *smallness or decay* of the local mean oscillation — which is exactly the
quantity $\mu_\xi(r)$ that the campaign measures, and exactly C5's hypothesis.
Running this through the $L^{3/2}$ vorticity balance gives regularity whenever a
dimensionless **geometric Reynolds number**

$$\mathrm{Re}_\xi:=\|\xi\|_{\rm BMO}\,\|\omega\|_{L^{3/2}}/\nu$$

is small. This is the classical critical small-data theorem with the constant
replaced by the direction field's oscillation. **As a reduction it is sharp and
worth stating; as a derivation it dies** — analyticity gives
$\mathrm{osc}_{\ell_\nu}\xi=O(1)$, not $o(1)$, and the natural energy
$\iint|\omega|^2|\nabla\xi|^2\le\int Y$ is not a priori bounded (B3).

### E5. Why modulus-of-continuity propagation does not transfer

The direction field obeys a harmonic-map heat flow with transport, cross-diffusion
and a **zeroth-order forcing**:

$$\partial_t\xi+u\cdot\nabla\xi=P^\perp_\xi S\xi+\nu(\Delta\xi+|\nabla\xi|^2\xi)
+2\nu(\nabla\log|\omega|\cdot\nabla)\xi .$$

Kiselev–Nazarov–Volberg's critical-SQG machinery needs (a) pure transport of the
quantity carrying the modulus, (b) a maximum principle, (c) *nonlocal*
dissipation of the same order as the velocity gain. **All three fail here.** At a
breakthrough pair at separation $r$: forcing $\lesssim M\,\omega(r)$, local
dissipation $\lesssim\nu\omega(r)/r^2$. Dissipation wins **iff $r\lesssim\ell_\nu$**
(and only up to $\ell_\nu/|\log\ell_\nu|^{1/2}$ for a log modulus). Above that
scale the equation is effectively inviscid and moduli grow like
$\exp\int\|S\|_\infty$ — BKM, circular. Below it, analyticity already gives
$\mathrm{osc}_r\xi\lesssim r/\ell_\nu$, which is $O(1)$ at $r=\ell_\nu$: **the
propagated modulus can be no better than its seed, and the seed is $O(1)$.** This
is Proposition B2 in harmonic-map form.

**Independent confirmation, and a qualification of this certificate.** Grujić,
arXiv:2609.05720 (4 Sept 2026), isolates the same geometric PDE — "the Harmonic
Map Heat Flow into the sphere supplemented with the fluid transport,
cross-diffusion and tangential strain" — for exactly the critical point
singularity $|\omega|\sim|x|^{-2}$ in $L^{3/2,\infty}$, and asks precisely the E5
question: *whether the NSE mechanics can propagate* the logarithmic decay
$\xi\in\mathrm{bmo}_{1/|\log r|}$ that the companion paper (A16) shows prevents
blow-up. That the equation above was derived here independently and agrees is
reassuring; that the question is posed as open in Sept 2026 is the relevant
status.

Two of its stated observations bear directly on the certificate above and
**partially undercut it**: (i) the $O(|x|^{-2})$ concentration, factored out of
the viscous cross-diffusion, generates an **outward radial drift
$4\nu x/|x|^2$ at the core**; (ii) the HMHF nonlinearity is **harmless for
$\tfrac12|\xi-e|^2$**. Neither is accounted for in the balance above, which
treated the dissipation as purely local and the nonlinearity as a liability. An
outward drift transports oscillation *away* from the core, which is a mechanism
of exactly the kind the certificate claims does not exist. **E5 is therefore
downgraded from DEAD to DEAD-FOR-GENERIC-TRANSFER: the KNV machinery does not
transfer as such, but the critical-concentration setting has structure our
argument does not cover, and the question is live in the literature.** We have
not read the paper body (arXiv is egress-blocked).

Sharper still: a scale-invariant profile $|\omega|=|x|^{-2}\Phi(\hat x)$,
$\xi=\Xi(\hat x)$ has $\mathrm{osc}_{B_r}\xi$ *independent of $r$*. So **C5's
hypothesis is precisely "the direction field has no scale-invariant angular
structure at the singular point"** — and deriving it means proving that.
Elgindi's $C^{1,\alpha}$ Euler blow-up has $\xi=\pm e_\theta$, $O(1)$ oscillation
on every ball meeting the axis, at every scale.

### E6. The Tao filter

| route | NS-specific structure used | survives averaging? |
|---|---|---|
| B1/B2 codimension | none (Hölder, Sobolev, CZ only) | **fails** — consistent with its death |
| a priori sparseness from energy-level bounds | none | **fails** — so Tao's blow-up is non-sparse, and no derivation from shared bounds can work |
| E3 alignment ODE | $DA/Dt=-A^2-\nabla^2p$, $W\xi=0$ | passes; dies elsewhere |
| E4 commutator | $\det(\hat y,\xi,\xi)=0$ — parallel vorticity induces no axial strain | passes; dies elsewhere |
| E5 modulus propagation | transport–stretching form | passes; dies elsewhere |

**The picture is consistent: every route that fails the filter is already dead by
exponent count, and every route that passes it dies at the same $O(1)$ unit
problem at scale $\ell_\nu$ and time $1/M$.** The kernel structure is *necessary*;
nothing in this ledger makes it *sufficient*.

### E7. Summary and what is *not* established

| sub-lemma | status | dies at |
|---|---|---|
| propagation of sparseness | **DEAD** | false for pure diffusion; thinness at $\delta cR_{\rm an}$ demands $\lvert S\rvert/\lvert\omega\rvert\ge(\delta c)^{-2}$, excluding tubes |
| a priori sparseness closing the gap | **DEAD** as a closable route | Tao filter (modulo a flagged transfer step) |
| "Lemma L" pointwise | **DEAD** ✔ | Burgers degeneracy + $t=0$ counterexample |
| "Lemma L" asymptotic | **OPEN-HARD** | no mechanism; $\xi\cdot H\xi$ sign-indefinite, $O(M^2)$ |
| restricted Euler | **DEAD** | Vieillefosse blow-up is a known artefact |
| KNV modulus propagation | **DEAD for generic transfer**; live in the critical setting | order-0 forcing vs local order-2 dissipation; but see the outward radial drift $4\nu x/|x|^2$ in arXiv:2609.05720, which our balance omits |
| C4 as a derivation target | **DEAD** | subcritical (E1) |
| commutator route | **DEAD as derivation**, valid as reduction | $\mathrm{osc}_{\ell_\nu}\xi=O(1)$ |
| Γ-criterion as a theorem | **OPEN-HARD** | $E_\lambda=\emptyset$ gives only $M'\le\lambda M^2$, blow-up-consistent |
| Lagrangian/Cauchy representation | **DEAD** | needs $\exp\int\|\nabla u\|_\infty$ — circular |

**Not established.** No route was found. Three items are OPEN-HARD and none is
recommended. The literature items were read only through search summaries
(arXiv is egress-blocked): in particular arXiv:2609.05720 (Sept 2026), on decay
of local mean oscillations of the vorticity direction, **is exactly the E5
sub-lemma and may change its status** — our certificate is derived from the
equation, not from that paper. The claim that Tao's averaged equation supports
the same analyticity local theory is plausible but unchecked. The non-vacuity
window for C3 (whether the analyticity radius already forces a core too fat to be
$\delta_0$-sparse) is a finite computation with two published constants that we
could not obtain, and it is the single most valuable next calculation: it would
determine whether C3 is a genuine criterion or is vacuous for tube-like profiles.

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

**A specific concern about the strongest prior candidate.** The originating run's
best case, ADV023, came from the `multiscale_random` family. Running that same
family here at $64^3$ with $\nu=0.02$ gives $k_{\max}\eta=1.33$, which **fails**
the resolution gate `HIGH_K_TAIL_RESOLVED` — the field is under-resolved from the
first step. A broadband initial spectrum puts energy at the grid scale by
construction, so it is the family most likely to be under-resolved at a given
$(N,\nu)$, and it is precisely the family whose apparent signal was strongest.

This can be made quantitative. At fixed initial enstrophy $\Omega_0$,
$\varepsilon=2\nu\Omega_0$ so $\eta=(\nu^2/2\Omega_0)^{1/4}$, and the viscosity
needed for the *initial condition* to satisfy $k_{\max}\eta\ge2$ is
$\nu_{\rm req}=(2/k_{\max})^2\sqrt{2\Omega_0}$. Evaluated on our seeded
broadband IC (`diagnostics.ic_resolution_check`):

| $N$ | 16 | 24 | 32 | 48 | 64 | 128 | 256 |
|---|---|---|---|---|---|---|---|
| $\nu_{\rm req}$ | 0.295 | 0.196 | 0.134 | 0.093 | **0.068** | 0.033 | 0.017 |

($\nu_{\rm req}$ falls more slowly than $k_{\max}^{-2}$ because a $-5/3$ initial
spectrum also carries more enstrophy as more high-$k$ modes are admitted.)

The prior rounds ran this family at $16^3$–$64^3$. Unless they used
$\nu\gtrsim0.07$ — an order of magnitude above typical DNS values, and a flow so
viscous that nothing interesting can happen — **the multiscale family was
under-resolved at every rung of that ladder**, most severely at the $16^3$–$24^3$
discovery stage where the candidate was selected.

We do not know the viscosity actually used and so do not assert this. But it is
exactly what a preregistered resolution gate exists to catch, and it should be
checked before $1.537\times$ growth or $\mathcal R\approx2.857$ from that case is
given any weight. Under-resolution inflates precisely these diagnostics: energy
piled at the grid scale raises $\|\nabla\omega\|$ and the grid-maximum
$\|\omega\|_\infty$ without corresponding physical structure.

Every other family we tested is comfortably resolved at $64^3$ for
$\nu\in\{0.02,0.05\}$ — the antiparallel tubes need only $\nu_{\rm req}=4\times10^{-4}$.
The defect is specific to broadband initial spectra, which is why a search
distribution built on them is the wrong instrument, independently of how
"adversarial" it looks.

### 6.2 The multifractal branch is downgraded, not abandoned

`NS-F01` (structure-function exponents $\to$ regularity) has **no theorem behind
it** and is marked `DIAGNOSTIC_ONLY`. It is replaced by `NS-F02`: the
Cheskidov–Shvydkoy number $\sup_q 2^{-q}\|\Delta_q u\|_\infty/\nu$ over the
dissipation range, which is connected to an actual regularity criterion (A15).
Measure that instead of $\zeta_p$.

### 6.3 What a resolution ladder can and cannot establish

**This subsection previously asserted that the geometric hypotheses are
unmeasurable except at $k_{\max}\eta\sim10$–$20$, and that the geometry
measurement and the blow-up search therefore compete for the same grid by an
order of magnitude in $\nu$. That assertion rested on evaluating sparseness at
$\ell_\nu$, which is not the scale A13 states. It is withdrawn.** The corrected
account is in `spec/gamma_measurement.md` §4; the short version:

- A13's scale is the radius of spatial analyticity $R_{\rm an}$, measurable from
  the exponential decay $E(k)\sim e^{-2R_{\rm an}k}$.
- Measured at $64^3$, $R_{\rm an}\approx4$–$7\,\ell_\nu$ and
  $R_{\rm an}/\Delta x\approx3$–$12$. The theorem's scale is resolved at standard
  resolution; $\ell_\nu$ is not.
- With the six component sets (§A13 correction 1) and this scale, the first
  informative measurement is $\delta_{\rm worst,p95}\approx0.28$–$0.32$ for
  antiparallel tubes, against the $1.000$ — i.e. no information — that the
  earlier setup returned.

What remains true, and is the real constraint:

At the resolutions reachable on four GPUs with $k_{\max}\eta\ge2$ enforced,
$512^3$ admits $\nu\approx9\times10^{-5}$ at unit volume-averaged enstrophy
(`code/scripts/resolution_budget.py`). **Viscous Navier–Stokes at a viscosity a
grid of that size can resolve will not blow up.** So a campaign whose primary
endpoint is "find a surviving candidate" is designed to return null whatever the
truth, and its null result is uninformative — the defect of both prior gauntlet
rounds. The endpoints must be reordered:

| endpoint | content | expected result |
|---|---|---|
| **primary** | $\delta_{\rm worst}(r/R_{\rm an})$, $\beta$, $\mu_\xi\log(L/r)$, $\mathcal R_E$, and their trends as $\nu\downarrow$ and under refinement | informative either way |
| **secondary** | blow-up candidate search through the gate ladder | **preregistered expected count: 0** |

The secondary endpoint is retained only so the gates are exercised and a surprise
would be caught.

A measured $\delta$ that stays below the A13 threshold and strengthens as
$\nu\to0$ is evidence for the depletion mechanism. A measured failure, stable
under refinement, localises where a singularity would have to hide. Either is
publishable; "no candidate found" is not.

**Note on this correction.** The withdrawn claim was not a slip of exposition: it
was a wrong scale propagated into a resolution budget, a campaign design and a
stated conclusion, and it was caught only by verifying the citation it came from.
That is the argument for `prereg/FREEZE_BLOCKERS.md` being worked through before
any GPU time is spent.

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
