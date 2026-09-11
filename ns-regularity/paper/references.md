# References

**Verification status.** arXiv, nature.com and par.nsf.gov were unreachable from
the runtime in which this document was drafted (egress-blocked). Entries marked
**[VERIFIED]** were confirmed; **[VERIFY]** entries
are stated from secondary sources (search summaries, abstracts) and their exact
hypotheses, exponents and constants **must be checked against the published
papers before any conditional theorem in `ns_gamma_codimension.md` is relied
upon**. This is not a formality: Theorems C3 and C5 are reductions *to* these
statements, so an error in the hypothesis of A13 or A16 propagates directly into
what we claim to have reduced.

## Primary source obtained in full

- **[VERIFIED — read directly]** Fefferman, C. L. *Existence and smoothness of the
  Navier–Stokes equation.* Official Clay Mathematics Institute problem statement
  (with errata). Supplied to this project as a PDF and read in full. Fixes the
  statements (A)–(D), the growth conditions (4)–(11), the weak-solution
  formulation (12)–(13), the CKN theorem in its parabolic form, and the
  bibliography below. Its closing assessment — "Standard methods from PDE appear
  inadequate to settle the problem. Instead, we probably need some deep, new
  ideas." — is the reason this project's Part B is devoted to ruling routes out.

## Foundational

- **[VERIFIED]** Leray, J. (1934). Sur le mouvement d'un liquide visqueux
  emplissant l'espace. *Acta Math.* **63**, 193–248. — A1, A12.
- **[VERIFIED]** Beale, J.T., Kato, T., Majda, A. (1984). Remarks on the breakdown
  of smooth solutions for the 3-D Euler equations. *Comm. Math. Phys.* **94**,
  61–66. — A6, and the log-inequality behind A5.
- **[VERIFIED]** Caffarelli, L., Kohn, R., Nirenberg, L. (1982). Partial regularity
  of suitable weak solutions of the Navier–Stokes equations. *CPAM* **35**,
  771–831. — A10.
- **[VERIFIED — Clay statement ref. 6]** Lin, F. (1998). A new proof of the
  Caffarelli–Kohn–Nirenberg theorem. *Comm. Pure & Appl. Math.* **51**, 241–257.
- **[VERIFIED]** Escauriaza, L., Seregin, G., Šverák, V. (2003). $L_{3,\infty}$-
  solutions of Navier–Stokes equations and backward uniqueness. *Russian Math.
  Surveys* **58**, 211–250. — A11.
- **[VERIFY]** Seregin, G. (2012). A certain necessary condition of potential blow
  up for Navier–Stokes equations. *Comm. Math. Phys.* **312**, 833–845.

## Geometric depletion (the mechanism — A7, A8, A9)

- **[VERIFIED]** Constantin, P., Fefferman, C. (1993). Direction of vorticity and
  the problem of global regularity for the Navier–Stokes equations. *Indiana
  Univ. Math. J.* **42**, 775–789. — A8.
- **[VERIFIED bibliographically; KERNEL BOUND NOT VERIFIED]** Constantin, P.
  (1994). Geometric statistics in turbulence. *SIAM Review* **36**, 73–98. — A7,
  the singular-integral representation of $\alpha$ with the $|\sin\angle|$ kernel
  bound. **The single most load-bearing citation in this program**: it is the
  mechanism every surviving route in Part C depends on, and it is what Tao's
  averaging destroys. The exact kernel bound could not be obtained from any
  reachable source and must be checked verbatim. Until then, A7 is used only
  qualitatively (stretching is suppressed where the direction field is locally
  coherent) and no quantitative constant from it enters any argument here.
- **[VERIFIED bibliographically; threshold CONFIRMED, full hypothesis NOT]**
  Beirão da Veiga, H., Berselli, L.C. (2002). On the regularizing effect of the
  vorticity direction in incompressible viscous flows. *Differential Integral
  Equations* **15**, 345–356. — A9. Confirmed: Lipschitz continuity of the
  vorticity direction can be replaced by $1/2$-Hölder continuity, and at
  $\beta=1/2$ one obtains $\omega\in L^\infty_tL^2_x$, hence smoothness. This is
  the origin of $\theta_c=1/2$ in Theorem C4. **Still unverified:** the precise
  function space and the region over which the Hölder bound is required (a
  neighbourhood of the high-vorticity set, with what uniformity). Theorem C4
  states it as uniform Hölder on an $r$-neighbourhood of $E_\lambda$; that form
  should be checked against the published hypothesis before C4 is relied upon.
- **[VERIFY]** Beirão da Veiga, H. (2016). Open problems concerning the Hölder
  continuity of the direction of vorticity for the Navier–Stokes equations.
  arXiv:1604.08083.

## Sparseness and the scaling gap (A13, A14 — the correct anchor)

- **[VERIFIED bibliographically; hypotheses PARTLY verified]** Grujić, Z. (2013).
  A geometric measure-type regularity criterion for solutions to the 3D
  Navier–Stokes equations. *Nonlinearity* **26**(1), 289 (DOI
  10.1088/0951-7715/26/1/289; arXiv:1111.0217). — A13.
  **Confirmed:** the framework is built on the scale of sparseness of the
  super-level sets of *the positive and negative parts of the vorticity
  components*; linear $\delta$-sparseness of a set $S$ around $x_0$ at scale $r$
  means there is a unit vector $d$ with the ratio of $|S\cap{\rm segment}|$ to
  the segment length bounded by $\delta$; the scale is comparable to the uniform
  lower bound on the radius of spatial analyticity; the mechanism is the
  harmonic-measure majorisation principle.
  **Still unverified:** the numerical threshold $\delta_0$ and the constant tying
  the sparseness scale to the analyticity radius. H2's $\delta_0$ stays UNSET.
- **[VERIFIED bibliographically; CLAIM CONTESTED]** Bradshaw, Z., Farhat, A.,
  Grujić, Z. (2019). An algebraic reduction of the 'scaling gap' in the
  Navier–Stokes regularity problem. *Arch. Ration. Mech. Anal.* **231**,
  1983–2005 (arXiv:1704.05546). — A14. See the Albritton–Bradshaw entry below:
  the scaling-gap claim is critically examined in a later paper co-authored by
  one of these authors. **Not relied upon anywhere in this project.**
- **[VERIFY]** Farhat, A., Grujić, Z., Leitmeyer, K. (2017). The space
  $B^{-1}_{\infty,\infty}$, volumetric sparseness, and 3D NSE.
- **[VERIFY]** Grujić, Z., Xu, L. Asymptotic criticality of the Navier–Stokes
  regularity problem.
- **[VERIFY]** Grujić, Z., Guberović, R. Localization of analytic regularity
  criteria on the vorticity and balance between the vorticity magnitude and
  coherence of the vorticity direction in the 3D NSE.
- **[VERIFY]** Bradshaw, Z., Grujić, Z. (2013). Blow-up scenarios for 3D NSE
  exhibiting sub-criticality with respect to the scaling of one-dimensional local
  sparseness. arXiv:1303.0257.
- **[VERIFIED bibliographically]** Albritton, D., Bradshaw, Z. (2022). Remarks on
  sparseness and regularity of Navier–Stokes solutions. *Nonlinearity* **35**,
  2858 (arXiv:2110.02187). **Two stated goals:** (i) a simple proof that
  sufficiently sparse Navier–Stokes solutions do not develop singularities — an
  alternative to Grujić's analyticity/harmonic-measure route, so the sufficiency
  direction of A13 has two independent proofs; (ii) an analysis of *the claims*
  that a priori sparseness estimates reduce the scaling gap. **Obtain and read
  this before freezing the preregistration** — it is the single most important
  unread source for this project.

## Log-weighted bmo depletion (A16 — the sharpest target)

- **[VERIFY]** Grujić, Z. (2026). Logarithmic Depletion of Vortex Stretching and
  Singularity Evasion in the 3D Navier–Stokes Equations. arXiv:2607.08866
  (13 July 2026). — A16, Theorem C5. Stretching eigenvalue recast as a
  singular-integral commutator; localised Coifman–Rochberg–Weiss estimate; direction
  in $\mathrm{bmo}_{1/|\log r|}$ (fails Dini, permits oscillatory defects) depletes
  stretching; $|\omega|$ forced into a subcritical Lorentz–Zygmund space under
  $L^{3/2,\infty}$ concentration.
  **Chain confirmed from the abstract:** commutator recast of the stretching
  eigenvalue → localised Coifman–Rochberg–Weiss with dyadic BMO tails →
  logarithmic envelope on shrinking super-level sets → interpolated De Giorgi →
  subcritical Lorentz–Zygmund for $|\omega|$ → **logarithmic gain transferred to
  the velocity, forcing the local 1D sparseness scale below the uniform radius of
  spatial analyticity** → no blow-up. That last step means A16 *terminates in
  A13's hypothesis*, so C5 is upstream of C3 rather than parallel to it.
  **Still to verify against the full text:** the precise definition of
  $\mathrm{bmo}_{1/|\log r|}$, the exact $L^{3/2,\infty}$ concentration
  hypothesis, and what "critical point singularities" restricts to. This paper
  post-dates the originating gauntlet runs and was not accounted for in them.
- **[VERIFY]** *On Decay of the Local Mean Oscillations of the Vorticity Direction
  in Critical Navier–Stokes Flows* (2026), arXiv:2609.05720. Directly relevant to
  the $\mu_\xi(r)$ measurement; obtain and read before freezing the preregistration.

## Spectral / Besov criteria (A15)

- **[VERIFIED bibliographically]** Cheskidov, A., Shvydkoy, R. (2010). The
  regularity of weak solutions of the 3D Navier–Stokes equations in
  $B^{-1}_{\infty,\infty}$. *Arch. Ration. Mech. Anal.* **195**(1), 159–169. —
  A15. And (2014) A unified approach to regularity problems for the 3D
  Navier–Stokes and Euler equations: the use of Kolmogorov's dissipation range.
  *J. Math. Fluid Mech.* **16**(2), 263–273. Basis for replacing the $\zeta_p$
  branch with the Cheskidov–Shvydkoy number.
- **[VERIFY]** Cheskidov, A., Dai, M. The determining wavenumber / dissipation
  wavenumber for Navier–Stokes.

## Enstrophy bounds and extreme states

- **[VERIFY]** Doering, C.R., Gibbon, J.D. (1995). *Applied Analysis of the
  Navier–Stokes Equations.* CUP. — A4.
- **[VERIFY]** Lu, L., Doering, C.R. (2008). Limits on enstrophy growth for
  solutions of the three-dimensional Navier–Stokes equations. *Indiana Univ.
  Math. J.* **57**, 2693–2727.
- **[VERIFIED]** Ayala, D., Protas, B. (2017). Extreme vortex states and the growth
  of enstrophy in 3D incompressible flows. *J. Fluid Mech.* **818**, 772–806.
- **[VERIFIED]** Kang, D., Yun, D., Protas, B. (2020). Maximum amplification of
  enstrophy in three-dimensional Navier–Stokes flows. *J. Fluid Mech.* **893**, A22.

## Numerical near-singular scenarios and DNS geometry

- **[VERIFY]** Kerr, R.M. (1993). Evidence for a singularity of the 3D
  incompressible Euler equations. *Phys. Fluids A* **5**, 1725–1746.
- **[VERIFY]** Kerr, R.M. (2018). Trefoil knot timescales for reconnection and
  helicity. — the $\sqrt{\nu}$ circulation scaling used as a solver-validation target.
- **[VERIFY]** Hou, T.Y., Li, R. (2006). Dynamic stability of the 3D axi-symmetric
  Navier–Stokes equations with swirl. — the $T^*$ drift test.
- **[VERIFY]** Luo, G., Hou, T.Y. (2014). Potentially singular solutions of the 3D
  axisymmetric Euler equations. *PNAS* **111**, 12968–12973.
- **[VERIFIED]** Rafner, J., Grujić, Z., et al. (2021). Geometry of turbulent
  dissipation and the Navier–Stokes regularity problem. *Scientific Reports* **11**,
  8824. — the only prior DNS measurement of the sparseness scale (Kida vortex and
  homogeneous isotropic turbulence). Reported that in the Kida case the sparseness
  scale was actualised "well beyond the guaranteed a priori bound and just beyond
  the critical bound sufficient for diffusion to fully engage".
  **Obtain the reported numerical values of $\delta$ and the sparseness scale**;
  they are the natural prior for H2 and the calibration target for our estimator.
- **[VERIFY]** Ashurst, W.T., Kerstein, A.R., Kerr, R.M., Gibson, C.H. (1987).
  Alignment of vorticity and scalar gradient with strain rate. *Phys. Fluids* **30**,
  2343. — the $e_2$ alignment expectation in `gamma.alignment_pdf`.
- **[VERIFY]** Meneveau, C., Sreenivasan, K.R. (1991). The multifractal nature of
  turbulent energy dissipation. *JFM* **224**, 429–484. — $D_q$ methodology.

## Obstructions (A17)

- **[VERIFIED]** Tao, T. (2016). Finite time blowup for an averaged
  three-dimensional Navier–Stokes equation. *J. Amer. Math. Soc.* **29**, 601–674.
- **[VERIFIED]** Elgindi, T. (2021). Finite-time singularity formation for
  $C^{1,\alpha}$ solutions to the incompressible Euler equations on $\mathbb R^3$.
  *Annals of Math.* **194**, 647–727.
- **[VERIFY]** Chen, J., Hou, T.Y. Stable nearly self-similar blowup of the 2D
  Boussinesq and 3D Euler equations with smooth data, I: Analysis and II: Rigorous
  Numerics (arXiv:2305.05660); *Singularity formation in 3D Euler equations with
  smooth initial data and boundary*, **PNAS (2025)**. Verify which statement is the
  published one before citing.

## Solvers

- **[VERIFY]** Mortensen, M., Langtangen, H.P. (2016). High performance Python for
  direct numerical simulations of turbulent flows. *Comput. Phys. Commun.* **203**,
  53–65. — spectralDNS, the proposed second solver for `MULTIPLE_SOLVERS_AGREE`.
