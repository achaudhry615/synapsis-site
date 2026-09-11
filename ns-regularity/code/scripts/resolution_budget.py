#!/usr/bin/env python3
"""Resolution required to MEASURE the Theorem-C3 hypothesis, vs. to resolve the flow.

The sparseness hypothesis lives at r ~ ell_nu = (nu/||omega||_inf)^{1/2}. To
estimate delta(r) at that scale the segment must span several grid cells, so the
binding constraint is ell_nu/dx, NOT the usual k_max*eta >= 2 DNS criterion.

  eta/dx    = (k_max eta)/pi                      [k_max = pi N / L]
  ell_nu/dx = (ell_nu/eta) * (k_max eta)/pi

ell_nu/eta is an O(1) flow-dependent ratio, measured at runtime (see results/).
"""
import numpy as np

L = 2 * np.pi
RATIO = 0.63   # ell_nu/eta measured on the 64^3 Taylor-Green run at its enstrophy peak


def table(Omega_v=1.0):
    print(f"assumed volume-averaged enstrophy Omega_v = {Omega_v}, L = 2*pi, "
          f"ell_nu/eta = {RATIO}\n")
    print(f"{'N':>5} {'k_max':>6} | {'target kmax*eta':>15} {'required nu':>13} "
          f"{'eta/dx':>7} {'ell_nu/dx':>10}  verdict")
    print("-" * 88)
    for N in (128, 256, 512, 1024):
        kmax = np.pi * N / L
        for target in (2.0, 4.0, 10.0):
            eta = target / kmax
            # eta = (nu^3/eps)^{1/4}, eps = 2 nu Omega_v  =>  nu = eta^2 sqrt(2 Omega_v)
            nu = eta**2 * np.sqrt(2 * Omega_v)
            eta_dx = target / np.pi
            ell_dx = RATIO * eta_dx
            if ell_dx >= 4:
                v = "delta(r~ell_nu) measurable"
            elif ell_dx >= 2:
                v = "marginal (5-point segment)"
            else:
                v = "NOT measurable: r~ell_nu is sub-grid"
            print(f"{N:>5} {kmax:6.0f} | {target:>15.1f} {nu:13.3e} "
                  f"{eta_dx:7.2f} {ell_dx:10.2f}  {v}")
        print()


if __name__ == "__main__":
    table()
    print("Conclusions (these set the campaign's nu schedule):")
    print("  * k_max*eta >= 2 -- the standard DNS criterion -- leaves ell_nu at ~0.4 dx.")
    print("    At that resolution delta(r ~ ell_nu) is identically 1 and carries NO")
    print("    information. This was confirmed empirically on the 64^3 run.")
    print("  * Measuring the Theorem-C3 hypothesis needs k_max*eta ~ 10, i.e. a factor")
    print("    ~5 in linear resolution beyond a 'well-resolved' DNS at the same nu,")
    print("    equivalently ~125x the cost, or a correspondingly larger nu at fixed N.")
    print("  * Consequence: a 512^3 run CAN measure the geometric hypotheses, but only")
    print("    at nu ~ 1e-3, i.e. a low-Reynolds flow that will certainly not blow up.")
    print("    The geometric measurement and the blow-up search are therefore competing")
    print("    for the same grid, and the campaign must choose the geometry (see")
    print("    paper section 6.3). Pretending one run does both is how this goes wrong.")
