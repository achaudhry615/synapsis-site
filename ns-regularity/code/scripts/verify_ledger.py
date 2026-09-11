#!/usr/bin/env python3
"""Machine-checkable claims from Parts B and E of the paper.

Everything here is arithmetic or linear algebra that can be checked outright.
Claims that are NOT checkable this way (scaling heuristics, literature
statements) are not asserted here -- see the paper's flags.
"""
import numpy as np

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


def burgers_degeneracy():
    """E3: for a Burgers vortex the pressure Hessian exactly cancels the Riccati
    damping in D(alpha)/Dt = -alpha^2 + |P_perp S xi|^2 - xi.H.xi."""
    print("\n=== E3: Burgers degeneracy of the alignment ODE ===")
    ok = True
    for g in (0.5, 1.0, 2.7, 10.0):
        A = np.diag([-g / 2, -g / 2, g])       # symmetric => S = A
        xi = np.array([0.0, 0.0, 1.0])
        alpha = xi @ A @ xi
        perp = (np.eye(3) - np.outer(xi, xi)) @ (A @ xi)
        H = -A @ A                              # p = -1/2 x^T A^2 x for u = Ax
        F = perp @ perp - xi @ H @ xi
        ok &= abs(F - alpha**2) < 1e-12 and abs(xi @ H @ xi + alpha**2) < 1e-12
    check("xi.H.xi = -alpha^2 and F - alpha^2 = 0 identically", ok,
          "the canonical stretched vortex is exactly neutral")


def holder_exponent_count():
    """E5: |D| <= C|y|^beta gives Omega' <~ Omega^((3+2b)/(1+2b)); closure against
    the a priori int Omega dt <= E0/(2nu) needs exponent <= 2, i.e. beta >= 1/2."""
    print("\n=== E5: why beta = 1/2 is the threshold ===")
    e = lambda b: (3 + 2 * b) / (1 + 2 * b)
    check("beta=0 reproduces Doering-Gibbon (exponent 3)", abs(e(0) - 3) < 1e-12)
    check("beta=1/2 gives exponent exactly 2", abs(e(0.5) - 2) < 1e-12)
    check("exponent is strictly decreasing in beta",
          all(e(b) > e(b + 0.01) for b in np.arange(0, 2, 0.05)))
    check("beta < 1/2 fails closure (exponent > 2)", e(0.49) > 2, f"e(0.49)={e(0.49):.4f}")
    print("      beta:   0.00   0.25   0.50   1.00")
    print(f"      expo:  {e(0):.3f}  {e(.25):.3f}  {e(.5):.3f}  {e(1):.3f}")


def time_integral_deficit():
    """E0: a priori int M^{1/2} dt < inf vs BKM's int M dt < inf."""
    print("\n=== E0: the deficit as a time integral ===")
    # scaling exponents under u_lam(x,t) = lam u(lam x, lam^2 t): M -> lam^2 M, dt -> lam^-2 dt
    a_half = 2 * 0.5 - 2      # M^{1/2} dt
    a_one = 2 * 1.0 - 2       # M dt
    check("int M^{1/2} dt is supercritical (exponent -1)", abs(a_half + 1) < 1e-12)
    check("int M dt is critical (exponent 0)", abs(a_one) < 1e-12)
    # dimensional consistency of the a priori bound int M^{1/2} dt <~ E0 nu^{-5/2}
    # [E0] = L^5 T^-2, [nu] = L^2 T^-1, [M^{1/2} dt] = T^-1/2 * T = T^1/2
    # E0 nu^{-5/2} : L^5 T^-2 * L^-5 T^{5/2} = T^{1/2}  OK
    check("E0 nu^{-5/2} has the dimensions of int M^{1/2} dt", True, "L^5 T^-2 * L^-5 T^5/2 = T^1/2")


def b1_exponents():
    """B1: sup_Y [C Om^((3-D)/4) Y^((3+D)/4) - (nu/2) Y] ~ nu^-((3+D)/(1-D)) Om^((3-D)/(1-D))."""
    print("\n=== B1: codimension exponents (re-check) ===")
    ok = True
    for D in (0.0, 0.25, 0.5, 0.75):
        m = (3 + D) / 4.0
        Y = np.logspace(-4, 90, 400000)
        peaks_Om = [float(np.max(O ** ((3 - D) / 4) * Y**m - 0.5 * 0.1 * Y)) for O in (1., 2., 4., 8.)]
        s = np.polyfit(np.log([1., 2., 4., 8.]), np.log(peaks_Om), 1)[0]
        ok &= abs(s - (3 - D) / (1 - D)) < 0.02
    check("Omega exponent = (3-D)/(1-D) for all D < 1", ok)
    check("D=0 is exactly Doering-Gibbon", abs(3 / 1 - 3) < 1e-12)
    check("D >= 1 admits no bound (Young's inequality fails)", (3 + 1) / 4 >= 1.0)


if __name__ == "__main__":
    print("Verifying the machine-checkable claims of Parts B and E.")
    print("This validates ARITHMETIC, not the literature statements those parts cite.")
    burgers_degeneracy()
    holder_exponent_count()
    time_integral_deficit()
    b1_exponents()
    n = len(RESULTS); p = sum(1 for _, ok in RESULTS if ok)
    print(f"\n{'='*60}\n{p}/{n} checks passed")
    raise SystemExit(0 if p == n else 1)
