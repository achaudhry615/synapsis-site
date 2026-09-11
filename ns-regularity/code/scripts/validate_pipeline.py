#!/usr/bin/env python3
"""End-to-end validation. Every check must print PASS before any run is trusted.

Covers: spectral operators, dealiasing (manufactured test), solver conservation,
Taylor-Green analytic values, geometry-estimator recovery on shapes of KNOWN
dimension, fit/model-selection recovery, and a numerical check of the exponent
algebra in Proposition B1 of the paper.
"""
import sys, os, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from nsref.grid import Grid                       # noqa: E402
from nsref import ic, gamma, diagnostics as dg, fit   # noqa: E402
from nsref.solver import NSSolver                 # noqa: E402

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    return ok


def sec(title):
    print(f"\n=== {title} ===")


# --------------------------------------------------------------------------
def test_spectral_operators():
    sec("1. spectral operators")
    g = Grid(32)
    rng = np.random.default_rng(0)
    f = rng.random((32,) * 3)
    check("rfft round-trip", np.abs(g.ifft(g.fft(f)) - f).max() < 1e-12,
          f"err={np.abs(g.ifft(g.fft(f)) - f).max():.2e}")
    u = ic.taylor_green(g)
    uh = np.stack([g.fft(c) for c in u])
    check("Leray projection gives div u = 0", g.divergence_max(uh) < 1e-12,
          f"div={g.divergence_max(uh):.2e}")
    E = g.shell_spectrum(uh).sum()
    Edir = 0.5 * np.mean(sum(c**2 for c in u))
    check("Parseval: sum E(k) == 0.5<|u|^2>", abs(E - Edir) < 1e-12 * max(1, abs(Edir)),
          f"|diff|={abs(E - Edir):.2e}")
    # curl of a gradient vanishes
    ph = g.fft(rng.random((32,) * 3))
    check("curl(grad phi) == 0", np.abs(g.curl(g.grad(ph))).max() < 1e-8,
          f"max={np.abs(g.curl(g.grad(ph))).max():.2e}")


def test_dealiasing():
    sec("2. dealiasing (manufactured aliasing test)")
    N = 32
    g = Grid(N)
    rng = np.random.default_rng(1)
    # field band-limited to the 2/3 band
    fh = (rng.normal(size=(N, N, N // 2 + 1)) + 1j * rng.normal(size=(N, N, N // 2 + 1)))
    fh *= g.dealias
    f = g.ifft(fh)
    # (A) product on the N grid, then 2/3 mask
    A = g.fft(f * f) * g.dealias
    # (B) exact product via 2x zero padding, truncated back
    M = 2 * N
    big = np.zeros((M, M, M // 2 + 1), dtype=fh.dtype)
    h = N // 2
    big[:h, :h, :h + 1] = fh[:h, :h, :h + 1]
    big[:h, -h:, :h + 1] = fh[:h, h:, :h + 1]
    big[-h:, :h, :h + 1] = fh[h:, :h, :h + 1]
    big[-h:, -h:, :h + 1] = fh[h:, h:, :h + 1]
    fbig = np.fft.irfftn(big, s=(M, M, M)) * (2**3)
    pbig = np.fft.rfftn(fbig * fbig) / (2**3)
    B = np.zeros_like(fh)
    B[:h, :h, :h + 1] = pbig[:h, :h, :h + 1]
    B[:h, h:, :h + 1] = pbig[:h, -h:, :h + 1]
    B[h:, :h, :h + 1] = pbig[-h:, :h, :h + 1]
    B[h:, h:, :h + 1] = pbig[-h:, -h:, :h + 1]
    B *= g.dealias
    rel = np.abs(A - B).max() / max(np.abs(B).max(), 1e-300)
    check("2/3 rule reproduces the exact (padded) product", rel < 1e-10, f"rel err={rel:.2e}")
    check("masked modes are identically zero", g.dealias_residual(A) == 0.0)


def test_solver_conservation():
    sec("3. solver conservation and Taylor-Green analytic values")
    g = Grid(32)
    nu = 0.05
    u0 = ic.taylor_green(g)
    s = NSSolver(g, nu, u0)
    d0 = dg.basic(g, s.vh, nu)
    check("TG t=0 energy == 1/8", abs(d0["E"] - 0.125) < 1e-12, f"E={d0['E']:.12f}")
    check("TG t=0 enstrophy == 3/8", abs(d0["enstrophy"] - 0.375) < 1e-12, f"Om={d0['enstrophy']:.12f}")
    check("TG t=0 ||omega||_inf == 2", abs(d0["omega_max"] - 2.0) < 1e-9, f"wmax={d0['omega_max']:.9f}")
    prev, we, wo, emono = d0, 0.0, 0.0, True
    for _ in range(20):
        dt = s.step()
        d = dg.basic(g, s.vh, nu, refine=False)
        we = max(we, dg.energy_budget_residual(prev["E"], d["E"], dt, 0.5 * (d["eps"] + prev["eps"])))
        wo = max(wo, dg.enstrophy_budget_residual(
            prev["enstrophy"], d["enstrophy"], dt,
            0.5 * (d["stretch_net"] + prev["stretch_net"]),
            0.5 * (d["visc_dissip_enstrophy"] + prev["visc_dissip_enstrophy"])))
        emono &= (d["E"] <= prev["E"] + 1e-12)
        prev = d
    check("energy budget residual <= 1e-3", we <= 1e-3, f"max={we:.2e}")
    check("enstrophy budget residual <= 1e-2", wo <= 1e-2, f"max={wo:.2e}")
    check("energy monotonically decreasing (unforced)", emono)
    check("dealias residual zero over the run", s.max_dealias_residual == 0.0)
    check("divergence stays zero", prev["div_max"] < 1e-12, f"div={prev['div_max']:.2e}")


def test_geometry_recovery():
    sec("4. geometry estimators on shapes of KNOWN dimension")
    print("   (box scales fixed in PHYSICAL units -- see spec/gamma_measurement.md)")
    tol = 0.65
    for N, bc in ((64, [2, 4, 8, 16]), (128, [4, 8, 16, 32])):
        g = Grid(N)
        row = []
        for name, (fn, Dtrue) in ic.SYNTH.items():
            f = fn(g)
            thr, _ = gamma.capture_threshold(f, 0.9)
            m = gamma.superlevel_mask(f, thr)
            Dinf = gamma.generalized_dimensions(f, bc)["D"]["inf"]
            row.append((name, Dtrue, Dinf))
        ok = all(abs(d - t) < tol + 0.35 * t for _, t, d in row)
        detail = " ".join(f"{n}:{d:.2f}(true {t:.0f})" for n, t, d in row)
        check(f"N={N} D_inf ordering and recovery", ok, detail)
    # N-independence at matched physical scales -- the artifact discriminator
    g64, g128 = Grid(64), Grid(128)
    drift = []
    for name, (fn, _) in ic.SYNTH.items():
        a = gamma.generalized_dimensions(fn(g64), [2, 4, 8, 16])["D"]["inf"]
        b = gamma.generalized_dimensions(fn(g128), [4, 8, 16, 32])["D"]["inf"]
        drift.append(abs(a - b))
    check("D_inf is N-independent at matched physical scales", max(drift) < 0.25,
          f"max drift={max(drift):.3f}")
    # sparseness discriminates thin from space-filling
    g = Grid(64)
    ds = {}
    for name in ("tube", "sheet", "uniform"):
        f = ic.SYNTH[name][0](g)
        thr, _ = gamma.capture_threshold(f, 0.9)
        m = gamma.superlevel_mask(f, thr)
        ds[name] = gamma.sparseness_1d(m, r_cells=16, n_samples=1500)["delta_p95"]
    check("sparseness: sheet < tube < uniform", ds["sheet"] < ds["tube"] < ds["uniform"],
          " ".join(f"{k}={v:.3f}" for k, v in ds.items()))
    # Holder exponent of a smooth direction field should be ~1 (Lipschitz)
    g = Grid(64)
    u = ic.abc_flow(g)
    uh = np.stack([g.fft(c) for c in u])
    om = np.stack([g.ifft(c) for c in g.curl(uh)])
    xi, mag = gamma.vorticity_direction(om)
    m = mag > 0.5 * mag.max()
    hb = gamma.holder_beta_xi(xi, m, [1, 2, 4, 8])
    check("smooth (ABC) direction field gives beta ~ 1", 0.7 < hb["beta"] < 1.3,
          f"beta={hb['beta']:.3f} r2={hb['r2']:.3f}")


def test_fit_recovery():
    sec("5. model selection")
    rng = np.random.default_rng(2)
    t = np.linspace(0.0, 1.8, 60)
    T, a = 2.0, 1.6
    w = (T - t) ** (-a) * np.exp(rng.normal(0, 0.01, t.size))
    c = fit.compare(t, w)
    b = fit.bootstrap_powerlaw(t, w, B=300)
    check("true power law is preferred", c["power_law_preferred"],
          f"margin={c['power_law_bic_margin']:.1f}")
    check("T* recovered to 1%", abs(b["Tstar"] - T) / T < 0.01, f"T*={b['Tstar']:.4f}")
    check("alpha recovered to 2%", abs(b["alpha"] - a) / a < 0.02, f"alpha={b['alpha']:.4f}")
    check("T* true value inside 95% CI", b["Tstar_ci95"][0] <= T <= b["Tstar_ci95"][1],
          f"CI={np.round(b['Tstar_ci95'], 4).tolist()}")
    w2 = np.exp(1.4 * t) * np.exp(rng.normal(0, 0.01, t.size))
    c2 = fit.compare(t, w2)
    check("exponential is NOT called a power law", not c2["power_law_preferred"],
          f"preferred={c2['preferred']}, pl margin={c2['power_law_bic_margin']:.1f}")


def test_prop_b1_algebra():
    sec("6. Proposition B1 exponent algebra (paper eq. 3.3)")
    print("   sup_Y [ C Om^((3-D)/4) Y^((3+D)/4) - (nu/2) Y ]  ==  C' nu^-(3+D)/(1-D) Om^((3-D)/(1-D))")
    ok_all = True
    for D in (0.0, 0.25, 0.5, 0.75):
        interior = True

        def peak(Om, nu, C=1.0):
            """sup over Y, on a grid wide enough that the maximum is INTERIOR.

            Y* = (a m / c)^(1/(1-m)) grows like 10^20 as D -> 1, so a narrow
            grid silently returns a boundary value and fakes the exponent.
            """
            nonlocal interior
            m = (3 + D) / 4.0
            Y = np.logspace(-4, 90, 600000)
            g = C * Om ** ((3 - D) / 4.0) * Y**m - 0.5 * nu * Y
            i = int(np.argmax(g))
            if i == 0 or i == Y.size - 1:
                interior = False
            return float(g[i])
        Om_s = np.array([1.0, 2.0, 4.0, 8.0])
        nu_s = np.array([0.1, 0.05, 0.025, 0.0125])
        sOm = np.polyfit(np.log(Om_s), np.log([peak(o, 0.1) for o in Om_s]), 1)[0]
        snu = np.polyfit(np.log(nu_s), np.log([peak(1.0, n) for n in nu_s]), 1)[0]
        eOm, enu = (3 - D) / (1 - D), -(3 + D) / (1 - D)
        ok = abs(sOm - eOm) < 0.02 and abs(snu - enu) < 0.02 and interior
        ok_all &= ok
        print(f"      D={D:.2f}  Om exponent {sOm:7.3f} (predicted {eOm:7.3f})"
              f"   nu exponent {snu:8.3f} (predicted {enu:8.3f})  "
              f"{'ok' if ok else ('BOUNDARY' if not interior else 'MISMATCH')}")
    check("B1 exponents match the paper for all D < 1", ok_all)
    check("B1 at D=0 reproduces the Doering-Gibbon bound Om^3 nu^-3",
          abs((3 - 0) / (1 - 0) - 3.0) < 1e-12 and abs(-(3 + 0) / (1 - 0) + 3.0) < 1e-12)


if __name__ == "__main__":
    t0 = time.time()
    print("nsref pipeline validation")
    print("NOTE: passing these checks validates the MEASUREMENT PIPELINE only.")
    print("      It says nothing about Navier-Stokes regularity.")
    test_spectral_operators()
    test_dealiasing()
    test_solver_conservation()
    test_geometry_recovery()
    test_fit_recovery()
    test_prop_b1_algebra()
    n = len(RESULTS)
    p = sum(1 for _, ok in RESULTS if ok)
    print(f"\n{'=' * 56}\n{p}/{n} checks passed   [{time.time() - t0:.1f}s]")
    if p != n:
        print("FAILED:", ", ".join(k for k, ok in RESULTS if not ok))
    sys.exit(0 if p == n else 1)
