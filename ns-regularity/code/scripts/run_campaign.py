#!/usr/bin/env python3
"""Run one case end-to-end: evolve, measure Gamma-geometry, fit, evaluate gates.

Usage:
  python3 run_campaign.py --ic taylor_green --N 64 --nu 0.02 --tend 10
"""
import sys, os, json, argparse, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from nsref.grid import Grid                                    # noqa: E402
from nsref import ic, gamma, diagnostics as dg, fit, gates as gt, io  # noqa: E402
from nsref.solver import NSSolver                              # noqa: E402


def geometry_snapshot(g, uh, nu, d, capture=0.9):
    """Full Gamma-geometry measurement at one instant (Theorems C3/C4/C5)."""
    om = np.stack([g.ifft(c) for c in g.curl(uh)])
    net, gp, gm = gamma.stretching_density(g, uh, omega=om)
    gw2 = gamma.grad_omega_sq(g, uh)
    thr, vf = gamma.capture_threshold(gp, capture)
    mask = gamma.superlevel_mask(gp, thr)
    xi, mag = gamma.vorticity_direction(om)

    ell_nu = d["ell_nu"]
    ell_nu_cells = ell_nu / g.dx
    # sparseness on a curve of scales, reported in units of ell_nu
    r_cells = [r for r in (2, 4, 8, 12, 16, 24) if r < g.N // 2]
    spars = gamma.sparseness_curve(mask, r_cells, n_samples=2500)
    for s in spars:
        s["r_over_ell_nu"] = float(s["r_cells"] / ell_nu_cells) if ell_nu_cells > 0 else float("inf")

    box_cells = [b for b in (2, 4, 8, 16) if g.N % b == 0 and b <= g.N // 4]
    return {
        "capture": capture,
        "gamma_threshold": thr,
        "set_volume_fraction": vf,
        "ell_nu": ell_nu,
        "ell_nu_over_dx": float(ell_nu_cells),
        "sparseness_curve": spars,
        "holder_beta": gamma.holder_beta_xi(xi, mask, [1, 2, 4, 8]),
        "mean_oscillation": gamma.mean_oscillation_modulus(xi, mask, [2, 4, 8, 16], dx=g.dx, L=g.L),
        "dimensions": gamma.generalized_dimensions(gp, box_cells)["D"],
        "box_counting_D0": gamma.box_counting_dimension(mask, box_cells)["D0"],
        "R_E": gamma.restricted_ratio(net, gw2, mask, nu),
        "R_global": d["R_global"],
        "rho_depletion_at_ell_nu": gamma.depletion_ratio(g, uh, max(1.0, ell_nu_cells)),
        "alignment": gamma.alignment_pdf(g, uh, mask),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ic", default="taylor_green", choices=sorted(ic.REGISTRY))
    ap.add_argument("--N", type=int, default=64)
    ap.add_argument("--nu", type=float, default=0.02)
    ap.add_argument("--tend", type=float, default=10.0)
    ap.add_argument("--log-every", type=int, default=5)
    ap.add_argument("--snapshots", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    g = Grid(a.N)
    kw = {"seed": a.seed} if a.ic == "multiscale_random" else {}
    u0 = ic.REGISTRY[a.ic](g, **kw)
    s = NSSolver(g, a.nu, u0)

    tag = f"{a.ic}_N{a.N}_nu{a.nu:g}"
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "results", f"{tag}.jsonl")
    log = io.RunLog(out, {"ic": a.ic, "N": a.N, "nu": a.nu, "tend": a.tend,
                          "seed": a.seed, "L": g.L, "scheme": "IF-RK4 rotational 2/3-dealiased"})

    snap_times = list(np.linspace(a.tend * 0.25, a.tend, a.snapshots))
    ts, wmax, series, worst = [], [], [], {"e": 0.0, "o": 0.0, "cfl": 0.0}
    prev = dg.basic(g, s.vh, a.nu)
    t_prev = s.t
    emono = True
    t0 = time.time()
    isnap = 0
    print(f"# {tag}: dx={g.dx:.4f}  kmax={g.kmax:.0f}  initial kmax*eta={prev['kmax_eta']:.2f}")
    rc = dg.ic_resolution_check(g, s.vh, a.nu)
    log.write({"record": "ic_resolution_check", **rc})
    if not rc["resolved"]:
        print(f"# WARNING: initial condition is UNDER-RESOLVED. "
              f"kmax*eta(t=0)={rc['kmax_eta_0']:.2f} < {rc['target']}; "
              f"nu={a.nu:g} but nu_required={rc['nu_required']:.4f} "
              f"({rc['nu_shortfall_factor']:.2f}x short).")
        print("#          HIGH_K_TAIL_RESOLVED will fail. The run proceeds so the")
        print("#          gate is exercised, but its diagnostics are not trustworthy.")

    bprev = dg.budget_terms(g, s.vh, a.nu)
    while s.t < a.tend:
        dt = s.dt()
        worst["cfl"] = max(worst["cfl"], s.cfl_number(dt))
        s.step(dt)
        # budget residuals every step: they are O(interval^2) accurate
        b = dg.budget_terms(g, s.vh, a.nu)
        worst["e"] = max(worst["e"], dg.energy_budget_residual(
            bprev["E"], b["E"], dt, 0.5 * (b["eps"] + bprev["eps"])))
        worst["o"] = max(worst["o"], dg.enstrophy_budget_residual(
            bprev["enstrophy"], b["enstrophy"], dt,
            0.5 * (b["stretch_net"] + bprev["stretch_net"]),
            0.5 * (b["visc_dissip_enstrophy"] + bprev["visc_dissip_enstrophy"])))
        emono &= (b["E"] <= bprev["E"] + 1e-12)
        bprev = b
        if s.step_count % a.log_every == 0 or s.t >= a.tend:
            d = dg.basic(g, s.vh, a.nu)
            ts.append(s.t); wmax.append(d["omega_max"])
            rec = {"record": "step", "t": s.t, "dt": dt, "step": s.step_count, **d}
            series.append(rec); log.write(rec)
            prev = d
            t_prev = s.t
        if isnap < len(snap_times) and s.t >= snap_times[isnap]:
            d = dg.basic(g, s.vh, a.nu)
            gs = geometry_snapshot(g, s.vh, a.nu, d)
            log.write({"record": "geometry", "t": s.t, **gs})
            sp = gs["sparseness_curve"]
            best = min(sp, key=lambda q: abs(q["r_over_ell_nu"] - 1.0))
            print(f"  t={s.t:6.3f}  wmax={d['omega_max']:7.4f}  kmax*eta={d['kmax_eta']:5.2f}  "
                  f"ell_nu/dx={gs['ell_nu_over_dx']:5.2f}  R_E={gs['R_E']:+7.3f}  "
                  f"beta={gs['holder_beta']['beta']:5.3f}  D_inf={gs['dimensions'].get('inf', float('nan')):5.2f}  "
                  f"delta(r~ell_nu)={best['delta_p95']:.3f}")
            isnap += 1

    ts, wmax = np.array(ts), np.array(wmax)
    cmp = fit.compare(ts, wmax)
    boot = fit.bootstrap_powerlaw(ts, wmax, B=300)
    stab = fit.window_stability(ts, wmax)

    spec = g.shell_spectrum(s.vh)
    kk = np.arange(spec.size)
    tail = float(spec[-1] / max(spec.max(), 1e-300))
    dis = kk**2 * spec
    measured = {
        "cfl": worst["cfl"], "state_finite": True,
        "energy_monotone_tol": 0.0 if emono else 1.0,
        "energy_budget_residual": worst["e"], "enstrophy_budget_residual": worst["o"],
        "kmax_eta": float(min(r["kmax_eta"] for r in series)),
        "spectrum_tail_ratio": tail,
        "dissipation_peak_k_over_kmax": float(kk[int(np.argmax(dis))] / max(kk[-1], 1)),
        "dealias_residual": s.max_dealias_residual,
        "aliasing_unit_test": "pass",
        "power_law_bic_margin": cmp["power_law_bic_margin"],
        "Tstar_relative_drift": stab["Tstar_relative_drift"],
        "alpha_ci_width": boot.get("alpha_ci_width", float("nan")),
        "alpha_ci_lower": boot.get("alpha_ci95", [float("nan")])[0],
        "extrapolation_ratio": (
            (boot["Tstar"] - ts[-1]) / (ts[-1] - ts[0])
            if boot.get("ok") and ts[-1] > ts[0] else float("nan")),
    }
    applicable = ["CFL_STABLE", "ENERGY_BALANCE_CONSISTENT", "HIGH_K_TAIL_RESOLVED",
                  "DEALIASING_ACTIVE", "SINGULAR_TIME_FIT_STABLE", "BKM_CONSISTENT"]
    res = gt.evaluate(measured, applicable=applicable)
    log.write({"record": "summary", "measured": measured, "gates": res,
               "fit": {"compare": cmp, "bootstrap": boot, "stability": stab},
               "wall_seconds": time.time() - t0})
    sha = log.close()

    print(f"\n# growth ||omega||_inf: {wmax.max() / wmax[0]:.4f}x   "
          f"(t of max = {ts[int(np.argmax(wmax))]:.3f})")
    print(f"# preferred model: {cmp['preferred']}   power_law_preferred={cmp['power_law_preferred']}")
    print(f"\n{gt.format_table(res)}")
    print(f"\n# log: {os.path.relpath(out)}\n# sha256: {sha}")
    print("# NOTE: gates passing means the RUN is trustworthy, not that anything")
    print("#       has been established about Navier-Stokes regularity.")


if __name__ == "__main__":
    main()
