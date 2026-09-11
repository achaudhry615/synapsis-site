#!/usr/bin/env python3
"""Build a cross-case table from results/*.jsonl.

Reports the geometric hypotheses (C3/C4/C5), the decisive ratio R_E (B2), the
gate outcome, and -- crucially -- whether r ~ ell_nu was RESOLVABLE at all. A
delta of 1.000 with ell_nu/dx ~ 1 is not a measurement; the table marks it.
"""
import sys, os, glob, json

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from nsref import io  # noqa: E402

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")


def _c5_bounded(mo):
    """Theorem C5: is mu_xi(r)*log(L/r) bounded as r -> 0?

    Entries are ordered small-r first, so 'bounded' means the product SHRINKS
    toward small r. Returns 'yes' / 'NO' / '?'.
    """
    vals = [o.get("mu_times_logL_over_r") for o in mo]
    vals = [v for v in vals if v is not None and v == v]
    if len(vals) < 2:
        return "?"
    return "yes" if vals[0] < vals[-1] else "NO"


def summarize(path):
    recs = io.read_jsonl(path)
    meta = next((r for r in recs if r.get("record") == "meta"), {})
    geo = [r for r in recs if r.get("record") == "geometry"]
    summ = next((r for r in recs if r.get("record") == "summary"), {})
    steps = [r for r in recs if r.get("record") == "step"]
    if not steps:
        return None
    wmax = [r["omega_max"] for r in steps]
    g = geo[len(geo) // 2] if geo else {}
    sp = g.get("sparseness_curve", [])
    near = min(sp, key=lambda q: abs(q["r_over_ell_nu"] - 1.0)) if sp else {}
    gates = summ.get("gates", {})
    applicable = [k for k, v in gates.get("gates", {}).items()
                  if v["status"] != "NOT_APPLICABLE"]
    npass = sum(1 for k in applicable if gates["gates"][k]["status"] == "PASS")
    return {
        "case": f"{meta.get('ic')}_N{meta.get('N')}",
        "nu": meta.get("nu"),
        "growth": max(wmax) / wmax[0],
        "ell_nu_dx": g.get("ell_nu_over_dx", float("nan")),
        "delta": near.get("delta_p95", float("nan")),
        "r_over_ell": near.get("r_over_ell_nu", float("nan")),
        "beta": g.get("holder_beta", {}).get("beta", float("nan")),
        "D_inf": g.get("dimensions", {}).get("inf", float("nan")),
        "R_E": g.get("R_E", float("nan")),
        "R_global": g.get("R_global", float("nan")),
        "mu_slope": (g.get("mean_oscillation") or [{}])[0].get("mu_loglog_slope", float("nan")),
        "c5_bounded": _c5_bounded(g.get("mean_oscillation") or []),
        "preferred": summ.get("fit", {}).get("compare", {}).get("preferred"),
        "pl": summ.get("fit", {}).get("compare", {}).get("power_law_preferred"),
        "gates": f"{npass}/{len(applicable)}",
    }


def main():
    paths = sorted(glob.glob(os.path.join(RESULTS, "*.jsonl")))
    rows = [r for r in (summarize(p) for p in paths) if r]
    if not rows:
        print("no results found"); return
    hdr = (f"{'case':<24}{'growth':>8}{'l_nu/dx':>9}{'delta':>7}"
           f"{'beta':>6}{'mu_sl':>7}{'C5':>5}{'D_inf':>7}{'R_E':>9}{'R_glob':>8}"
           f"{'model':>14}{'gates':>7}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['case']:<24}{r['growth']:>8.4f}{r['ell_nu_dx']:>9.2f}"
              f"{r['delta']:>7.3f}{r['beta']:>6.2f}{r['mu_slope']:>7.2f}"
              f"{r['c5_bounded']:>5}{r['D_inf']:>7.2f}{r['R_E']:>9.4f}"
              f"{r['R_global']:>8.3f}{str(r['preferred']):>14}{r['gates']:>7}")
    print("-" * len(hdr))
    unres = [r for r in rows if r["ell_nu_dx"] < 2.0]
    if unres:
        print(f"\nWARNING: {len(unres)}/{len(rows)} cases have ell_nu < 2 dx.")
        print("For these, delta(r ~ ell_nu) is SUB-GRID and its value carries no")
        print("information about the Theorem C3 hypothesis (spec/gamma_measurement.md §4).")
        print("The sparseness column must NOT be read as a measurement for these rows.")
    print("\nColumns: delta = 1D sparseness p95 near r ~ ell_nu (Thm C3);")
    print("  beta = Holder exponent of xi on the set (Thm C4, threshold 0.5;")
    print("         values > 1 mean a locally CONSTANT direction field, which")
    print("         satisfies C4 trivially -- the Holder range is (0,1]);")
    print("  mu_sl = log-log slope of the mean-oscillation modulus (~1 Lipschitz, ~0 none);")
    print("  C5 = is mu_xi(r)*log(L/r) bounded as r -> 0 (Thm C5 hypothesis);")
    print("  R_E = set-restricted stretching/dissipation (Prop. B2 -- the decisive one);")
    print("  R_glob = the global ratio, which is > 1 iff enstrophy is rising (Cor. B2.1,")
    print("           uninformative, shown only to make that point concrete).")
    print("\nNo row of this table is evidence for or against global regularity.")


if __name__ == "__main__":
    main()
