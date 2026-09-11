"""Model selection for ||omega||_inf(t): power-law blow-up vs the alternatives.

Fits are performed on y = log ||omega||_inf(t). A power law is preferred only if
it beats EVERY alternative by a preregistered BIC margin AND returns a >= 1
(Leray, A12) AND T* is stable across fit windows and resolutions.

The point of fitting M2-M4 is that a power law will fit almost any convex
increasing curve over a short window; the alternatives are what make the
comparison informative.
"""
import numpy as np

try:
    from scipy.optimize import curve_fit
    _HAVE_SCIPY = True
except Exception:                                     # pragma: no cover
    _HAVE_SCIPY = False


def _powerlaw(t, logA, a, Tstar):
    return logA - a * np.log(np.maximum(Tstar - t, 1e-12))


def _exponential(t, logA, b):
    return logA + b * t


def _stretched(t, logA, b, c):
    return logA + b * np.power(np.maximum(t, 1e-12), c)


def _double_exp(t, logA, B, c):
    return logA + B * np.exp(np.clip(c * t, -50, 50))


MODELS = {
    "power_law":    (_powerlaw,    3),
    "exponential":  (_exponential, 2),
    "stretched_exp": (_stretched,  3),
    "double_exp":   (_double_exp,  3),
}


def _ic_scores(resid, k):
    n = resid.size
    rss = float(np.sum(resid**2))
    rss = max(rss, 1e-300)
    ll_term = n * np.log(rss / n)
    aic = ll_term + 2 * k
    aicc = aic + (2 * k * (k + 1) / (n - k - 1)) if n - k - 1 > 0 else np.inf
    bic = ll_term + k * np.log(n)
    return {"rss": rss, "aic": float(aic), "aicc": float(aicc), "bic": float(bic)}


def _p0(name, t, y):
    tmax = float(t[-1])
    span = float(t[-1] - t[0]) or 1.0
    if name == "power_law":
        return [float(y[0]), 1.5, tmax + 0.25 * span]
    if name == "exponential":
        return [float(y[0]), max((y[-1] - y[0]) / span, 1e-3)]
    if name == "stretched_exp":
        return [float(y[0]), max((y[-1] - y[0]) / span, 1e-3), 1.0]
    return [float(y[0]), 0.1, 1.0 / span]


def fit_model(name, t, y):
    f, k = MODELS[name]
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    if t.size <= k + 1:
        return {"model": name, "ok": False, "reason": "too few points"}
    bounds = (-np.inf, np.inf)
    if name == "power_law":
        bounds = ([-np.inf, 0.0, t[-1] + 1e-9], [np.inf, 20.0, t[-1] + 100.0])
    try:
        if _HAVE_SCIPY:
            p, _ = curve_fit(f, t, y, p0=_p0(name, t, y), bounds=bounds, maxfev=40000)
        else:
            p = _nelder_mead(lambda q: np.sum((f(t, *q) - y) ** 2), _p0(name, t, y))
        resid = f(t, *p) - y
    except Exception as e:
        return {"model": name, "ok": False, "reason": str(e)[:120]}
    out = {"model": name, "ok": True, "params": [float(v) for v in p], "k": k}
    out.update(_ic_scores(resid, k))
    if name == "power_law":
        out["Tstar"], out["alpha"] = float(p[2]), float(p[1])
    return out


def _nelder_mead(fun, x0, iters=4000):               # pragma: no cover
    x0 = np.asarray(x0, float)
    n = x0.size
    sim = np.vstack([x0] + [x0 + (0.05 * abs(x0[i]) + 0.05) * np.eye(n)[i] for i in range(n)])
    fv = np.array([fun(s) for s in sim])
    for _ in range(iters):
        o = np.argsort(fv); sim, fv = sim[o], fv[o]
        cen = sim[:-1].mean(0)
        xr = cen + (cen - sim[-1]); fr = fun(xr)
        if fr < fv[0]:
            xe = cen + 2 * (cen - sim[-1]); fe = fun(xe)
            sim[-1], fv[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < fv[-2]:
            sim[-1], fv[-1] = xr, fr
        else:
            xc = cen + 0.5 * (sim[-1] - cen); fc = fun(xc)
            if fc < fv[-1]:
                sim[-1], fv[-1] = xc, fc
            else:
                sim[1:] = sim[0] + 0.5 * (sim[1:] - sim[0])
                fv[1:] = [fun(s) for s in sim[1:]]
    return sim[np.argmin(fv)]


def compare(t, omega_max, bic_margin=10.0):
    """Fit all models; report the preferred one and the margin over the runner-up."""
    y = np.log(np.asarray(omega_max, float))
    fits = {n: fit_model(n, t, y) for n in MODELS}
    ok = {n: f for n, f in fits.items() if f.get("ok")}
    if not ok:
        return {"fits": fits, "preferred": None, "reason": "no model converged"}
    ranked = sorted(ok.values(), key=lambda f: f["bic"])
    best = ranked[0]
    margin = (ranked[1]["bic"] - best["bic"]) if len(ranked) > 1 else np.inf
    pl = ok.get("power_law")
    pl_margin = (min(f["bic"] for n, f in ok.items() if n != "power_law") - pl["bic"]) if pl and len(ok) > 1 else -np.inf
    return {
        "fits": fits,
        "preferred": best["model"],
        "bic_margin_over_runner_up": float(margin),
        "power_law_bic_margin": float(pl_margin),
        "power_law_preferred": bool(pl and best["model"] == "power_law" and pl_margin >= bic_margin),
        "bic_margin_required": bic_margin,
    }


def bootstrap_powerlaw(t, omega_max, B=1000, seed=0):
    """Wild (Rademacher) residual bootstrap CI on T* and alpha."""
    t = np.asarray(t, float)
    y = np.log(np.asarray(omega_max, float))
    base = fit_model("power_law", t, y)
    if not base.get("ok"):
        return {"ok": False, "reason": base.get("reason")}
    f = MODELS["power_law"][0]
    yhat = f(t, *base["params"])
    resid = y - yhat
    rng = np.random.default_rng(seed)
    Ts, As = [], []
    for _ in range(B):
        yb = yhat + resid * rng.choice([-1.0, 1.0], size=resid.size)
        fb = fit_model("power_law", t, yb)
        if fb.get("ok"):
            Ts.append(fb["Tstar"]); As.append(fb["alpha"])
    if len(Ts) < B // 4:
        return {"ok": False, "reason": "bootstrap did not converge"}
    Ts, As = np.array(Ts), np.array(As)
    return {
        "ok": True, "B": len(Ts),
        "Tstar": base["Tstar"], "alpha": base["alpha"],
        "Tstar_ci95": [float(np.percentile(Ts, 2.5)), float(np.percentile(Ts, 97.5))],
        "alpha_ci95": [float(np.percentile(As, 2.5)), float(np.percentile(As, 97.5))],
        "alpha_ci_width": float(np.percentile(As, 97.5) - np.percentile(As, 2.5)),
    }


def window_stability(t, omega_max, fracs=(1.0, 0.9, 0.8)):
    """Refit on truncated windows. A real T* does not drift (Hou-Li drift test)."""
    t = np.asarray(t, float)
    w = np.asarray(omega_max, float)
    out = []
    for fr in fracs:
        n = max(6, int(len(t) * fr))
        f = fit_model("power_law", t[:n], np.log(w[:n]))
        out.append({"frac": fr, "n": n,
                    "Tstar": f.get("Tstar", float("nan")),
                    "alpha": f.get("alpha", float("nan")), "ok": f.get("ok", False)})
    Ts = [o["Tstar"] for o in out if o["ok"] and np.isfinite(o["Tstar"])]
    drift = (max(Ts) - min(Ts)) / abs(np.mean(Ts)) if len(Ts) >= 2 and np.mean(Ts) != 0 else float("nan")
    return {"windows": out, "Tstar_relative_drift": float(drift)}
