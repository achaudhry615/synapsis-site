"""Gamma-field geometry: super-level sets, sparseness, dimensions, direction regularity.

Implements the measurable hypotheses of Theorems C3 / C4 / C5 of
`paper/ns_gamma_codimension.md`:

  C3  sparseness_1d(...)            -> delta(r) at r ~ ell_nu           [NS-G02]
  C4  holder_beta_xi(...)           -> beta, threshold 1/2              [NS-G04]
  C5  mean_oscillation_modulus(...) -> mu_xi(r) vs 1/|log r|            [NS-G06]

Secondary diagnostics (Prop. B4: dimension is NOT the theorem-relevant object):
  generalized_dimensions(...), volume_slope_dimension(...)
"""
import numpy as np

# 13 canonical directions of the cubic lattice (axes, face diagonals, body diagonals)
_DIRS13 = np.array([
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1),
    (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
], dtype=np.float64)
_DIRS13 /= np.linalg.norm(_DIRS13, axis=1, keepdims=True)


def stretching_density(grid, uh, omega=None):
    """Return (net, gamma_plus, gamma_minus) where net = omega . S omega.

    Uses omega_i omega_j d_j u_i, which equals omega.S.omega because omega_i omega_j
    is symmetric and annihilates the antisymmetric part of the velocity gradient.
    Avoids materialising the full (3,3,N,N,N) strain tensor.
    """
    if omega is None:
        omega = np.stack([grid.ifft(c) for c in grid.curl(uh)])
    k = (grid.kx, grid.ky, grid.kz)
    net = np.zeros((grid.N,) * 3, dtype=grid.dtype)
    for i in range(3):
        for j in range(3):
            duij = grid.ifft(1j * k[j] * uh[i])       # d u_i / d x_j
            net += omega[i] * omega[j] * duij
    return net, np.maximum(net, 0.0), np.maximum(-net, 0.0)


def strain_frobenius(grid, uh):
    """|S| = sqrt(S_ij S_ij), computed without storing the full tensor."""
    k = (grid.kx, grid.ky, grid.kz)
    du = [[grid.ifft(1j * k[j] * uh[i]) for j in range(3)] for i in range(3)]
    acc = np.zeros((grid.N,) * 3, dtype=grid.dtype)
    for i in range(3):
        for j in range(3):
            s = 0.5 * (du[i][j] + du[j][i])
            acc += s * s
    return np.sqrt(acc)


def grad_omega_sq(grid, uh):
    """|grad omega|^2 pointwise."""
    wh = grid.curl(uh)
    k = (grid.kx, grid.ky, grid.kz)
    acc = np.zeros((grid.N,) * 3, dtype=grid.dtype)
    for i in range(3):
        for j in range(3):
            d = grid.ifft(1j * k[j] * wh[i])
            acc += d * d
    return acc


# --------------------------------------------------------------------------
# super-level sets
# --------------------------------------------------------------------------
def capture_threshold(field, capture=0.9):
    """Smallest-volume super-level threshold carrying `capture` of sum(field).

    This is the primary thresholding rule: it is scale-free and does not require
    guessing an absolute level. Returns (threshold, volume_fraction).
    """
    f = field.ravel()
    f = f[f > 0]
    if f.size == 0:
        return np.inf, 0.0
    s = np.sort(f)[::-1]
    c = np.cumsum(s)
    total = c[-1]
    idx = int(np.searchsorted(c, capture * total))
    idx = min(idx, s.size - 1)
    return float(s[idx]), float(idx + 1) / field.size


def superlevel_mask(field, threshold):
    return field > threshold


# --------------------------------------------------------------------------
# C3: 1D sparseness  (the theorem-relevant geometric quantity)
# --------------------------------------------------------------------------
def sparseness_1d(mask, r_cells, n_samples=4000, n_along=None, dirs=_DIRS13, rng=None):
    """1D sparseness ratio delta of `mask` at scale r_cells (in grid cells).

    For each sampled point of the set and each direction n, compute the occupied
    fraction of the segment (x0 - r n, x0 + r n). Take the MINIMUM over
    directions (a set is delta-thin if SOME direction is thin -- Lemma C1).

    Returns dict with the distribution of that per-point minimum. Small delta
    means thin, i.e. favourable to Theorem C3.
    """
    rng = np.random.default_rng(0 if rng is None else rng)
    N = mask.shape[0]
    pts = np.argwhere(mask)
    if pts.shape[0] == 0:
        return {"n_points": 0, "delta_mean": float("nan"), "delta_p50": float("nan"),
                "delta_p95": float("nan"), "delta_max": float("nan"), "r_cells": r_cells}
    if pts.shape[0] > n_samples:
        pts = pts[rng.choice(pts.shape[0], n_samples, replace=False)]
    if n_along is None:
        n_along = max(5, int(2 * r_cells) + 1)

    s = np.linspace(-r_cells, r_cells, n_along)          # (A,)
    # positions: (P, D, A, 3)
    offs = s[None, :, None] * dirs[:, None, :]           # (D, A, 3)
    idx = np.rint(pts[:, None, None, :] + offs[None]).astype(np.int64) % N
    occ = mask[idx[..., 0], idx[..., 1], idx[..., 2]]    # (P, D, A)
    frac = occ.mean(axis=2)                              # (P, D)
    delta = frac.min(axis=1)                             # per-point min over directions
    return {
        "n_points": int(pts.shape[0]),
        "r_cells": float(r_cells),
        "delta_mean": float(delta.mean()),
        "delta_p50": float(np.percentile(delta, 50)),
        "delta_p95": float(np.percentile(delta, 95)),
        "delta_max": float(delta.max()),
    }


# --------------------------------------------------------------------------
# C4: Holder exponent of the vorticity direction on the set
# --------------------------------------------------------------------------
def vorticity_direction(omega, eps_rel=1e-6):
    """xi = omega/|omega|, with a relative floor to avoid 0/0 at vorticity nulls."""
    mag = np.sqrt(sum(omega[i] ** 2 for i in range(3)))
    floor = eps_rel * float(mag.max()) if mag.max() > 0 else 1.0
    safe = np.maximum(mag, floor)
    return np.stack([omega[i] / safe for i in range(3)]), mag


def holder_beta_xi(xi, mask, h_cells, dirs=_DIRS13, n_samples=4000, rng=None):
    """Estimate beta from <|xi(x+h)-xi(x)|^2>_{x in mask}^{1/2} ~ h^beta.

    Returns (beta, intercept, h_cells, values, r2). Theorem C4 needs beta >= 1/2.
    """
    rng = np.random.default_rng(0 if rng is None else rng)
    N = xi.shape[1]
    pts = np.argwhere(mask)
    if pts.shape[0] < 16:
        return {"beta": float("nan"), "r2": float("nan"), "h_cells": list(map(float, h_cells)),
                "struct": [], "n_points": int(pts.shape[0])}
    if pts.shape[0] > n_samples:
        pts = pts[rng.choice(pts.shape[0], n_samples, replace=False)]
    xi0 = xi[:, pts[:, 0], pts[:, 1], pts[:, 2]]         # (3, P)
    vals = []
    for h in h_cells:
        acc = []
        for d in dirs:
            q = np.rint(pts + h * d).astype(np.int64) % N
            xih = xi[:, q[:, 0], q[:, 1], q[:, 2]]
            acc.append(np.mean(sum((xih[i] - xi0[i]) ** 2 for i in range(3))))
        vals.append(float(np.sqrt(np.mean(acc))))
    lx, ly = np.log(np.asarray(h_cells, float)), np.log(np.maximum(vals, 1e-300))
    A = np.vstack([lx, np.ones_like(lx)]).T
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    resid = ly - A @ coef
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = float(1 - np.sum(resid**2) / ss_tot) if ss_tot > 0 else float("nan")
    return {"beta": float(coef[0]), "intercept": float(coef[1]),
            "h_cells": list(map(float, h_cells)), "struct": vals, "r2": r2,
            "n_points": int(pts.shape[0])}


# --------------------------------------------------------------------------
# C5: local mean oscillation of xi, and its modulus
# --------------------------------------------------------------------------
def mean_oscillation_modulus(xi, mask, box_cells, dx=1.0, L=None):
    """mu_xi(r) = < (1/|B_r|) int_{B_r} |xi - <xi>_{B_r}| > over boxes meeting `mask`.

    Theorem C5 (A16) asks whether mu_xi(r) decays like 1/log(L/r) as r -> 0, i.e.
    whether mu_xi(r) * log(L/r) stays BOUNDED.

    The weight MUST be referenced to the outer scale L. Using |log r| instead is
    a trap: it vanishes at r = 1 and is non-monotone across it, so the product
    reports a smooth field as "growing" and a white-noise field as "decaying" --
    exactly backwards. log(L/r) is positive and monotone for r < L.

    Also returns the log-log slope of mu(r), which is the robust discriminator:
    slope ~ +1 for a Lipschitz direction field, ~0 for one with no regularity.
    """
    N = xi.shape[1]
    if L is None:
        L = N * dx
    out = []
    for b in box_cells:
        if N % b:
            continue
        m = N // b
        # block-reshape to (m,b,m,b,m,b) -> (m,m,m,b^3)
        def blocks(a):
            return (a.reshape(m, b, m, b, m, b)
                     .transpose(0, 2, 4, 1, 3, 5)
                     .reshape(m, m, m, b**3))
        xb = np.stack([blocks(xi[i]) for i in range(3)])       # (3,m,m,m,b^3)
        mb = blocks(mask)                                       # (m,m,m,b^3)
        sel = mb.any(axis=3)
        if not sel.any():
            continue
        mean = xb.mean(axis=4, keepdims=True)
        dev = np.sqrt(sum((xb[i] - mean[i]) ** 2 for i in range(3)))  # (m,m,m,b^3)
        osc = dev.mean(axis=3)[sel]
        r = b * dx
        mu = float(osc.mean())
        w = np.log(L / r) if 0 < r < L else float("nan")
        out.append({"box_cells": int(b), "r": float(r), "mu": mu,
                    "log_outer_over_r": float(w),
                    "mu_times_logL_over_r": float(mu * w),
                    "n_boxes": int(sel.sum())})
    if len(out) >= 3:
        rs = np.log([o["r"] for o in out])
        ms = np.log([max(o["mu"], 1e-300) for o in out])
        slope = float(np.polyfit(rs, ms, 1)[0])
        for o in out:
            o["mu_loglog_slope"] = slope
    return out


# --------------------------------------------------------------------------
# Secondary: dimensions (Prop. B4 -- diagnostic only, NOT theorem-relevant)
# --------------------------------------------------------------------------
def generalized_dimensions(field, box_cells, qs=(0, 1, 2, 4)):
    """D_q of the measure mu = field/sum(field) via box partition, plus D_inf.

    chi_q(l) = sum_B mu(B)^q ~ l^{(q-1) D_q}.
    """
    N = field.shape[0]
    tot = float(field.sum())
    if tot <= 0:
        return {"D": {}, "scales": []}
    scales, chi, mx = [], {q: [] for q in qs}, []
    for b in box_cells:
        if N % b:
            continue
        m = N // b
        mu = (field.reshape(m, b, m, b, m, b).sum(axis=(1, 3, 5)) / tot).ravel()
        mu = mu[mu > 0]
        scales.append(b / N)
        for q in qs:
            if q == 1:
                chi[q].append(float(-np.sum(mu * np.log(mu))))   # entropy
            else:
                chi[q].append(float(np.sum(mu**q)))
        mx.append(float(mu.max()))
    if len(scales) < 3:
        return {"D": {}, "scales": scales}
    ls = np.log(np.asarray(scales))
    D = {}
    for q in qs:
        y = np.asarray(chi[q])
        if q == 1:
            slope = np.polyfit(ls, y, 1)[0]
            D[1] = float(slope)                                  # S(l) ~ D_1 log l
        else:
            slope = np.polyfit(ls, np.log(np.maximum(y, 1e-300)), 1)[0]
            D[q] = float(slope / (q - 1))
    D["inf"] = float(np.polyfit(ls, np.log(np.maximum(mx, 1e-300)), 1)[0])
    return {"D": D, "scales": scales, "chi": {str(k): v for k, v in chi.items()},
            "max_mu": mx}


def volume_slope_dimension(field, capture, box_cells):
    """D from |E(l)| ~ l^{3-D} where E is the capture-threshold set of the
    l-coarse-grained field. Diagnostic only (Prop. B1)."""
    N = field.shape[0]
    vols, scales = [], []
    for b in box_cells:
        if N % b:
            continue
        m = N // b
        cg = field.reshape(m, b, m, b, m, b).mean(axis=(1, 3, 5))
        thr, vf = capture_threshold(cg, capture)
        if vf <= 0:
            continue
        vols.append(vf)
        scales.append(b / N)
    if len(scales) < 3:
        return {"D": float("nan"), "scales": scales, "vol_frac": vols}
    slope = np.polyfit(np.log(scales), np.log(vols), 1)[0]
    return {"D": float(3.0 - slope), "slope": float(slope),
            "scales": scales, "vol_frac": vols}


# --------------------------------------------------------------------------
# B2: the ratio that actually decides the balance
# --------------------------------------------------------------------------
def restricted_ratio(net_stretch, gradw2, mask, nu):
    """R_E = int_E omega.S.omega / (nu int_E |grad omega|^2).

    Unlike the global ratio (Cor. B2.1) this is NOT a restatement of Omega' > 0.
    """
    num = float(net_stretch[mask].sum())
    den = float(nu * gradw2[mask].sum())
    return num / den if den > 0 else float("inf")


def depletion_ratio(grid, uh, ell_cells):
    """rho(l) = int Gamma[u_l] / int |omega_l|^2 |S_l|  for the filtered field."""
    ell = ell_cells * grid.dx
    G = grid.gaussian_filter_hat(ell)
    uhf = np.stack([uh[i] * G for i in range(3)])
    om = np.stack([grid.ifft(c) for c in grid.curl(uhf)])
    net, gp, _ = stretching_density(grid, uhf, omega=om)
    Sf = strain_frobenius(grid, uhf)
    mag2 = sum(om[i] ** 2 for i in range(3))
    den = float((mag2 * Sf).sum())
    return float(gp.sum()) / den if den > 0 else float("nan")


def alignment_pdf(grid, uh, mask, nbins=40):
    """cos angle between xi and each strain eigenvector, conditioned on `mask`."""
    S = grid.strain(uh)
    om = np.stack([grid.ifft(c) for c in grid.curl(uh)])
    xi, _ = vorticity_direction(om)
    pts = np.argwhere(mask)
    if pts.shape[0] == 0:
        return {}
    if pts.shape[0] > 20000:
        pts = pts[np.random.default_rng(0).choice(pts.shape[0], 20000, replace=False)]
    Sp = S[:, :, pts[:, 0], pts[:, 1], pts[:, 2]].transpose(2, 0, 1)  # (P,3,3)
    xp = xi[:, pts[:, 0], pts[:, 1], pts[:, 2]].T                      # (P,3)
    evals, evecs = np.linalg.eigh(Sp)          # ascending: e1 <= e2 <= e3
    out = {}
    for j, name in enumerate(("e1_compress", "e2_intermediate", "e3_extensional")):
        c = np.abs(np.einsum("pi,pi->p", xp, evecs[:, :, j]))
        h, edges = np.histogram(c, bins=nbins, range=(0, 1), density=True)
        out[name] = {"hist": h.tolist(), "edges": edges.tolist(),
                     "mean_abs_cos": float(c.mean())}
    out["eigenvalue_means"] = [float(evals[:, j].mean()) for j in range(3)]
    return out


def sparseness_curve(mask, r_cells_list, n_samples=2000, rng=None):
    """delta(r) over a range of scales.

    The theorem-relevant value is delta at r ~ ell_nu. Sweeping r is mandatory:
    delta(r) is meaningless for r below the transverse feature size (every
    direction is then fully occupied) and saturates at the global volume
    fraction for r >> the object. Report the curve, never a single r.
    """
    return [sparseness_1d(mask, r, n_samples=n_samples, rng=rng) for r in r_cells_list]


def box_counting_dimension(mask, box_cells):
    """D_0 of a SET (indicator), via box counting.

    Must be applied to a thresholded mask, never to a smooth field: a smooth
    positive field has full support and yields D_0 = 3 identically.
    """
    N = mask.shape[0]
    counts, scales = [], []
    for b in box_cells:
        if N % b:
            continue
        m = N // b
        occ = mask.reshape(m, b, m, b, m, b).any(axis=(1, 3, 5))
        n = int(occ.sum())
        if n == 0:
            continue
        counts.append(n)
        scales.append(b / N)
    if len(scales) < 3:
        return {"D0": float("nan"), "scales": scales, "counts": counts}
    slope = np.polyfit(np.log(scales), np.log(counts), 1)[0]
    return {"D0": float(-slope), "scales": scales, "counts": counts}


# --------------------------------------------------------------------------
# C3, corrected: Grujic's criterion is stated on the super-level sets of the
# POSITIVE AND NEGATIVE PARTS OF THE VORTICITY COMPONENTS, not of |omega|.
# --------------------------------------------------------------------------
def component_superlevel_sets(omega, lam=0.5):
    """The six sets {omega_i^{+/-} > lam * ||omega_i^{+/-}||_inf}, i = 1,2,3.

    Grujic (2013) and Bradshaw-Farhat-Grujic (2019) define the scale of
    sparseness on the super-level sets of the positive and negative parts of the
    vorticity COMPONENTS. Using {|omega| > lam ||omega||_inf} instead is a
    different (and in general strictly larger, less anisotropic) family of sets,
    and a sparseness statement about it is NOT the hypothesis of the theorem.

    Returns a list of (label, mask, sup) triples, one per component and sign.
    """
    out = []
    for i, axis in enumerate("xyz"):
        for sgn, tag in ((+1.0, "+"), (-1.0, "-")):
            part = np.maximum(sgn * omega[i], 0.0)
            sup = float(part.max())
            if sup <= 0:
                continue
            out.append((f"w{axis}{tag}", part > lam * sup, sup))
    return out


def sparseness_grujic(omega, r_cells, lam=0.5, n_samples=2000, rng=None):
    """1D sparseness of ALL SIX component super-level sets (the C3 hypothesis).

    The criterion requires every one of the six sets to be thin, so the
    reportable quantity is the WORST (largest) delta over them. Reporting only
    the mean would hide the set that fails.
    """
    sets = component_superlevel_sets(omega, lam=lam)
    per = {}
    for label, mask, sup in sets:
        st = sparseness_1d(mask, r_cells, n_samples=n_samples, rng=rng)
        st["sup"] = sup
        st["volume_fraction"] = float(mask.mean())
        per[label] = st
    if not per:
        return {"per_set": {}, "delta_worst_p95": float("nan"), "lam": lam,
                "r_cells": float(r_cells)}
    worst = max(per.values(), key=lambda q: (q["delta_p95"] if q["delta_p95"] == q["delta_p95"] else -1))
    worst_label = [k for k, v in per.items() if v is worst][0]
    return {
        "per_set": per,
        "lam": lam,
        "r_cells": float(r_cells),
        "delta_worst_p95": float(worst["delta_p95"]),
        "delta_worst_set": worst_label,
        "delta_mean_over_sets": float(np.mean([v["delta_p95"] for v in per.values()])),
    }
