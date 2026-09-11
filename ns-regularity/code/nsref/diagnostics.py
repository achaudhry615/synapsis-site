"""Scalar diagnostics. All volume-averaged (per unit volume), not integrated."""
import numpy as np
from . import gamma as gmod


def omega_max_refined(grid, uh, pad=2):
    """||omega||_inf via zero-padded spectral interpolation.

    The grid maximum hops between cells as a vortex translates, injecting
    O(1%) noise that destroys nonlinear T* fits. Padding by 2 removes it at
    the cost of one transform on a (pad*N)^3 grid.
    """
    wh = grid.curl(uh)
    N, M = grid.N, pad * grid.N
    big = np.zeros((3, M, M, M // 2 + 1), dtype=wh.dtype)
    h = N // 2
    big[:, :h, :h, :h + 1] = wh[:, :h, :h, :h + 1]
    big[:, :h, -h:, :h + 1] = wh[:, :h, h:, :h + 1]
    big[:, -h:, :h, :h + 1] = wh[:, h:, :h, :h + 1]
    big[:, -h:, -h:, :h + 1] = wh[:, h:, h:, :h + 1]
    w = np.stack([np.fft.irfftn(big[i], s=(M, M, M), axes=(0, 1, 2)) for i in range(3)])
    w *= pad**3                                     # rfftn normalisation on the larger grid
    return float(np.sqrt(sum(w[i] ** 2 for i in range(3))).max())


def basic(grid, uh, nu, refine=True):
    u = np.stack([grid.ifft(c) for c in uh])
    w = np.stack([grid.ifft(c) for c in grid.curl(uh)])
    u2 = sum(u[i] ** 2 for i in range(3))
    w2 = sum(w[i] ** 2 for i in range(3))
    E = 0.5 * float(u2.mean())
    Om = 0.5 * float(w2.mean())
    eps = 2.0 * nu * Om                              # dissipation rate per unit mass
    eta = (nu**3 / eps) ** 0.25 if eps > 0 else np.inf
    wmax_grid = float(np.sqrt(w2).max())
    wmax = omega_max_refined(grid, uh) if refine else wmax_grid
    Y = float(gmod.grad_omega_sq(grid, uh).mean())
    net, gp, gm = gmod.stretching_density(grid, uh, omega=w)
    L3 = float((u2 ** 1.5).mean() ** (1.0 / 3.0))
    return {
        "E": E, "enstrophy": Om, "palinstrophy": Y,
        "eps": eps, "eta": eta, "kmax_eta": float(grid.kmax * eta),
        "omega_max": wmax, "omega_max_grid": wmax_grid,
        "omega_max_refine_gain": wmax / wmax_grid if wmax_grid > 0 else float("nan"),
        "ell_nu": float(np.sqrt(nu / wmax)) if wmax > 0 else float("inf"),
        "stretch_net": float(net.mean()), "gamma_plus": float(gp.mean()),
        "gamma_minus": float(gm.mean()),
        "visc_dissip_enstrophy": float(nu * Y),
        "R_global": float(net.mean() / (nu * Y)) if Y > 0 else float("inf"),
        "L3_norm": L3,
        "div_max": grid.divergence_max(uh),
    }


def budget_terms(grid, uh, nu):
    """Cheap per-step budget terms: E, enstrophy, stretching, viscous dissipation.

    Budget residuals are O(interval^2) accurate, so they must be evaluated on
    EVERY step, not at the logging cadence -- a midpoint difference over 5 steps
    reports ~2% error that is an artefact of the measurement, not the solver.
    Skips the padded omega_max refinement and the L3 norm.
    """
    u = np.stack([grid.ifft(c) for c in uh])
    w = np.stack([grid.ifft(c) for c in grid.curl(uh)])
    E = 0.5 * float(sum(u[i] ** 2 for i in range(3)).mean())
    Om = 0.5 * float(sum(w[i] ** 2 for i in range(3)).mean())
    net, _, _ = gmod.stretching_density(grid, uh, omega=w)
    Y = float(gmod.grad_omega_sq(grid, uh).mean())
    return {"E": E, "enstrophy": Om, "eps": 2.0 * nu * Om,
            "stretch_net": float(net.mean()), "visc_dissip_enstrophy": float(nu * Y)}


def energy_budget_residual(E_prev, E_now, dt, eps_mid):
    """|dE/dt + eps| / eps  -- the ENERGY_BALANCE_CONSISTENT gate quantity."""
    if eps_mid <= 0 or dt <= 0:
        return float("nan")
    return abs((E_now - E_prev) / dt + eps_mid) / eps_mid


def enstrophy_budget_residual(Om_prev, Om_now, dt, stretch_mid, visc_mid):
    """Residual of Omega' = <omega.S.omega> - nu<|grad omega|^2> (identity A3)."""
    if dt <= 0:
        return float("nan")
    lhs = (Om_now - Om_prev) / dt
    rhs = stretch_mid - visc_mid
    # Normalise by the magnitude of the individual terms, NOT by |lhs| or |rhs|:
    # Omega' crosses zero at the enstrophy peak, where a difference-normalised
    # residual diverges for purely algebraic reasons and reports a solver error
    # that does not exist.
    scale = max(abs(stretch_mid) + abs(visc_mid), 1e-300)
    return abs(lhs - rhs) / scale


def cheskidov_shvydkoy_number(grid, uh, nu):
    """sup_q 2^{-q} ||Delta_q u||_inf / nu over dyadic bands (proxy for A15).

    Replaces the structure-function exponents: unlike zeta_p, this quantity is
    attached to an actual regularity criterion.
    """
    kmag = grid.kmag
    out, best = [], 0.0
    q = 0
    while 2**q <= grid.kmax:
        band = (kmag >= 2**q) & (kmag < 2 ** (q + 1))
        if band.any():
            vb = np.stack([uh[i] * band for i in range(3)])
            ub = np.stack([grid.ifft(c) for c in vb])
            linf = float(np.sqrt(sum(ub[i] ** 2 for i in range(3))).max())
            val = (2.0 ** (-q)) * linf / nu
            out.append({"q": q, "linf": linf, "value": val})
            best = max(best, val)
        q += 1
    return {"cs_number": best, "bands": out}
