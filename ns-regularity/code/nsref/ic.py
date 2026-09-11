"""Initial conditions, and synthetic fields of KNOWN geometry for validating
the estimators in gamma.py.

Dynamic ICs return a divergence-free velocity in physical space, shape (3,N,N,N).
Synthetic analysis fields return a scalar density with known dimension D.
"""
import numpy as np


def _mesh(grid):
    x = np.arange(grid.N) * grid.dx
    return np.meshgrid(x, x, x, indexing="ij")


def _project(grid, u):
    uh = np.stack([grid.fft(c) for c in u])
    uh = grid.apply_dealias(grid.project(uh))
    return np.stack([grid.ifft(c) for c in uh])


# --------------------------------------------------------------------------
# dynamic initial conditions
# --------------------------------------------------------------------------
def taylor_green(grid, amp=1.0):
    X, Y, Z = _mesh(grid)
    u = np.stack([amp * np.sin(X) * np.cos(Y) * np.cos(Z),
                  -amp * np.cos(X) * np.sin(Y) * np.cos(Z),
                  np.zeros_like(X)])
    return _project(grid, u)


def abc_flow(grid, A=1.0, B=1.0, C=1.0):
    """Arnold-Beltrami-Childress: an exact Beltrami (fully helical) field."""
    X, Y, Z = _mesh(grid)
    u = np.stack([A * np.sin(Z) + C * np.cos(Y),
                  B * np.sin(X) + A * np.cos(Z),
                  C * np.sin(Y) + B * np.cos(X)])
    return _project(grid, u)


def antiparallel_tubes(grid, sep=1.6, core=0.35, amp=1.0, perturb=0.12, kpert=1.0):
    """Two counter-rotating Gaussian vortex tubes along x, sinusoidally perturbed.

    The Kerr-type geometry. Vorticity is imposed and the velocity recovered by
    Biot-Savart (curl of the inverse-Laplacian of omega), then Leray-projected.
    """
    X, Y, Z = _mesh(grid)
    L = grid.L
    yc, zc = L / 2, L / 2
    dz = perturb * np.sin(kpert * X)
    w = np.zeros((3,) + X.shape)
    for sgn, off in ((+1.0, +sep / 2), (-1.0, -sep / 2)):
        dy = (Y - (yc + off) + L / 2) % L - L / 2
        dzz = (Z - (zc + sgn * dz) + L / 2) % L - L / 2
        w[0] += sgn * amp * np.exp(-(dy**2 + dzz**2) / core**2)
    wh = np.stack([grid.fft(c) for c in w])
    wh = grid.project(wh)                      # enforce div omega = 0
    # u_hat = i k x omega_hat / k^2
    uh = np.empty_like(wh)
    kx, ky, kz = grid.kx, grid.ky, grid.kz
    uh[0] = 1j * (ky * wh[2] - kz * wh[1]) / grid.k2_nz
    uh[1] = 1j * (kz * wh[0] - kx * wh[2]) / grid.k2_nz
    uh[2] = 1j * (kx * wh[1] - ky * wh[0]) / grid.k2_nz
    uh[:, 0, 0, 0] = 0
    uh = grid.apply_dealias(grid.project(uh))
    return np.stack([grid.ifft(c) for c in uh])


def vortex_sheet(grid, thickness=0.3, amp=1.0, perturb=0.1):
    """Smoothed shear layer (tanh profile) with a transverse perturbation."""
    X, Y, Z = _mesh(grid)
    L = grid.L
    yy = (Y - L / 2) / thickness
    u = np.stack([amp * np.tanh(yy),
                  perturb * np.sin(X) * np.exp(-yy**2),
                  perturb * np.sin(Z) * np.exp(-yy**2)])
    return _project(grid, u)


def multiscale_random(grid, seed=0, amp=1.0, kpeak=3.0, slope=-5.0 / 3.0, kmin=1):
    """Divergence-free random field with a prescribed multiscale spectrum.

    Regression case for the originating run's strongest candidate (ADV023).
    """
    rng = np.random.default_rng(seed)
    N = grid.N
    shape = (3, N, N, N // 2 + 1)
    ph = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    k = np.sqrt(grid.k2)
    env = np.where(k >= kmin, (np.maximum(k, 1e-12) / kpeak) ** (slope / 2.0), 0.0)
    env = np.where(k == 0, 0.0, env)
    uh = grid.apply_dealias(grid.project(ph * env))
    u = np.stack([grid.ifft(c) for c in uh])
    u *= amp / np.sqrt(np.mean(sum(c**2 for c in u)))
    return _project(grid, u)


REGISTRY = {
    "taylor_green": taylor_green,
    "abc": abc_flow,
    "antiparallel_tubes": antiparallel_tubes,
    "vortex_sheet": vortex_sheet,
    "multiscale_random": multiscale_random,
}


# --------------------------------------------------------------------------
# synthetic fields of KNOWN geometry (analysis validation only -- not flows)
# --------------------------------------------------------------------------
def synth_blob(grid, width=0.25):
    """Isolated Gaussian blob: box-counting dimension D = 0."""
    X, Y, Z = _mesh(grid)
    L = grid.L
    d2 = sum(((q - L / 2 + L / 2) % L - L / 2) ** 2 for q in (X, Y, Z))
    return np.exp(-d2 / width**2)


def synth_tube(grid, radius=0.18):
    """Straight tube along x: D = 1, and 1D-sparse transversally."""
    X, Y, Z = _mesh(grid)
    L = grid.L
    dy = (Y - L / 2 + L / 2) % L - L / 2
    dz = (Z - L / 2 + L / 2) % L - L / 2
    return np.exp(-(dy**2 + dz**2) / radius**2)


def synth_sheet(grid, thickness=0.18):
    """Plane sheet normal to y: D = 2, still 1D-sparse along y."""
    _, Y, _ = _mesh(grid)
    L = grid.L
    dy = (Y - L / 2 + L / 2) % L - L / 2
    return np.exp(-(dy / thickness) ** 2)


def synth_uniform(grid, seed=0):
    """Space-filling random field: D = 3, NOT sparse in any direction."""
    rng = np.random.default_rng(seed)
    return rng.random((grid.N,) * 3)


SYNTH = {"blob": (synth_blob, 0.0), "tube": (synth_tube, 1.0),
         "sheet": (synth_sheet, 2.0), "uniform": (synth_uniform, 3.0)}
