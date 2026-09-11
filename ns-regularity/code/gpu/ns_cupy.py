"""CuPy port of the reference solver for the 4-GPU campaign.

Algorithmically IDENTICAL to nsref/solver.py -- rotational form, 2/3 dealiasing,
integrating-factor RK4, FP64 -- so that a CPU run at 64^3/128^3 and a GPU run at
the same parameters must agree to round-off. That equivalence is the first thing
to check on the workstation, before any production run:

    python3 -m gpu.ns_cupy --selftest --N 64 --nu 0.02 --steps 20

Not runnable in the drafting environment (no CUDA). Import is guarded.

Memory at 512^3 FP64: one real field 1.07 GB, one complex spectral field 1.08 GB.
The state is 3 spectral fields; RK4 holds 4 stage arrays plus 2 scratch; the
nonlinear evaluation holds 6 real fields. Budget ~22-30 GB including the cuFFT
workspace. Single GPU if VRAM >= 40 GB; otherwise run 512^3 in FP32 and make
PRECISION_ESCALATION_CONSISTENT mandatory at 256^3.
"""
from __future__ import annotations
import argparse

try:
    import cupy as xp
    import cupyx.scipy.fft as _cufft          # noqa: F401
    HAVE_CUPY = True
except Exception:                              # pragma: no cover
    import numpy as xp
    HAVE_CUPY = False


class GPUGrid:
    def __init__(self, N, L=2 * xp.pi if HAVE_CUPY else 6.283185307179586, dtype="float64"):
        self.N, self.L = N, float(L)
        self.dtype = xp.float64 if dtype == "float64" else xp.float32
        self.cdtype = xp.complex128 if dtype == "float64" else xp.complex64
        self.dx = self.L / N
        k1 = xp.fft.fftfreq(N, d=1.0 / N) * (2 * xp.pi / self.L)
        kz = xp.fft.rfftfreq(N, d=1.0 / N) * (2 * xp.pi / self.L)
        self.kx, self.ky, self.kz = k1[:, None, None], k1[None, :, None], kz[None, None, :]
        self.k2 = self.kx**2 + self.ky**2 + self.kz**2
        self.k2_nz = xp.where(self.k2 == 0, 1.0, self.k2)
        self.kmax = N // 2 * (2 * xp.pi / self.L)
        kcut = (N // 3) * (2 * xp.pi / self.L)
        self.dealias = ((abs(self.kx) <= kcut) & (abs(self.ky) <= kcut) & (abs(self.kz) <= kcut))

    def fft(self, f):
        return xp.fft.rfftn(f, axes=(0, 1, 2))

    def ifft(self, fh):
        return xp.fft.irfftn(fh, s=(self.N,) * 3, axes=(0, 1, 2))

    def project(self, vh):
        f = (self.kx * vh[0] + self.ky * vh[1] + self.kz * vh[2]) / self.k2_nz
        return xp.stack([vh[0] - self.kx * f, vh[1] - self.ky * f, vh[2] - self.kz * f])

    def curl(self, vh):
        return xp.stack([1j * (self.ky * vh[2] - self.kz * vh[1]),
                         1j * (self.kz * vh[0] - self.kx * vh[2]),
                         1j * (self.kx * vh[1] - self.ky * vh[0])])

    def dealias_residual(self, vh):
        return float(abs(vh * (~self.dealias)).max())


class GPUSolver:
    def __init__(self, grid, nu, u0, cfl=0.4):
        self.g, self.nu, self.cfl = grid, float(nu), float(cfl)
        self.t, self.step_count, self.max_dealias_residual = 0.0, 0, 0.0
        vh = xp.stack([grid.fft(c) for c in u0])
        self.vh = grid.project(vh) * grid.dealias
        self._L = -self.nu * grid.k2

    def _nonlinear(self, vh):
        g = self.g
        u = xp.stack([g.ifft(c) for c in vh])
        w = xp.stack([g.ifft(c) for c in g.curl(vh)])
        cx = xp.stack([u[1] * w[2] - u[2] * w[1],
                       u[2] * w[0] - u[0] * w[2],
                       u[0] * w[1] - u[1] * w[0]])
        return g.project(xp.stack([g.fft(c) for c in cx]) * g.dealias)

    def dt(self):
        u = xp.stack([self.g.ifft(c) for c in self.vh])
        return float(self.cfl * self.g.dx / (3**0.5 * (float(abs(u).max()) + 1e-30)))

    def step(self, dt=None):
        g = self.g
        dt = self.dt() if dt is None else dt
        E, E2 = xp.exp(self._L * dt), xp.exp(self._L * dt / 2)
        v = self.vh
        k1 = dt * self._nonlinear(v)
        k2 = dt * self._nonlinear(E2 * (v + 0.5 * k1))
        k3 = dt * self._nonlinear(E2 * v + 0.5 * k2)
        k4 = dt * self._nonlinear(E * v + E2 * k3)
        self.vh = g.project((E * v + (E * k1 + 2 * E2 * (k2 + k3) + k4) / 6) * g.dealias)
        self.max_dealias_residual = max(self.max_dealias_residual, g.dealias_residual(self.vh))
        self.t += dt
        self.step_count += 1
        return dt


def selftest(N=64, nu=0.02, steps=20):
    """CPU/GPU equivalence: identical algorithm must give identical trajectories."""
    import os, sys
    import numpy as np
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    from nsref.grid import Grid
    from nsref import ic
    from nsref.solver import NSSolver

    gc = Grid(N)
    u0 = ic.taylor_green(gc)
    sc = NSSolver(gc, nu, u0, cfl=0.4)
    gg = GPUGrid(N)
    sg = GPUSolver(gg, nu, xp.asarray(u0) if HAVE_CUPY else u0, cfl=0.4)
    for _ in range(steps):
        dt = sc.dt()
        sc.step(dt)
        sg.step(dt)
    a = np.asarray(sc.vh)
    b = xp.asnumpy(sg.vh) if HAVE_CUPY else np.asarray(sg.vh)
    rel = np.abs(a - b).max() / max(np.abs(a).max(), 1e-300)
    print(f"backend={'cupy' if HAVE_CUPY else 'numpy(fallback)'}  N={N} steps={steps}")
    print(f"max relative difference CPU vs GPU: {rel:.3e}   "
          f"{'PASS' if rel < 1e-10 else 'FAIL'}")
    return rel < 1e-10


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--selftest", action="store_true")
    p.add_argument("--N", type=int, default=64)
    p.add_argument("--nu", type=float, default=0.02)
    p.add_argument("--steps", type=int, default=20)
    a = p.parse_args()
    if a.selftest:
        raise SystemExit(0 if selftest(a.N, a.nu, a.steps) else 1)
    print(__doc__)
