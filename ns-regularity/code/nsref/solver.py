"""Pseudospectral Navier-Stokes solver, rotational form, integrating-factor RK4.

Deliberate choices (see spec/campaign_4gpu.md):
  * rotational form P[u x omega]  -- conserves energy better than the divergence form
  * 2/3 dealiasing applied to every nonlinear product, asserted each step
  * exact integrating factor for the viscous term (no viscous time-step limit)
  * RK4, not RK2: the RK2 used in the originating runs has O(dt^2) phase error
    that corrupts ||omega||_inf timeseries and therefore the T* fits
  * NO hyperviscosity, ever.
"""
import numpy as np


class NSSolver:
    def __init__(self, grid, nu, u0, cfl=0.4, dt_max=None):
        self.g, self.nu, self.cfl = grid, float(nu), float(cfl)
        self.dt_max = dt_max
        self.t = 0.0
        self.step_count = 0
        uh = np.stack([grid.fft(c) for c in u0])
        self.vh = grid.apply_dealias(grid.project(uh))
        self._Lop = -self.nu * grid.k2
        self._dt_cached = None
        self.max_dealias_residual = 0.0

    # ---- right-hand side -------------------------------------------------
    def _nonlinear(self, vh):
        g = self.g
        u = np.stack([g.ifft(c) for c in vh])
        w = np.stack([g.ifft(c) for c in g.curl(vh)])
        cx = np.stack([u[1] * w[2] - u[2] * w[1],
                       u[2] * w[0] - u[0] * w[2],
                       u[0] * w[1] - u[1] * w[0]])
        nh = np.stack([g.fft(c) for c in cx])
        nh = g.apply_dealias(nh)
        return g.project(nh)

    def velocity(self):
        return np.stack([self.g.ifft(c) for c in self.vh])

    def vorticity(self):
        return np.stack([self.g.ifft(c) for c in self.g.curl(self.vh)])

    # ---- time step -------------------------------------------------------
    def dt(self):
        u = self.velocity()
        umax = float(np.max(np.abs(u))) + 1e-30
        dt = self.cfl * self.g.dx / (np.sqrt(3.0) * umax)
        if self.dt_max is not None:
            dt = min(dt, self.dt_max)
        return dt

    def cfl_number(self, dt):
        u = self.velocity()
        return float(np.sqrt(3.0) * np.max(np.abs(u)) * dt / self.g.dx)

    def step(self, dt=None):
        """One integrating-factor RK4 step. Returns dt actually taken."""
        g = self.g
        dt = self.dt() if dt is None else dt
        E = np.exp(self._Lop * dt)
        E2 = np.exp(self._Lop * dt / 2.0)
        v = self.vh
        k1 = dt * self._nonlinear(v)
        k2 = dt * self._nonlinear(E2 * (v + 0.5 * k1))
        k3 = dt * self._nonlinear(E2 * v + 0.5 * k2)
        k4 = dt * self._nonlinear(E * v + E2 * k3)
        self.vh = E * v + (E * k1 + 2.0 * E2 * (k2 + k3) + k4) / 6.0
        self.vh = g.apply_dealias(g.project(self.vh))
        # DEALIASING_ACTIVE evidence
        r = max(g.dealias_residual(self.vh[i]) for i in range(3))
        self.max_dealias_residual = max(self.max_dealias_residual, r)
        self.t += dt
        self.step_count += 1
        if not np.isfinite(np.abs(self.vh).max()):
            raise FloatingPointError(f"non-finite state at t={self.t}, step {self.step_count}")
        return dt
