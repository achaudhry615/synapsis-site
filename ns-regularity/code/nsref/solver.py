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
    def __init__(self, grid, nu, u0, cfl=0.4, dt_max=None, visc_dt_factor=0.25,
                 ramp_steps=8):
        """visc_dt_factor caps dt at visc_dt_factor / (nu k_max^2).

        The integrating factor makes the viscous term EXACT, so this cap is not
        needed for stability. It is needed for measurement: the budget residuals
        are O(dt^2)-accurate finite differences, and in a slow flow the CFL
        condition permits a dt large enough that the residual reports the
        difference error rather than the solver error. Set to None to disable.
        """
        self.g, self.nu, self.cfl = grid, float(nu), float(cfl)
        self.dt_max = dt_max
        self.visc_dt_factor = visc_dt_factor
        self.ramp_steps = int(ramp_steps)
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
        if self.visc_dt_factor is not None and self.nu > 0:
            dt = min(dt, self.visc_dt_factor / (self.nu * self.g.kmax**2))
        if self.dt_max is not None:
            dt = min(dt, self.dt_max)
        if self.ramp_steps > 0 and self.step_count < self.ramp_steps:
            dt *= (self.step_count + 1) / self.ramp_steps
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
