"""unittest wrappers over the core invariants. `python3 -m unittest discover tests`"""
import os, sys, unittest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from nsref.grid import Grid                                   # noqa: E402
from nsref import ic, gamma, diagnostics as dg, fit           # noqa: E402
from nsref.solver import NSSolver                             # noqa: E402


class TestSpectral(unittest.TestCase):
    def setUp(self):
        self.g = Grid(32)

    def test_roundtrip(self):
        f = np.random.default_rng(0).random((32,) * 3)
        self.assertLess(np.abs(self.g.ifft(self.g.fft(f)) - f).max(), 1e-12)

    def test_divergence_free(self):
        uh = np.stack([self.g.fft(c) for c in ic.taylor_green(self.g)])
        self.assertLess(self.g.divergence_max(uh), 1e-12)

    def test_parseval(self):
        u = ic.taylor_green(self.g)
        uh = np.stack([self.g.fft(c) for c in u])
        direct = 0.5 * np.mean(sum(c**2 for c in u))
        self.assertAlmostEqual(self.g.shell_spectrum(uh).sum(), direct, places=12)

    def test_dealias_mask_zeroes(self):
        uh = np.stack([self.g.fft(c) for c in ic.abc_flow(self.g)])
        self.assertEqual(self.g.dealias_residual(self.g.apply_dealias(uh[0])), 0.0)


class TestSolver(unittest.TestCase):
    def test_taylor_green_analytic_initial_values(self):
        g = Grid(32)
        d = dg.basic(g, np.stack([g.fft(c) for c in ic.taylor_green(g)]), 0.05)
        self.assertAlmostEqual(d["E"], 0.125, places=12)
        self.assertAlmostEqual(d["enstrophy"], 0.375, places=12)
        self.assertAlmostEqual(d["omega_max"], 2.0, places=9)

    def test_budgets_and_monotone_energy(self):
        g = Grid(32)
        nu = 0.05
        s = NSSolver(g, nu, ic.taylor_green(g))
        prev = dg.budget_terms(g, s.vh, nu)
        we = wo = 0.0
        for _ in range(10):
            dt = s.step()
            b = dg.budget_terms(g, s.vh, nu)
            we = max(we, dg.energy_budget_residual(prev["E"], b["E"], dt,
                                                   0.5 * (b["eps"] + prev["eps"])))
            wo = max(wo, dg.enstrophy_budget_residual(
                prev["enstrophy"], b["enstrophy"], dt,
                0.5 * (b["stretch_net"] + prev["stretch_net"]),
                0.5 * (b["visc_dissip_enstrophy"] + prev["visc_dissip_enstrophy"])))
            self.assertLessEqual(b["E"], prev["E"] + 1e-12)
            prev = b
        self.assertLess(we, 1e-3)
        self.assertLess(wo, 1e-2)
        self.assertEqual(s.max_dealias_residual, 0.0)


class TestGamma(unittest.TestCase):
    def test_known_geometry_ordering(self):
        g = Grid(64)
        d = {}
        for name, (fn, _) in ic.SYNTH.items():
            d[name] = gamma.generalized_dimensions(fn(g), [2, 4, 8, 16])["D"]["inf"]
        self.assertLess(d["blob"], d["tube"])
        self.assertLess(d["tube"], d["sheet"])
        self.assertLess(d["sheet"], d["uniform"])

    def test_sparseness_discriminates(self):
        g = Grid(64)
        out = {}
        for name in ("sheet", "tube", "uniform"):
            f = ic.SYNTH[name][0](g)
            thr, _ = gamma.capture_threshold(f, 0.9)
            out[name] = gamma.sparseness_1d(gamma.superlevel_mask(f, thr), 16,
                                            n_samples=800)["delta_p95"]
        self.assertLess(out["tube"], out["uniform"])
        self.assertLess(out["sheet"], out["uniform"])

    def test_D0_requires_thresholding(self):
        """A smooth positive field has full support: D_0 on it is 3 identically."""
        g = Grid(64)
        f = ic.SYNTH["tube"][0](g)
        full = gamma.box_counting_dimension(f > 0, [2, 4, 8, 16])["D0"]
        self.assertAlmostEqual(full, 3.0, places=2)
        thr, _ = gamma.capture_threshold(f, 0.9)
        masked = gamma.box_counting_dimension(f > thr, [2, 4, 8, 16])["D0"]
        self.assertLess(masked, 2.0)

    def test_stretching_identity(self):
        """omega_i omega_j d_j u_i must equal omega . S omega."""
        g = Grid(32)
        uh = np.stack([g.fft(c) for c in ic.abc_flow(g)])
        net, _, _ = gamma.stretching_density(g, uh)
        S = g.strain(uh)
        om = np.stack([g.ifft(c) for c in g.curl(uh)])
        direct = sum(om[i] * S[i, j] * om[j] for i in range(3) for j in range(3))
        self.assertLess(np.abs(net - direct).max() / max(np.abs(direct).max(), 1e-30), 1e-10)


class TestFit(unittest.TestCase):
    def test_power_law_recovered(self):
        rng = np.random.default_rng(3)
        t = np.linspace(0, 1.8, 60)
        w = (2.0 - t) ** (-1.6) * np.exp(rng.normal(0, 0.01, t.size))
        self.assertTrue(fit.compare(t, w)["power_law_preferred"])
        b = fit.bootstrap_powerlaw(t, w, B=200)
        self.assertLess(abs(b["Tstar"] - 2.0), 0.02)

    def test_exponential_not_called_power_law(self):
        rng = np.random.default_rng(4)
        t = np.linspace(0, 1.8, 60)
        w = np.exp(1.4 * t) * np.exp(rng.normal(0, 0.01, t.size))
        self.assertFalse(fit.compare(t, w)["power_law_preferred"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
