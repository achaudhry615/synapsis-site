"""Spectral grid: wavenumbers, dealiasing mask, Leray projector."""
import numpy as np


class Grid:
    """Uniform periodic grid on [0, L)^3 with rfft layout (N, N, N//2+1)."""

    def __init__(self, N, L=2 * np.pi, dtype=np.float64):
        if N % 2:
            raise ValueError("N must be even")
        self.N, self.L, self.dtype = N, float(L), dtype
        self.cdtype = np.complex128 if dtype == np.float64 else np.complex64
        self.dx = self.L / N

        k1 = np.fft.fftfreq(N, d=1.0 / N) * (2 * np.pi / self.L)
        kz = np.fft.rfftfreq(N, d=1.0 / N) * (2 * np.pi / self.L)
        self.kx = k1[:, None, None]
        self.ky = k1[None, :, None]
        self.kz = kz[None, None, :]
        self.k2 = self.kx**2 + self.ky**2 + self.kz**2
        self.k2_nz = np.where(self.k2 == 0, 1.0, self.k2)      # safe divisor
        self.kmag = np.sqrt(self.k2)

        # 2/3 rule: retain |k_i| <= kmax_dealias along every axis
        self.kmax = N // 2 * (2 * np.pi / self.L)
        kcut = (N // 3) * (2 * np.pi / self.L)
        self.dealias = ((np.abs(self.kx) <= kcut)
                        & (np.abs(self.ky) <= kcut)
                        & (np.abs(self.kz) <= kcut))

        # Hermitian multiplicity for correct rfft-space sums (Parseval)
        w = np.full((N, N, N // 2 + 1), 2.0)
        w[:, :, 0] = 1.0
        if N % 2 == 0:
            w[:, :, -1] = 1.0
        self.herm_w = w

        # integer shell index for spectra
        self.shell = np.rint(self.kmag * self.L / (2 * np.pi)).astype(np.int32)
        self.nshell = int(self.shell.max()) + 1

    # --- transforms -------------------------------------------------------
    def fft(self, f):
        return np.fft.rfftn(f, axes=(0, 1, 2))

    def ifft(self, fh):
        return np.fft.irfftn(fh, s=(self.N,) * 3, axes=(0, 1, 2))

    def apply_dealias(self, vh):
        return vh * self.dealias

    def dealias_residual(self, vh):
        """Max modulus of any mode that the 2/3 rule should have zeroed.

        Evidence for the DEALIASING_ACTIVE gate.
        """
        masked = vh * (~self.dealias)
        return float(np.abs(masked).max()) if masked.size else 0.0

    # --- vector calculus in spectral space --------------------------------
    def project(self, vh):
        """Leray projection onto divergence-free fields."""
        kdotv = self.kx * vh[0] + self.ky * vh[1] + self.kz * vh[2]
        f = kdotv / self.k2_nz
        out = np.empty_like(vh)
        out[0] = vh[0] - self.kx * f
        out[1] = vh[1] - self.ky * f
        out[2] = vh[2] - self.kz * f
        return out

    def divergence_max(self, vh):
        kdotv = 1j * (self.kx * vh[0] + self.ky * vh[1] + self.kz * vh[2])
        return float(np.abs(self.ifft(kdotv)).max())

    def curl(self, vh):
        out = np.empty_like(vh)
        out[0] = 1j * (self.ky * vh[2] - self.kz * vh[1])
        out[1] = 1j * (self.kz * vh[0] - self.kx * vh[2])
        out[2] = 1j * (self.kx * vh[1] - self.ky * vh[0])
        return out

    def grad(self, fh):
        return np.stack([1j * self.kx * fh, 1j * self.ky * fh, 1j * self.kz * fh])

    def strain(self, vh):
        """Symmetric strain-rate tensor S_ij in physical space, shape (3,3,N,N,N)."""
        k = (self.kx, self.ky, self.kz)
        S = np.empty((3, 3) + (self.N,) * 3, dtype=self.dtype)
        du = [[self.ifft(1j * k[j] * vh[i]) for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3):
                S[i, j] = 0.5 * (du[i][j] + du[j][i])
        return S

    def gaussian_filter_hat(self, ell):
        """Spectral Gaussian filter kernel at scale ell (Germano convention)."""
        return np.exp(-self.k2 * ell**2 / 24.0)

    def shell_spectrum(self, vh):
        """E(k) = 0.5 * sum_{|k| in shell} |v_hat|^2, normalised so sum E(k) = 0.5<|v|^2>."""
        e = 0.5 * self.herm_w * sum(np.abs(vh[i]) ** 2 for i in range(3))
        e = e / self.N**6
        return np.bincount(self.shell.ravel(), weights=e.ravel(), minlength=self.nshell)

    def mean(self, f):
        return float(f.mean())
