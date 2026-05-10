"""Discrete-population PFL model: replaces the integral with a sum over N neurons."""

from __future__ import annotations

import numpy as np

from steering.models.continuous import ContinuousPFLModel
from steering.nonlinearities import get_nonlinearity
from steering.params import ModelParams


class DiscretePFLModel(ContinuousPFLModel):
    """Population sum over ``N`` uniformly-spaced preferred headings.

    ``N`` is taken from ``params.N_neurons``. With ``N -> infinity`` this
    recovers :class:`ContinuousPFLModel`.
    """

    def __init__(self, N_default: int = 8):
        super().__init__(n_quad=2)  # quadrature unused, but keep parent happy
        self.N_default = N_default

    def _phi_grid(self, params: ModelParams) -> np.ndarray:
        N = params.N_neurons or self.N_default
        # Symmetric grid with N points on [-pi, pi). Note: the grid -pi + 2*pi*j/N
        # is symmetric about 0 only when N is even (it includes -pi and 0).
        return -np.pi + 2.0 * np.pi * np.arange(N) / N

    def _population_sum(self, theta, params: ModelParams, shift: float):
        f = get_nonlinearity(params.nonlinearity, params.nonlinearity_params)
        theta_arr = np.atleast_1d(np.asarray(theta, dtype=float))
        phi = self._phi_grid(params)
        N = phi.size
        # Riemann-sum approximation of int_{-pi}^{pi} f(h) dphi: weight 2 pi / N.
        phi_b = phi[None, :]
        theta_b = theta_arr[:, None]
        h = params.S * np.exp(params.kappa_h * np.cos(theta_b - phi_b + shift)) + params.A * (
            np.exp(params.kappa_g * np.cos(params.delta - phi_b))
            + np.exp(params.kappa_g * np.cos(-params.delta - phi_b))
        )
        fh = f(h)
        result = (2.0 * np.pi / N) * fh.sum(axis=1)
        if np.isscalar(theta) or (np.ndim(theta) == 0):
            return float(result[0])
        return result
