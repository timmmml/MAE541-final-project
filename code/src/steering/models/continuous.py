"""Continuous-population PFL model with arbitrary nonlinearity ``f``."""

from __future__ import annotations

import numpy as np
from scipy.integrate import fixed_quad

from steering.models.base import SteeringModel
from steering.nonlinearities import get_nonlinearity
from steering.params import ModelParams


def _heading_drive(phi, theta, kappa_h, shift):
    """``g_h(phi - shift, theta) = exp(kappa_h cos(theta - phi + shift))``."""
    return np.exp(kappa_h * np.cos(theta - phi + shift))


def _goal_drive(phi, delta, kappa_g):
    """Two-bump goal profile at ``+/- delta``."""
    return np.exp(kappa_g * np.cos(delta - phi)) + np.exp(kappa_g * np.cos(-delta - phi))


def _h_population(phi, theta, params: ModelParams, shift: float) -> np.ndarray:
    """Total drive into a PFL cell with goal preference ``phi``.

    ``shift = -Delta_pop`` for the L subpopulation (heading neuron preference is
    ``phi - Delta_pop``, so the heading cosine becomes ``theta - phi + Delta``)
    and ``shift = +Delta_pop`` for the R subpopulation.
    """
    return params.S * _heading_drive(phi, theta, params.kappa_h, shift) + params.A * _goal_drive(
        phi, params.delta, params.kappa_g
    )


class ContinuousPFLModel(SteeringModel):
    r"""Numerical evaluation of population sums via Gauss-Legendre quadrature.

    ``G_{P3L/R}(theta, delta) = int_{-pi}^{pi} f(h_{P3L/R}(theta, phi, delta)) dphi``,
    ``U(theta) = G_{P3L} - G_{P3R}``.

    With ``f = quadratic`` and equal ``kappa_h, kappa_g`` this matches
    :class:`BesselSteeringModel` exactly (consistency test).
    """

    def __init__(self, n_quad: int = 256):
        # Fixed Gauss-Legendre order; 256 nodes is overkill for typical kappa
        # but keeps the integral exact even for sharp bumps (kappa ~ 30).
        self.n_quad = int(n_quad)
        self._nodes, self._weights = np.polynomial.legendre.leggauss(self.n_quad)
        # Map nodes from [-1, 1] to [-pi, pi]
        self._phi_nodes = np.pi * self._nodes
        self._phi_weights = np.pi * self._weights

    # ------------------------------------------------------------------
    def _population_sum(self, theta, params: ModelParams, shift: float) -> np.ndarray:
        f = get_nonlinearity(params.nonlinearity, params.nonlinearity_params)
        theta_arr = np.atleast_1d(np.asarray(theta, dtype=float))
        # Build (n_theta, n_phi) array of h, evaluate f, integrate.
        # phi shape (n_phi,), theta shape (n_theta,)
        phi = self._phi_nodes[None, :]
        theta_b = theta_arr[:, None]
        h = params.S * np.exp(params.kappa_h * np.cos(theta_b - phi + shift)) + params.A * (
            np.exp(params.kappa_g * np.cos(params.delta - phi))
            + np.exp(params.kappa_g * np.cos(-params.delta - phi))
        )
        fh = f(h)
        result = np.einsum("ij,j->i", fh, self._phi_weights)
        if np.isscalar(theta) or (np.ndim(theta) == 0):
            return float(result[0])
        return result

    def steering_drive(self, theta, params: ModelParams):
        # PFL3L heading shift: phi - Delta_pop -> heading cos(theta - phi + Delta).
        G_L = self._population_sum(theta, params, shift=+params.Delta_pop)
        G_R = self._population_sum(theta, params, shift=-params.Delta_pop)
        return G_L - G_R

    # PFL2 has the same goal preference but heading shift of pi (180 deg) by
    # default. Expose it for the full circuit model.
    def pfl2_population_sum(self, theta, params: ModelParams) -> np.ndarray:
        return self._population_sum(theta, params, shift=params.Delta_pop_pfl2)

    def pfl3_population_sums(self, theta, params: ModelParams):
        """Return ``(G_L, G_R)`` so the full circuit model can reuse them."""
        G_L = self._population_sum(theta, params, shift=+params.Delta_pop)
        G_R = self._population_sum(theta, params, shift=-params.Delta_pop)
        return G_L, G_R
