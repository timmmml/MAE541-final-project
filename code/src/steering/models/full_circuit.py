"""Full circuit with PFL2 indirect-pathway modulation."""

from __future__ import annotations

import numpy as np

from steering.models.continuous import ContinuousPFLModel
from steering.models.discrete import DiscretePFLModel
from steering.nonlinearities import get_nonlinearity
from steering.params import ModelParams


class FullCircuitModel:
    """PFL3 direct + PFL2 indirect (DNa03) pathway.

    ``U(theta) = (G_{P3L} - G_{P3R}) + W_{D3} [f_d(G_{P3L} + G_{P2}) - f_d(G_{P3R} + G_{P2})]``.

    Reduces to the ``ContinuousPFLModel`` (or ``DiscretePFLModel`` when
    ``N_neurons`` is set) when ``W_{D3} = 0``.
    """

    def __init__(self, n_quad: int = 256):
        # Two underlying engines: continuous (integral) or discrete (sum).
        self._continuous = ContinuousPFLModel(n_quad=n_quad)
        self._discrete = DiscretePFLModel()

    def _engine(self, params: ModelParams) -> ContinuousPFLModel:
        return self._discrete if params.N_neurons else self._continuous

    def steering_drive(self, theta, params: ModelParams):
        engine = self._engine(params)
        G_L, G_R = engine.pfl3_population_sums(theta, params)
        if params.W_D3 == 0.0:
            return G_L - G_R
        # Use a (possibly distinct) DNa03 nonlinearity for the indirect path.
        f_d_spec = params.nonlinearity_pfl2 or params.nonlinearity
        f_d = get_nonlinearity(f_d_spec, params.nonlinearity_params)
        G_P2 = engine.pfl2_population_sum(theta, params)
        indirect = f_d(G_L + G_P2) - f_d(G_R + G_P2)
        return (G_L - G_R) + params.W_D3 * indirect

    def steering_drive_derivatives(self, theta, params: ModelParams, order: int = 3, h: float = 1e-3):
        # Inherit finite-difference fallback from base SteeringModel.
        from steering.models.base import SteeringModel

        return SteeringModel.steering_drive_derivatives(self, theta, params, order=order, h=h)

    def steering_potential(self, theta, params: ModelParams):
        from steering.models.base import SteeringModel

        return SteeringModel.steering_potential(self, theta, params)
