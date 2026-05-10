"""2D acceleration-control dynamics: ddot theta + gamma dot theta = U + F cos(omega t)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from steering.dynamics.base import Dynamics
from steering.models.base import SteeringModel
from steering.params import ForcingParams, ModelParams


@dataclass
class AccelerationDynamics(Dynamics):
    """Damped, forced second-order steering dynamics.

    State is ``(theta, v)``. Equations of motion::

        theta_dot = v
        v_dot     = -gamma * v + U(theta, delta(t)) + F cos(omega t)

    With ``forcing_type='parametric_delta'`` the goal half-separation itself
    oscillates: ``delta(t) = delta_0 + F cos(omega t)``, and the additive
    ``F cos(omega t)`` term is dropped.
    """

    model: SteeringModel
    gamma: float = 0.0
    topology: Literal["planar", "cylindrical"] = "planar"

    def state_dim(self) -> int:
        return 2

    def rhs(
        self,
        t: float,
        state: np.ndarray,
        params: ModelParams,
        forcing: ForcingParams | None = None,
    ) -> np.ndarray:
        theta = float(state[0])
        v = float(state[1])
        if forcing is not None and forcing.forcing_type == "parametric_delta":
            delta_t = params.delta + forcing.F * np.cos(forcing.omega * t)
            params = params.replace(delta=delta_t)
            U = float(np.asarray(self.model.steering_drive(theta, params)))
        else:
            U = float(np.asarray(self.model.steering_drive(theta, params)))
            if forcing is not None and forcing.forcing_type == "additive":
                U += forcing.F * np.cos(forcing.omega * t)
        return np.array([v, -self.gamma * v + U])

    # Convenience wrappers used by analysis tools.
    def jacobian(
        self,
        theta: float,
        v: float,
        params: ModelParams,
        h: float = 1e-4,
    ) -> np.ndarray:
        """Jacobian of the autonomous (unforced) RHS at ``(theta, v)``."""
        # Use analytic U' if available.
        try:
            _, Up = self.model.steering_drive_derivatives(theta, params, order=1)
            Up = float(np.asarray(Up))
        except Exception:
            Up = (
                float(np.asarray(self.model.steering_drive(theta + h, params)))
                - float(np.asarray(self.model.steering_drive(theta - h, params)))
            ) / (2.0 * h)
        return np.array([[0.0, 1.0], [Up, -self.gamma]])
