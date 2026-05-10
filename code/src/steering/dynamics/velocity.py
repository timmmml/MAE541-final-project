"""1D velocity-control dynamics: theta_dot = U(theta)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from steering.dynamics.base import Dynamics
from steering.models.base import SteeringModel
from steering.params import ForcingParams, ModelParams


@dataclass
class VelocityDynamics(Dynamics):
    """``theta_dot = U(theta, delta(t))``.

    Pure additive forcing is degenerate in 1D (the system stays bounded,
    cannot be chaotic), but parametric forcing through ``delta(t)`` produces a
    non-autonomous 1.5D system that can be interesting.
    """

    model: SteeringModel
    topology: Literal["planar", "cylindrical"] = "cylindrical"

    def state_dim(self) -> int:
        return 1

    def rhs(
        self,
        t: float,
        state: np.ndarray,
        params: ModelParams,
        forcing: ForcingParams | None = None,
    ) -> np.ndarray:
        theta = float(state[0])
        if forcing is not None and forcing.forcing_type == "parametric_delta":
            delta_t = params.delta + forcing.F * np.cos(forcing.omega * t)
            params = params.replace(delta=delta_t)
        U = float(np.asarray(self.model.steering_drive(theta, params)))
        if forcing is not None and forcing.forcing_type == "additive":
            U += forcing.F * np.cos(forcing.omega * t)
        return np.array([U])
