"""Abstract base class for dynamics."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Literal

import numpy as np

from steering.params import ForcingParams, ModelParams


class Dynamics(ABC):
    """Wraps a :class:`SteeringModel` to produce an ODE right-hand side.

    Subclasses choose the control law (velocity vs acceleration) and the state
    space topology (planar ``R^n`` or cylindrical ``S^1 x R^(n-1)``).
    """

    topology: Literal["planar", "cylindrical"] = "planar"

    @abstractmethod
    def rhs(
        self,
        t: float,
        state: np.ndarray,
        params: ModelParams,
        forcing: ForcingParams | None = None,
    ) -> np.ndarray: ...

    @abstractmethod
    def state_dim(self) -> int: ...

    def wrap_state(self, state: np.ndarray) -> np.ndarray:
        """Apply topology to a state vector. Default: identity (planar)."""
        if self.topology == "cylindrical":
            out = np.array(state, dtype=float, copy=True)
            out[..., 0] = (out[..., 0] + np.pi) % (2.0 * np.pi) - np.pi
            return out
        return np.asarray(state, dtype=float)
