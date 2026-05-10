"""Abstract base class for all steering models."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from scipy.integrate import quad

from steering.params import ModelParams


class SteeringModel(ABC):
    """Common interface for steering drives ``U(theta, params)``."""

    @abstractmethod
    def steering_drive(self, theta, params: ModelParams):
        """Return ``U(theta)`` for scalar or array ``theta``."""
        ...

    def steering_drive_derivatives(
        self, theta, params: ModelParams, order: int = 3, h: float = 1e-3
    ):
        """Return ``(U, U', U'', U''')`` at ``theta`` via central differences.

        Subclasses with closed-form derivatives may override.
        """
        theta = np.asarray(theta, dtype=float)

        def U(t):
            return self.steering_drive(t, params)

        Uv = U(theta)
        out = [Uv]
        if order >= 1:
            out.append((U(theta + h) - U(theta - h)) / (2.0 * h))
        if order >= 2:
            out.append((U(theta + h) - 2.0 * U(theta) + U(theta - h)) / (h * h))
        if order >= 3:
            out.append(
                (U(theta + 2 * h) - 2 * U(theta + h) + 2 * U(theta - h) - U(theta - 2 * h))
                / (2.0 * h**3)
            )
        return tuple(out[: order + 1])

    def steering_potential(self, theta, params: ModelParams) -> np.ndarray:
        """Return ``V(theta) = -int_0^theta U(t) dt`` by quadrature."""
        theta = np.atleast_1d(np.asarray(theta, dtype=float))
        out = np.empty_like(theta)
        for i, th in enumerate(theta):
            val, _ = quad(
                lambda t: float(self.steering_drive(t, params)),
                0.0,
                float(th),
                limit=200,
            )
            out[i] = -val
        return out if theta.ndim else float(out[0])
