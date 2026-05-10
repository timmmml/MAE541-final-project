"""Cubic (Duffing) normal-form steering model: U = c1*theta + c3*theta^3."""

from __future__ import annotations

from dataclasses import dataclass
from math import cosh, pi, sqrt

import numpy as np

from steering.models.base import SteeringModel
from steering.params import ModelParams


@dataclass
class DuffingModel(SteeringModel):
    """Cubic truncation ``U(theta) = c1 * theta + c3 * theta**3``.

    Either supply ``c1, c3`` directly or build from ``ModelParams`` via
    :meth:`from_params`. The Bessel-derived coefficients are computed by
    delegating to :class:`BesselSteeringModel`.
    """

    c1: float
    c3: float

    @classmethod
    def from_params(cls, params: ModelParams) -> "DuffingModel":
        # Local import to avoid circular dependency.
        from steering.models.bessel import BesselSteeringModel

        c1, c3 = BesselSteeringModel().taylor_coefficients(params)
        return cls(c1=float(c1), c3=float(c3))

    # ---- SteeringModel interface ----------------------------------------

    def steering_drive(self, theta, params: ModelParams | None = None):
        theta = np.asarray(theta, dtype=float)
        return self.c1 * theta + self.c3 * theta * theta * theta

    def steering_drive_derivatives(
        self, theta, params: ModelParams | None = None, order: int = 3, h: float = 0.0
    ):
        theta = np.asarray(theta, dtype=float)
        U0 = self.c1 * theta + self.c3 * theta**3
        out = [U0]
        if order >= 1:
            out.append(self.c1 + 3.0 * self.c3 * theta * theta)
        if order >= 2:
            out.append(6.0 * self.c3 * theta)
        if order >= 3:
            out.append(np.full_like(theta, 6.0 * self.c3))
        return tuple(out[: order + 1])

    def steering_potential(self, theta, params: ModelParams | None = None):
        theta = np.asarray(theta, dtype=float)
        return -0.5 * self.c1 * theta * theta - 0.25 * self.c3 * theta**4

    # ---- Duffing-specific tools -----------------------------------------

    def fixed_points(self) -> list[float]:
        """Return the real roots of ``U(theta) = 0``."""
        if self.c3 == 0.0:
            return [0.0] if self.c1 != 0.0 else []
        roots = [0.0]
        ratio = -self.c1 / self.c3
        if ratio > 0.0:
            r = sqrt(ratio)
            roots += [-r, r]
        return sorted(roots)

    def hamiltonian(self, theta, v) -> np.ndarray:
        theta = np.asarray(theta, dtype=float)
        v = np.asarray(v, dtype=float)
        return 0.5 * v * v - 0.5 * self.c1 * theta * theta - 0.25 * self.c3 * theta**4

    def homoclinic_orbit(
        self, t_array: np.ndarray, branch: str = "positive"
    ) -> tuple[np.ndarray, np.ndarray]:
        """Closed-form sech orbit. Requires ``c1 > 0`` and ``c3 < 0``."""
        if not (self.c1 > 0.0 and self.c3 < 0.0):
            raise ValueError(
                "Homoclinic orbit requires c1 > 0 (saddle at origin) "
                "and c3 < 0 (closing wells)."
            )
        sign = 1.0 if branch == "positive" else -1.0
        amp = sqrt(2.0 * self.c1 / abs(self.c3))
        s = np.sqrt(self.c1) * t_array
        sech = 1.0 / np.cosh(s)
        theta = sign * amp * sech
        v = -sign * amp * np.sqrt(self.c1) * sech * np.tanh(s)
        return theta, v

    def melnikov_threshold(self, gamma: float, omega: float) -> float:
        """Critical forcing amplitude ``F_crit(omega)`` for chaos.

        ``F_crit = 4 gamma c1^{3/2} / (3 pi omega sqrt(2 |c3|)) * cosh(pi omega / (2 sqrt c1))``.

        Note: ``SPEC.md`` writes ``c1`` instead of ``c1^{3/2}``; that is a typo
        (verified by direct integration of the Melnikov function on the
        analytical sech homoclinic). At ``c1 = 1`` the formulas coincide.
        """
        if not (self.c1 > 0.0 and self.c3 < 0.0):
            raise ValueError("Melnikov threshold requires c1 > 0, c3 < 0.")
        c1, c3 = self.c1, self.c3
        return (
            (4.0 * gamma * c1**1.5)
            / (3.0 * pi * omega * sqrt(2.0 * abs(c3)))
            * cosh(pi * omega / (2.0 * sqrt(c1)))
        )
