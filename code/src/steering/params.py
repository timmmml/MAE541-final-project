"""Parameter dataclasses for steering models and forcing."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from math import pi
from typing import Any, Callable

import numpy as np


@dataclass(frozen=True)
class ModelParams:
    """Parameters shared by every steering model.

    All angles are in radians. The PFL3 population phase shift ``Delta_pop``
    defaults to 67.5 deg = 3 pi / 8, which is the value used throughout the
    spec; it can be overridden when exploring perturbations of the population
    geometry.
    """

    kappa_h: float
    kappa_g: float
    delta: float
    Delta_pop: float = 3 * pi / 8
    S: float = 1.0
    A: float = 1.0
    W_D3: float = 0.0
    Delta_pop_pfl2: float = pi
    N_neurons: int | None = None
    nonlinearity: str | Callable[[np.ndarray], np.ndarray] = "quadratic"
    nonlinearity_params: dict[str, Any] = field(default_factory=dict)
    nonlinearity_pfl2: str | Callable | None = None

    @property
    def kappa_equal(self) -> float | None:
        """Return ``kappa`` if ``kappa_h == kappa_g``, else ``None``."""
        if np.isclose(self.kappa_h, self.kappa_g):
            return float(self.kappa_h)
        return None

    def replace(self, **changes: Any) -> "ModelParams":
        """Return a new ``ModelParams`` with the listed fields replaced."""
        return replace(self, **changes)


@dataclass(frozen=True)
class ForcingParams:
    """External forcing applied to ``AccelerationDynamics``.

    ``forcing_type='additive'`` adds ``F cos(omega t)`` directly to ``v_dot``.
    ``forcing_type='parametric_delta'`` modulates the goal half-separation as
    ``delta(t) = delta_0 + F cos(omega t)`` so that ``U`` itself is time
    dependent.
    """

    F: float = 0.0
    omega: float = 1.0
    forcing_type: str = "additive"

    def __post_init__(self) -> None:
        if self.forcing_type not in ("additive", "parametric_delta"):
            raise ValueError(
                f"Unknown forcing_type {self.forcing_type!r}; "
                "use 'additive' or 'parametric_delta'."
            )
