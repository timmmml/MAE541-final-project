"""Thin wrapper around ``scipy.integrate.solve_ivp`` for steering dynamics."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Literal

import numpy as np
from scipy.integrate import solve_ivp

from steering.dynamics.base import Dynamics
from steering.params import ForcingParams, ModelParams


@dataclass
class SimulationResult:
    t: np.ndarray
    states: np.ndarray  # shape (n_times, state_dim)
    dynamics: Dynamics
    params: ModelParams
    forcing: ForcingParams | None = None
    sol: object | None = None  # the underlying solve_ivp OdeResult, if needed
    success: bool = True
    message: str = ""

    @property
    def theta(self) -> np.ndarray:
        return self.states[:, 0]

    @property
    def v(self) -> np.ndarray | None:
        return self.states[:, 1] if self.states.shape[1] >= 2 else None


@dataclass
class Simulation:
    """Bundle a ``Dynamics`` + parameters into a runnable integrator.

    The forcing is part of the problem specification, so it lives on the
    ``Simulation`` (not the ``Dynamics``). This means a single dynamics object
    can be reused for many force sweeps.
    """

    dynamics: Dynamics
    params: ModelParams
    forcing: ForcingParams | None = None
    method: str = "RK45"
    rtol: float = 1e-9
    atol: float = 1e-11
    max_step: float | None = None

    def _rhs(self):
        dyn = self.dynamics
        params = self.params
        forcing = self.forcing

        def f(t, y):
            return dyn.rhs(t, y, params, forcing)

        return f

    def run(
        self,
        state0: np.ndarray,
        t_span: tuple[float, float],
        t_eval: np.ndarray | None = None,
        events: list[Callable] | None = None,
        dense_output: bool = False,
    ) -> SimulationResult:
        rhs = self._rhs()
        kwargs = dict(
            method=self.method,
            rtol=self.rtol,
            atol=self.atol,
            t_eval=t_eval,
            events=events,
            dense_output=dense_output,
        )
        if self.max_step is not None:
            kwargs["max_step"] = self.max_step
        sol = solve_ivp(rhs, t_span, np.asarray(state0, dtype=float), **kwargs)
        states = sol.y.T
        # Apply topology wrapping after integration.
        states = self.dynamics.wrap_state(states)
        return SimulationResult(
            t=sol.t,
            states=states,
            dynamics=self.dynamics,
            params=self.params,
            forcing=self.forcing,
            sol=sol,
            success=sol.success,
            message=sol.message,
        )
