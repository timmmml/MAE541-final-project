"""Stroboscopic Poincaré sections and bifurcation diagrams."""

from __future__ import annotations

from math import pi

import numpy as np

from steering.dynamics.acceleration import AccelerationDynamics
from steering.integrator import Simulation, SimulationResult
from steering.models.base import SteeringModel
from steering.params import ForcingParams, ModelParams


def _normalize_initial_states(state0: np.ndarray | None) -> np.ndarray:
    """Return initial conditions as an array of shape ``(n_states, state_dim)``."""
    if state0 is None:
        return np.array([[0.1, 0.0]], dtype=float)

    states0 = np.asarray(state0, dtype=float)
    if states0.ndim == 1:
        return states0[None, :]
    if states0.ndim == 2:
        return states0.copy()
    raise ValueError("state0 must be a 1D state vector or a 2D array of states")


def stroboscopic_section(
    result: SimulationResult,
    omega: float,
    transient_periods: int = 100,
) -> np.ndarray:
    """Sample state at ``t = 2 pi n / omega`` (after a transient).

    Uses ``solve_ivp``'s dense output if available, else linear interpolation.
    """
    T = 2.0 * pi / omega
    t_start = result.t[0] + transient_periods * T
    t_end = result.t[-1]
    n_strobe = max(int((t_end - t_start) / T), 0)
    strobe_times = t_start + np.arange(n_strobe) * T
    if result.sol is not None and getattr(result.sol, "sol", None) is not None:
        states = result.sol.sol(strobe_times).T
    else:
        # Linear interpolation in each component.
        n_dim = result.states.shape[1]
        states = np.column_stack([
            np.interp(strobe_times, result.t, result.states[:, k]) for k in range(n_dim)
        ])
    if result.dynamics.topology == "cylindrical":
        states = result.dynamics.wrap_state(states)
    return states


def bifurcation_diagram_poincare(
    dynamics: AccelerationDynamics,
    model: SteeringModel,
    params: ModelParams,
    gamma: float,
    omega: float,
    sweep_param: str,
    sweep_values: np.ndarray,
    F: float | None = None,
    state0: np.ndarray | None = None,
    n_periods_transient: int = 200,
    n_periods_record: int = 100,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    continuation: bool = True,
) -> dict:
    """Sweep ``sweep_param`` over ``sweep_values`` and return strobe ``theta`` arrays.

    For each sweep value we run ``(transient + record)`` forcing periods and
    record ``theta`` at stroboscopic times. ``sweep_param`` can be any field of
    ``ModelParams`` or ``F``, ``omega``, ``forcing_type`` (the latter on the
    forcing). ``state0`` may be a single initial condition of shape
    ``(state_dim,)`` or multiple initial conditions of shape
    ``(n_initial_conditions, state_dim)``. With ``continuation=True`` each run
    starts from the previous run's final state for each initial condition,
    smoothing the bifurcation tracking.
    """
    initial_states = _normalize_initial_states(state0)
    dyn_gamma = AccelerationDynamics(
        model=dynamics.model, gamma=gamma, topology=dynamics.topology
    )
    T = 2.0 * pi / omega
    t_total = (n_periods_transient + n_periods_record) * T
    forcing_keys = {"F", "omega", "forcing_type"}

    out_thetas: list[np.ndarray] = []
    out_vs: list[np.ndarray] = []
    cur_states = initial_states.copy()
    for v in np.asarray(sweep_values):
        if sweep_param in forcing_keys:
            f_kwargs = {"omega": omega}
            if F is not None:
                f_kwargs["F"] = F
            f_kwargs[sweep_param] = float(v) if sweep_param != "forcing_type" else str(v)
            forcing = ForcingParams(**f_kwargs)
            cur_params = params
        else:
            cur_params = params.replace(**{sweep_param: float(v)})
            forcing = ForcingParams(F=F if F is not None else 0.0, omega=omega)
        sim = Simulation(
            dyn_gamma, cur_params, forcing, rtol=rtol, atol=atol
        )

        thetas_for_value: list[np.ndarray] = []
        vs_for_value: list[np.ndarray] = []
        next_states: list[np.ndarray] = []
        for cur_state in cur_states:
            result = sim.run(cur_state, (0.0, t_total), dense_output=True)
            section = stroboscopic_section(result, omega, transient_periods=n_periods_transient)
            if section.size > 0:
                thetas_for_value.append(section[:, 0])
                vs_for_value.append(section[:, 1])
                next_states.append(section[-1].copy())
            else:
                next_states.append(np.array(cur_state, copy=True))

        out_thetas.append(
            np.concatenate(thetas_for_value) if thetas_for_value else np.array([], dtype=float)
        )
        out_vs.append(
            np.concatenate(vs_for_value) if vs_for_value else np.array([], dtype=float)
        )
        if continuation:
            cur_states = np.asarray(next_states, dtype=float)
    return {
        "sweep_param": sweep_param,
        "sweep_values": np.asarray(sweep_values),
        "thetas": out_thetas,
        "vs": out_vs,
    }
