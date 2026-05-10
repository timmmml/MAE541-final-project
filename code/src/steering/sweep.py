"""Parameter sweep infrastructure (multiprocessing-friendly)."""

from __future__ import annotations

import itertools
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from steering.dynamics.base import Dynamics
from steering.integrator import Simulation, SimulationResult
from steering.params import ForcingParams, ModelParams


def _run_one(args):
    dynamics, base_params, base_forcing, axes, idx, state0, t_span, sim_kwargs, analysis_fn = args
    overrides = {axis: vals[i] for axis, (vals, i) in zip(axes.keys(), zip(axes.values(), idx))}
    # Some overrides may be on ForcingParams (F, omega, forcing_type).
    forcing_keys = {"F", "omega", "forcing_type"}
    p_over = {k: v for k, v in overrides.items() if k not in forcing_keys}
    f_over = {k: v for k, v in overrides.items() if k in forcing_keys}
    params = base_params.replace(**p_over) if p_over else base_params
    forcing = base_forcing
    if f_over:
        if forcing is None:
            forcing = ForcingParams()
        from dataclasses import replace as _replace

        forcing = _replace(forcing, **f_over)
    sim = Simulation(dynamics, params, forcing, **sim_kwargs)
    result = sim.run(state0, t_span)
    return idx, analysis_fn(result)


@dataclass
class ParameterSweep:
    """Run ``Simulation`` over a Cartesian product of parameter values.

    Parameters
    ----------
    dynamics : Dynamics
    base_params : ModelParams
        Defaults; overridden by ``sweep_axes``.
    sweep_axes : dict[str, np.ndarray]
        Mapping ``parameter name -> array of values``. Recognized names are any
        field of ``ModelParams`` plus ``F, omega, forcing_type`` (which set
        ``ForcingParams``).
    analysis_fn : Callable[[SimulationResult], Any]
        Reduces a simulation to a small piece of data (e.g. final state,
        fixed-point classification, max Lyapunov exponent).
    base_forcing : ForcingParams | None
    sim_kwargs : dict
        Forwarded to :class:`Simulation` (e.g. ``method``, ``rtol``).
    """

    dynamics: Dynamics
    base_params: ModelParams
    sweep_axes: dict[str, np.ndarray]
    analysis_fn: Callable[[SimulationResult], Any]
    base_forcing: ForcingParams | None = None
    sim_kwargs: dict | None = None

    def run(
        self,
        state0: np.ndarray,
        t_span: tuple[float, float],
        n_workers: int | None = None,
    ):
        axes = {k: np.asarray(v) for k, v in self.sweep_axes.items()}
        shape = tuple(len(v) for v in axes.values())
        # Pre-allocate as object array to allow analysis_fn to return anything.
        out = np.empty(shape, dtype=object)

        sim_kwargs = self.sim_kwargs or {}
        all_indices = list(itertools.product(*[range(s) for s in shape]))
        jobs = [
            (
                self.dynamics,
                self.base_params,
                self.base_forcing,
                axes,
                idx,
                state0,
                t_span,
                sim_kwargs,
                self.analysis_fn,
            )
            for idx in all_indices
        ]

        if n_workers is None or n_workers == 1:
            for j in jobs:
                idx, val = _run_one(j)
                out[idx] = val
        else:
            with ProcessPoolExecutor(max_workers=n_workers) as ex:
                for idx, val in ex.map(_run_one, jobs):
                    out[idx] = val
        return SweepResult(values=out, axes=axes)


@dataclass
class SweepResult:
    values: np.ndarray  # object array
    axes: dict[str, np.ndarray]

    def to_xarray(self):
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError("xarray is required for SweepResult.to_xarray") from e

        # If all values are scalar-like, stack into a numeric array.
        try:
            data = np.asarray(self.values.tolist(), dtype=float)
        except (TypeError, ValueError):
            data = self.values
        return xr.DataArray(
            data,
            coords={k: v for k, v in self.axes.items()},
            dims=list(self.axes.keys()),
        )
