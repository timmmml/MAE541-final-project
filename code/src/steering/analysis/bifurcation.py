"""Pitchfork bifurcation analysis in the (kappa, delta) plane (and others)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from steering.analysis.fixed_points import find_fixed_points_1d
from steering.models.base import SteeringModel
from steering.models.bessel import BesselSteeringModel
from steering.params import ModelParams


@dataclass
class BifurcationDiagram:
    param1_name: str
    param2_name: str
    param1_values: np.ndarray
    param2_values: np.ndarray
    c1_grid: np.ndarray
    c3_grid: np.ndarray
    metadata: dict = field(default_factory=dict)


def _taylor(model: SteeringModel, params: ModelParams) -> tuple[float, float]:
    """Return ``(c1, c3)``: analytical for Bessel, finite difference otherwise."""
    if isinstance(model, BesselSteeringModel):
        return model.taylor_coefficients(params)
    h = 5e-4
    U = lambda t: float(np.asarray(model.steering_drive(t, params)))
    c1 = (U(h) - U(-h)) / (2.0 * h)
    c3 = (U(2 * h) - 2.0 * U(h) + 2.0 * U(-h) - U(-2 * h)) / (12.0 * h**3)
    # Note: U(2h) - 2U(h) + 2U(-h) - U(-2h) = 2 h^3 U'''(0) + O(h^5), so
    # divide by 2 h^3 to get U'''(0), then divide by 6 for c3.
    c3 = (U(2 * h) - 2 * U(h) + 2 * U(-h) - U(-2 * h)) / (2.0 * h**3) / 6.0
    return float(c1), float(c3)


def pitchfork_bifurcation_curve(
    model: SteeringModel,
    param1_name: str,
    param1_range: np.ndarray,
    param2_name: str,
    param2_range: np.ndarray,
    base_params: ModelParams,
) -> BifurcationDiagram:
    """Compute ``c1, c3`` on a 2D grid in parameter space.

    The pitchfork locus is the contour ``c1 = 0``; ``sgn(c3)`` along it
    distinguishes supercritical (``c3 < 0``) from subcritical (``c3 > 0``).
    """
    p1 = np.asarray(param1_range)
    p2 = np.asarray(param2_range)
    c1_grid = np.empty((len(p1), len(p2)))
    c3_grid = np.empty_like(c1_grid)
    for i, v1 in enumerate(p1):
        for j, v2 in enumerate(p2):
            params = base_params.replace(**{param1_name: float(v1), param2_name: float(v2)})
            c1, c3 = _taylor(model, params)
            c1_grid[i, j] = c1
            c3_grid[i, j] = c3
    return BifurcationDiagram(
        param1_name=param1_name,
        param2_name=param2_name,
        param1_values=p1,
        param2_values=p2,
        c1_grid=c1_grid,
        c3_grid=c3_grid,
    )


def numerical_bifurcation_diagram_1d(
    model: SteeringModel,
    param_name: str,
    param_range: np.ndarray,
    base_params: ModelParams,
    theta_range: tuple = (-np.pi, np.pi),
) -> dict:
    """For each value of ``param_name``, list all fixed points (theta, stability).

    Returns ``{'param': vals, 'thetas': list[list[float]], 'stabilities': list[list[str]]}``.
    """
    vals = np.asarray(param_range)
    theta_lists, stab_lists = [], []
    for v in vals:
        p = base_params.replace(**{param_name: float(v)})
        fps = find_fixed_points_1d(model, p, theta_range)
        theta_lists.append([fp.theta for fp in fps])
        stab_lists.append([fp.stability for fp in fps])
    return {"param": vals, "thetas": theta_lists, "stabilities": stab_lists}
