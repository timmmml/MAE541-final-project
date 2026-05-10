"""Fixed-point detection and stability classification."""

from __future__ import annotations

from dataclasses import dataclass
from math import pi
from typing import Literal

import numpy as np
from scipy.optimize import brentq

from steering.dynamics.acceleration import AccelerationDynamics
from steering.models.base import SteeringModel
from steering.params import ModelParams


@dataclass
class FixedPoint:
    theta: float
    U_prime: float
    stability: Literal["stable", "unstable", "marginal"]


@dataclass
class FixedPoint2D:
    theta: float
    eigenvalues: tuple[complex, complex]
    classification: str


def _bracket_zeros(theta_grid: np.ndarray, vals: np.ndarray):
    """Yield ``(a, b)`` brackets where ``vals`` changes sign on consecutive nodes."""
    for i in range(len(theta_grid) - 1):
        if vals[i] == 0.0:
            yield theta_grid[i], theta_grid[i]
        elif vals[i] * vals[i + 1] < 0.0:
            yield theta_grid[i], theta_grid[i + 1]


def find_fixed_points_1d(
    model: SteeringModel,
    params: ModelParams,
    theta_range: tuple[float, float] = (-pi, pi),
    n_initial: int = 400,
    tol: float = 1e-10,
) -> list[FixedPoint]:
    """Find zeros of ``U(theta)`` on ``theta_range`` and classify by ``sgn U'``.

    For 1D autonomous flow ``theta_dot = U(theta)``, a fixed point is stable
    iff ``U' < 0`` there.
    """
    grid = np.linspace(*theta_range, n_initial)
    Uvals = np.asarray(model.steering_drive(grid, params))
    fps: list[FixedPoint] = []
    seen: list[float] = []
    for a, b in _bracket_zeros(grid, Uvals):
        if a == b:
            theta_star = float(a)
        else:
            try:
                theta_star = brentq(
                    lambda t: float(np.asarray(model.steering_drive(t, params))),
                    a,
                    b,
                    xtol=tol,
                )
            except ValueError:
                continue
        if any(abs(theta_star - s) < 1e-6 for s in seen):
            continue
        seen.append(theta_star)
        try:
            _, Up = model.steering_drive_derivatives(theta_star, params, order=1)
            Up = float(np.asarray(Up))
        except Exception:
            h = 1e-5
            Up = (
                float(np.asarray(model.steering_drive(theta_star + h, params)))
                - float(np.asarray(model.steering_drive(theta_star - h, params)))
            ) / (2.0 * h)
        if abs(Up) < 1e-10:
            stability = "marginal"
        elif Up < 0:
            stability = "stable"
        else:
            stability = "unstable"
        fps.append(FixedPoint(theta=theta_star, U_prime=Up, stability=stability))
    return sorted(fps, key=lambda fp: fp.theta)


def _classify_2d(eigs: tuple[complex, complex]) -> str:
    e1, e2 = eigs
    real1, real2 = e1.real, e2.real
    imag1, imag2 = abs(e1.imag), abs(e2.imag)
    is_complex = (imag1 > 1e-10) or (imag2 > 1e-10)
    # In a 2x2 real matrix the two eigenvalues are either both real or
    # complex conjugates, so ``imag1 == imag2``.
    if is_complex:
        if real1 > 1e-10:
            return "unstable_spiral"
        if real1 < -1e-10:
            return "stable_spiral"
        return "center"
    # Both real.
    if real1 * real2 < 0.0:
        return "saddle"
    if real1 < 0.0 and real2 < 0.0:
        return "stable_node"
    if real1 > 0.0 and real2 > 0.0:
        return "unstable_node"
    return "degenerate"


def find_fixed_points_2d(
    dynamics: AccelerationDynamics,
    params: ModelParams,
    theta_range: tuple[float, float] = (-pi, pi),
    n_initial: int = 400,
) -> list[FixedPoint2D]:
    """Fixed points of the 2D system ``(theta_dot=v, v_dot=-gamma v + U(theta))``.

    These coincide with the zeros of ``U(theta)`` (with ``v = 0``); only the
    stability classification differs from the 1D case.
    """
    fps_1d = find_fixed_points_1d(dynamics.model, params, theta_range, n_initial)
    out: list[FixedPoint2D] = []
    for fp in fps_1d:
        J = dynamics.jacobian(fp.theta, 0.0, params)
        eigs = np.linalg.eigvals(J)
        eig_tuple = (complex(eigs[0]), complex(eigs[1]))
        out.append(
            FixedPoint2D(
                theta=fp.theta,
                eigenvalues=eig_tuple,
                classification=_classify_2d(eig_tuple),
            )
        )
    return out
