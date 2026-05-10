"""Homoclinic-orbit construction (analytical and numerical)."""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from steering.analysis.fixed_points import FixedPoint2D, find_fixed_points_2d
from steering.dynamics.acceleration import AccelerationDynamics
from steering.models.base import SteeringModel
from steering.models.duffing import DuffingModel
from steering.params import ModelParams


def duffing_homoclinic(
    c1: float, c3: float, t_array: np.ndarray, branch: str = "positive"
) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form sech homoclinic orbit of the cubic Hamiltonian.

    Requires ``c1 > 0`` and ``c3 < 0``.
    """
    return DuffingModel(c1=c1, c3=c3).homoclinic_orbit(t_array, branch=branch)


def numerical_homoclinic(
    model: SteeringModel,
    params: ModelParams,
    saddle: FixedPoint2D | None = None,
    t_max: float = 50.0,
    eps: float = 1e-8,
    branch: str = "positive",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Numerically integrate the homoclinic orbit of the undamped system.

    Procedure: find the saddle at ``theta = 0``; perturb along the unstable
    eigenvector; integrate forward until the orbit reaches the apex (``v = 0``)
    and reflect by time-reversal symmetry to get the second half. This gives
    an exponentially-accurate homoclinic without depending on the orbit
    "returning" close to the saddle (which is asymptotic).
    """
    dyn = AccelerationDynamics(model=model, gamma=0.0, topology="planar")
    if saddle is None:
        fps = find_fixed_points_2d(dyn, params)
        saddle_candidates = [fp for fp in fps if "saddle" in fp.classification]
        if not saddle_candidates:
            raise RuntimeError("No saddle found; cannot construct homoclinic orbit.")
        saddle = min(saddle_candidates, key=lambda fp: abs(fp.theta))

    J = dyn.jacobian(saddle.theta, 0.0, params)
    eigvals, eigvecs = np.linalg.eig(J)
    iu = int(np.argmax(eigvals.real))
    eu = np.real(eigvecs[:, iu])
    eu = eu / np.linalg.norm(eu)
    sign = 1.0 if branch == "positive" else -1.0
    # Perturb in the direction that makes v positive (so we integrate the
    # "going up" half).
    if eu[1] * sign < 0:
        eu = -eu
    base = np.array([saddle.theta, 0.0])
    init = base + sign * eps * eu

    def rhs(t, y):
        return dyn.rhs(t, y, params)

    # Apex: v = 0 going from positive to negative (or negative to positive on
    # the negative branch). Use direction depending on branch sign.
    def apex_event(t, y):
        return y[1]

    apex_event.terminal = True
    apex_event.direction = -1.0 * sign  # we leave with v > 0 on positive branch.

    sol_f = solve_ivp(
        rhs,
        (0.0, t_max),
        init,
        rtol=1e-11,
        atol=1e-13,
        events=apex_event,
        dense_output=False,
        max_step=0.05,
    )
    if sol_f.t_events[0].size == 0:
        # Apex not found — orbit didn't close. Could be heteroclinic or
        # the parameter set is outside the homoclinic regime.
        raise RuntimeError(
            "Numerical homoclinic: apex (v=0) not reached within t_max. "
            "Orbit may be heteroclinic (saddle to saddle) or unbounded."
        )

    t_half = sol_f.t
    y_half = sol_f.y
    # Mirror by time-reversal: the second half is the first half reversed,
    # with v negated.
    t_full = np.concatenate([-t_half[::-1][:-1], t_half])
    theta_full = np.concatenate([y_half[0, ::-1][:-1], y_half[0]])
    v_full = np.concatenate([-y_half[1, ::-1][:-1], y_half[1]])
    return t_full, theta_full, v_full
