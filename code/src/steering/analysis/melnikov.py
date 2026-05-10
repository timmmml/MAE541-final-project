"""Melnikov analysis for additive and parametric forcing."""

from __future__ import annotations

from math import cosh, pi, sqrt

import numpy as np
from scipy.integrate import simpson

from steering.analysis.homoclinic import duffing_homoclinic, numerical_homoclinic
from steering.models.base import SteeringModel
from steering.params import ModelParams


def melnikov_analytical(
    c1: float,
    c3: float,
    gamma: float,
    omega: float,
    F: float,
    t0_array: np.ndarray,
) -> np.ndarray:
    """Closed-form ``M(t_0)`` for the Duffing reduction with additive forcing.

    ``M(t_0) = -4 gamma c_1^{3/2}/(3 |c_3|) + F pi omega sqrt(2/|c_3|) / cosh(pi omega / (2 sqrt c_1)) * sin(omega t_0)``.

    Note: SPEC.md writes the forcing amplitude as ``sqrt(2 c_1/|c_3|)`` and
    ``F_crit`` with a bare ``c_1`` factor; both contain the same typo (a
    spurious ``sqrt(c_1)`` in the forcing-amplitude expression). The forms
    here are verified by direct integration along the analytical sech orbit.
    """
    if not (c1 > 0.0 and c3 < 0.0):
        raise ValueError("Analytical Melnikov requires c1 > 0, c3 < 0.")
    damping_part = -4.0 * gamma * c1 ** 1.5 / (3.0 * abs(c3))
    forcing_amp = F * pi * omega * sqrt(2.0 / abs(c3)) / cosh(pi * omega / (2.0 * sqrt(c1)))
    return damping_part + forcing_amp * np.sin(omega * np.asarray(t0_array))


def melnikov_critical_forcing(
    c1: float, c3: float, gamma: float, omega: float
) -> float:
    """``F_crit(omega) = 4 gamma c_1^{3/2} / (3 pi omega sqrt(2 |c_3|)) * cosh(pi omega / (2 sqrt c_1))``.

    See note on :func:`melnikov_analytical` for the discrepancy with SPEC.md.
    """
    if not (c1 > 0.0 and c3 < 0.0):
        raise ValueError("Melnikov threshold requires c1 > 0, c3 < 0.")
    return (
        (4.0 * gamma * c1 ** 1.5)
        / (3.0 * pi * omega * sqrt(2.0 * abs(c3)))
        * cosh(pi * omega / (2.0 * sqrt(c1)))
    )


def melnikov_numerical(
    model: SteeringModel,
    params: ModelParams,
    gamma: float,
    omega: float,
    F: float,
    t0_array: np.ndarray,
    homoclinic: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> np.ndarray:
    """Compute ``M(t_0)`` along the numerical homoclinic of ``model``.

    ``M(t_0) = -gamma int v_0(t)^2 dt + F int v_0(t) cos(omega(t + t_0)) dt``.
    """
    if homoclinic is None:
        homoclinic = numerical_homoclinic(model, params)
    t, theta_h, v_h = homoclinic
    damping_int = simpson(v_h * v_h, x=t)
    cos_term = simpson(v_h * np.cos(omega * t), x=t)
    sin_term = simpson(v_h * np.sin(omega * t), x=t)
    t0 = np.asarray(t0_array)
    # F * int v_0 cos(omega(t + t_0)) dt = F [cos(omega t_0) cos_term - sin(omega t_0) sin_term]
    forcing = F * (np.cos(omega * t0) * cos_term - np.sin(omega * t0) * sin_term)
    return -gamma * damping_int + forcing


def melnikov_critical_forcing_numerical(
    model: SteeringModel,
    params: ModelParams,
    gamma: float,
    omega: float,
    homoclinic: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> float:
    """Numerical ``F_crit`` from the homoclinic of ``model``.

    Set ``M(t_0)`` = damping_term + F * sinusoid; find F at which the sinusoid
    amplitude equals the damping bias.
    """
    if homoclinic is None:
        homoclinic = numerical_homoclinic(model, params)
    t, _, v_h = homoclinic
    damping_int = simpson(v_h * v_h, x=t)
    cos_term = simpson(v_h * np.cos(omega * t), x=t)
    sin_term = simpson(v_h * np.sin(omega * t), x=t)
    R = np.hypot(cos_term, sin_term)
    if R <= 0:
        return np.inf
    return float(gamma * damping_int / R)
