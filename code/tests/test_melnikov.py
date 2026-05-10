"""Melnikov tests: analytical vs numerical, formula consistency."""

from __future__ import annotations

import numpy as np
from scipy.integrate import simpson

from steering.analysis.melnikov import (
    melnikov_analytical,
    melnikov_critical_forcing,
    melnikov_critical_forcing_numerical,
    melnikov_numerical,
)
from steering.models import BesselSteeringModel, DuffingModel
from steering.params import ModelParams


def test_analytical_formula_matches_direct_integration():
    """F_crit formula matches direct integration along the analytical sech orbit."""
    c1, c3, gamma, omega = 4.0, -1.0, 0.1, 1.0
    duff = DuffingModel(c1=c1, c3=c3)
    t = np.linspace(-30.0, 30.0, 200_001)
    th, v = duff.homoclinic_orbit(t)
    I1 = simpson(v * v, x=t)
    a = simpson(v * np.cos(omega * t), x=t)
    b = simpson(v * np.sin(omega * t), x=t)
    R = np.hypot(a, b)
    F_direct = gamma * I1 / R
    F_formula = melnikov_critical_forcing(c1, c3, gamma, omega)
    assert abs(F_direct - F_formula) / F_formula < 1e-3


def test_analytical_M_matches_direct():
    """Analytical M(t_0) matches integration of the same quantity."""
    c1, c3, gamma, omega, F = 4.0, -1.0, 0.1, 1.0, 0.5
    duff = DuffingModel(c1=c1, c3=c3)
    t = np.linspace(-30.0, 30.0, 200_001)
    _, v = duff.homoclinic_orbit(t)
    t0_arr = np.linspace(0.0, 2 * np.pi / omega, 17)
    M_formula = melnikov_analytical(c1, c3, gamma, omega, F, t0_arr)
    # Direct: M(t0) = -gamma int v^2 dt + F int v cos(omega(t+t0)) dt.
    I1 = simpson(v * v, x=t)
    a = simpson(v * np.cos(omega * t), x=t)
    b = simpson(v * np.sin(omega * t), x=t)
    M_direct = -gamma * I1 + F * (np.cos(omega * t0_arr) * a - np.sin(omega * t0_arr) * b)
    assert np.max(np.abs(M_formula - M_direct)) < 1e-3


def test_bessel_numerical_close_to_duffing_near_bifurcation():
    """Near the pitchfork bifurcation Bessel F_crit ~ Duffing F_crit."""
    p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.34)
    bessel = BesselSteeringModel()
    c1, c3 = bessel.taylor_coefficients(p)
    F_a = melnikov_critical_forcing(c1, c3, 0.05, 1.0)
    F_n = melnikov_critical_forcing_numerical(bessel, p, 0.05, 1.0)
    # theta_max(Duffing) ~ 0.53 here; Bessel correction lifts F_crit by O(20%).
    assert 0.7 < F_n / F_a < 2.0
