"""Detailed checks on the analytical Taylor coefficients."""

from __future__ import annotations

import numpy as np
import pytest

from steering.models import BesselSteeringModel
from steering.params import ModelParams


@pytest.mark.parametrize("delta", [0.2, 0.5, 1.0, 1.3, 1.5])
@pytest.mark.parametrize("kappa", [1.0, 2.0, 4.0])
def test_c1_matches_finite_difference(delta, kappa):
    bessel = BesselSteeringModel()
    p = ModelParams(kappa_h=kappa, kappa_g=kappa, delta=delta)
    c1_a, _ = bessel.taylor_coefficients(p)
    h = 1e-4
    Up = float(np.asarray(bessel.steering_drive(h, p)))
    Um = float(np.asarray(bessel.steering_drive(-h, p)))
    c1_fd = (Up - Um) / (2.0 * h)
    rel = abs(c1_a - c1_fd) / (abs(c1_fd) + 1e-9)
    assert rel < 1e-4, f"delta={delta}, kappa={kappa}: c1 a={c1_a}, fd={c1_fd}"


@pytest.mark.parametrize("delta", [0.5, 1.0, 1.3])
@pytest.mark.parametrize("kappa", [1.0, 2.0])
def test_c3_matches_finite_difference(delta, kappa):
    bessel = BesselSteeringModel()
    p = ModelParams(kappa_h=kappa, kappa_g=kappa, delta=delta)
    _, c3_a = bessel.taylor_coefficients(p)
    h = 5e-3
    U = lambda t: float(np.asarray(bessel.steering_drive(t, p)))
    c3_fd = (U(2 * h) - 2 * U(h) + 2 * U(-h) - U(-2 * h)) / (2.0 * h**3) / 6.0
    rel = abs(c3_a - c3_fd) / (abs(c3_fd) + 1e-9)
    assert rel < 5e-3, f"delta={delta}, kappa={kappa}: c3 a={c3_a}, fd={c3_fd}, rel={rel}"


def test_pitchfork_locus_sign_change():
    """c1 changes sign as delta crosses the pitchfork at fixed kappa."""
    bessel = BesselSteeringModel()
    deltas = np.linspace(1.2, 1.4, 21)
    c1s = []
    for d in deltas:
        c1, _ = bessel.taylor_coefficients(ModelParams(kappa_h=2.0, kappa_g=2.0, delta=d))
        c1s.append(c1)
    c1s = np.array(c1s)
    assert np.any(c1s < 0) and np.any(c1s > 0)
