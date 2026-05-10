"""Cross-model consistency tests (spec §9.3)."""

from __future__ import annotations

import numpy as np
import pytest

from steering.models import (
    BesselSteeringModel,
    ContinuousPFLModel,
    DiscretePFLModel,
    DuffingModel,
    FullCircuitModel,
)
from steering.params import ModelParams


@pytest.fixture
def theta_grid():
    return np.linspace(-np.pi + 0.01, np.pi - 0.01, 33)


def test_bessel_equals_continuous_quadratic(params_monostable, theta_grid):
    """§9.3 (1): Bessel == Continuous PFL with quadratic f."""
    bessel = BesselSteeringModel()
    cont = ContinuousPFLModel(n_quad=128)
    Ub = bessel.steering_drive(theta_grid, params_monostable)
    Uc = cont.steering_drive(theta_grid, params_monostable)
    assert np.max(np.abs(Ub - Uc)) < 1e-10


def test_duffing_matches_bessel_near_origin(params_close_to_bif):
    """§9.3 (2): Duffing.from_params matches Bessel for |theta| << 1."""
    bessel = BesselSteeringModel()
    duff = DuffingModel.from_params(params_close_to_bif)
    theta_small = np.linspace(-0.1, 0.1, 21)
    Ub = bessel.steering_drive(theta_small, params_close_to_bif)
    Ud = duff.steering_drive(theta_small)
    rel = np.max(np.abs(Ub - Ud)) / (np.max(np.abs(Ub)) + 1e-12)
    assert rel < 1e-3


def test_discrete_converges_to_continuous(params_monostable, theta_grid):
    """§9.3 (3): DiscretePFL with large N converges to ContinuousPFL."""
    cont = ContinuousPFLModel(n_quad=128)
    disc = DiscretePFLModel()
    Uc = cont.steering_drive(theta_grid, params_monostable)
    p_big_N = params_monostable.replace(N_neurons=512)
    Ud = disc.steering_drive(theta_grid, p_big_N)
    # 512 trapezoid nodes for a smooth periodic integrand should be near
    # exact (exponentially converging trapezoid for analytic periodic
    # integrands).
    assert np.max(np.abs(Uc - Ud)) < 1e-8


def test_full_circuit_reduces_with_W_D3_zero(params_monostable, theta_grid):
    """§9.3 (4): FullCircuit(W_D3=0) reduces to ContinuousPFL/DiscretePFL."""
    full = FullCircuitModel()
    cont = ContinuousPFLModel(n_quad=256)
    Uf = full.steering_drive(theta_grid, params_monostable)
    Uc = cont.steering_drive(theta_grid, params_monostable)
    assert np.allclose(Uf, Uc, atol=1e-10)
    # And in discrete mode.
    p_disc = params_monostable.replace(N_neurons=8)
    Uf_d = full.steering_drive(theta_grid, p_disc)
    disc = DiscretePFLModel()
    Ud = disc.steering_drive(theta_grid, p_disc)
    assert np.allclose(Uf_d, Ud, atol=1e-12)


def test_taylor_coefficients_finite_difference(params_close_to_bif):
    """Bessel.taylor_coefficients matches a finite-difference computation."""
    bessel = BesselSteeringModel()
    c1, c3 = bessel.taylor_coefficients(params_close_to_bif)
    h = 1e-3
    # Higher-order central FD for the first and third derivatives.
    U = lambda t: float(np.asarray(bessel.steering_drive(t, params_close_to_bif)))
    c1_fd = (U(h) - U(-h)) / (2.0 * h)
    c3_fd = (U(2 * h) - 2 * U(h) + 2 * U(-h) - U(-2 * h)) / (2.0 * h**3) / 6.0
    assert abs(c1 - c1_fd) / abs(c1_fd + 1e-12) < 1e-4
    assert abs(c3 - c3_fd) / abs(c3_fd + 1e-12) < 1e-2
