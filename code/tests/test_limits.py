"""Limiting-case tests: kappa -> 0, N -> infinity, etc."""

from __future__ import annotations

import numpy as np
import pytest

from steering.models import BesselSteeringModel, ContinuousPFLModel, DiscretePFLModel
from steering.params import ModelParams


def test_small_kappa_limit():
    """For kappa -> 0, c1 ~ -8 pi S A * 2 sin(Delta) cos(delta) (broad-bump limit)."""
    bessel = BesselSteeringModel()
    Delta = 3 * np.pi / 8
    delta = np.pi / 6
    p = ModelParams(kappa_h=0.05, kappa_g=0.05, delta=delta, Delta_pop=Delta)
    c1, _ = bessel.taylor_coefficients(p)
    # In the broad-bump limit: I_1(K)/K ~ 1/2 + K^2/16 + ...
    # so c1 -> -8 pi S A kappa_h kappa_g * (1/2)(sin(Delta - delta) + sin(Delta + delta))
    #       = -8 pi kappa_h kappa_g * sin(Delta) cos(delta)
    #       = -8 pi (0.05)^2 * sin(3pi/8) cos(pi/6).
    expected = -8.0 * np.pi * 0.05 * 0.05 * np.sin(Delta) * np.cos(delta)
    assert abs(c1 - expected) / abs(expected) < 5e-3


def test_discrete_to_continuous_convergence():
    """L^infinity error -> 0 as N -> infinity."""
    cont = ContinuousPFLModel(n_quad=256)
    disc = DiscretePFLModel()
    p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=np.pi / 4)
    theta = np.linspace(-np.pi + 0.05, np.pi - 0.05, 21)
    errs = []
    for N in [4, 8, 16, 32, 64, 128]:
        p_N = p.replace(N_neurons=N)
        err = np.max(np.abs(cont.steering_drive(theta, p) - disc.steering_drive(theta, p_N)))
        errs.append(err)
    # Error should monotonically decrease (or be near machine precision at the end).
    assert errs[-1] < 1e-8
    # Strict monotone decrease in the rapid-convergence regime (N up to 32).
    assert all(errs[i] >= errs[i + 1] - 1e-12 for i in range(3))
