"""Shared pytest fixtures."""

from __future__ import annotations

import numpy as np
import pytest

from steering.params import ModelParams


@pytest.fixture
def params_monostable():
    return ModelParams(kappa_h=2.0, kappa_g=2.0, delta=np.pi / 4)


@pytest.fixture
def params_weakly_bistable():
    # delta = 1.4 places kappa = 2 just past the pitchfork (c1 ~ 22).
    return ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)


@pytest.fixture
def params_close_to_bif():
    # Just barely bistable: theta_max ~ 0.5.
    return ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.34)
