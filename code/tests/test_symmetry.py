"""Symmetry tests: U(-theta) = -U(theta), U(0) = 0."""

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


@pytest.mark.parametrize(
    "model_factory",
    [
        BesselSteeringModel,
        ContinuousPFLModel,
        lambda: DuffingModel(c1=10.0, c3=-5.0),
        FullCircuitModel,
    ],
)
def test_odd_symmetry(model_factory, params_weakly_bistable):
    model = model_factory()
    if isinstance(model, DuffingModel):
        Up = model.steering_drive(0.5)
        Um = model.steering_drive(-0.5)
        assert abs(Up + Um) < 1e-12
        assert abs(model.steering_drive(0.0)) < 1e-12
    else:
        theta = np.array([-0.7, -0.3, 0.0, 0.3, 0.7, 1.1])
        Up = np.asarray(model.steering_drive(theta, params_weakly_bistable))
        Um = np.asarray(model.steering_drive(-theta, params_weakly_bistable))
        assert np.max(np.abs(Up + Um)) < 1e-9


def test_discrete_symmetry_for_even_N(params_weakly_bistable):
    """Even-N discrete grid preserves theta -> -theta symmetry exactly."""
    disc = DiscretePFLModel()
    theta = np.array([-1.0, -0.4, 0.0, 0.4, 1.0])
    p = params_weakly_bistable.replace(N_neurons=8)
    Up = disc.steering_drive(theta, p)
    Um = disc.steering_drive(-theta, p)
    assert np.max(np.abs(Up + Um)) < 1e-10
