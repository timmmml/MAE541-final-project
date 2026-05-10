"""Topology-handling tests for cylindrical S^1 wrap."""

from __future__ import annotations

import numpy as np

from steering.dynamics import AccelerationDynamics, VelocityDynamics
from steering.integrator import Simulation
from steering.models import DuffingModel
from steering.params import ModelParams


def test_velocity_cylindrical_wraps():
    """A trajectory leaving (-pi, pi) is wrapped post-step."""
    duff = DuffingModel(c1=10.0, c3=-1.0)  # bistable, large drive.
    dyn = VelocityDynamics(model=duff, topology="cylindrical")
    sim = Simulation(dyn, ModelParams(kappa_h=1.0, kappa_g=1.0, delta=0.0))
    # Start near the unstable saddle at theta=pi - tiny epsilon and watch
    # it flow back to a well; check that values stay in [-pi, pi).
    res = sim.run(np.array([np.pi - 0.01]), (0.0, 5.0))
    assert np.all(res.states[:, 0] >= -np.pi)
    assert np.all(res.states[:, 0] < np.pi)


def test_acceleration_planar_no_wrap():
    """Planar topology leaves the state as-is (no wrap)."""
    duff = DuffingModel(c1=-1.0, c3=-1.0)
    dyn = AccelerationDynamics(model=duff, gamma=0.0, topology="planar")
    state = np.array([5.0, 1.0])
    out = dyn.wrap_state(state)
    assert np.allclose(out, state)
