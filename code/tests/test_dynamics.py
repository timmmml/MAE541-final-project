"""Dynamics and integrator tests."""

from __future__ import annotations

import numpy as np

from steering.analysis.hamiltonian import hamiltonian_from_model
from steering.dynamics import AccelerationDynamics, VelocityDynamics
from steering.integrator import Simulation
from steering.models import BesselSteeringModel, DuffingModel
from steering.params import ForcingParams, ModelParams


def test_hamiltonian_conservation_undamped(params_weakly_bistable):
    """§9.3 (5): undamped, unforced AccelerationDynamics conserves H."""
    bessel = BesselSteeringModel()
    dyn = AccelerationDynamics(model=bessel, gamma=0.0)
    sim = Simulation(dyn, params_weakly_bistable, rtol=1e-11, atol=1e-13)
    res = sim.run(np.array([0.5, 0.5]), (0.0, 30.0))
    H = hamiltonian_from_model(bessel, params_weakly_bistable, res.states[:, 0], res.states[:, 1])
    drift = np.max(np.abs(H - H[0]))
    assert drift / (abs(H[0]) + 1e-3) < 1e-5


def test_velocity_dynamics_basic(params_monostable):
    """1D velocity dynamics drives theta -> 0 (the stable monostable fixed point)."""
    bessel = BesselSteeringModel()
    dyn = VelocityDynamics(model=bessel, topology="planar")
    sim = Simulation(dyn, params_monostable, rtol=1e-9, atol=1e-11)
    res = sim.run(np.array([0.3]), (0.0, 5.0))
    assert abs(res.states[-1, 0]) < 1e-3


def test_topology_wrapping():
    """Cylindrical topology wraps theta to [-pi, pi)."""
    duff = DuffingModel(c1=1.0, c3=-1.0)  # has saddles at 0, +-1.
    dyn = AccelerationDynamics(model=duff, gamma=0.0, topology="cylindrical")
    state = np.array([3.5, 0.0])  # theta > pi.
    wrapped = dyn.wrap_state(state)
    assert -np.pi <= wrapped[0] < np.pi
    assert abs(wrapped[1] - state[1]) < 1e-12


def test_forcing_additive_acceleration():
    """Additive forcing adds F cos(omega t) to v_dot."""
    duff = DuffingModel(c1=-1.0, c3=-1.0)  # monostable
    dyn = AccelerationDynamics(model=duff, gamma=0.5)
    forc = ForcingParams(F=0.3, omega=1.0)
    rhs = dyn.rhs(0.0, np.array([0.0, 0.0]), ModelParams(kappa_h=1.0, kappa_g=1.0, delta=0.0), forc)
    # At theta=0, v=0: theta_dot = 0; v_dot = -gamma*0 + U(0) + F*1 = 0 + 0 + 0.3.
    assert abs(rhs[0]) < 1e-12
    assert abs(rhs[1] - 0.3) < 1e-12
