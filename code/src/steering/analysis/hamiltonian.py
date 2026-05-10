"""Hamiltonian, potential, and separatrix tools for the undamped system."""

from __future__ import annotations

import numpy as np

from steering.analysis.fixed_points import FixedPoint2D
from steering.dynamics.acceleration import AccelerationDynamics
from steering.integrator import Simulation
from steering.models.base import SteeringModel
from steering.models.duffing import DuffingModel
from steering.params import ModelParams


def hamiltonian_from_model(
    model: SteeringModel,
    params: ModelParams,
    theta: np.ndarray,
    v: np.ndarray,
) -> np.ndarray:
    """Return ``H = 1/2 v^2 + V(theta)`` evaluated on a meshgrid.

    For :class:`DuffingModel` the potential is closed-form; otherwise we
    integrate ``-int_0^theta U(t) dt`` numerically.
    """
    theta = np.asarray(theta, dtype=float)
    v = np.asarray(v, dtype=float)
    if isinstance(model, DuffingModel):
        V = -0.5 * model.c1 * theta * theta - 0.25 * model.c3 * theta**4
    else:
        # Vectorize potential evaluation: integrate along a 1D grid in theta,
        # then interpolate. The grid spans the requested range.
        flat = np.atleast_1d(theta).ravel()
        if flat.size == 0:
            V = np.zeros_like(theta)
        else:
            t_min, t_max = float(flat.min()), float(flat.max())
            t_lo, t_hi = min(t_min, 0.0), max(t_max, 0.0)
            n = max(801, flat.size * 4)
            grid = np.linspace(t_lo, t_hi, n)
            U_grid = np.asarray(model.steering_drive(grid, params))
            # cumulative integral of U from grid[0] to grid; then shift so V(0) = 0.
            from scipy.integrate import cumulative_trapezoid

            int_U = np.concatenate([[0.0], cumulative_trapezoid(U_grid, grid)])
            # Find offset so V(0) = -int_0^0 U dt = 0.
            i0 = np.argmin(np.abs(grid))
            int_U = int_U - int_U[i0]
            V_flat = -np.interp(flat, grid, int_U)
            V = V_flat.reshape(theta.shape)
    return 0.5 * v * v + V


def compute_separatrix(
    dynamics: AccelerationDynamics,
    params: ModelParams,
    saddle: FixedPoint2D,
    t_span: float = 60.0,
    eps: float = 1e-7,
    direction: str = "both",
) -> dict[str, np.ndarray]:
    """Trace the stable / unstable manifolds of a 2D saddle.

    Integrates from ``saddle +- eps * e_unstable`` (forward, for the unstable
    manifold) and from ``saddle +- eps * e_stable`` (backward, for the stable
    manifold). Returns a dict of (theta, v) arrays.
    """
    if "saddle" not in saddle.classification:
        raise ValueError(f"compute_separatrix expects a saddle, got {saddle.classification!r}")

    J = dynamics.jacobian(saddle.theta, 0.0, params)
    eigvals, eigvecs = np.linalg.eig(J)
    # Choose the unstable (positive real) eigenvector.
    iu = int(np.argmax(eigvals.real))
    is_ = int(np.argmin(eigvals.real))
    eu = np.real(eigvecs[:, iu])
    eu = eu / np.linalg.norm(eu)
    es = np.real(eigvecs[:, is_])
    es = es / np.linalg.norm(es)

    sim = Simulation(dynamics, params, forcing=None, rtol=1e-10, atol=1e-12)
    out: dict[str, np.ndarray] = {}
    base = np.array([saddle.theta, 0.0])

    def _run(state0, span):
        return sim.run(state0, span)

    if direction in ("unstable", "both"):
        r1 = _run(base + eps * eu, (0.0, t_span))
        r2 = _run(base - eps * eu, (0.0, t_span))
        out["unstable_plus"] = r1.states
        out["unstable_minus"] = r2.states
        out["unstable_t"] = r1.t
    if direction in ("stable", "both"):
        # Stable manifold = unstable manifold in reverse time.
        r3 = _run(base + eps * es, (0.0, -t_span))
        r4 = _run(base - eps * es, (0.0, -t_span))
        out["stable_plus"] = r3.states
        out["stable_minus"] = r4.states
        out["stable_t"] = r3.t
    return out
