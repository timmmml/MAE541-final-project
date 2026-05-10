"""Maximal Lyapunov exponent via the standard tangent-vector renormalization."""

from __future__ import annotations

from math import log, pi

import numpy as np
from scipy.integrate import solve_ivp

from steering.dynamics.acceleration import AccelerationDynamics
from steering.models.base import SteeringModel
from steering.params import ForcingParams, ModelParams


def _augmented_rhs(dyn: AccelerationDynamics, params: ModelParams, forcing: ForcingParams | None):
    """Combined RHS: (theta, v) for the trajectory and (delta_theta, delta_v) for the variational."""

    def f(t, z):
        x = z[:2]
        dx = z[2:]
        dxdt = dyn.rhs(t, x, params, forcing)
        # Variational equation: d(delta_x)/dt = J(x, t) delta_x.
        J = dyn.jacobian(float(x[0]), float(x[1]), params)
        ddxdt = J @ dx
        return np.concatenate([dxdt, ddxdt])

    return f


def max_lyapunov_exponent(
    dynamics: AccelerationDynamics,
    params: ModelParams,
    gamma: float,
    F: float,
    omega: float,
    state0: np.ndarray,
    t_total: float = 1000.0,
    renorm_interval: float = 1.0,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    forcing_type: str = "additive",
) -> float:
    """Estimate ``lambda_max`` of the forced damped system.

    Implements Benettin's algorithm: integrate trajectory + tangent vector,
    renormalize the tangent every ``renorm_interval`` units of time and
    accumulate ``log(|delta|)``.
    """
    dyn = AccelerationDynamics(model=dynamics.model, gamma=gamma, topology=dynamics.topology)
    forcing = ForcingParams(F=F, omega=omega, forcing_type=forcing_type)
    f = _augmented_rhs(dyn, params, forcing)

    delta = np.array([1.0, 0.0])
    delta = delta / np.linalg.norm(delta)
    z = np.concatenate([np.asarray(state0, dtype=float), delta])
    t = 0.0
    log_sum = 0.0
    n_renorm = 0
    while t < t_total:
        t_next = min(t + renorm_interval, t_total)
        sol = solve_ivp(f, (t, t_next), z, rtol=rtol, atol=atol, method="RK45")
        z = sol.y[:, -1]
        d = z[2:]
        norm = np.linalg.norm(d)
        if norm <= 0:
            break
        log_sum += log(norm)
        z[2:] = d / norm
        t = t_next
        n_renorm += 1
    if n_renorm == 0 or t == 0.0:
        return 0.0
    return log_sum / t
