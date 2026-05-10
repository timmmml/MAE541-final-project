"""Phase portraits with optional Hamiltonian level sets and separatrices."""

from __future__ import annotations

import numpy as np

from steering.analysis.fixed_points import find_fixed_points_2d
from steering.analysis.hamiltonian import compute_separatrix, hamiltonian_from_model
from steering.dynamics.acceleration import AccelerationDynamics
from steering.models.base import SteeringModel
from steering.params import ModelParams


def plot_phase_portrait(
    dynamics: AccelerationDynamics,
    params: ModelParams,
    theta_range: tuple[float, float] = (-np.pi, np.pi),
    v_range: tuple[float, float] = (-3.0, 3.0),
    n_grid: int = 25,
    overlay_hamiltonian: bool = True,
    overlay_separatrix: bool = True,
    overlay_fixed_points: bool = True,
    n_contours: int = 18,
    trajectories: list[np.ndarray] | None = None,
    ax=None,
):
    """Quiver plot of the (theta, v) vector field with optional overlays."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(5.5, 4.5))

    theta = np.linspace(*theta_range, n_grid)
    v = np.linspace(*v_range, n_grid)
    TH, VV = np.meshgrid(theta, v)
    U = np.asarray(dynamics.model.steering_drive(TH, params))
    dTH = VV
    dVV = -dynamics.gamma * VV + U
    speed = np.sqrt(dTH**2 + dVV**2) + 1e-12
    ax.streamplot(theta, v, dTH, dVV, color=np.log(speed), cmap="Greys", density=1.2, linewidth=0.7)

    if overlay_hamiltonian:
        TH_f = np.linspace(*theta_range, 200)
        VV_f = np.linspace(*v_range, 200)
        TGRID, VGRID = np.meshgrid(TH_f, VV_f)
        H = hamiltonian_from_model(dynamics.model, params, TGRID, VGRID)
        ax.contour(TGRID, VGRID, H, levels=n_contours, colors="C0", alpha=0.4, linewidths=0.7)

    if overlay_fixed_points or overlay_separatrix:
        fps = find_fixed_points_2d(dynamics, params, theta_range=theta_range)
        for fp in fps:
            color = {
                "stable_node": "g",
                "stable_spiral": "g",
                "saddle": "k",
                "unstable_node": "r",
                "unstable_spiral": "r",
                "center": "b",
            }.get(fp.classification, "0.5")
            marker = "o" if "stable" in fp.classification else (
                "x" if "saddle" in fp.classification else "s"
            )
            ax.plot(fp.theta, 0.0, marker=marker, mec=color, mfc="white", ms=8, mew=1.6, zorder=5)
        if overlay_separatrix:
            for fp in fps:
                if "saddle" in fp.classification and dynamics.gamma == 0.0:
                    sep = compute_separatrix(dynamics, params, fp, t_span=20.0)
                    for key in ("unstable_plus", "unstable_minus", "stable_plus", "stable_minus"):
                        if key in sep:
                            arr = sep[key]
                            ax.plot(arr[:, 0], arr[:, 1], "k-", lw=1.0, alpha=0.7, zorder=4)

    if trajectories is not None:
        for traj in trajectories:
            ax.plot(traj[:, 0], traj[:, 1], "C3-", lw=0.9, alpha=0.7)

    ax.set_xlim(theta_range)
    ax.set_ylim(v_range)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$v = \dot\theta$")
    return ax
