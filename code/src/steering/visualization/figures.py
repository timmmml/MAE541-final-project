"""Composite multi-panel figure builders for the report."""

from __future__ import annotations

from math import pi

import numpy as np

from steering.analysis.bifurcation import pitchfork_bifurcation_curve
from steering.analysis.fixed_points import find_fixed_points_2d
from steering.analysis.hamiltonian import hamiltonian_from_model
from steering.analysis.homoclinic import numerical_homoclinic
from steering.analysis.melnikov import (
    melnikov_critical_forcing,
    melnikov_critical_forcing_numerical,
)
from steering.analysis.poincare import bifurcation_diagram_poincare, stroboscopic_section
from steering.dynamics import AccelerationDynamics
from steering.integrator import Simulation
from steering.models import BesselSteeringModel, DuffingModel
from steering.params import ForcingParams, ModelParams
from steering.visualization.bifurcation_plot import plot_bifurcation_2d
from steering.visualization.phase_portrait import plot_phase_portrait
from steering.visualization.profiles import plot_steering_drive
from steering.visualization.style import use_paper_style


# ---------------------------------------------------------------------------
# Figure 1: pitchfork bifurcation in (kappa, delta) plane
# ---------------------------------------------------------------------------


def figure_1_bifurcation(
    base_params: ModelParams | None = None,
    kappa_range: tuple[float, float] = (0.5, 6.0),
    delta_range: tuple[float, float] = (0.0, pi / 2),
    n_kappa: int = 41,
    n_delta: int = 41,
    inset_kappa: float = 3.0,
    inset_deltas: tuple[float, float] = (0.5, 1.4),
):
    """Pitchfork in (kappa, delta) plane with two U(theta) insets."""
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    use_paper_style()
    if base_params is None:
        base_params = ModelParams(kappa_h=1.0, kappa_g=1.0, delta=0.5)

    bessel = BesselSteeringModel()
    kappas = np.linspace(*kappa_range, n_kappa)
    deltas = np.linspace(*delta_range, n_delta)

    # Use kappa_h as the sweep axis but tie kappa_g = kappa_h.
    # The diagram tool sweeps two named fields; we sweep kappa_h vs delta and
    # post-process by overwriting kappa_g = kappa_h via a wrapper params.
    c1_grid = np.empty((n_kappa, n_delta))
    c3_grid = np.empty_like(c1_grid)
    for i, k in enumerate(kappas):
        for j, d in enumerate(deltas):
            p = base_params.replace(kappa_h=float(k), kappa_g=float(k), delta=float(d))
            c1, c3 = bessel.taylor_coefficients(p)
            c1_grid[i, j] = c1
            c3_grid[i, j] = c3

    from steering.analysis.bifurcation import BifurcationDiagram

    diag = BifurcationDiagram(
        param1_name="kappa", param2_name="delta",
        param1_values=kappas, param2_values=deltas,
        c1_grid=c1_grid, c3_grid=c3_grid,
    )

    fig = plt.figure(figsize=(8.5, 5))
    gs = GridSpec(2, 3, width_ratios=[2.0, 0.05, 1.4], height_ratios=[1, 1], wspace=0.3, hspace=0.4)
    ax_main = fig.add_subplot(gs[:, 0])
    ax_a = fig.add_subplot(gs[0, 2])
    ax_b = fig.add_subplot(gs[1, 2])

    plot_bifurcation_2d(diag, param1_label=r"$\kappa$", param2_label=r"$\delta$", ax=ax_main)
    ax_main.set_title("Pitchfork locus in $(\\kappa,\\delta)$ plane")
    # Mark the insets' parameter points.
    for d in inset_deltas:
        ax_main.plot(inset_kappa, d, marker="*", color="k", ms=10, zorder=6)

    # Insets: U(theta) profiles.
    duffing_color = "#4477AA"
    bessel_color = "#EE6677"
    theta = np.linspace(-pi, pi, 401)
    for ax, d in zip([ax_a, ax_b], inset_deltas):
        p = base_params.replace(kappa_h=inset_kappa, kappa_g=inset_kappa, delta=d)
        c1, c3 = bessel.taylor_coefficients(p)
        Ub = np.asarray(bessel.steering_drive(theta, p))
        Ud = c1 * theta + c3 * theta ** 3
        ax.plot(theta, Ub, color=bessel_color, label="Bessel")
        ax.plot(theta, Ud, color=duffing_color, ls="--", label="Duffing")
        ax.axhline(0, color="0.6", lw=0.5)
        regime = "monostable" if c1 < 0 else "bistable"
        ax.set_title(rf"$\delta={d:.2f}$ ({regime})", fontsize=10)
        ax.set_xlabel(r"$\theta$")
        ax.set_ylabel(r"$U$")
        ax.set_xlim(-pi, pi)
    ax_a.legend(fontsize=8, loc="best")
    return fig, (ax_main, ax_a, ax_b)


# ---------------------------------------------------------------------------
# Figure 2: phase portraits Duffing vs Bessel, undamped vs damped
# ---------------------------------------------------------------------------


def figure_2_phase_portraits(
    params: ModelParams | None = None,
    gamma: float = 0.1,
    theta_range: tuple[float, float] = (-pi, pi),
    v_range: tuple[float, float] = (-6.0, 6.0),
):
    """2x2 grid: rows = (undamped, damped), cols = (Duffing, Bessel)."""
    import matplotlib.pyplot as plt

    use_paper_style()
    if params is None:
        params = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)

    bessel = BesselSteeringModel()
    duff = DuffingModel.from_params(params)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    titles = [
        (r"Duffing, $\gamma=0$", duff, 0.0),
        (r"Bessel, $\gamma=0$", bessel, 0.0),
        (rf"Duffing, $\gamma={gamma}$", duff, gamma),
        (rf"Bessel,  $\gamma={gamma}$", bessel, gamma),
    ]
    for ax, (title, model, g) in zip(axes.flat, titles):
        dyn = AccelerationDynamics(model=model, gamma=g, topology="planar")
        plot_phase_portrait(
            dyn,
            params,
            theta_range=theta_range,
            v_range=v_range,
            n_grid=22,
            n_contours=14,
            ax=ax,
        )
        ax.set_title(title)

    fig.tight_layout()
    return fig, axes


# ---------------------------------------------------------------------------
# Figure 3: Melnikov F_crit curves
# ---------------------------------------------------------------------------


def figure_3_melnikov(
    base_params: ModelParams | None = None,
    gamma: float = 0.05,
    omega_range: tuple[float, float] = (0.2, 4.0),
    n_omega: int = 40,
    delta_range: tuple[float, float] = (1.32, 1.55),
    n_delta: int = 18,
    omega_fixed: float = 1.0,
):
    """Two panels: F_crit(omega) at fixed delta, F_crit(delta) at fixed omega."""
    import matplotlib.pyplot as plt

    use_paper_style()
    if base_params is None:
        base_params = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)

    bessel = BesselSteeringModel()

    # Panel A: F_crit vs omega.
    omegas = np.linspace(*omega_range, n_omega)
    c1, c3 = bessel.taylor_coefficients(base_params)
    F_a_duff = np.array([melnikov_critical_forcing(c1, c3, gamma, om) for om in omegas])
    homo = numerical_homoclinic(bessel, base_params)
    F_a_bess = np.array(
        [melnikov_critical_forcing_numerical(bessel, base_params, gamma, om, homoclinic=homo) for om in omegas]
    )

    # Panel B: F_crit vs delta at omega_fixed.
    deltas = np.linspace(*delta_range, n_delta)
    F_b_duff = np.empty_like(deltas)
    F_b_bess = np.empty_like(deltas)
    for i, d in enumerate(deltas):
        p = base_params.replace(delta=float(d))
        c1d, c3d = bessel.taylor_coefficients(p)
        if c1d > 0 and c3d < 0:
            F_b_duff[i] = melnikov_critical_forcing(c1d, c3d, gamma, omega_fixed)
            try:
                F_b_bess[i] = melnikov_critical_forcing_numerical(bessel, p, gamma, omega_fixed)
            except Exception:
                F_b_bess[i] = np.nan
        else:
            F_b_duff[i] = np.nan
            F_b_bess[i] = np.nan

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(10, 4))
    axA.plot(omegas, F_a_duff, color="#4477AA", lw=1.8, label="Duffing (analytical)")
    axA.plot(omegas, F_a_bess, color="#EE6677", lw=1.8, ls="--", label="Bessel (numerical)")
    axA.set_xlabel(r"$\omega$")
    axA.set_ylabel(r"$F_{\mathrm{crit}}$")
    axA.set_yscale("log")
    axA.legend()
    axA.set_title(rf"$F_{{\mathrm{{crit}}}}(\omega)$ at $\delta={base_params.delta:.2f}$")

    axB.plot(deltas, F_b_duff, color="#4477AA", lw=1.8, label="Duffing")
    axB.plot(deltas, F_b_bess, color="#EE6677", lw=1.8, ls="--", label="Bessel")
    axB.set_xlabel(r"$\delta$")
    axB.set_ylabel(r"$F_{\mathrm{crit}}$")
    axB.legend()
    axB.set_title(rf"$F_{{\mathrm{{crit}}}}(\delta)$ at $\omega={omega_fixed}$")

    fig.tight_layout()
    return fig, (axA, axB)


# ---------------------------------------------------------------------------
# Figure 4: Poincaré bifurcation + sample sections
# ---------------------------------------------------------------------------


def figure_4_poincare(
    base_params: ModelParams | None = None,
    gamma: float = 0.05,
    omega: float = 1.0,
    F_range: tuple[float, float] = (0.0, 1.5),
    n_F: int = 80,
    F_section_examples: tuple[float, ...] = (0.05, 0.5, 1.0),
    n_periods_transient: int = 200,
    n_periods_record: int = 80,
):
    """Top: bifurcation diagram theta_strobe vs F. Bottom: example sections."""
    import matplotlib.pyplot as plt

    use_paper_style()
    if base_params is None:
        base_params = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)

    bessel = BesselSteeringModel()
    dyn = AccelerationDynamics(model=bessel, gamma=gamma, topology="planar")
    F_vals = np.linspace(*F_range, n_F)
    diag = bifurcation_diagram_poincare(
        dyn, bessel, base_params, gamma=gamma, omega=omega,
        sweep_param="F", sweep_values=F_vals,
        F=None, state0=np.array([0.5, 0.0]),
        n_periods_transient=n_periods_transient,
        n_periods_record=n_periods_record,
    )

    n_examples = len(F_section_examples)
    fig = plt.figure(figsize=(10, 7))
    gs = fig.add_gridspec(2, n_examples, height_ratios=[1.2, 1.0], hspace=0.4, wspace=0.3)

    ax_top = fig.add_subplot(gs[0, :])
    for v, vals in zip(diag["sweep_values"], diag["thetas"]):
        ax_top.plot(np.full_like(vals, v), vals, "k.", ms=0.7, alpha=0.5)
    ax_top.set_xlabel(r"forcing amplitude $F$")
    ax_top.set_ylabel(r"$\theta$ (stroboscopic)")
    ax_top.set_title(r"Poincaré bifurcation diagram, Bessel model")

    for j, F in enumerate(F_section_examples):
        ax = fig.add_subplot(gs[1, j])
        forc = ForcingParams(F=F, omega=omega)
        sim = Simulation(dyn, base_params, forc, rtol=1e-8, atol=1e-10)
        # Long enough to populate the section.
        T = 2 * pi / omega
        n_record_section = max(300, n_periods_record * 4)
        result = sim.run(
            np.array([0.5, 0.0]),
            (0.0, (n_periods_transient + n_record_section) * T),
            dense_output=True,
        )
        sec = stroboscopic_section(result, omega=omega, transient_periods=n_periods_transient)
        ax.plot(sec[:, 0], sec[:, 1], "k.", ms=1.0, alpha=0.6)
        ax.set_xlabel(r"$\theta$")
        if j == 0:
            ax.set_ylabel(r"$v$")
        ax.set_title(rf"$F={F}$")

    return fig, ax_top
