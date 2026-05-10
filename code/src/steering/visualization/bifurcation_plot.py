"""Bifurcation plots."""

from __future__ import annotations

import numpy as np

from steering.analysis.bifurcation import BifurcationDiagram


def plot_bifurcation_2d(
    diagram: BifurcationDiagram,
    param1_label: str | None = None,
    param2_label: str | None = None,
    overlay_numerical: dict | None = None,
    ax=None,
):
    """Plot the pitchfork locus ``c_1 = 0`` in a 2D parameter plane.

    The locus is drawn as a contour, and the regions where ``c_3 > 0`` (
    subcritical) vs ``c_3 < 0`` (supercritical) are shaded distinctly.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    p1 = diagram.param1_values
    p2 = diagram.param2_values
    P1, P2 = np.meshgrid(p1, p2, indexing="ij")
    # Shading by sgn(c_3): -1 for supercritical (c3<0), +1 for subcritical (c3>0).
    sign = np.where(diagram.c3_grid < 0.0, -1.0, 1.0)
    cmap = ListedColormap(["#cce6cc", "#fde0e0"])
    ax.pcolormesh(P1, P2, sign, cmap=cmap, alpha=0.8, shading="auto", vmin=-1.0, vmax=1.0)
    cs = ax.contour(P1, P2, diagram.c1_grid, levels=[0.0], colors="k", linewidths=2.0)
    ax.clabel(cs, fmt={0.0: r"$c_1 = 0$"}, fontsize=9)
    if overlay_numerical is not None:
        for key, (xs, ys) in overlay_numerical.items():
            ax.plot(xs, ys, label=key, ls="--")
        ax.legend()
    ax.set_xlabel(param1_label or diagram.param1_name)
    ax.set_ylabel(param2_label or diagram.param2_name)
    return ax


def plot_poincare_bifurcation(
    data: dict,
    sweep_label: str | None = None,
    component: str = "theta",
    ax=None,
    s: float = 1.0,
):
    """Plot the stroboscopic Poincaré bifurcation diagram (theta vs sweep param)."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    series = data["thetas"] if component == "theta" else data["vs"]
    for v, vals in zip(data["sweep_values"], series):
        ax.plot(np.full_like(vals, v), vals, "k.", ms=1, alpha=0.5)
    ax.set_xlabel(sweep_label or data["sweep_param"])
    ax.set_ylabel(r"$\theta_{\mathrm{strobe}}$" if component == "theta" else r"$v_{\mathrm{strobe}}$")
    return ax
