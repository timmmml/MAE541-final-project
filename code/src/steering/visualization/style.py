"""Shared style for plots: color palette, line styles, and small helpers."""

from __future__ import annotations

import matplotlib as mpl

# Colorblind-friendly palette (Tol's bright). Shared across all figures so the
# same model has the same color in every panel.
PALETTE = {
    "duffing": "#4477AA",      # blue
    "bessel": "#EE6677",       # red
    "continuous": "#228833",   # green
    "discrete": "#CCBB44",     # yellow
    "full_circuit": "#AA3377", # purple
    "stable": "#228833",
    "unstable": "#EE6677",
    "saddle": "#000000",
    "separatrix": "#222222",
    "trajectory": "#4477AA",
}


LINE_STYLES = {
    "duffing": "--",
    "bessel": "-",
    "continuous": "-",
    "discrete": ":",
    "full_circuit": "-.",
}


def use_paper_style() -> None:
    """Apply a tidy paper-figure default to matplotlib's rcParams."""
    mpl.rcParams.update({
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "lines.linewidth": 1.6,
    })


def color_for(model_key: str) -> str:
    return PALETTE.get(model_key, "#444444")


def linestyle_for(model_key: str) -> str:
    return LINE_STYLES.get(model_key, "-")
