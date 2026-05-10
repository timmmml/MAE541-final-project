"""Stroboscopic Poincaré section plots."""

from __future__ import annotations

import numpy as np


def plot_poincare_section(section: np.ndarray, xlabel: str = r"$\theta$", ylabel: str = r"$v$", ax=None, **kwargs):
    """Scatter plot of ``(theta, v)`` strobed at the forcing period.

    ``section`` is the array returned by ``stroboscopic_section`` (shape
    ``(n_crossings, state_dim)``).
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))
    s = kwargs.pop("s", 1.0)
    ax.plot(section[:, 0], section[:, 1], "k.", ms=2.0, alpha=0.5, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax
