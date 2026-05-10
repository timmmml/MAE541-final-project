"""Plot ``U(theta)`` and ``V(theta)`` for one or more models."""

from __future__ import annotations

import numpy as np

from steering.models.base import SteeringModel
from steering.params import ModelParams


def plot_steering_drive(
    models: list[SteeringModel],
    params_list: list[ModelParams] | ModelParams,
    theta_range: tuple[float, float] = (-np.pi, np.pi),
    n: int = 401,
    labels: list[str] | None = None,
    ax=None,
):
    """Overlay ``U(theta)`` curves for several models on one set of axes."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    if isinstance(params_list, ModelParams):
        params_list = [params_list] * len(models)
    theta = np.linspace(*theta_range, n)
    for i, (model, params) in enumerate(zip(models, params_list)):
        label = labels[i] if labels is not None else model.__class__.__name__
        U = np.asarray(model.steering_drive(theta, params))
        ax.plot(theta, U, label=label)
    ax.axhline(0.0, color="0.5", lw=0.5)
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(r"$U(\theta)$")
    ax.legend()
    return ax


def plot_potential(
    models: list[SteeringModel],
    params_list: list[ModelParams] | ModelParams,
    theta_range: tuple[float, float] = (-np.pi, np.pi),
    n: int = 401,
    labels: list[str] | None = None,
    ax=None,
):
    """Overlay ``V(theta)`` curves."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    if isinstance(params_list, ModelParams):
        params_list = [params_list] * len(models)
    theta = np.linspace(*theta_range, n)
    for i, (model, params) in enumerate(zip(models, params_list)):
        label = labels[i] if labels is not None else model.__class__.__name__
        V = np.asarray(model.steering_potential(theta, params))
        ax.plot(theta, V, label=label)
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(r"$V(\theta)$")
    ax.legend()
    return ax
