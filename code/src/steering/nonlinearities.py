"""Nonlinearity registry for PFL population activity ``f(h)``."""

from __future__ import annotations

from typing import Callable

import numpy as np


def quadratic(x: np.ndarray) -> np.ndarray:
    return x * x


def relu(x: np.ndarray) -> np.ndarray:
    return np.maximum(0.0, x)


def elu(x: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    x = np.asarray(x)
    return np.where(x >= 0, x, alpha * (np.exp(x) - 1.0))


def softplus(x: np.ndarray, beta: float = 1.0) -> np.ndarray:
    x = np.asarray(x)
    # Numerically-stable softplus.
    z = beta * x
    return np.where(z > 30.0, z, np.log1p(np.exp(z))) / beta


_REGISTRY: dict[str, Callable[..., np.ndarray]] = {
    "quadratic": quadratic,
    "relu": relu,
    "elu": elu,
    "softplus": softplus,
}


def get_nonlinearity(
    name_or_callable: str | Callable, params: dict | None = None
) -> Callable[[np.ndarray], np.ndarray]:
    """Resolve a nonlinearity spec to a unary function ``f(x)``.

    ``name_or_callable`` may be a registry key, or a callable taking an array
    and any extra keyword arguments listed in ``params``.
    """
    if callable(name_or_callable):
        f = name_or_callable
    else:
        try:
            f = _REGISTRY[name_or_callable]
        except KeyError as exc:
            raise ValueError(
                f"Unknown nonlinearity {name_or_callable!r}. "
                f"Known: {sorted(_REGISTRY)}"
            ) from exc

    extra = dict(params or {})
    if not extra:
        return f

    def wrapped(x: np.ndarray, _f=f, _kw=extra) -> np.ndarray:
        return _f(x, **_kw)

    return wrapped


def register(name: str, fn: Callable[..., np.ndarray]) -> None:
    """Add a new nonlinearity to the registry."""
    _REGISTRY[name] = fn
