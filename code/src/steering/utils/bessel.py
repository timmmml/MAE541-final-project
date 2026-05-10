"""Overflow-safe wrappers around modified Bessel functions of the first kind.

The steering drive involves differences of ``I_0`` evaluated at similar
arguments. For ``z`` of order a few hundred ``scipy.special.i0(z)`` overflows;
the exponentially scaled ``i0e(z) = e^{-|z|} I_0(z)`` does not. Since the four
arguments entering ``U(theta)`` are typically of comparable size, the leading
``e^{|z|}`` prefactors approximately cancel in the difference. We expose
``i0_safe`` etc. that fall back on the scaled forms above a threshold.
"""

from __future__ import annotations

import numpy as np
from scipy import special

# Threshold above which we switch to exp-scaled Bessel functions.
_SCALE_THRESHOLD = 500.0


def _safe(z, fn, fn_e):
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z, dtype=float)
    small = np.abs(z) < _SCALE_THRESHOLD
    out[small] = fn(z[small])
    big = ~small
    if np.any(big):
        out[big] = fn_e(z[big]) * np.exp(np.abs(z[big]))
    return out


def i0_safe(z):
    """``I_0(z)`` with overflow protection."""
    z = np.asarray(z, dtype=float)
    if z.ndim == 0:
        return float(_safe(z[None], special.i0, special.i0e)[0])
    return _safe(z, special.i0, special.i0e)


def i1_safe(z):
    """``I_1(z)`` with overflow protection."""
    z = np.asarray(z, dtype=float)
    if z.ndim == 0:
        return float(_safe(z[None], special.i1, special.i1e)[0])
    return _safe(z, special.i1, special.i1e)


def iv_safe(n, z):
    """``I_n(z)`` for integer ``n`` with overflow protection."""
    z = np.asarray(z, dtype=float)

    def fn(zz):
        return special.iv(n, zz)

    def fn_e(zz):
        return special.ive(n, zz)

    if z.ndim == 0:
        return float(_safe(z[None], fn, fn_e)[0])
    return _safe(z, fn, fn_e)


def i_n_ratio(n: int, z):
    """``I_n(z) / I_{n-1}(z)`` computed with the scaled forms (well behaved).

    Useful when explicit ``I_n`` would overflow but only the ratio is needed.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    z = np.asarray(z, dtype=float)
    num = special.ive(n, z)
    den = special.ive(n - 1, z)
    return num / den
