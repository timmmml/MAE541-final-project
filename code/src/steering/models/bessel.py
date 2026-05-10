"""Exact closed-form steering drive in terms of modified Bessel functions.

Valid for the continuous-population PFL model with quadratic nonlinearity.
"""

from __future__ import annotations

from math import pi

import numpy as np

from steering.models.base import SteeringModel
from steering.params import ModelParams
from steering.utils.bessel import i0_safe, i1_safe, iv_safe


def _kappa(a, b, kh: float, kg: float):
    """Combined concentration ``sqrt(kh^2 + kg^2 + 2 kh kg cos(a - b))``."""
    return np.sqrt(kh * kh + kg * kg + 2.0 * kh * kg * np.cos(a - b))


class BesselSteeringModel(SteeringModel):
    r"""Exact ``U(theta, delta)`` from the four-Bessel closed form.

    ``U(theta, delta) = 4 pi S A * [ I_0(kappa(theta+Delta, delta))
                                  + I_0(kappa(theta+Delta, -delta))
                                  - I_0(kappa(theta-Delta, delta))
                                  - I_0(kappa(theta-Delta, -delta)) ]``
    """

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _four_kappas(theta, params: ModelParams):
        kh, kg, delta, Delta = (
            params.kappa_h,
            params.kappa_g,
            params.delta,
            params.Delta_pop,
        )
        theta = np.asarray(theta, dtype=float)
        a_p = theta + Delta
        a_m = theta - Delta
        K_pp = _kappa(a_p, delta, kh, kg)
        K_pm = _kappa(a_p, -delta, kh, kg)
        K_mp = _kappa(a_m, delta, kh, kg)
        K_mm = _kappa(a_m, -delta, kh, kg)
        return a_p, a_m, K_pp, K_pm, K_mp, K_mm

    # ------------------------------------------------------------------
    # SteeringModel interface
    # ------------------------------------------------------------------
    def steering_drive(self, theta, params: ModelParams):
        _, _, K_pp, K_pm, K_mp, K_mm = self._four_kappas(theta, params)
        kh, kg = params.kappa_h, params.kappa_g
        prefactor = 4.0 * pi * params.S * params.A
        # Pure-product-of-bumps terms cancel between L and R; what remains is
        # the cross-term in (S g_h + A g_g)^2 integrated over phi -> 2pi I_0.
        # The factor 2 pi got absorbed into the 4 pi prefactor (factor of 2
        # is from the (S g_h + A g_g)^2 cross term).
        out = (
            i0_safe(K_pp) + i0_safe(K_pm) - i0_safe(K_mp) - i0_safe(K_mm)
        )
        return prefactor * out

    def steering_drive_derivatives(
        self, theta, params: ModelParams, order: int = 3, h: float = 1e-3
    ):
        """Analytical first derivative; higher orders fall back to finite diff."""
        if order >= 4:
            raise ValueError("order must be <= 3")
        theta = np.asarray(theta, dtype=float)

        # First derivative via chain rule: d I_0(K)/dx = I_1(K) dK/dx.
        # K(a, b) = sqrt(kh^2 + kg^2 + 2 kh kg cos(a - b))
        # dK/d theta = - kh kg sin(a - b) / K  (since da/dtheta = +-1)
        kh, kg, delta, Delta = (
            params.kappa_h,
            params.kappa_g,
            params.delta,
            params.Delta_pop,
        )
        a_p = theta + Delta
        a_m = theta - Delta
        # arguments of cos / sin
        d1 = a_p - delta
        d2 = a_p + delta
        d3 = a_m - delta
        d4 = a_m + delta
        K1 = np.sqrt(kh * kh + kg * kg + 2 * kh * kg * np.cos(d1))
        K2 = np.sqrt(kh * kh + kg * kg + 2 * kh * kg * np.cos(d2))
        K3 = np.sqrt(kh * kh + kg * kg + 2 * kh * kg * np.cos(d3))
        K4 = np.sqrt(kh * kh + kg * kg + 2 * kh * kg * np.cos(d4))

        prefactor = 4.0 * pi * params.S * params.A
        U0 = prefactor * (i0_safe(K1) + i0_safe(K2) - i0_safe(K3) - i0_safe(K4))

        # d/dtheta(I_0(K)) = I_1(K) * dK/dtheta = I_1(K) * (- kh*kg*sin(d)/K) * (da/dtheta)
        # for terms 1,2 da/dtheta = +1, for terms 3,4 da/dtheta = +1 too (a_m = theta - Delta has d/dtheta = +1)
        # but the - sign in front comes from the U formula (terms 3,4 are subtracted).
        def chain(K, d):
            # safeguard against K = 0 (only happens if kh = kg and cos(d) = -1)
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio = np.where(K > 0, i1_safe(K) / K, 0.5)
            return -kh * kg * np.sin(d) * ratio

        dU_dtheta = prefactor * (chain(K1, d1) + chain(K2, d2) - chain(K3, d3) - chain(K4, d4))

        out = [U0, dU_dtheta]
        if order >= 2:
            # Higher derivatives via finite difference of the analytic first
            # derivative for accuracy.
            def Up(t):
                _, dU = self.steering_drive_derivatives(t, params, order=1)
                return dU

            out.append((Up(theta + h) - Up(theta - h)) / (2.0 * h))
        if order >= 3:
            def Up(t):
                _, dU = self.steering_drive_derivatives(t, params, order=1)
                return dU

            out.append(
                (Up(theta + h) - 2 * Up(theta) + Up(theta - h)) / (h * h)
            )
        return tuple(out[: order + 1])

    # ------------------------------------------------------------------
    # Taylor coefficients at theta = 0
    # ------------------------------------------------------------------
    def taylor_coefficients(self, params: ModelParams) -> tuple[float, float]:
        """Return ``(c1, c3)`` so that ``U(theta) ~ c1 theta + c3 theta^3``."""
        kh, kg = params.kappa_h, params.kappa_g
        delta, Delta = params.delta, params.Delta_pop
        S, A = params.S, params.A

        def K(arg):
            return float(np.sqrt(kh * kh + kg * kg + 2 * kh * kg * np.cos(arg)))

        d_minus = Delta - delta
        d_plus = Delta + delta
        Km = K(d_minus)
        Kp = K(d_plus)

        # c1 = U'(0) (boxed formula in spec).
        def safe_div(a, b):
            return a / b if b > 0 else 0.0

        c1 = -8.0 * pi * S * A * kh * kg * (
            safe_div(i1_safe(Km), Km) * np.sin(d_minus)
            + safe_div(i1_safe(Kp), Kp) * np.sin(d_plus)
        )

        # c3 = U'''(0) / 6 with U'''(0) given by the boxed D(K, C, S) formula.
        def D(K_, Cv, Sv):
            if K_ <= 0:
                return 0.0
            I0 = float(i0_safe(K_))
            I1 = float(i1_safe(K_))
            I2 = float(iv_safe(2, K_))
            I3 = float(iv_safe(3, K_))
            term1 = I1 * (1.0 - 3.0 * Cv / K_**2 - 3.0 * Sv**2 / K_**4)
            term2 = 1.5 * (I0 + I2) * (Cv / K_ + Sv**2 / K_**3)
            term3 = -((3.0 * I1 + I3) / 4.0) * (Sv**2 / K_**2)
            return (term1 + term2 + term3) / K_

        Cm = kh * kg * np.cos(d_minus)
        Cp = kh * kg * np.cos(d_plus)
        Sm = kh * kg * np.sin(d_minus)
        Sp = kh * kg * np.sin(d_plus)

        Uppp = (
            8.0 * pi * S * A * kh * kg
            * (np.sin(d_minus) * D(Km, Cm, Sm) + np.sin(d_plus) * D(Kp, Cp, Sp))
        )
        c3 = Uppp / 6.0
        return float(c1), float(c3)
