"""Poincaré section for a grid of initial conditions, colored by IC.

Strobe at the forcing period T = 2π/ω. For each initial state on a grid in
(θ, v) ∈ [-π, π] × [-1, 1], simulate ``T_skip + T`` forcing periods, drop the
first ``T_skip`` for transient, and scatter the remaining stroboscopic samples.
Each point is colored by its initial condition via a 2-D HSV map
(hue ← θ₀, value ← v₀).
"""

from __future__ import annotations

import argparse
import colorsys
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from steering.params import ModelParams, ForcingParams
from steering.models import (
    BesselSteeringModel,
    DiscretePFLModel,
    DuffingModel,
    ContinuousPFLModel,
)
from steering.dynamics import AccelerationDynamics
from steering.integrator import Simulation
from steering.analysis.poincare import stroboscopic_section


REPO_ROOT = Path(__file__).resolve().parents[1]


def ic_color(theta0: float, v0: float, v_lim: float) -> tuple[float, float, float]:
    """Map (θ₀, v₀) to an RGB color: hue from θ₀, value from v₀."""
    h = (theta0 + np.pi) / (2.0 * np.pi)
    val = 0.35 + 0.65 * (v0 + v_lim) / (2.0 * v_lim)
    return colorsys.hsv_to_rgb(h % 1.0, 0.85, float(np.clip(val, 0.0, 1.0)))


# Worker-side state. Each worker process builds its own Simulation once,
# so per-IC work is just `.run(...)` and `stroboscopic_section(...)`.
_WORKER: dict = {}


def _build_sim(cfg: dict):
    p = ModelParams(kappa_h=cfg["kappa_h"], kappa_g=cfg["kappa_g"], delta=cfg["delta"])
    name = cfg["model"]
    if name == "duffing":
        bessel = BesselSteeringModel()
        c1, c3 = bessel.taylor_coefficients(p)
        model = DuffingModel(c1, c3)
    elif name == "discrete_pfl":
        model = DiscretePFLModel(N_default=cfg["N"])
    elif name == "continuous_pfl":
        model = ContinuousPFLModel()
    elif name == "bessel":
        model = BesselSteeringModel()
    else:
        raise ValueError(f"unknown model {name!r}")
    dyn = AccelerationDynamics(model=model, gamma=cfg["gamma"])
    forc = ForcingParams(F=cfg["F"], omega=cfg["omega"])
    return Simulation(dyn, p, forc, rtol=1e-9, atol=1e-11)


def _init_worker(cfg: dict) -> None:
    _WORKER["sim"] = _build_sim(cfg)
    _WORKER["omega"] = cfg["omega"]
    _WORKER["T_skip"] = cfg["T_skip"]
    _WORKER["t_total"] = cfg["t_total"]


def _simulate_ic(item):
    idx, theta0, v0 = item
    sim = _WORKER["sim"]
    res = sim.run(
        np.array([theta0, v0]), (0.0, _WORKER["t_total"]), dense_output=True
    )
    sec = stroboscopic_section(res, _WORKER["omega"], transient_periods=_WORKER["T_skip"])
    return idx, sec


def run(args) -> None:
    p = ModelParams(kappa_h=args.kappa_h, kappa_g=args.kappa_g, delta=args.delta)

    omega = args.omega
    if omega is None:
        bessel = BesselSteeringModel()
        c1, _ = bessel.taylor_coefficients(p)
        omega = float(np.sqrt(c1))
        print(f"[omega] auto: sqrt(c1) = {omega:.4f}")

    T_period = 2.0 * np.pi / omega
    t_total = (args.T_skip + args.T) * T_period

    cfg = dict(
        kappa_h=args.kappa_h, kappa_g=args.kappa_g, delta=args.delta,
        gamma=args.gamma, F=args.F, omega=omega,
        model=args.model, N=args.N,
        T_skip=args.T_skip, t_total=t_total,
    )

    thetas0 = np.linspace(-np.pi, np.pi, args.n_theta, endpoint=False)
    vs0 = np.linspace(-args.v_lim, args.v_lim, args.n_v)

    items = [
        (i * args.n_v + j, float(th), float(v))
        for i, th in enumerate(thetas0) for j, v in enumerate(vs0)
    ]
    n_ic = len(items)

    n_workers = args.n_workers or max(1, (os.cpu_count() or 2) - 1)
    n_workers = min(n_workers, n_ic)
    print(
        f"[grid] {args.n_theta}×{args.n_v} = {n_ic} ICs, "
        f"T_skip={args.T_skip} periods, T={args.T} periods, "
        f"ω={omega:.4f}, F={args.F:.4f}, γ={args.gamma:.4f}, "
        f"workers={n_workers}"
    )

    results: dict[int, np.ndarray] = {}
    with Pool(processes=n_workers, initializer=_init_worker, initargs=(cfg,)) as pool:
        for k, (idx, sec) in enumerate(
            pool.imap_unordered(_simulate_ic, items, chunksize=1), start=1
        ):
            results[idx] = sec
            if k % max(1, n_ic // 20) == 0 or k == n_ic:
                print(f"  {k}/{n_ic} ICs done")

    fig, ax = plt.subplots(figsize=(7, 6))
    for idx, th0, v0 in items:
        sec = results.get(idx)
        if sec is None or sec.size == 0:
            continue
        # Wrap θ to [-π, 3π) so the section displays two copies of the
        # fundamental domain side by side.
        theta_wrapped = (sec[:, 0] + np.pi) % (4.0 * np.pi) - np.pi
        c = ic_color(th0, v0, args.v_lim)
        ax.scatter(
            theta_wrapped, sec[:, 1],
            color=c, s=args.markersize, alpha=args.alpha, linewidths=0,
        )

    ax.set_xlim(-np.pi, 3.0 * np.pi)
    ax.set_ylim(-args.v_lim, args.v_lim)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$v$")
    ax.set_title(
        rf"Poincaré grid: $\kappa_h={args.kappa_h}, \kappa_g={args.kappa_g}, "
        rf"\delta={args.delta},\ \gamma={args.gamma},\ "
        rf"F={args.F},\ \omega={omega:.3f}$"
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi)
    print(f"[saved] {args.out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kappa_h", type=float, default=1.0)
    ap.add_argument("--kappa_g", type=float, default=1.0)
    ap.add_argument("--delta", type=float, default=1.49)
    ap.add_argument("--gamma", type=float, default=0.1)
    ap.add_argument("--F", type=float, default=0.018)
    ap.add_argument(
        "--omega", type=float, default=None,
        help="Forcing angular frequency. Default: sqrt(c1) from Bessel Taylor expansion.",
    )
    ap.add_argument(
        "--T", type=int, default=400,
        help="Number of forcing periods to record (after T_skip).",
    )
    ap.add_argument(
        "--T_skip", type=int, default=100,
        help="Number of forcing periods to discard as transient.",
    )
    ap.add_argument("--n_theta", type=int, default=12)
    ap.add_argument("--n_v", type=int, default=7)
    ap.add_argument("--v_lim", type=float, default=1.0)
    ap.add_argument(
        "--model",
        choices=["discrete_pfl", "continuous_pfl", "bessel", "duffing"],
        default="discrete_pfl",
    )
    ap.add_argument("--N", type=int, default=12, help="N for discrete PFL model.")
    ap.add_argument("--markersize", type=float, default=1.5)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument(
        "--n_workers", type=int, default=0,
        help="Number of parallel worker processes. 0 = cpu_count - 1.",
    )
    ap.add_argument(
        "--out", type=Path,
        default=REPO_ROOT / "figures" / "poincare_grid.png",
    )
    args = ap.parse_args()
    run(args)


if __name__ == "__main__":
    main()
