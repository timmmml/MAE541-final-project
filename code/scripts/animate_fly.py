"""Render mp4 animations of the fly cartoon tracking simulated heading θ(t).

Scenario: κ_g = κ_h = 1.0, δ = 1.49, DiscretePFLModel(N_default=12),
AccelerationDynamics with γ = 0.1, additive forcing F cos(ω t), ω = sqrt(c1).
One mp4 is written per forcing amplitude.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from PIL import Image

from steering.params import ModelParams, ForcingParams
from steering.models import BesselSteeringModel, DiscretePFLModel
from steering.dynamics import AccelerationDynamics
from steering.integrator import Simulation


REPO_ROOT = Path(__file__).resolve().parents[1]
FLY_PATH = REPO_ROOT / "cartoon" / "fly.png"


def simulate(F: float, T: float, fps: int, gamma: float = 0.1):
    p = ModelParams(kappa_h=1.0, kappa_g=1.0, delta=1.49)
    bessel = BesselSteeringModel()
    c1, _ = bessel.taylor_coefficients(p)
    omega = float(np.sqrt(c1))

    neural = DiscretePFLModel(N_default=12)
    dyn = AccelerationDynamics(model=neural, gamma=gamma)
    forc = ForcingParams(F=F, omega=omega)
    sim = Simulation(dyn, p, forc, rtol=1e-9, atol=1e-11)

    t_eval = np.linspace(0.0, T, int(T * fps) + 1)
    res = sim.run(np.array([0.01, 0.0]), (0.0, T), t_eval=t_eval)
    return res, omega


def render(
    res,
    omega: float,
    F: float,
    fps: int,
    out: Path,
    skip: float,
    box: bool,
    annotate: bool,
) -> None:
    fly = Image.open(FLY_PATH).convert("RGBA")
    # Pad to a transparent square with side = ceil(diagonal) so rotation
    # at any angle never clips the fly out of the canvas.
    diag = int(np.ceil(np.hypot(*fly.size)))
    canvas = Image.new("RGBA", (diag, diag), (0, 0, 0, 0))
    canvas.paste(fly, ((diag - fly.width) // 2, (diag - fly.height) // 2), fly)
    fly = canvas

    start_idx = int(np.searchsorted(res.t, skip))
    t = res.t[start_idx:]
    theta = res.theta[start_idx:]
    if len(t) == 0:
        raise ValueError(f"skip={skip}s leaves no frames (T must exceed skip).")

    fig, ax = plt.subplots(figsize=(5, 5))
    if annotate:
        fig.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.05)
    else:
        fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
    im_artist = ax.imshow(fly)
    ax.set_aspect("equal")
    if box:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.set_axis_off()

    title = ax.set_title("") if annotate else None

    def update(i: int):
        rot = fly.rotate(np.degrees(theta[i]), resample=Image.BICUBIC, expand=False)
        im_artist.set_data(rot)
        artists = [im_artist]
        if title is not None:
            title.set_text(
                f"F = {F:.4f}   ω = {omega:.3f}   "
                f"t = {t[i]:6.2f} s   θ = {theta[i]:+.2f} rad"
            )
            artists.append(title)
        return tuple(artists)

    ani = FuncAnimation(fig, update, frames=len(t), interval=1000 / fps, blit=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(out), writer=FFMpegWriter(fps=fps, bitrate=2400))
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--F", type=float, nargs="+", default=[0.018],
        help="Forcing amplitude(s). One mp4 is written per value.",
    )
    ap.add_argument(
        "--T", type=float, default=60.0,
        help="Simulated/real-time duration in seconds.",
    )
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--gamma", type=float, default=0.1)
    ap.add_argument(
        "--skip", type=float, default=0.0,
        help="Discard this many seconds from the start of the trajectory "
             "(transient). Annotated time begins at this offset.",
    )
    ap.add_argument(
        "--box", dest="box", action="store_true", default=True,
        help="Draw the axes frame around the fly (default).",
    )
    ap.add_argument(
        "--no-box", dest="box", action="store_false",
        help="Hide the axes frame.",
    )
    ap.add_argument(
        "--annotate", dest="annotate", action="store_true", default=True,
        help="Show the F / ω / t / θ title overlay (default).",
    )
    ap.add_argument(
        "--no-annotate", dest="annotate", action="store_false",
        help="Hide the title overlay.",
    )
    ap.add_argument(
        "--out-dir", type=Path, default=REPO_ROOT / "figures" / "fly_anim",
    )
    args = ap.parse_args()

    for F in args.F:
        res, omega = simulate(F=F, T=args.T, fps=args.fps, gamma=args.gamma)
        out = args.out_dir / f"fly_F{F:.4f}_T{args.T:g}s.mp4"
        print(f"[{F=}] rendering -> {out}")
        render(
            res, omega=omega, F=F, fps=args.fps, out=out,
            skip=args.skip, box=args.box, annotate=args.annotate,
        )
        print(f"  done: {out}")


if __name__ == "__main__":
    main()
