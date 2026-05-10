"""Figure 1: pitchfork bifurcation in (kappa, delta) plane."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless

from steering.visualization.figures import figure_1_bifurcation


def main(out: Path = Path("figures/fig1_bifurcation.png")) -> Path:
    fig, _ = figure_1_bifurcation(
        kappa_range=(0.5, 6.0),
        delta_range=(0.0, 1.55),
        n_kappa=60,
        n_delta=60,
        inset_kappa=2.0,
        inset_deltas=(0.5, 1.45),
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"wrote {out}")
    return out


if __name__ == "__main__":
    main()
