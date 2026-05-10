"""Figure 2: phase portraits Duffing vs Bessel, undamped vs damped."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from steering.params import ModelParams
from steering.visualization.figures import figure_2_phase_portraits


def main(out: Path = Path("figures/fig2_phase_portraits.png")) -> Path:
    params = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
    fig, _ = figure_2_phase_portraits(params=params, gamma=0.1)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"wrote {out}")
    return out


if __name__ == "__main__":
    main()
