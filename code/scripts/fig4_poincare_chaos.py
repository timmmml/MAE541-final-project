"""Figure 4: Poincaré bifurcation diagram + sample sections."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from steering.params import ModelParams
from steering.visualization.figures import figure_4_poincare


def main(out: Path = Path("figures/fig4_poincare_chaos.png"), fast: bool = False) -> Path:
    params = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
    fig, _ = figure_4_poincare(
        base_params=params,
        gamma=0.15,
        omega=1.0,
        F_range=(0.05, 1.6),
        n_F=40 if fast else 60,
        F_section_examples=(0.1, 0.6, 1.2),
        n_periods_transient=60 if fast else 120,
        n_periods_record=40 if fast else 60,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"wrote {out}")
    return out


if __name__ == "__main__":
    import sys
    main(fast="--fast" in sys.argv)
