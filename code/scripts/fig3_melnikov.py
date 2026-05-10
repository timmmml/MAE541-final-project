"""Figure 3: Melnikov F_crit curves vs omega and vs delta."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from steering.params import ModelParams
from steering.visualization.figures import figure_3_melnikov


def main(out: Path = Path("figures/fig3_melnikov.png")) -> Path:
    params = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
    fig, _ = figure_3_melnikov(
        base_params=params,
        gamma=0.05,
        omega_range=(0.2, 4.0),
        n_omega=40,
        delta_range=(1.33, 1.55),
        n_delta=14,
        omega_fixed=1.0,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"wrote {out}")
    return out


if __name__ == "__main__":
    main()
