# drosophila-steering

Reference implementation accompanying the MAE541 final project on dual-goal
steering dynamics in the *Drosophila* central complex. Provides a hierarchy of
models for the PFL3 steering drive (cubic Duffing normal form → exact Bessel
closed-form → continuous / discrete population integration → full circuit with
PFL2 indirect pathway), along with dynamics, integration, and analysis tools
(fixed-point and bifurcation analysis, Hamiltonian / homoclinic computation,
Melnikov function, Poincaré sections, Lyapunov exponents).

See `../SPEC.md` for the full design specification.

## Install

```
cd code
pip install -e .          # core
pip install -e .[dev]     # core + tests + jupyter + xarray
```

## Quick start

```python
import numpy as np
from steering.params import ModelParams
from steering.models import DuffingModel, BesselSteeringModel
from steering.dynamics import AccelerationDynamics
from steering.integrator import Simulation

params = ModelParams(kappa_h=3.0, kappa_g=3.0, delta=np.pi/4)
bessel = BesselSteeringModel()
duffing = DuffingModel.from_params(params)

dyn = AccelerationDynamics(model=bessel, gamma=0.1)
sim = Simulation(dyn, params)
result = sim.run(state0=np.array([0.1, 0.0]), t_span=(0.0, 200.0))
```

## Layout

```
src/steering/
  params.py            ModelParams, ForcingParams
  nonlinearities.py    quadratic / elu / relu / softplus registry
  utils/bessel.py      overflow-safe I_n
  models/              Steering drive U(theta) at five abstraction levels
  dynamics/            Velocity / acceleration RHS + topology
  integrator.py        Simulation wrapper around solve_ivp
  sweep.py             ParameterSweep (multiprocessing)
  analysis/            fixed_points, bifurcation, hamiltonian, homoclinic,
                       melnikov, poincare, lyapunov
  visualization/       style, profiles, phase_portrait, bifurcation_plot,
                       poincare_plot
tests/                 Consistency / symmetry / limits / dynamics / Melnikov
```

## Tests

```
cd code
pytest
```
