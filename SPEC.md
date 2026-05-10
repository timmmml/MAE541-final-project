# Drosophila Dual-Goal Steering Dynamics — Project Specification 

## Overview

This codebase investigates how two simultaneous goal representations in the *Drosophila* central complex steering circuit generate Duffing-like dynamics, pitchfork bifurcations, and chaos under periodic forcing. The system is studied at multiple levels of biological detail, from a cubic (Duffing) normal form through the exact Bessel closed-form steering drive to a full discretized neural circuit with PFL2 indirect-pathway modulation.

---

## 1. Mathematical Background

### 1.1 The Steering Drive

The steering signal is the difference in population activity between PFL3L and PFL3R subpopulations. Under a quadratic nonlinearity $f(x) = x^2$ and continuous population coverage, this has the closed form:

$$
U(\theta, \delta) = 4\pi S A \bigl[
  I_0(\kappa(\theta + \Delta, \delta)) + I_0(\kappa(\theta + \Delta, -\delta))
  - I_0(\kappa(\theta - \Delta, \delta)) - I_0(\kappa(\theta - \Delta, -\delta))
\bigr]
$$

where $\Delta \equiv \Delta_{P3L} = 67.5° = 3\pi/8$ is the PFL3 population phase shift, $\delta$ is the half-separation between goals, $S$ is synaptic scaling, $A$ is goal amplitude, and:

$$
\kappa(a, b) = \sqrt{\kappa_h^2 + \kappa_g^2 + 2\kappa_h \kappa_g \cos(a - b)}
$$

with $\kappa_h, \kappa_g$ the heading and goal bump concentrations respectively. Under equal concentrations $\kappa_h = \kappa_g = \kappa$:

$$
\kappa(a, b) = 2\kappa \left|\cos\!\left(\frac{a - b}{2}\right)\right|
$$

### 1.2 Taylor Coefficients

Expanding $U(\theta, \delta)$ around $\theta = 0$:

$$
U(\theta, \delta) = c_1(\delta, \kappa, \Delta)\,\theta + c_3(\delta, \kappa, \Delta)\,\theta^3 + O(\theta^5)
$$

where $c_0 = 0$ (odd symmetry $U(-\theta) = -U(\theta)$), $c_2 = 0$ (same), and $c_1, c_3$ are given here, involving $I_0, I_1, I_2, I_3$ evaluated at the combined concentrations.

$$ \boxed{U'(0, \delta) = -8\pi SA, \kappa_h\kappa_g \left[ \frac{I_1(\kappa(\Delta_{P3L}, \delta))}{\kappa(\Delta_{P3L}, \delta)} \sin(\Delta_{P3L} - \delta) + \frac{I_1(\kappa(\Delta_{P3L}, -\delta))}{\kappa(\Delta_{P3L}, -\delta)} \sin(\Delta_{P3L} + \delta) \right]} $$

the third order coefficient is more complicated
$$
\boxed{U'''(0, \delta) = 8\pi SA\kappa_{h}\kappa_{g}[\sin(\Delta_{P3L}-\delta)\cdot D(\kappa_{-},C_{-},S_{-})+\sin (\Delta_{P_{3}L}+\delta)\cdot D(\kappa_{+},C_{+},S_{+})]}
$$
with
$$
D(K,C,S)=\frac{1}{K}​\left[ I_{1}​(K)\left( 1-\frac{3C}{K^{2}}​−\frac{3S^{2}}{K^{4}} \right)+\frac{3(I_{0}(K)+I_{2}​(K))}{2}​\left( \frac{C}{K} + \frac{S^{2}}{K^{3}}​ \right)−\frac{3I_{1}(K)+I_{3}(K)}{4}\cdot \frac{S^{2}}{K^{2}} \right]
$$
and
$$
\begin{align}
\kappa_- &\equiv \kappa(\Delta_{P3L}, \delta) = \sqrt{\kappa_h^2 + \kappa_g^2 + 2\kappa_h\kappa_g\cos(\Delta_{P3L} - \delta)} \\
\kappa_+ &\equiv \kappa(\Delta_{P3L}, -\delta) = \sqrt{\kappa_h^2 + \kappa_g^2 + 2\kappa_h\kappa_g\cos(\Delta_{P3L} + \delta)} \\
C_- &\equiv \kappa_h\kappa_g\cos(\Delta_{P3L} - \delta), \qquad C_+ \equiv \kappa_h\kappa_g\cos(\Delta_{P3L} + \delta) \\
S_- &\equiv \kappa_h\kappa_g\sin(\Delta_{P3L} - \delta), \qquad S_+ \equiv \kappa_h\kappa_g\sin(\Delta_{P3L} + \delta)
\end{align}
$$


### 1.3 Dynamical Systems

**Velocity control** (1D gradient):
$$
\dot{\theta} = U(\theta, \delta)
$$

**Acceleration control** (Duffing-like, 2D):
$$
\ddot{\theta} + \gamma \dot{\theta} = U(\theta, \delta)
$$

**Forced acceleration control**:
$$
\ddot{\theta} + \gamma \dot{\theta} = U(\theta, \delta) + F \cos(\omega t)
$$

Each can be studied on $\mathbb{R}$ (or $\mathbb{R}^2$) for local analysis, or on $S^1$ (or $S^1 \times \mathbb{R}$) for global analysis.

### 1.4 Homoclinic Orbit (Reduced System, $\mathbb{R}^2$)

For the Duffing truncation $U \approx c_1 \theta + c_3 \theta^3$ with $c_1 > 0$ and $c_3 < 0$, the unforced undamped system is Hamiltonian with:

$$
H = \frac{1}{2}v^2 - \frac{1}{2}c_1 \theta^2 - \frac{1}{4}c_3 \theta^4
$$

Homoclinic orbit to the origin:


$$
\theta_0(t) = \pm\sqrt{\frac{2c_1}{|c_3|}} \operatorname{sech}(\sqrt{c_1}\, t), \qquad
v_0(t) = \dot{\theta}_0(t)
$$

### 1.5 Melnikov Function

For the forced damped system with $\gamma, F = O(\epsilon)$:

$$
M(t_0) = -\frac{4\gamma\, c_1^{3/2}}{3|c_3|}
+ \frac{F \pi \omega \sqrt{2c_1/|c_3|}}{\cosh\!\left(\frac{\pi\omega}{2\sqrt{c_1}}\right)} \sin(\omega t_0)
$$

Simple zeros (and Smale horseshoe) when:

$$
F > F_{\mathrm{crit}}(\omega) = \frac{4\gamma\, c_1}{3\pi\omega\sqrt{2|c_3|}}\,\cosh\!\left(\frac{\pi\omega}{2\sqrt{c_1}}\right)
$$

---

## 2. Model Hierarchy

All models expose the same interface: given state and parameters, return the steering drive $U(\theta)$ (and optionally its derivatives). This lets all analysis and integration tools work against any model interchangeably.

### 2.1 Abstract Base: `SteeringModel`

```
class SteeringModel(ABC):
    """
    Base class for all steering models.
    All models compute U(θ, params) — the net steering drive.
    """
    @abstractmethod
    def steering_drive(self, theta: float | ndarray, params: ModelParams) -> float | ndarray:
        """Return U(θ)."""

    def steering_drive_derivatives(self, theta, params, order=3):
        """Return (U, U', U'', U''') at θ. Default: finite differences. Subclasses may override analytically."""

    def steering_potential(self, theta, params):
        """Return V(θ) = -∫U dθ via numerical quadrature. Subclasses may override."""
```

`ModelParams` is a frozen dataclass:

```
@dataclass(frozen=True)
class ModelParams:
    kappa_h: float          # heading bump concentration
    kappa_g: float          # goal bump concentration
    delta: float            # half-separation between goals (radians)
    Delta_pop: float = 3*pi/8   # PFL3 population phase shift (67.5°)
    S: float = 1.0          # synaptic scaling
    A: float = 1.0          # goal amplitude
    W_D3: float = 0.0       # indirect pathway weight (PFL2→DNa03→DNa02)
    N_neurons: int | None = None   # None = continuous limit
    nonlinearity: str = "quadratic"  # "quadratic" | "elu" | "relu" | "softplus" | callable
    nonlinearity_params: dict = field(default_factory=dict)  # e.g. {"alpha": 1.0} for ELU

    @property
    def kappa_equal(self) -> float | None:
        """Return κ if κ_h == κ_g, else None."""
        if np.isclose(self.kappa_h, self.kappa_g):
            return self.kappa_h
        return None
```

### 2.2 Level 0: `DuffingModel`

Cubic truncation: $U(\theta) = c_1 \theta + c_3 \theta^3$.

- Stores $c_1, c_3$ either directly or computed from `ModelParams` using the Bessel coefficient formulas (eqs. 14, 16).
- Analytical derivatives trivially available.
- Method `from_params(params) -> DuffingModel` computes $c_1, c_3$ from $(\kappa, \delta, \Delta, S, A)$.
- Method `fixed_points() -> list[float]` returns $\{0\}$ or $\{0, \pm\sqrt{-c_1/c_3}\}$.
- Method `homoclinic_orbit(t_array) -> (theta, v)` returns the sech solution.
- Method `melnikov_threshold(gamma, omega) -> F_crit`.
- Method `hamiltonian(theta, v) -> H`.

### 2.3 Level 1: `BesselSteeringModel`

Exact closed-form from eq. (9): evaluates $U(\theta)$ using $I_0$ at the four combined concentrations. This is exact for quadratic nonlinearity and continuous population, with no Taylor truncation.

- Uses `scipy.special.i0` (and `i1`, `iv` for derivatives).
- Analytical first derivative $U'(\theta)$ available via eq. (13); higher derivatives via the recurrence in eq. (15). Implement these analytically, with finite-difference fallback for validation.
- Method `taylor_coefficients(order=3) -> (c1, c3)` computes the Taylor coefficients at $\theta = 0$ and returns a `DuffingModel`.
- This is the primary "full-but-tractable" model for comparing against the Duffing reduction.

### 2.4 Level 2: `ContinuousPFLModel`

Numerically integrates the population activity (eq. 4) with an arbitrary nonlinearity $f$:

$$
G_{P3L}(\theta, \delta) = \int_{-\pi}^{\pi} f\bigl(h_{P3L}(\theta, \phi, \delta)\bigr) \, d\phi
$$

- Uses `scipy.integrate.quad` (or fixed-point Gauss-Legendre quadrature for speed).
- The nonlinearity $f$ is selected from a registry or passed as a callable.
- The PFL3L/R drives are computed from the von Mises heading and goal profiles directly: 
	Taking a reference frame where the two goals, $\theta_{\text{goal}}^{(0)}$ and $\theta_{\text{goal}}^{(1)}$, are located at $\delta$ and $-\delta$ respectively, we can write the von Mises bump representation, under this dual-goal, for a goal neuron with preference $\phi$ as:
	
	$$ g_{\text{goal}}(\phi, \delta) = \exp(\kappa_{\text{goal}} \cos(\delta - \phi)) + \exp(\kappa_{\text{goal}} \cos(-\delta - \phi)) \tag{1} $$
	
	similarly write down the heading neurons with preference $\phi$ firing to the current heading $\theta$ as:
	
	$$ g_{\text{heading}}(\phi, \theta) = \exp(\kappa_{\text{heading}} \cos(\theta - \phi)) \tag{2} $$
	
	Considering a PFL cell with goal preference $\phi$, we abstract that each cell receives a direct synaptic connection from the associated goal neuron, as well as a synaptic connection from a heading neuron tuned to $\Delta_{\text{pop}}$ away from $\phi$. Showing the PFL3L subpopulation as an example, we can write down the total drive to a PFL3L neuron with preference $\phi$, $h_{P3L}$, and the neural activity after the nonlinearity $f$, $g_{P3L}$, as:
	
	$$ \begin{aligned} h_{P3L}(\theta, \phi, \delta) &= S \cdot g_{\text{heading}}(\phi - \Delta_{P3L}, \theta) + A \cdot g_{\text{goal}}(\phi, \delta) \ &= S \cdot \exp(\kappa_{\text{heading}} \cos(\theta - (\phi - \Delta_{P3L}))) \ &\quad + A \cdot (\exp(\kappa_{\text{goal}} \cos(\delta - \phi)) + \exp(\kappa_{\text{goal}} \cos(-\delta - \phi))) \ g_{P3L}(\theta, \phi, \delta) &= f(h_{P3L}(\theta, \phi, \delta)) \end{aligned} \tag{3} $$
	
	As downstream neurons to the PFLs (DNa02 and DNa03) effectively sum inputs from all neurons of each population (PFL3L/R, PFL2), it behooves us to consider the population activity as a whole, integrating eq. (3) over all $\phi$. We can write down the population activity of the PFL3L subpopulation:
	
	$$ \begin{aligned} G_{P3L}(\theta, \delta) &= \int_{-\pi}^{\pi} g_{P3L}(\theta, \phi, \delta), d\phi \ &= \int_{-\pi}^{\pi} f(h_{P3L}(\theta, \phi, \delta)), d\phi \end{aligned} \tag{4} $$
	
	We choose the quadratic nonlinearity $f(x) = x^2$, as a common choice in neuroscience and one that allows conjunctive tuning for heading and goal states on the level of PFL3 neurons. This is important to allow differences between the L and R populations after integration.
	
	Hence we have the following:
	
	$$ \begin{aligned} G_{P3L}(\theta, \delta) &= \int_{-\pi}^{\pi} h_{P3L}(\theta, \phi, \delta)^2 , d\phi \ &= \int_{-\pi}^{\pi} (S \cdot g_{\text{heading}}(\phi - \Delta_{P3L}, \theta) + A \cdot g_{\text{goal}}(\phi, \delta))^2 , d\phi \ &= \int_{-\pi}^{\pi} \big[ S^2 \cdot g_{\text{heading}}(\phi - \Delta_{P3L}, \theta)^2 + A^2 \cdot g_{\text{goal}}(\phi, \delta)^2 \ &\quad + 2SA \cdot g_{\text{heading}}(\phi - \Delta_{P3L}, \theta) \cdot g_{\text{goal}}(\phi, \delta) \big] d\phi \end{aligned} \tag{5} $$
	

- Agrees with `BesselSteeringModel` when $f(x) = x^2$ (use this as a consistency check).

**Nonlinearity registry** (module `steering.nonlinearities`):
```
def quadratic(x): return x**2
def elu(x, alpha=1.0): return np.where(x >= 0, x, alpha*(np.exp(x) - 1))
def relu(x): return np.maximum(0, x)
def softplus(x, beta=1.0): return np.log1p(np.exp(beta * x)) / beta
```

### 2.5 Level 3: `DiscretePFLModel`

Replaces the integral with a sum over $N$ neurons at uniformly-spaced preferred headings $\phi_j = -\pi + 2\pi j / N$ for $j = 0, \ldots, N-1$:

$$
G_{P3L}(\theta, \delta) = \frac{2\pi}{N} \sum_{j=0}^{N-1} f\bigl(h_{P3L}(\theta, \phi_j, \delta)\bigr)
$$

- Inherits from `ContinuousPFLModel`, overriding the integration with a sum.
- $N$ is set via `ModelParams.N_neurons`. The biological value for *Drosophila* is $N = 8$ (8 glomeruli per hemisphere in the protocerebral bridge, though this should be a parameter).
- Recovers `ContinuousPFLModel` in the limit $N \to \infty$.

### 2.6 Level 4: `FullCircuitModel`

Adds the PFL2 population and the indirect pathway (eq. 6 of the proposal):

$$
U(\theta) = \bigl[\Sigma_{P3R} - \Sigma_{P3L}\bigr]
+ W_{D3} \bigl[f(\Sigma_{P3R} + \Sigma_{P2}) - f(\Sigma_{P3L} + \Sigma_{P2})\bigr]
$$

where $\Sigma_{P2}$ is the PFL2 population sum with its own phase shift ($\Delta_{P2} = \pi$, i.e., 180°).

- Can operate in continuous or discrete mode (controlled by `ModelParams.N_neurons`).
- When $W_{D3} = 0$, reduces to `DiscretePFLModel` or `ContinuousPFLModel`.
- The PFL2 pathway uses the same nonlinearity as PFL3 by default, but can use a different one (add `nonlinearity_pfl2` field to `ModelParams` if needed).

---

## 3. Dynamics Layer

The dynamics layer wraps a `SteeringModel` with a control law, topology, and (optionally) forcing. It produces the right-hand side for ODE integration.

### 3.1 Abstract Base: `Dynamics`

```
class Dynamics(ABC):
    model: SteeringModel
    topology: Literal["planar", "cylindrical"]  # R^n vs S^1 × R^(n-1)

    @abstractmethod
    def rhs(self, t: float, state: ndarray, params: ModelParams, forcing_params: ForcingParams | None) -> ndarray:
        """Return dstate/dt."""

    @abstractmethod
    def state_dim(self) -> int:
        """Dimension of state space."""

    def wrap_state(self, state: ndarray) -> ndarray:
        """Apply topology. For cylindrical: wrap θ to [-π, π)."""
```

```
@dataclass(frozen=True)
class ForcingParams:
    F: float = 0.0          # forcing amplitude
    omega: float = 1.0      # forcing frequency
    forcing_type: str = "additive"  # "additive" | "parametric_delta"
    # For parametric: δ(t) = δ_0 + F·cos(ωt), so U depends on time through δ
```

### 3.2 `VelocityDynamics`

State: $\theta$ (1D).

$$
\dot{\theta} = U(\theta, \delta)
$$

No forcing in the standard sense (1D autonomous systems can't be chaotic). But parametric forcing through $\delta(t)$ makes it a 1.5D nonautonomous system — implementable but not the primary focus.

### 3.3 `AccelerationDynamics`

State: $(\theta, v)$ (2D).

$$
\dot{\theta} = v, \qquad \dot{v} = -\gamma v + U(\theta, \delta) + F\cos(\omega t)
$$

Parameters: $\gamma$ (damping), plus `ForcingParams`.

For parametric forcing through $\delta$: $U$ is evaluated at $\delta(t) = \delta_0 + F\cos(\omega t)$.

---

## 4. Integration and Simulation

### 4.1 Module: `steering.integrator`

Wraps `scipy.integrate.solve_ivp` with sensible defaults for stiff and non-stiff systems.

```
class Simulation:
    """
    Integrates a Dynamics object over a time span.
    """
    def __init__(self, dynamics: Dynamics, params: ModelParams,
                 forcing: ForcingParams | None = None,
                 gamma: float = 0.0,
                 method: str = "RK45",
                 rtol: float = 1e-10, atol: float = 1e-12):
        ...

    def run(self, state0: ndarray, t_span: tuple[float, float],
            t_eval: ndarray | None = None,
            events: list[callable] | None = None) -> SimulationResult:
        ...

@dataclass
class SimulationResult:
    t: ndarray
    states: ndarray            # shape (n_times, state_dim)
    model: SteeringModel
    params: ModelParams
    forcing: ForcingParams | None
```

### 4.2 Module: `steering.sweep`

Parameter sweep infrastructure for batch simulations:

```
class ParameterSweep:
    """
    Runs simulations over a grid of parameter values.
    Supports parallelism via multiprocessing.
    """
    def __init__(self, dynamics: Dynamics, base_params: ModelParams,
                 sweep_axes: dict[str, ndarray],
                 # e.g. {"delta": linspace(0, pi/2, 50), "kappa_h": linspace(0.5, 5, 50)}
                 analysis_fn: Callable[[SimulationResult], Any]):
        ...

    def run(self, state0, t_span, n_workers=None) -> xarray.DataArray:
        ...
```

---

## 5. Analysis Tools

### 5.1 Module: `steering.analysis.fixed_points`

```
def find_fixed_points_1d(model: SteeringModel, params: ModelParams,
                          theta_range: tuple = (-pi, pi),
                          n_initial: int = 200) -> list[FixedPoint]:
    """
    Find zeros of U(θ) on the given range via rootfinding from a grid of initial guesses.
    Returns FixedPoint objects with location, stability (sign of U'), and index.
    """

def find_fixed_points_2d(dynamics: AccelerationDynamics, params: ModelParams,
                          gamma: float,
                          theta_range: tuple = (-pi, pi)) -> list[FixedPoint2D]:
    """
    Find fixed points of (v=0, -γv + U(θ) = 0), i.e., same θ* as 1D but with
    2D stability classification (node/spiral/saddle) from eigenvalues of the Jacobian.
    """

@dataclass
class FixedPoint:
    theta: float
    U_prime: float       # dU/dθ at fixed point
    stability: str       # "stable" | "unstable"

@dataclass
class FixedPoint2D:
    theta: float
    eigenvalues: tuple[complex, complex]
    classification: str  # "stable_node" | "stable_spiral" | "unstable_node" | "unstable_spiral" | "saddle" | "center"
```

### 5.2 Module: `steering.analysis.bifurcation`

```
def pitchfork_bifurcation_curve(model: SteeringModel,
                                  param1_name: str, param1_range: ndarray,
                                  param2_name: str, param2_range: ndarray,
                                  base_params: ModelParams) -> BifurcationDiagram:
    """
    Compute the bifurcation curve c1(param1, param2) = 0 in a 2D parameter plane.
    Also evaluates sgn(c3) to classify supercritical vs subcritical.

    For BesselSteeringModel: uses analytical Taylor coefficients.
    For other models: uses finite-difference derivatives of U at θ=0.

    Returns:
        BifurcationDiagram with:
        - c1_grid: 2D array of c1 values
        - c3_grid: 2D array of c3 values
        - bifurcation_curve: contour at c1 = 0
        - type_along_curve: "supercritical" (c3 < 0) or "subcritical" (c3 > 0)
    """

def numerical_bifurcation_diagram_1d(model: SteeringModel,
                                      param_name: str, param_range: ndarray,
                                      base_params: ModelParams,
                                      theta_range=(-pi, pi)) -> dict:
    """
    For each parameter value, find all fixed points and their stability.
    Returns arrays suitable for plotting stable/unstable branches.
    Used for overlaying on the analytical pitchfork locus.
    """
```

### 5.3 Module: `steering.analysis.hamiltonian`

```
def hamiltonian_from_model(model: SteeringModel, params: ModelParams,
                            theta: ndarray, v: ndarray) -> ndarray:
    """
    Compute H(θ, v) = ½v² + V(θ) where V(θ) = -∫₀^θ U(θ') dθ'.
    Works for any model (numerical quadrature for V).
    For DuffingModel, uses exact V(θ) = -½c₁θ² - ¼c₃θ⁴.
    Returns 2D array on the (θ, v) meshgrid.
    """

def compute_separatrix(model: SteeringModel, params: ModelParams,
                        saddle: FixedPoint2D,
                        t_span: float = 50.0,
                        direction: str = "both") -> ndarray:
    """
    Integrate along the stable/unstable manifold of a saddle point
    to trace the separatrix (homoclinic or heteroclinic orbit).
    Uses eigenvector perturbation from the saddle.
    """
```

### 5.4 Module: `steering.analysis.homoclinic`

```
def duffing_homoclinic(c1: float, c3: float, t_array: ndarray,
                        branch: str = "positive") -> tuple[ndarray, ndarray]:
    """
    Analytical homoclinic orbit for the Duffing truncation.
    θ₀(t) = ±√(2c₁/|c₃|) sech(√c₁ t)
    v₀(t) = ∓√(2c₁/|c₃|) √c₁ sech(√c₁ t) tanh(√c₁ t)
    Requires c₁ > 0, c₃ < 0.
    """

def numerical_homoclinic(model: SteeringModel, params: ModelParams,
                          gamma: float = 0.0,
                          saddle: FixedPoint2D | None = None,
                          **kwargs) -> tuple[ndarray, ndarray, ndarray]:
    """
    Numerically compute the homoclinic orbit for any model by integrating
    from a small perturbation along the unstable eigenvector of the saddle,
    in the undamped (γ=0) Hamiltonian system.
    Returns (t, theta, v).
    """
```

### 5.5 Module: `steering.analysis.melnikov`

```
def melnikov_analytical(c1: float, c3: float, gamma: float, omega: float,
                         F: float, t0_array: ndarray) -> ndarray:
    """
    Evaluate the Melnikov function M(t₀) for the Duffing reduction
    with additive forcing F cos(ωt).
    """

def melnikov_critical_forcing(c1: float, c3: float, gamma: float,
                                omega: float) -> float:
    """
    Return F_crit(ω) = (4γc₁)/(3πω√(2|c₃|)) · cosh(πω/(2√c₁)).
    """

def melnikov_numerical(model: SteeringModel, params: ModelParams,
                        gamma: float, omega: float, F: float,
                        t0_array: ndarray,
                        homoclinic_t: ndarray | None = None) -> ndarray:
    """
    Numerically compute the Melnikov integral along the (numerical) homoclinic
    orbit of the given model. For the Bessel model, this captures higher-order
    corrections beyond the Duffing truncation.

    M(t₀) = -γ ∫ v₀(t)² dt + F ∫ v₀(t) cos(ω(t + t₀)) dt

    where (θ₀(t), v₀(t)) is the unperturbed homoclinic orbit.
    """
```

### 5.6 Module: `steering.analysis.poincare`

```
def stroboscopic_section(result: SimulationResult, omega: float,
                          transient_periods: int = 100) -> ndarray:
    """
    Extract stroboscopic Poincaré section: sample state at t = 2πn/ω.
    Discard the first `transient_periods` crossings.
    Returns array of shape (n_crossings, state_dim).
    """

def bifurcation_diagram_poincare(dynamics: AccelerationDynamics,
                                   model: SteeringModel,
                                   params: ModelParams,
                                   gamma: float, omega: float,
                                   sweep_param: str,
                                   sweep_values: ndarray,
                                   F: float | None = None,
                                   state0: ndarray | None = None,
                                   n_periods_transient: int = 200,
                                   n_periods_record: int = 100) -> dict:
    """
    For each value of the sweep parameter (e.g., F, gamma, omega, or delta),
    run a long simulation, extract the stroboscopic Poincaré section,
    and record the θ (or v) values. This produces the standard
    period-doubling / chaos bifurcation diagram.

    Supports continuation from final state of previous run for smooth tracking.
    """
```

### 5.7 Module: `steering.analysis.lyapunov`

```
def max_lyapunov_exponent(dynamics: AccelerationDynamics,
                           model: SteeringModel,
                           params: ModelParams,
                           gamma: float, F: float, omega: float,
                           state0: ndarray,
                           t_total: float = 1000.0,
                           renorm_interval: float = 1.0) -> float:
    """
    Compute the maximal Lyapunov exponent via the standard algorithm:
    integrate the variational equation alongside the trajectory,
    periodically renormalizing the tangent vector.
    """
```

---

## 6. Visualization

### 6.1 Module: `steering.visualization.profiles`

```
def plot_steering_drive(models: list[SteeringModel],
                         params_list: list[ModelParams],
                         theta_range: tuple = (-pi, pi),
                         labels: list[str] | None = None,
                         ax=None) -> matplotlib.axes.Axes:
    """
    Plot U(θ) for one or more models on the same axes.
    Used for comparing Duffing vs Bessel vs full model.
    """

def plot_potential(models, params_list, **kwargs):
    """Same but for V(θ) = -∫U dθ."""
```

### 6.2 Module: `steering.visualization.phase_portrait`

```
def plot_phase_portrait(dynamics: AccelerationDynamics,
                         model: SteeringModel,
                         params: ModelParams,
                         gamma: float = 0.0,
                         theta_range: tuple = (-pi, pi),
                         v_range: tuple = (-3, 3),
                         n_grid: int = 20,
                         overlay_hamiltonian: bool = True,
                         overlay_separatrix: bool = True,
                         overlay_fixed_points: bool = True,
                         trajectories: list[ndarray] | None = None,
                         ax=None) -> matplotlib.axes.Axes:
    """
    Quiver plot of the vector field with optional overlays.
    """
```

### 6.3 Module: `steering.visualization.bifurcation_plot`

```
def plot_bifurcation_2d(diagram: BifurcationDiagram,
                         param1_label: str, param2_label: str,
                         overlay_numerical: dict | None = None,
                         ax=None):
    """
    Contour plot of c₁ = 0 in the 2D parameter plane,
    shaded by sgn(c₃) (supercritical vs subcritical).
    Optionally overlay numerical fixed point branches.
    """

def plot_poincare_bifurcation(data: dict,
                                sweep_label: str,
                                ax=None):
    """
    Standard Poincaré bifurcation diagram (e.g. θ* vs F).
    """
```

### 6.4 Module: `steering.visualization.poincare_plot`

```
def plot_poincare_section(section: ndarray,
                           xlabel=r"$\theta$", ylabel=r"$v$",
                           ax=None):
    """
    Scatter plot of stroboscopic Poincaré section points.
    For strange attractors vs period-n orbits.
    """
```

### 6.5 Module: `steering.visualization.figures`

Composite figure builders that assemble multi-panel figures for the final report:

```
def figure_1_bifurcation(bessel_model, duffing_model, params, ...):
    """
    Main panel: bifurcation curve in (κ, δ) plane, two models overlaid.
    Inset panels: example U(θ) in each regime.
    """

def figure_2_phase_portraits(models, params_list, gamma, ...):
    """
    Side-by-side phase portraits: Duffing vs Bessel vs full circuit.
    With Hamiltonian level sets and separatrices.
    """

def figure_3_melnikov(c1, c3, gamma, omega_range, delta_range, ...):
    """
    F_crit vs δ at fixed ω and κ.
    F_crit vs ω at fixed δ and κ.
    Comparison between analytical (Duffing) and numerical (Bessel) Melnikov.
    """

def figure_4_poincare(results_duffing, results_bessel, results_full, ...):
    """
    Bifurcation diagram (Poincaré section values vs F or ω).
    Example Poincaré sections at representative parameter values.
    Three-model comparison.
    """
```

---

## 7. Repository Structure

```
drosophila-steering/
├── SPEC.md                          # This file
├── pyproject.toml                   # Project metadata, dependencies
├── README.md
│
├── src/
│   └── steering/
│       ├── __init__.py
│       ├── params.py                # ModelParams, ForcingParams dataclasses
│       ├── nonlinearities.py        # Nonlinearity registry
│       │
│       ├── models/
│       │   ├── __init__.py          # Re-exports all models
│       │   ├── base.py              # SteeringModel ABC
│       │   ├── duffing.py           # DuffingModel
│       │   ├── bessel.py            # BesselSteeringModel
│       │   ├── continuous.py        # ContinuousPFLModel
│       │   ├── discrete.py          # DiscretePFLModel
│       │   └── full_circuit.py      # FullCircuitModel (with PFL2)
│       │
│       ├── dynamics/
│       │   ├── __init__.py
│       │   ├── base.py              # Dynamics ABC
│       │   ├── velocity.py          # VelocityDynamics
│       │   └── acceleration.py      # AccelerationDynamics (handles forcing + topology)
│       │
│       ├── integrator.py            # Simulation, SimulationResult
│       ├── sweep.py                 # ParameterSweep
│       │
│       ├── analysis/
│       │   ├── __init__.py
│       │   ├── fixed_points.py
│       │   ├── bifurcation.py
│       │   ├── hamiltonian.py
│       │   ├── homoclinic.py
│       │   ├── melnikov.py
│       │   ├── poincare.py
│       │   └── lyapunov.py
│       │
│       └── visualization/
│           ├── __init__.py
│           ├── style.py             # matplotlib style, color palettes, shared config
│           ├── profiles.py          # U(θ) and V(θ) plots
│           ├── phase_portrait.py
│           ├── bifurcation_plot.py
│           ├── poincare_plot.py
│           └── figures.py           # Composite multi-panel figures
│
├── notebooks/
│   ├── 00_model_walkthrough.ipynb           # Model hierarchy tour + consistency checks
│   ├── 01_steering_profiles.ipynb           # U(θ) and V(θ) across parameter space
│   ├── 02_pitchfork_bifurcation.ipynb       # Bifurcation analysis (velocity + acceleration)
│   ├── 03_phase_portraits_hamiltonian.ipynb  # Phase portraits, Hamiltonian level sets, separatrices
│   ├── 04_homoclinic_orbits.ipynb           # Analytical vs numerical homoclinics, validity
│   ├── 05_melnikov_analysis.ipynb           # Melnikov function, chaos thresholds
│   ├── 06_poincare_chaos.ipynb              # Poincaré sections, bifurcation diagrams, Lyapunov
│   ├── 07_topology_comparison.ipynb         # R² vs S¹×R: when does topology matter?
│   ├── 08_discrete_neurons.ipynb            # Finite-N effects, biological N=8
│   ├── 09_pfl2_indirect_pathway.ipynb       # Full circuit with PFL2 modulation
│   └── 10_extensions.ipynb                  # Asymmetric goals, parametric forcing, noise
│
├── scripts/
│   ├── fig1_bifurcation.py          # Generate Figure 1
│   ├── fig2_phase_portraits.py      # Generate Figure 2
│   ├── fig3_melnikov.py             # Generate Figure 3
│   ├── fig4_poincare_chaos.py       # Generate Figure 4
│   └── exploratory/
│       ├── sweep_S_A.py             # Effect of S, A on bifurcation structure
│       ├── discrete_vs_continuous.py # Finite-N effects
│       ├── pfl2_effects.py          # Indirect pathway exploration
│       └── parametric_forcing.py    # δ(t) forcing variant
│
├── tests/
│   ├── test_models.py               # Consistency checks across model levels
│   ├── test_bessel_coefficients.py  # Taylor coefficient formulas vs finite diff
│   ├── test_symmetry.py             # U(-θ) = -U(θ), U(0) = 0
│   ├── test_limits.py               # κ→0, κ→∞, N→∞ limits
│   ├── test_dynamics.py             # Conservation of H in undamped system
│   ├── test_melnikov.py             # Analytical vs numerical Melnikov
│   └── test_topology.py            # S^1 wrapping correctness
│
└── data/
    ├── .gitkeep                     # For saved simulation results
    └── README.md                    # Documents caching conventions for notebooks
```

---

## 8. Dependencies

```toml
[project]
name = "drosophila-steering"
requires-python = ">=3.11"
dependencies = [
    "numpy>=1.26",
    "scipy>=1.12",
    "matplotlib>=3.8",
    "xarray>=2024.1",     # For labeled parameter sweep results
]

[project.optional-dependencies]
dev = [
    "pytest>=8.0",
    "pytest-xdist",       # Parallel test execution
    "ipython",
    "jupyter",
    "ipywidgets>=8.0",    # Interactive sliders in notebooks
    "jupyterlab>=4.0",
]
```

---

## 9. Key Implementation Notes

### 9.1 Bessel Function Numerics

`scipy.special.i0`, `i1`, `iv` can overflow for large arguments ($\kappa \gtrsim 300$). Use the **exponentially scaled** versions `i0e`, `i1e`, `ive` (which return $e^{-|z|}I_n(z)$) and track the exponential prefactors analytically. Since $U$ involves differences of $I_0$ at similar arguments, the exponential factors approximately cancel and the result is well-conditioned even for large $\kappa$. Implement a utility:

```python
# steering/utils/bessel.py
def i0_safe(z):
    """I₀(z) using exponentially scaled form for large z."""
    z = np.asarray(z, dtype=float)
    return np.where(z < 500, scipy.special.i0(z), scipy.special.i0e(z) * np.exp(z))
```

### 9.2 Topology Handling

In `AccelerationDynamics`, topology is enforced at the integration level via a wrapper around the ODE solver's step output, not inside the RHS function. This avoids discontinuities in the derivative computation.

For `S^1` analysis: after each integration step, wrap $\theta$ to $[-\pi, \pi)$. The velocity $v$ is not wrapped.

For Poincaré section comparisons: when comparing $S^1 \times \mathbb{R}$ vs $\mathbb{R}^2$ results, fold the $\mathbb{R}^2$ trajectory's $\theta$ values modulo $2\pi$ before computing section statistics.

### 9.3 Consistency Tests

The following identity chains must hold numerically (up to integration tolerance):

1. `BesselSteeringModel` with `params` $\equiv$ `ContinuousPFLModel` with `nonlinearity="quadratic"` and same `params`, for all $\theta$.
2. `DuffingModel.from_params(params).steering_drive(theta)` $\approx$ `BesselSteeringModel.steering_drive(theta)` for $|\theta| \ll 1$.
3. `DiscretePFLModel(N=1000)` $\approx$ `ContinuousPFLModel` for all $\theta$.
4. `FullCircuitModel(W_D3=0)` $\equiv$ `DiscretePFLModel` or `ContinuousPFLModel` (depending on `N_neurons`).
5. Hamiltonian is conserved ($|H(t) - H(0)| < \epsilon$) under undamped, unforced `AccelerationDynamics`.

### 9.4 Performance Considerations

- The `DiscretePFLModel` and `FullCircuitModel` are called $O(10^6)$ times in a typical parameter sweep. Vectorize over $\theta$ arrays and pre-compute the neuron preferred headings.
- For Poincaré bifurcation diagrams: use the dense output of `solve_ivp` to interpolate exactly at stroboscopic times $t = 2\pi n/\omega$, rather than using a fixed `t_eval` grid. This is both more accurate and faster.
- Parameter sweeps should use `multiprocessing.Pool` (or `concurrent.futures.ProcessPoolExecutor`) since each parameter point is independent.

### 9.5 Numerical Homoclinic Computation

For models beyond Duffing where no closed-form homoclinic exists:

1. Find the saddle at $\theta = 0$ and compute its unstable eigenvector $\mathbf{e}_u$ (with eigenvalue $\lambda_u > 0$).
2. Initialize at $\mathbf{x}_0 = (0, 0) + \epsilon \mathbf{e}_u$ with $\epsilon \sim 10^{-8}$.
3. Integrate forward and backward in time until the trajectory returns close to the saddle (distance $< \epsilon_{\text{tol}}$).
4. The forward branch traces one wing; the backward branch (with $-\epsilon$) traces the other.
5. For the Melnikov integral, the orbit only needs to be accurate where $|v_0(t)|$ is appreciable — the exponential decay near the saddle means truncation errors in the tails are harmless.

### 9.6 Forcing Types

Implement two forcing entry points in `AccelerationDynamics`:

- **Additive**: $\ddot{\theta} + \gamma\dot{\theta} = U(\theta, \delta_0) + F\cos\omega t$. The forcing is a direct torque.
- **Parametric (through $\delta$)**: $\ddot{\theta} + \gamma\dot{\theta} = U(\theta, \delta_0 + F\cos\omega t)$. The goals oscillate. This requires re-evaluating $U$ at each timestep with a time-dependent $\delta$. The Melnikov analysis for this case involves $\partial U/\partial\delta$ evaluated along the homoclinic, which the code should also compute.

### 9.7 On the Validity Regime

Every analysis function should accept and propagate an optional `validity_check: bool = True` flag. When enabled:

- For `DuffingModel` analyses: warn if the homoclinic orbit amplitude $\theta_{\max} = \sqrt{2c_1/|c_3|}$ exceeds a threshold (default: $\pi/4$), indicating the cubic truncation and/or the $\mathbb{R}^2$ approximation may be unreliable.
- For Melnikov results: warn if $\gamma$ or $F$ are not small relative to $|c_1|$ (the unperturbed energy scale), violating the $O(\epsilon)$ assumption.
- Attach these warnings to result objects so they can be displayed or logged.

---

## 10. Figure Specifications

### Figure 1: Pitchfork Bifurcation in the $(\kappa, \delta)$ Plane

**Layout**: 1 main panel + 2 inset panels.

- **Main panel**: Heatmap or contour plot of $c_1(\kappa, \delta)$ with the bifurcation curve $c_1 = 0$ as a solid contour. Shade by $\mathrm{sgn}(c_3)$: one color for supercritical, another for subcritical. Overlay the numerical bifurcation curve from `BesselSteeringModel` (dotted) and optionally from `DiscretePFLModel` (dashed, at biological $N$).
- **Inset A**: $U(\theta)$ for parameter values in the monostable regime ($c_1 < 0$). Show Duffing (dashed) and Bessel (solid).
- **Inset B**: $U(\theta)$ for parameter values in the bistable regime ($c_1 > 0$). Show double-well structure, mark stable/unstable fixed points.
- **Axes**: $\kappa \in [0.5, 10]$ (x), $\delta \in [0, \pi/2]$ (y). Fix $S = 1, A = 1, \Delta = 3\pi/8$.

### Figure 2: Phase Portraits and Hamiltonian Structure

**Layout**: 2 × 2 panels (or 1 × 2 if space is tight).

- **Top left**: Phase portrait for Duffing model, undamped. Hamiltonian level sets as contours, separatrix (homoclinic) in bold, fixed points marked.
- **Top right**: Same for `BesselSteeringModel` (full $U(\theta)$ on $[-\pi, \pi]$).
- **Bottom left**: Phase portrait for Duffing model, lightly damped ($\gamma$ small). Show a few representative trajectories spiraling into the stable foci.
- **Bottom right**: Same for `BesselSteeringModel`.
- Mark centers, saddles with standard symbols. Note the extra saddle at $\theta = \pi$ in the Bessel model that is absent in Duffing.

### Figure 3: Melnikov Analysis and Chaos Threshold

**Layout**: 2 panels.

- **Left**: $F_{\mathrm{crit}}$ vs $\omega$ at fixed $\delta, \kappa, \gamma$. Show analytical (Duffing) curve and numerical (Bessel homoclinic) curve. The $\cosh$ envelope should be visible.
- **Right**: $F_{\mathrm{crit}}$ vs $\delta$ at fixed $\omega, \kappa, \gamma$. This is the biologically interesting plot: which goal separations are most susceptible to chaos?

### Figure 4: Poincaré Sections and Chaos

**Layout**: 2 rows.

- **Top row**: Poincaré bifurcation diagram ($\theta$ at stroboscopic times vs $F$) for (a) Duffing on $\mathbb{R}^2$, (b) Bessel on $S^1 \times \mathbb{R}$, (c) `FullCircuitModel`.
- **Bottom row**: Example Poincaré sections ($\theta$ vs $v$) at representative $F$ values showing (a) period-1, (b) period-doubling, (c) strange attractor.

---

## 11. Extension Points

The following are not core deliverables but should be easy to add given the architecture:

- **Unequal goal amplitudes**: Replace the symmetric two-bump goal profile with $A_1 \exp(\kappa\cos(\delta - \phi)) + A_2 \exp(\kappa\cos(-\delta - \phi))$. Breaks the $\theta \to -\theta$ symmetry, turning the pitchfork into an imperfect bifurcation (cusp unfolding).
- **Three or more goals**: Generalize the goal profile to $\sum_i A_i \exp(\kappa\cos(\delta_i - \phi))$.
- **Noise**: Add Wiener process to `AccelerationDynamics` for Kramers escape rate analysis between wells.
- **Slow goal dynamics**: Let $\delta(t)$ evolve on a slow timescale, creating a slow-fast system for adiabatic analysis.

---

## 12. Notebook Specifications

Each notebook is a self-contained exploration of one aspect of the project. They are ordered to follow the logical flow of the analysis: build intuition for the models, then bifurcation structure, then Hamiltonian/phase portrait structure, then chaos. The later notebooks (08–10) branch into extensions. All notebooks share a common preamble:

```python
import numpy as np
import matplotlib.pyplot as plt
from steering.params import ModelParams, ForcingParams
from steering.models import DuffingModel, BesselSteeringModel, ContinuousPFLModel, DiscretePFLModel, FullCircuitModel
from steering.dynamics import VelocityDynamics, AccelerationDynamics
from steering.integrator import Simulation
# + analysis/visualization imports as needed
```

Notebooks should use `%matplotlib widget` for interactive figures where useful (phase portraits, parameter sliders), with `%matplotlib inline` as fallback.

---

### `00_model_walkthrough.ipynb` — Model Hierarchy Tour

**Purpose**: Verify that the codebase is working and build intuition for what each model level does.

**Cells**:

1. **Instantiate all five models** at a common set of parameters (e.g., $\kappa = 3$, $\delta = \pi/4$, $S = 1$, $A = 1$, $\Delta = 3\pi/8$). For the discrete model, use $N = 8$ (biological) and $N = 64$ (near-continuous).

2. **Plot $U(\theta)$ for all models on one axis**. This is the first sanity check: Duffing should match near $\theta = 0$; Bessel and ContinuousPFL (with $f = x^2$) should be identical; DiscretePFL($N = 64$) should be close to continuous; DiscretePFL($N = 8$) will show visible discretization artifacts.

3. **Consistency test table**: For a grid of $\theta$ values, print $|U_{\text{Bessel}}(\theta) - U_{\text{ContinuousPFL}}(\theta)|$ and $|U_{\text{Duffing}}(\theta) - U_{\text{Bessel}}(\theta)|/|U_{\text{Bessel}}(\theta)|$. Verify the former is $< 10^{-10}$ and the latter grows with $|\theta|$.

4. **Nonlinearity comparison**: Fix one parameter set in the bistable regime. Plot $U(\theta)$ for $f = x^2$ (Bessel-matchable), $f = \text{ELU}$, $f = \text{ReLU}$, $f = \text{softplus}$ using `ContinuousPFLModel`. Observe how the nonlinearity choice changes the shape of $U$ — in particular, whether the locations and depths of the wells shift. This motivates why the Bessel closed form (which assumes $f = x^2$) is an approximation even within the continuous-population model when the true nonlinearity is ELU.

5. **FullCircuitModel at $W_{D3} = 0$**: Verify it reproduces the DiscretePFL or ContinuousPFL result (consistency test §9.3 item 4).

---

### `01_steering_profiles.ipynb` — $U(\theta)$ and $V(\theta)$ Landscape

**Purpose**: Map out how the steering drive and potential landscape change across parameter space, building geometric intuition for the bifurcation analysis.

**Cells**:

1. **$U(\theta)$ heatmap**: On a $(\theta, \delta)$ grid with $\theta \in [-\pi, \pi]$ and $\delta \in [0, \pi/2]$, plot $U(\theta, \delta)$ as a heatmap (fixed $\kappa$). The zero contours of this heatmap are the fixed points as a function of $\delta$ — this is a direct visualization of the pitchfork.

2. **$V(\theta)$ cross-sections**: For $\delta$ values below, at, and above bifurcation, plot $V(\theta) = -\int_0^\theta U(\theta')\,d\theta'$. Animate or use a slider over $\delta$. Watch the single well split into a double well.

3. **Same plots for varying $\kappa$** at fixed $\delta$. The sharpening of bumps ($\kappa$ increasing) deepens the wells and narrows them — this is the concentration parameter's role.

4. **Bessel vs Duffing overlay**: For each cross-section, show the Duffing cubic approximation overlaid. Annotate where the approximation breaks down (quantify relative error).

5. **Effect of $S/A$ ratio**: Hold $SA$ constant and vary the ratio. This changes how much the steering drive is dominated by the cross-term (heading × goal) vs the pure goal term. Explore whether the bifurcation point shifts.

---

### `02_pitchfork_bifurcation.ipynb` — Bifurcation Analysis

**Purpose**: Core analysis — compute and visualize the pitchfork bifurcation in both velocity and acceleration control.

**Cells**:

1. **Compute $c_1(\kappa, \delta)$ and $c_3(\kappa, \delta)$** over a grid using `BesselSteeringModel.taylor_coefficients()`. Store as 2D arrays.

2. **Bifurcation diagram in $(\kappa, \delta)$ plane**: Contour $c_1 = 0$. Shade by $\operatorname{sgn}(c_3)$. This is the main panel of Figure 1. Annotate the supercritical and subcritical regions. Verify: in the broad-bump limit ($\kappa \to 0$), the bifurcation condition reduces to $\sin\Delta\cos\delta = 0$ (analytical anchor from previous conversation).

3. **Overlay numerical bifurcation locus**: For each $\kappa$ on the grid, numerically find the $\delta$ at which the origin changes stability by rootfinding on $U'(0, \delta) = 0$ using the `BesselSteeringModel`. Plot as a dotted line — should overlay exactly on the $c_1 = 0$ contour. Then do the same with `DiscretePFLModel(N=8)` and `ContinuousPFLModel(f=ELU)`. These will deviate, and the deviations are scientifically interesting.

4. **1D bifurcation diagrams**: Fix $\kappa$ and sweep $\delta$. For each $\delta$, find all fixed points of $U(\theta) = 0$ and classify stability. Plot the classic pitchfork diagram (stable branches as solid, unstable as dashed). Do this for Duffing, Bessel, and DiscretePFL on the same axes.

5. **Velocity vs acceleration control**: For the acceleration system $\ddot\theta + \gamma\dot\theta = U(\theta)$, the fixed points are the same but stability classification uses eigenvalues of the $2\times 2$ Jacobian. Show that the pitchfork locus is identical but the dynamics near it differ (nodes vs spirals depending on $\gamma$). Compute the critical damping $\gamma_c = 2\sqrt{|c_1|}$ at which the approach to the stable fixed points transitions from oscillatory to overdamped.

6. **Bifurcation theorem verification**: At a point on the bifurcation curve, verify the nondegeneracy conditions for the pitchfork: (i) $c_1 = 0$, (ii) $\partial c_1 / \partial \delta \neq 0$ (transversality), (iii) $c_3 \neq 0$. Compute these numerically and print a pass/fail table.

---

### `03_phase_portraits_hamiltonian.ipynb` — Phase Portraits and Hamiltonian Structure

**Purpose**: Visualize the 2D phase space of the acceleration system, identify the qualitative regimes, and compare models.

**Cells**:

1. **Parameter selection**: Choose 3 representative parameter sets: (A) monostable ($c_1 < 0$), (B) weakly bistable ($c_1 > 0$ small, near bifurcation), (C) strongly bistable ($c_1 > 0$ large). Print $(c_1, c_3, \theta_{\max})$ for each.

2. **Duffing phase portraits (undamped)**: For each parameter set, plot the Hamiltonian $H(\theta, v)$ as filled contours. Overlay the separatrix (homoclinic orbit) in bold. Mark centers and saddles. For (A): single center at origin, no separatrix. For (B, C): saddle at origin, two centers at $\pm\sqrt{-c_1/c_3}$, figure-eight homoclinic.

3. **Bessel phase portraits (undamped)**: Same layout, but using the numerically-computed Hamiltonian $H = \frac{1}{2}v^2 + V(\theta)$ where $V$ is obtained from integrating the full Bessel $U(\theta)$. Key difference: the Bessel model has additional saddles/centers at $\theta \approx \pm\pi$ that the Duffing model misses. For regime (C), the separatrix topology may differ globally (heteroclinic connections between $\theta = 0$ and $\theta = \pi$ saddles).

4. **Side-by-side comparison**: 2×3 grid — Duffing on top, Bessel on bottom, one column per regime. This is Figure 2 material.

5. **Damped trajectories**: Add $\gamma = 0.1$ (or appropriate value). Integrate several initial conditions and overlay trajectories on the Hamiltonian contours. Show spiral-in to the stable foci. Compare Duffing vs Bessel — do trajectories converge to the same attractors?

6. **$S^1 \times \mathbb{R}$ phase cylinder**: For the Bessel model, plot the phase portrait on the cylinder by tiling $[-\pi, \pi]$ periodically. Visualize the phase cylinder as a rectangle with identified left/right edges. The periodicity makes the $\theta = \pm\pi$ saddle visible and shows how the two homoclinic loops are globally connected.

---

### `04_homoclinic_orbits.ipynb` — Homoclinic Orbit Construction and Validation

**Purpose**: Compute the homoclinic orbits that enter the Melnikov integral, analytically (Duffing) and numerically (Bessel/full), and quantify the validity of the Duffing approximation along the orbit.

**Cells**:

1. **Analytical Duffing homoclinic**: Plot $\theta_0(t) = \sqrt{2c_1/|c_3|}\operatorname{sech}(\sqrt{c_1}\,t)$ and $v_0(t)$ for a representative parameter set. Mark $\theta_{\max}$ and the characteristic timescale $1/\sqrt{c_1}$.

2. **Numerical Bessel homoclinic**: Using the saddle eigenvector perturbation method (§9.5), compute the homoclinic orbit of the undamped Bessel system. Overlay on the Duffing analytical orbit. Quantify the discrepancy in $\theta_{\max}$, orbit shape, and period (time to traverse the loop, formally $\infty$ for a true homoclinic, but measured as time between threshold crossings).

3. **Pointwise error along the orbit**: Plot $|U_{\text{Duffing}}(\theta_0(t)) - U_{\text{Bessel}}(\theta_0(t))|$ as a function of $t$. Verify that the error is peaked at $t = 0$ (the orbit's apex) and exponentially small in the tails. This directly shows why the Melnikov integral is robust to the Taylor truncation error.

4. **Validity vs parameters**: Sweep over $c_1$ (distance from bifurcation). For each value, compute $\theta_{\max}$ and the relative $L^2$ error between the Duffing and Bessel homoclinics:
   $$
   \varepsilon_{\text{rel}} = \frac{\|\theta_{0,\text{Duffing}}(t) - \theta_{0,\text{Bessel}}(t)\|_{L^2}}{\|\theta_{0,\text{Bessel}}(t)\|_{L^2}}
   $$
   Plot $\varepsilon_{\text{rel}}$ vs $\theta_{\max}/\pi$. This gives a quantitative "validity boundary" for the Duffing reduction.

5. **Topology check**: For parameter values where $\theta_{\max}$ approaches $O(1)$, compare the $\mathbb{R}^2$ homoclinic with the $S^1 \times \mathbb{R}$ separatrix. When do they differ? Does the $S^1$ separatrix connect to the $\theta = \pi$ saddle instead of returning to $\theta = 0$? This is the topological transition where the local Duffing picture qualitatively breaks down.

---

### `05_melnikov_analysis.ipynb` — Melnikov Function and Chaos Thresholds

**Purpose**: Compute the Melnikov function analytically (Duffing) and numerically (Bessel), determine $F_{\text{crit}}$ curves, and compare with direct simulation.

**Cells**:

1. **Analytical Melnikov for Duffing**: Plot $M(t_0)$ for several $(F, \omega, \gamma)$ values. Show cases where $M$ has simple zeros (chaos) and where it doesn't (no horseshoe). Annotate the damping and forcing contributions separately.

2. **$F_{\text{crit}}(\omega)$ curve (Duffing)**: Plot the analytical threshold $F_{\text{crit}} = \frac{4\gamma c_1}{3\pi\omega\sqrt{2|c_3|}}\cosh\!\left(\frac{\pi\omega}{2\sqrt{c_1}}\right)$. Discuss the $\cosh$ envelope: high-frequency forcing is exponentially ineffective. Mark the minimum of $F_{\text{crit}}(\omega)$ — this is the "most dangerous" frequency.

3. **Numerical Melnikov for Bessel**: Compute the Melnikov integral numerically using the Bessel homoclinic orbit from notebook 04. Compare $F_{\text{crit}}$ curves between Duffing-analytical and Bessel-numerical. The discrepancy quantifies how much the higher-order terms in $U(\theta)$ shift the chaos threshold.

4. **$F_{\text{crit}}(\delta)$ at fixed $\omega$**: The biologically interesting plot. As goal separation $\delta$ increases from 0, $c_1$ increases (the double well deepens), which lowers $F_{\text{crit}}$ — widely separated goals are more susceptible to chaotic steering. But at very large $\delta$, the cubic truncation breaks down. Show this non-monotonic structure.

5. **$F_{\text{crit}}(\delta)$ at fixed $\omega$, multiple $\kappa$**: Overlay curves for $\kappa = 1, 2, 3, 5$. Show that sharper bumps (higher $\kappa$) lead to deeper wells and lower chaos thresholds — tighter goal representations are paradoxically more chaos-prone.

6. **Parametric forcing variant**: Compute the Melnikov integral for forcing through $\delta(t) = \delta_0 + F\cos\omega t$. The perturbation Hamiltonian is $h_1(\theta, v, t) = -\gamma v^2 + \frac{\partial U}{\partial \delta}\big|_{\delta_0} F\cos(\omega t) \cdot \theta$ (to leading order). The resulting Melnikov integral involves different moments of the homoclinic. Compare $F_{\text{crit}}$ between additive and parametric forcing — which is more effective at inducing chaos?

7. **Validation via direct simulation**: Pick 2-3 $(F, \omega)$ points: one below $F_{\text{crit}}$, one near, one above. Run long simulations and check whether the trajectory is periodic or chaotic (by Lyapunov exponent or visual inspection of Poincaré section). Does the Melnikov prediction hold?

---

### `06_poincare_chaos.ipynb` — Poincaré Sections and Route to Chaos

**Purpose**: Construct stroboscopic Poincaré maps, visualize the strange attractor, and trace the route to chaos via period-doubling cascades.

**Cells**:

1. **Poincaré bifurcation diagram vs $F$**: Fix $\omega, \gamma, \delta, \kappa$. Sweep $F$ from 0 to well above $F_{\text{crit}}$. For each $F$, run a long simulation ($\sim 500$ forcing periods), discard transients, record $\theta$ at stroboscopic times. Plot $\theta_{\text{strobe}}$ vs $F$. Expect: period-1 → period-2 → period-4 → ... → chaos, possibly with periodic windows.

2. **Three-model comparison**: Same sweep for (a) Duffing on $\mathbb{R}^2$, (b) Bessel on $\mathbb{R}^2$, (c) Bessel on $S^1 \times \mathbb{R}$. Side-by-side panels. Key question: does the topology ($S^1$) change the bifurcation sequence? (It may introduce additional periodic orbits winding around the cylinder.)

3. **Example Poincaré sections**: At 3-4 representative $F$ values (period-1, period-4, strange attractor, periodic window), plot the full $(\theta, v)$ Poincaré section. For the strange attractor, use $\sim 10^4$ stroboscopic points to resolve the fractal structure.

4. **Lyapunov exponent vs $F$**: Compute $\lambda_{\max}(F)$ across the same sweep. Overlay on the bifurcation diagram. $\lambda_{\max} > 0$ regions should correspond exactly to the chaotic bands. Check consistency with the Melnikov prediction for the onset.

5. **Poincaré bifurcation diagram vs $\omega$**: Fix $F$ above threshold and sweep $\omega$. This reveals resonance structure — certain frequencies produce larger chaotic windows. Compare with the $F_{\text{crit}}(\omega)$ curve from notebook 05.

6. **Poincaré bifurcation diagram vs $\delta$**: Fix $F, \omega$ and sweep goal separation. This is the most biologically interpretable: which goal separations produce chaotic steering? Compare the onset with the Melnikov $F_{\text{crit}}(\delta)$ prediction.

7. **Full circuit comparison**: Repeat the $F$-sweep bifurcation diagram for `FullCircuitModel` with $W_{D3} > 0$. Does the indirect pathway (PFL2 gain modulation) suppress or enhance chaos? The hypothesis is that PFL2's adaptive gain control regularizes the dynamics, raising $F_{\text{crit}}$.

---

### `07_topology_comparison.ipynb` — $\mathbb{R}^2$ vs $S^1 \times \mathbb{R}$

**Purpose**: Systematically investigate when the planar approximation breaks down and the cylindrical topology produces qualitatively different dynamics.

**Cells**:

1. **Trajectory comparison**: For parameters well inside the Duffing validity regime ($\theta_{\max} \ll \pi$), simulate the forced damped system on both $\mathbb{R}^2$ and $S^1 \times \mathbb{R}$. Overlay trajectories. They should be indistinguishable — this is the "boring" (but important) baseline.

2. **Push toward breakdown**: Increase $c_1$ (move far from bifurcation) until $\theta_{\max} \sim \pi/2$. Repeat the comparison. Where do trajectories first diverge? Is it during chaotic transients (where small differences amplify) or even in the periodic regime?

3. **Winding number**: In the chaotic regime on $S^1 \times \mathbb{R}$, the trajectory can wind around the cylinder. Define the winding number $W(T) = \frac{1}{2\pi}\int_0^T \dot\theta\,dt$ (net rotations). Plot $W(T)/T$ vs $T$. For $\mathbb{R}^2$ this is trivially near zero (bounded orbits); for $S^1 \times \mathbb{R}$ it can be nonzero, indicating net rotational drift — the fly spins persistently while switching between goals. This is a qualitatively new phenomenon absent in the planar picture.

4. **Phase portrait on the cylinder**: Visualize the $S^1 \times \mathbb{R}$ phase space as a rectangle $[-\pi, \pi] \times [-v_{\max}, v_{\max}]$ with identified edges. Plot the Poincaré section on this. For chaotic trajectories, does the attractor wrap around the cylinder or stay local?

5. **Quantitative divergence metric**: For a grid of $(\kappa, \delta)$ values, compute the Hausdorff distance (or simpler: symmetric difference of Poincaré section point clouds) between $\mathbb{R}^2$ and $S^1 \times \mathbb{R}$ at fixed $(F, \omega, \gamma)$. Plot as a heatmap. This gives a "validity map" for the planar approximation, complementing the analytical $\theta_{\max}/\pi$ criterion.

---

### `08_discrete_neurons.ipynb` — Finite-$N$ Effects

**Purpose**: Explore how discretizing the neural population (from continuous integral to $N$ neurons) changes the dynamics.

**Cells**:

1. **$U(\theta)$ convergence**: Plot $U(\theta)$ for $N = 4, 8, 16, 32, 64, 256$ and the continuous limit. Quantify $L^\infty$ error vs $N$. At biological $N = 8$, how large is the discretization error relative to the Duffing truncation error?

2. **Bifurcation locus shift**: Compute the pitchfork bifurcation locus $c_1 = 0$ in $(\kappa, \delta)$ for $N = 8, 16, 64, \infty$. How much does the bifurcation curve shift at biological $N$? Is the shift larger or smaller than the shift from changing nonlinearity ($x^2$ vs ELU)?

3. **Broken symmetry at finite $N$**: The continuous model has exact $U(-\theta) = -U(\theta)$ symmetry. Does this hold at finite $N$? If the preferred headings $\phi_j$ are not symmetric about 0, the antisymmetry breaks and the pitchfork unfolds into a saddle-node + transcritical. Check by computing $U(0)$ (should be exactly 0 by symmetry) and $U''(0)$ (should be 0 for the pitchfork). Explore: for which $N$ and neuron arrangements is the symmetry preserved vs broken?

4. **Phase portrait at $N = 8$**: Show the phase portrait for the biological neuron count. How does it compare to the continuous Bessel model? Are the wells at slightly different locations?

5. **Poincaré section at $N = 8$**: Does the strange attractor look qualitatively similar to the continuous case, or does discretization introduce visible artifacts (e.g., preferred directions due to the finite neuron grid)?

---

### `09_pfl2_indirect_pathway.ipynb` — Full Circuit with PFL2 Modulation

**Purpose**: Explore how the indirect pathway (PFL2 → DNa03 → DNa02) modulates the steering dynamics.

**Cells**:

1. **PFL2 population activity**: Plot $\Sigma_{P2}(\theta, \delta)$ — the PFL2 population sum. Since PFL2 has a 180° phase shift, it peaks at the anti-goal. Verify this.

2. **Effect of $W_{D3}$ on $U(\theta)$**: Plot $U(\theta)$ for $W_{D3} = 0, 0.1, 0.5, 1.0, 2.0$ at a fixed bistable parameter set. The indirect pathway should amplify the steering drive at large heading errors (where PFL2 is active) and have little effect near the goals.

3. **Bifurcation locus with PFL2**: Compute the pitchfork locus in $(\kappa, \delta)$ for several $W_{D3}$ values. Does PFL2 shift the bifurcation? The hypothesis: PFL2 primarily affects the gain (slope of $U$ near fixed points) rather than the fixed point locations, so the bifurcation locus shifts only weakly but the stability properties (overdamped vs underdamped approach) change more.

4. **Phase portrait with PFL2**: Compare undamped phase portraits with and without PFL2. The potential $V(\theta)$ should be modified — does PFL2 make the wells deeper or shallower? Does it change the barrier height?

5. **Melnikov analysis with PFL2**: Compute the numerical Melnikov integral for the full circuit model. Does PFL2 raise or lower $F_{\text{crit}}$? If PFL2 acts as adaptive gain (increasing effective damping at large errors), it should increase $F_{\text{crit}}$, making chaos harder to achieve — a stabilizing role.

6. **Poincaré comparison**: Side-by-side Poincaré bifurcation diagrams (vs $F$) with $W_{D3} = 0$ and $W_{D3} > 0$. Does the period-doubling cascade shift? Are there new periodic windows?

7. **Speed-accuracy tradeoff**: The biological role of PFL2 is managing the speed-accuracy tradeoff. Quantify this: for a step change in goal ($\delta \to \delta + \epsilon$), measure the settling time and overshoot as a function of $W_{D3}$. The indirect pathway should reduce overshoot (higher accuracy) at the cost of slower convergence (lower speed), or vice versa depending on the regime.

---

### `10_extensions.ipynb` — Asymmetric Goals, Parametric Forcing, and Noise

**Purpose**: Sandbox notebook for extension investigations beyond the core deliverables.

**Cells**:

1. **Asymmetric goal amplitudes**: Set $A_1 \neq A_2$ in the goal profile. This breaks the $\theta \to -\theta$ symmetry. The pitchfork unfolds into an imperfect bifurcation (a saddle-node plus a persisting branch). Plot the bifurcation diagram in $(\delta, A_1/A_2)$ — this is a cusp (codimension-2) unfolding. Show the cusp point where the two saddle-node branches meet.

2. **Parametric forcing detailed study**: $\delta(t) = \delta_0 + F\cos\omega t$. Run the full forced system and compare Poincaré sections with the additive forcing case. Parametric forcing can produce parametric resonance at $\omega \approx 2\omega_0$ (where $\omega_0 = \sqrt{|c_1|}$ is the natural frequency of the wells) — check for this.

3. **Noise**: Add small Gaussian white noise $\sqrt{2D}\,\xi(t)$ to $\dot{v}$. Use Euler-Maruyama integration. In the bistable regime without forcing, observe noise-induced switching between wells (Kramers escape). Compute mean escape time vs noise intensity $D$ and compare with the Kramers rate $\sim \exp(\Delta V / D)$ where $\Delta V$ is the barrier height.

4. **Three goals**: Superpose three von Mises bumps at $\delta_1, \delta_2, \delta_3$. Explore the fixed point structure — does the system still reduce to a Duffing-like normal form near any of the bifurcations, or is the codimension higher? Plot $V(\theta)$ and the phase portrait for representative configurations.

5. **Slow $\delta$ dynamics**: Let $\dot\delta = \epsilon g(\delta, \theta)$ for some slow function $g$ (e.g., $g = -\delta + \delta_{\text{input}}$ with slow timescale $1/\epsilon$). This creates a slow-fast system. Explore whether the system tracks a moving bifurcation point (adiabatic regime) or exhibits delayed bifurcation / canard-like transitions.

---

### Notebook Conventions

- **Naming**: Notebooks are numbered `00`–`10` for ordering. The number reflects dependency: notebook $n$ may use results/intuitions from notebooks $< n$, but not vice versa.
- **Parameters**: Each notebook defines its parameters at the top in a clearly marked "Parameters" cell. Use a dictionary or `ModelParams` instance so it's easy to re-run with different values.
- **Figures**: Notebooks produce exploratory figures (not publication-quality). The `scripts/fig*` files produce the final versions. However, notebooks should use `steering.visualization.style` for consistent color schemes.
- **Saving results**: Expensive computations (parameter sweeps, long simulations) should cache results to `data/` as `.npz` or `.pkl` files with a descriptive name, and check for cached results before recomputing. Use a pattern like:

```python
cache_path = Path("data/melnikov_sweep_kappa3.npz")
if cache_path.exists():
    data = np.load(cache_path)
    F_crit_array = data["F_crit"]
else:
    F_crit_array = compute_melnikov_sweep(...)
    np.savez(cache_path, F_crit=F_crit_array, omega=omega_array, ...)
```

- **Markdown cells**: Each notebook should have substantial markdown narrative explaining the mathematical motivation for each computation, not just code. The notebooks double as a lab notebook / analysis record for the project.
- **Interactive widgets**: Where appropriate (especially notebooks 01, 03, 04), use `ipywidgets` sliders for parameter exploration. Example: a slider over $\delta$ that updates the $V(\theta)$ plot in real time, showing the well split as it happens. List `ipywidgets` as an optional dependency.
