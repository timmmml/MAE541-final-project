# Von Mises Drosophila Dual-Goal Dynamics

**Melissa Tian, Tim Lin**
May 6, 2026

## Abstract

Hello world.

---
## 1. Introduction

Full circuit model
$$
\begin{align}
U(\theta) & =f\left( \sum g_{PFL3L}+W_{DNa03}\cdot f\left( \sum g_{PFL3L}+\sum g_{PFL2L} \right) \right) \\
 & \quad -f\left( \sum g_{PFL3R}+W_{DNa03}\cdot f\left( \sum g_{PFL3R}+\sum g_{PFL2R} \right) \right)
\end{align}
$$

## 2. Derivation of the Reduced System
### 1.1 Analytical Form of the Steering Drive 

Taking a reference frame where the two goals, $\theta_{\text{goal}}^{(0)}$ and $\theta_{\text{goal}}^{(1)}$, are located at $\delta$ and $-\delta$ respectively, we can write the von Mises bump representation, under this dual-goal, for a goal neuron with preference $\phi$ as:

$$
g_{\text{goal}}(\phi, \delta) = \exp(\kappa_{\text{goal}} \cos(\delta - \phi)) + \exp(\kappa_{\text{goal}} \cos(-\delta - \phi)) \tag{1}
$$

similarly write down the heading neurons with preference $\phi$ firing to the current heading $\theta$ as:

$$
g_{\text{heading}}(\phi, \theta) = \exp(\kappa_{\text{heading}} \cos(\theta - \phi)) \tag{2}
$$

Considering a PFL cell with goal preference $\phi$, we abstract that each cell receives a direct synaptic connection from the associated goal neuron, as well as a synaptic connection from a heading neuron tuned to $\Delta_{\text{pop}}$ away from $\phi$. Showing the PFL3L subpopulation as an example, we can write down the total drive to a PFL3L neuron with preference $\phi$, $h_{PFL3}$, and the neural activity after the nonlinearity $f$, $g_{PFL3}$, as:

$$
\begin{aligned}
h_{PFL3L}(\theta, \phi, \delta) &= S \cdot g_{\text{heading}}(\phi - \Delta_{PFL3}, \theta) + A \cdot g_{\text{goal}}(\phi, \delta) \\
&= S \cdot \exp(\kappa_{\text{heading}} \cos(\theta - (\phi - \Delta_{PFL3}))) \\
&\quad + A \cdot (\exp(\kappa_{\text{goal}} \cos(\delta - \phi)) + \exp(\kappa_{\text{goal}} \cos(-\delta - \phi))) \\
g_{PFL3L}(\theta, \phi, \delta) &= f(h_{PFL3}(\theta, \phi, \delta))
\end{aligned} \tag{3}
$$

As downstream neurons to the PFLs (DNa02 and DNa03) effectively sum inputs from all neurons of each population (PFL3L/R, PFL2), it behooves us to consider the population activity as a whole, integrating eq. (3) over all $\phi$. We can write down the population activity of the PFL3L subpopulation:

$$
\begin{aligned}
G_{PFL3L}(\theta, \delta) &= \int_{-\pi}^{\pi} g_{PFL3L}(\theta, \phi, \delta)\, d\phi \\
&= \int_{-\pi}^{\pi} f(h_{PFL3L}(\theta, \phi, \delta))\, d\phi
\end{aligned} \tag{4}
$$

We choose the quadratic nonlinearity $f(x) = x^2$, as a common choice in neuroscience and one that allows conjunctive tuning for heading and goal states on the level of PFL3 neurons. This is important to allow differences between the L and R populations after integration.

Hence we have the following:

$$
\begin{aligned}
G_{PFL3L}(\theta, \delta) &= \int_{-\pi}^{\pi} h_{PFL3}(\theta, \phi, \delta)^2 \, d\phi \\
&= \int_{-\pi}^{\pi} (S \cdot g_{\text{heading}}(\phi - \Delta_{PFL3}, \theta) + A \cdot g_{\text{goal}}(\phi, \delta))^2 \, d\phi \\
&= \int_{-\pi}^{\pi} \big[ S^2 \cdot g_{\text{heading}}(\phi - \Delta_{PFL3}, \theta)^2 + A^2 \cdot g_{\text{goal}}(\phi, \delta)^2 \\
&\quad + 2SA \cdot g_{\text{heading}}(\phi - \Delta_{PFL3}, \theta) \cdot g_{\text{goal}}(\phi, \delta) \big] d\phi
\end{aligned} \tag{5}
$$

As $\int_{-\pi}^{\pi} g_{\text{heading}}(\phi - \Delta_{PFL3}, \theta)^2 \, d\phi = \int_{-\pi}^{\pi} g_{\text{heading}}(\phi + \Delta_{PFL3}, \theta)^2 \, d\phi$, we can simplify the difference in population activity between the PFL3L and PFL3R subpopulations as:

$$
\begin{aligned}
G_{PFL3L}(\theta, \delta) - G_{PFL3R}(\theta, \delta) &= 2SA \int_{-\pi}^{\pi} \big( g_{\text{heading}}(\phi - \Delta_{PFL3}, \theta) - g_{\text{heading}}(\phi + \Delta_{PFL3}, \theta) \big) g_{\text{goal}}(\phi, \delta)\, d\phi \\
&= 2SA \int_{-\pi}^{\pi} \big( \exp(\kappa_{\text{heading}} \cos(\theta - (\phi - \Delta_{PFL3}))) \\
&\quad - \exp(\kappa_{\text{heading}} \cos(\theta - (\phi + \Delta_{PFL3}))) \big) \\
&\quad \cdot \big( \exp(\kappa_{\text{goal}} \cos(\delta - \phi)) + \exp(\kappa_{\text{goal}} \cos(-\delta - \phi)) \big) d\phi
\end{aligned} \tag{6}
$$

Let $\xi_\phi(a, b): S^1 \times S^1 \to \mathbb{R} = \exp(\kappa_{\text{heading}} \cos(a - \phi) + \kappa_{\text{goal}} \cos(b - \phi))$, we can rewrite the above as:

$$
\begin{aligned}
G_{PFL3L}(\theta, \delta) - G_{PFL3R}(\theta, \delta) = 2SA \int_{-\pi}^{\pi} \big( &\xi_\phi(\theta + \Delta_{PFL3}, \delta) + \xi_\phi(\theta + \Delta_{PFL3}, -\delta) \\
- &\xi_\phi(\theta - \Delta_{PFL3}, \delta) - \xi_\phi(\theta - \Delta_{PFL3}, -\delta) \big) d\phi
\end{aligned} \tag{7}
$$

Noting that:

$$
\kappa_{\text{heading}} \cos(\phi - a) + \kappa_{\text{goal}} \cos(\phi - b) = \kappa(a, b) \cos(\phi - \mu(a,b))
$$

where the effective concentration is

$$
\kappa(a, b) = \sqrt{\kappa_{\text{heading}}^2 + \kappa_{\text{goal}}^2 + 2\kappa_{\text{heading}}\kappa_{\text{goal}} \cos(a - b)}
$$

and the effective phase is

$$
\mu(a, b) = \operatorname{atan2}\big(\kappa_{\text{heading}}\sin a + \kappa_{\text{goal}}\sin b,\; \kappa_{\text{heading}}\cos a + \kappa_{\text{goal}}\cos b\big).
$$

This is the standard phasor sum: two cosines with the same argument variable $\phi$ but different phases $a, b$ combine into a single cosine with concentration $\kappa(a,b)$ and phase $\mu(a,b)$. Note that $\mu$ depends on $a, b$ but, crucially, drops out after integration over the full period in $\phi$ (see below).

We can evaluate each piece of the integral as:

$$
\begin{aligned}
\int_{-\pi}^{\pi} \xi_\phi(a, b)\, d\phi &= \int_{-\pi}^{\pi} \exp(\kappa_{\text{heading}} \cos(a - \phi) + \kappa_{\text{goal}} \cos(b - \phi))\, d\phi \\
&= \int_{-\pi}^{\pi} \exp(\kappa(a, b) \cos(\phi - \mu(a,b)))\, d\phi \\
&= 2\pi I_0(\kappa(a, b))
\end{aligned} \tag{8}
$$

The last equality uses the standard identity $\int_{-\pi}^{\pi} \exp(\kappa \cos(\phi - \mu))\, d\phi = 2\pi I_0(\kappa)$, which is independent of $\mu$ by translation invariance of the integration measure on $S^1$. Here $I_0$ is the modified Bessel function of the first kind, order 0.

Hence we have the following closed form solution for the difference in population activity between the PFL3L and PFL3R subpopulations:

$$
\begin{aligned}
G_{PFL3L}(\theta, \delta) - G_{PFL3R}(\theta, \delta) = 4\pi SA \big( &I_0(\kappa(\theta + \Delta_{PFL3}, \delta)) + I_0(\kappa(\theta + \Delta_{PFL3}, -\delta)) \\
- &I_0(\kappa(\theta - \Delta_{PFL3}, \delta)) - I_0(\kappa(\theta - \Delta_{PFL3}, -\delta)) \big)
\end{aligned} \tag{9}
$$

Now, in the simplified system where there is no PFL2 input and the DNa02s and DNa03s simply sum the PFL3L and PFL3R population activity respectively, we can write down the steering drive as:

$$
U(\theta, \delta) = G_{PFL3}(\theta, \delta) - G_{PFL3R}(\theta, \delta) \tag{10}
$$

This form of $U$ will be referred to as the **Bessel model** in the following text and figures. 

Numerical calculation from the Bessel model, the discrete PFL3 model with 12 neurons on each side (corresponding to the biological circuit), and the full circuit model (with PFL2 neurons considered) yield very similarly shaped $U(\theta)$, supporting the validity of our simplifications. For the following analysis, we will continue with the PFL3 only circuit, using the 12 neuron discrete model as the ground-truth comparison.

![[Pasted image 20260510184649.png]]
**Figure 1.** Steering drive $U(\theta)$ from different models.


### 1.2 Taylor Expansion Approximation of the Steering Drive
#### 1.2.1 Symmetry of $U(\theta, \delta)$

Before computing the Taylor expansion, we record two symmetry properties that will simplify the derivatives at $\theta = 0$.

**Property (a): $\kappa(a, b)$ depends on $a, b$ only through $\cos(a - b)$.** Direct from the definition of $\kappa(a,b)$.

**Property (b): $U$ is odd in $\theta$ at fixed $\delta$.** Apply the substitution $\theta \to -\theta$ to the four arguments in eq. (9):

| Original                                  | After $\theta \to -\theta$                 |
| ----------------------------------------- | ------------------------------------------ |
| $\kappa(\theta + \Delta_{PFL3}, \delta)$  | $\kappa(-\theta + \Delta_{PFL3}, \delta)$  |
| $\kappa(\theta + \Delta_{PFL3}, -\delta)$ | $\kappa(-\theta + \Delta_{PFL3}, -\delta)$ |
| $\kappa(\theta - \Delta_{PFL3}, \delta)$  | $\kappa(-\theta - \Delta_{PFL3}, \delta)$  |
| $\kappa(\theta - \Delta_{PFL3}, -\delta)$ | $\kappa(-\theta - \Delta_{PFL3}, -\delta)$ |

Using property (a) — that $\kappa$ depends on the *difference* of its arguments via cosine — we have $\kappa(-\theta + \Delta_{PFL3}, \delta) = \kappa(\theta - \Delta_{PFL3}, -\delta)$ (both have $\cos$ of $-\theta + \Delta_{PFL3} - \delta$ inside). Working through all four, the four terms after $\theta \to -\theta$ are precisely the four original terms with their signs swapped. Therefore

$$
U(-\theta, \delta) = -U(\theta, \delta),
$$

i.e. $U$ is odd in $\theta$. As corollaries:

- $U(0, \delta) = 0$ for all $\delta$.
- All even-order derivatives vanish at $\theta = 0$: $U^{(2k)}(0, \delta) = 0$.

This is the systematic version of "the function is odd around $\theta = 0$" used below.

#### 1.2.2 Taylor expansion of $U(\theta, \delta)$ around $\theta = 0$

We expand $U(\theta, \delta)$ around $\theta = 0$ for our dreamed Duffing connection:

$$
U(\theta, \delta) = U(0, \delta) + U'(0, \delta)\,\theta + \frac{1}{2} U''(0, \delta)\,\theta^2 + \frac{1}{6} U'''(0, \delta)\,\theta^3 + O(\theta^4) \tag{11}
$$

By the symmetry result of §1.1:

$$
U(0, \delta) = 0, \qquad U''(0, \delta) = 0. \tag{12}
$$

Only the odd-order terms survive at this order, giving the Duffing-form

$$
U(\theta, \delta) \approx U'(0, \delta)\,\theta + \frac{1}{6} U'''(0, \delta)\,\theta^3 + O(\theta^5).
$$

#### 1.2.3 First-order coefficient

Define $K(\theta; a, b) \equiv \kappa(\theta + a, b)$. By the chain rule:

$$
\frac{d}{d\theta} I_0(K(\theta; a, b)) = I_0'(K) \cdot \frac{dK}{d\theta} = I_1(K) \cdot \frac{dK}{d\theta}, \tag{13a}
$$

using the identity $I_0'(x) = I_1(x)$. The derivative of $K$ is computed from $K^2 = \kappa_h^2 + \kappa_g^2 + 2\kappa_h\kappa_g\cos(\theta + a - b)$, giving $2K K' = -2\kappa_h\kappa_g\sin(\theta + a - b)$, hence

$$
\frac{dK}{d\theta} = -\frac{\kappa_h\kappa_g \sin(\theta + a - b)}{K(\theta; a, b)}. \tag{13b}
$$

Combining:

$$
\frac{d}{d\theta} I_0(K(\theta; a, b)) = -I_1(K) \cdot \frac{\kappa_h\kappa_g \sin(\theta + a - b)}{K}. \tag{13c}
$$

Applying this to each of the four terms in eq. (9) and evaluating at $\theta = 0$:

$$
\begin{aligned}
U'(0, \delta) = -4\pi SA\, \kappa_h\kappa_g \bigg[ &\frac{I_1(\kappa(\Delta_{PFL3}, \delta))}{\kappa(\Delta_{PFL3}, \delta)} \sin(\Delta_{PFL3} - \delta) \\
+\, &\frac{I_1(\kappa(\Delta_{PFL3}, -\delta))}{\kappa(\Delta_{PFL3}, -\delta)} \sin(\Delta_{PFL3} + \delta) \\
-\, &\frac{I_1(\kappa(-\Delta_{PFL3}, \delta))}{\kappa(-\Delta_{PFL3}, \delta)} \sin(-\Delta_{PFL3} - \delta) \\
-\, &\frac{I_1(\kappa(-\Delta_{PFL3}, -\delta))}{\kappa(-\Delta_{PFL3}, -\delta)} \sin(-\Delta_{PFL3} + \delta) \bigg].
\end{aligned} \tag{14}
$$

Using $\kappa(-\Delta_{PFL3}, \delta) = \kappa(\Delta_{PFL3}, -\delta)$ and $\kappa(-\Delta_{PFL3}, -\delta) = \kappa(\Delta_{PFL3}, \delta)$ (property (a)), and $\sin(-x) = -\sin(x)$, the four terms collapse pairwise:

$$
\boxed{\;U'(0, \delta) = -8\pi SA\, \kappa_h\kappa_g \left[ \frac{I_1(\kappa(\Delta_{PFL3}, \delta))}{\kappa(\Delta_{PFL3}, \delta)} \sin(\Delta_{PFL3} - \delta) + \frac{I_1(\kappa(\Delta_{PFL3}, -\delta))}{\kappa(\Delta_{PFL3}, -\delta)} \sin(\Delta_{PFL3} + \delta) \right]\;}
$$

#### 1.2.4 Third-order coefficient

We need $\frac{d^3}{d\theta^3} I_0(K(\theta; a, b))$ at $\theta = 0, with derivatives $K', K'', K'''$ taken with respect to $\theta$. From $K^2 = \kappa_h^2 + \kappa_g^2 + 2\kappa_h\kappa_g\cos(\theta + a - b)$, define

$$
C(\theta) \equiv \kappa_h\kappa_g\cos(\theta + a - b), \qquad S(\theta) \equiv \kappa_h\kappa_g\sin(\theta + a - b).
$$

Then $K^2 = \kappa_h^2 + \kappa_g^2 + 2C$, and successive differentiation gives:

$$
\begin{aligned}
K' &= -\frac{S}{K} \\
K'' &= -\frac{C}{K} - \frac{(K')^2}{K} = -\frac{C}{K} - \frac{S^2}{K^3} \\
K''' &= \frac{d}{d\theta}\left[ -\frac{C}{K} \right]+\frac{d}{d\theta}\left[ -\frac{S^{2}}{K^{3}} \right]=\frac{S}{K}-\frac{CS}{K^{3}}-\frac{2SC}{K^{3}}-\frac{3S^{3}}{K^{5}}=\frac{S}{K}-\frac{3SC}{K^{3}}-\frac{3S^{3}}{K^{5}}
\end{aligned}
$$

For the Bessel side, use the recurrences $I_0' = I_1$, $I_1' = \tfrac{1}{2}(I_0 + I_2)$, $I_1'' = \tfrac{1}{4}(3I_1 + I_3)$.

Applying the chain rule three times:

$$
\begin{align}
I_{0}' & =I_{1}\\
I_{0}'' & =I_{1}'=\frac{1}{2}[I_{0}(x)+I_{2}(x)] \\
 I_{0}'''& =I_{1}''=\frac{1}{2}[I_{0}(x)+I_{2}(x)]'=\frac{1}{2}\left[ I_{1}(x)+\frac{1}{2}(I_{1})(x)+I_{3}(x) \right]=\frac{1}{4}[3I_{1}(x)+I_{3}(x)] \\
\end{align}
$$
Taken together, we have

$$
\begin{align}
\frac{d^3}{d\theta^3} I_0(K) & = I_0'''(K) (K')^3 + 3 I_0''(K) K' K'' + I_0'(K) K'''  \\
 & = \frac{1}{4}[3I_{1}(K)+I_{3}(K)] \cdot\left( -\frac{S^{3}}{K^{3}} \right)\\
  & \quad +​ \underbrace{3 \cdot \frac{1}{2}[I_0(K) + I_2(K)] \cdot \left(\frac{SC}{K^2} + \frac{S^3}{K^4}\right)}_{\text{from } 3 I_0''(K) K' K''} \\
   & \quad +\underbrace{I_1(K) \cdot \left(\frac{S}{K} - \frac{3CS}{K^3} - \frac{3 S^3}{K^5}\right)}_{\text{from } I_0'(K) K'''}​
\end{align}
$$
With $\frac{d^3}{d\theta^3} I_0(K)$ in hand, the third-order coefficient of $U$ is the sum over the four arguments of eq. (9):

$$
U'''(0, \delta) = 4\pi SA \sum_{(a, b) \in \mathcal{P}} \sigma(a) \cdot \frac{d^3}{d\theta^3} I_0(K(\theta; a, b)) \bigg|_{\theta = 0} \tag{16}
$$

where $\mathcal{P} = \{(\Delta_{PFL3}, \delta), (\Delta_{PFL3}, -\delta), (-\Delta_{PFL3}, \delta), (-\Delta_{PFL3}, -\delta)\}$ and $\sigma(a) = +1$ for $a = +\Delta_{PFL3}$, $-1$ for $a = -\Delta_{PFL3}$ (matching the signs in eq. 9).

By property (a), the four-term sum simplifies to two terms after pairing, which gives us

$$
\boxed{U'''(0, \delta) = 8\pi SA\kappa_{h}\kappa_{g}[\sin(\Delta_{PFL3}-\delta)\cdot D(\kappa_{-},C_{-},S_{-})+\sin (\Delta_{P_{3}L}+\delta)\cdot D(\kappa_{+},C_{+},S_{+})]} \tag{17}
$$

with
$$
D(K,C,S)=\frac{1}{K}​\left[ I_{1}​(K)\left( 1-\frac{3C}{K^{2}}​−\frac{3S^{2}}{K^{4}} \right)+\frac{3(I_{0}(K)+I_{2}​(K))}{2}​\left( \frac{C}{K} + \frac{S^{2}}{K^{3}}​ \right)−\frac{3I_{1}(K)+I_{3}(K)}{4}\cdot \frac{S^{2}}{K^{2}} \right]
$$
and
$$
\begin{align}
\kappa_- &\equiv \kappa(\Delta_{PFL3}, \delta) = \sqrt{\kappa_h^2 + \kappa_g^2 + 2\kappa_h\kappa_g\cos(\Delta_{PFL3} - \delta)} \\
\kappa_+ &\equiv \kappa(\Delta_{PFL3}, -\delta) = \sqrt{\kappa_h^2 + \kappa_g^2 + 2\kappa_h\kappa_g\cos(\Delta_{PFL3} + \delta)} \\
C_- &\equiv \kappa_h\kappa_g\cos(\Delta_{PFL3} - \delta), \qquad C_+ \equiv \kappa_h\kappa_g\cos(\Delta_{PFL3} + \delta) \\
S_- &\equiv \kappa_h\kappa_g\sin(\Delta_{PFL3} - \delta), \qquad S_+ \equiv \kappa_h\kappa_g\sin(\Delta_{PFL3} + \delta)
\end{align}
$$

#### 1.2.5 Local cubic expansion of the steering drive

Combining §1.3–§1.5, the local expansion of the steering drive is

$$
\dot{\theta}=U(\theta, \delta) \approx A(\delta)\,\theta - B(\delta)\,\theta^3 + O(\theta^5),
$$
with
$$
A(\delta) \equiv U'(0, \delta), \qquad B(\delta) \equiv -\frac{1}{6} U'''(0, \delta).
$$
in order the match the conventional expression of the Duffing equation.

### 1.3 Body Mechanics and Second-Order Dynamics

Consider the fly as a rigid body with rotational inertia $I$ and rotational drag $\gamma$. 

The reduction up to §1.2.5 yields a **first-order** equation, in which the circuit's output is interpreted (following Westeinde et al.) as a direct angular-velocity command. This is the natural overdamped limit of the fly's heading dynamics: velocity is directly proportional to the steering signal with no inertial transients.

If we restore the inertia term and instead interpret the steering signal $U(\theta,\delta)$ as a **torque command** rather than a angular velocity command, Newton's second law for then gives

$$
I \ddot{\theta} + \gamma \dot{\theta} = U(\theta,\delta)
$$

This contains the first-order model as the overdamped limit $I \to 0$: when inertia is negligible relative to drag, $\ddot\theta$ drops out and one recovers $\gamma\dot\theta = U$.

We note that this is a modeling reinterpretation, not a derivation: the same $U(\theta, \delta)$ that previously carried units of velocity now carries units of torque. The constants $I, \gamma$ enter as phenomenological parameters describing the body dynamics rather than being computed from the circuit itself. An alternative interpretation is that $I$ represents an *effective* inertia arising from circuit and motor lags.

Substituting the local cubic expansion of §1.2.5,

$$

I\ddot\theta + c\dot\theta = A(\delta, \kappa)\,\theta - B(\delta, \kappa)\,\theta^3,

\tag{19}

$$

which is the canonical **unforced double-well Duffing oscillator**.


## 3. Bifurcation Analysis
### 2.1 Supercritical Pitchfork bifurcation and double-well criterion

The first-order system that we directly derived from the circuit model shares **equilibrium structure** with the unforced double-well Duffing oscillator $\ddot{\theta}+\gamma \dot{\theta}=A\theta-B\theta^{3}$. In both systems, equilibria are reached at $\theta=0$ and $\theta=\pm \sqrt{ \frac{A}{B} }$ if $\frac{A}{B}>0$.

Linearize at the origin, $\dot{\theta}=A\theta$, so the fixed point at $\theta = 0$ is **unstable** (i.e., a double-well structure exists in the first-order system) $\text{iff.}\, A(\delta,\kappa) = U'(0, \delta) > 0$. Substituting from eq. (14), the bifurcation curve $A(\delta) = 0$ is given by

$$
-8\pi SA\, \kappa_h\kappa_g \left[ \frac{I_1(\kappa(\Delta_{PFL3}, \delta))}{\kappa(\Delta_{PFL3}, \delta)} \sin(\Delta_{PFL3} - \delta) + \frac{I_1(\kappa(\Delta_{PFL3}, -\delta))}{\kappa(\Delta_{PFL3}, -\delta)} \sin(\Delta_{PFL3} + \delta) \right] =0
$$

When in addition $B(\delta,\kappa)>0$ (i.e. $U'''(0,\delta)<0$), two new equilibria emerge at $\theta=\pm \sqrt{ \frac{A}{B} }$​, giving the double-well potential structure. Substituting from eq. (17), the condition in which the pitchfork bifurcation is supercritical is given by 

$$
8\pi SA\kappa_{h}\kappa_{g}[\sin(\Delta_{PFL3}-\delta)\cdot D(\kappa_{-},C_{-},S_{-})+\sin (\Delta_{P_{3}L}+\delta)\cdot D(\kappa_{+},C_{+},S_{+})]<0
$$

These are transcendental equations in $(\delta, \kappa_h, \kappa_g)$ at fixed $\Delta_{PFL3} = 67.5°$, and setting $\kappa_h = \kappa_g \equiv \kappa$ collapses the parameter space to $(\delta, \kappa)$ for presentation. Numerical solutions yield the following bifurcation diagram.

![[fig1_bifurcation2.png]]
**Figure 2.** **Pitchfork bifurcation in $(\kappa,\delta)$ parameter space.** **Left**, the black curve marks the pitchfork bifurcation curve, where the origin changes stability. Above the curve, $A>0$ and the origin is unstable; below the curve, $A<0$ and the origin is the unique stable fixed point. Green shading indicates where the cubic coefficient $B(\delta,\kappa)>0$, corresponding to a supercritical pitchfork; red shading indicates $B(\delta,\kappa)<0$. **Right**, the exact steering drive $U(\theta)$ (red, from the Bessel model) with its local Duffing truncation (blue dashed) at two sample parameter points.

We can see that the bifurcation boundary lies in the area where $B(\delta,\kappa)>0$ (green shadow), indicating that the bifurcation is a **supercritical pitchfork** everywhere. This implies the transition from monostable to bistable is smooth — no jumps or hysteresis as parameters cross the boundary.

The top right panel shows a sample parameter point from the green area above the curve, which is the double-well regime that we are interested in. The approximation is good near $\theta=0$, capturing the saddle and two flanking stable equilibria. At large $\theta$, it diverges from the exact model, where the Bessel curve has additional zero-crossings near $\theta=\pm \pi$ that the Duffing truncation misses entirely. These correspond to fixed points associated with the anti-midpoint of the two goals, which lie outside the local domain of validity of the Taylor expansion. The Duffing reduction therefore is reliable for **small-amplitude dynamics**, but global features of the $S^{1}$ phase space require the exact steering drive.

## 4. Chaotic Behavior Under Periodic Forcing

### 4.1 Phase portrait and homoclinic orbit


## 5. Discussion

It's worth noting that the full circuit model also yield surprisingly similar $U(\theta)$, diminishing the necessity of having PFL2 in the circuit, which is not the case according to experimental observations. This inconsistency might result from our choice of quadratic non-linearity as opposed to ELU in the original paper. The asymmetric property of ELU might act as a gating mechanism to select a subpopulation of active control neurons, which cannot be captured by our model.


## Notation and conventions

For ease of cross-reference:

- $\theta$: fly's heading.
- $\delta$: half-separation of the two goals (goals at $\pm\delta$).
- $\phi$: preferred direction of a PFL/goal/heading neuron (integration variable).
- $\Delta_{PFL3}$: PFL3L heading-input phase shift ($+67.5°$ for L, $-67.5°$ for R, by convention here).
- $\kappa_{\text{heading}}, \kappa_{\text{goal}}$: von Mises concentrations of heading and goal bumps. Often taken equal to $\kappa$.
- $\kappa(a, b) = \sqrt{\kappa_h^2 + \kappa_g^2 + 2\kappa_h\kappa_g\cos(a-b)}$: effective concentration after phasor sum.
- $\mu(a, b)$: effective phase (drops out after $\phi$ integration).
- $S, A$: synaptic gain on heading and goal inputs (both positive, set by biology).
- $f(x) = x^2$: quadratic nonlinearity (conjunctive tuning).
- $U(\theta, \delta)$: total steering drive, $G_{PFL3} - G_{PFL3R}$.

