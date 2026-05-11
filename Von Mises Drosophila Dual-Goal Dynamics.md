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

We have established that the steering circuit, when reduced to its local cubic normal form (eq. 19), is the unforced double-well Duffing oscillator, and that the bistable regime is reached through a supercritical pitchfork in $(\delta, \kappa)$ parameter space. We now ask what happens when this system is driven by periodic forcing — biologically, this models oscillatory perturbations to the steering command (e.g. visual or proprioceptive input modulated by wingbeat-locked sensory streams, or a periodic exogenous torque). The forced system reads

$$
\ddot\theta + \gamma\dot\theta = U(\theta, \delta) + F\cos(\omega t), \tag{20}
$$

where we have set $I = 1$ for notational economy (this fixes a time-scale; all reported frequencies should be read in units of $\sqrt{I^{-1}}$). The forced Duffing oscillator is the canonical example of a low-dimensional smooth system that exhibits Smale-horseshoe chaos for sufficiently strong forcing [4]. Our goal in this section is twofold: (i) to derive an analytical condition on $(F, \omega, \gamma, \delta, \kappa)$ at which chaotic dynamics first becomes possible in the reduced system, via Melnikov's method, and (ii) to verify this condition numerically in both the reduced (Duffing) and Bessel models, and to characterize a second, qualitatively distinct chaotic regime that emerges in the global $S^1$ system but is invisible to the local Duffing reduction.

### 4.1 Phase portrait and homoclinic orbit

The Melnikov analysis is built around the homoclinic orbit of the unforced, undamped reduced system

$$
\ddot\theta = A\theta - B\theta^3, \tag{21}
$$

with $A, B > 0$ (the supercritical regime of §3.1). This is a planar Hamiltonian system with conserved energy

$$
H(\theta, v) = \frac{1}{2}v^2 - \frac{1}{2}A\theta^2 + \frac{1}{4}B\theta^4, \qquad v \equiv \dot\theta. \tag{22}
$$

The level sets of $H$ foliate the phase plane. The origin is a hyperbolic saddle with $H(0,0) = 0$ and unstable/stable eigenvalues $\pm\sqrt{A}$; the two centers at $\theta_\pm = \pm\sqrt{A/B}$ have $H(\theta_\pm, 0) = -A^2/(4B) < 0$ and linearized frequency $\omega_0 = \sqrt{2A}$. The zero-energy level set $H = 0$ consists of two homoclinic loops joining the saddle to itself, one in each half-plane. Solving $\frac{1}{2}\dot\theta^2 = \frac{1}{2}A\theta^2 - \frac{1}{4}B\theta^4$ on the right branch yields the closed-form

$$
\boxed{\;\theta_0(t) = \sqrt{\frac{2A}{B}}\operatorname{sech}(\sqrt{A}\,t), \qquad v_0(t) = -\sqrt{\frac{2A}{B}}\sqrt{A}\,\operatorname{sech}(\sqrt{A}\,t)\tanh(\sqrt{A}\,t).\;} \tag{23}
$$

This orbit reaches its apex $\theta_{\max} = \sqrt{2A/B}$ at $t = 0$ and approaches the saddle exponentially as $t \to \pm\infty$ at rate $\sqrt{A}$. The corresponding loop in the left half-plane is obtained by $\theta_0 \to -\theta_0$; the union forms the figure-eight separatrix that divides the phase plane into the two wells (oscillation around either center) and the region of unbounded growth (excluded in our problem by the global structure of the Bessel $U$, but present in the cubic truncation).

The same construction does **not** carry over to the full Bessel system without modification. The Bessel steering drive $U(\theta, \delta)$ on $S^1$ has additional fixed points near $\theta = \pm\pi$ (the "anti-midpoint" of the two goals; cf. §3.1 and Figure 2), giving rise to a richer phase-cylinder structure with multiple saddles and the possibility of heteroclinic connections between them. We will return to this point in §4.4. For the local Melnikov analysis we adopt a **second approximation**, in addition to the Taylor truncation already made in §1.2.5: we treat the phase cylinder $S^1 \times \mathbb{R}$ as the plane $\mathbb{R}^2$, on which the homoclinic (23) is well-defined. This planar approximation is internally consistent with the cubic truncation, since both are controlled by the same small parameter $\theta_{\max} \ll \pi$.

### 4.2 The Melnikov criterion

Following Guckenheimer & Holmes [4, §4.5], we cast the forced damped system (20) as an $O(\epsilon)$ perturbation of the Hamiltonian system (21):

$$
\dot\theta = v, \qquad \dot v = A\theta - B\theta^3 + \epsilon\bigl[-\gamma v + F\cos\omega t\bigr], \tag{24}
$$

with the understanding that $\gamma$ and $F$ are formally $O(1)$ and the perturbation strength is controlled by $\epsilon \ll 1$. The Melnikov function $M(t_0)$ measures, to leading order in $\epsilon$, the signed distance between the stable and unstable manifolds of the perturbed saddle (now a hyperbolic periodic orbit of period $2\pi/\omega$ in the extended phase space) at a fixed phase $\omega t_0$ of the forcing cycle. Simple zeros of $M(t_0)$ correspond to transverse intersections of these manifolds, and by the Smale-Birkhoff theorem imply the existence of a Smale horseshoe — hence symbolic dynamics, positive topological entropy, and sensitive dependence on initial conditions on the invariant set [4, Theorem 4.5.3].

The Melnikov function for a planar Hamiltonian system $\dot{\mathbf{x}} = J\nabla H + \epsilon\mathbf{g}(\mathbf{x}, t)$ with homoclinic orbit $\boldsymbol{\gamma}_0(t)$ is the integral

$$
M(t_0) = \int_{-\infty}^{\infty} \nabla H(\boldsymbol{\gamma}_0(\tau)) \cdot \mathbf{g}(\boldsymbol{\gamma}_0(\tau), \tau + t_0)\,d\tau, \tag{25}
$$

evaluated along the unperturbed orbit. In our case $\mathbf{g} = (0, -\gamma v + F\cos\omega t)^T$ and $\nabla H = (-A\theta - B\theta^3, v)^T = (-\dot v_0, v_0)$ (the second equality using the unperturbed equation of motion), so $\nabla H \cdot \mathbf{g} = v_0(\tau)\bigl[-\gamma v_0(\tau) + F\cos\omega(\tau + t_0)\bigr]$ and the integral splits into a damping term and a forcing term:

$$
M(t_0) = -\gamma \underbrace{\int_{-\infty}^{\infty} v_0(\tau)^2\,d\tau}_{\equiv\,D} + F\int_{-\infty}^{\infty} v_0(\tau)\cos\omega(\tau + t_0)\,d\tau. \tag{26}
$$

The damping integral $D$ is shift-invariant. Substituting (23) and $u = \sqrt{A}\,\tau$,

$$
D = \frac{2A^2}{B}\int_{-\infty}^{\infty}\operatorname{sech}^2 u\,\tanh^2 u\,du \cdot \frac{1}{\sqrt{A}} = \frac{4A^{3/2}}{3B}, \tag{27}
$$

using the elementary integral $\int_{-\infty}^{\infty}\operatorname{sech}^2 u\,\tanh^2 u\,du = 2/3$. For the forcing integral, expand $\cos\omega(\tau + t_0) = \cos\omega\tau\cos\omega t_0 - \sin\omega\tau\sin\omega t_0$. Since $v_0(\tau)$ is odd in $\tau$, the $\cos\omega\tau$ piece (odd $\times$ even) integrates to zero, leaving

$$
\int_{-\infty}^{\infty} v_0(\tau)\cos\omega(\tau + t_0)\,d\tau = -\sin\omega t_0\int_{-\infty}^{\infty} v_0(\tau)\sin\omega\tau\,d\tau.
$$

The remaining integral is the sine transform of $-\sqrt{2A/B}\cdot\sqrt{A}\operatorname{sech}(\sqrt{A}\,\tau)\tanh(\sqrt{A}\,\tau)$. Using the identity

$$
\int_{-\infty}^{\infty} \operatorname{sech}(\alpha\tau)\tanh(\alpha\tau)\sin\omega\tau\,d\tau = -\frac{\pi\omega}{\alpha^2}\operatorname{sech}\!\left(\frac{\pi\omega}{2\alpha}\right),
$$

which follows from contour integration over the poles of $\operatorname{sech}$ at $\tau = i\pi(n + \frac{1}{2})/\alpha$, with $\alpha = \sqrt{A}$ we obtain

$$
\int_{-\infty}^{\infty} v_0(\tau)\sin\omega\tau\,d\tau = \pi\omega\sqrt{\frac{2}{B}}\,\operatorname{sech}\!\left(\frac{\pi\omega}{2\sqrt{A}}\right). \tag{28}
$$

Combining (26)–(28),

$$
\boxed{\;M(t_0) = -\frac{4\gamma A^{3/2}}{3B} - F\pi\omega\sqrt{\frac{2}{B}}\,\operatorname{sech}\!\left(\frac{\pi\omega}{2\sqrt{A}}\right)\sin\omega t_0.\;} \tag{29}
$$

$M(t_0)$ has simple zeros if and only if the amplitude of the sinusoidal part exceeds the damping bias, which yields the **Melnikov threshold**

$$
\boxed{\;F_{\text{crit}}(\omega) = \frac{4\gamma A}{3\pi\omega\sqrt{2B}}\cosh\!\left(\frac{\pi\omega}{2\sqrt{A}}\right).\;} \tag{30}
$$

Two features of (30) deserve emphasis. First, $F_{\text{crit}} \to \infty$ as $\omega \to 0$ (linear divergence) and as $\omega \to \infty$ (exponential divergence through $\cosh$): both very-low-frequency forcing (which shifts the equilibrium adiabatically rather than exciting the homoclinic) and very-high-frequency forcing (which averages out before the orbit can respond) are inefficient at inducing chaos. The threshold is minimized at an intermediate "homoclinic-resonance" frequency $\omega^* \approx 0.76\sqrt{A}$, obtained from $\tanh(\pi\omega^*/2\sqrt{A}) = 2\sqrt{A}/(\pi\omega^*)$. Second, the dependence on the well structure enters only through the two coefficients $A, B$ — but these in turn depend nontrivially on $(\delta, \kappa, S, A_{\text{syn}})$ via eqs. (14) and (17), so (30) implicitly maps the chaos-onset surface into circuit parameter space.

### 4.3 Numerical verification near the homoclinic-resonance frequency

We fix the parameters $\kappa = 1.0$, $\delta = 1.49$ (well inside the bistable regime of Figure 2), and $\gamma = 0.1$. This yields $A \approx 0.18$, $B \approx 1.31$, giving a saddle eigenvalue $\sqrt{A} \approx 0.42$ and well frequency $\omega_0 = \sqrt{2A} \approx 0.60$. We choose the forcing frequency

$$
\omega = \sqrt{A} \approx 0.42,
$$

which is the convention used in Guckenheimer & Holmes [4] for the canonical twin-well Duffing. This $\omega$ is close to (though not exactly at) the homoclinic-resonance optimum, with $\cosh(\pi\omega/2\sqrt{A}) = \cosh(\pi/2) \approx 2.51$. Substituting into (30) gives $F_{\text{crit}} \approx 0.012$.

We then perform a loose grid search over the forcing amplitude,

$$
F \in \{1.0,\, 1.5,\, 2.0,\, 2.5,\, 3.0\} \times F_{\text{crit}},
$$

integrating (20) for the reduced Duffing model and for the 12-neuron discrete circuit model from initial conditions in the right well ($\theta_0 = \sqrt{A/B}, v_0 = 0$). For each run we discard the first 200 forcing periods as transient and record the next 200 stroboscopic samples $(\theta_n, v_n) = (\theta(t_n), v(t_n))$ at $t_n = 2\pi n/\omega$. We then sweep $F$ over a denser grid to construct the Poincaré bifurcation diagram $\theta_n$ vs $F$.

**Figure 4** presents three rows: (A) reduced Duffing, (B) 12-neuron full circuit at $\delta = 1.49$ (the bistable regime), and (C) 6-neuron full circuit at the same bistable regime, but as previously shown only has a single attractor (hence belonging to the monostable regime). Each row contains four panels: (i) stroboscopic Poincaré sections at the five grid-search values of $F$; (ii) the Poincaré bifurcation diagram on a narrow zoom $F \in [0, 0.1]$ revealing the route to chaos; (iii) full phase-space trajectories at the same five $F$ values; (iv) the bifurcation diagram on a wide zoom $F \in [0, 10]$ exposing the regime of large forcing studied in §4.4.

Across both the reduced Duffing (row A) and the 12-neuron circuit (row B), the same qualitative scenario unfolds. At $F = F_{\text{crit}}$ (column 1 of panels i and iii) the attractor remains a period-1 orbit confined to the right well — the horseshoe predicted by Melnikov has been born, but it is a chaotic *saddle*, not an attractor. By $F \approx 1.5\,F_{\text{crit}}$ (column 2) the periodic orbit has lost stability and the trajectory has filled out a strange attractor that visibly explores both wells, crossing the unperturbed separatrix infinitely often. The bifurcation diagrams in panel (ii) show the characteristic interleaving of chaotic bands and narrow periodic windows in this $F$ range. By $F \approx 2.5\text{–}3\,F_{\text{crit}}$ (columns 4–5) the trajectory has re-collapsed onto a periodic orbit of large amplitude, encircling both centers.

The factor-of-$1.5$ gap between the Melnikov threshold and the visible onset of chaos is consistent with the topological-vs-metric distinction noted above and with values reported by Holmes for the canonical Duffing system at comparable damping and frequency ratios [4, §4.6]. The agreement between the reduced model (A) and the 12-neuron circuit (B) is quantitatively close, with the chaos band of (B) shifted slightly to higher $F$. This validates the Duffing reduction as a faithful predictor of chaos onset in the biological circuit, at least in this parameter regime.

Row (C), where the system sits below the pitchfork curve, provides a useful control: here the unforced system has a single stable fixed point and no homoclinic orbit, so the Melnikov machinery does not apply. Correspondingly, the trajectories in panel (C-iii) remain periodic across the entire grid — no chaos emerges at these forcing amplitudes. The strange attractor of (A) and (B) is genuinely a consequence of the double-well structure produced by the pitchfork bifurcation.

We restrict attention to initial conditions in the right well throughout. By the $\theta \to -\theta$ symmetry of $U$, there exists a mirror-image strange attractor for initial conditions in the left well, with the same Poincaré bifurcation structure under $\theta \to -\theta$. For sufficiently large $F$ the two attractors merge into a single symmetric one (visible in the wide-zoom panel iv), but their structure for $F \sim F_{\text{crit}}$ is independent.

### 4.4 A second chaotic regime: phase-slipping and the pendulum analogy

The wide-zoom bifurcation diagrams in panels (A-iv) and (B-iv) of Figure 4, together with the snapshots in **Figure 5**, reveal a striking divergence between the reduced and full systems at large forcing amplitudes ($F \gtrsim 2$, two orders of magnitude above $F_{\text{crit}}$).

In the reduced Duffing model (Figure 5A), increasing $F$ eventually drives the trajectory out of the chaotic regime into a large-amplitude periodic orbit that encircles both wells. The Poincaré section in this regime is a small number of isolated points (Figure 5A-i, $F = 2.25$ through $3.0$) and the phase-space trajectory (Figure 5A-ii) is a smooth closed curve oscillating between $\theta \approx \pm 3$. This is the "global" periodic attractor familiar from the standard forced Duffing oscillator.

The full circuit (Figure 5B) behaves qualitatively differently in the same parameter range. Rather than settling onto a closed periodic orbit, the trajectory undergoes **phase slipping**: $\theta$ executes net rotations around the cylinder $S^1$ rather than remaining bounded. The Poincaré sections (Figure 5B-i) are wide diffuse clouds extending across the full range of $\theta \in [-\pi, \pi]$ and reaching velocities $|v| \sim 10$, an order of magnitude larger than in the Duffing-type chaos of Figure 4. The phase trajectories (Figure 5B-ii) wrap repeatedly around the cylinder, with the heading drifting persistently in one direction while the velocity executes large excursions.

This phenomenon is invisible to the Duffing reduction because the cubic truncation does not "know" about the global $S^1$ topology of the heading variable. It is, however, naturally understood through analogy with the forced damped pendulum,

$$
\ddot\varphi + \gamma\dot\varphi + \sin\varphi = F\cos\omega t,
$$

which has been studied extensively as a model system for chaos on the phase cylinder [4, §2.2 and §4.8]. The pendulum has two equivalent saddle points at $\varphi = \pm\pi$ (identified on $S^1$) connected by two heteroclinic orbits — the upper and lower branches of the separatrix that bound the regions of librating versus rotating motion. Under periodic forcing, transverse intersections of the stable and unstable manifolds of these heteroclinic connections generate a chaotic invariant set on which the trajectory irregularly switches between librating and rotating behavior, with the rotating component manifesting macroscopically as phase slip [4, §4.8].

Our full Bessel system inherits this pendulum-like structure: the steering drive $U(\theta, \delta)$ vanishes at $\theta = 0$ (the unstable saddle between the two goals, when bistable) and at $\theta = \pm\pi$ (the saddles associated with the anti-midpoint of the two goals; cf. Figure 2 right panel). When the trajectory has enough energy to traverse the global separatrix, it can transit between adjacent saddles, producing phase slipping analogous to the rotating regime of the forced pendulum. The chaotic dynamics of Figure 5B is most naturally described as transverse intersection of the heteroclinic manifolds connecting the $\theta = 0$ and $\theta = \pm\pi$ saddles, rather than the homoclinic intersection that governs Figure 4.

A revealing comparison is provided by row C of Figure 5, which uses the same large-$F$ protocol but where the two goals effectively merge (the monostable regime of Figure 4C). In this configuration the cubic reduction predicts no chaos at any $F$ — there is no double-well structure and no homoclinic orbit. Yet the full system (Figure 5C) still exhibits phase-slipping chaos at large $F$, qualitatively indistinguishable from Figure 5B. This is direct evidence that the second chaotic regime is **not** a consequence of the goal-induced double well, but of the global $S^1$ structure of the steering drive. The reduced Duffing model is silent on this regime precisely because it is local: the cubic truncation captures the saddle at $\theta = 0$ but discards the saddles at $\theta = \pm\pi$.

A detailed Melnikov analysis of phase-slipping chaos would require: (i) characterizing the heteroclinic orbits connecting $\theta = 0$ and $\theta = \pm\pi$ in the Bessel system, which do not admit a simple closed form; (ii) computing the Melnikov integral along each heteroclinic branch using the numerical orbit; (iii) accounting for the codimension of the heteroclinic chain on the cylinder. Each step is in principle straightforward but operationally involved, and the resulting expression for the chaos threshold would depend on $(\delta, \kappa, \gamma, \omega)$ through circuit-specific integrals analogous to (28). We leave this as a direction for future work. For the present, we content ourselves with the numerical demonstration that two distinct chaotic regimes exist in the steering circuit — one captured faithfully by the Duffing reduction near $F_{\text{crit}}$, and a second, global regime accessible only through the full Bessel or neural circuit model.


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

<<<<<<< HEAD
![[Pasted image 20260510223250.png]]
![[Pasted image 20260510223801.png]]
=======
>>>>>>> 73d40a9ee4322a3933c2344f71e95ae81b2e519a
