I would like to develop a code base for the final project (as specced in the proposal). Here is a list of analyses I would like to run and the associated simulations i have in mind. 

![[Pasted image 20260509163401.png]]
## section 1: pitchfork bifurcation

recall the form 
$$
\begin{align}
U(\theta, \delta) =  \dots
\end{align}
$$

we consider two cases: 

$$
\begin{align}
\text{velocity control: } &  \\ 
\dot{\theta} \propto U(\theta) \\ 
\text{acceleration control: } &  \\ 
I \ddot{\theta} + \gamma \dot{\theta} \propto U(\theta) 
\end{align}
$$

A2

we have access to the full system as well as the reduced system (in the sense of the functional form $U(\theta)$, where the reduced system is literally Duffing with $\alpha, \beta$ as functions of $\kappa$, $\delta$, $S$, and $A$. we will probably focus on the situation here holding $S, A constant and plot out the bifurcation diagram in the $\kappa, \delta$ plane. 

we would do analytical derivation of the bifurcation condition and classify the bifurcation type (supercritical/subcritical), with bifurcation theorem. The end product will be a bifurcation diagram in the $\kappa, \delta$ plane, with the bifurcation curve plotted out and the type of bifurcation indicated. We will compare the analytical results (i.e. boundary conditions) with numerical simulations of the full system. 

This will be figure 1, which shows a contour diagram of the bifurcation curve in the $\kappa, \delta$ plane, with the bifurcation type indicated, clearly indicating the differences between the original system and the reduced system. The figure will have two smaller panels showing example $U(\theta)$ in these regimes. 

B1, B2: (these and below will all be acceleration control)
Hamiltonian structure and phase portrait
- for each parameter regime, draw the phase portrait and overlay level sets of the hamiltonian in the reduced system. Use two panels to plot out the phase portrait of the original system and the reduced system. identify the centers.saddles, POs and homoclinic orbits
- we will write down the equation for the homoclinic orbit in the local space (treating $S^{1}$ as $\mathbb{R}^{1}$). 

Then we are going to investigate chaotic dynamics under periodic forcing. we will find conditions $\gamma, \omega, \delta$ for which the melnikov function has simple zeros and argue the existence of the smale's horseshoe. for this, the Melnikov analysis lives locally (we take a second approximation, in that we treat the phase cylinder as $\mathbb{R}^{2}$, to write down the homoclinic orbit). We will then perform numerical simulations for the reduced system in $\mathbb{R}^{2}$ and the reduced system in $S^{1}\times \mathbb{R}^{1}$, and then the full system, to compare the behaviors (in the form of bifurcation diagrams in $\gamma$ or $\omega$ where we plot the asymptotic poincare section values, as well as poincare sections themselves by constructing stroboscopic poincare maps (sampled at $t = \frac{2\pi n}{\omega}$) to see strange atractor vs. period n orbit)  

also plot characteristic stuff like $F_{crit}$ vs $\delta$ at fixed $\omega$ and $\kappa$. 