I would like to develop a code base for the final project (as specced in the proposal). Here is a list of analyses I would like to run and the associated simulations i have in mind. 


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

we have access to the full system as well as the reduced system (in the sense of the functional form $U(\theta)$, where the reduced system is literally Duffing with $\alpha, \beta$ as functions of $\kappa$ and $\delta$