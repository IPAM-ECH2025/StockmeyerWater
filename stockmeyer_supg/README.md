**Stockmeyer PDE - Implicit integration with SUPG stabilization (stockmeyer_supg)**

Obtains the time volution of the Stockmeyer equation reduced to 3 dimensions ($y$, $\theta$, and $\omega$):
$$
\partial f_\tau  = -\omega\partial f_\theta - (\nu_3 \mathcal{G}(f)-\nu_4\vec{E})\cdot n^\bot(\theta)\partial_\omega f. \tag{1}
$$

We can write Eq. $(1)$ as
$$
\partial f_\tau  = \vec{u}\cdot\nabla f, \tag{2}
$$
where
$$
\begin{aligned}
\vec{u} &= \left[ u_y,u_{\theta}, u_\omega \right] \\
&= \left [0, -\omega, - (\nu_3 \mathcal{G}(f)-\nu_4\vec{E})\cdot n^\bot(\theta) \right].
\end{aligned}
$$

The SUPG stabilization term is included by replacing the test funcion for the weak form, $\phi$, with 
$$
\phi_S=\phi +\tau_{S}\vec{u}\cdot\nabla f, \tag{3}
$$ 
where
$$
\nabla f=(\partial_y f, \partial_\theta f, \partial_\omega f).
$$
For a rectangular finite element mesh with isotropic elements,  $\tau_S$ in Eq. (3) is given by
$$
\tau_S=\frac{\sqrt{3}h}{2(\xi+1)|\vec{u}|}, \tag{4}
$$ 
where $h$ is the linear element length and $\xi$ is the element degree. 

Currently, only the first term of the right-hand side is implemented:
$$
\partial f_\tau  = -\omega\partial f_\theta.
$$
