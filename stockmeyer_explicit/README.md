**Stockmeyer PDE - Explicit integration (stockmeyer_explicit)**

Obtains the time volution of the Stockmeyer equation reduced to 3 dimensions ($y$, $\theta$, and $\omega$):
$$
\partial f_\tau  = -\omega\partial f_\theta - (\nu_3 \mathcal{G}(f)-\nu_4\vec{E})\cdot n^\bot(\theta)\partial_\omega f
$$

Uses an explicit Euler time integration scheme *without* stabilization. The numerical solution eventually becomes unstable.

To include only the advection term, $-\omega\partial f_\theta$, set Ex0 = 0 and Ey0 = 0 in  ```parameters.prm```.