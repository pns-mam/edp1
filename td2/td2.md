![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1
# 2025-26
# TD 2 - Études de consistance


Considérons l'équation de la chaleur en une dimension d'espace dans le domaine borné $\Omega = (0,1)$ :

$$ \begin{cases}
\displaystyle\frac{\partial u}{\partial t}(t,x) - \nu \frac{\partial^2 u}{\partial x^2}(t,x) = 0,\quad t > 0,\quad x \in (0,1)\\
u(x,0)=u_0(x),\quad x\in (0,1).
\end{cases} $$

Les conditions aux limites sont des conditions de Dirichlet homogène : 

$$ u(t,0)=u(t,1)=0,\quad t > 0. $$

On discrétise le domaine en utilisant un maillage régulier

$$ (t_n,x_j) := (n\Delta t,j\Delta x),\quad n \ge 0,\quad j = 0,\dots,N, $$

avec $\Delta t>0$ et $\Delta x := 1/N$.

## Exercice 1

Étudier la consistance du $\theta$-schéma, $\theta \in [0,1]$ :

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t}-\theta\nu \frac{u_{j+1}^{n}-2u_j^n+u_{j-1}^{n}}{\Delta x^2}-(1-\theta)\nu \frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2}=0. $$

En déduire la consistance des schémas d'Euler explicite, Euler implicite, et Crank-Nicolson.

## Exercice 2

Étudier la consistance du schéma de Gear :

$$ \frac{3u_j^{n+1}-4u_j^n+u_j^{n-1}}{2\Delta t}-\nu \frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2}=0. $$

## Exercice 3

Étudier la consistance du schéma de Du Fort-Frankel : 

$$ \frac{u_j^{n+1}-u_j^{n-1}}{2\Delta t}-\nu \frac{u_{j+1}^{n}-u_j^{n+1}-u_j^{n-1}+u_{j-1}^{n}}{\Delta x^2}=0. $$
