![PNS](https://raw.githubusercontent.com/pns-mam/mi1/master/logo-pns.png)

## MAM4

# EDP1 
# 2024-25
# Exam CC no. 1

**Durée 2H00. Documents non autorisés. Tous les exercices sont indépendants.
Le barème prévisionnel est indiqué pour chaque exercice.**

## Exercice 1 (6 points)

On considère l'équation de la chaleur sur le domaine $\Omega = ]0,1[$ avec conditions de Dirichlet :

$$ \begin{aligned}
  & \frac{\partial u}{\partial t}(t, x) - \frac{\partial^2 u}{\partial x^2}(t, x) = 0,\quad t>0,\quad x \in \Omega,\\
  & u(t, 0) = u(t, 1) = 0,\quad t > 0,\\
  & u(0, x) = u_0(x),\quad x \in \Omega.
\end{aligned} $$

### 1.1

On définit l'énergie

$$ E(t) := \frac{1}{2} \int_0^1 u^2(t,x)\,\mathrm{d}x. $$

En supposant $u$ solution aussi régulière que nécessaire, montrer que

$$ E'(t) =- \int_0^1 \left( \frac{\partial u}{\partial x}(t,x) \right)^2\,\mathrm{d}x. $$

### 1.2

En déduire que $E$ est décroissante et qu'on a unicité de solution pour l'EDP considérée. 

## Exercice 2 (5 points)

Montrer que le schéma implicite centré ci-dessous

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + V \frac{u_{j+1}^{n+1}-u_{j-1}^{n+1}}{2\Delta x} = 0 $$

est consistant  pour l'équation d'advection avec $V > 0$ (on précisera les ordres en $t$ en $x$, par exemple en faisant un développement en $(t,x) = (t_{n+1},x_j)$) :

$$ \frac{\partial u}{\partial t}(t, x) + V\frac{\partial u}{\partial x}(t, x) = 0. $$

## Exercice 3 (4 points)

Montrer que le schéma décentré aval ci-dessous, avec $V < 0$ tel que $\sigma := |V|\Delta t/\Delta x \leq 1$,

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + V \frac{u_{j+1}^n-u_j^n}{\Delta x} = 0 $$

vérifie le principe du maximum discret : si $u^0_j \in [m, M]$ pour tout $j$, alors $u^n_j \in [m, M]$ pour tout $j$ et tout $n > 0$.

Que peut-on en déduire ?

## Exercice 4 (5 points)

Étudier la stabilité au sens $L^2$ du schéma implicite suivant ($\nu > 0$) :

$$ \frac{2u_j^{n+1} - 3u_j^n + u_j^{n-1}}{\Delta t} - \nu \frac{u_{j+1}^{n+1} - 2u_j^{n+1} + u_{j-1}^{n+1}}{\Delta x^2} = 0. $$
