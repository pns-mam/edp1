![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1
# 2026-27
# Exam CC no. 1 (corrigé)

**Durée 2H00. Documents et objets connectés non autorisés. Tous les exercices sont indépendants.
Le barème prévisionnel est indiqué pour chaque exercice. Rendre sur deux copies séparées les exercices 1 et 2 d'une part, 3 d'autre part.**

## Exercice 1 (6 points)

On considère l'équation de convection-diffusion sur le domaine $\Omega = ]0,1[$ avec conditions aux limites de Dirichlet homogènes :

$$ \begin{aligned}
  & u_t(t, x) + V u_x(t, x) - \nu u_{xx}(t, x) = 0,\quad t>0,\quad x \in \Omega,\\
  & u(t, 0) = u(t, 1) = 0,\quad t > 0,\\
  & u(0, x) = u_0(x), \quad x \in \Omega,
\end{aligned} $$

où $V \in \mathbf{R}$ (vitesse d'advection) et $\nu > 0$ (coefficient de diffusion) sont des constantes. On suppose $u$ solution aussi régulière que nécessaire.

### 1.1

On définit l'énergie

$$ E(t) := \frac{1}{2} \int_0^1 u^2(t,x)\,\mathrm{d}x. $$

Montrer que

$$ E'(t) = -\nu \int_0^1 u_x^2(t,x)\,\mathrm{d}x. $$

**Réponse.** En dérivant sous le signe intégral et en utilisant l'équation :

$$ E'(t) = \int_0^1 u\,u_t\,\mathrm{d}x = -V \int_0^1 u\,u_x\,\mathrm{d}x + \nu \int_0^1 u\,u_{xx}\,\mathrm{d}x. $$

Le terme de convection est nul :

$$ \int_0^1 u\,u_x\,\mathrm{d}x = \left[ \frac{u^2}{2} \right]_0^1 = 0 $$

grâce aux conditions de Dirichlet. Pour le terme de diffusion, une intégration par parties donne

$$ \int_0^1 u\,u_{xx}\,\mathrm{d}x = \left[ u\,u_x \right]_0^1 - \int_0^1 u_x^2\,\mathrm{d}x = - \int_0^1 u_x^2\,\mathrm{d}x, $$

d'où le résultat. La convection transporte l'énergie sans la modifier ; seule la diffusion la dissipe.

### 1.2

En déduire qu'il y a unicité de solution.

**Réponse.** Si $u$ et $v$ sont deux solutions, $w := u - v$ est solution du même problème (linéarité) avec donnée initiale nulle. L'énergie $E_w$ associée à $w$ vérifie $E_w' \le 0$ d'après 1.1, et $E_w(0) = 0$, donc $0 \le E_w(t) \le E_w(0) = 0$ pour tout $t \ge 0$. Ainsi $\int_0^1 w^2(t,x)\,\mathrm{d}x = 0$ et, $w$ étant continue, $w(t, \cdot) = 0$ pour tout $t$, *i.e.* $u = v$.

### 1.3

On admet l'inégalité de Poincaré : pour toute fonction $v$ continûment dérivable sur $[0,1]$ telle que $v(0) = 0$,

$$ \int_0^1 v^2(x)\,\mathrm{d}x \le \int_0^1 |v'(x)|^2\,\mathrm{d}x. $$

Montrer que $E(t) \le E(0)\,e^{-2\nu t}$ pour tout $t \ge 0$.

**Réponse.** Pour $t > 0$ fixé, $v := u(t,\cdot)$ vérifie $v(0) = 0$, donc d'après Poincaré

$$ E'(t) = -\nu \int_0^1 u_x^2\,\mathrm{d}x \le -\nu \int_0^1 u^2\,\mathrm{d}x = -2\nu E(t). $$

Alors $\left( e^{2\nu t} E(t) \right)' = e^{2\nu t}\left( E'(t) + 2\nu E(t) \right) \le 0$, donc $e^{2\nu t}E(t) \le E(0)$.

## Exercice 2 (6 points)

On considère l'équation de convection-diffusion de l'exercice 1 (constantes $V \in \mathbf{R}$, $\nu > 0$), sans se préoccuper des conditions aux limites, et le schéma explicite centré

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + V \frac{u_{j+1}^n-u_{j-1}^n}{2\Delta x} -
   \nu \frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{\Delta x^2} = 0. $$

On note $\lambda := V\Delta t/\Delta x$ et $\sigma := \nu\Delta t/\Delta x^2$.

### 2.1

Montrer que le schéma est consistant, d'ordre $1$ en temps et $2$ en espace.

**Réponse.** Développements de Taylor au point $(t_n, x_j)$ :

$$ \frac{u(t_{n+1},x_j)-u(t_n,x_j)}{\Delta t} = u_t + \frac{\Delta t}{2} u_{tt} + O(\Delta t^2), $$

$$ \frac{u(t_n,x_{j+1})-u(t_n,x_{j-1})}{2\Delta x} = u_x + \frac{\Delta x^2}{6} u_{xxx} + O(\Delta x^4), $$

$$ \frac{u(t_n,x_{j+1})-2u(t_n,x_j)+u(t_n,x_{j-1})}{\Delta x^2} = u_{xx} + \frac{\Delta x^2}{12} u_{xxxx} + O(\Delta x^4). $$

Comme $u_t + Vu_x - \nu u_{xx} = 0$, l'erreur de consistance vaut

$$ \varepsilon_j^n = \frac{\Delta t}{2} u_{tt} + \Delta x^2 \left( \frac{V}{6} u_{xxx} - \frac{\nu}{12} u_{xxxx} \right) + O(\Delta t^2) + O(\Delta x^4). $$

Le schéma est donc consistant d'ordre $1$ en temps et $2$ en espace, le coefficient $\frac{1}{2}u_{tt}$ de $\Delta t$ n'étant pas nul en général. (On peut l'exprimer en espace : $u_{tt} = V^2 u_{xx} - 2V\nu\, u_{xxx} + \nu^2 u_{xxxx}$, mais ce n'est pas nécessaire.)

### 2.2

Étudier la stabilité au sens $L^2$ du schéma.

**Réponse.** On injecte le mode de Fourier $u_j^n = a^n(\xi)\,e^{2i\pi\xi j\Delta x}$ et on pose $\varphi := 2\pi\xi\Delta x$. En divisant par $a^n(\xi)e^{2i\pi\xi j\Delta x}$ :

$$ a(\xi) = 1 - \frac{\lambda}{2}\left( e^{i\varphi} - e^{-i\varphi} \right) + \sigma\left( e^{i\varphi} - 2 + e^{-i\varphi} \right) = 1 - 4\sigma\sin^2(\varphi/2) - i\lambda\sin\varphi, $$

car $1 - \cos\varphi = 2\sin^2(\varphi/2)$. En notant $s := \sin^2(\varphi/2) \in [0,1]$, avec $\sin^2\varphi = 4s(1-s)$ :

$$ |a(\xi)|^2 = (1 - 4\sigma s)^2 + 4\lambda^2 s(1-s) = 1 + 4s\left[ (\lambda^2 - 2\sigma) + s\,(4\sigma^2 - \lambda^2) \right]. $$

Ainsi $|a(\xi)| \le 1$ pour tout $\xi$ ssi la fonction affine $g(s) := (\lambda^2 - 2\sigma) + s\,(4\sigma^2 - \lambda^2)$ est $\le 0$ sur $[0,1]$, *i.e.* ssi $g(0) \le 0$ et $g(1) \le 0$ :

- $g(0) = \lambda^2 - 2\sigma \le 0$, *i.e.* $V^2\Delta t \le 2\nu$ ;
- $g(1) = 4\sigma^2 - 2\sigma = 2\sigma(2\sigma - 1) \le 0$, *i.e.* $\sigma \le 1/2$.

Sous la condition $\lambda^2 \le 2\sigma \le 1$, le schéma est donc stable au sens $L^2$. La première condition porte sur $\Delta t$ seul : elle est contraignante quand la convection domine la diffusion. En particulier, pour $\nu = 0$ (advection pure), $|a(\xi)|^2 = 1 + \lambda^2\sin^2\varphi > 1$ pour $\lambda \neq 0$ et $\sin\varphi \neq 0$ : c'est la diffusion qui stabilise le schéma.

## Exercice 3 (8 points)

On considère l'équation d'advection

$$ u_t(t, x) + V u_x(t, x) = 0, $$

avec $V \in \mathbf{R}$, $V \neq 0$, et le schéma *saute-mouton* (Leapfrog), centré en temps et en espace :

$$ \frac{u_j^{n+1}-u_j^{n-1}}{2\Delta t} + V \frac{u_{j+1}^n-u_{j-1}^n}{2\Delta x} = 0. $$

On note $\lambda := V\Delta t/\Delta x$.

### 3.1

Montrer que le schéma est consistant d'ordre $2$ en temps et en espace.

**Réponse.** Développements de Taylor au point $(t_n, x_j)$ (les termes pairs se compensent par symétrie) :

$$ \frac{u(t_{n+1},x_j)-u(t_{n-1},x_j)}{2\Delta t} = u_t + \frac{\Delta t^2}{6} u_{ttt} + O(\Delta t^4), $$

$$ \frac{u(t_n,x_{j+1})-u(t_n,x_{j-1})}{2\Delta x} = u_x + \frac{\Delta x^2}{6} u_{xxx} + O(\Delta x^4). $$

Comme $u_t + Vu_x = 0$, l'erreur de consistance vaut

$$ \varepsilon_j^n = \frac{\Delta t^2}{6} u_{ttt} + \frac{V\Delta x^2}{6} u_{xxx} + O(\Delta t^4) + O(\Delta x^4). $$

Le schéma est consistant d'ordre $2$ en temps et en espace, les coefficients de $\Delta t^2$ et $\Delta x^2$ n'étant pas nuls en général. (Comme $u_t = -V u_x$ avec $V$ constant, on a $u_{tt} = -V u_{xt} = V^2 u_{xx}$ puis $u_{ttt} = -V^3 u_{xxx}$, et $\varepsilon_j^n = \frac{V}{6}\left( \Delta x^2 - V^2\Delta t^2 \right) u_{xxx} + O(\Delta t^4) + O(\Delta x^4)$.)

### 3.2

Montrer que le facteur d'amplification est solution d'une équation de degré deux que l'on précisera.

**Réponse.** En injectant le mode de Fourier $u_j^n = a^n(\xi)\,e^{2i\pi\xi j\Delta x}$ et en divisant par $a^{n-1}(\xi)e^{2i\pi\xi j\Delta x}$ :

$$ \frac{a^2 - 1}{2\Delta t} + V\, a\, \frac{e^{i\varphi} - e^{-i\varphi}}{2\Delta x} = 0 \iff a^2 + 2ia\lambda\sin\varphi - 1 = 0, \quad \varphi := 2\pi\xi\Delta x. $$

### 3.3

Discuter la stabilité au sens $L^2$ en fonction de $|\lambda|$.

**Réponse.** Le discriminant réduit vaut $\delta := 1 - \lambda^2\sin^2\varphi$ et les racines sont $a_{1,2}(\xi) = -i\lambda\sin\varphi \pm \sqrt{\delta}$ ; leur produit vaut $a_1 a_2 = -1$.

- Si $|\lambda| < 1$ : $\delta > 0$ pour tout $\xi$, les racines sont distinctes et $|a_{1,2}|^2 = \lambda^2\sin^2\varphi + \delta = 1$. La condition $|a_{1,2}(\xi)| \le 1$ est satisfaite pour tout $\xi$ : le schéma est stable au sens $L^2$.
- Si $|\lambda| = 1$ : pour $\varphi = \pi/2$ (*i.e.* $\xi = 1/(4\Delta x)$), $\delta = 0$ et $a_1 = a_2 = -i\lambda$ est racine double de module $1$. La suite associée contient alors un terme de la forme $n(-i\lambda)^n$, de module $n$, qui croît linéairement : la condition (racine double de module $< 1$) n'est pas satisfaite.
- Si $|\lambda| > 1$ : il existe $\xi$ tel que $\lambda^2\sin^2\varphi > 1$ ; les racines $a_{1,2} = -i\left( \lambda\sin\varphi \mp \sqrt{\lambda^2\sin^2\varphi - 1} \right)$ sont distinctes, imaginaires pures, de modules différents et de produit de module $1$ : l'une est de module $> 1$.

La condition de Von Neumann est donc satisfaite pour tout $\xi$ ssi $|\lambda| < 1$, *i.e.* $|V|\Delta t < \Delta x$.
