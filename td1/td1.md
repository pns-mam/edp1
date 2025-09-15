![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1
# 2025-26
# TD 1 - Propriétés qualitatives des solutions d'EDP

Le but de cette série d'exercices est de montrer quelques propriétés qualitatives des solutions régulières de certaines équations aux dérivées partielles de nature différente (chaleur, équation des ondes, Schrödinger) qui font notamment intervenir la notion d'*énergie*.

## Exercice 1

On considère le problème de la chaleur en une dimension d'espace posé dans le domaine Ω=(0,1) :

$$\left\{\begin{array}{rcl}
\frac{\partial u}{\partial t} -  \frac{\partial^2 u}{\partial x^2}&=
&0,\quad t > 0,\quad x \in (0,1),\\
u(t,x)&=&0,\quad t > 0,\quad x \in \{0,1\},\\
u(0,x)&=&u_0(x),\, x\in (0,1).
\end{array}\right.$$

### 1.1

Montrer d'abord que toute fonction $v(x)$ continûment dérivable sur $[0,1]$ t.q. $v(0)=0$, vérifie l'inégalité de Poincaré

$$\int_0^1 v^2(x) \mathrm{d}x \le \int_0^1\left(\frac{\mathrm{d}v}{\mathrm{d}x}\right)^2 \mathrm{d}x.$$

(On écrira d'abord $v(x)$ comme intégrale de sa dérivée sur l'intervalle $(0,x)$ et ensuite on appliquera l'inégalité de Cauchy-Schwarz).

### 1.2

On notera dans ce qui suit par 
$$E(t) = \int_0^1 u^2(t,x) \mathrm{d}x$$
l'énergie à l'instant $t$. En multipliant l'équation de la chaleur par $u$ et en intégrant par rapport à $x$, établir l'égalité vérifiée par cette quantité:
$$\frac{1}{2}\frac{\mathrm{d}E(t)}{\mathrm{d}t} =-\int_0^1\left(\frac{\partial u}{\partial x}\right)^2 \mathrm{d}x.$$

### 1.3

En appliquant l'inégalité de Poincaré à $v(x)= u(t,x)$ en déduire que
$$E(t) \le E(0)e^{-2t}.$$
On dit que l'énergie *décroît exponentiellement en temps*.

## Exercice 2

On se propose de caractériser la solution régulière du problème des ondes en une dimension d'espace dans le domaine $\Omega\subset \mathbf{R}$ : 

$$\left\{\begin{array}{rcl}
\frac{\partial^2 u}{\partial t^2} -  \frac{\partial^2 u}{\partial x^2} &= &0,\quad t > 0,\quad
x \in\Omega,\\
u(t,x)&=&0,\quad t > 0,\quad x \in\partial\Omega,\\
u(0,x)&=&u_0(x),\quad x\in\Omega,\\
\frac{\partial u}{\partial t}(0,x)& = & u_1(x),\quad x\in \Omega,
\end{array}\right.$$

où $u_0$ et $u_1$ sont des fonctions régulières et $U_1$ une primitive de $u_1$.

### 2.1

Supposons maintenant que $\Omega= (0,1)$. En multipliant l'équation par $\frac{\partial u}{\partial t}$ et en intégrant par rapport à $x$, établir *l'égalité de l'énergie*:

$$\frac{\mathrm{d}}{\mathrm{d}t} \int_0^1 \left(\frac{\partial u}{\partial t}(t,x)\right)^2 + \left(\frac{\partial u}{\partial x}(t,x)\right)^2 \mathrm{d}x = 0.$$

### 2.2

En utilisant l'égalité de l'énergie montrer que cette solution est unique. (On choisira $u_0=u_1=0$ et on montrera que l'unique solution est celle nulle.)

### 2.3

Vérifier que la fonction suivante est solution du problème des ondes si $\Omega=\mathbf{R}$

$$u(t,x) = \frac{1}{2}(u_0(x+t)+u_0(x-t))+\frac{1}{2}(U_1(x+t)-U_1(x-t)).$$

Cette dernière relation porte le nom de *formule de d'Alembert*.

## Exercice 3

L'équation de Schrödinger décrit l'évolution de la fonction d'onde $u:\mathbf{R}^n \times \mathbf{R} \rightarrow \mathbf{C}$ d'une particule soumise à un potentiel $V:\mathbf{R}^n \rightarrow \mathbf{R}$. La quantité $|u|^2$ s'interprète comme la densité de probabilité pour détecter que la particule se trouve au point $(t,x)$.

On se propose de montrer les principes de conservation de l'énergie pour une solution régulière de l'équation de Schrödinger uni-dimensionnelle:

$$\left\{\begin{array}{rcl}
i\frac{\partial u}{\partial t}+\frac{\partial^2 u}{\partial x^2} -V u &=&0,\quad t > 0, \quad x \in\mathbf{R},\\
\lim_{|x|\rightarrow\infty}u(t,x)  &=& 0,\, t > 0,\\
\lim_{|x|\rightarrow\infty}\frac{\partial u}{\partial x}(t,x) &=& 0,\quad t > 0,\\
u(0,x) & = & u_0(x),\,x\in \mathbf{R}.
\end{array}\right.$$

### 3.1

Montrer que pour toute fonction dérivable $v : \mathbf{R}^2 \to \mathbf{R}$ on a

$$\mathrm{Re}\left(\frac{\partial v}{\partial t}\bar v\right)(t,x) = \frac{1}{2}\frac{\partial |v|^2}{\partial t}(t,x),$$

où $\mathrm{Re}(v)$ désigne la partie réelle de la fonction $v$, et $\bar v$ sa conjuguée.

### 3.2 

En multipliant l'équation par $\bar u$ et en intégrant par rapport à $x$, établir l'égalité de l'énergie:

$$\int_{\mathbf{R}}|u(t,x)|^2 \mathrm{d}x = \int_{\mathbf{R}}|u_0(x)|^2 \mathrm{d}x.$$

### 3.3

En multipliant l'équation par $\frac{\partial \bar u}{\partial t}$, montrer que

$$\int_{\mathbf{R}}\left(\left|\frac{\partial u}{\partial x}(x)\right|^2 +V(x) |u(x,t)|^2\right) \mathrm{d}x = \int_{\mathbf{R}}\left(\left|\frac{\partial u_0}{\partial x}(x)\right|^2 +V(x) |u_0(x)|^2\right) \mathrm{d}x.$$
