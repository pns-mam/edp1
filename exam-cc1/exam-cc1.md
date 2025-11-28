![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1 
# 2025-26
# Exam CC no. 1

**Durée 2H00. Documents non autorisés. Tous les exercices sont indépendants.
Le barème prévisionnel est indiqué pour chaque exercice. Rendre sur deux copies séparées les exercices 1 et 2 d'une part, 3 d'autre part.**

## Exercice 1 (6 points)

On considère l'équation des ondes sur le domaine $\Omega = ]0,1[$ avec conditions aux limites périodiques :

$$ \begin{aligned}
  & \frac{\partial^2 u}{\partial t^2}(t, x) - \frac{\partial^2 u}{\partial x^2}(t, x) = 0,\quad t>0,\quad x \in \Omega,\\
  & u(t, 0) = u(t, 1),\quad 
    \frac{\partial u}{\partial x}(t, 0) = \frac{\partial u}{\partial x}(t, 1),\quad t > 0,\\
  & u(0, x) = u_0(x),\quad \frac{\partial u}{\partial t}(0, x) = u_1(x), \quad x \in \Omega.
\end{aligned} $$

### 1.1

On définit l'énergie

$$ E(t) := \frac{1}{2} \int_0^1 \left( \frac{\partial u}{\partial t}(t,x) \right)^2 +
   \left( \frac{\partial u}{\partial x}(t,x) \right)^2\,\mathrm{d}x. $$

En supposant $u$ solution aussi régulière que nécessaire, montrer que l'énergie est constante.

**Réponse.** 

$$ 
\frac{dE}{dt}=\int_0^1\frac{\partial^2 u}{\partial t^2}\frac{\partial u}{\partial t} dx +\int_0^1\frac{\partial^2 u}{\partial t \partial x}\frac{\partial u}{\partial x} dx.
$$

En intégrant par partie: 

$$ 
\int_0^1\frac{\partial^2 u}{\partial t \partial x}\frac{\partial u}{\partial x} dx= \left[\frac{\partial u}{\partial t}\frac{\partial u}{\partial x}\right]_0^1-
\int_0^1\frac{\partial u}{\partial t}\frac{\partial^2 u}{\partial x^2}.
$$

Or $\partial_x u(t,0)=\partial_x u(t,1)$ et $u(t,0)=u(t,1)$ pour tout $t\geq 0$. En particulier on a alors en dérivatn par rapport à $t$: $\partial_t u(t,0)=\partial_t u(t,1)$.
Comme de plus $\partial^2_t u(t, x) = \partial^2_x u(t, x)$, on en déduit 

$$ 
\frac{dE}{dt}=\int_0^1\frac{\partial^2 u}{\partial t^2}\frac{\partial u}{\partial t} dx -\int_0^1\frac{\partial^2 u}{ \partial t}\frac{\partial^2 u}{\partial t^2} dx=0.
$$

### 1.2

En déduire qu'on a unicité de solution pour l'EDP considérée. 

**Réponse.** Si $u$ et $v$ sont deux solutions, alors $w=u-v$ est slution de l'équation des ondes (sur domaine périodique) avec conditions initiales nulles. En particulier, on en déduit que l'énergie associée à $w$ est constamment nulle ($u(0,x)=0$ implique $\partial_x u(0,x)=0$. De plus, on en déduit (E étant l'intégrale d'une somme de carrés) que $\partial_t w(t,x)=0$ et $\partial_x w(t,x)=0$ pour tout $x\in[0,1]$, $t\geq 0$. De cela on conclut que $w$ est la constante nulle.

## Exercice 2 (6 points)

On considère le schéma explicite ci-dessous

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + \frac{u_{j+1}^n-u_{j-1}^n}{2\Delta x} -
   \Delta t \frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{2\Delta x^2} = 0 $$

pour l'équation d'advection

$$ \frac{\partial u}{\partial t}(t, x) + \frac{\partial u}{\partial x}(t, x) = 0. $$

### 2.1

En supposant $u$ solution aussi régulière que nécessaire, montrer que

$$ \frac{\partial^2 u}{\partial t^2}(t, x) = \frac{\partial^2 u}{\partial x^2}(t, x). $$

**Réponse.** 

$$ \frac{\partial^2 u}{\partial t^2}=-\frac{\partial }{\partial t}\frac{\partial u}{\partial x}=-\frac{\partial }{\partial x}\frac{\partial u}{\partial t}=\frac{\partial^2 u}{\partial x^2}. $$

### 2.2

En déduire que le schéma proposé est consistant d'ordre $2$ en temps et en espace. (On pourra utiliser des développements à l'ordre $4$ en temps et en espace.)

**Réponse.** Taylor à l'ordre 4 :

$$ u_j^{n+1} = u + \Delta t \ u_t + \frac{\Delta t^2}{2} u_{tt} + \frac{\Delta t^3}{6} u_{ttt} + \frac{\Delta t^4}{24} u_{tttt} + O(\Delta t^5),$$

$$u_{j\pm 1}^n = u \pm \Delta x\ u_x + \frac{\Delta x^2}{2} u_{xx} \pm \frac{\Delta x^3}{6} u_{xxx} + \frac{\Delta x^4}{24} u_{xxxx} + O(\Delta x^5).$$

En remplaçant par les développements de Taylor et en utilisant $u_{tt} = u_{xx}$ :

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} = u_t + \frac{\Delta t}{2} u_{xx} + \frac{\Delta t^2}{6} u_{ttt} + O(\Delta t^3),$$

$$ \frac{u_{j+1}^n-u_{j-1}^n}{2\Delta x} = u_x + \frac{\Delta x^2}{6} u_{xxx} + O(\Delta x^4), $$

$$ \Delta t \frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{2\Delta x^2} = \frac{\Delta t}{2} u_{xx} + \frac{\Delta t \Delta x^2}{24} u_{xxxx} + O(\Delta t \Delta x^4). $$

Le résidu du schéma, i.e. l’erreur de consistance $E_j^n$, est donc

$$ E_j^n = \frac{u_j^{n+1}-u_j^n}{\Delta t} + \frac{u_{j+1}^n-u_{j-1}^n}{2\Delta x} - \Delta t \frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{2\Delta x^2}
= \frac{1}{6}(\Delta x^2 - \Delta t^2) u_{xxx} - \frac{\Delta t \Delta x^2}{24} u_{xxxx} + O(\Delta t^3) + O(\Delta x^3) + O(\Delta t \Delta x^4).
$$

Le schéma est donc consistant d'ordre 2 en temps et en espace.




## Exercice 3 (8 points)

### 3.1

Montrer que le schéma explicite ci-dessous

$$ \frac{u_j^{n+1}-2u_j^n+u_j^{n-1}}{\Delta t^2} - \frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{\Delta x^2} = 0 $$

pour l'équation des ondes (voir Exercice 1) est consistant d'ordre $2$ en temps et en espace. (On pourra utiliser des développements à l'ordre $4$ en temps et en espace.)

**Réponse.** 

### 3.2

Montrer que ce schéma est stable au sens de Von Neumann sous une condition que l'on précisera.

**Réponse.** 

### 3.3

Comment pourrait-on améliorer la stabilité ?

**Réponse.** En implicitant le schéma, par exemple selon

$$ \frac{u_j^{n+1}-2u_j^n+u_j^{n-1}}{\Delta t^2} - \frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2} = 0. $$
