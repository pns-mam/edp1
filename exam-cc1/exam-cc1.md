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

### 1.2

En déduire qu'on a unicité de solution pour l'EDP considérée. 

## Exercice 2 (6 points)

On considère le schéma explicite ci-dessous

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + \frac{u_{j+1}^n-u_{j-1}^n}{2\Delta x} -
   \Delta t \frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{2\Delta x^2} = 0 $$

pour l'équation d'advection

$$ \frac{\partial u}{\partial t}(t, x) + \frac{\partial u}{\partial x}(t, x) = 0. $$

### 2.1

En supposant $u$ solution aussi régulière que nécessaire, montrer que

$$ \frac{\partial^2 u}{\partial t^2}(t, x) = \frac{\partial^2 u}{\partial x^2}(t, x). $$

### 2.2

En déduire que le schéma proposé est consistant d'ordre $2$ en temps et en espace. (On pourra utiliser des développements à l'ordre $4$ en temps et en espace.)

## Exercice 3 (8 points)

### 3.1

Montrer que le schéma explicite ci-dessous

$$ \frac{u_j^{n+1}-2u_j^n+u_j^{n-1}}{\Delta t^2} - \frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{\Delta x^2} = 0 $$

pour l'équation des ondes (voir Exercice 1) est consistant d'ordre $2$ en temps et en espace. (On pourra utiliser des développements à l'ordre $4$ en temps et en espace.)

### 3.2

Montrer que ce schéma est stable au sens de Von Neumann sous une condition que l'on précisera.

### 3.3

Comment pourrait-on améliorer la stabilité ?
