---
header-includes:
  - \usepackage{mathrsfs}
---
![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1 
# 2024-25
# Exam CC no. 2

**Durée 2H00. Documents non autorisés. Tous les exercices sont indépendants.
Le barème prévisionnel est indiqué pour chaque exercice.**

## Exercice 1 (5 points)

On considère le schéma de Gear

$$ \frac{3u_j^{n+1} - 4u_j^n + u_j^{n-1}}{2\Delta t} - \nu\frac{u_{j+1}^{n+1} - 2u_j^{n+1} + u_{j-1}^{n+1}}{\Delta x^2} = 0 $$

pour l'équation de la chaleur 

$$ \frac{\partial u}{\partial t} - \nu\frac{\partial^2 u}{\partial x^2} = 0. $$

Montrer que le schéma est consistant à l'ordre $2$ à la fois en temps et en espace.

## Exercice 2 (4 points)

### 2.1

On considère le schéma de Dufort-Frankel (avec $\nu > 0$) :

$$ \frac{u_j^{n+1} - u_j^{n-1}}{2\Delta t} - \nu \frac{u_{j+1}^n - u_j^{n+1} - u_j^{n-1} + u_{j-1}^n}{\Delta x^2} = 0. $$

Montrer que, sous une condition CFL que l'on précisera, $u_j^{n+1}$ s'écrit comme combinaison convexe des $(u_j^n)_j$.

### 2.2

En déduire la stabilité $L^\infty$ de ce schéma sous cette même condition.

## Exercice 3 (4 points)

### 3.1

On considère l'équation des ondes sur un domaine $\Omega \subset \mathbf{R}$

$$ \frac{\partial^2 u}{\partial t^2}(t,x) - \frac{\partial^2 u}{\partial x^2}(t,x) = 0,\quad t > 0,\quad x \in \Omega. $$

Montrer qu'en posant $v := \partial u/\partial t$, $w := \partial u/\partial x$ et $U : = (v,w)$ on peut mettre l'équation sous la forme

$$ \frac{\partial U}{\partial t}(t,x) - J \frac{\partial U}{\partial x}(t,x) = 0,\quad t > 0,\quad x \in \Omega, $$

où $J$ est une matrice $2 \times 2$ que l'on précisera.

### 3.2

On considère le schéma vectoriel associé ci-dessous :

$$ \frac{U_j^{n+1} - U_j^n}{\Delta t} - J\frac{U_{j+1}^{n+1} - U_{j-1}^{n+1}}{2\Delta x} = 0. $$

Montrer que la stabilité $L^2$ de ce schéma revient à étudier, pour chaque mode $k \in \mathbf{Z}$, les valeurs propres d'une matrice $A(k)$ que l'on précisera. (On ne demande ni de déterminer ces valeurs propres, ni de faire l'étude de stabilité.)

## Exercice 4 (4 points)

On considère la portion de code ci-dessous pour l'équation d'advection : 
```julia
A = spdiagm(-1 => -.5σ^2 * ones(J),
             0 => (1 + σ^2) * ones(J + 1),
             1 => -.5σ^2 * ones(J))
A[1, end] = -.5σ^2
A[end, 1] = -.5σ^2

F = cholesky(A)

w = zeros(J + 1)

for n ∈ 1:N
    w[2:end-1] = u[2:end-1] - .5σ * (u[3:end] - u[1:end-2])
    w[1] = u0(x[1] - V * n * Δt)
    w[end] = w[1]
    u = F \ w
end
```

### 4.1

Quel est le schéma numérique utilisé ? (On donnera son expression.)

### 4.2

Quel est l'intérêt d'utiliser la fonction `spdiagm` pour construire la matrice `A` ?

### 4.3

Expliquer pourquoi on calcule `w[1]` *après* `w[2:end-1]`.

### 4.4

Expliquer le rôle de la ligne `u = F \ w`.

## Exercice 5 (3 points)

On cherche $u$ dans $\mathscr{C}^2(\overline{\Omega})$ (où $\Omega$ est un ouvert de $\mathbf{R}^n$ à bord régulier)
telle que

$$ -\Delta u(x) + u(x) = f(x),\quad x \in \Omega, $$

avec conditions aux limites de Dirichlet, $u = 0$ sur $\partial\Omega$, et $f$ une fonction continue sur
$\overline\Omega$. On admet que le sous-espace

$$ H_0^1(\Omega) := \lbrace v \in H^1(\Omega)\ |\ v = 0 \text{ sur } \partial\Omega \rbrace $$

muni du produit scalaire de $H^1(\Omega)$ est encore un espace de Hilbert. Montrer que toute solution du problème précédent, dite *solution forte*, est encore solution d'un nouveau problème que l'on précisera.
