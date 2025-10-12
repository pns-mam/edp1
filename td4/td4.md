![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1
# 2025-26
# TD 4 - Stabilité (2/2)


## Exercice 1

### 1.1

Montrer que sous une condition CFL que l'on précisera, le schéma d'Euler explicite pour la chaleur (avec coefficient de diffusion $\nu > 0$) ci-dessous est tel que $u^{n+1}_j$ s'écrit comme combinaison convexe de $u^n_{j+1}$, $u^n_j$, $u^n_{j-1}$ : 

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} - \nu \frac{u_{j+1}^{n}-2u_j^n+u_{j-1}^{n}}{\Delta x^2} = 0,\quad j \in \mathbf{Z},\quad n \in \mathbf{N}. $$

### 1.2

En déduire que, sous cette condition, ce schéma est stable au sens $L^\infty$, c'est-à-dire que, si $u^0_j \in [m, M]$ pour tout $j \in \mathbf{Z}$, on a aussi

$$ u^n_j \in [m, M],\quad j \in \mathbf{Z},\quad n \in \mathbf{N}. $$

### 1.3

Supposons que la condition CFL précédente n'est pas vérifiée ; montrer que, pour $u^0_j = (-1)^j$, $j \in \mathbf{Z}$, les $u^n_j$ ne sont pas bornés. 

## Exercice 2

Montrer que le schéma d'Euler implicite pour la chaleur (avec coefficient de diffusion $\nu > 0$)

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} - \nu \frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2} = 0,\quad j=1,\dots,J-1,\quad n \in \mathbf{N}, $$

et $u^{n+1}_0 = u^{n+1}_J = 0$ (Dirichlet), se met sous la forme $U^{n+1} = AU^n$, avec $U_n := (u^n_j)_{j=0,\dots,J}$ et $A$ une matrice inversible que l'on précisera.

### 2.2

En déduire que ce schéma est (inconditionnellement) stable au sens $L^\infty$, c'est-à-dire que, si $u^0_j \in [m, M]$ pour tout $j=0,\dots,J$, on a aussi

$$ u^n_j \in [m, M],\quad j=0,\dots,J,\quad n \in \mathbf{N}. $$

## Exercice 3

### 3.1

Montrer que sous une condition CFL que l'on précisera, le schéma décentré amont pour l'équation d'advection avec $V > 0$ ci-dessous est tel que $u^{n+1}_j$ s'écrit comme combinaison convexe de $u^n_{j+1}$, $u^n_j$, $u^n_{j-1}$ : 

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + V \frac{u_{j}^{n}-u_{j-1}^n}{\Delta x} = 0,\quad j \in \mathbf{Z},\quad n \in \mathbf{N}. $$

### 3.2

En déduire que, sous cette condition, ce schéma est stable au sens $L^\infty$, c'est-à-dire que, si $u^0_j \in [m, M]$ pour tout $j \in \mathbf{Z}$, on a aussi

$$ u^n_j \in [m, M],\quad j \in \mathbf{Z},\quad n \in \mathbf{N}. $$

### 3.2

Supposons que la condition CFL précédente n'est pas vérifiée ; montrer que, pour $u^0_j = (-1)^j$, $j \in \mathbf{Z}$, les $u^n_j$ ne sont pas bornés.
