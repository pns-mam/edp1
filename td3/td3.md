![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1
# 2025-26
# TD 3 - Stabilité (1/2)


## Exercice 1

Étudier la stabilité au sens de Von Neumann du $\theta$-schéma, $\theta \in [0,1]$, avec coefficient de diffusion $\nu > 0$ :

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t}-(1-\theta)\nu \frac{u_{j+1}^{n}-2u_j^n+u_{j-1}^{n}}{\Delta x^2}-\theta\nu \frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2}=0. $$

## Exercice 2

Étudier la stabilité au sens de Von Neumann du schéma de Crank-Nicolson ci-dessous, avec vitesse d'advection $V$ :

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + V \frac{u_{j+1}^{n}-u_{j-1}^{n}}{4\Delta x} + V \frac{u_{j+1}^{n+1}-u_{j-1}^{n+1}}{4\Delta x} = 0. $$

## Exercice 3

### 3.1

Étudier la stabilité au sens de Von Neumann du schéma Leapfrog (avec $\nu > 0$) :

$$ \frac{u_j^{n+1}-u_j^{n-1}}{2\Delta t}-\nu \frac{u_{j+1}^{n}-2u_j^{n}-+u_{j-1}^{n}}{\Delta x^2}=0. $$

### 3.2

Étudier la stabilité au sens de Von Neumann du schéma de Du Fort-Frankel : 

$$ \frac{u_j^{n+1}-u_j^{n-1}}{2\Delta t}-\nu \frac{u_{j+1}^{n}-u_j^{n+1}-u_j^{n-1}+u_{j-1}^{n}}{\Delta x^2}=0. $$
