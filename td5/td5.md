![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1
# 2025-26
# TD 4 - Ondes


## Exercice 1

Étudier la stabilité au sens de Von Neumann du $\theta$-schéma ci-dessous, $\theta \in [0,1/2]$ :

$$ \frac{u_j^{n+1}-2u_j^n+u^{n-1}_j}{\Delta t^2} - \theta \frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2} - (1-2\theta) \frac{u_{j+1}^{n}-2u_j^{n}+u_{j-1}^{n}}{\Delta x^2} - \theta \frac{u_{j+1}^{n-1}-2u_j^{n-1}+u_{j-1}^{n-1}}{\Delta x^2}= 0. $$

## Exercice 2

Étudier la stabilité au sens de Von Neumann du schéma vectoriel de Lax-Friedrichs ci-dessous, où $A$ est une matrice que l'on précisera :

$$ \frac{2U_j^{n+1}-U_{j+1}^{n}-U_{j-1}^n}{2\Delta t} + A \frac{U_{j+1}^{n}-U_{j-1}^{n}}{2\Delta x} = 0. $$

## Exercice 3

Étudier la stabilité au sens de Von Neumann du schéma vectoriel de Lax-Wendroff ci-dessous :

$$ \frac{2U_j^{n+1}-U_{j+1}^{n}-U_{j-1}^n}{2\Delta t} + A \frac{U_{j+1}^{n}-U_{j-1}^{n}}{2\Delta x} + A^2\Delta t \frac{U_{j+1}^{n}-2U_{j}^{n}+U_{j-1}^{n}}{2\Delta x^2}= 0. $$
