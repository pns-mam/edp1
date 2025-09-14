![PNS](https://raw.githubusercontent.com/pns-mam/edp1/main/logo-pns.png)

## MAM4

# Équations aux dérivées partielles 1
# 2025-26

## 1. Équation de la chaleur

- équation $-u'' = f$ en dimension $1$ et conditions aux limites associées (Dirichlet, Neumann, Robin)
- équation de Poisson-Laplace en dimension deux
- régime instationnaire : équation de la chaleur (conditions initiale et aux limites associées)
- équation de convection (diffusion, advection)
- classification des EDP linéaires du deuxième ordre
- noyau de la chaleur en dimension $1$, principe du maximum, vitesse de propagation infinie

## 2. Consistance

- équation de la chaleur 1D sur un ouvert borné, CL de Dirichlet
- exemples de schémas : explicite, implicite, $\theta$-schéma
- déf. schéma numérique (et support / stencil associé)
- déf. erreur de troncature, consistance et ordre (en temps et en espace) d'un schéma
- équation équivalente associée à un schéma
- exemple : Lax-Friedrichs (advection) et diffusion numérique

## 3. Stabilité
- stabilité au sens $L^\infty$
- stabilité au sens $L^2$ (et interprétation par Fourier)
- exemple : Lax-Friedrichs (advection)
- théorème de Lax
- cas multiniveaux

## 4. Équation des ondes
- équations des ondes 1D sur un domaine borné, CL périodiques
- condition de moyenne nulle (solutions bornées)
- motivation (corde, son, lumière / électromagnétisme)
- autres CL (Dirichlet, Neumann...)
- conservation de l'énergie et unicité
- $\theta$-schéma
- écriture à l'ordre 1 en vectoriel et schéma de Lax-Friedrichs associé

## 5. Introduction aux formulations variationnelles
- équation $-\Delta u + \alpha u = f$ avec CL de Neumann
- déf. espace $H^1(\Omega)$
- A. toute solution forte est solution faible
- Th. de Lax-Milgram
- B. existence et unicité de solution faible
- D. toute solution faible *assez régulière (C.)* est solution forte
