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

**Réponse.** Étant donnés les ordres recherchés en $t$ en $x$, on fait un DL en $(t_{n+1},x_j)$ à l'ordre $3$ en temps, $4$ en espace. En écrivant, notamment, que

$$ u(t_{n-1},x_j) = u + \frac{\partial u}{\partial t}(-2\Delta t) + \frac{1}{2}\frac{\partial^2 u}{\partial t^2}(-2\Delta t)^2 + \frac{1}{6}\frac{\partial^3 u}{\partial t^3}(-2\Delta t)^3 + O(\Delta t^4) $$
   
où le membre de droite est évalué en $(t_{n+1},x_j)$, on obtient que l'erreur de troncature vérifie

$$ E^{n+1}_j = -\frac{1}{3}\frac{\partial^3 u}{\partial t^3}\Delta t^2 + O(\Delta t^3)
   -\frac{\nu}{12}\frac{\partial^4 u}{\partial x^4}\Delta x^2 + O(\Delta x^3). $$

## Exercice 2 (4 points)

### 2.1

On considère le schéma de Dufort-Frankel (avec $\nu > 0$) :

$$ \frac{u_j^{n+1} - u_j^{n-1}}{2\Delta t} - \nu \frac{u_{j+1}^n - u_j^{n+1} - u_j^{n-1} + u_{j-1}^n}{\Delta x^2} = 0. $$

Montrer que, sous une condition CFL que l'on précisera, $u_j^{n+1}$ s'écrit comme combinaison convexe des $(u_j^n)_j$ et des $(u_j^{n-1})_j$.


**Réponse.** Avec $\sigma := \nu\Delta t/\Delta x$, on a

$$ u^{n+1}_j = \frac{1}{1+2\sigma}( (1-2\sigma)u^{n-1}_j + 2\sigma u^n_{j-1} + 2\sigma u^n_{j+1} ) $$

qui est bien une combinaison convexe si $(0 <)\ \sigma \leq 1/2$ (CFL).

### 2.2

En déduire la stabilité $L^\infty$ de ce schéma sous cette même condition.

**Réponse.** Sous cette condition, si les $(u^0_j)_j$ ainsi que les $(u^1_j)_j$ sont dans le convexe $[m,M]$, on vérifie par une récurrence immédiate que tous les $(u^n)_j$ aussi pour $n \geq 2$, d'où la stabilité $L^\infty$.

## Exercice 3 (4 points)

### 3.1

On considère l'équation des ondes sur un domaine $\Omega \subset \mathbf{R}$

$$ \frac{\partial^2 u}{\partial t^2}(t,x) - \frac{\partial^2 u}{\partial x^2}(t,x) = 0,\quad t > 0,\quad x \in \Omega. $$

Montrer qu'en posant $v := \partial u/\partial t$, $w := \partial u/\partial x$ et $U : = (v,w)$ on peut mettre l'équation sous la forme

$$ \frac{\partial U}{\partial t}(t,x) - J \frac{\partial U}{\partial x}(t,x) = 0,\quad t > 0,\quad x \in \Omega, $$

où $J$ est une matrice $2 \times 2$ que l'on précisera.

**Réponse.** En supposant $u$ assez régulière et en utilisant le théorème de Schwarz (égalité des dérivées partielles croisées secondes), on a la forme voulue avec


$$ J = \left[ \begin{array}{cc} 0 & 1\\ 1 & 0 \end{array} \right]. $$

### 3.2

On considère le schéma vectoriel associé ci-dessous :

$$ \frac{U_j^{n+1} - U_j^n}{\Delta t} - J\frac{U_{j+1}^{n+1} - U_{j-1}^{n+1}}{2\Delta x} = 0. $$

Montrer que la stabilité $L^2$ de ce schéma revient à étudier, pour chaque mode $k \in \mathbf{Z}$, les valeurs propres d'une matrice $A(k)$ que l'on précisera. (On ne demande ni de déterminer ces valeurs propres, ni de faire l'étude de stabilité.)

**Réponse.** En injectant un mode de Fourier (vectoriel) $k$ sous la forme $U^n_j = e^{2i\pi kj\Delta x}A(k)$, on obtient

$$ (I - i\sigma\sin\alpha_k\,J)A(k) = I $$

où $\sigma := \Delta t/\Delta x$, et où $\alpha_k := 2\pi k\Delta x$. Le déterminant du facteur de gauche valant $1 + \sigma^2\sin^2\alpha_k \geq 1 > 0$, on en déduit que

$$ A(k) = (I - i\sigma\sin\alpha_k\,J)^{-1}. $$

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

**Réponse.** Il s'agit du schéma de Lax-Wendroff implicite : 

$$ u_j^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n)
             + \frac{\sigma^2}{2}(u_{j+1}^{n+1} - 2u_j^{n+1} + u_{j-1}^{n+1}). $$  

### 4.2

Quel est l'intérêt d'utiliser la fonction `spdiagm` pour construire la matrice `A` ?

**Réponse.** La matrice construite est creuse (*sparse*) ce qui occasione un stockage réduit et permet une factorisation plus efficace.

### 4.3

Expliquer pourquoi on calcule `w[1]` *après* `w[2:end-1]`.

**Réponse.** La valeur mise à jour de `w[1]` serait sinon utilisée (à tort) pour la mise à jour de `w[2:end]`.

### 4.4

Expliquer le rôle de la ligne `u = F \ w`.

**Réponse.** On réutilise la factorisation (creuse) de Cholesky, faite une fois pour toute avant la boucle d'itération sur le pas de temps, pour faire une résolution du système linéaire `A * u = w` par simple remontée triangulaire.

## Exercice 5 (3 points)

On cherche $u$ dans $\mathscr{C}^2(\overline{\Omega})$ (où $\Omega$ est un ouvert de $\mathbf{R}^n$ à bord régulier)
telle que

$$ -\Delta u(x) + u(x) = f(x),\quad x \in \Omega, $$

avec conditions aux limites de Dirichlet, $u = 0$ sur $\partial\Omega$, et $f$ une fonction continue sur
$\overline\Omega$. On admet que le sous-espace

$$ H_0^1(\Omega) := \lbrace v \in H^1(\Omega)\ |\ v = 0 \text{ sur } \partial\Omega \rbrace $$

muni du produit scalaire de $H^1(\Omega)$ est encore un espace de Hilbert. Montrer que toute solution du problème précédent, dite *solution forte*, est encore solution d'un nouveau problème que l'on précisera.

**Réponse.** Si $u$ est solution forte, pour tout $v$ dans $H^1_0(\Omega)$ on a

$$ \int_\Omega (-\Delta u + u)v\,\mathrm{d}x = \int_\Omega fv\,\mathrm{d}x. $$

En utilisant le théorème de Stockes, on en déduit ($\Gamma := \partial\Omega$)

$$ -\int_\Gamma \frac{\partial u}{\partial n}v\,\mathrm{d}\sigma + (u|v)_{H^1} = \int_\Omega fv\,\mathrm{d}x. $$

L'intégrale sur le bord étant nulle grâce à l'appartenance de $v$ à $H^1_0$, on déduit la formulation faible attendue.
