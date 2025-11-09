![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1 
# 2025-26
# TD 6 - Formulations variationnelles

## Exercice 1

On cherche $u$ dans $\mathscr{C}^2(\overline{\Omega})$ (où $\Omega$ est un ouvert borné de $\mathbf{R}^n$ à bord régulier - de classe $\mathscr{C}^1$) telle que

$$ -\Delta u(x) + \alpha u(x) = f(x),\quad x \in \Omega, $$

avec $\alpha > 0$ et conditions aux limites de Dirichlet, $u = 0$ sur $\partial\Omega$, et $f$ une fonction continue sur $\overline\Omega$.

### 1.1

On admet que l'application trace $v \mapsto v_{|\partial \Omega}$ de $\mathscr{C}^1(\overline{\Omega})$ dans $L^2(\Omega)$ se prolonge continûment à $H^1(\Omega)$. En déduire que le sous-espace

$$ H_0^1(\Omega) := \lbrace v \in H^1(\Omega)\ |\ v_{|\partial\Omega} = 0 \rbrace $$

muni du produit scalaire de $H^1(\Omega)$ est encore un espace de Hilbert. 

### 1.2


Montrer que toute solution du problème précédent, dite *solution forte*, est encore solution (dite *solution faible*) d'un nouveau problème que l'on précisera.

### 1.3

Montrer qu'on a existence et unicité de solution faible. Interpréter de façon variationnelle.

### 1.4

Montrer que toute solution faible dans $\mathscr{C}^2(\overline{\Omega})$ est solution forte.

## Exercice 2

On cherche $u$ dans $\mathscr{C}^2(\overline{\Omega})$ (où $\Omega$ est un ouvert borné de $\mathbf{R}^n$ à bord régulier - de classe $\mathscr{C}^1$)
telle que (problème de Poisson)

$$ -\Delta u(x)= f(x),\quad x \in \Omega, $$

avec conditions aux limites de Robin

$$ u + \frac{\partial u}{\partial n} = g,\quad x \in \partial\Omega, $$

et $f, g$ des fonctions continues sur $\overline\Omega$ et $\partial\Omega$, respectivement.

### 2.1

Montrer que toute solution solution forte est encore solution d'un nouveau problème que l'on précisera.

### 2.2

À l'aide de l'inégalité de Poincaré ci-dessous,

$$ \|v\|_2 \leq C(\|v_{|\partial \Omega}\|_2 + \|\nabla v\|_2),\quad v \in H^1(\Omega), $$

montrer qu'on a existence et unicité de solution faible. Interpréter de façon variationnelle.

### 2.3

Montrer que toute solution faible dans $\mathscr{C}^2(\overline{\Omega})$ est solution forte.
