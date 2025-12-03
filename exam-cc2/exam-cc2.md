![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4

# EDP1 
# 2025-26
# Exam CC no. 2

**Durée 2H00. Documents non autorisés. Tous les exercices sont indépendants. Le barème prévisionnel est indiqué pour chaque exercice. Rendre sur deux copies séparées les exercices 1 et 2 d'une part, 3 et 4 d'autre part.**

## Exercice 1 (5 points)

On considère, pour l'équation d'advection avec $V > 0$

$$ u_t + V u_x = 0, $$

le schéma *décentré amont* :

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + V\frac{u_j^n-u_{j-1}^n}{\Delta x} = 0. $$

### 1.1

Montrer que, sous une condition CFL que l'on précisera, ce schéma est stable au sens $L^\infty$.

### 1.2

Montrer que ce schéma est consistant d'ordre $1$ en temps et en espace.

### 1.3

Montrer que l'équation équivalente de ce schéma est

$$ u_t + V u_x - \nu u_{xx} = 0 $$

avec $\nu$ une constante (dépendant de $\Delta t$ et $\Delta x$) que l'on précisera. Sous la condition CFL précédente, quel comportement numérique de la méthode observe-t-on quand $n$ croît ?

## Exercice 2 (5 points)

Étudier la stabilité au sens de Von Neumann du schéma de Gear ($\nu > 0$) :

$$ \frac{3u_j^{n+1}-4u_j^n+u_j^{n-1}}{2\Delta t}-\nu \frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2}=0. $$

## Exercice 3 (5 points)

### 3.1

On utilise le schéma de Lax-Wendroff ci-dessous ($\sigma := V\Delta t/\Delta x$) pour résoudre l'équation d'advection avec $V > 0$ et CL périodiques : 

$$ u_j^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n)
             + \frac{\sigma^2}{2}(u_{j+1}^{n+1} - 2u_j^{n+1} + u_{j-1}^{n+1}). $$ 

On note $U^n := (u_j^n)_j$, montrer que ce schéma s'écrit $AU^{n+1} = W^n$ avec $A$ une matrice et $W$ un vecteur que l'on précisera. Montrer que cette matrice est inversible.

### 3.2

Dans la portion de code associé ci-dessous, quel est le rôle de l'expression `u = F \ w` et que représente `F` ?

```julia=
for n ∈ 1:N
    w = u - ...
    u = F \ w
end
```

### 3.3

Dans le code ci-dessous, à quoi sert la commande `spzeros` ? Compléter la ligne 11.

```julia=
function laplacian(G; Δx = 1.0, Δy = 1.0)
    ij = findall(G .≠ 0)
    K = length(ij)
    p = sortperm(G[ij])
    ij = ij[p] # ij is reordered so that G[ij[k]] = k
    Δ = spzeros(K, K)
    for k ∈ 1:K
        Δ[k, k] = -2 / Δx^2 - 2 / Δy^2
        ijₖ = ij[k]
        kW = G[ijₖ + CartesianIndex(0, -1)]
        if kW ≠ 0 Δ[k, kW] = ... end # à compléter
        kS = G[ijₖ + CartesianIndex(-1, 0)]
        if kS ≠ 0 Δ[k, kS] = 1 / Δy^2 end
        kE = G[ijₖ + CartesianIndex(0,  1)]
        if kE ≠ 0 Δ[k, kE] = 1 / Δx^2 end
        kN = G[ijₖ + CartesianIndex( 1, 0)]
        if kN ≠ 0 Δ[k, kN] = 1 / Δy^2 end
    end
    return Δ
end
```

## Exercice 4 (5 points)

On cherche $u$ dans $\mathscr{C}^2(\overline{\Omega})$ (où $\Omega$ est un ouvert borné de $\mathbf{R}^n$ à bord $\Gamma := \partial\Omega$ régulier) telle que

$$ -\Delta u(x) + u(x) = f(x),\quad x \in \Omega, $$

avec conditions aux limites de Neumann

$$ \frac{\partial u}{\partial n} = g,\quad x \in \Gamma, $$

et $f, g$ des fonctions continues sur $\overline\Omega$ et $\Gamma$, respectivement.

### 4.1

Montrer que toute solution forte $u$ est solution (dite *faible*) d'un nouveau problème de la forme

$$ (\forall v \in H^1(\Omega)) : a(u, v) = \varphi \cdot v, $$

avec $a$ bilinéaire et $\varphi$ linéaire que l'on précisera.

### 4.2

Montrer qu'on a existence et unicité de solution faible.

### 4.3

Montrer que toute solution faible dans $\mathscr{C}^2(\overline{\Omega})$ est solution forte.
