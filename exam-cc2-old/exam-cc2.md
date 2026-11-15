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

**Correction.** En posant $\lambda=\frac{V \Delta t}{\Delta x}$, on a $u_j^{n+1} = (1-\lambda)u_j^{n} +\lambda u_{j-1}^{n}$.

Lorsque $0\leq \lambda \leq 1$, $u_j^{n+1}$ est une combinaison convexe de $(u_j^{n},u_{j-1}^{n})$. En particulier

$$\|u^{n+1}\|_ {\infty}\leq (1-\lambda) \| u^{n} \|_\infty + \lambda  \| u^{n} \| _\infty = \| u^{n} \| _\infty.$$

La condition CFL est donc $\frac{V \Delta t}{\Delta x}\leq 1$.

### 1.2

Montrer que ce schéma est consistant d'ordre $1$ en temps et en espace.

**Correction.** On a

$$ \frac{u_j^{n+1}- u_j^n}{\Delta t}  =  u_t(x_j,t_n) + \frac{\Delta t}{2} u_{tt} + O(\Delta t^2),$$

$$ \frac{u_j^n - u_{j-1}^n}{\Delta x} =   u_x(x_j,t_n) - \frac{\Delta x}{2} u_{xx} + O(\Delta x^2).$$

En substituant dans le schéma, comme $u_t(x_j,t_n) + V u_x(x_j,t_n) = 0$, on a l'erreur de troncature :

$$E_j^n=\frac{u_j^{n+1}- u_j^n}{\Delta t} +V \frac{u_j^n - u_{j-1}^n}{\Delta x} = \frac{1}{2} \left(\Delta t \ u_{tt} - V \Delta x \  u_{xx}\right)+O(\Delta t^2)+ O(\Delta x^2) $$

Le schéma est donc consistant d'ordre 1 en temps et en espace.

### 1.3

Montrer que l'équation équivalente de ce schéma est

$$ u_t + V u_x - \nu u_{xx} = 0 $$

avec $\nu$ une constante (dépendant de $\Delta t$ et $\Delta x$) que l'on précisera. Sous la condition CFL précédente, quel comportement numérique de la méthode observe-t-on quand $n$ croît ?

**Correction.** En remarquant que $u_{t}= -Vu_x$ implique $u_{tt} = V^2 u_{xx}$, on obtient 

$$ \frac{1}{2} \left(\Delta t\  u_{tt} - V \Delta x u_{xx}\right) = \frac{V}{2} \left(V \Delta t - \Delta x \right)u_{xx} =   \frac{V\Delta x}{2} \left(\lambda - 1 \right)u_{xx}$$

Sous la condition CFL $\lambda \leq 1$ et en posant $\nu = \frac{V\Delta x}{2} \left(1-\lambda \right)\geq 0$, on obtient donc l'équation équivalente du schéma 

$$ u_t + V u_x - \nu u_{xx} = 0. $$

L'équation équivalente du schéma est une équation de convection-diffusion, on s'attend donc à observer de la diffusion numérique.

## Exercice 2 (5 points)

Étudier la stabilité au sens de Von Neumann du schéma de Gear ($\nu > 0$) :

$$ \frac{3u_j^{n+1}-4u_j^n+u_j^{n-1}}{2\Delta t}-\nu \frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2}=0. $$

**Correction.** On pose $u_j^n = \hat{u}^n(\xi) e^{2i\pi \xi j\Delta x}$. En substituant dans le schéma et en simplifiant par $e^{2i\pi\xi j\Delta x}$, on obtient

$$\frac{3\hat{u}^{n+1}(\xi)-4\hat{u}^n(\xi)+\hat{u}^{n-1}(\xi)}{2\Delta t}-\nu \hat{u}^{n+1}(\xi) \frac{e^{2i\pi\xi\Delta x}-2+e^{-2i\pi\xi\Delta x}}{\Delta x^2}=0$$

On pose $\mu = \frac{2\nu \Delta t}{\Delta x^2}$ et on utilise $e^{2i\pi\xi\Delta x}+e^{-2i\pi\xi\Delta x}-2 = 2\cos(2\pi\xi\Delta x)-2 = -4 \sin^2(\pi\xi\Delta x)$ :

$$(3+4\mu \sin^2(\pi\xi\Delta x))\hat{u}^{n+1}(\xi)-4\hat{u}^n(\xi)+\hat{u}^{n-1}(\xi)=0$$

En posant $\theta = \sin^2(\pi\xi\Delta x) \in [0,1]$ et en cherchant un facteur d'accroissmeent $a(\xi)$, on torouve qu'il vérifie l'équation :

$$(3+4\mu\theta)a^2(\xi)-4a(\xi)+1=0$$

Le discriminant vaut : $\Delta = 16 - 4 (3+4\mu\theta) = 4(1- 4\mu \theta)$.

Cas 1 : si $4\mu\theta \leq 1$, les racines sont réelles, $a_\pm(\xi) = \frac{2 \pm \sqrt{1-4\mu\theta}}{3+4\mu\theta}$

On a $a_+(\xi) a_-(\xi) = \frac{1}{3+4\mu\theta} \in [0,1/3]$. Comme $a_+(\xi)>0$, on en déduit $0<a_-(\xi)\leq a_+(\xi)$.

Pour $a_+(\xi)$ : On a $2 + \sqrt{1-4\mu\theta} \leq 3 < 3+4\mu\theta$, donc $0<a_-(\xi)\leq a_+(\xi) < 1$.


Cas 2 : Si $4\mu\theta > 1$, les racines sont complexes conjuguées $a_\pm(\xi) = \frac{2 \pm i \sqrt{4\mu\theta-1}}{3+4\mu\theta}$. On a :
$$|a_\pm(\xi)|^2 = \frac{4+(4\mu \theta-1)}{(3+4\mu\theta)^2} = \frac{1}{3+4\mu\theta} < 1$$

donc $|a_\pm(\xi)| < 1$.

Conclusion : Le schéma est inconditionnellement stable pour tout $\mu > 0$.


## Exercice 3 (5 points)

### 3.1

On utilise le schéma de Lax-Wendroff ci-dessous ($\sigma := V\Delta t/\Delta x$) pour résoudre l'équation d'advection avec $V > 0$ et CL périodiques : 

$$ u_j^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n)
             + \frac{\sigma^2}{2}(u_{j+1}^{n+1} - 2u_j^{n+1} + u_{j-1}^{n+1}). $$ 

On note $U^n := (u_j^n)_j$, montrer que ce schéma s'écrit $AU^{n+1} = W^n$ avec $A$ une matrice et $W$ un vecteur que l'on précisera. Montrer que cette matrice est inversible.

**Correction.** Grâce à la périodicité ($u_{-1} = u_{J-1}$, $u_J = u_0$), on a $J$ inconnues $U^n := (u_j^n)_{0 \le j \le J-1}$ ($J \ge 3$). En passant les termes en $n+1$ à gauche, le schéma s'écrit pour tout $j$

$$ -\frac{\sigma^2}{2} u_{j-1}^{n+1} + (1+\sigma^2)\, u_j^{n+1} - \frac{\sigma^2}{2} u_{j+1}^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n), $$

c'est-à-dire $AU^{n+1} = W^n$ avec

$$ A = \begin{pmatrix} 1+\sigma^2 & -\sigma^2/2 & & & -\sigma^2/2 \\ -\sigma^2/2 & 1+\sigma^2 & -\sigma^2/2 & & \\ & \ddots & \ddots & \ddots & \\ & & -\sigma^2/2 & 1+\sigma^2 & -\sigma^2/2 \\ -\sigma^2/2 & & & -\sigma^2/2 & 1+\sigma^2 \end{pmatrix}, \qquad W^n = \left( u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n) \right)_{j}. $$

La matrice $A$ est tridiagonale, avec en plus les deux coins $A_{0,J-1} = A_{J-1,0} = -\sigma^2/2$ dus à la périodicité.

*Inversibilité.* Soit $U$ tel que $AU = 0$ et $j_0$ un indice tel que $|u_{j_0}| = \max_j |u_j|$. La ligne $j_0$ donne $(1+\sigma^2)\, u_{j_0} = \frac{\sigma^2}{2}(u_{j_0+1} + u_{j_0-1})$, donc

$$ (1+\sigma^2)\,|u_{j_0}| \le \frac{\sigma^2}{2}\left( |u_{j_0+1}| + |u_{j_0-1}| \right) \le \sigma^2 |u_{j_0}|, $$

d'où $|u_{j_0}| \le 0$ et $U = 0$. La matrice carrée $A$ est donc injective, donc inversible, pour tout $\sigma$. (On peut aussi remarquer que $A$ est symétrique à diagonale strictement dominante, donc définie positive, ses valeurs propres étant $1 + \sigma^2(1 - \cos\varphi) \ge 1$.)

### 3.2

Dans la portion de code associé ci-dessous, quel est le rôle de l'expression `u = F \ w` et que représente `F` ?

```julia=
for n ∈ 1:N
    w = u - ...
    u = F \ w
end
```

**Correction.** `F` est la factorisation de Cholesky de la matrice $A$ (symétrique définie positive), calculée une seule fois avant la boucle en temps. Dans la boucle, `w` représente $W^n$ et l'expression `u = F \ w` résout le système linéaire $AU^{n+1} = W^n$ à l'aide de la factorisation (deux résolutions triangulaires, peu coûteuses grâce à la structure creuse) : `u` contient ensuite $U^{n+1}$. On évite ainsi d'inverser $A$ et de la refactoriser à chaque pas de temps.

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

**Correction.** La commande `spzeros(K, K)` crée une matrice creuse (*sparse*) $K \times K$ nulle : seuls les coefficients non nuls affectés ensuite (au plus $5$ par ligne pour le laplacien discret) sont stockés, ce qui représente une mémoire en $O(K)$ au lieu de $O(K^2)$ et permet d'utiliser des algorithmes d'algèbre linéaire creuse.

La ligne 11 concerne le voisin de gauche (ouest), dans la direction $x$ comme le voisin de droite (est) de la ligne 15 : on complète par

```julia
if kW ≠ 0 Δ[k, kW] = 1 / Δx^2 end
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

**Correction.** Soit $u \in \mathscr{C}^2(\overline{\Omega})$ solution forte et $v \in H^1(\Omega)$. On multiplie l'équation par $v$ et on intègre sur $\Omega$ :

$$ -\int_\Omega \Delta u\, v\,\mathrm{d}x + \int_\Omega u\, v\,\mathrm{d}x = \int_\Omega f v\,\mathrm{d}x. $$

La formule de Green

$$ \int_\Omega \Delta u\, v\,\mathrm{d}x = -\int_\Omega \nabla u \cdot \nabla v\,\mathrm{d}x + \int_\Gamma \frac{\partial u}{\partial n}\, v\,\mathrm{d}\sigma, $$

valable pour $u \in \mathscr{C}^2(\overline{\Omega})$ et $v \in \mathscr{C}^1(\overline{\Omega})$ (le bord est régulier), s'étend à $v \in H^1(\Omega)$ par densité de $\mathscr{C}^1(\overline{\Omega})$ dans $H^1(\Omega)$ et continuité de l'application trace $H^1(\Omega) \to L^2(\Gamma)$. Avec la condition de Neumann $\partial u/\partial n = g$ sur $\Gamma$, on obtient

$$ \int_\Omega \nabla u \cdot \nabla v\,\mathrm{d}x + \int_\Omega u\, v\,\mathrm{d}x = \int_\Omega f v\,\mathrm{d}x + \int_\Gamma g\, v\,\mathrm{d}\sigma. $$

Ainsi $u \in H^1(\Omega)$ (car $u \in \mathscr{C}^2(\overline{\Omega})$ et $\Omega$ est borné) est solution du problème faible $a(u,v) = \varphi \cdot v$ pour tout $v \in H^1(\Omega)$, avec

$$ a(u,v) = \int_\Omega \left( \nabla u \cdot \nabla v + u\, v \right) \mathrm{d}x, \qquad \varphi \cdot v = \int_\Omega f v\,\mathrm{d}x + \int_\Gamma g\, v\,\mathrm{d}\sigma. $$

La condition de Neumann, non imposée dans l'espace, apparaît dans la forme linéaire (condition *naturelle*).

### 4.2

Montrer qu'on a existence et unicité de solution faible.

**Correction.** On applique le théorème de Lax-Milgram dans l'espace de Hilbert $H^1(\Omega)$, muni de la norme $\|v\|_{H^1}^2 = \int_\Omega (|\nabla v|^2 + v^2)\,\mathrm{d}x$.

- *Continuité de $a$.* Par l'inégalité de Cauchy-Schwarz, $|a(u,v)| \le \|\nabla u\|_{L^2}\|\nabla v\|_{L^2} + \|u\|_{L^2}\|v\|_{L^2} \le \|u\|_{H^1}\|v\|_{H^1}$.
- *Coercivité de $a$.* $a(v,v) = \|v\|_{H^1}^2$ pour tout $v \in H^1(\Omega)$, avec constante $1$.
- *Continuité de $\varphi$.* Comme $\Omega$ est borné et $f$ continue sur $\overline{\Omega}$, $f \in L^2(\Omega)$ ; comme $\Gamma$ est compact et $g$ continue, $g \in L^2(\Gamma)$. Par le théorème de trace, il existe $C > 0$ tel que $\|v\|_{L^2(\Gamma)} \le C \|v\|_{H^1}$, donc
$$ |\varphi \cdot v| \le \|f\|_{L^2(\Omega)} \|v\|_{L^2(\Omega)} + \|g\|_{L^2(\Gamma)} \|v\|_{L^2(\Gamma)} \le \left( \|f\|_{L^2(\Omega)} + C \|g\|_{L^2(\Gamma)} \right) \|v\|_{H^1}. $$

Le théorème de Lax-Milgram donne l'existence et l'unicité de $u \in H^1(\Omega)$ solution faible. (La forme $a$ étant le produit scalaire de $H^1(\Omega)$, c'est aussi une conséquence du théorème de représentation de Riesz.)

### 4.3

Montrer que toute solution faible dans $\mathscr{C}^2(\overline{\Omega})$ est solution forte.

**Correction.** Soit $u \in \mathscr{C}^2(\overline{\Omega})$ solution faible.

*Équation dans $\Omega$.* Pour $v \in \mathscr{C}_c^\infty(\Omega) \subset H^1(\Omega)$ (nulle au bord), la formule de Green donne $\int_\Omega \nabla u \cdot \nabla v\,\mathrm{d}x = -\int_\Omega \Delta u\, v\,\mathrm{d}x$, et la formulation faible devient

$$ \int_\Omega \left( -\Delta u + u - f \right) v\,\mathrm{d}x = 0 \quad \text{pour tout } v \in \mathscr{C}_c^\infty(\Omega). $$

La fonction $-\Delta u + u - f$ étant continue, on en déduit $-\Delta u + u = f$ en tout point de $\Omega$ (lemme fondamental du calcul des variations), puis sur $\overline{\Omega}$ par continuité.

*Condition au bord.* Soit maintenant $v \in H^1(\Omega)$ quelconque. La formule de Green et l'équation précédente donnent

$$ a(u,v) = \int_\Omega \left( -\Delta u + u \right) v\,\mathrm{d}x + \int_\Gamma \frac{\partial u}{\partial n}\, v\,\mathrm{d}\sigma = \int_\Omega f v\,\mathrm{d}x + \int_\Gamma \frac{\partial u}{\partial n}\, v\,\mathrm{d}\sigma. $$

Comme $u$ est solution faible, $a(u,v) = \int_\Omega f v + \int_\Gamma g\, v$, d'où, avec $h := \partial u/\partial n - g$,

$$ \int_\Gamma h\, v\,\mathrm{d}\sigma = 0 \quad \text{pour tout } v \in H^1(\Omega). $$

La fonction $h$ est continue sur $\Gamma$ (car $u \in \mathscr{C}^2(\overline{\Omega})$, $g$ est continue et le bord est régulier). En prenant pour $v$ des polynômes, dont les restrictions à $\Gamma$ (compact) sont denses dans $\mathscr{C}(\Gamma)$ pour la norme uniforme (Stone-Weierstrass), on obtient $\int_\Gamma h\,\psi\,\mathrm{d}\sigma = 0$ pour tout $\psi \in \mathscr{C}(\Gamma)$, et $\psi = h$ donne $\int_\Gamma h^2\,\mathrm{d}\sigma = 0$, soit $h = 0$ : $\partial u/\partial n = g$ sur $\Gamma$. La fonction $u$ est donc solution forte.
