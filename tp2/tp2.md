![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4 - EDP1
# 2026-27
# **TP 2 - Équation de la chaleur**

On souhaite simuler l'évolution de la température dans une pièce en forme de L
- dont le chauffage se situe au centre de la pièce, à droite du coin rentrant,
- dont les parois sont isolantes,
- à l'exception d'une fenêtre située sur le mur du bas.

On rappelle que cette température suit l'EDP de la chaleur 2D :

$$ \begin{aligned}
  & \frac{\partial u}{\partial t} - \nu\Delta u = f,\quad t > 0,\quad (x,y) \in \Omega,\\
  & u(0,x,y) = u_0(x,y),\quad  (x,y) \in \Omega,
\end{aligned} $$

avec des conditions aux limites appropriées sur le bord du domaine.

**Le problème.** La pièce est $\Omega := ]0,1[^2 \setminus [0,\tfrac12]\times[\tfrac12,1]$ (un carré dont on a retiré le coin supérieur gauche). Le chauffage est la source $f = f_0 = 250$ sur le carré $[0.6,0.8]\times[0.4,0.6]$ et $f = 0$ ailleurs, avec $\nu = 1$. Sur les parois, la condition est de Neumann homogène (parois isolantes) :
$$ \frac{\partial u}{\partial n} = 0, $$
à l'exception de la fenêtre $\{(x,0)\,;\ \tfrac14 \le x \le \tfrac34\}$ sur le mur du bas, où la température est imposée (Dirichlet inhomogène) : $u = u_w = -1$. À l'instant initial, la pièce est à la température de la fenêtre : $u_0 \equiv u_w$.

![](chaleur_piece.png)

**Conventions numériques.**

- *Maillage en cellules.* On découpe le carré $]0,1[^2$ en $N \times N$ cellules carrées de côté $h = 1/N$, et l'on garde les cellules contenues dans $\Omega$ ; $N$ est un multiple de $20$ (pour que le coin retiré, le chauffage et la fenêtre soient formés de cellules entières). L'inconnue $u_k$ approche la température au centre de la cellule $k$. En Julia, la cellule $(j,i)$ est repérée par sa ligne $j$ (en $y$) et sa colonne $i$ (en $x$).
- *Volumes finis.* On intègre $-\nu\Delta u = f$ sur chaque cellule : le flux à travers une face commune à deux cellules voisines $k$ et $l$ est approché par $-\nu(u_l - u_k)/h$ ; un mur isolant ne donne aucun flux (le terme correspondant est simplement absent) ; une face de la fenêtre est à la température $u_w$, à distance $h/2$ du centre de la cellule, d'où le flux $\nu(u_k - u_w)/(h/2)$.
- *Pas de temps.* On pose $\sigma := \nu\Delta t/h^2$ (c'est $\sigma_x = \sigma_y$ du cours). Le temps final $t_f$ est fixé et $\Delta t = t_f/N_t$ avec $N_t$ le plus petit entier tel que $\nu\Delta t/h^2 \le \sigma$ : faire varier $\sigma$ ne change pas $t_f$.
- *Distance à l'état stationnaire.* On mesure $\|u\| := \big(h^2\sum_k u_k^2\big)^{1/2}$, la norme $L^2$ discrète.

Le code ci-dessous est à exécuter une fois en début de session.

```julia
using SparseArrays, LinearAlgebra, Plots

# Problem data

ν  = 1.0      # diffusion coefficient
uw = -1.0     # window temperature
f0 = 250.0    # heating power density
tf = 4.0      # final time

heater0 = ((0.6, 0.8), (0.4, 0.6))    # heater: x-range, y-range
umax = 10.0   # upper bound of the colour scale

# Geometry: Ω = ]0,1[² minus [0,1/2]×[1/2,1], cells of size h = 1/N.
# Cell (j, i): row j ↔ y, column i ↔ x. N must be a multiple of 20.
# Returns: mask (cells in Ω), win (cells touching the window), f (heating)

function lroom(N; heater = heater0)
    @assert N % 20 == 0
    mask = trues(N, N)
    mask[N÷2+1:N, 1:N÷2] .= false               # removed square
    win = falses(N, N)
    win[1, N÷4+1:3*N÷4] .= true                 # window: south wall, 1/4 ≤ x ≤ 3/4
    c = ((1:N) .- 0.5) / N                      # cell centres (in x and in y)
    (x0, x1), (y0, y1) = heater
    f = zeros(N, N)
    f[y0 .< c .< y1, x0 .< c .< x1] .= f0       # rows ↔ y, columns ↔ x
    @assert all(mask[f .> 0])                   # the heater must lie in Ω
    return mask, win, f
end

# Cell values (vector) → matrix with NaN outside Ω (NaN cells are not drawn)

field(u, mask) = (U = fill(NaN, size(mask)); U[mask] = u; U)

# Filled contour plot of the cell values u, with the heater drawn as a red square

function heatplot(u, mask; title = "", heater = heater0)
    N = size(mask, 1)
    xc = ((1:N) .- 0.5) / N
    p = contourf(xc, xc, field(u, mask); aspect_ratio = :equal, c = :inferno, levels = range(uw, umax, length = 33),
                 clims = (uw, umax), linewidth = 0, xlims = (0, 1), ylims = (0, 1), title = title)
    (x0, x1), (y0, y1) = heater
    plot!(p, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]; color = :red, linewidth = 2, label = false)
end
```

## Exercice 1

On commence par le problème stationnaire (Poisson) : trouver $u$ telle que $-\nu\Delta u = f$ dans $\Omega$ avec les conditions aux limites ci-dessus. Le code suivant assemble la matrice $A$ et le second membre $b$ du système linéaire $Au = b$.

```julia
# Finite volume assembly of A ≃ -νΔ and b (heating + Dirichlet data).
# Cells are numbered columnwise: idx[mask] = 1:K

function assemble(mask, win, f, h; ν = 1.0, uw = 0.0)
    K = count(mask)
    idx = zeros(Int, size(mask))
    idx[mask] = 1:K
    rows, cols, vals = Int[], Int[], Float64[]
    b = f[mask]                                 # new array
    for c ∈ findall(mask)                       # c = CartesianIndex(j, i)
        k = idx[c]
        d = 0.0
        for δ ∈ (CartesianIndex(0, -1), CartesianIndex(0, 1), CartesianIndex(-1, 0), CartesianIndex(1, 0))
            n = c + δ
            if checkbounds(Bool, mask, n) && mask[n]    # neighbour in Ω
                push!(rows, k); push!(cols, idx[n]); push!(vals, -ν / h^2)
                d += ν / h^2
            end                                          # else: insulating wall, no flux
        end
        if win[c]                                        # Dirichlet face at distance h/2
            d += 2ν / h^2
            b[k] += 2ν / h^2 * uw
        end
        push!(rows, k); push!(cols, k); push!(vals, d)
    end
    return sparse(rows, cols, vals, K, K), b, idx
end

N = 40
h = 1 / N
mask, win, f = lroom(N)
A, b, idx = assemble(mask, win, f, h; ν = ν, uw = uw)
K = size(A, 1)

# Stationary solution: -νΔu = f (A is symmetric positive definite)

ustat = cholesky(Symmetric(A)) \ b

# Plot (move the heater with lroom(N; heater = ...), or the window in `lroom`)

heatplot(ustat, mask)
```

**1.1** Montrer que, pour une cellule $k$, la ligne $k$ du système s'écrit
$$ \frac{\nu}{h^2}\Big( n_k u_k - \sum_{l \sim k} u_l \Big) + \frac{2\nu}{h^2}\,u_k\,\mathbf{1}_{k \in W} = f_k + \frac{2\nu}{h^2}\,u_w\,\mathbf{1}_{k \in W}, $$
où $l \sim k$ désigne les cellules voisines de $k$ *dans* $\Omega$, $n_k$ leur nombre et $W$ l'ensemble des cellules de la fenêtre. Écrire $A$ à l'aide du laplacien du graphe des cellules, et vérifier que $A$ est symétrique.

**1.2** Montrer que $A$ est symétrique définie positive. Quel est le rôle de la fenêtre ? Que se passerait-il sans fenêtre (Neumann sur tout le bord) ?

**1.3** Exécuter le code et tracer la solution stationnaire. Vérifier numériquement le bilan d'énergie : la puissance de chauffe $h^2\sum_k f_k$ est égale au flux sortant par la fenêtre $h^2\sum_{k\in W}\frac{2\nu}{h^2}(u_k - u_w)$. Démontrer ce bilan au niveau discret.

```julia
# Energy balance: heating power = flux through the window

power = h^2 * sum(f)
flux  = h^2 * sum(2ν / h^2 * (ustat[idx[win]] .- uw))
@show power flux
```

**1.4** Déplacer le chauffage (argument `heater` de `lroom`) et la fenêtre, et commenter. Par exemple, comparer le chauffage en $[0.6,0.8]\times[0.4,0.6]$ et en $[0.6,0.8]\times[0,0.2]$ : que devient la température maximale, à puissance de chauffe égale ? Le bilan d'énergie change-t-il ?

**1.5** *Validation.* Pour vérifier le code d'assemblage, on prend le carré $]0,1[^2$ avec Dirichlet $u = 0$ sur tout le mur du bas et Neumann ailleurs, et $f(x,y) = \nu\frac{\pi^2}{4}\sin\big(\frac{\pi y}{2}\big)$. Donner la solution exacte, puis calculer l'erreur maximale aux centres des cellules pour $N = 20, 40, 80, 160$ et en déduire l'ordre observé.

```julia
# Validation: unit square, u = 0 on the whole south wall, Neumann elsewhere,
# f = (π²/4) sin(πy/2) (ν = 1); exact solution u = sin(πy/2)

function test_square(N)
    h = 1 / N
    mask = trues(N, N)
    win = falses(N, N)
    win[1, :] .= true
    y = ((1:N) .- 0.5) * h
    f = repeat(π^2 / 4 * sin.(π * y / 2), 1, N)
    A, b, _ = assemble(mask, win, f, h)         # ν = 1, uw = 0
    u = cholesky(Symmetric(A)) \ b
    return maximum(abs.(u - vec(repeat(sin.(π * y / 2), 1, N))))
end

for N ∈ (20, 40, 80, 160)
    println("N = $N: max error = $(test_square(N))")
end
```

**1.6** *Convergence sur la pièce en L.* Calculer $\int_\Omega u \simeq h^2\sum_k u_k$ et $\max u$ pour $N = 20, 40, 80, 160$. Quel ordre observe-t-on ? Pourquoi est-il différent de celui de la question précédente ?

```julia
# Convergence on the L-shaped room: ∫u ≃ h² Σ u_k and max u

function steady(N)
    h = 1 / N
    mask, win, f = lroom(N)
    A, b, _ = assemble(mask, win, f, h; ν = ν, uw = uw)
    u = cholesky(Symmetric(A)) \ b
    return h^2 * sum(u), maximum(u)
end

for N ∈ (20, 40, 80, 160)
    I, m = steady(N)
    println("N = $N: ∫u = $I, max u = $m")
end
```

## Exercice 2

On calcule maintenant l'évolution temporelle de la température. On utilise le $\theta$-schéma
$$ \frac{u^{n+1} - u^n}{\Delta t} + A\big(\theta u^{n+1} + (1-\theta)u^n\big) = b, $$
c'est-à-dire
$$ (I + \theta\Delta t\,A)\,u^{n+1} = \big(I - (1-\theta)\Delta t\,A\big)u^n + \Delta t\,b, $$
où $A$ et $b$ sont ceux de l'exercice 1 ($b$ contient à la fois le chauffage et la donnée de Dirichlet). Le cas $\theta = 0$ est le schéma explicite, $\theta = 1$ le schéma d'Euler implicite et $\theta = \tfrac12$ le schéma de Crank-Nicolson.

**2.1** Pour une cellule intérieure (ses quatre voisines sont dans $\Omega$, pas de terme source), écrire le schéma de Crank-Nicolson ($\theta = \frac12$) et le comparer à
$$ \begin{aligned}
  & u_k^{n+1} - \frac{\sigma_x}{2}(u_{kW}^{n+1} - 2u_k^{n+1} + u_{kE}^{n+1})
              - \frac{\sigma_y}{2}(u_{kN}^{n+1} - 2u_k^{n+1} + u_{kS}^{n+1})\\
  & = u^n_k   + \frac{\sigma_x}{2}(u_{kW}^n - 2u_k^n + u_{kE}^{n})
              + \frac{\sigma_y}{2}(u_{kN}^n - 2u_k^n + u_{kS}^n) + \Delta t\, f_k,
\end{aligned} $$
où $\sigma_x = \nu\Delta t/\Delta x^2$ et $\sigma_y = \nu\Delta t/\Delta y^2$. Que deviennent les termes de bord (mur isolant, fenêtre) dans l'écriture matricielle ?

**2.2** *Stabilité.* Les valeurs propres de $A$ sont notées $0 < \lambda_1 \le \dots \le \lambda_K$. Exprimer le facteur d'amplification du $\theta$-schéma sur chaque mode propre. Montrer que le schéma de Crank-Nicolson est stable pour tout $\Delta t > 0$. Montrer que $\lambda_K \le 8\nu/h^2$ (lemme de Gershgorin) et en déduire une condition suffisante de stabilité pour le schéma explicite en fonction de $\sigma$.

**2.3** Compléter la fonction `evolve` (`# TBC`) pour programmer le $\theta$-schéma, puis exécuter Crank-Nicolson avec $\sigma = 2$ et $t_f = 4$. Que se passe-t-il au bout d'un temps long ? Comparer avec la solution stationnaire de l'exercice 1 en traçant $\|u(t) - u_{\mathrm{stat}}\|$ en fonction de $t$ (échelle logarithmique), pour Crank-Nicolson avec $\sigma = 0.2,\ 2,\ 20,\ 200$ et Euler implicite avec $\sigma = 20$. Quelle est la vitesse de convergence observée ? La relier à $\lambda_1$ (obtenue par `eigvals(Symmetric(Matrix(A)), 1:1)`).

```julia
# θ-scheme: (I + θΔt A) uⁿ⁺¹ = (I - (1-θ)Δt A) uⁿ + Δt b.
# θ = 0: explicit, θ = 1/2: Crank-Nicolson, θ = 1: implicit Euler.
# tf is fixed; the number of steps is deduced from σ = νΔt/h²

Id = sparse(1.0I, K, K)
nsteps(σ, T) = ceil(Int, T * ν / (σ * h^2) - 1e-9)
dist(u) = h * norm(u - ustat)                   # distance to the stationary solution

function evolve(θ, σ; T = tf, callback = (n, t, u) -> nothing)
    Nt = nsteps(σ, T)
    Δt = T / Nt
    F = nothing                                 # TBC: factorisation of I + θΔt A
    E = nothing                                 # TBC: I - (1-θ)Δt A
    u = fill(uw, K)                             # u0 = uw
    err = [dist(u)]
    callback(0, 0.0, u)
    for n ∈ 1:Nt
        u = u                                   # TBC: uⁿ⁺¹ from uⁿ
        push!(err, dist(u))
        callback(n, n * Δt, u)
    end
    return u, err, Δt
end

# Crank-Nicolson, σ = 2: distance to the stationary solution along time

u, err, Δt = evolve(0.5, 2.0)
heatplot(u, mask)

p = plot(; yscale = :log10, xlabel = "t", ylabel = "‖u(t) - u_stat‖")
for (θ, σ, name) ∈ ((0.5, 0.2, "CN σ=0.2"), (0.5, 2.0, "CN σ=2"), (0.5, 20.0, "CN σ=20"), (0.5, 200.0, "CN σ=200"), (1.0, 20.0, "IE σ=20"))
    _, e, dt = evolve(θ, σ)
    plot!(p, (0:length(e)-1) * dt, e; label = name)
end
display(p)

# Slowest mode: smallest eigenvalue of A (dense, K ≃ 1200)

λ1 = eigvals(Symmetric(Matrix(A)), 1:1)[1]
```

**2.4** *Schéma explicite.* Exécuter `evolve(0.0, σ)` pour $\sigma = 0.2$ et $\sigma = 0.26$ (avec un temps final plus court si besoin). Combien de pas de temps faut-il à $t_f = 4$ pour le schéma explicite à la limite de stabilité, pour Crank-Nicolson avec $\sigma = 2$ et avec $\sigma = 20$ ? Commenter le coût.

**2.5** *Ordre en temps.* Calculer une solution de référence à $t = 0.5$ avec Crank-Nicolson et $\sigma = 0.02$, puis l'erreur de Crank-Nicolson et d'Euler implicite pour $\sigma = 0.2,\ 0.4,\ \dots,\ 6.4$. Quels ordres observe-t-on ?

```julia
# Order in time at t = 0.5: reference = Crank-Nicolson with σ = 0.02

uref, _, _ = evolve(0.5, 0.02; T = 0.5)
for θ ∈ (0.5, 1.0), σ ∈ (0.2, 0.4, 0.8, 1.6, 3.2, 6.4)
    uθ, _, _ = evolve(θ, σ; T = 0.5)
    println("θ = $θ, σ = $σ: error = $(h * norm(uθ - uref))")
end
```

**2.6** *Pas de temps très grand.* Pour $\sigma = 200$, comparer Crank-Nicolson et Euler implicite : erreur à $t_f = 4$, et monotonie en temps. Pourquoi la température doit-elle être croissante en tout point (rappel : $f \ge 0$ et $u_0 \equiv u_w$) ? Ce que l'on observe est-il compatible avec cette propriété ? Expliquer avec $a(\lambda)$ quand $\Delta t\lambda \to \infty$.

```julia
# Monotonicity: smallest increment uⁿ⁺¹_k - uⁿ_k over all cells and steps

function min_increment(θ, σ)
    uprev = fill(uw, K)
    m = Inf
    evolve(θ, σ; callback = (n, t, u) -> begin
        n > 0 && (m = min(m, minimum(u - uprev)))
        uprev = copy(u)
    end)
    return m
end

for (θ, name) ∈ ((0.5, "Crank-Nicolson"), (1.0, "implicit Euler"))
    _, e, _ = evolve(θ, 200.0)
    println("$name, σ = 200: error at tf = $(e[end]), min increment = $(min_increment(θ, 200.0))")
end

# Explicit scheme just above and just below σ = 1/4 (T shortened to 400 steps)

for σ ∈ (0.2, 0.26)
    ux, _, _ = evolve(0.0, σ; T = 400 * σ * h^2 / ν)
    println("explicit, σ = $σ: max|u| after 400 steps = $(maximum(abs, ux))")
end
```
