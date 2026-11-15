![PNS](https://raw.githubusercontent.com/pns-mam/edp1/master/logo-pns.png)

## MAM4 - EDP1
# 2026-27
# **TP 1 - Équation d'advection**

On résout numériquement l'équation d'advection sur $\Omega := ]0,1[$ avec conditions aux limites périodiques ($V > 0$) :

$$ \left\lbrace \begin{aligned}
& \frac{\partial u}{\partial t}(t,x) + V\frac{\partial u}{\partial x}(t,x) = 0,\quad x \in \Omega,\quad t > 0,\\
& u(t, 0) = u(t, 1),\quad t \ge 0,\\
& u(0, x) = u_0(x),\quad x \in \Omega.
\end{aligned} \right. $$

On suppose la donnée initiale $u_0$ $1$-périodique. On vérifie aisément que la solution exacte est donnée par $u(t,x) = u_0(x-Vt)$. On pose $t_n = n\Delta t$ et $x_j = j\Delta x$, et on note $u_j^n \simeq u(t_n,x_j)$.

**Conventions numériques.**

- *Maillage périodique.* On utilise $J$ points $x_j = j\Delta x$, $j = 0,\dots,J-1$, avec $\Delta x = 1/J$ (le point $x = 1$ n'est pas stocké : il est identifié à $x = 0$). La périodicité s'écrit $u_{-1} := u_{J-1}$ et $u_J := u_0$ ; en Julia, `circshift(u, 1)` donne $(u_{j-1})_j$ et `circshift(u, -1)` donne $(u_{j+1})_j$.
- *Temps final fixé.* On fixe le temps final $t_f$ et le nombre de pas de temps $N$ est déduit du nombre de Courant $\sigma = V\Delta t/\Delta x$ : on prend $\Delta t = t_f/N$ avec $N$ le plus petit entier tel que $\sigma$ effectif $\le \sigma$ demandé. Ainsi, **faire varier $\sigma$ ne change pas $t_f$** et les résultats sont comparables. Avec $V = 1$ et $t_f = 1$, la solution exacte à l'instant final est $u_0$ (une révolution).
- *Erreur.* On mesure l'erreur dans la norme $L^2$ discrète $\|u\|_2 = \big(\frac{1}{J}\sum_j u_j^2\big)^{1/2}$.
- *Deux données initiales* : une donnée lisse `u0_smooth` et une donnée discontinue `u0_disc` (créneau).

Le code ci-dessous (utilisé dans tous les exercices) est fourni ; il est à exécuter une fois en début de session. Les fonctions `solve`, `exact`, `err_L2` et `plot_result` n'ont pas à être modifiées. Dans les exercices, les parties marquées `# TBC` (*to be completed*) sont à compléter.

```julia
using LinearAlgebra, SparseArrays, Plots

# Problem data

V  = 1.0                                   # advection speed
tf = 1.0                                   # final time (one revolution, so u(tf, ·) = u0)

u0_smooth(x) = exp(-10 * (1 - cos(2π * (x - 0.25))))   # smooth, 1-periodic
u0_disc(x)   = Float64(0.2 < mod(x, 1) < 0.3)          # discontinuous, 1-periodic

datasets = [(u0_smooth, "smooth data"), (u0_disc, "discontinuous data")]

# Periodic grid: J points x_j = jΔx, j = 0, …, J-1, with Δx = 1/J

grid(J) = (0:J-1) / J

# Number of time steps N for a Courant number ≤ σ (tf is fixed, Δt = tf / N)

nsteps(J, σ) = ceil(Int, V * tf * J / σ - 1e-9)

# A scheme is a function (J, σ) -> (u -> next u).
# `solve` returns the numerical solution at tf, N and the effective Courant number

function solve(scheme, u0, J, σ)
    N = nsteps(J, σ)
    σe = V * (tf / N) * J
    step = scheme(J, σe)
    u = u0.(grid(J))
    for n in 1:N
        u = step(u)
    end
    return u, N, σe
end

exact(u0, J) = u0.(grid(J) .- V * tf)
err_L2(u, uex) = norm(u - uex) / sqrt(length(u))       # discrete L² norm on [0, 1)

function plot_result(x, u, uex; title = "")
    p1 = plot(x, u; color = :red, lw = 4, label = "numerical", ylabel = "u", title = title)
    plot!(p1, x, uex; color = :black, lw = 2, label = "exact")
    p2 = plot(x, u - uex; legend = false, xlabel = "x", ylabel = "error")
    return plot(p1, p2; layout = (2, 1), size = (700, 600))
end

# Run one case: plot the result and print the L² error

function run_case(scheme, u0, J, σ; name = "")
    u, N, σe = solve(scheme, u0, J, σ)
    display(plot_result(grid(J), u, exact(u0, J); title = "$name, σ = $(round(σe, digits = 3)), N = $N"))
    println("$name: J = $J, N = $N, σ = $σe, L² error = ", err_L2(u, exact(u0, J)))
end
```

## Exercice 1. Schéma décentré amont

On cherche à approcher la solution par le schéma *décentré amont* (*cf.* $V>0$) :

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + V\frac{u_j^n-u_{j-1}^n}{\Delta x} = 0,\quad j = 0,\dots,J-1,\quad n \geq 0, $$

soit, en réarrangeant les termes,

$$ u_j^{n+1} = u_j^n - \sigma(u_j^n-u_{j-1}^n), $$

où $\sigma = V\Delta t/\Delta x$ est le [*nombre de Courant*](https://fr.wikipedia.org/wiki/Nombre_de_Courant), avec $u_j^0 = u_0(x_j)$.

```julia
# Upwind scheme (periodic): u_j^{n+1} = u_j^n - σ (u_j^n - u_{j-1}^n)

upwind(J, σ) = u -> u - σ * (u - circshift(u, 1))      # circshift(u, 1)[j] = u[j-1]

for (u0, data) in datasets
    run_case(upwind, u0, 200, 0.8; name = "upwind, $data")
end
```

**1.1** Exécuter le code, qui teste les deux données initiales (`u0_smooth` et `u0_disc`). Comparer les deux solutions numériques à la solution exacte.

**1.2** Augmenter progressivement $\sigma$ (par exemple $\sigma \in \{0.5, 0.9, 1, 1.05, 1.1\}$) et observer le résultat. Quelle est la valeur critique ? Que se passe-t-il pour $\sigma = 1$ ? Retrouver cette valeur par l'analyse de stabilité de Von Neumann (facteur d'amplification du schéma) et par un argument de domaine de dépendance (condition CFL).

**1.3** À $\sigma$ fixé, diminuer $J$ (par exemple $J \in \{50, 100, 200, 400\}$) et constater que la solution numérique est amortie au fil des itérations en temps (phénomène de *diffusion numérique*). Comment cet amortissement dépend-il de $J$ et de $\sigma$ ?

## Exercice 2. Schéma de Lax-Wendroff

En suivant le modèle précédent, implémenter le schéma

$$ u_j^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n) + \frac{\sigma^2}{2}(u_{j+1}^n - 2u_j^n + u_{j-1}^n),\quad j = 0,\dots,J-1,\quad n \geq 0. $$

```julia
# Explicit Lax-Wendroff scheme

lax_wendroff(J, σ) = u -> zero(u)                      # TBC (circshift(u, -1)[j] = u[j+1], circshift(u, 1)[j] = u[j-1])

for (u0, data) in datasets
    run_case(lax_wendroff, u0, 200, 0.8; name = "Lax-Wendroff, $data")
end
```

**2.1** Tester le schéma pour les deux données initiales, avec $\sigma = 0.8$. Quelles différences qualitatives observe-t-on avec le schéma décentré amont ? Quelle donnée initiale met le mieux en évidence ces différences ?

**2.2** Augmenter progressivement $\sigma$ et observer. Quelle est la valeur critique ?

**2.3** Que peut-on dire de la diffusion numérique observée dans le cas du schéma décentré amont ? Le schéma de Lax-Wendroff est-il diffusif ? Quelle est la nature de l'erreur dominante ?

## Exercice 3. Schéma de Lax-Wendroff implicite

On modifie le schéma précédent en évaluant le terme de diffusion au temps $n+1$ :

$$ u_j^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n) + \frac{\sigma^2}{2}(u_{j+1}^{n+1} - 2u_j^{n+1} + u_{j-1}^{n+1}),\quad j = 0,\dots,J-1,\quad n \geq 0. $$

En introduisant $W^n := (u_j^n - (\sigma/2)(u_{j+1}^n - u_{j-1}^n))_j$ et $U^n := (u_j^n)_j$, le schéma s'écrit $AU^{n+1} = W^n$. La détermination de $U^{n+1}$ à partir de $U^n$ nécessite la résolution d'un système linéaire : on exploitera la structure creuse de $A$ et on la factorisera une seule fois, en dehors de la boucle en temps.

```julia
# Implicit Lax-Wendroff scheme: A U^{n+1} = W^n

function lax_wendroff_implicit(J, σ)
    A = sparse(1.0I, J, J)                  # TBC: sparse periodic matrix A (tridiagonal + corners)
    F = cholesky(A)                         # factorisation, computed once, outside of the time loop
    return u -> zero(u)                     # TBC: solve A U^{n+1} = W^n with F (compute W^n first)
end

J = 400
for (u0, data) in datasets, σ in [0.5, 1, 2, 5, 10, 20]
    u, N, σe = solve(lax_wendroff_implicit, u0, J, σ)
    println("$data, σ = $σ, N = $N, L² error: ", err_L2(u, exact(u0, J)))
end

for (u0, data) in datasets
    run_case(lax_wendroff_implicit, u0, J, 5.0; name = "implicit Lax-Wendroff, $data")
end
```

**3.1** Écrire la matrice $A$ ; justifier qu'elle est symétrique définie positive et préciser ses coins (périodicité). Compléter le code.

**3.2** Que constate-t-on lorsque $\sigma$ dépasse $1$ ? Calculer le facteur d'amplification du schéma et conclure sur sa stabilité.

**3.3** Comment évolue l'erreur lorsque $\sigma$ augmente ? Pourquoi peut-on avoir à la fois stabilité et mauvaise précision ?

## Exercice 4. Convergence et coût

**4.1 Ordre de convergence.** Compléter le code ci-dessous pour calculer l'erreur $L^2$ des trois schémas sur les maillages $J = 50, 100,\dots, 1600$ à $\sigma = 0.8$, et l'ordre observé (pente de la droite de régression de $\log(\text{erreur})$ en fonction de $\log \Delta x$). Commenter les ordres obtenus pour la donnée lisse puis pour la donnée discontinue.

```julia
Js = [50, 100, 200, 400, 800, 1600]
schemes = [("upwind", upwind), ("Lax-Wendroff", lax_wendroff), ("implicit LW", lax_wendroff_implicit)]

# L² errors on the successive grids (Courant number σ fixed)

function errors(scheme, u0; σ = 0.8)
    errs = Float64[]
    for J in Js
        u, = solve(scheme, u0, J, σ)
        push!(errs, 0.0)                     # TBC: L² error between u and the exact solution
    end
    return errs
end

# Observed order: slope of the least squares line through the points (log Δx, log error)

function order(errs)
    A = [ones(length(Js)) log.(1 ./ Js)]
    return 0.0                            # TBC: second coefficient of A \ log.(errs)
end

for (u0, name) in datasets
    p = plot(; xscale = :log10, yscale = :log10, xlabel = "Δx", ylabel = "L² error", title = name, legend = :bottomright)
    for (label, scheme) in schemes
        errs = errors(scheme, u0)
        plot!(p, 1 ./ Js, errs; marker = :o, label = "$label (order $(round(order(errs), digits = 2)))")
        println(name, ", ", label, ": observed order ", round(order(errs), digits = 2))
    end
    display(p)
end
```

**4.2 Coût.** Le code suivant mesure, à $J = 400$ et pour la donnée lisse, le temps de calcul, le nombre de pas de temps et l'erreur. Discuter les performances relatives des schémas explicite et implicite. Pour atteindre une erreur $L^2$ donnée (par exemple $10^{-2}$), quel schéma est le moins coûteux ? On pourra adapter le code en faisant varier $J$.

```julia
# Computing time, number of steps and error

function timing(scheme, u0, J, σ)
    solve(scheme, u0, J, σ)                                        # compilation
    t = minimum(@elapsed(solve(scheme, u0, J, σ)) for _ in 1:3)
    u, N, σe = solve(scheme, u0, J, σ)
    return t, N, err_L2(u, exact(u0, J))
end

J = 400
for (label, scheme, σ) in [("LW", lax_wendroff, 0.9), ("implicit LW", lax_wendroff_implicit, 0.5),
                           ("implicit LW", lax_wendroff_implicit, 1.0), ("implicit LW", lax_wendroff_implicit, 2.0),
                           ("implicit LW", lax_wendroff_implicit, 5.0)]
    t, N, e = timing(scheme, u0_smooth, J, σ)
    println(rpad(label, 12), "σ = $σ, N = $N, time = $(round(1e3 * t, digits = 2)) ms, L² error = $e")
end
```

## Exercice 5 (bonus). Équation équivalente

**5.1** Implémenter le schéma explicite centré $u_j^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n)$ et le tester avec les deux données initiales pour plusieurs valeurs de $\sigma$. Que constate-t-on ? Expliquer par l'analyse de Von Neumann.

```julia
# Explicit centered scheme

centered(J, σ) = u -> zero(u)                        # TBC: explicit Euler in time, centered in space

J = 200
for (u0, data) in datasets, σ in [0.1, 0.5, 0.9]
    u, N, σe = solve(centered, u0, J, σ)
    println("$data, σ = $σ, N = $N, max |u| = ", maximum(abs.(u)))
end
```

**5.2** Par développement de Taylor, déterminer l'équation équivalente du schéma décentré amont, sous la forme $u_t + Vu_x = \nu\, u_{xx} + \dots$ Pour $u_0(x) = \sin(2\pi x)$, la solution de $u_t = \nu u_{xx}$ est amortie d'un facteur $e^{-\nu(2\pi)^2 t}$. Compléter le code pour comparer l'amortissement observé à l'amortissement prédit.

```julia
# Damping of u0(x) = sin(2πx) by the upwind scheme

u0_sin(x) = sin(2π * x)

for J in [50, 100, 200], σ in [0.5, 0.9]
    u, N, σe = solve(upwind, u0_sin, J, σ)
    ν = 0.0                                    # TBC: numerical diffusion coefficient of the upwind scheme
    println("J = $J, σ = $σ: max u = ", maximum(u), ", predicted: ", exp(-ν * (2π)^2 * tf))
end
```

**5.3** Écrire l'équation équivalente du schéma de Lax-Wendroff. Quel terme domine ? Interpréter les oscillations observées à l'exercice 2.

## Annexe. Animations

La fonction suivante produit un GIF de l'évolution de la solution numérique (en rouge) et de la solution exacte (en noir), un panneau par schéma. Elle peut servir à visualiser les phénomènes observés.

```julia
# Animation: `cases` is a vector of (name, scheme, initial data), one panel per case

function animate(cases, J, σ; file = "evolution.gif", every = 4, fps = 15)
    x = grid(J)
    N = nsteps(J, σ)
    σe = V * (tf / N) * J
    steps = [scheme(J, σe) for (_, scheme, _) in cases]
    us = [u0.(x) for (_, _, u0) in cases]
    anim = Animation()
    for n in 0:N
        if n % every == 0
            panels = map(enumerate(cases)) do (k, (name, _, u0))
                p = plot(x, us[k]; color = :red, lw = 4, label = "numerical", ylims = (-0.4, 1.4), title = "$name, σ = $(round(σe, digits = 2))")
                plot!(p, x, u0.(x .- V * n * tf / N); color = :black, lw = 2, label = "exact")
            end
            frame(anim, plot(panels...; layout = (length(cases), 1), size = (700, 300 * length(cases))))
        end
        us = [steps[k](us[k]) for k in eachindex(cases)]
    end
    gif(anim, file; fps = fps)
end
```

```julia
animate([("upwind", upwind, u0_disc)], 200, 0.8; file = "evolution.gif")
```
