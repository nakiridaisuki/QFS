#set page(
  paper: "a4",
)
#show title: set text(size: 24pt)
#show title: set align(center)
#show title: set block(below: 1em)
#title[
  CA2026 Final Project\
  Adaptive Grid Fluid Simulator
]

#set text(
  font: ("Times New Roman", "New SimSun"),
  size: 12pt,
  lang: "en",
  region: "us",
)
#show heading.where(level: 1): set text(size: 18pt, weight: "bold")
#show heading.where(level: 2): set text(size: 16pt, weight: "bold")
#show heading.where(level: 3): set text(size: 14pt, weight: "bold")
#show heading.where(level: 4): set text(size: 12pt, weight: "bold")
#show selector.or(
  heading.where(level: 1),
  heading.where(level: 2),
  heading.where(level: 4),
): set heading(numbering: (..nums) => {
  let n = nums.pos()
  if n.len() <= 2 {
    numbering("I.1.", ..n)
  } else {
    numbering("a.", n.last())
  }
})
#show raw.where(block: true): it => block(
  fill: luma(240),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
  it,
)


#set align(center)

Author: 113550058 蔡昀呈

#set align(left)


Outline
- Introduction
  - Motivation
  - Fluid Simulation
  - What does this paper done
  - What did I do

- Implementation
  - MAC / FLIP
  - Quad Tree Simulator
  - Optimizations

- Conclude

= Introduction

== Motivation
A few month ago, I watched some captivating fluid simulator videos on Youtube (#link("https://youtu.be/rSKMYc1CQHE?si=LpY_0M06A26NAMnJ")[like this guy])
and some #link("https://youtu.be/Mh2y2Z6Iy0U?si=DDcJg1_ybWIZeHaR")[paper introduction] about fluid simulation.
I found the simulations incredibly cool and want to build one myself some day. This project is the result of that aspiration.

== Fluid Simulation
To simulate fluid, we need to solve the Navier-Stokes equations:
$
  rho ((partial arrow(u)) / (partial t) + (arrow(u) dot nabla) arrow(u)) = -nabla p + rho arrow(g) + mu nabla^2 arrow(u)
$
where $rho$, $arrow(u)$, $p$, and $arrow(g)$ are density, velocity, pressure, and gravity respectively.

In this project, we want to simulate water, so we set $rho=1$ and ignore the viscosity term $mu nabla^2 arrow(u)$.
We also enforce the incompressibility constraint (divergence-free condition). Under these assumptions, the equations simplify to:
$
  (partial arrow(u)) / (partial t) = -nabla p + arrow(g)\
  nabla dot arrow(u) = 0
$
To solve these equations numerically, we can advance the fluid by a tiny time step and solve for the pressure that makes the velocity divergence-free in the next step.

So, how do we solve for the pressure? Assume the intermediate velocity before pressure projection (which already incorporates advection and gravity) is $arrow(u)^*$. The projection step is defined as:
$
      & (arrow(u)_"next" - arrow(u)^*)/ (Delta t) = - nabla p \
  ==> & arrow(u)_"next" = arrow(u)^* - Delta t nabla p
$
Then, by taking the divergence on both sides:
$
      & nabla dot arrow(u)_"next" = nabla dot arrow(u)^* - Delta t nabla^2 p \
  ==> & 0 = nabla dot arrow(u)^* - Delta t nabla^2 p \
  ==> & Delta t nabla^2 p = nabla dot arrow(u)^*
$

To solve the pressure on a discrete grid like a MAC grid, we can turn this Poisson equation into a system of linear equations. More implementation details will be provided in the following sections.

However, the traditional MAC method requires discretizing the entire simulation domain---including deep water regions---using a uniform grid.
In 3D scenarios, this leads to a memory footprint that scales cubically as $O(N^3)$ with the grid resolution N,
posing a severe bottleneck as the simulation scale expands or the resolution increases.

In gaming or visual effects, we don't really need to simulate deep water regions in high detail;
instead, we only focus on surface phenomena like splashes or waves.
In this situation, using an adaptive grid near the surface can help us resolve this bottleneck.

== What does this paper done
To simulate details of surface phenomena, there are three things we need to do:
1. Find the surface
2. Subdivide grid around surface
3. Calculate the pressure on adaptive grid correctly
Let's talk about how to do that.

=== 1. Find the surface
To record the data of surface, author record a value $Phi$ at cell center.
Represent the level set, meaning the minimum distance from cell center to surface.
This is a SDF (signed distance field), less 0 means under the water.

It's easy to calculate precise $Phi$ when you add water drop like a circle into simulation zone.
But how to keep track $Phi$ every frame?
After semi-Lagrangian advect, $Phi$ will become unpricise.
We can use following two process to correct $Phi$: \
a. Get estimate position of $Phi = 0$ \
b. Use this estimate value to solve Eikonal equation $|nabla Phi| = 1$ \

For process a, we can check all adjacent cell if their $Phi$ have difference sign, then the surface must cross them.
To calculate estimate $Phi$, use this formula:
$
  Phi_"i_new" = - Phi_i / (Phi_i - Phi_j)
$
and then mark this cell as updated.

For process b, we can start from updated cell with minimum $Phi_"i_new"$,
and update it's neighbors using Eikonal equation to get corrected $Phi$.
Just work like Dijkstra algorithm, and this method called Fast Marching.
After this two process, we have new, precise estimate $Phi$ SDF.

=== 2. Subdivide grid
Author use two metrics to decide if a grid need to be subdivide:\
a. $Phi$ < node size\
b. Size function > 1 / node size\
If $Phi$ < node size means surface across this node, so it might need to be subdivide. \
And author use following size function to calculate how many details inside a node:
$
  S = gamma_Phi|nabla dot nabla (Phi + Phi_"solid")| + gamma_u sqrt(Sigma_i ((partial u_i )/ (partial x_i))^2)
$
where $gamma_Phi = 4$ and $gamma_u = 3$ and I don's have solid in simulation zone, therefore
$
  S = 4|nabla^2 Phi| + 2 sqrt(((partial u)/ (partial x))^2 + ((partial v)/ (partial y))^2)
$
for my 2D simulation. This size function represent a depth status and move status. Faster moving area should has move details.

=== 3. Solve pressure
For adaptive grid, how to solve pressure is a hard problem for a long time.
It need to build boundary-conforming grid, Vornonoi diagram or Power diagram in the past. (But I don't know what is that.)
This use a modified MLS interpolate and discretize pressure, produce a simple way to solve pressure on adaptive grid.
Following equation is simple discretized pressure Laplacian equation:
$
  nabla^T [V A] [F nabla] p = nabla^T [V A] u
$
where $V_"2D" = 3 Delta x^2$, $0 < A < 1$ is non-solid fraction of a face.
Let $nabla_k$ be the kth row of $nabla$ where k denotes a reference to a face
and $delta x$ be the smaller cell's size.
$
  nabla_k & = "sign"("face") [ -1, 1/2, 1/2] / (1.5 Delta x) "for T junction" \
  nabla_k & = [ -1, 1 ] / (Delta x) "for regular cells"
$
and $F nabla_k = max(W_k, 0) nabla_k$ where
$
  W_k = (Sigma_(j in Q^*_k) c_"kj" Phi_j) / (Sigma_(i in Q_k) c_"ki" Phi_i)
$
where $Q^*_k$ donate all cells adjacent face k and $Q_k$ donate cells under water.
Author gives two way to calculate $F$, and this way I use keep the final Laplacian matrix be positive definition and still don't look difference.
So it can use efficient solver like PCG to solve.
