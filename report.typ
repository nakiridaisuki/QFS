#set page(
  paper: "a4",
  numbering: "1",
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
    // numbering("i.", n.last())
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
To simulate the fine details of surface phenomena, three key tasks must be addressed:
1. Representing and tracking the free surface.
2. Subdivide the grid dynamically around the surface.
3. Solving the pressure Poisson equation accurately on the adaptive grid.
Below is an overview of the methodologies proposed to solve these challenges.

=== 1. Tracking the Free Surface
To represent the liquid interface, the simulation utilizes a level-set method,
storing a signed distance field (SDF) value, $Phi$, at each cell center.
This value represents the minimum signed distance from the cell center to the surface,
where a negative value ($Phi < 0$) indicates that the cell is inside (under) the liquid.

While it's easy to initialize $Phi$ precisely for simple geometric shapes (such as circular water drop),
maintaining an accurate SDF over subsequent simulation frames is more challenging.
After semi-Lagrangian advection, the field $Phi$ typically deviates from true SDF and becomes imprecise.
To restore its signed-distance property, a two-step redistancing process is applied:\
Step a: Estimate the location of the zero-level set ($Phi = 0$). \
Step b: Use this estimate as a boundary condition to solve Eikonal equation, $|nabla Phi| = 1$. \

For step a, we can check adjacent cells that have opposite signs for $Phi$, indicating that the surface crosses them.
To estimate new $Phi$, we can use linear interpolation:
$
  Phi_"i_new" = (- Phi_i) / (Phi_i - Phi_j)
$
and then mark the cells as updated.

For step b, starting from updated cells with the smallest $Phi_"i_new"$,
the Fast Marching Method--which operates similarly to Dijkstra's algorithm--is
used to propagate the values outward and solve the Eikonal equation for neighboring cells.
This two-step approach reconstructs a clean, accurate SDF $Phi$ across the domain.

=== 2. Subdivide Grid
The framework use two main criteria to determine whether a grid cell should be subdivide: \
a. $Phi$ < node size\
b. Size function > 1 / node size\
The first criterion ensure that cells intersection the surface are targeted for subdivision.

The second criterion uses a size function, S, to measure local details density based on geometric curvature and flow kinematics:
$
  S = gamma_Phi|nabla dot nabla (Phi + Phi_"solid")| + gamma_u sqrt(Sigma_i ((partial u_i )/ (partial x_i))^2)
$
where $gamma_Phi = 4$ and $gamma_u = 3$ and.
Since there are no solid boundaries in the simulation zone,
the sizing function simplifies for a 2D simulation to:
$
  S = 4|nabla^2 Phi| + 2 sqrt(((partial u)/ (partial x))^2 + ((partial v)/ (partial y))^2)
$
This function represents both geometric complexity (surface curvature) and kinematic activity (velocity gradients).
Consequently, regions with high curvature or rapid flow undergo finer subdivision.

=== 3. Solve pressure
Solving the pressure projection equation on adaptive, non-uniform grids has historically been a significant challenge.
Traditional methods often required constructing boundary-conforming grids, Voronoi diagrams, or Power diagrams,
which are geometrically complex and difficult to implement.

To address this, the authors utilize a modified Moving Least Squares (MLS) interpolation to discretize pressure.
This approach offers a simpler yet accurate formulation to solve for pressure on adaptive grids.
The discretized pressure Laplacian equation is expressed as:
$
  nabla^T [V A] [F nabla] p = nabla^T [V A] u
$
where $V_"2D" = 3 Delta x^2$ represents the control volume for a 2D cell, and $0 < A < 1$ is non-solid fraction of a face.

Let $nabla_k$ denote the discrete gradient operator of face k,
and $delta x$ represent the smaller cell's size. The discrete gradient is defined as:
$
  nabla_k & = "sign"("face") [ -1, 1/2, 1/2] / (1.5 Delta x) &    "for T junction" \
  nabla_k & = [ -1, 1 ] / (Delta x)                          & "for regular cells"
$
and $F nabla_k = max(W_k, 0) nabla_k$ where $W_k$ represents the interpolation weights derived from the level-set values:
$
  W_k = (Sigma_(j in Q^*_k) c_"kj" Phi_j) / (Sigma_(i in Q_k) c_"ki" Phi_i)
$
where $Q^*_k$ denotes the set of all cells adjacent to face k
and $Q_k$ denotes the subset of cells under water.

Although the authors propose two methods to compute the scaling factor $F$,
the chosen formulation ensures that the resulting discrete Laplacian matrix remains symmetric positive-definite (SPD)
while yielding indistinguishable visual results. Maintaining an SPD system is highly advantageous,
as it allows the pressure equation to be solved efficiently using linear solvers such as the Preconditioned Conjugate Gradient (PCG) method.

== What did I do
In this project, I implement the Eulerian/FLIP MAC simulator and an Eulerian Quadtree simulator using C++.
Here's the environment requirement:
- C++ 17+
- CMake 3.14+
- OpenMP 5.0+
and used libraries:
- raylib (for window display)
- raygui (for GUI)
- Eigen  (for PCG solving)
All codes are written by myself and AI since I didn't find the source code.

In Quadtree simulator, I implement pressure solver, surface tension, MLS interpolation and redistancing in this paper.
I didn't implement the smooth method in paper because I didn't understand how to do that.
So I implement simpler but hard to parallelize one.

= Implementation
== MAC grid
On MAC grid, we save velocity on edge, and pressure, density at cell center.
The simulation pipeline has 4 main steps:\
1. Update grid data
2. Add external force
3. Solve pressure to divergence free
4. Post-process

The source data of step 1 can be data from last frame or particles.
And the post-process contain velocity extrapolate and some particles updates.

=== Eulerian Method
In Eulerian method, we store density at cell center and consider a cell be water by a threshold of density.
The simulation pipeline of Eulerian method is:
1. Advect density, velocity from last frame
2. Apply gravity and surface tension
3. Solve pressure and update velocity
4. Velocity extrapolate

==== 1. Advect Datas
To use semi-Lagrangian method advect data, we first calculate position of cell center in last frame.
Then use bilinear interpolation to get value from that position. \
The pseudo code look like:
```c
bilerp: bilinear interpolation function
field_old <- old field from last frame
field <- new field need to be advected

[x, y] <- position of data point

[v_x, v_y] <- velocity at position [x, y]

x_prev = x - v_x * dt
y_prev = y - v_y * dt

field[(x, y)] = bilerp(field_old, x_prev, y_prev);
```

==== 2. Apply Gravity and Surface Tension
Apply gravity is easy, just add gravity force to all horizontal edges (store vertical velocity).
```py
G: gravity constant

for all horizontal edges e:
  if (any cell adjacent e has water):
    e.velo += G * dt
```
We use CSF (Continuum Surface Force) model to calculate surface tension.
Since this is just a additional part of our simulator, I just briefly introduce it.
It contain 4 steps:
1. Initialize color field by density and Jacobi smooth it.
2. Calculate unit normal.
3. Calculate curvature using divergence of normal.
4. Calculate surface tension and apply it.

==== 3. Solve pressure
Recall our equation. If we merge $Delta t$ to $p$, it will become:
$
  nabla^2 p = nabla dot arrow(u)^*
$
Let u, v be horizontal, vertical velocity respectively, the discrete divergence is:
$
  (nabla dot arrow(u)^*)_"i, j" approx (arrow(u)^*_"i+1, j" - arrow(u)^*_"i, j")
  + (arrow(v)^*_"i, j+1" - arrow(v)^*_"i, j")
$
and discrete Laplacian operator is:
$
  (nabla^2 p)_"i, j" approx p_"i-1,j" + p_"i+1,j" + p_"i,j-1" + p_"i,j+1" - 4p_"i,j"
$
Combine them together, we get
$
  p_"i-1,j" + p_"i+1,j" + p_"i,j-1" + p_"i,j+1" - 4p_"i,j"
  = (arrow(u)^*_"i+1, j" - arrow(u)^*_"i, j")
  + (arrow(v)^*_"i, j+1" - arrow(v)^*_"i, j")
$
to solve $p$, it become a huge linear system. To solve it, we can use tools like PCG.
In order to use PCG, the last thing we need to do is switch the sign, let it become positive-definite.
$
  4p_"i,j" - p_"i-1,j" - p_"i+1,j" - p_"i,j-1" - p_"i,j+1"
  = - ((arrow(u)^*_"i+1, j" - arrow(u)^*_"i, j")
    + (arrow(v)^*_"i, j+1" - arrow(v)^*_"i, j"))
$

==== 4. Velocity Extrapolate
Since we advect data using semi-Lagrangian method, some sample points may in the air cell around the surface.
If the air cell's edges don't have velocity, the interpolate result will be affected.
To solve this problem, we can propagate velocity from surface to air for a few grids.
Use BFS can complete this task.

=== FLIP
A major drawback of Eulerian method is its inherent numerical diffusion, which washes away fine structures over time.
This motivates our transition to FLIP, which carries velocity on lagrangian particles to reduce dissipation.

FLIP advect velocities using particles instead of grid interpolation,
successfully preserves fluid energy and fine surface structures.
The simulation pipeline become:
1. Transfer velocity from particles to grid
2. Apply gravity and surface tension
3. Solve pressure and update velocity
4. Velocity extrapolate
5. Transfer velocity from grid to particles
6. Advect particles
7. Maintaining particles

We can notice that the 2, 3, 4 step are the same with Eulerian method.
The only difference is about how to manipulate particles.
==== 1. Transfer velocity between grid and particles
We simply use bilinear interpolation and distribution to transfer velocity between grid and particles.
==== 2. Advect particles
Since we have particle now, advect just like particle simulation. We update particles' position by it's velocity.
```
p.x += p.u * dt
p.y += p.v * dt
```
==== 3. Maintaining particles
Now, a cell has water or not depends on if there are particles in this cell.
To prevent an air cell appear under water, we need to resample particles after particle advection.
We consider a cell is under water if it's around by water cells.

There are two condition for a water cell:\
a. Too many particles\
b. Too few particles

We wish the number of particles in a single cell between 4 to 8.
If there are more then 8 particles in one cell, we simply delete the extra.
If there are less then 4 particles in a under water cell,
then we random sample needed particles number and bilinear interpolate velocity to it.

== Quadtree Grid
