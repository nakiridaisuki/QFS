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
#set heading(numbering: (..nums) => {
  let n = nums.pos()
  if n.len() <= 3 {
    numbering("I.1.a", ..n)
  }
})
#show raw.where(block: true): it => block(
  fill: luma(240),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
  it,
)
#set align(left)

#set page(
  paper: "a4",
  margin: (x: 2.5cm, top: 3cm, bottom: 2.5cm),
  header: align(right)[
    #text(8pt, fill: luma(120))[CA2026 Final Project: Adaptive Grid Fluid Simulator]
  ],
  numbering: "1",
)

// 美化標題與作者區塊
#align(center)[
  #v(2cm)
  #text(24pt, weight: "bold")[CA2026 Final Project] \
  #v(0.5em)
  #text(18pt, fill: rgb("1c3d5a"))[Adaptive Grid Fluid Simulator]
  #v(2cm)

  #grid(
    columns: 1fr,
    gutter: 1.2em,
    text(12pt)[*Author:* 蔡昀呈 (113550058)],
    text(11pt, fill: luma(100))[Department of Computer Science \ National Yang Ming Chiao Tung University]
  )
  #v(2cm)
]


#outline(indent: 1.5em, depth: 2)
#pagebreak()

= Introduction

This project implements the method described in
Ryoichi Ando and Christopher Batty. 2020. "A Practical Octree Liquid Simulator with Adaptive Surface Resolution." ACM Transactions on Graphics (TOG).

== Motivation
A few months ago, I watched some captivating fluid simulator videos on Youtube (#link("https://youtu.be/rSKMYc1CQHE?si=LpY_0M06A26NAMnJ")[e.g., this channel])
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

#pagebreak()
However, the traditional MAC method requires discretizing the entire simulation domain---including deep water regions---using a uniform grid.
In 3D scenarios, this leads to a memory footprint that scales cubically as $O(N^3)$ with the grid resolution N,
posing a severe bottleneck as the simulation scale expands or the resolution increases.

In computer graphics and real-time visual effects, high-resolution details are typically unnecessary in deep water regions;
instead, we only focus on surface phenomena like splashes or waves.
In this situation, using an adaptive grid near the surface can help us resolve this bottleneck.

== Overview of the Reference Method
To simulate the fine details of surface phenomena, three key tasks must be addressed:
1. Representing and tracking the free surface.
2. Subdividing the grid dynamically around the surface.
3. Solving the pressure Poisson equation accurately on the adaptive grid.
Below is an overview of the methodologies proposed to solve these challenges.

=== Tracking the Free Surface <sec-surface>
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
  Phi_"i_new" = (|Phi_i|) / (|Phi_i| + |Phi_j|)
$
and then mark the cells as updated.

For step b, starting from updated cells with the smallest $Phi_"i_new"$,
the Fast Marching Method--which operates similarly to Dijkstra's algorithm--is
used to propagate the values outward and solve the Eikonal equation for neighboring cells.
This two-step approach reconstructs a clean, accurate SDF $Phi$ across the domain.

#pagebreak()
=== Grid Subdivide <sec-subdivide>
The framework uses two main criteria to determine whether a grid cell should be subdivided: \
a. $Phi$ < node size\
b. Size function > 1 / node size\
The first criterion ensures that cells intersecting the surface are targeted for subdivision.

The second criterion uses a size function, S, to measure local details density based on geometric curvature and flow kinematics:
$
  S = gamma_Phi|nabla dot nabla (Phi + Phi_"solid")| + gamma_u sqrt(Sigma_i ((partial u_i )/ (partial x_i))^2)
$
where $gamma_Phi = 4$ and $gamma_u = 3$ and.
Since there are no solid boundaries in the simulation zone,
the sizing function simplifies for a 2D simulation to:
$
  S = 4|nabla^2 Phi| + 3 sqrt(((partial u)/ (partial x))^2 + ((partial v)/ (partial y))^2)
$
This function represents both geometric complexity (surface curvature) and kinematic activity (velocity gradients).
Consequently, regions with high curvature or rapid flow undergo finer subdivision.

=== Solving the Pressure <sec-pressure>
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

== Project Scope
In this project, I implemented the Eulerian/FLIP MAC simulator and an Eulerian Quadtree simulator using C++.
Here's the environment requirement:
- C++ 17+
- CMake 3.14+
- OpenMP 5.0+
and used libraries:
- raylib (for window display)
- raygui (for GUI)
- Eigen  (for PCG solving)

Because the original source code was not publicly available, the entire codebase was developed from scratch, with assistance from AI tools to accelerate the implementation.

In the Quadtree simulator, I successfully implemented the pressure solver, surface tension, Moving Least Squares (MLS) interpolation, and the level-set redistancing algorithm described in the paper. However, the surface smoothing method proposed by the authors was omitted due to its complexity. Instead, I implemented a simpler alternative, though it presents greater challenges for parallelization.

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

==== 1. Data Advection
To advect data using the semi-Lagrangian method, we first calculate position of cell center in last frame.
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
Apply gravity is straightforward, just add gravity force to all horizontal edges (store vertical velocity).
```py
G: gravity constant

for all horizontal edges e:
  if (any cell adjacent e has water):
    e.velo += G * dt
```
We use CSF (Continuum Surface Force) model to calculate surface tension.
Since this is just an additional feature of our simulator, I just briefly introduce it.
It contains 4 steps:
1. Initialize color field by density and Jacobi smooth it.
2. Calculate unit normal.
3. Calculate curvature using divergence of normal.
4. Calculate surface tension and apply it.

#pagebreak()
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
Since we advect data using semi-Lagrangian method, some sample points may be in the air cell around the surface.
If the air cell's edges don't have velocity, the interpolate result will be affected.
To solve this problem, we can propagate velocity from surface to air for a few grids.
A Breadth-First Search (BFS) was implemented to complete this velocity propagation task.

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
Since we have particle now, advect just like particle simulation. We update particles' position by its velocity.
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

The QT(Quadtree) simulator rebuilds the entire quadtree from scratch each frame.
The simulation pipeline is:
1. Build a new quadtree refined around the surface and high-velocity regions.
2. Apply gravity and surface tension on the new edges.
3. Solve the pressure Poisson equation (projection).
4. Extrapolate velocity from fluid cells to nearby air cells (same as MAC, on adaptive edges).
5. Reconstruct $Phi$ via redistancing.

Steps 2 and 4 follow the same principles as in MAC but operate on the adaptive edge and leaf structures.
The key implementation differences lie in step 1 (tree construction), step 3 (pressure discretization),
and step 5 (redistance).

=== Building a New Tree Each Frame

The tree reconstruction consists of six sub-steps, starting from the existing tree
and producing a fully refined new tree with all data advected.

==== 1. Computing the Size Function

We compute $S$ for each leaf using the formula from @sec-subdivide. The Laplacian of $Phi$ uses
non-uniform finite differences to account for varying neighbor cell sizes ($L_"dis"$, $R_"dis"$, etc.
are the distances to each neighbor). The velocity gradients
$(partial u) / (partial x)$ and $(partial v) / (partial y)$ are estimated from the leaf's adjacent u and v edges.
```c
// For each leaf. L/R/U/D_dis and L/R/U/D_phi is neighbor's data
d2f_dx2 = 2 * (L_phi / (L_dis * (L_dis + R_dis))
            - leaf.phi / (L_dis * R_dis)
            + R_phi / (R_dis * (L_dis + R_dis)))
d2f_dy2 = ... // computed similarly with D/U
S_geo = gamma_phi * |d2f_dx2 + d2f_dy2|

u_diff <- difference of horizontal velocity
v_diff <- difference of vertical velocity
S_vel = gamma_u * sqrt(u_diff^2 + v_diff^2)
        / leaf.size

S_new = S_geo + S_vel
S_advected = MLS interpolated S from old tree

S = max(R_t * S_advected, S_new)  // R_t = 0.9^(dt / 0.01)
```

==== 2. Propagating the Size Function

To avoid abrupt refinement changes, we propagate $S$ over 5 iterations via area-weighted
averaging with neighbors:
```c
for iter in 1..5:
  for each leaf:
    total = leaf.S * leaf.area
    total_area = leaf.area
    for each cached neighbor:
        total += max(neighbor.S, leaf.S) * neighbor.area
        total_area += neighbor.area
    leaf.S = total / total_area
```

==== 3. Recursive Tree Subdivision

We allocate a new root and recursively subdivide using the two criteria from @sec-subdivide
($|Phi| < "size"$ and $S > 1 / "size"$).
My code can only handle square simulation zone with power of 2 side length, so the minimum cell size will be 1.
$Phi$ and $S$ at child nodes are obtained via MLS interpolation from the old tree.
User water add/delete interactions merge a circle SDF into $Phi$ during this traversal.
```c
recursiveBuildTree(node):
    if node.size < 1.5: return

    exp_phi = advect(node.x, node.y)    // semi-Lagrangian from old tree
    exp_S   = MLSinterpolate(node.x, node.y)
    if user_interaction: exp_phi = merge(circleSDF, exp_phi)

    if |exp_phi| < node.size and exp_S > 1/node.size:
        subdivide(node)  // allocate 4 children, inherit parent edge ids
        for child in children: recursiveBuildTree(child)
```

==== 4. Smoothing the Tree

We utilize a min-heap based on node sizes to enforce the 2:1 depth ratio constraint,
ensuring that adjacent leaves satisfy $|"depth"_a - "depth"_b| <= 1$.
When a leaf detects a neighbor more than one level coarser,
that neighbor is subdivided. This bounds T-junctions to at most two small cells meeting
one large cell, keeping the pressure stencil tractable.

==== 5. Finding All Edges

This step translates the cell-centered tree into staggered (MAC-style) velocity edges.
For each leaf, we iterate over its four faces (left, up, right, down) and decide whether
to create a new edge or reuse an existing one from a neighbor.

First, we adopt a *face ownership rule*: when two equal-depth cells share a face,
only the right/bottom cell creates the edge; the left/top cell will reuse that edge.
This avoids creating duplicate edges at regular interfaces.

For each leaf face from large to small, we examine its neighbor:
- *No neighbor*: The face is on the domain boundary. Create an edge with
  `solid_fraction = 1`, the leaf itself as the sole adjacent cell, and gradient
  coefficient ${0}$ (the face carries zero velocity regardless of pressure).
- *Two smaller neighbors (T-junction)*: The face spans this leaf on one side and
  two smaller leaves on the other. Create an edge with three adjacent cells---this
  leaf plus the two smaller neighbors---using the T-junction gradient stencil
  $"sign" dot [-1, 1/2, 1/2] / (1.5 Delta x)$ from @sec-pressure The smaller
  leaves' coefficients are both $0.5$, and the larger leaf's coefficient is $-1$.
- *One neighbor*:
  - If the neighbor has the same depth and we are on the right or bottom face
    (`dir > 1`): create a regular edge with two adjacent cells and the
    $[-1, 1] / Delta x$ stencil from @sec-pressure (negated for left/top).
  - If the neighbor has the same depth and we are on the left or top face
    (`dir <= 1`): this is an overlay case---the neighbor is the right/bottom
    cell that already owns this face. Reuse its existing edge id.
  - If the neighbor is larger (coarser): overlay case---our face maps to the
    larger neighbor's corresponding edge. Reuse that edge id.

In the code, for parallelize, we first count how many edges each leaf will create,
then allocate the u and v edge arrays and fill them.
The overlay case is done in a separate pass after all edges are created.

==== 6. Advecting Data to the New Tree

Cell-centered data ($Phi$, $S$) and edge velocities are advected via semi-Lagrangian
backward trace + MLS interpolation from the old tree. Newly created edges skip this step
and retain their initial zero velocity.

#pagebreak()
=== Pressure Solve on the Adaptive Grid

Recall the discrete system from @sec-pressure $nabla^T [V A] [F nabla] p = nabla^T [V A] u^*$,
where $V = 3 Delta x^2$, $A = 1 - "solid_fraction"$, and $F = max(W_k, 0.01)$.

We first compute a `fluid_id` mapping (via prefix sum) to index only the $Phi < 0$ cells
into the linear system. Then we iterate over all u and v edges:
```c
for each face in QTu, QTv:
  W_k = computeWk(face)  // use the formula
  W_prime = max(W_k, 0.01)
  face_weight = 3 * face.length^2 * (1 - solid_frac) * W_prime  // VA * F

  for each adjacent fluid cell_i:
    div[fluid_id[i]] += VA * face.val * grad_coeff[i] // RHS (VA * u^*)
    for each adjacent fluid cell_j:
      a = face_weight * grad_coeff[i] * grad_coeff[j]
      Add a into correspond position in Laplacian matrix
```

Each OpenMP thread accumulates triplets locally; after the face loop, thread-local vectors
are merged. The resulting SPD matrix is solved via Eigen's PCG with tolerance $10^(-3)$.
After solving, velocity is updated:
```c
for each face:
  grad_p = sum of (grad_coeff * pressure) for each adjacent fluid cell
  face.val -= W_prime * grad_p
```

=== Level Set Maintenance

The redistancing follows the two-stage approach from @sec-surface, but adapted to the
non-uniform quadtree.

==== Stage 1: Surface Cell Detection

As described in @sec-surface step a, we check cached neighbors for opposite-sign $Phi$.
The distance to interface is estimated by linear interpolation (the formula from @sec-surface),
scaled by the average of the two cell sizes $(L_i + L_j) / 2$ to account for non-uniform resolution:
$
  Phi_"i_new" = (|Phi_i|) / (|Phi_i| + |Phi_j|) dot (L_i + L_j) / 2
$
Known surface cells are collected into a contiguous array,
with each known cell's signed $Phi$ serving as the initial boundary condition for FMM.

==== Stage 2: Fast Marching Method

This implements step b from @sec-surface A min-heap initialized with known surface cells
propagates $Phi$ outward. The per-neighbor Eikonal solve differs from the uniform case:
each axis may have a different effective grid spacing $h$ (average of the two cell sizes).
```c
// For each unvisited neighbor, collect the best known |phi| per axis:
for each known neighbor n_i of target cell:
  axis = (|dx| > |dy|) ? X_AXIS : Y_AXIS
  if |n_i.phi| < axis_data[axis].u:
    axis_data[axis] = {|n_i.phi|, 0.5*(target.size + n_i.size)}

// Sort valid axes by u ascending, then solve:
u_1D = axes[0].u + axes[0].h
if len(axes) >= 2 and u_1D > axes[1].u:
  // 2D quadratic: A*u^2 + B*u + C = 0
  // from (u - u1)^2 / h1^2 + (u - u2)^2 / h2^2 = 1
  u_new = (-B + sqrt(disc)) / (2*A)
else:
  u_new = u_1D
```
After FMM completes, `phi_new` replaces `phi` for all leaves.


== Optimizations

In both the uniform MAC and Quadtree simulators, several optimization strategies
were employed to achieve real-time performance at moderate grid resolutions.

=== Reducing Memory Allocation Overhead

In both MAC and Quadtree simulator, all vectors utilized during simulation will be
declared as a member variable to prevent memory reallocation.

In Quadtree simulation. The tree nodes are stored in `node_pool`---a double buffer where `root_list` toggles between 0 and 1 each frame.
At the start of building new tree, the `node_pool[1 - root_list]` is cleared and reused, avoiding repeated deallocation.
And in MLS interpolation, edge mirroring uses an epoch-based visited-faces
scheme ($O(1)$ clear via incrementing a `uint64_t` counter) to deduplicate samples without resetting arrays.

=== OpenMP Parallelization

Both simulators use `#pragma omp parallel for` extensively for all simple for loop parallelization.

=== Caching Leaves and Neighbors

Leaves are collected into a flat array.
A leaf table maps any integer grid coordinate to its containing leaf in $O(1)$.
For each leaf, eight probe points (two per direction) query the table to build cached neighbor and
per-direction neighbor count. This precomputation amortizes neighbor lookups used throughout
the pipeline.

#pagebreak()
=== Parallel Prefix Sum Optimization

The Quadtree simulator frequently needs to build compact arrays from sparsely marked
elements. A recurring pattern is:

1. Mark valid elements with 1 in parallel.
2. Compute a sequential prefix sum to determine output offsets.
3. Scatter data in parallel using the precomputed offsets.

This pattern is used in `collectLeafNodes()` (gathering all leaf indices), `findAllEdges()`
(computing u/v edge counts and write offsets), `project()` (building the fluid cell ID
mapping for the pressure system), and `redistancing()` (collecting known surface cells
for FMM initialization). While the sequential prefix sum is $O(N)$ and not parallelized,
the surrounding marking and scattering steps are fully parallel, yielding a significant
speedup over atomic-push alternatives.

= Result

== MAC Simulator
Based on testing, the PCG solver's tolerance error needs to be less than 0.1 to produce visually correct behavior.
The main bottleneck in the MAC simulator is the PCG solver. It often requires hundreds of iterations, depending on the size of the simulation area.
The simulator runs at approximately 30 FPS at a 512x512 resolution.

== Quadtree Simulator
Unlike the MAC simulator, the PCG solver is not the main bottleneck. It only needs about 20 iterations to reach an error of $10^{-3}$ (which is better than the MAC simulator).
However, it only runs at about 15 FPS at a 512x512 resolution.
After profiling, I found that the main costs come from parallelization and MLS interpolation.
Since this simulator has more serial steps, we need to use different parallel algorithms.
Also, MLS interpolation is a heavy calculation. Since the system automatically falls back to bilinear/trilinear interpolation anyway, we should implement a faster calculation for these fallback cases instead of using MLS.

#align(center)[
  #table(
    columns: (auto, auto, auto, auto, auto),
    inset: 10pt,
    align: center + horizon,
    [*Simulator*], [*Resolution*], [*PCG Iterations*], [*Frame Rate*], [*Main Bottleneck*],
    [Uniform MAC], [512x512], [100+], [~30 FPS], [PCG Linear Solver],
    [Adaptive QT], [512x512], [~20], [~15 FPS], [Tree Rebuild & MLS],
  )
]

= Limitations & Future Work

While the Eulerian solvers perform robustly, the Lagrangian particle extension (EXNBFLIP) remains incomplete. The primary challenge lies in reconstructing a smooth, continuous signed distance field (SDF) $Phi$ from discrete particle distributions.

Simply calculating the minimum distance to the nearest particle in water cells does not yield a correct level set, as particles often cluster unevenly. Although assigning a predefined radius to each particle is a common alternative, finding an optimal kernel radius that prevents both artificial surface volume loss and numerical noise proved difficult during development.

Furthermore, particle resampling introduces stability challenges. Since the surface is reconstructed dynamically from the particle layout, any noise that accidentally places a resampled particle in the air zone ($Phi > 0$) causes the surface to continuously expand outward. This artifact propagates across frames, eventually filling the entire simulation domain with phantom fluid.

In future work, I plan to address these limitations by:
1. Implementing a more robust narrow-band particle-to-level-set projection method.
2. Exploring alternative surface-smoothing algorithms on the adaptive grid.
3. Optimizing the MLS interpolation step to reduce the serial processing overhead in the Quadtree simulator.

= Conclusion

In this project, I implemented and evaluated both a uniform MAC-grid fluid simulator and an adaptive Quadtree-grid fluid simulator based on the state-of-the-art formulation by Ando and Batty (2020).

The uniform MAC simulator achieves interactive framerates (30 FPS at $512 times 512$) but suffers from high iteration counts in the PCG solver due to the large grid size. Conversely, the Quadtree simulator demonstrates exceptional efficiency in pressure projection, requiring only about 20 iterations to reach $10^(-3)$ tolerance. However, its overall performance is bounded at approximately 15 FPS due to the computational overhead of dynamic tree reconstruction, parallelization bottlenecks in serial steps, and heavy MLS interpolations.

This project provided valuable hands-on experience in solving the Navier-Stokes equations, managing complex data structures on adaptive grids, and identifying system bottlenecks through profiling.
