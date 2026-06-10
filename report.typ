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
  ==> & nabla^2 p = (nabla dot arrow(u)^*) / (Delta t)
$

To solve the pressure on a discrete grid like a MAC grid, we can turn this Poisson equation into a system of linear equations. More implementation details will be provided in the following sections.

However, the traditional MAC method requires discretizing the entire simulation domain---including deep water regions---using a uniform grid.
In 3D scenarios, this leads to a memory footprint that scales cubically as $O(N^3)$ with the grid resolution N,
posing a severe bottleneck as the simulation scale expands or the resolution increases.

In gaming or visual effects, we don't really need to simulate deep water regions in high detail;
instead, we only focus on surface phenomena like splashes or waves.
In this situation, using an adaptive grid near the surface can help us resolve this bottleneck.
