- Make scalar fields act on Euclidean vectors
- Consider making CircularInterval::hull non-commutative
- Consider not providing concrete default template arguments
- Make our own precondition macro

## Direction

Ordered so that each stage is verifiable on its own and derisks the
next. The limits motivating the work are in
`docs/exactness_frontier.md`.

Where two items look equally valuable, prefer the one that carries over
to the plane, the torus and meshes; `docs/geometry_portability.md`
records what already does. That consideration retired the
spherical-harmonic basis: it addressed conditioning in high-degree
global fits, and the global path is the one that does not generalise.

### 1. Flat C1 elements
The flat densities are C0: quadratics through samples on a periodic
grid, refined to a requested tolerance by probing. The planar
Powell-Sabin element is the textbook case and the projection
machinery is already structured for it, so the sphere's certified,
least-squares C1 path has a direct flat analogue waiting.

*Verified by:* a certified positive floor and a telescoped residual
on a flat density, as on the sphere.

### 2. Walls
Deferred by design: every edge look except the true conformed one
comes from cropping the torus, so the bounded domain and its terms in
the exact Hessian wait until a piece demands what a cut cannot give.
The design is in `docs/geometry_portability.md`.

*Verified by:* a walled rectangle driving the shared optimizer stack.
