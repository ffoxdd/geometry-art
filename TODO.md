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

### 1. The flat family
Plane, cylinder and torus are one implementation: each direction of
the domain is either wrapped or walled, and the three are the three
settings of that switch. The second instance, and the one that lets a
`Geometry` abstraction be extracted from two cases rather than one.

The torus comes first: both directions wrapped, so no boundary
mathematics at all -- a periodic Delaunay diagram, an identity site
manifold that will expose wherever the spherical one is silently
load-bearing, and `Field::integrals` taking the region by concept
rather than by the spherical type. `geometry/planar` already has the
region types with divergence-theorem moments, so the polynomial layer
is shared. Walls follow, bringing the bounded domain and the wall
terms in the exact Hessian, and unlocking the cylinder and the plane.
Edge treatments -- cut, conformed, fade -- are configuration on top;
the design is in `docs/geometry_portability.md`.

*Verified by:* the optimizer stack driving the torus unchanged.
