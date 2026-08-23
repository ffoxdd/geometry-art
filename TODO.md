- Make scalar fields act on Euclidean vectors
- Consider making CircularInterval::hull non-commutative
- Consider not providing concrete default template arguments
- Make our own precondition macro

## Direction

Ordered so that each stage is verifiable on its own and derisks the
next. The limits motivating the work are in
`docs/exactness_frontier.md`.

Speed comes before representation: the system already handles the
fields we care about, and handles them slowly.

Where two items look equally valuable, prefer the one that carries over
to the plane, the torus and meshes; `docs/geometry_portability.md`
records what already does. That consideration retired the
spherical-harmonic basis: it addressed conditioning in high-degree
global fits, and the global path is the one that does not generalise.

### 1. Exact constraint curvature
The inner Newton model is Gauss-Newton: exact for the energy, slope-only
for the penalty, and missing the multiplier-weighted curvature of the
constraints altogether. That term does not vanish at a solution, because
mass has a price wherever the density varies. Reading the true curvature
by finite differences (`--newton-curvature finite-difference`) cuts the
inner iterations on the noise field by two to three times. Deriving it
needs the Voronoi vertex velocities -- the implicit function theorem on
the three equidistance equations -- and the density gradient along the
bisectors. Both are geometry-portable in the way the sweep rate is.

*Verified by:* the analytic curvature matching the finite-difference
one to the difference error, and matching its iteration counts.

### 2. Project onto the space rather than fitting each vertex separately
Vertex readings now come from a least-squares quadratic fit over a
neighbourhood the size of the mesh rather than from samples at the vertex
(`LocalQuadraticFit`). Structure finer than the mesh lands in the residual
instead of being aliased into the result, and no gradient is ever read as
a difference quotient. Measured on the noise at 200 sites, worst relative
mass error over cell-sized caps, and the certified lower bound:

| mesh | cell-scale error | lower bound |
|---|---|---|
| level 3, about the cell scale | 5.8% -> 3.9% | 0.039 -> 0.145 |
| level 4, a quarter of it | 0.97% -> 0.46% | 0.066 -> 0.198 |

The lower bound matters as much as the error: the field now keeps the 0.2
floor it was given, where before it dipped to a third of that and so posed
a harder capacity problem than the one asked for. Gradient damping no
longer triggers at any resolution.

What remains is that each vertex is still fitted on its own. The best
approximation in the space is its L2 projection, which is one sparse
symmetric solve over the whole mesh rather than a fit per vertex, and it
would close the rest of the gap to the Lagrange field. Positivity would
then need the damping backstop that currently never fires.

*Verified by:* cell-scale mass error at a mesh matched to the cell scale
falling to the level a mesh four times finer reaches today.

### 3. Choose the resolution from the site count
Only worth doing once the operator above stops aliasing, because until
then the answer is "several times finer than the cell", which is the
answer that makes the choice not worth computing.

*Verified by:* a requested cell-scale tolerance being met without the
caller naming a subdivision level.

### 4. Accuracy as the input
Take a requested accuracy and choose the representation to meet it,
rather than taking a mesh level and a degree and reporting what
accuracy came out. Sits on the two stages above, which are what make raising
degree and refining reliable.

*Verified by:* a requested tolerance being met without the caller
naming a degree or a subdivision level.

### 5. Image densities
Load an image, map it onto the sphere and sample it onto the C¹ field
at the cell scale. A C¹ field cannot hold a discontinuity, and the
bandwidth rule says the tessellation cannot express one either, so
alignment to image edges is unnecessary: a step becomes a ramp one cell
wide.

*Verified by:* a synthetic step edge producing a clean size transition
with no ringing in the fitted density.

### 6. The flat family
Plane, cylinder and torus are one implementation: all three are
intrinsically flat, differing only in which directions wrap. The
cheapest second instance, and the one that would let a `Geometry`
abstraction be extracted from two cases rather than one.

`geometry/planar` now has the region types, with moments by the
divergence theorem reported in the same table the spherical ones use,
so the polynomial layer is already shared. What remains: a bounded
domain, which the sphere never needed; a diagram from
`Delaunay_triangulation_2`; an identity site manifold, which will
expose wherever the spherical one is silently load-bearing; and
`Field::integrals` taking the region by concept rather than by the
spherical type. Sequencing against C¹ elements is open -- see
`docs/geometry_portability.md`.

*Verified by:* the optimizer stack driving it unchanged.
