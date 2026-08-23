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

### 2. C¹ elements
Spherical Powell-Sabin elements in Bernstein-Bézier form: C¹ quadratics
on a six-way split of each mesh triangle, determined by value and
gradient at the vertices, sampled from any callable. Positivity is a
sign check on the coefficients; mesh resolution comes from the site
count, not the field. Needed so the constraint curvature above is
continuous, which is what lets the second-order rate hold on piecewise
fields; today Newton needs about twice the inner iterations on the
piecewise noise field that it needs on a polynomial one.

*Verified by:* C¹ across every seam, exact reproduction of quadratics,
positivity on a floored input, and noise-field inner iteration counts
matching the polynomial fields'.

### 3. Accuracy as the input
Take a requested accuracy and choose the representation to meet it,
rather than taking a mesh level and a degree and reporting what
accuracy came out. Sits on the two stages above, which are what make raising
degree and refining reliable.

*Verified by:* a requested tolerance being met without the caller
naming a degree or a subdivision level.

### 4. Image densities
Load an image, map it onto the sphere and sample it onto the C¹ field
at the cell scale. A C¹ field cannot hold a discontinuity, and the
bandwidth rule says the tessellation cannot express one either, so
alignment to image edges is unnecessary: a step becomes a ramp one cell
wide.

*Verified by:* a synthetic step edge producing a clean size transition
with no ringing in the fitted density.

### 5. The flat family
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
