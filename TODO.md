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
records what already does. That consideration is why the
spherical-harmonic item sits last rather than in the middle.

### 1. C¹ elements
Spherical Clough–Tocher or Powell–Sabin elements, or a smooth
partition-of-unity blend. Now the first thing to do rather than a
refinement: on a piecewise field both inner solvers stall on some
seeds and the second-order one stalls on more, because that energy is
C¹ but not C². Until it is C², the faster solver cannot be trusted on
exactly the fields that matter for images.

*Verified by:* a piecewise field reaching the polynomial path's
precision floor across seeds that stall today.

### 2. Accuracy as the input
Take a requested accuracy and choose the representation to meet it,
rather than taking a mesh level and a degree and reporting what
accuracy came out. Sits on the two stages above, which are what make raising
degree and refining reliable.

*Verified by:* a requested tolerance being met without the caller
naming a degree or a subdivision level.

### 3. Image densities
Load an image, map it onto the sphere and sample it onto a piecewise
field. Choose mesh resolution from the cell scale rather than the image
resolution, and align mesh edges to the image's discontinuities: an
edge approximated by a great-circle polyline stays exactly integrable,
while refining a uniform mesh toward it does not.

*Verified by:* cells aligning along a synthetic step edge, with no
ringing in the fitted density.

### 4. The flat family
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

### 5. Spherical-harmonic basis
Projection onto an orthonormal basis in place of the least-squares fit
onto the homogeneous `(d, d-1)` monomial basis, which has the right
dimension but is not orthogonal on the sphere, so conditioning is what
caps degree today. Needs analytic `∫_P Y_lm dA` over a spherical
polygon.

Last rather than middle: it is the least transferable item on the list.
The plane and torus would want Fourier, a mesh would want
Laplace-Beltrami eigenfunctions with no analytic polygon integral, and
its motivation is the global path, which is the one that does not
generalise.

*Verified by:* fitting at degrees the monomial path cannot condition,
with residual falling monotonically in degree.
