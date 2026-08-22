- Make scalar fields act on Euclidean vectors
- Consider making CircularInterval::hull non-commutative
- Consider not providing concrete default template arguments
- Make our own precondition macro

## Direction

Ordered by leverage. The reasoning behind the ordering, and the limits
that motivate each item, are in `docs/exactness_frontier.md`.

### Accuracy as the input
Take a requested accuracy and choose the representation to meet it,
rather than taking a mesh level and a degree and reporting what
accuracy came out. Raise degree and refine until the representation
error is below what the capacity tolerance needs.

### Spherical-harmonic basis for the global path
Replace the least-squares fit onto the homogeneous `(d, d−1)` monomial
basis with projection onto an orthonormal spherical-harmonic basis.
The monomial basis has the right dimension but is not orthogonal on the
sphere, so conditioning is what currently caps degree. Needs analytic
`∫_P Y_lm dA` over a spherical polygon, and the `l → l ± 1` recurrence
for multiplication by a coordinate to keep the first moments closed.

### Parallel per-cell integration
Cell integrals within one gradient evaluation are independent. The
Delaunay rebuild stays serial.

### Second-order step on the constraint manifold
Assemble the exact CVT Hessian and constraint Jacobian — both sparse,
both built from arc integrals the gradient already computes — and take
Newton or SQP steps instead of growing a penalty. Needs a trust region
or regularisation, since the energy is not convex.

### C¹ elements
Spherical Clough–Tocher or Powell–Sabin elements, or a smooth
partition-of-unity blend, to make the piecewise energy C² and let it
reach the same floor as the polynomial path.

### Image densities
Load an image, map it onto the sphere and sample it onto a
PiecewisePolynomialField. Choose mesh resolution from the cell scale
rather than the image resolution, and align mesh edges to the image's
discontinuities: an edge approximated by a great-circle polyline stays
exactly integrable, while refining a uniform mesh toward it does not.
