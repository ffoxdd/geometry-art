- Make scalar fields act on Euclidean vectors
- Consider making CircularInterval::hull non-commutative
- Consider not providing concrete default template arguments
- Make our own precondition macro

## Direction

Ordered so that each stage is verifiable on its own and derisks the
next. The limits motivating the work are in
`docs/exactness_frontier.md`.

Speed comes before representation: the system already handles the
fields we care about, and handles them slowly. Within speed, the
second-order work is split so the matrix is proved correct before any
optimizer depends on it.

### 1. Benchmark harness
Record iteration counts and wall time for a fixed set of
field-and-site-count configurations, so later stages are measured
rather than asserted.

*Verified by:* rerunning it reproduces the same iteration counts.

### 2. Parallel per-cell integration
Cell integrals within one gradient evaluation are independent. The
Delaunay rebuild stays serial. A constant factor, but it multiplies
every stage after it.

*Verified by:* identical results to the serial path, bitwise where the
summation order is fixed.

### 3. Exact Hessian assembly
`SphereState` already stores the per-edge `RegionIntegrals` the Hessian
needs: the off-diagonal block for a Delaunay neighbour is the shared
bisector's arc integral weighted by the normal velocity that the
capacity gradient already uses. Assemble the sparse matrix and the
constraint Jacobian; change no optimizer yet.

*Verified by:* agreement with finite differences of the exact gradient,
to the same tolerance the gradient tests use.

### 4. Trust-region Newton on the unconstrained energy
Replace the Lloyd warm start. Self-contained, and it exercises the
Hessian and the Riemannian correction without the constraints in play.
Needs a trust region rather than a raw Newton step, since the energy is
not convex.

*Verified by:* reaching a lower energy than Lloyd in fewer evaluations
on the polynomial fields.

### 5. SQP on the constraint manifold
Sparse KKT solves in place of a growing penalty, whose inner problems
stiffen by construction.

*Verified by:* the capacity floor reached in an order of magnitude
fewer inner iterations.

### 6. C¹ elements
Spherical Clough–Tocher or Powell–Sabin elements, or a smooth
partition-of-unity blend. Stages 4 and 5 only pay their full return
where the energy is C², so until this lands the piecewise and image
paths gain a constant factor rather than a rate.

*Verified by:* a piecewise field reaching the polynomial path's
precision floor.

### 7. Spherical-harmonic basis
Replace the least-squares fit onto the homogeneous `(d, d−1)` monomial
basis with projection onto an orthonormal basis. The monomial basis has
the right dimension but is not orthogonal on the sphere, so
conditioning is what caps degree today. Needs analytic `∫_P Y_lm dA`
over a spherical polygon, and the `l → l ± 1` recurrence for
multiplication by a coordinate to keep the first moments closed.

*Verified by:* fitting at degrees the monomial path cannot condition,
with residual falling monotonically in degree.

### 8. Accuracy as the input
Take a requested accuracy and choose the representation to meet it,
rather than taking a mesh level and a degree and reporting what
accuracy came out. Sits on stages 6 and 7, which are what make raising
degree and refining reliable.

*Verified by:* a requested tolerance being met without the caller
naming a degree or a subdivision level.

### 9. Image densities
Load an image, map it onto the sphere and sample it onto a piecewise
field. Choose mesh resolution from the cell scale rather than the image
resolution, and align mesh edges to the image's discontinuities: an
edge approximated by a great-circle polyline stays exactly integrable,
while refining a uniform mesh toward it does not.

*Verified by:* cells aligning along a synthetic step edge, with no
ringing in the fitted density.
