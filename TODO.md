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

### 1. SQP on the constraint manifold
Sparse KKT solves in place of a growing penalty, whose inner problems
stiffen by construction. The CVT Hessian is exact and the constraint
Jacobian already is; the constraints' own second derivatives need the
Voronoi vertices' velocities, so start from the Gauss-Newton form,
where the penalty contributes `penalty * J^T J` and needs no new
derivation.

*Verified by:* the capacity floor reached in an order of magnitude
fewer inner iterations.

### 2. C¹ elements
Spherical Clough–Tocher or Powell–Sabin elements, or a smooth
partition-of-unity blend. The second-order work only pays its full return
where the energy is C², so until this lands the piecewise and image
paths gain a constant factor rather than a rate.

*Verified by:* a piecewise field reaching the polynomial path's
precision floor.

### 3. Spherical-harmonic basis
Replace the least-squares fit onto the homogeneous `(d, d−1)` monomial
basis with projection onto an orthonormal basis. The monomial basis has
the right dimension but is not orthogonal on the sphere, so
conditioning is what caps degree today. Needs analytic `∫_P Y_lm dA`
over a spherical polygon, and the `l → l ± 1` recurrence for
multiplication by a coordinate to keep the first moments closed.

*Verified by:* fitting at degrees the monomial path cannot condition,
with residual falling monotonically in degree.

### 4. Accuracy as the input
Take a requested accuracy and choose the representation to meet it,
rather than taking a mesh level and a degree and reporting what
accuracy came out. Sits on the two stages above, which are what make raising
degree and refining reliable.

*Verified by:* a requested tolerance being met without the caller
naming a degree or a subdivision level.

### 5. Image densities
Load an image, map it onto the sphere and sample it onto a piecewise
field. Choose mesh resolution from the cell scale rather than the image
resolution, and align mesh edges to the image's discontinuities: an
edge approximated by a great-circle polyline stays exactly integrable,
while refining a uniform mesh toward it does not.

*Verified by:* cells aligning along a synthetic step edge, with no
ringing in the fitted density.
