# The exactness frontier

What the symbolic path can integrate exactly, what bounds the work, and
which limits are worth spending effort on.

## The closure invariant

Exact integration survives exactly as long as the density is
**piecewise polynomial over a partition whose boundaries are the same
kind of curve the cells are made of.**

Great-circle-bounded pieces stay great-circle-bounded under clipping,
and polynomials stay polynomial under restriction, so the Stokes
recursion always terminates. The icosphere mesh is the simplest
instance of this class, not the idea itself: the partition need not be
uniform, need not be an icosphere, and need not be fine.

Consequences:

- A density discontinuity is representable exactly when it lies *along*
  partition boundaries. Approximating an edge by a great-circle
  polyline and making it a mesh edge restores exactness; refining a
  uniform mesh toward the edge does not, and pays Gibbs ringing plus a
  triangle explosion for the failure.
- Small-circle boundaries (radial densities, cap-supported features)
  extend the same machinery — the parametrisation is still
  trigonometric — but the outward conormal is no longer constant along
  an edge, so the arc term carries more algebra.
- A density whose natural description is neither is outside the exact
  path, and must be projected into it before anything else happens.

## Bandwidth matching bounds the work

The constraint is an integral over a cell, so density structure finer
than a cell cannot influence the solution except through its local
average. A cell for `N` sites has angular radius about `2/sqrt(N)`, so
there is a matched bandwidth: represent the density to cell resolution
and no further.

Low-passing to that scale is the correct representation rather than a
concession. It also caps the whole problem — 1000 sites needs harmonic
degree around 50, a few thousand coefficients.

Cost scales as `sites × bandwidth²`. Below cell scale the question is
not well posed, so no method recovers that detail and none should.

## Which limits are real

**Arithmetic floor — not worth chasing.** The 1e-8 relative capacity
RMS is `sqrt(eps)`, the signature of cancellation in the Stokes
differences. Compensated summation or double-double would move it, and
1e-8 relative area is already sub-micron on a metre globe.

**Geometric tolerances — real, and the first thing to rule out.** Any
"close enough counts as equal" rule in the geometry is a discontinuity
in the objective of that size, and a trust region aiming at 1e-8
cannot tell a 1e-6 jump from a lying model. Tolerances stay at rounding
scale, and geometry that is known exactly — the great circle an arc
lies on, the bisector normal from two sites, the clip circle that
closes a clipped polygon — is carried rather than rederived from
rounded points, so near-degenerate pieces contribute in proportion to
their length and nothing else. The diagnostic that separates a
representation fault from a smoothness one is to store a smooth
polynomial piecewise and compare every cell and edge integral against
the global field; `test/piecewise_precision_test.cpp` holds that to
1e-11.

**Smoothness — real, and it is a rate limit rather than a correctness
one.** Degree-`d` Lagrange elements share nodes, so the piecewise
density is C⁰. The CVT energy's slope and curvature both see the
density only through boundary integrals, so they stay continuous; the
constraints' curvature sees the density's gradient on the bisectors,
so it jumps wherever a bisector crosses a mesh edge, and the augmented
Lagrangian is C¹ but not C². With the representation exact, both inner
solvers converge on the piecewise noise field across the seeds tried;
the cost shows up as iteration count, with Newton needing roughly
twice the inner iterations it needs on a polynomial field. Restoring
C² means making the density C¹ — C¹ elements.

**Positivity — real, and now certified rather than sampled for.** The
Bernstein-Bezier form bounds a piece below by its coefficients, so
`PowellSabinInterpolant` reports a lower bound on the density it built
and no search over the sphere can miss a dip. Where the mesh is too
coarse to carry the sampled gradients the interpolant would overshoot
below zero; damping the gradients per vertex pulls it back and keeps the
field C1, and the damping is reported because it measures whether the
mesh resolves the field.

**Positivity of a global fit — real, and it is what rules out the global fit.** Least
squares constrains no value, so fitting a function with a floor
overshoots below it: the degree-8 fit of the floored noise reaches
-0.047 against a floor of 0.2, with an RMS residual of 0.12 on a field
whose range is 0.8. A density that reaches zero cannot be balanced
against the others and one that goes negative is not a density, so the
tessellation distorts around those regions. The fitter reports its
lowest value and the factory warns when it is not positive.
`test/field_positivity_test.cpp` pins both halves: the piecewise field
holds its floor, the global fit does not. A Bernstein-Bezier
representation makes positivity a sign check on the coefficients, which
is the durable fix and another reason the piecewise path is the one
that generalises.

**Conditioning — real, and it is what caps degree.** The homogeneous
`(d, d−1)` basis has exactly `(d+1)²` functions, which is the true
dimension of polynomials restricted to the sphere, so nothing is
wasted. But those monomials are not orthogonal on the sphere, so the
least-squares normal equations degrade as degree rises. An orthonormal
spherical-harmonic basis turns fitting into projection, removes the
solve, and makes truncation by coefficient energy meaningful.

**Non-convexity — real, and inherent to the aesthetic constraint.**
Unweighted cells is what makes this hard; power diagrams would make the
problem semi-discrete optimal transport, and convex. The dimension
count is friendly: `2N` site freedoms against `N−1` independent
constraints, since the capacities sum to the total mass automatically,
leaving an `(N+1)`-dimensional feasible manifold. The question is which
local minimum, not whether a solution exists.

**Combinatorial cost — real, and the least explored.** Per-cell
integration is embarrassingly parallel and currently serial. The
augmented Lagrangian's inner problems stiffen as the penalty grows,
which is inherent to penalty methods. Both second-order ingredients are
already within reach: the constraint Jacobian is exact and sparse, and
the CVT Hessian blocks are built from the same shared-bisector arc
integrals of `rho` and `rho * x` that the capacity gradient already
computes, nonzero only between Delaunay neighbours. A Newton or SQP
step on the constraint manifold therefore costs little beyond a sparse
KKT solve, and attacks the iteration count directly.

## Direction

Ordered by leverage per unit of effort. See `TODO.md` for status.

1. Exact constraint curvature. The inner Newton model is Gauss–Newton:
   it omits the multiplier-weighted curvature of the constraints, which
   does not vanish at a solution because mass has a price wherever the
   density varies. Reading the true curvature by finite differences cuts
   the inner iterations on the noise field by two to three times, so
   the term is worth deriving; it needs the Voronoi vertex velocities
   and the density gradient on the bisectors.
2. C¹ elements, so that curvature is continuous and the second-order
   rate holds on piecewise fields. Spherical Powell–Sabin in
   Bernstein–Bézier form: value and gradient sampled at mesh vertices
   from any callable, positivity by coefficient sign, mesh resolution
   from the site count.
3. Make accuracy the input and the discretisation the output: request a
   tolerance, raise degree and refine the mesh until the representation
   error is below what the capacity tolerance needs.
4. Image densities, sampled onto the C¹ field at the cell scale; a C¹
   field cannot hold a discontinuity, so alignment to image edges is
   unnecessary under the bandwidth rule.

Done: per-cell integration runs in parallel; the exact CVT Hessian and
the constraint Jacobian are assembled and drive both an unconstrained
trust-region relaxation and a second-order inner solve for the
constrained problem; the piecewise representation agrees with the
global one to 1e-11 on every cell and edge, and every seed that used to
stall on it converges.
