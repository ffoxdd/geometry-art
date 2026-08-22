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

**Arithmetic floor — not worth chasing.** The 1e-9 relative capacity
RMS is `sqrt(eps)`, the signature of cancellation in the Stokes
differences. Compensated summation or double-double would move it, and
1e-9 relative area is already sub-micron on a metre globe.

**Smoothness — real, and it is a rate limit rather than a correctness
one.** Degree-`d` Lagrange elements share nodes, so the piecewise
density is C⁰; integrating a C⁰ density over a moving cell gives an
energy that is C¹ but not C². The gradient is continuous everywhere,
so the piecewise stall near 1e-6 is L-BFGS losing curvature
information at the kinks, not a broken derivative. Restoring C² means
making the density C¹ — a global smooth basis, or C¹ elements.

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

1. C¹ elements, so the piecewise path reaches the same floor as the
   polynomial one. This moved to the front once the second-order solve
   landed: on a piecewise field both inner solvers stall on some seeds
   and the second-order one stalls on more, while neither does on the
   polynomial ones. That is the C¹ limit above showing up as a
   measurement, and it caps the solver work until it is lifted.
2. Move the global path to a spherical-harmonic basis.
3. Make accuracy the input and the discretisation the output: request a
   tolerance, raise degree and refine the mesh until the representation
   error is below what the capacity tolerance needs.
4. Image densities on a mesh aligned to their discontinuities.

Done: per-cell integration runs in parallel; the exact CVT Hessian and
the constraint Jacobian are assembled and drive both an unconstrained
trust-region relaxation and a second-order inner solve for the
constrained problem.
