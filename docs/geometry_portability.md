# Portability across geometries

The target is capacity-constrained tessellation on any domain: the
sphere, the plane, a flat torus, an arbitrary mesh. The sphere is the
instance that exists. This note records what already carries over, what
does not, and how to use that as a tiebreak when two pieces of work
look equally valuable.

## What is already general

**The sweep rate.** Perturbing a site moves the bisector it shares with
a neighbour at an outward normal speed of

```
(x - s_k) . delta / |s_j - s_k|
```

Deriving this from the planar bisector and from the spherical one gives
the same expression. `CapacityJacobian`, the Lagrangian's site
gradient, and the block structure of the Hessian are therefore already
written in a form that holds on the plane and the flat torus as well.

**The optimizer stack.** `TrustRegionStep`, `CapacityJacobian` and
`CapacityConstrainedHessian` name no region or boundary type at all,
and the two trust-region loops name one spherical vector each. Whatever
produces the per-cell and per-bisector integrals is the boundary; above
it, nothing is spherical except the vector type.

## What is not

Five seams carry all the domain knowledge:

1. Region and boundary types, and their `moments(degree)`. Great-circle
   arcs and the Stokes recursion on S² here; segments and Green's
   theorem on the plane.
2. Diagram construction — which triangulation, and whether it is
   periodic.
3. The site manifold: its constraint, a retraction onto it, and the
   tangent projection. `Normalization` is this for the sphere; on the
   plane it is the identity, which makes the plane a good test of
   whether the abstraction leaks.
4. Field integrals over those regions.
5. `CvtHessian`'s block, which uses the tangentially reduced
   `integral of rho x x^T`. The general form is
   `integral of rho x (x - s_k)^T`; the two differ by a rank-one term
   that the sphere's tangent projection removes.

## The flat family is one implementation

A cylinder and a flat torus are intrinsically flat: both are the plane
under a translation quotient, and their geodesics unroll to straight
lines. So a Voronoi diagram on a cylinder is a planar one periodic in
one direction, and on a torus a planar one periodic in both. Cones and
the other developables join them.

That collapses the target list. There are three classes, not five:

- **Flat** -- plane, cylinder, torus, developables. One implementation
  with a periodicity parameter. Straight edges, moments by the
  divergence theorem, an unconstrained site manifold.
- **Constant curvature** -- the sphere. The implementation that exists.
- **Variable curvature** -- an arbitrary mesh, where curvature
  concentrates at vertices and geodesic bisectors bend as they cross
  them. A different problem rather than a port.

Two consequences worth carrying:

- The flat family needs a **bounded domain**, which the sphere never
  did because it is compact. The plane needs bounds in both directions,
  a cylinder in one, a torus in neither -- so the torus is structurally
  the closest of the three to the sphere.
- On a cylinder embedded in space, chord distance and geodesic distance
  disagree, unlike the sphere where both order the same way. The
  intrinsic one is both the tessellation people want and the one that
  unrolls, so the flat family means intrinsic throughout.

## Using it as a tiebreak

Prefer work that lands in the general part, or in a seam whose other
implementations are well understood.

- A **global spectral basis** is the least transferable choice
  available, and this retired the spherical-harmonic plan. Spherical
  harmonics have an analytic polygon integral; the plane and torus
  would want Fourier instead, and a mesh would want Laplace-Beltrami
  eigenfunctions, which have no analytic polygon integral at all. The
  piecewise representation is what every geometry shares, so the
  conditioning limit on high-degree global fits is a limit on a path
  not being taken.
- **C¹ macro-elements** are worth deriving where the classical
  constructions apply. Powell-Sabin and Clough-Tocher are textbook on
  the plane, and per-face polynomials on a mesh are the planar case, so
  the planar element serves three geometries and the spherical one
  serves one.
- **A mesh is a different difficulty class**, not a port. Geodesic
  bisectors on a polyhedral surface are not straight in the unfolded
  picture, so the closure invariant in `exactness_frontier.md` may not
  hold; that geometry likely needs its own answer.

## The open sequencing question

Abstracting a `Geometry` concept from the sphere alone risks putting
the seams in the wrong places. The plane is the cheapest second
instance — easier moments, a mature triangulation, an identity site
manifold — and adding it before extracting anything would let the
abstraction come from two instances rather than one. Against that, it
delays the piecewise and image work on the sphere, and C¹ elements are
what that work is waiting on.

Both orders are defensible. The planar region types now exist, which
settles nothing about the order but makes the second instance cheaper
to reach.
