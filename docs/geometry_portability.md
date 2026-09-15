# Portability across geometries

The target is capacity-constrained tessellation on any domain: the
sphere, the plane, a flat torus, an arbitrary mesh. The sphere, the
flat torus and the walled plane are the instances that exist, sharing
one optimizer stack.
This note records what carries over, what does not, and how to use
that as a tiebreak when two pieces of work look equally valuable.

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

**The optimizer stack.** The trust-region Newton descent, the
augmented Lagrangian outer loop, `TrustRegionStep`, `CapacityJacobian`
and the Hessian block structure exist once, in `geometry_art::voronoi`, and
drive both geometries through small model classes. The diagram state
is stored in chart-invariant relative moments -- each bisector's
moments about its two sites, with the separation -- so one bisector
serves cells whose charts sit a period apart, and the Jacobian needs
no site coordinates at all.

**The contracts.** `testing/contracts` states the laws any diagram and
any assembled state must satisfy -- partition of the domain, neighbor
symmetry, equidistance, the separation identity on relative moments --
and both geometries instantiate the same typed suites.

**The skeleton.** Thickening the edge graph into a printable solid is
one construction, in `geometry_art::skeleton`: the region within half a
bar of the graph is the domain less every cell's inset, so each cell
hands over the ring between its boundary and its inset and the builder
triangulates, offsets and walls it the same way everywhere. A geometry
takes part by outlining its cells -- the sphere with cap polygons whose
rims are small circles, the flat family through the embedding its
closure names -- and by nothing else.

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

- Walled directions need a **bounded domain**, which the sphere never
  did because it is compact. The plane needs bounds in both directions,
  a cylinder in one, a torus in neither -- so the torus is structurally
  the closest of the three to the sphere, and the default domain below.
- On a cylinder embedded in space, chord distance and geodesic distance
  disagree, unlike the sphere where both order the same way. The
  intrinsic one is both the tessellation people want and the one that
  unrolls, so the flat family means intrinsic throughout.

## Boundaries and edge treatments

The flat family's domain is a per-axis choice: each direction is
either **wrapped** or **walled**. Plane (walled, walled), cylinder
(wrapped, walled), torus (wrapped, wrapped). Wrapping is exact
periodicity in the diagram -- a cylinder's pattern must close
seamlessly around its circumference, and no margin can fake that.
A wall clips cells inside the optimizer and contributes its own
derivative terms: a cell edge on the wall does not sweep when sites
move, and a Voronoi vertex on the wall is constrained by the wall
plus one equidistance instead of three-site equidistance -- the slot
the sphere constraint occupies in the endpoint system holds the wall
equation instead. The flat diagram builds every cell the same way
whichever closure its axes have: the site's fundamental region --
the strip nearer the site than its own period images along a wrapped
axis, the domain itself along a walled one -- clipped by the
bisectors of its Delaunay neighbours, each edge remembering the cut
it lies on. A wall and a self-image seam are then the same kind of
edge, one that shapes the cell without sweeping, and the exact
constraint curvature reads the cut on the far side of each bisector
endpoint to know whether a third site or a wall pins it.

**The torus is the default domain.** It is the only member with no
boundary mathematics at all, so the optimizer stack ports unchanged;
it is statistically homogeneous, with no wall anywhere to nucleate
order; and every cell is usable pattern. The walled members extend
it rather than precede it. The plane exists for the densities a
wrapped axis cannot carry: a gradient rising across the width is not
periodic, so it needs a wall at each end.

The edge of a finished piece is a separate, artistic choice, and the
treatments decouple from the domain:

- **Cut.** The default. Clip the tessellation at the frame at render
  time. Border cells are sliced mid-cell -- partial cells are part of
  the look -- and the piece reads as a window onto a pattern that
  continues past the edge. From the torus this is free;
  from a walled domain it needs a margin between wall and frame. A
  cylinder whose axial period equals the frame height wastes nothing:
  each cell's clipped top fragment reappears as its bottom fragment,
  so the visible band holds every cell exactly once.
- **Conformed.** The frame is the wall. Every cell, border cells
  included, is a complete equal-capacity cell pressed flat against
  the edge, and the piece reads as a self-contained object.
- **Scalloped.** Keep only the cells lying wholly inside the frame,
  so the silhouette is the cells' own walls: an organic edge, and
  nothing but a filter on which cells to draw.
- **Closed cut.** Slice the border cells and seal them along the
  frame line, so every visible cell reads as complete with one flat
  side. This imitates the conformed edge, except the sealed cells
  hold less than a full capacity; only the true conformed edge needs
  the optimizer to know the wall.
- **Fade.** The density dips toward the line where the frame will
  cut, cells swell in the band, and the cut runs through the sparse
  zone. On the torus this is pure density design -- the domain is
  compact, so no wall is needed to keep masses finite -- and it
  wastes almost nothing, since the cells in the band are as few as
  the fade makes them. The dip bottoms out at a positive floor,
  never zero: the capacity constraint needs every cell to hold
  mass, and the positivity certification is precisely a proof the
  density has a floor. The band must also be at least a cell wide,
  or the bandwidth rule says the tessellation cannot read it.

The renderer owns each frame edge independently, so the treatments
mix per edge -- a cylinder cut at the top and scalloped at the
bottom. Every look except the true conformed edge comes from the
torus, in the density or the renderer: the optimizer never sees a
boundary, the density sees the frame only as a mask, and the frame
itself is a crop. The plane is the conformed edge on all four sides
and the cylinder is conformed at both rims, so a density can rise
up its height; the cut look survives as a way of viewing a torus,
rolled into a tube with the rims slicing through its cells.

One caution for margins: a straight wall nucleates rows of cells
aligned parallel to it, the way crystallization starts at the flat
side of a container, and the layering can propagate several rows
deep at near-constant density. A margin is therefore measured in
cell layers, not length, and its required depth should be measured
-- count the layers before cell-shape statistics match the interior
-- rather than assumed. Noise and image densities disrupt the
alignment naturally.

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

## What the second instance taught

The flat torus was built against the sphere and the seams settled
where the two instances agreed:

- The energy formula unified: the CVT energy reads the squared-norm
  moment, which on the unit sphere is the mass itself, so one formula
  serves both. The gradient's position term stayed per geometry -- the
  sphere's tangent projection discards what the flat domain keeps.
- The diagram is the geometry object: it owns construction, charts,
  rebuild-with-retraction, and the per-edge opposite-site data the
  exact curvature needs.
- Curvature blocks stayed per geometry, as they should: rotation of
  great circles against a spherical constraint on one side, straight
  lines with planar equidistance systems on the other, both held to
  the same finite-difference referee and both handing the shared
  loops the same block structure.
- Wrapped domains name bisectors by site pair AND period offset: the
  same two sites can share two bisectors, and a cell can border
  itself, in which case its two edge appearances cancel exactly.
