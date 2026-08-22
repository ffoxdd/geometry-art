# Capacity-Constrained Centroidal Voronoi Tessellation on the Sphere

The generator places `n` sites on the unit sphere so that every ordinary
(unweighted) Voronoi cell holds the same mass of a density field ρ, while
each site sits at the ρ-weighted centroid of its cell as far as the
constraint allows. Cells are unweighted Voronoi cells rather than power
cells for aesthetic reasons, which is what makes the problem non-convex
and the capacity constraints non-trivial.

Reference: Balzer, Schlömer, Deussen, *Capacity-constrained point
distributions: a variant of Lloyd's method* (SIGGRAPH 2009); de Goes et
al., *Blue noise through optimal transport* (SIGGRAPH Asia 2012) for the
boundary-integral gradients.

## Formulation

With sites `s_i ∈ S²`, Voronoi cells `V_i`, density ρ and
`t = (∫_{S²} ρ) / n`:

```
minimize   E_cvt(s) = Σ_i ∫_{V_i} ρ(x) |x − s_i|² dA
subject to M_i(s) = ∫_{V_i} ρ(x) dA = t   for all i
```

Solved as an augmented Lagrangian

```
L(s, λ, μ) = E_cvt(s) + Σ_i λ_i (M_i − t) + (μ/2) Σ_i (M_i − t)²
```

with L-BFGS inner solves over unconstrained `x ∈ R^{3n}`, `s_i = x_i/|x_i|`,
multiplier updates `λ_i ← λ_i + μ (M_i − t)` and penalty growth when the
maximum violation fails to shrink by the required factor. A few passes
of density-weighted Lloyd relaxation (`s_i ← normalize(∫_{V_i} ρ x dA)`)
warm-start the solve.

## Gradients

All gradients are exact. With `F_i = ∫_{V_i} ρ x dA`:

- `∂E_cvt/∂s_i = −2 F_i` (the cell-boundary terms cancel because the
  integrand is continuous across every bisector).
- Moving `s_k` moves only the bisectors between `k` and its neighbours
  `j`, with normal velocity `(x·d)/|s_j − s_k|` along the shared arc, so

  ```
  ∂M_k/∂s_k = Σ_j (∫_{e_kj} ρ x ds − s_k ∫_{e_kj} ρ ds) / |s_j − s_k|
  ∂M_j/∂s_k = −(that edge's term)
  ```

  and the constraint part of `∇_{s_k} L` is
  `Σ_j (w_k − w_j) · edge_term(k, j)` with `w_i = λ_i + μ (M_i − t)`.
- The ambient gradient is mapped to the optimisation variables by the
  chain rule through normalisation: `(I − s sᵀ) g / |x|`.

## Exact integration

Every quantity above is an integral of a polynomial over a spherical
polygon or a great-circle arc. `Polygon::moments(d)` and
`Arc::moments(d)` return all monomial moments `∫ x^a y^b z^c` up to degree
`d` in closed form:

- Arc: parametrise `x(t) = cos t · u + sin t · n`, expand binomially and
  reduce `∫ cos^p t sin^q t dt` by recurrence.
- Polygon: area by Gauss–Bonnet (`2π − Σ turning angles`), then the
  Stokes recursion on S² with `ν` the outward conormal (constant along a
  great-circle arc):

  ```
  ∫_P x_k f dA = (∫_P ∂_k f dA − Σ_e ν_{e,k} ∫_e f ds) / (deg f + 2)
  ```

  which reduces a degree-`k` polygon moment to a degree-`(k−2)` polygon
  moment and degree-`(k−1)` arc moments.

A `PolynomialField` holds monomial coefficients; its mass and first
moment over a region are dot products with that region's moments, and its
total mass is the closed-form unit-sphere integral. Arbitrary scalar
fields are projected onto the `(degree, degree − 1)` homogeneous monomial
basis — the span of spherical harmonics up to that degree — by least
squares on Fibonacci sample points (`PolynomialFieldFitter`). High
frequency content beyond the chosen degree is not representable this
way; that is what the piecewise representation is for.

`PiecewisePolynomialField` carries a homogeneous degree-`d` polynomial
on every triangle of a spherical `TriangleMesh` (an icosphere), fitted
to the degree-`d` Lagrange nodes of the triangle — corners, edge points,
interior points projected to the sphere. A homogeneous degree-`d`
polynomial has exactly as many monomials as that lattice has nodes, its
trace along an edge is determined by the edge's nodes so the field is
continuous, and for even `d` it reproduces constants exactly. Cell and
arc integrals clip the region against each overlapping triangle
(Sutherland–Hodgman against great-circle half-spaces, candidates from a
kd-tree over triangle centres) and sum the per-triangle polynomial
integrals, so the integration stays exact for the interpolant.

## Geometry conventions

- An `Arc` is traversed counter-clockwise about its `normal()`; the
  two-argument constructor takes the minor arc, the three-argument form
  allows antipodal and major arcs.
- A `Polygon` is a closed loop of arcs oriented counter-clockwise seen
  from outside the sphere; its interior is to the left. An inside-out
  loop describes the complement, consistently for area, moments and
  `contains`.
- Voronoi vertices are normalised Delaunay face normals computed in
  double precision; the Delaunay combinatorics come from CGAL's exact
  predicates.

## Precision floor

The constraint terms compete with `E_cvt` inside one double-precision
objective, so the relative RMS capacity error bottoms out around 1e-9;
the optimizer reports `stalled` there. The default tolerance is 1e-7.
