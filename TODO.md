- Make scalar fields act on Euclidean vectors
- Consider making CircularInterval::hull non-commutative
- Consider not providing concrete default template arguments
- Make our own precondition macro

## Direction

Ordered so that each stage is verifiable on its own and derisks the
next. The limits motivating the work are in
`docs/exactness_frontier.md`.

Where two items look equally valuable, prefer the one that carries over
to the plane, the torus and meshes; `docs/geometry_portability.md`
records what already does. That consideration retired the
spherical-harmonic basis: it addressed conditioning in high-degree
global fits, and the global path is the one that does not generalise.

### 1. Flat C1 elements
The flat densities are C0: quadratics through samples on a periodic
grid, refined to a requested tolerance by probing. The planar
Powell-Sabin element is the textbook case and the projection
machinery is already structured for it, so the sphere's certified,
least-squares C1 path has a direct flat analogue waiting.

*Verified by:* a certified positive floor and a telescoped residual
on a flat density, as on the sphere.

### 2. Walls
Deferred by design: every edge look except the true conformed one
comes from cropping the torus, so the bounded domain and its terms in
the exact Hessian wait until a piece demands what a cut cannot give.
The design is in `docs/geometry_portability.md`.

*Verified by:* a walled rectangle driving the shared optimizer stack.

### 3. Torus-aware density
A field on the closed minor dimension v, built from two independent
factors that multiply: the embedding correction — the area element
R + r cos phi of the ring torus, configured by the tube aspect r/R,
which makes cells look uniform on the embedded surface — and an
artistic profile in phi, mirror-symmetric (f(v) = f(h - v)) so it is
specified on half the circle and seamless by construction. Dense at
both equators and sparse between is one such profile, not a special
case. Parameters are in user-meaningful units
transformed internally: the aspect as r/R, the profile strength as
a densest-to-sparsest cell area ratio, phase in turns of phi --
never raw coefficients in chart coordinates.

*Verified by:* a snapshot on the embedded torus whose cell areas are
uniform under the identity profile, and match the profile otherwise.

### 4. Field parameters in the studio
Each field declares its parameters in user-meaningful units --
contrast as a densest-to-sparsest ratio, level, phase in turns --
and the applet renders controls from that declaration, so adding a
field never touches the viewer.

*Verified by:* a field with a parameter driven from the applet end
to end.
