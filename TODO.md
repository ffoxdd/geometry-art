- Make scalar fields act on Euclidean vectors
- Consider making CircularInterval::hull non-commutative
- Consider not providing concrete default template arguments
- Make our own precondition macro

## Future Features

### Piecewise density fields
Represent image-like densities exactly by clipping each Voronoi cell against a
piecewise-linear (or piecewise-constant) density defined on a spherical mesh, and
integrating each piece with the exact moment machinery. Global polynomial fits
cannot carry high-frequency content.
