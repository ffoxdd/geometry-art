# What the solver spends its time on

Measured on 2026-08-25, on one machine, at the site counts named. The
proportions come from sampling profiles, which measure shares inside a
process and so survive a loaded machine; the iteration counts are
deterministic given a seed and survive it too. Wall-clock numbers are
deliberately absent, because that machine was never quiet enough for
them to mean anything.

## The two stages of a trial step

A trial step rebuilds the diagram and integrates every cell. Which of
those dominates was, for a while, the difference between the two
geometries.

The sphere triangulates one point per site. The flat torus triangulated
nine, because every site was copied into all the surrounding tiles, and
each copy's location walk crossed a structure nine times too large.
Insertion work grows like the vertex count to the three-halves power, so
that is twenty-seven times the point location at equal site count, and
the copies of one site sit a period apart, which leaves the walk no
locality to exploit either. At four thousand sites it was 43% of the
solver's serial time against nothing measurable on the sphere.

Replicating only the sites a seam can reach took it to 3%. What remains
is what the sphere already spent its trial time on: integrating cells,
which parallelises.

The integration itself is cheaper on the flat domain -- planar segment
moments are a double loop where spherical arc moments are a triple loop
that also wants trigonometric integrals -- so with the triangulation no
longer in the way the flat family should stay ahead of the sphere at
equal site count.

## The conditioning, not the site count

The inner descent is a trust-region Newton method whose steps come from
a truncated conjugate gradient. Unpreconditioned, half of every step's
conjugate gradients were hitting their hundred-iteration cap: the
descent was working from steps the model had not finished computing and
paying for the shortfall in more steps.

Preconditioning by the operator's block diagonal, penalty included, at
a thousand sites on the cylinder: 30,260 Hessian-vector products across
518 Newton steps became 7,404 across 232.

## Coarse-to-fine continuation does not pay

Solving at a coarse site count and carrying the answer up by splitting
the heaviest cells was tried and removed. Weighting iterations by site
count, on the cylinder at eight thousand sites: one level cost 4.50
million, three levels 8.25, four levels 8.26.

The premise is that long-wavelength error costs about the same number of
iterations whatever the site count, so relaxing it where sites are few
is the same work at a discount. That premise is false here. The
single-level solve took 562 iterations at eight thousand sites and 725
at four thousand -- the count does not grow with the site count at all,
so there is nothing to amortise, and the coarse levels are pure
addition. The conditioning that would have produced such growth is
what the preconditioner already removed.

A refinement that seeded the fine level better than a Lloyd relaxation
does might change this. Splitting the heaviest cells is not that:
relaxing after each split was needed just to stop it being actively
harmful.
