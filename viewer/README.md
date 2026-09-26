# Geometry Art studio

`./studio` serves the page and runs the programs behind it.

```
./studio            # http://localhost:8731/
./studio --help     # --port, --binary, --aggregate-binary, --skeletonize-binary, --runs
```

Standard library only; needs a built `tessellate` and/or `aggregate`
(defaults under `build-release/`). A program that is not built is not
offered.

Runs are executed one at a time, each kept under `runs/<id>/` with its
parameters, log and snapshot. The snapshot is rewritten every second while
a run is in flight, so the page redraws the structure as it grows; a run
can be cancelled while running and deleted afterwards.

## Parameters

Each program declares its parameters once, in `server.py`: the flag the
parameter becomes, how to present it, and the conditions under which it
applies at all -- `when` for the parameter, `choice_when` for one of its
choices. The page builds its form from that declaration and the server
validates against it, so a control that cannot apply is never offered, and
teaching the studio a new flag means adding it in one place.

## Exporting a model

A finished tessellation can leave the studio as a solid: `skeletonize`
thickens its edge graph into bars and the Export section, shown for a
selected `tessellate` run, downloads the result in STL, OBJ, PLY or OFF.
The bars are set in output units -- width along the surface, thickness off
it -- with `scale` saying how many of those units one model unit is: the
sphere's radius, or one unit of a flat rectangle. `resolution` bounds how
long a facet may run along a curved surface. A `window` cuts the model
down to a centred square of that width, laid flat, with a frame of
`frame width` around it; left blank, the whole domain is exported and the
frame defaults to the bar width. While a window is set, the drawing
outlines it on the tessellation: the inner line is the window's edge, where
the openings stop, and the outer one is the frame's, where the part ends.
Like a program, an export
declares its parameters once in `server.py`, and the page renders them.

## The page

`index.html` is markup and `studio.css` is the look; everything else is ES
modules beside them, split the way the screen is: the scene and the drawings
it holds on one side, the panels -- view, parameter form, runs -- on the
other. Both panels are built by the same renderer from a declaration, the
run form's coming from the server and each drawing's from the drawing
itself, so a new knob is declared rather than wired up.

The figure turns, not the camera: drag to roll it about any axis, and a
flick hands the idle spin that axis and speed, easing back to a slow
drift. Scroll to zoom -- the camera goes all the way inside -- and click
a rail's title to fold that panel away.

The view alone works without the server. Open the page from any static
server and drop a snapshot's `<name>.json` onto it, or address a reading
directly:

```
index.html?snapshot=<file>&colour=<name>&weight=<n>&ramp=<stops>&run=<id>
```

## Drawing an aggregate as a tree

Every particle but the seed froze onto exactly one parent, so an aggregate
is a tree and can be drawn as one: a truncated cone per edge, running from
its parent's thickness down to its own so that branches meet flush.
Thickness follows descendant count between two ends you set directly: a
tip is exactly `tip width` and the root exactly `root width`, both in
particle radii, however large the aggregate grows, and `taper` only bends
the path between those pinned ends -- above one the tree stays slender
and swells late, below one it thickens straight away.

Every node ends in one of three fates -- staying in the trunk, gathered
into a bud, or culled outright -- chosen by a pipeline of small rules
over the node's own attributes. A wish rule asks who leaves the trunk,
by one of two measures: `bud depth` gates on links from the tips, and
`bud under` gates on descendant count, gathering exactly the crowns
lighter than it says; `depth by age` sways either gate by growth age, so
the young or the old side of the aggregate is pruned harder.
Connectivity then makes any wish honest -- a node leaves only when
everything below it leaves, so what is drawn is always one tree. Last, a
fate rule says which pruned clusters stand as buds; the rest are culled
and simply gone. It can bud everywhere, cull everything, keep only
clusters above a size bar as buds, or split by place: a node with at
least `trunk from` descendants is trunk, and clusters carried within
`within links` of one are culled while the rest bud.

`by rules` replaces the named fates with a program: an ordered list of
assignments in a tiny closed language -- arithmetic, comparisons, logic,
a short function table, nothing that can run code -- evaluated pointwise
at every node. The sources are the tree's structural quantities, one
pass each and normalised where a natural scale exists: `mass` (subtree
size) and `mass_fraction` of the whole, `attachment_ratio` (a branch's
share of its parent, scale-free tininess), `tip_links`, `root_links`,
`spine_links` (distance to the maximal-mass path from the root), `age`,
and `cluster_mass` (a carrier's pruned cluster size -- 0 while `pruned`
is being decided, since clusters do not exist yet). The sinks are the
two fate decisions: assign `pruned` to say who leaves the trunk and
`budded`, read at each carrying node, to say which clusters stand as
buds; unassigned, nothing is pruned and every cluster buds. Any other
name is a free intermediate, so a program is the linearisation of a
dependency graph. Connectivity still binds every program, and one that
does not compile prunes nothing and says why in the readout.

`thinnest branch` is the drawable floor: below it a twig dissolves to
sub-pixel width under antialiasing and everything it carries appears to
float, so no branch renders thinner.

A bud is sized to hold its cluster's volume and stretched so that its
own centre of mass lands where the cluster's does.

A cone opens outwards the way a spray widens, a spike points the other way,
a ball forgets the axis and keeps only the mass. `bud aspect` stretches a
bud along its axis at constant volume, and `bud base` sets how far the
dome over a pointed bud's flat face rises, from nearly flat to a full
half-ball.

The tree takes one of two skins. `parts` assembles instanced primitives --
fast enough to redraw at any size, and it holds one invariant everywhere:
a flat rim is either flush with a ball's equator or tucked inside a ball,
never proud of one. Stems recess into the joint that covers their parent
end rather than stand out of it, and a dome caps each pointed bud's flat
face -- the base a spike stands on, the mouth a cone opens to.

`lofted rings` sweeps the skin through the structure itself. The tree
decomposes into paths -- at every branch point the heaviest child carries
the path on, so the trunk's skin runs unbroken through its junctions,
and each other child starts a path of its own -- and each path is swept:
a ring at every node, perpendicular to the path and sized to it, joined
to its neighbours, frames carried by parallel transport so tubes do not
twist. A side path enters at the carrier's surface and its rings inside
the carrier are cut away -- the union of the tubes, taken by
construction -- and a twig that never leaves its carrier's inside is not
drawn at all, its mass already in the carrier's thickness. `smoothing`
relaxes centres and radii along each path, rounding off the
particle-scale jitter of the bonds.

`smooth envelope` gives up primitives for an implicit surface: every part
of the trunk reduces to a tapered capsule, each capsule splats a soft
kernel into a scalar field, and marching cubes pulls the level set out as
one continuous skin -- welded where parts meet, swelling where they
crowd, with exact per-edge taper and vertex colours blended from whatever
welded there. The grid decides only the topology: every vertex is then
projected onto the exact isosurface and shaded by the field's own
gradient, so neither the geometry nor the lighting carries the lattice.
Buds stay out of the envelope, drawn as their crisp instanced shapes on
top of it. `blend` sets the kernel's reach past each capsule's own
radius: how far apart two parts can sit and still fuse, and the most the
skin may inflate where they pile up. `skin grid` trades resolution for
time, and sets the floor of the drawable: a part thinner than a cell is
drawn at the thinnest size the grid resolves, so fine twigs thicken to
the floor rather than shatter. A running solver draws parts and skins on
settling.

## Colour

Every part of a drawing is painted from its own source, and each source is
either a colour picked for it or the colour map read at the value that part
stands for. Branches and buds are separate parts, so a bud can be its own
colour, follow its branches, or be read off the map while the branches stay
picked. Finish -- satin, matte or metal -- is a choice separate from colour
and each part takes its own, so the buds can be gold metal on a tree of
satin bone.

A source that reads a value -- growth age, capacity error, cell area --
draws it through an editable map: stops on a bar, each a position and a colour, read linearly between neighbours. Click
the bar to add a stop, drag to move it, click it to recolour, right-click
to remove. Stops interpolate in linear light, and the bar is painted with
the same sampler the drawing uses, so it shows what it will get. The map
travels in the address as `position-rrggbb` pairs.
