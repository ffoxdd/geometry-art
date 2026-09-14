import * as THREE from 'three';

// A skin lofted through rings. The tree decomposes into paths -- at every
// branch the heaviest child carries the path on, and each other child
// starts a path of its own there, entering at its own thickness and
// plunging into the branch that carries it -- and each path is swept: a
// ring at each node, perpendicular to the path, sized to the node, joined
// to its neighbours by quads. Ring frames are carried along the path by
// parallel transport, so tubes do not twist; `relax` rounds off the
// particle-scale jitter of the bond directions by smoothing centres and
// radii along each path, ends pinned. Rings crowding closer than their
// own thickness are dropped, so the triangles follow the structure.

const SPACING = 0.6;
const SEGMENTS = { least: 8, most: 24 };

export function skin(pruned, tree, radius, at, appearance, shade, anchors = new Map()) {
  const built = build(rails(pruned, tree, radius, at, appearance, anchors), appearance, shade);

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.Float32BufferAttribute(built.positions, 3));
  geometry.setAttribute('color', new THREE.Float32BufferAttribute(built.colors, 3));
  geometry.setIndex(built.kept);
  geometry.computeVertexNormals();

  return {
    geometry,
    triangles: built.kept.length / 3,
    hidden: built.dropped.length / 3,
  };
}

// The whole skin, with its interior faces separated out. A face is
// interior when its three corners sit strictly inside one witness
// capsule: a capsule is convex, so three inside corners put the whole
// face inside -- no point of it can lie on the union's boundary, and
// dropping it can never open a hole. Each witness shrinks to stay inside
// the solid it stands for -- by the tilt of its rings, the bend at its
// joints, and pulled in from a path's capped ends -- so the cut is
// conservative: what is dropped is provably hidden, and what is merely
// hidden may survive.
export function build(railList, appearance, shade) {
  const out = { positions: [], colors: [], indices: [], owners: [], witnesses: [] };
  const finest = railList.reduce((least, rail) => Math.min(least, ...rail.radii), Infinity);

  railList.forEach((rail, railIndex) => {
    const segments = Math.min(
      SEGMENTS.most,
      Math.max(SEGMENTS.least, Math.ceil(Math.PI * Math.max(...rail.radii) / finest))
    );

    sweep(rail, railIndex, segments, appearance, shade, out);
  });

  const sifted = culled(out);

  return {
    positions: out.positions,
    colors: out.colors,
    indices: out.indices,
    kept: sifted.kept,
    dropped: sifted.dropped,
    interior: sifted.interior,
    owners: out.owners,
  };
}

// The rings the skin passes through, path by path. A carrier is always
// railed before the paths that hang from it, and what it records in
// `anchors` is its centreline as smoothing left it -- which is what a
// side path must attach to, not where the raw particle sat.
export function rails(pruned, tree, radius, at, appearance, anchors = new Map()) {
  return paths(pruned, tree, radius)
    .map(trail => railed(trail, at, radius, appearance.relax ?? 2, anchors))
    .filter(rail => rail !== null && rail.nodes.length >= 2);
}

// The trunk decomposition: every shown node lies on exactly one path, and
// a path begins either at a root or at the branch node it hangs from.
export function paths(pruned, tree, radius) {
  const children = new Map();

  for (const edge of pruned.edges) {
    if (children.has(edge.parent)) children.get(edge.parent).push(edge.child);
    else children.set(edge.parent, [edge.child]);
  }

  for (const list of children.values()) {
    list.sort((first, second) => radius[second] - radius[first]);
  }

  const trails = [];
  const pending = pruned.shown
    .filter(node => tree.particles[node].parent == null)
    .map(root => ({ nodes: [root], side: false }));

  while (pending.length > 0) {
    const trail = pending.pop();
    let node = trail.nodes[trail.nodes.length - 1];

    for (let kids = children.get(node); kids !== undefined; kids = children.get(node)) {
      for (let extra = 1; extra < kids.length; extra++) {
        pending.push({ nodes: [node, kids[extra]], side: true });
      }

      trail.nodes.push(kids[0]);
      node = kids[0];
    }

    trails.push(trail);
  }

  return trails;
}

function railed(trail, at, radius, relax, anchors) {
  let nodes = trail.nodes;
  let centres = nodes.map(node => {
    const place = at(node);
    return [place.x, place.y, place.z];
  });

  let radii = nodes.map(node => radius[node]);

  // A side path enters at its own thickness, from the carrier as it is
  // actually drawn: its start snaps to the carrier's smoothed centreline,
  // and only rings past the carrier's surface survive -- the buried ones
  // are the union's interior and are cut away. A path that never emerges
  // is not drawn at all, its mass already in the carrier's own thickness.
  // The entry ring sits just under the surface, so the weld never gapes.
  if (trail.side) {
    const anchor = anchors.get(nodes[0]) ?? { centre: centres[0], radius: radii[0] };
    centres[0] = [...anchor.centre];
    radii[0] = radii[1];

    const buried = 0.85 * anchor.radius;

    let out = 1;
    while (out < centres.length && distance(centres[out], centres[0]) <= buried) out++;

    // Drawn or cut away as the union's interior, a node still records
    // where things hanging from it attach: a buried node's branches and
    // buds belong to the carrier's own centreline.
    if (out === centres.length) {
      for (let i = 1; i < nodes.length; i++) {
        if (!anchors.has(nodes[i])) anchors.set(nodes[i], { centre: anchor.centre, radius: radii[i] });
      }

      return null;
    }

    if (out > 1) {
      const entry = touched(centres[0], centres[out - 1], centres[out], buried);

      for (let i = 1; i < out - 1; i++) {
        if (!anchors.has(nodes[i])) anchors.set(nodes[i], { centre: entry, radius: radii[i] });
      }

      nodes = nodes.slice(out - 1);
      centres = centres.slice(out - 1);
      radii = radii.slice(out - 1);
      centres[0] = entry;
      radii[0] = radii[1];
    }
  }

  for (let pass = 0; pass < relax; pass++) {
    for (let i = 1; i < centres.length - 1; i++) {
      for (let axis = 0; axis < 3; axis++) {
        centres[i][axis] = 0.5 * centres[i][axis] + 0.25 * (centres[i - 1][axis] + centres[i + 1][axis]);
      }

      radii[i] = 0.5 * radii[i] + 0.25 * (radii[i - 1] + radii[i + 1]);
    }
  }

  // Recorded after smoothing and before thinning: later branches attach
  // to this, the path as drawn.
  for (let i = 0; i < nodes.length; i++) {
    if (!anchors.has(nodes[i])) {
      anchors.set(nodes[i], { centre: centres[i], radius: radii[i] });
    }
  }

  return thinned(nodes, centres, radii);
}

// Where the walk from `from` to `to` first reaches `reach` from the
// origin, taken linearly along the step.
function touched(origin, from, to, reach) {
  const near = distance(from, origin);
  const far = distance(to, origin);
  const along = far - near < 1e-12 ? 0 : Math.min(1, Math.max(0, (reach - near) / (far - near)));

  return [0, 1, 2].map(axis => from[axis] + along * (to[axis] - from[axis]));
}

function thinned(nodes, centres, radii) {
  const keep = [0];

  for (let i = 1; i < centres.length; i++) {
    const last = keep[keep.length - 1];

    if (i === centres.length - 1 || distance(centres[i], centres[last]) >= SPACING * radii[last]) {
      keep.push(i);
    }
  }

  return {
    nodes: keep.map(i => nodes[i]),
    centres: keep.map(i => centres[i]),
    radii: keep.map(i => radii[i]),
  };
}

function sweep(rail, railIndex, segments, appearance, shade, out) {
  const { positions, colors, indices, owners } = out;
  const { nodes, centres, radii } = rail;
  const tangents = centres.map((centre, i) => direction(
    centres[Math.max(0, i - 1)],
    centres[Math.min(centres.length - 1, i + 1)]
  ));

  let normal = perpendicular(tangents[0]);
  const first = positions.length / 3;

  for (let i = 0; i < centres.length; i++) {
    if (i > 0) normal = transported(centres[i - 1], centres[i], tangents[i - 1], tangents[i], normal);

    const binormal = cross(tangents[i], normal);
    const colour = shade(nodes[i]);

    for (let s = 0; s < segments; s++) {
      const angle = (2 * Math.PI * s) / segments;
      const sway = Math.cos(angle);
      const lean = Math.sin(angle);

      for (let axis = 0; axis < 3; axis++) {
        positions.push(centres[i][axis] + radii[i] * (sway * normal[axis] + lean * binormal[axis]));
      }

      colors.push(colour.r, colour.g, colour.b);
    }

    if (i === 0) continue;

    const below = first + (i - 1) * segments;
    const here = first + i * segments;

    for (let s = 0; s < segments; s++) {
      const next = (s + 1) % segments;
      indices.push(below + s, below + next, here + next);
      indices.push(below + s, here + next, here + s);
      owners.push(railIndex, i - 1, railIndex, i - 1);
    }
  }

  cap(rail, 0, tangents[0], -1, segments, first, appearance, shade, out, railIndex);
  cap(rail, centres.length - 1, tangents[centres.length - 1], 1, segments, first, appearance, shade, out, railIndex);
  witnessed(rail, tangents, railIndex, out.witnesses);
}

// A path's open ends are closed by a fan to an apex: flat on the ring for
// a flat cap and every path's start (which hides inside what it grew
// from), one radius out for a round one, `TIP_CONE` radii for a cone.
function cap(rail, ring, tangent, facing, segments, first, appearance, shade, out, railIndex) {
  const { positions, colors, indices, owners } = out;
  const reach = facing < 0 || appearance.cap === 'flat' ? 0 : (appearance.cap === 'cone' ? 2.5 : 1);
  const centre = rail.centres[ring];
  const apex = positions.length / 3;
  const colour = shade(rail.nodes[ring]);
  const segment = ring === 0 ? 0 : rail.centres.length - 2;

  for (let axis = 0; axis < 3; axis++) {
    positions.push(centre[axis] + facing * reach * rail.radii[ring] * tangent[axis]);
  }

  colors.push(colour.r, colour.g, colour.b);

  const base = first + ring * segments;

  for (let s = 0; s < segments; s++) {
    const next = (s + 1) % segments;

    if (facing > 0) indices.push(base + s, base + next, apex);
    else indices.push(base + next, base + s, apex);

    owners.push(railIndex, segment);
  }
}

// The witnesses: one conservative frustum per drawn segment, provably
// inside the solid the segment draws. A point counts only strictly
// between the segment's own ring planes, so a witness can never claim
// ground past its rings, where its rail's wall has bent away; within the
// slab, the radii shrink by the tilt of each ring against the chord and
// by the bend at each joint, which is the tangency condition against the
// dented wall.
//
// A tube thinner than a few pixels at the standard framing witnesses
// nothing: a face hidden inside one is beneath notice, skipping a witness
// only ever keeps more, and leaving the hairline canopy out keeps the cut
// affordable at full detail.
const FAINT = 0.006;

function witnessed(rail, tangents, railIndex, witnesses) {
  const { centres, radii } = rail;
  const last = centres.length - 2;
  const chords = [];

  for (let i = 0; i <= last; i++) chords.push(direction(centres[i], centres[i + 1]));

  for (let i = 0; i <= last; i++) {
    const chord = chords[i];
    const bendA = i === 0 ? 1 : Math.max(0, dot(chords[i - 1], chord));
    const bendB = i === last ? 1 : Math.max(0, dot(chord, chords[i + 1]));
    const ra = radii[i] * Math.max(0, dot(tangents[i], chord)) * bendA * 0.999;
    const rb = radii[i + 1] * Math.max(0, dot(tangents[i + 1], chord)) * bendB * 0.999;

    if (ra < 1e-12 || rb < 1e-12) continue;
    if (Math.max(ra, rb) < FAINT) continue;

    witnesses.push({ a: centres[i], b: centres[i + 1], ra, rb, rail: railIndex, segment: i });
  }
}

// The cut itself: a face whose three corners all sit strictly inside one
// witness of ANOTHER path is interior, and goes. A path never witnesses
// against itself -- a tube cannot cut its own wall -- so a fold's
// self-overlap keeps its hidden faces, which no outside eye sees anyway.
function culled({ positions, indices, owners, witnesses }) {
  const buckets = crated(witnesses);
  const kept = [];
  const dropped = [];
  const interior = new Uint8Array(indices.length / 3);

  for (let triangle = 0; triangle < indices.length / 3; triangle++) {
    const a = indices[triangle * 3];
    const b = indices[triangle * 3 + 1];
    const c = indices[triangle * 3 + 2];
    const listed = buckets.get(crateOf(positions[a * 3], positions[a * 3 + 1], positions[a * 3 + 2]));

    let hidden = false;

    if (listed) {
      const rail = owners[triangle * 2];

      for (const witness of listed) {
        if (witness.rail === rail) continue;

        if (within(witness, positions, a) && within(witness, positions, b) && within(witness, positions, c)) {
          hidden = true;
          break;
        }
      }
    }

    interior[triangle] = hidden ? 1 : 0;
    (hidden ? dropped : kept).push(a, b, c);
  }

  return { kept, dropped, interior };
}

function within(witness, positions, vertex) {
  const x = positions[vertex * 3];
  const y = positions[vertex * 3 + 1];
  const z = positions[vertex * 3 + 2];
  const spanX = witness.b[0] - witness.a[0];
  const spanY = witness.b[1] - witness.a[1];
  const spanZ = witness.b[2] - witness.a[2];
  const squared = spanX * spanX + spanY * spanY + spanZ * spanZ;

  if (squared < 1e-18) return false;

  const along =
    ((x - witness.a[0]) * spanX + (y - witness.a[1]) * spanY + (z - witness.a[2]) * spanZ) / squared;

  if (along <= 0 || along >= 1) return false;

  const nearX = x - (witness.a[0] + along * spanX);
  const nearY = y - (witness.a[1] + along * spanY);
  const nearZ = z - (witness.a[2] + along * spanZ);
  const reach = witness.ra + (witness.rb - witness.ra) * along;

  return nearX * nearX + nearY * nearY + nearZ * nearZ < reach * reach;
}

const CRATES = 48;

function crated(witnesses) {
  const buckets = new Map();

  for (const witness of witnesses) {
    const widest = Math.max(witness.ra, witness.rb);
    const low = axis => crate(Math.min(witness.a[axis], witness.b[axis]) - widest);
    const high = axis => crate(Math.max(witness.a[axis], witness.b[axis]) + widest);

    for (let z = low(2); z <= high(2); z++) {
      for (let y = low(1); y <= high(1); y++) {
        for (let x = low(0); x <= high(0); x++) {
          const key = (z * CRATES + y) * CRATES + x;

          if (buckets.has(key)) buckets.get(key).push(witness);
          else buckets.set(key, [witness]);
        }
      }
    }
  }

  return buckets;
}

function crateOf(x, y, z) {
  return (crate(z) * CRATES + crate(y)) * CRATES + crate(x);
}

function crate(value) {
  return Math.max(0, Math.min(CRATES - 1, Math.floor((value + 1) * CRATES / 2)));
}

function transported(from, to, tangentFrom, tangentTo, normal) {
  const lead = [to[0] - from[0], to[1] - from[1], to[2] - from[2]];
  const led = dot(lead, lead);

  if (led < 1e-18) return normal;

  const reflectedNormal = reflected(normal, lead, led);
  const reflectedTangent = reflected(tangentFrom, lead, led);
  const settle = [
    tangentTo[0] - reflectedTangent[0],
    tangentTo[1] - reflectedTangent[1],
    tangentTo[2] - reflectedTangent[2],
  ];
  const settled = dot(settle, settle);

  if (settled < 1e-18) return reflectedNormal;

  return normalised(reflected(reflectedNormal, settle, settled));
}

function reflected(vector, mirror, squared) {
  const lean = (2 * dot(vector, mirror)) / squared;

  return [
    vector[0] - lean * mirror[0],
    vector[1] - lean * mirror[1],
    vector[2] - lean * mirror[2],
  ];
}

function perpendicular(tangent) {
  const shy = [Math.abs(tangent[0]), Math.abs(tangent[1]), Math.abs(tangent[2])];
  const axis = shy[0] <= shy[1] && shy[0] <= shy[2] ? [1, 0, 0] : (shy[1] <= shy[2] ? [0, 1, 0] : [0, 0, 1]);
  const held = dot(axis, tangent);

  return normalised([
    axis[0] - held * tangent[0],
    axis[1] - held * tangent[1],
    axis[2] - held * tangent[2],
  ]);
}

function direction(from, to) {
  return normalised([to[0] - from[0], to[1] - from[1], to[2] - from[2]]);
}

function normalised(vector) {
  const length = Math.hypot(vector[0], vector[1], vector[2]) || 1;

  return [vector[0] / length, vector[1] / length, vector[2] / length];
}

function cross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

function dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function distance(a, b) {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}