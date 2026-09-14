import * as THREE from 'three';
import { group } from './scene.js';
import { ballDetail } from './detail.js';
import * as envelope from './envelope.js';
import * as finish from './finish.js';
import * as loft from './loft.js';
import * as rules from './rules.js';

// An aggregate is a tree: every particle but the seed froze onto exactly
// one parent. Drawn as one, a branch gets the thickness its subtree has
// earned, and the last generations before each tip can be gathered into a
// bud, leaving the structure rather than the spray. Nothing is thrown away
// by gathering: budded particles still count towards the thickness of the
// branch that carries them, and towards the bud that stands for them.

const UP = new THREE.Vector3(0, 1, 0);
const AROUND = 7;
const TIP_CONE = 2.5;
const RATIO_STEP = Math.log(1.1);

export function draw(snapshot, appearance, paint) {
  const scale = 0.98 / snapshot.reach;
  const tree = read(snapshot.particles);
  const pruned = prune(tree, appearance);
  const radius = floored(
    radii(tree.mass, appearance, snapshot.particleRadius * scale),
    appearance
  );
  const at = index => placeOf(snapshot.particles[index].center, scale);
  const grain = (4 / 3) * Math.PI * Math.pow(snapshot.particleRadius * scale, 3);
  const smooth = appearance.skin === 'smooth';
  const drawn = { tree, pruned, radius, at, scale, grain };

  if (appearance.skin === 'rings') {
    return ringed(snapshot, appearance, paint, drawn);
  }

  if (smooth && !appearance.live) {
    return skinned(snapshot, appearance, paint, drawn);
  }

  const branch = {
    material: new THREE.MeshPhysicalMaterial({ color: 0xffffff, ...finish.of(appearance.finish) }),
    shade: paint.branch,
  };

  const parts = [
    ...stems(pruned.edges, radius, at).map(part => ({ ...part, ...branch })),
    { placements: joints(pruned, radius, at, appearance.cap), geometry: () => ball(pruned.shown.length), ...branch },
    { placements: spikes(pruned, tree, radius, at, appearance.cap), geometry: () => spike(), ...branch },
  ];

  const resources = [branch.material];

  for (const part of parts.filter(part => part.placements.length > 0)) {
    resources.push(...instances(part.geometry(), part.material, part.placements, part.shade));
  }

  resources.push(...blooming(pruned, at, scale, grain, appearance, paint));

  return {
    resources,
    readout: [
      ['particles', snapshot.particles.length],
      ['branches', pruned.edges.length],
      ['budded', pruned.budded],
      ['culled', pruned.culled],
      ['reach', snapshot.reach.toFixed(1)],
    ],
    note: [pruned.note, smooth ? 'parts while running' : ''].filter(Boolean).join('; '),
  };
}

// The ringed skin lofts one mesh through cross-section rings along every
// path of the trunk; buds stay their crisp instanced selves, hung on the
// smoothed centreline the loft actually drew rather than the raw node it
// moved away from.
function ringed(snapshot, appearance, paint, drawn) {
  const { tree, pruned, radius, at, scale, grain } = drawn;

  const anchors = new Map();
  const lofted = loft.skin(pruned, tree, radius, at, appearance, paint.branch, anchors);
  const material = new THREE.MeshPhysicalMaterial({ vertexColors: true, ...finish.of(appearance.finish) });
  group.add(new THREE.Mesh(lofted.geometry, material));

  const anchored = node => {
    const held = anchors.get(node);
    return held ? new THREE.Vector3(held.centre[0], held.centre[1], held.centre[2]) : at(node);
  };

  return {
    resources: [
      lofted.geometry,
      material,
      ...blooming(pruned, anchored, scale, grain, appearance, paint),
    ],
    readout: [
      ['particles', snapshot.particles.length],
      ['branches', pruned.edges.length],
      ['budded', pruned.budded],
      ['culled', pruned.culled],
      ['triangles', lofted.triangles],
      ['interior cut', lofted.hidden],
      ['reach', snapshot.reach.toFixed(1)],
    ],
    note: pruned.note,
  };
}

// The envelope reduces the trunk to one primitive, the tapered capsule:
// stems taper exactly from parent to child with no quantised buckets, and
// capsule ends round every joint and tip on their own. Buds stay out of
// the field, drawn as their crisp instanced shapes on top of the skin.
function skinned(snapshot, appearance, paint, drawn) {
  const { tree, pruned, radius, at, scale, grain } = drawn;

  const capsules = [
    ...pruned.edges.map(edge => ({
      a: at(edge.parent),
      b: at(edge.child),
      ra: radius[edge.parent],
      rb: radius[edge.child],
      colour: paint.branch(edge.child),
    })),

    ...tipCones(pruned, tree, radius, at, appearance).map(cone => ({ ...cone, colour: paint.branch(cone.node) })),
  ];

  const skinned = envelope.skin(capsules, appearance);

  return {
    resources: [
      ...skinned.resources,
      ...blooming(pruned, at, scale, grain, appearance, paint),
    ],
    readout: [
      ['particles', snapshot.particles.length],
      ['branches', pruned.edges.length],
      ['budded', pruned.budded],
      ['culled', pruned.culled],
      ['triangles', skinned.triangles],
      ['reach', snapshot.reach.toFixed(1)],
    ],
    note: [pruned.note, skinned.truncated ? 'skin truncated -- coarsen the grid' : ''].filter(Boolean).join('; '),
  };
}

function tipCones(pruned, tree, radius, at, appearance) {
  if (appearance.cap !== 'cone') return [];

  return pruned.tips
    .filter(node => tree.particles[node].parent != null)
    .map(node => {
      const base = at(node);
      const heading = direction(at(tree.particles[node].parent), base);

      return {
        node,
        a: base,
        b: new THREE.Vector3().copy(base).addScaledVector(heading, TIP_CONE * radius[node]),
        ra: radius[node],
        rb: 0,
      };
    });
}

// The buds as both skins draw them: the shapes and their covering domes,
// instanced, in the bud's own finish and paint.
function blooming(pruned, at, scale, grain, appearance, paint) {
  const bloom = budding(blooms(pruned, at, scale, grain, appearance));
  const material = new THREE.MeshPhysicalMaterial({ color: 0xffffff, ...finish.of(appearance.bud_finish) });
  const resources = [material];

  const parts = [
    { placements: bloom.shapes, geometry: () => budShape(pruned, appearance) },
    { placements: bloom.caps, geometry: () => ball(pruned.buds.length) },
  ];

  for (const part of parts.filter(part => part.placements.length > 0)) {
    resources.push(...instances(part.geometry(), material, part.placements, paint.bud));
  }

  return resources;
}

function instances(geometry, material, placements, shade) {
  const mesh = new THREE.InstancedMesh(geometry, material, placements.length);

  placements.forEach((placement, index) => {
    mesh.setMatrixAt(index, placement.matrix);
    mesh.setColorAt(index, shade(placement.node));
  });

  group.add(mesh);

  return [geometry, mesh];
}

// One truncated cone per edge, running from its parent's thickness at one
// end down to its own at the other, so branches meet flush instead of
// stepping. Instances of a shared geometry can only be scaled uniformly
// across, so taper ratios are quantised to 10% steps and edges near the
// same ratio share a geometry -- rounded down, keeping the invariant that
// makes the skin read as one surface: every rim is flush with a ball's
// equator or inside a ball, never proud of one. A recess of at most one
// step hides inside the joint that covers it.
function stems(edges, radius, at) {
  const buckets = new Map();

  for (const edge of edges) {
    const ratio = quantised(radius[edge.parent] / radius[edge.child]);

    const placement = {
      node: edge.child,
      matrix: spanning(at(edge.parent), at(edge.child), radius[edge.child]),
    };

    if (buckets.has(ratio)) buckets.get(ratio).push(placement);
    else buckets.set(ratio, [placement]);
  }

  return [...buckets].map(([ratio, placements]) => ({ placements, geometry: () => stem(ratio) }));
}

function quantised(ratio) {
  return Math.exp(Math.floor(Math.log(ratio) / RATIO_STEP) * RATIO_STEP);
}

// A ball at every node rounds the corner between two stems. Whether the
// tips get one is the cap: round keeps them, flat and cone do not.
function joints(pruned, radius, at, cap) {
  const tips = new Set(pruned.tips);
  const capped = cap === 'round' ? pruned.shown : pruned.shown.filter(node => !tips.has(node));

  return capped.map(node => ({ node, matrix: sitting(at(node), radius[node]) }));
}

function spikes(pruned, tree, radius, at, cap) {
  if (cap !== 'cone') return [];

  return pruned.tips
    .filter(node => tree.particles[node].parent != null)
    .map(node => ({
      node,
      matrix: pointing(at(node), direction(at(tree.particles[node].parent), at(node)), radius[node], TIP_CONE * radius[node]),
    }));
}

// A cluster of budded leaves is drawn as the shape it would fill: hung off
// the branch carrying them, along the axis to their centre of mass, sized
// to hold the same volume as the particles it stands for and stretched so
// that its own centre of mass lands where theirs does.
//
//   spike  points away from the branch; a cone's centre of mass sits a
//          quarter of the way up from its base, so it runs to 4 times the
//          reach
//   cone   the same volume opening outwards instead, meeting the branch at
//          a point, so it runs to 4/3 of the reach
//   ball   no axis, just the volume centred on the mass -- stretched into
//          an ellipsoid only when the sphere would drift clear of its
//          branch, just far enough to keep touching it
//
// `bud_aspect` stretches a spike or cone along its axis at constant volume,
// letting the centre of mass slide in exchange.
//
function blooms(pruned, at, scale, grain, appearance) {
  const forms = [];

  if (appearance.bud === 'none' || appearance.bud === undefined) return forms;

  for (const bud of pruned.buds) {
    const root = at(bud.node);
    const centre = placeOf(bud.centroid, scale);
    const axis = new THREE.Vector3().subVectors(centre, root);
    const reach = axis.length();
    const volume = appearance.bulk * bud.count * grain;
    const round = Math.cbrt(3 * volume / (4 * Math.PI));

    if (appearance.bud === 'ball' && reach <= round) {
      forms.push({ node: bud.node, kind: 'ball', centre, round });
      continue;
    }

    if (reach < 1e-9) continue;

    const heading = axis.divideScalar(reach);

    if (appearance.bud === 'ball') {
      forms.push({
        node: bud.node,
        kind: 'stretched',
        centre,
        heading,
        along: reach,
        across: Math.sqrt(3 * volume / (4 * Math.PI * reach)),
      });
      continue;
    }

    const length = reach * (appearance.bud === 'cone' ? 4 / 3 : 4) * (appearance.bud_aspect ?? 1);
    const width = Math.sqrt(3 * volume / (Math.PI * length));

    forms.push({
      node: bud.node,
      kind: appearance.bud,
      root,
      heading,
      length,
      width,
      dome: Math.max(appearance.bud_base ?? 0.35, 0.02) * width,
    });
  }

  return forms;
}

// A dome caps each flat face -- the base a spike stands on, the mouth a
// cone opens to -- so a bud meets its branch and the world with curvature
// everywhere. `bud base` sets how far the dome rises over the face, from
// nearly flat to a full half-ball. The domes ride outside the volume
// match, small against the shapes they cap.
function budding(forms) {
  const shapes = [];
  const caps = [];

  for (const form of forms) {
    if (form.kind === 'ball') {
      shapes.push({ node: form.node, matrix: sitting(form.centre, form.round) });
      continue;
    }

    if (form.kind === 'stretched') {
      shapes.push({ node: form.node, matrix: holding(form.centre, form.heading, form.across, form.along) });
      continue;
    }

    if (form.kind === 'cone') {
      shapes.push({ node: form.node, matrix: opening(form.root, form.heading, form.width, form.length) });
      caps.push({
        node: form.node,
        matrix: holding(
          new THREE.Vector3().copy(form.root).addScaledVector(form.heading, form.length),
          form.heading, form.width, form.dome
        ),
      });
      continue;
    }

    shapes.push({ node: form.node, matrix: pointing(form.root, form.heading, form.width, form.length) });
    caps.push({ node: form.node, matrix: holding(form.root, form.heading, form.width, form.dome) });
  }

  return { shapes, caps };
}

const scratch = {
  middle: new THREE.Vector3(),
  axis: new THREE.Vector3(),
  size: new THREE.Vector3(),
  turn: new THREE.Quaternion(),
};

const UNTURNED = new THREE.Quaternion();

// A unit stem stands along +Y from -0.5 to 0.5, so an edge is one instance
// of it, turned onto the edge and stretched to its length.
function spanning(from, to, radius) {
  const axis = scratch.axis.subVectors(to, from);
  const length = axis.length();

  scratch.turn.setFromUnitVectors(UP, axis.divideScalar(Math.max(length, 1e-12)));
  scratch.middle.addVectors(from, to).multiplyScalar(0.5);

  return new THREE.Matrix4().compose(scratch.middle, scratch.turn, scratch.size.set(radius, length, radius));
}

function sitting(centre, radius) {
  return new THREE.Matrix4().compose(centre, UNTURNED, scratch.size.setScalar(radius));
}

// An ellipsoid of revolution about the heading: `along` is the semi-axis on
// it, `across` the waist.
function holding(centre, heading, across, along) {
  scratch.turn.setFromUnitVectors(UP, heading);

  return new THREE.Matrix4().compose(centre, scratch.turn, scratch.size.set(across, along, across));
}

// A unit spike stands along +Y with its point at the top, so a tip cap sits
// base-down on the node it grows from.
function pointing(base, heading, radius, length) {
  scratch.turn.setFromUnitVectors(UP, heading);
  scratch.middle.copy(base).addScaledVector(heading, length / 2);

  return new THREE.Matrix4().compose(scratch.middle, scratch.turn, scratch.size.set(radius, length, radius));
}

// The same spike inverted: the point is pinned to the node and the mouth
// opens outwards, which is the way a spray of leaves actually widens.
function opening(apex, heading, radius, length) {
  scratch.turn.setFromUnitVectors(UP, scratch.axis.copy(heading).negate());
  scratch.middle.copy(apex).addScaledVector(heading, length / 2);

  return new THREE.Matrix4().compose(scratch.middle, scratch.turn, scratch.size.set(radius, length, radius));
}

function direction(from, to) {
  return new THREE.Vector3().subVectors(to, from).normalize();
}

function placeOf(centre, scale) {
  return new THREE.Vector3(centre[0] * scale, centre[1] * scale, centre[2] * scale);
}

// A unit stem stands along +Y with the parent's end down, so its
// cross-section is 1 at the top and `ratio` at the bottom. Closed, so a
// flat-capped tip is a face rather than a hole; every other end's face
// hides inside a joint.
function stem(ratio) {
  return new THREE.CylinderGeometry(1, ratio, 1, AROUND, 1);
}

function spike() {
  return new THREE.ConeGeometry(1, 1, AROUND);
}

function ball(count) {
  return new THREE.IcosahedronGeometry(1, ballDetail(count));
}

function budShape(pruned, appearance) {
  return appearance.bud === 'ball' ? ball(pruned.buds.length) : spike();
}

// The screen cannot show a branch thinner than a pixel: below that a twig
// dissolves under antialiasing, and everything it carried appears to
// float. `thinnest branch` is the drawable floor, in scene units where
// the whole structure spans two.
function floored(radius, appearance) {
  const thinnest = appearance.floor ?? 0.001;

  for (let index = 0; index < radius.length; index++) {
    radius[index] = Math.max(radius[index], thinnest);
  }

  return radius;
}

// A tip is exactly `tip` wide and the root exactly `root` wide, both in
// particle radii, however large the aggregate grows; `taper` only bends
// the path between those two pinned ends, read off each branch's
// descendant count -- above one the tree stays slender and swells late,
// below one it thickens straight away.
function radii(carried, appearance, particleRadius) {
  const total = Math.max(carried.reduce((most, mass) => Math.max(most, mass), 1), 2);
  const across = Math.log(total);
  const tip = particleRadius * (appearance.tip ?? 0.2);
  const spread = Math.log((appearance.root ?? 1.2) / (appearance.tip ?? 0.2));
  const bend = appearance.taper ?? 1;

  return Float64Array.from(carried, mass =>
    tip * Math.exp(spread * Math.pow(Math.log(mass) / across, bend))
  );
}

// A particle is always frozen after the one it stuck to, so a single pass
// backwards over the array visits every child before its parent.
export function read(particles) {
  const mass = new Float64Array(particles.length).fill(1);
  const tipDistance = new Int32Array(particles.length);

  for (let index = particles.length - 1; index >= 0; index--) {
    const parent = particles[index].parent;

    if (parent == null) continue;

    mass[parent] += mass[index];
    tipDistance[parent] = Math.max(tipDistance[parent], tipDistance[index] + 1);
  }

  return {
    particles,
    mass,
    tipDistance,
    height: tipDistance.reduce((most, distance) => Math.max(most, distance), 0),
  };
}

// Every node ends in one of three fates -- staying in the trunk, gathered
// into a bud, or culled outright -- chosen by a pipeline
// of small rules over the node's own attributes. A wish rule asks who
// leaves the trunk; connectivity then makes any wish honest, a node
// leaving only when everything below it leaves, so what is drawn is always
// one tree and every pruned particle hangs below exactly one drawn node;
// and a fate rule says what each pruned cluster becomes.
export function prune(tree, appearance) {
  const { particles } = tree;
  const program = programmed(tree, appearance);
  const pruned = normalised(particles, program ? program.wish : wishes(tree, appearance));

  const shown = [];
  const edges = [];
  const carries = new Uint8Array(particles.length);

  for (let index = 0; index < particles.length; index++) {
    if (pruned[index]) continue;

    shown.push(index);
    const parent = particles[index].parent;

    if (parent == null) continue;

    edges.push({ parent, child: index });
    carries[parent] = 1;
  }

  const gathered = accumulate(particles, pruned);

  const clusters = shown.filter(node => gathered.count[node] > 0).map(node => ({
    node,
    count: gathered.count[node],
    centroid: [0, 1, 2].map(axis => gathered.sum[node * 3 + axis] / gathered.count[node]),
  }));

  const keep = program ? program.keep(gathered.count) : keeping(appearance, tree);
  const buds = clusters.filter(keep);
  const budded = buds.reduce((total, bud) => total + bud.count, 0);

  return {
    shown,
    edges,
    tips: shown.filter(node => !carries[node]),
    buds,
    budded,
    culled: clusters.reduce((total, cluster) => total + cluster.count, 0) - budded,
    note: program ? program.note : '',
  };
}

// `by rules` hands both fate decisions to a compiled program: `pruned`
// says who leaves the trunk, and `budded`, read at each carrying node,
// says whether the cluster it carries stands as a bud. `cluster_mass` is
// the one two-stage source -- it does not exist until pruning settles,
// so it reads 0 while `pruned` is decided and holds each carrier's
// cluster size when `budded` is. A program that does not compile prunes
// nothing and says why in the readout.
function programmed(tree, appearance) {
  if (appearance.fate !== 'rules') return null;

  const count = tree.particles.length;
  const environment = sources(tree);

  let program;

  try {
    program = rules.compile(appearance.rules ?? '', Object.keys(environment));
  } catch (error) {
    return { wish: new Uint8Array(count), keep: () => () => true, note: `rules: ${error.message}` };
  }

  const first = program.run({ ...environment }, count);
  const wish = Uint8Array.from({ length: count }, (nothing, index) =>
    (first.pruned ? first.pruned[index] : 0) !== 0 ? 1 : 0
  );

  return {
    wish,
    keep: clusterMass => {
      environment.cluster_mass = Float64Array.from(clusterMass);
      const second = program.run({ ...environment }, count);

      return cluster => (second.budded ? second.budded[cluster.node] : 1) !== 0;
    },
    note: '',
  };
}

// The program's sources: every structural quantity a rule can read, each
// one pass over the tree, normalised where a natural scale exists so a
// rule keeps its meaning at any aggregate size.
function sources(tree) {
  const { particles, mass, tipDistance } = tree;
  const count = particles.length;
  const span = Math.max(1, count - 1);
  const onSpine = spine(particles, mass);

  const rootLinks = new Float64Array(count);
  const spineLinks = new Float64Array(count);

  for (let index = 0; index < count; index++) {
    const parent = particles[index].parent;

    if (parent == null) continue;

    rootLinks[index] = rootLinks[parent] + 1;
    spineLinks[index] = onSpine[index] ? 0 : spineLinks[parent] + 1;
  }

  return {
    mass: Float64Array.from(mass),
    mass_fraction: Float64Array.from(mass, carried => carried / count),

    attachment_ratio: Float64Array.from(particles, (particle, index) =>
      particle.parent == null ? 1 : mass[index] / mass[particle.parent]
    ),

    tip_links: Float64Array.from(tipDistance),
    root_links: rootLinks,
    spine_links: spineLinks,
    age: Float64Array.from(particles, (particle, index) => index / span),
    cluster_mass: new Float64Array(count),
  };
}

// The spine is the maximal-mass path from the root: at every branch point
// the heaviest child carries it on. Every other subtree is a side branch,
// and its distance to the spine runs through its ancestors.
function spine(particles, mass) {
  const heaviest = new Int32Array(particles.length).fill(-1);

  for (let index = 0; index < particles.length; index++) {
    const parent = particles[index].parent;

    if (parent == null) continue;

    if (heaviest[parent] < 0 || mass[index] > mass[heaviest[parent]]) heaviest[parent] = index;
  }

  const onSpine = new Uint8Array(particles.length);
  const root = particles.findIndex(particle => particle.parent == null);

  for (let node = root; node >= 0; node = heaviest[node]) onSpine[node] = 1;

  return onSpine;
}

// The wish gate reads one of two growth measures against a threshold: how
// many links a node stands from its farthest tip, or how many descendants
// it carries -- zero descendants is a leaf, so `bud under` gathers exactly
// the crowns lighter than it says. `sway` bends either gate by growth
// age, so the young or the old side of the aggregate is pruned harder,
// and clamping to the largest possible value keeps the root.
function wishes(tree, appearance) {
  const { tipDistance, mass, height, particles } = tree;
  const sway = appearance.sway ?? 0;
  const span = Math.max(1, particles.length - 1);
  const byMass = appearance.measure === 'mass';
  const base = byMass ? (appearance.under ?? 0) : appearance.depth;
  const most = byMass ? particles.length - 1 : height;

  return Uint8Array.from(particles, (particle, index) => {
    const gate = Math.min(base * Math.exp(sway * (index / span - 0.5)), most);
    const grown = byMass ? mass[index] - 1 : tipDistance[index];

    return grown < gate ? 1 : 0;
  });
}

// A parent always precedes its children, so one pass backwards hears from
// every child before deciding the parent.
function normalised(particles, wish) {
  const pruned = new Uint8Array(particles.length);
  const anchored = new Uint8Array(particles.length);

  for (let index = particles.length - 1; index >= 0; index--) {
    pruned[index] = wish[index] && !anchored[index] ? 1 : 0;
    const parent = particles[index].parent;

    if (parent != null && !pruned[index]) anchored[parent] = 1;
  }

  return pruned;
}

// The fate rule: which pruned clusters stand as buds. The rest are
// culled: nothing is drawn for them, and `root width` is the dial for any
// thickening they might have suggested.
function keeping(appearance, tree) {
  if (appearance.fate === 'swell') return () => false;
  if (appearance.fate === 'sift') return cluster => cluster.count >= (appearance.least ?? 20);

  if (appearance.fate === 'near') {
    const away = awayFromTrunk(tree, appearance.pivot ?? 500);
    const reach = appearance.reach ?? 4;

    return cluster => away[cluster.node] > reach;
  }

  return () => true;
}

// How many links up to the nearest node with at least `pivot` descendants
// -- a node that big is trunk, and it shows itself by what hangs below it.
// Parents always precede children, so one forward pass answers.
function awayFromTrunk(tree, pivot) {
  const { particles, mass } = tree;
  const away = new Int32Array(particles.length);

  for (let index = 0; index < particles.length; index++) {
    if (mass[index] - 1 >= pivot) {
      away[index] = 0;
      continue;
    }

    const parent = particles[index].parent;
    away[index] = parent == null ? 1 << 30 : away[parent] + 1;
  }

  return away;
}

// How many pruned particles hang below each node and where their mass sits,
// gathered by the same pass backwards: a pruned particle hands its own
// tally up, a drawn one keeps what it has been handed.
function accumulate(particles, pruned) {
  const count = new Float64Array(particles.length);
  const sum = new Float64Array(particles.length * 3);

  for (let index = particles.length - 1; index >= 0; index--) {
    const parent = particles[index].parent;

    if (!pruned[index] || parent == null) continue;

    count[parent] += count[index] + 1;

    for (let axis = 0; axis < 3; axis++) {
      sum[parent * 3 + axis] += sum[index * 3 + axis] + particles[index].center[axis];
    }
  }

  return { count, sum };
}
