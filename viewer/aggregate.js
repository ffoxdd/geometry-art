import * as THREE from 'three';
import { scene, group, BACKGROUND } from './scene.js';
import { ballDetail } from './detail.js';
import * as envelope from './envelope.js';
import * as finish from './finish.js';
import * as tree from './tree.js';

export const CONTROLS = {
  colour: {
    label: 'colour',
    type: 'text',
    choices: { plain: 'plain', age: 'growth age' },
    default: 'plain',
  },
  tint: {
    label: 'tint',
    type: 'colour',
    default: '#d9d2c5',
    when: { colour: ['plain'] },
  },
  finish: {
    label: 'finish',
    type: 'text',
    choices: finish.CHOICES,
    default: 'satin',
  },
  form: {
    label: 'draw as',
    type: 'text',
    wide: true,
    choices: { particles: 'particles', tree: 'tree' },
    default: 'particles',
  },
  weight: {
    label: 'particle size',
    type: 'number',
    low: 1,
    high: 30,
    default: 9,
    range: true,
    wide: true,
    when: { form: ['particles'] },
  },
  tip: {
    label: 'tip width',
    type: 'number',
    low: 0.02,
    high: 2,
    default: 0.2,
    scale: 'log',
    range: true,
    wide: true,
    when: { form: ['tree'] },
  },
  root: {
    label: 'root width',
    type: 'number',
    low: 0.1,
    high: 60,
    default: 1.2,
    scale: 'log',
    range: true,
    wide: true,
    when: { form: ['tree'] },
  },
  taper: {
    label: 'taper',
    type: 'number',
    low: 0.3,
    high: 3,
    step: 0.05,
    default: 1,
    range: true,
    wide: true,
    when: { form: ['tree'] },
  },
  floor: {
    label: 'thinnest branch',
    type: 'number',
    low: 0.0002,
    high: 0.02,
    default: 0.001,
    scale: 'log',
    range: true,
    wide: true,
    when: { form: ['tree'] },
  },
  skin: {
    label: 'skin',
    type: 'text',
    wide: true,
    choices: { parts: 'parts', rings: 'lofted rings', smooth: 'smooth envelope' },
    default: 'parts',
    when: { form: ['tree'] },
  },
  relax: {
    label: 'smoothing',
    type: 'integer',
    low: 0,
    high: 10,
    default: 2,
    range: true,
    wide: true,
    when: { skin: ['rings'] },
  },
  blend: {
    label: 'blend',
    type: 'number',
    low: 0.1,
    high: 1.5,
    step: 0.05,
    default: 0.35,
    range: true,
    wide: true,
    when: { skin: ['smooth'] },
  },
  grid: {
    label: 'skin grid',
    type: 'text',
    choices: { 96: 'coarse', 128: 'standard', 160: 'fine' },
    default: '128',
    when: { skin: ['smooth'] },
  },
  fate: {
    label: 'pruned become',
    type: 'text',
    wide: true,
    choices: {
      bud: 'buds',
      sift: 'buds, small ones culled',
      near: 'culled near the trunk',
      swell: 'all culled',
      rules: 'by rules',
    },
    default: 'bud',
    when: { form: ['tree'] },
  },
  rules: {
    label: 'rules',
    type: 'code',
    default: 'pruned = tip_links < 4\nbudded = cluster_mass >= 20',
    wide: true,
    when: { fate: ['rules'] },
  },
  measure: {
    label: 'prune by',
    type: 'text',
    wide: true,
    choices: { links: 'links from a tip', mass: 'descendant count' },
    default: 'links',
    when: { form: ['tree'], fate: ['bud', 'sift', 'near', 'swell'] },
  },
  depth: {
    label: 'bud depth',
    type: 'integer',
    low: 0,
    high: 40,
    default: 0,
    when: { form: ['tree'], measure: ['links'] },
  },
  under: {
    label: 'bud under',
    type: 'integer',
    low: 0,
    high: 1000000,
    default: 0,
    when: { form: ['tree'], measure: ['mass'] },
  },
  sway: {
    label: 'depth by age',
    type: 'number',
    low: -3,
    high: 3,
    step: 0.1,
    default: 0,
    range: true,
    wide: true,
    when: { form: ['tree'], fate: ['bud', 'sift', 'near', 'swell'] },
  },
  least: {
    label: 'bud from',
    type: 'integer',
    low: 1,
    high: 100000,
    default: 20,
    when: { fate: ['sift'] },
  },
  pivot: {
    label: 'trunk from',
    type: 'integer',
    low: 1,
    high: 1000000,
    default: 500,
    when: { fate: ['near'] },
  },
  reach: {
    label: 'within links',
    type: 'integer',
    low: 0,
    high: 100,
    default: 4,
    when: { fate: ['near'] },
  },
  cap: {
    label: 'tip cap',
    type: 'text',
    choices: { round: 'round', flat: 'flat', cone: 'cone' },
    default: 'round',
    when: { form: ['tree'] },
  },
  bud: {
    label: 'bud shape',
    type: 'text',
    choices: { spike: 'point out', cone: 'point in', ball: 'ball', none: 'none' },
    default: 'spike',
    when: { form: ['tree'], fate: ['bud', 'sift', 'near', 'rules'] },
  },
  bud_colour: {
    label: 'bud colour',
    type: 'text',
    choices: { tint: 'its own', branch: 'as branches', age: 'growth age' },
    default: 'tint',
    when: { form: ['tree'], bud: ['spike', 'cone', 'ball'] },
  },
  bud_tint: {
    label: 'bud tint',
    type: 'colour',
    default: '#c98a3f',
    when: { bud_colour: ['tint'] },
  },
  bud_finish: {
    label: 'bud finish',
    type: 'text',
    choices: finish.CHOICES,
    default: 'satin',
    when: { form: ['tree'], bud: ['spike', 'cone', 'ball'] },
  },
  bud_aspect: {
    label: 'bud aspect',
    type: 'number',
    low: 0.25,
    high: 4,
    default: 1,
    scale: 'log',
    range: true,
    wide: true,
    when: { form: ['tree'], bud: ['spike', 'cone'] },
  },
  bud_base: {
    label: 'bud base',
    type: 'number',
    low: 0,
    high: 1,
    step: 0.05,
    default: 0.35,
    range: true,
    wide: true,
    when: { form: ['tree'], bud: ['spike', 'cone'] },
  },
  bulk: {
    label: 'bud bulk',
    type: 'number',
    low: 0.05,
    high: 3,
    step: 0.05,
    default: 1,
    range: true,
    wide: true,
    when: { form: ['tree'], bud: ['spike', 'cone', 'ball'] },
  },
};

export function controls() {
  return CONTROLS;
}

export function panelOf() {
  return 'aggregate';
}

export const RAMPED = ['age'];

// Depth comes from two cues: fog fades the far side toward the background
// as the camera moves, and -- for the particle cloud, where lighting alone
// leaves a washed-out ball -- crowded particles are darkened as if occluded.
export function draw(snapshot, appearance) {
  scene.fog = scene.fog || new THREE.Fog(BACKGROUND, 3, 6);

  if (!(appearance.form === 'tree' && appearance.skin === 'smooth')) {
    envelope.drop();
  }

  const paint = painting(snapshot, appearance);

  return appearance.form === 'tree'
    ? tree.draw(snapshot, appearance, paint)
    : cloud(snapshot, appearance, occluded(snapshot, paint.branch));
}

export function clear() {
  scene.fog = null;
  envelope.drop();
}

// A bag of equal spheres: one instanced mesh draws them all, and one line
// buffer holds the parent bonds -- hidden inside the spheres until the size
// slider shrinks them.
function cloud(snapshot, appearance, shade) {
  const particles = snapshot.particles;
  const scale = 0.98 / snapshot.reach;
  const size = snapshot.particleRadius * scale * appearance.weight / 9;

  return {
    resources: [...balls(particles, scale, size, shade, appearance), ...bonds(particles, scale, shade)],
    readout: [
      ['particles', particles.length],
      ['reach', snapshot.reach.toFixed(1)],
    ],
    note: '',
  };
}

function balls(particles, scale, size, shade, appearance) {
  const geometry = new THREE.IcosahedronGeometry(size, ballDetail(particles.length));
  const material = new THREE.MeshPhysicalMaterial({ color: 0xffffff, ...finish.of(appearance.finish) });
  const spheres = new THREE.InstancedMesh(geometry, material, particles.length);
  const matrix = new THREE.Matrix4();

  particles.forEach((particle, index) => {
    matrix.setPosition(...particle.center.map(value => value * scale));
    spheres.setMatrixAt(index, matrix);
    spheres.setColorAt(index, shade(index));
  });

  group.add(spheres);

  return [geometry, material, spheres];
}

function bonds(particles, scale, shade) {
  const positions = [];
  const colours = [];

  particles.forEach((particle, index) => {
    if (particle.parent == null) return;

    const parent = particles[particle.parent];
    const colour = shade(index);

    positions.push(
      ...particle.center.map(value => value * scale),
      ...parent.center.map(value => value * scale)
    );
    colours.push(colour.r, colour.g, colour.b, colour.r, colour.g, colour.b);
  });

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  geometry.setAttribute('color', new THREE.Float32BufferAttribute(colours, 3));

  const material = new THREE.LineBasicMaterial({ vertexColors: true });
  group.add(new THREE.LineSegments(geometry, material));

  return [geometry, material];
}

// Every part of a drawing is painted from its own source: a colour picked
// for it, or the colour map read at the value that part stands for. Buds can
// also simply take what the branches took.
function painting(snapshot, appearance) {
  const branch = tone(snapshot, appearance, appearance.colour, appearance.tint);

  if (appearance.bud_colour === 'branch') {
    return { branch, bud: branch };
  }

  return { branch, bud: tone(snapshot, appearance, appearance.bud_colour, appearance.bud_tint) };
}

function tone(snapshot, appearance, source, tint) {
  const particles = snapshot.particles;

  return source === 'age'
    ? index => appearance.ramp(index / Math.max(1, particles.length - 1))
    : () => new THREE.Color(tint);
}

// Branches are lit geometry with a silhouette of their own, so they need no
// help reading as solid; a cloud of equal balls does.
function occluded(snapshot, shade) {
  const exposure = exposures(snapshot.particles, snapshot.particleRadius);

  return index => shade(index).multiplyScalar(exposure[index]);
}

let measured = { particles: null, exposure: null };

// Crowding is the expensive part of drawing a cloud and it does not change
// while one is on screen, so recolouring or resizing reuses it.
function exposures(particles, radius) {
  if (measured.particles === particles) {
    return measured.exposure;
  }

  measured = { particles, exposure: measure(particles, radius) };

  return measured.exposure;
}

// How buried each particle is, as a stand-in for ambient occlusion: branch
// tips stand nearly alone while trunk particles sit in crowds, so a
// neighbour count on a hash grid separates them for a fraction of the cost
// of real shadowing. Crowdedness is judged against the 90th percentile, so
// the shading adapts to any aggregate size.
function measure(particles, radius) {
  const neighbourhood = 5 * radius;
  const squared = neighbourhood * neighbourhood;
  const cellOf = value => Math.floor(value / neighbourhood) + 2048;
  const keyOf = (cellX, cellY, cellZ) => (cellX * 4096 + cellY) * 4096 + cellZ;

  const grid = new Map();

  particles.forEach((particle, index) => {
    const key = keyOf(cellOf(particle.center[0]), cellOf(particle.center[1]), cellOf(particle.center[2]));
    const bucket = grid.get(key);
    if (bucket) bucket.push(index); else grid.set(key, [index]);
  });

  const counts = new Float32Array(particles.length);

  particles.forEach((particle, index) => {
    const [x, y, z] = particle.center;
    let crowd = 0;

    for (let dx = -1; dx <= 1; dx++)
      for (let dy = -1; dy <= 1; dy++)
        for (let dz = -1; dz <= 1; dz++) {
          const bucket = grid.get(keyOf(cellOf(x) + dx, cellOf(y) + dy, cellOf(z) + dz));
          if (!bucket) continue;

          for (const other of bucket) {
            const center = particles[other].center;
            const offsetX = center[0] - x;
            const offsetY = center[1] - y;
            const offsetZ = center[2] - z;
            if (offsetX * offsetX + offsetY * offsetY + offsetZ * offsetZ <= squared) crowd++;
          }
        }

    counts[index] = crowd;
  });

  const crowded = Math.max(2, Float32Array.from(counts).sort()[Math.floor(counts.length * 0.9)]);

  return Float32Array.from(counts, count => 1 - 0.82 * Math.min(count / crowded, 1));
}
