import * as THREE from 'three';
import { mergeGeometries } from 'three/addons/utils/BufferGeometryUtils.js';
import { scene, group } from './scene.js';
import { ballDetail } from './detail.js';
import * as finish from './finish.js';

export const CONTROLS = {
  colour: {
    label: 'colour',
    type: 'text',
    choices: { plain: 'plain', capacity: 'capacity error', area: 'cell area' },
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
    wide: true,
    choices: finish.CHOICES,
    default: 'metal',
  },
  weight: {
    label: 'edge weight',
    type: 'number',
    low: 1,
    high: 30,
    default: 9,
    range: true,
    wide: true,
  },
};

export const RAMPED = ['capacity', 'area'];

// A running solver rewrites its snapshot faster than tubes can be built for
// a large tessellation, so it is drawn as bare lines until it settles.
const OUTLINE_CELLS = 400;

let snapshot = null;

export function draw(next, appearance) {
  snapshot = next;
  enterScene();

  const cells = snapshot.cells;
  const shade = shading(cells, appearance);
  const edges = uniqueEdges(cells);
  const outline = appearance.live && cells.length > OUTLINE_CELLS;

  return {
    resources: outline ? buildLines(edges, shade) : buildTubes(edges, cells, shade, appearance),
    readout: readout(cells),
    note: outline ? 'outline while running' : '',
  };
}

export function clear() {
  shell.visible = false;
}

// The shell turns with the figure, so it joins the drawing's group -- and
// rejoins it every draw, since clearing the group detaches it.
function enterScene() {
  shell.visible = true;
  scene.fog = null;
  reshapeShell();
  group.add(shell);
}

function readout(cells) {
  const target = snapshot.targetMass;
  const masses = cells.map(cell => cell.mass);
  const spread = (Math.max(...masses) - Math.min(...masses)) / target;

  return [
    ['cells', cells.length],
    ['capacity rms', snapshot.relativeRmsCapacityError.toExponential(2)],
    ['mass spread', spread.toExponential(2)],
  ];
}

// Capacity error is scaled to the errors actually present rather than to a
// fixed multiple of the target: at convergence a fixed scale reports one
// flat colour, which says the solver worked but shows nothing about where
// the residue sits.
function shading(cells, appearance) {
  const target = snapshot.targetMass;
  const worstError = Math.max(1e-300, ...cells.map(cell => Math.abs(cell.mass - target)));
  const areas = cells.map(cell => cell.area);
  const lowArea = Math.min(...areas);
  const highArea = Math.max(...areas);

  if (appearance.colour === 'capacity') {
    return cell => appearance.ramp(Math.abs(cell.mass - target) / worstError);
  }

  if (appearance.colour === 'area') {
    return cell => appearance.ramp((cell.area - lowArea) / Math.max(1e-12, highArea - lowArea));
  }

  return () => new THREE.Color(appearance.tint);
}

// Two vertices per segment, one buffer, one draw call: cheap enough to
// rebuild every time a running solver reports a new tessellation.
function buildLines(edges, shade) {
  const positions = [];
  const colours = [];

  for (const [from, to, cell] of edges) {
    const colour = shade(cell);

    for (const points of edgePaths(from, to, 4)) {
      for (let i = 0; i < points.length - 1; i++) {
        positions.push(points[i].x, points[i].y, points[i].z, points[i + 1].x, points[i + 1].y, points[i + 1].z);
        colours.push(colour.r, colour.g, colour.b, colour.r, colour.g, colour.b);
      }
    }
  }

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  geometry.setAttribute('color', new THREE.Float32BufferAttribute(colours, 3));

  const material = new THREE.LineBasicMaterial({ vertexColors: true });
  group.add(new THREE.LineSegments(geometry, material));

  return [geometry, material];
}

function buildTubes(edges, cells, shade, appearance) {
  const radius = appearance.weight / 1200;

  const material = new THREE.MeshPhysicalMaterial({ vertexColors: true, ...finish.of(appearance.finish) });

  const radial = cells.length > 400 ? 5 : 8;
  const tubes = [];

  for (const [from, to, cell] of edges) {
    const colour = shade(cell);

    for (const points of edgePaths(from, to, 6)) {
      const curve = new THREE.CatmullRomCurve3(points);
      const tube = new THREE.TubeGeometry(curve, 6, radius, radial, false);
      tube.deleteAttribute('uv');

      const count = tube.attributes.position.count;
      const colours = new Float32Array(count * 3);

      for (let i = 0; i < count; i++) {
        colours[i * 3] = colour.r;
        colours[i * 3 + 1] = colour.g;
        colours[i * 3 + 2] = colour.b;
      }

      tube.setAttribute('color', new THREE.BufferAttribute(colours, 3));
      tubes.push(tube);
    }
  }

  const merged = mergeGeometries(tubes);
  for (const tube of tubes) tube.dispose();
  group.add(new THREE.Mesh(merged, material));

  const corners = cells.reduce((total, cell) => total + cell.boundary.length, 0);
  const joint = new THREE.IcosahedronGeometry(radius, ballDetail(corners));
  const joints = new THREE.InstancedMesh(joint, material, corners);
  let index = 0;
  const matrix = new THREE.Matrix4();

  for (const cell of cells) {
    for (const point of cell.boundary) {
      matrix.setPosition(...jointPosition(point));
      joints.setMatrixAt(index++, matrix);
    }
  }

  joints.count = index;
  group.add(joints);

  // The instanced mesh owns a GPU buffer of its own, separate from the
  // geometry it draws, and it is not freed by removing it from the scene.
  return [merged, joint, material, joints];
}

// Each bisector belongs to two cells, so it would otherwise be built twice
// and the two copies would fight for the same pixels.
function uniqueEdges(cells) {
  const seen = new Set();
  const edges = [];

  for (const cell of cells) {
    const boundary = cell.boundary;

    for (let i = 0; i < boundary.length; i++) {
      const from = boundary[i];
      const to = boundary[(i + 1) % boundary.length];
      const a = from.map(value => value.toFixed(6)).join();
      const b = to.map(value => value.toFixed(6)).join();
      const id = a < b ? a + '|' + b : b + '|' + a;

      if (!seen.has(id)) {
        seen.add(id);
        edges.push([from, to, cell]);
      }
    }
  }

  return edges;
}

function edgePaths(from, to, steps) {
  if (!flatGeometry()) {
    return [greatCircle(from, to, steps)];
  }

  return flatPieces(from, to).map(([a, b]) => {
    const points = [];

    for (let i = 0; i <= steps; i++) {
      const t = i / steps;
      points.push(embed(a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])));
    }

    return points;
  });
}

function greatCircle(from, to, steps) {
  const a = new THREE.Vector3(...from).normalize();
  const b = new THREE.Vector3(...to).normalize();
  const points = [];

  for (let i = 0; i <= steps; i++) {
    points.push(new THREE.Vector3().lerpVectors(a, b, i / steps).normalize());
  }

  return points;
}

function jointPosition(point) {
  if (!flatGeometry()) {
    return point;
  }

  const wrapped = ((point[1] % snapshot.height) + snapshot.height) % snapshot.height;
  const position = embed(point[0], wrapped);

  return [position.x, position.y, position.z];
}

// The flat snapshots carry cells in their own charts over a rectangle of
// periods; the viewer rolls the rectangle into space. The cylinder's rims
// cut straight across the cells -- partial cells are part of the look --
// while the seam stays periodic, so a protruding piece reappears from the
// other rim; the torus needs no cutting at all.
function flatGeometry() {
  return snapshot && snapshot.geometry && snapshot.geometry !== 'sphere';
}

function flatScale() {
  const radius = snapshot.width / (2 * Math.PI);

  if (snapshot.geometry === 'cylinder') {
    return 0.98 / Math.hypot(radius, snapshot.height / 2);
  }

  const tube = Math.min(snapshot.height / (2 * Math.PI), 0.75 * radius);
  return 0.98 / (radius + tube);
}

function embed(x, y) {
  const scale = flatScale();
  const radius = snapshot.width / (2 * Math.PI);
  const around = 2 * Math.PI * x / snapshot.width;

  if (snapshot.geometry === 'cylinder') {
    return new THREE.Vector3(
      radius * Math.cos(around) * scale,
      (y - snapshot.height / 2) * scale,
      radius * Math.sin(around) * scale
    );
  }

  const tube = Math.min(snapshot.height / (2 * Math.PI), 0.75 * radius);
  const along = 2 * Math.PI * y / snapshot.height;
  const ring = radius + tube * Math.cos(along);

  return new THREE.Vector3(
    ring * Math.cos(around) * scale,
    tube * Math.sin(along) * scale,
    ring * Math.sin(around) * scale
  );
}

// The chart-space pieces of one edge that survive the frame: the whole
// edge on the torus; on the cylinder, each period image clipped to the
// band between the rims.
function flatPieces(from, to) {
  if (snapshot.geometry !== 'cylinder') {
    return [[from, to]];
  }

  const pieces = [];

  for (const tile of [-1, 0, 1]) {
    const a = [from[0], from[1] + tile * snapshot.height];
    const b = [to[0], to[1] + tile * snapshot.height];
    let low = 0;
    let high = 1;
    const span = b[1] - a[1];

    for (const bound of [0, snapshot.height]) {
      const fromInside = bound === 0 ? a[1] >= 0 : a[1] <= bound;
      const toInside = bound === 0 ? b[1] >= 0 : b[1] <= bound;
      if (fromInside && toInside) continue;
      if (!fromInside && !toInside) { low = 1; high = 0; break; }
      const crossing = (bound - a[1]) / span;
      if (fromInside) high = Math.min(high, crossing); else low = Math.max(low, crossing);
    }

    if (low < high) {
      pieces.push([
        [a[0] + low * (b[0] - a[0]), a[1] + low * span],
        [a[0] + high * (b[0] - a[0]), a[1] + high * span],
      ]);
    }
  }

  return pieces;
}

// The shell sits just inside the cell network so the far side reads as
// occluded rather than as a second layer of lines; it takes the shape of
// whatever surface the snapshot lives on.
const shellMaterial = new THREE.MeshStandardMaterial({ color: 0x14151a, roughness: 0.85, metalness: 0.0 });

let shell = new THREE.Mesh(new THREE.IcosahedronGeometry(0.985, 5), shellMaterial);
let shellSurface = 'sphere';
shell.visible = false;

function reshapeShell() {
  const wanted = surface();

  if (wanted === shellSurface) return;

  shell.removeFromParent();
  shell.geometry.dispose();
  shell = new THREE.Mesh(shellGeometry(snapshot.geometry || 'sphere'), shellMaterial);
  shell.visible = true;
  shellSurface = wanted;
}

// A flat shell is cut to the run's own rectangle of periods, so the shape
// alone does not say whether the current shell still fits.
function surface() {
  if (!flatGeometry()) {
    return 'sphere';
  }

  return `${snapshot.geometry}:${snapshot.width}:${snapshot.height}`;
}

function shellGeometry(shape) {
  if (shape === 'sphere') {
    return new THREE.IcosahedronGeometry(0.985, 5);
  }

  const scale = flatScale();
  const radius = snapshot.width / (2 * Math.PI);

  if (shape === 'cylinder') {
    return new THREE.CylinderGeometry(
      0.985 * radius * scale,
      0.985 * radius * scale,
      snapshot.height * scale,
      128, 1, false
    );
  }

  const tube = Math.min(snapshot.height / (2 * Math.PI), 0.75 * radius);
  const doughnut = new THREE.TorusGeometry(radius * scale, 0.985 * tube * scale, 96, 192);
  doughnut.rotateX(Math.PI / 2);

  return doughnut;
}
