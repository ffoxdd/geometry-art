import * as THREE from 'three';
import { MarchingCubes } from 'three/addons/objects/MarchingCubes.js';
import { group } from './scene.js';
import * as finish from './finish.js';

// A smooth envelope over a set of tapered capsules -- each a segment from
// `a` to `b` whose radius runs `ra` to `rb`, carrying a colour. Every
// capsule splats a soft kernel into a scalar field on a grid over the unit
// cube and marching cubes pulls the level set out as one skin: welded
// where parts meet, swelling where they crowd, rounded everywhere.
// `blend` is the kernel's reach past a capsule's own radius -- how far
// apart two parts can sit and still weld, and the most the skin can
// inflate where they pile up.
//
// The grid only decides the topology. Every vertex is then projected onto
// the exact isosurface and shaded by the field's own analytic gradient,
// so neither the geometry nor the lighting carries the lattice.

const COARSE = 16;

export function skin(capsules, appearance) {
  const size = Number(appearance.grid) || 128;
  const blend = appearance.blend ?? 0.35;

  const { mesh, material } = fetch(size);
  Object.assign(material, finish.of(appearance.finish));

  const thick = thickened(capsules, size);

  mesh.reset();
  const { isolation } = splat(thick, size, blend, { field: mesh.field, palette: mesh.palette });
  mesh.isolation = isolation;
  mesh.update();

  refine(mesh, thick, size, blend, isolation);

  group.add(mesh);

  return {
    resources: [],
    triangles: Math.floor(mesh.count / 3),
    truncated: mesh.count / 3 > mesh.budget,
  };
}

// The mesh keeps large buffers, so one is held and refilled rather than
// reallocated every redraw; leaving the envelope lets it go.
let kept = null;

export function drop() {
  if (kept === null) return;

  kept.mesh.geometry.dispose();
  kept.material.dispose();
  kept = null;
}

function fetch(size) {
  if (kept !== null && kept.size === size) return kept;

  drop();

  const material = new THREE.MeshPhysicalMaterial({ vertexColors: true });
  const budget = size * size * 24;
  const mesh = new MarchingCubes(size, material, false, true, budget);
  mesh.budget = budget;

  kept = { size, mesh, material };

  return kept;
}

// The field alone, pure: alongside each cell's kernel sum runs its
// kernel-weighted colour, so a vertex takes the average of what welded
// there. A lone capsule's surface sits exactly at its own radius, which is
// where the kernel passes the returned isolation.
export function splat(capsules, size, blend, out) {
  const field = out ? out.field : new Float32Array(size * size * size);
  const palette = out ? out.palette : new Float32Array(size * size * size * 3);
  const support = 1 + blend;
  const half = size / 2;

  for (const capsule of thickened(capsules, size)) {
    const ax = capsule.a.x;
    const ay = capsule.a.y;
    const az = capsule.a.z;
    const spanX = capsule.b.x - ax;
    const spanY = capsule.b.y - ay;
    const spanZ = capsule.b.z - az;
    const squared = spanX * spanX + spanY * spanY + spanZ * spanZ;
    const thickA = capsule.ra;
    const thickB = capsule.rb;
    const red = capsule.colour.r;
    const green = capsule.colour.g;
    const blue = capsule.colour.b;
    const widest = Math.max(thickA, thickB) * support;

    const lowX = Math.max(1, Math.floor((Math.min(ax, capsule.b.x) - widest) * half + half));
    const lowY = Math.max(1, Math.floor((Math.min(ay, capsule.b.y) - widest) * half + half));
    const lowZ = Math.max(1, Math.floor((Math.min(az, capsule.b.z) - widest) * half + half));
    const highX = Math.min(size - 2, Math.ceil((Math.max(ax, capsule.b.x) + widest) * half + half));
    const highY = Math.min(size - 2, Math.ceil((Math.max(ay, capsule.b.y) + widest) * half + half));
    const highZ = Math.min(size - 2, Math.ceil((Math.max(az, capsule.b.z) + widest) * half + half));

    for (let z = lowZ; z <= highZ; z++) {
      const wz = (z - half) / half;

      for (let y = lowY; y <= highY; y++) {
        const wy = (y - half) / half;
        const row = size * size * z + size * y;

        for (let x = lowX; x <= highX; x++) {
          const wx = (x - half) / half;

          const along = squared < 1e-18 ? 0 : Math.min(Math.max(
            ((wx - ax) * spanX + (wy - ay) * spanY + (wz - az) * spanZ) / squared,
            0), 1);

          const nearX = wx - (ax + along * spanX);
          const nearY = wy - (ay + along * spanY);
          const nearZ = wz - (az + along * spanZ);
          const radius = (thickA + (thickB - thickA) * along) * support;
          const spread = (nearX * nearX + nearY * nearY + nearZ * nearZ) / (radius * radius);

          if (spread >= 1) continue;

          const kernel = (1 - spread) * (1 - spread);
          const cell = row + x;

          field[cell] += kernel;
          palette[cell * 3] += red * kernel;
          palette[cell * 3 + 1] += green * kernel;
          palette[cell * 3 + 2] += blue * kernel;
        }
      }
    }
  }

  for (let cell = 0; cell < field.length; cell++) {
    if (field[cell] === 0) continue;

    palette[cell * 3] /= field[cell];
    palette[cell * 3 + 1] /= field[cell];
    palette[cell * 3 + 2] /= field[cell];
  }

  const edge = 1 / (support * support);

  return { field, palette, isolation: (1 - edge) * (1 - edge) };
}

// The grid cannot show anything thinner than a cell. A capsule below that
// is drawn at the thinnest size the grid resolves -- 0.9 of a cell covers
// a cell centre even when the axis threads a cell corner -- so a twig
// thickens to the floor rather than shattering into disconnected lumps; a
// finer grid lowers the floor.
function thickened(capsules, size) {
  const thinnest = 1.8 / size;

  return capsules.map(capsule => ({
    ...capsule,
    ra: Math.max(capsule.ra, thinnest),
    rb: Math.max(capsule.rb, thinnest),
  }));
}

// Newton steps carry each vertex from its grid-edge estimate onto the
// exact isosurface, and the field's own gradient replaces the grid's
// finite differences as the normal. Capsules are bucketed coarsely, each
// listed in every bucket its kernel can reach, so a vertex reads one
// bucket. A step is bounded to a cell, and a vertex only moves onto
// ground the field still covers, so the last evaluation always stands.
function refine(mesh, capsules, size, blend, isolation) {
  const support = 1 + blend;
  const buckets = bucketed(capsules, support);
  const positions = mesh.geometry.getAttribute('position');
  const normals = mesh.geometry.getAttribute('normal');
  const step = 2 / size;

  for (let vertex = 0; vertex < mesh.count; vertex++) {
    let wx = positions.getX(vertex);
    let wy = positions.getY(vertex);
    let wz = positions.getZ(vertex);

    if (!fieldAt(buckets, wx, wy, wz, support)) continue;

    let value = felt.value;
    let gradientX = felt.gradientX;
    let gradientY = felt.gradientY;
    let gradientZ = felt.gradientZ;

    for (let pass = 0; pass < 2; pass++) {
      const steep = gradientX * gradientX + gradientY * gradientY + gradientZ * gradientZ;

      if (steep < 1e-12) break;

      const move = (isolation - value) / steep;
      let deltaX = move * gradientX;
      let deltaY = move * gradientY;
      let deltaZ = move * gradientZ;
      const span = Math.hypot(deltaX, deltaY, deltaZ);

      if (span > step) {
        deltaX *= step / span;
        deltaY *= step / span;
        deltaZ *= step / span;
      }

      if (!fieldAt(buckets, wx + deltaX, wy + deltaY, wz + deltaZ, support)) break;

      wx += deltaX;
      wy += deltaY;
      wz += deltaZ;
      value = felt.value;
      gradientX = felt.gradientX;
      gradientY = felt.gradientY;
      gradientZ = felt.gradientZ;
    }

    const length = Math.hypot(gradientX, gradientY, gradientZ);

    if (length < 1e-9) continue;

    positions.setXYZ(vertex, wx, wy, wz);
    normals.setXYZ(vertex, -gradientX / length, -gradientY / length, -gradientZ / length);
  }

  positions.needsUpdate = true;
  normals.needsUpdate = true;
}

// The one scratch the evaluator fills, so a hot loop never allocates.
const felt = { value: 0, gradientX: 0, gradientY: 0, gradientZ: 0 };

function fieldAt(buckets, wx, wy, wz, support) {
  const listed = buckets.get(bucketOf(wx, wy, wz));

  if (!listed) return false;

  felt.value = 0;
  felt.gradientX = 0;
  felt.gradientY = 0;
  felt.gradientZ = 0;

  for (const capsule of listed) {
    const along = capsule.squared < 1e-18 ? 0 : Math.min(Math.max(
      ((wx - capsule.a.x) * capsule.spanX +
        (wy - capsule.a.y) * capsule.spanY +
        (wz - capsule.a.z) * capsule.spanZ) / capsule.squared,
      0), 1);

    const nearX = wx - (capsule.a.x + along * capsule.spanX);
    const nearY = wy - (capsule.a.y + along * capsule.spanY);
    const nearZ = wz - (capsule.a.z + along * capsule.spanZ);
    const radius = (capsule.ra + (capsule.rb - capsule.ra) * along) * support;
    const spread = (nearX * nearX + nearY * nearY + nearZ * nearZ) / (radius * radius);

    if (spread >= 1) continue;

    const fade = 1 - spread;
    const pull = 4 * fade / (radius * radius);

    felt.value += fade * fade;
    felt.gradientX -= pull * nearX;
    felt.gradientY -= pull * nearY;
    felt.gradientZ -= pull * nearZ;
  }

  return felt.value > 0;
}

function bucketed(capsules, support) {
  const buckets = new Map();

  for (const capsule of capsules) {
    const spanX = capsule.b.x - capsule.a.x;
    const spanY = capsule.b.y - capsule.a.y;
    const spanZ = capsule.b.z - capsule.a.z;
    const squared = spanX * spanX + spanY * spanY + spanZ * spanZ;
    const listed = { ...capsule, spanX, spanY, spanZ, squared };
    const widest = Math.max(capsule.ra, capsule.rb) * support;

    const low = axis => grid(Math.min(capsule.a[axis], capsule.b[axis]) - widest);
    const high = axis => grid(Math.max(capsule.a[axis], capsule.b[axis]) + widest);

    for (let z = low('z'); z <= high('z'); z++) {
      for (let y = low('y'); y <= high('y'); y++) {
        for (let x = low('x'); x <= high('x'); x++) {
          const key = (z * COARSE + y) * COARSE + x;
          const bucket = buckets.get(key);

          if (bucket) bucket.push(listed);
          else buckets.set(key, [listed]);
        }
      }
    }
  }

  return buckets;
}

function bucketOf(wx, wy, wz) {
  return (grid(wz) * COARSE + grid(wy)) * COARSE + grid(wx);
}

function grid(value) {
  return Math.max(0, Math.min(COARSE - 1, Math.floor((value + 1) * COARSE / 2)));
}
