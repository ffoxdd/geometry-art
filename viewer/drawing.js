import { group } from './scene.js';
import * as tessellation from './tessellation.js';
import * as aggregate from './aggregate.js';

const DRAWINGS = { tessellation, aggregate };

let kind = null;
let resources = [];
let milliseconds = 0;

// A snapshot says what it is by its shape: cells and their boundaries, or
// particles and their bonds.
export function kindOf(snapshot) {
  return Array.isArray(snapshot.particles) ? 'aggregate' : 'tessellation';
}

// What a drawing lets you set, and which of its colourings read a value off
// every cell or particle and so have a colour map to edit.
export function controls(snapshot) {
  return DRAWINGS[kindOf(snapshot)].CONTROLS;
}

export function ramped(snapshot) {
  return DRAWINGS[kindOf(snapshot)].RAMPED;
}

export function draw(snapshot, appearance) {
  const started = performance.now();
  const wanted = kindOf(snapshot);

  release(wanted);
  const drawn = DRAWINGS[wanted].draw(snapshot, appearance);

  resources = drawn.resources;
  kind = wanted;
  milliseconds = performance.now() - started;

  return drawn;
}

// What the last drawing cost, which is what a live run's poll interval is
// paced against.
export function cost() {
  return milliseconds;
}

// Geometries, materials and instanced meshes hold GPU buffers that removing
// them from the scene does not free; a drawing being replaced by another
// kind also takes its scene fixtures back down.
function release(wanted) {
  for (const resource of resources) resource.dispose();
  resources = [];
  group.clear();

  if (kind !== null && kind !== wanted) {
    DRAWINGS[kind].clear();
  }
}
