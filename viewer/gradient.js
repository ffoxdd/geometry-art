import * as THREE from 'three';

// A colour map is stops on [0, 1], each a position and an sRGB hex, read
// linearly between neighbours. Three converts a hex to the linear working
// space on the way in and back on the way out, so the interpolation happens
// in linear light and the strip in the panel shows exactly what is drawn.

export const PRESETS = {
  spectral: [
    { position: 0.0, colour: '#3963d0' },
    { position: 0.1667, colour: '#39c1d0' },
    { position: 0.3333, colour: '#39d082' },
    { position: 0.5, colour: '#4ed039' },
    { position: 0.6667, colour: '#acd039' },
    { position: 0.8333, colour: '#d09739' },
    { position: 1.0, colour: '#d03939' },
  ],
  ember: [
    { position: 0.0, colour: '#1b1035' },
    { position: 0.35, colour: '#7a1f3d' },
    { position: 0.65, colour: '#d1462f' },
    { position: 0.85, colour: '#f5a623' },
    { position: 1.0, colour: '#fff2c4' },
  ],
  ice: [
    { position: 0.0, colour: '#04131f' },
    { position: 0.4, colour: '#135e83' },
    { position: 0.75, colour: '#4fb8c4' },
    { position: 1.0, colour: '#cdeef5' },
  ],
  bone: [
    { position: 0.0, colour: '#1a1712' },
    { position: 1.0, colour: '#d9d2c5' },
  ],
};

export const DEFAULT = 'spectral';

export function sampler(stops) {
  const ordered = sorted(stops);
  const colours = ordered.map(stop => new THREE.Color(stop.colour));

  return value => between(ordered, colours, Math.min(Math.max(value, 0), 1));
}

export function sorted(stops) {
  return [...stops].sort((one, other) => one.position - other.position);
}

// Every reading of the view is addressable, this one included: stops travel
// as `position-rrggbb` pairs.
export function format(stops) {
  return sorted(stops)
    .map(stop => `${Number(stop.position.toFixed(4))}-${stop.colour.slice(1)}`)
    .join(',');
}

export function parse(text) {
  const stops = (text || '').split(',').map(readStop).filter(Boolean);

  return stops.length >= 2 ? stops : null;
}

// The strip is painted rather than handed to a CSS gradient, which would
// interpolate in its own space and disagree with the drawing.
export function paint(canvas, stops) {
  const context = canvas.getContext('2d');
  const sample = sampler(stops);
  const width = canvas.width;

  for (let x = 0; x < width; x++) {
    context.fillStyle = sample(x / (width - 1)).getStyle();
    context.fillRect(x, 0, 1, canvas.height);
  }
}

export function named(stops) {
  const written = format(stops);

  return Object.keys(PRESETS).find(name => format(PRESETS[name]) === written) || null;
}

function between(ordered, colours, value) {
  const upper = ordered.findIndex(stop => stop.position >= value);

  if (upper === -1) return colours.at(-1).clone();
  if (upper === 0) return colours[0].clone();

  const low = ordered[upper - 1];
  const high = ordered[upper];
  const span = high.position - low.position;

  return new THREE.Color().lerpColors(
    colours[upper - 1],
    colours[upper],
    span > 0 ? (value - low.position) / span : 0
  );
}

function readStop(text) {
  const [position, colour] = text.split('-');

  if (!/^[0-9a-f]{6}$/i.test(colour || '')) return null;

  const at = Number(position);

  return Number.isFinite(at) ? { position: Math.min(Math.max(at, 0), 1), colour: `#${colour.toLowerCase()}` } : null;
}
