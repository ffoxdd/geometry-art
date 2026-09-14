import { PRESETS, paint, sampler, sorted, named } from './gradient.js';

const strip = document.getElementById('strip');
const canvas = document.getElementById('rampCanvas');
const preset = document.getElementById('preset');
const picker = document.getElementById('picker');

const CLICK_SLOP = 3;
const LEAST_STOPS = 2;

let stops = [];
let handles = [];
let selected = 0;
let announce = () => {};

export function use(initial, onChange) {
  announce = onChange;
  offerPresets();
  replace(initial);
}

export function current() {
  return stops;
}

function replace(next) {
  stops = sorted(next);
  selected = Math.min(selected, stops.length - 1);
  render();
}

function render() {
  paint(canvas, stops);
  renderHandles();
  preset.value = named(stops) || 'custom';
}

function renderHandles() {
  handles = stops.map((stop, index) => handle(stop, index));
  strip.replaceChildren(canvas, picker, ...handles);
}

function handle(stop, index) {
  const element = document.createElement('button');
  element.type = 'button';
  element.className = index === selected ? 'stop selected' : 'stop';
  element.style.left = `${stop.position * 100}%`;
  element.style.background = stop.colour;
  element.title = 'drag to move, click to recolour, right-click to remove';

  element.addEventListener('pointerdown', event => grab(event, index));
  element.addEventListener('contextmenu', event => {
    event.preventDefault();
    remove(index);
  });

  return element;
}

// A press on the bar itself plants a stop of the colour already there and
// then behaves as though that stop had been grabbed, so one gesture adds
// and places it.
strip.addEventListener('pointerdown', event => {
  if (event.target !== canvas) return;

  const position = positionOf(event);
  const colour = `#${sampler(stops)(position).getHexString()}`;

  stops = sorted([...stops, { position, colour }]);
  grab(event, stops.findIndex(stop => stop.position === position));
  changed();
});

function grab(event, index) {
  event.preventDefault();
  selected = index;

  const origin = event.clientX;
  let moved = false;

  const move = pointer => {
    if (Math.abs(pointer.clientX - origin) > CLICK_SLOP) moved = true;
    if (!moved) return;

    stops[index] = { ...stops[index], position: positionOf(pointer) };
    render();
    announce(stops);
  };

  const release = pointer => {
    strip.removeEventListener('pointermove', move);
    strip.removeEventListener('pointerup', release);
    strip.releasePointerCapture(pointer.pointerId);

    if (moved) {
      replace(stops);
      changed();
    } else {
      recolour(index);
    }
  };

  strip.setPointerCapture(event.pointerId);
  strip.addEventListener('pointermove', move);
  strip.addEventListener('pointerup', release);
  render();
}

function recolour(index) {
  picker.style.left = `${stops[index].position * 100}%`;
  picker.value = stops[index].colour;

  picker.oninput = () => {
    stops[index] = { ...stops[index], colour: picker.value };
    render();
    announce(stops);
  };

  picker.onchange = changed;
  picker.click();
}

function remove(index) {
  if (stops.length <= LEAST_STOPS) return;

  stops = stops.filter((_, at) => at !== index);
  replace(stops);
  changed();
}

function offerPresets() {
  const custom = new Option('custom', 'custom');
  custom.disabled = true;

  preset.replaceChildren(...Object.keys(PRESETS).map(name => new Option(name, name)), custom);
  preset.addEventListener('input', () => {
    replace(PRESETS[preset.value]);
    changed();
  });
}

function positionOf(event) {
  const bounds = canvas.getBoundingClientRect();

  return Math.min(Math.max((event.clientX - bounds.left) / bounds.width, 0), 1);
}

function changed() {
  announce(stops, { settled: true });
}
