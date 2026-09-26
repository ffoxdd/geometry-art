import * as drawing from './drawing.js';
import * as controls from './controls.js';
import * as gradient from './gradient.js';
import * as editor from './gradient_editor.js';
import * as address from './address.js';

const appearance = document.getElementById('appearance');
const statistics = document.getElementById('statistics');
const rampSection = document.getElementById('ramp');
const readout = document.getElementById('readout');
const note = document.getElementById('note');
const drop = document.getElementById('drop');

const RECORD_DELAY = 400;

const remembered = address.all();
const stops = gradient.parse(remembered.ramp) || gradient.PRESETS[gradient.DEFAULT];

// Each kind of snapshot brings its own controls, and keeping the panel it
// built means a reading survives leaving that kind and coming back.
const panels = new Map();

let snapshot = null;
let live = false;
let ramp = gradient.sampler(stops);
let cut = null;

// A snapshot arriving from a run in flight is drawn cheaply; one the user
// is looking at is drawn properly.
export function show(next, options = {}) {
  snapshot = next;
  live = options.live === true;
  render();
}

export function load(text) {
  show(JSON.parse(text));
}

// Draws the snapshot already on screen again, which is how a run that has
// just settled trades its cheap drawing for the full one.
export function redraw(options = {}) {
  if (options.live !== undefined) live = options.live;
  render();
}

// The cut an export would make, drawn over the snapshot while it is set.
export function frame(next) {
  if (JSON.stringify(next) === JSON.stringify(cut)) return;

  cut = next;
  queueRender();
}

// What the last drawing cost, which paces how often a live run is polled.
export function cost() {
  return drawing.cost();
}

function render() {
  if (snapshot === null) return;

  const panel = panelFor(snapshot);
  const chosen = panel.values();

  appearance.hidden = false;
  statistics.hidden = false;
  rampSection.hidden = !Object.values(chosen).some(value => drawing.ramped(snapshot).includes(value));

  const drawn = drawing.draw(snapshot, { ...chosen, ramp, live, cut });

  renderReadout(drawn.readout);
  note.textContent = drawn.note;
}

function panelFor(next) {
  const kind = drawing.panelOf(next);

  if (!panels.has(kind)) {
    const declaration = drawing.controls(next);

    panels.set(kind, controls.create({
      parameters: declaration,
      initial: remembered,
      onChange: () => {
        queueRender();
        recordLater(declaration, panels.get(kind));
      },
    }));
  }

  const panel = panels.get(kind);

  if (appearance.firstChild !== panel.element) {
    appearance.replaceChildren(panel.element);
  }

  return panel;
}

function renderReadout(rows) {
  readout.replaceChildren(...rows.flatMap(([label, value]) => {
    const term = document.createElement('dt');
    term.textContent = label;

    const detail = document.createElement('dd');
    detail.textContent = value;

    return [term, detail];
  }));
}

// Slider drags request redraws faster than they are worth doing; coalesce
// to one per frame.
let queued = false;

function queueRender() {
  if (queued) return;
  queued = true;

  requestAnimationFrame(() => {
    queued = false;
    render();
  });
}

// The address is rewritten once a drag has stopped rather than at every
// step of it: a browser throttles history writes, and a half-finished
// setting is not a reading worth linking.
let recording = null;

function recordLater(declaration, panel) {
  clearTimeout(recording);
  recording = setTimeout(() => record(declaration, panel), RECORD_DELAY);
}

function record(declaration, panel) {
  for (const [name, value] of Object.entries(panel.values())) {
    const settled = value === null || String(value) === String(declaration[name].default);
    address.set(name, settled ? null : String(value));
  }
}

editor.use(stops, (edited, options = {}) => {
  ramp = gradient.sampler(edited);
  queueRender();

  if (options.settled) {
    address.set('ramp', gradient.format(edited));
  }
});

drop.addEventListener('click', () => {
  const picker = document.createElement('input');
  picker.type = 'file';
  picker.accept = '.json';
  picker.onchange = () => picker.files[0].text().then(load);
  picker.click();
});

document.addEventListener('dragover', event => {
  event.preventDefault();
  drop.classList.add('over');
});

document.addEventListener('dragleave', () => drop.classList.remove('over'));

document.addEventListener('drop', event => {
  event.preventDefault();
  drop.classList.remove('over');
  event.dataTransfer.files[0].text().then(load);
});
