import * as address from './address.js';
import * as form from './form.js';
import * as view from './view.js';
import * as exporter from './export.js';

const list = document.getElementById('runs');
const empty = document.getElementById('runsEmpty');
const logPanel = document.getElementById('logPanel');
const log = document.getElementById('log');
const failure = document.getElementById('error');
const runButton = document.getElementById('run');

let programs = {};
let runs = [];
let selected = null;
let version = null;
let live = false;
let drawing = false;

export async function start(declared, requested) {
  programs = declared;
  selected = requested;
  runButton.addEventListener('click', startRun);

  await refresh();
  if (selected) poll();

  tick();
}

export function describeButton(program) {
  runButton.textContent = `Run ${program.label}`;
}

async function startRun() {
  failure.textContent = '';

  try {
    const run = await api('/api/runs', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(form.values()),
    });

    await refresh();
    select(run.id);
  } catch (error) {
    failure.textContent = error.message;
  }
}

function select(id) {
  selected = id;
  version = null;
  address.set('run', id);
  render();
  poll();
}

// Polling costs a fraction of what drawing costs, so a run whose structure
// is expensive to draw simply updates less often.
async function tick() {
  await poll().catch(() => {});

  if (runs.some(run => run.status === 'running' || run.status === 'queued')) {
    await refresh().catch(() => {});
  }

  setTimeout(tick, Math.min(8000, Math.max(1000, 4 * view.cost())));
}

async function poll() {
  if (!selected) return;

  const run = await api(`/api/runs/${selected}`);
  log.textContent = run.log.slice(-12).join('\n');
  logPanel.hidden = run.log.length === 0;
  exporter.offer(run);

  const running = run.status === 'running' || run.status === 'queued';
  const settled = live && !running;
  live = running;

  await redraw(run, settled);

  const known = runs.find(other => other.id === run.id);
  if (known && known.status !== run.status) refresh();
}

// A drawing can take longer than the poll interval on a large run, so the
// next snapshot is only fetched once the previous one is on screen. The
// snapshot is fetched again when the run settles, to replace the cheap
// drawing of a run in flight with the full one.
async function redraw(run, settled) {
  const changed = run.snapshot && run.snapshot_version !== version;

  if (drawing || !(changed || settled)) return;

  drawing = true;

  try {
    if (changed) {
      version = run.snapshot_version;
      const response = await fetch(`${run.snapshot}?v=${run.snapshot_version}`);
      view.show(JSON.parse(await response.text()), { live });
    } else {
      view.redraw({ live });
    }
  } finally {
    drawing = false;
  }
}

async function refresh() {
  runs = await api('/api/runs');
  render();
}

function render() {
  list.replaceChildren(...runs.map(item));
  empty.hidden = runs.length > 0;
}

function item(run) {
  const element = document.createElement('li');
  element.className = run.id === selected ? 'selected' : '';
  element.addEventListener('click', () => select(run.id));
  element.append(status(run), text(run), action(run));

  return element;
}

function status(run) {
  const dot = document.createElement('span');
  dot.className = `dot ${run.status}`;

  return dot;
}

function text(run) {
  const element = document.createElement('div');

  const title = document.createElement('div');
  title.className = 'title';
  title.textContent = describe(run);

  const detail = document.createElement('div');
  detail.className = 'detail';
  const summary = run.summary ? run.summary.replace(/^\s*Final\s*/, '') : run.status;
  detail.textContent = `${summary} ${elapsed(run)}`.trim();

  element.append(title, detail);

  return element;
}

function action(run) {
  const button = document.createElement('button');
  const active = run.status === 'running' || run.status === 'queued';
  button.type = 'button';
  button.textContent = active ? 'stop' : '×';
  button.title = active ? 'cancel this run' : 'delete this run';

  button.addEventListener('click', async event => {
    event.stopPropagation();
    await api(`/api/runs/${run.id}` + (active ? '/cancel' : ''), { method: active ? 'POST' : 'DELETE' });

    if (!active && selected === run.id) {
      selected = null;
      exporter.offer(null);
    }

    refresh();
  });

  return button;
}

// A run is named by the parameters its program says identify it. A program
// with no binary behind it is not declared, so its old runs keep their name
// and nothing else.
function describe(run) {
  const program = programs[run.program];

  if (!program) return run.program;

  const headline = program.headline
    .map(name => run.parameters[name])
    .filter(value => value != null);

  return [program.label, ...headline].join(' · ');
}

function elapsed(run) {
  if (!run.started) return '';

  const seconds = (run.finished || Date.now() / 1000) - run.started;

  return seconds < 90 ? `${seconds.toFixed(0)}s` : `${(seconds / 60).toFixed(1)}m`;
}

async function api(path, options) {
  const response = await fetch(path, options);
  const payload = await response.json();

  if (!response.ok) throw new Error(payload.error || response.statusText);

  return payload;
}
