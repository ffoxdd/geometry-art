import * as controls from './controls.js';

const section = document.getElementById('export');
const summary = document.getElementById('exportSummary');
const root = document.getElementById('exportForm');
const button = document.getElementById('download');
const failure = document.getElementById('exportError');

let exports = {};
let chosen = null;
let panel = null;
let run = null;

// The exports are declared by the server like the programs are, so the
// form is that declaration rendered. The studio offers the first one.
export function build(declared) {
  exports = declared;
  chosen = Object.keys(exports)[0] || null;

  if (!chosen) return;

  const definition = exports[chosen];
  summary.textContent = definition.summary;
  panel = controls.create({ parameters: definition.parameters });
  root.replaceChildren(panel.element);
  button.addEventListener('click', download);
}

// An export is offered for the run on screen, once it has a snapshot to
// read, and only when the export applies to what that run's program made.
export function offer(current) {
  run = current;
  failure.textContent = '';

  const applicable = chosen && run && run.snapshot && exports[chosen].programs.includes(run.program);
  section.hidden = !applicable;
  button.textContent = applicable ? `Download ${exports[chosen].label}` : 'Download';
}

async function download() {
  failure.textContent = '';
  button.disabled = true;

  try {
    const query = new URLSearchParams(Object.entries(panel.values()).filter(([, value]) => value !== null));
    const response = await fetch(`/api/runs/${run.id}/export/${chosen}?${query}`);

    if (!response.ok) {
      const payload = await response.json();
      throw new Error(payload.error || response.statusText);
    }

    save(await response.blob(), filenameOf(response));
  } catch (error) {
    failure.textContent = error.message;
  } finally {
    button.disabled = false;
  }
}

function filenameOf(response) {
  const disposition = response.headers.get('Content-Disposition') || '';
  const match = disposition.match(/filename="([^"]+)"/);

  return match ? match[1] : `${run.id}.${panel.values().format}`;
}

function save(blob, filename) {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;
  anchor.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
