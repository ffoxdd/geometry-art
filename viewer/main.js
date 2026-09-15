import * as address from './address.js';
import * as view from './view.js';
import * as form from './form.js';
import * as runs from './runs.js';
import * as exporter from './export.js';

// Either rail folds away to its title, from its title.
for (const rail of document.querySelectorAll('.rail')) {
  rail.querySelector('h1, h2').addEventListener('click', () => rail.classList.toggle('folded'));
}

const requested = address.get('snapshot');

if (requested) {
  fetch(requested).then(response => response.text()).then(view.load).catch(() => {});
}

// The run panel is present only when the page is served by viewer/server.py,
// which is what knows the programs and can start them. Only that question is
// answered by failing quietly; what the panel then does is allowed to raise.
const programs = await declared('/api/programs');

if (programs) {
  document.getElementById('studio').hidden = false;
  form.build(programs, { onProgram: (name, program) => runs.describeButton(program) });
  exporter.build(await declared('/api/exports') || {});
  runs.start(programs, address.get('run'));
}

async function declared(path) {
  try {
    const response = await fetch(path);
    return response.ok ? await response.json() : null;
  } catch {
    return null;
  }
}
