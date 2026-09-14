import * as controls from './controls.js';

const tabs = document.getElementById('programs');
const summary = document.getElementById('programSummary');
const root = document.getElementById('parameters');

let programs = {};
let panels = new Map();
let current = null;
let announce = () => {};

// The parameters of every program are declared by the server; the form is
// that declaration rendered, so a parameter is described in one place only.
export function build(declared, { onProgram }) {
  programs = declared;
  announce = onProgram;

  tabs.replaceChildren(...Object.entries(programs).map(([name, program]) => tab(name, program)));
  tabs.hidden = Object.keys(programs).length < 2;

  select(Object.keys(programs)[0]);
}

export function values() {
  return { program: current.name, ...current.panel.values() };
}

function tab(name, program) {
  const button = document.createElement('button');
  button.type = 'button';
  button.textContent = program.label;
  button.addEventListener('click', () => select(name));

  return button;
}

function select(name) {
  current = panels.get(name) || remember(name);
  summary.textContent = programs[name].summary;
  root.replaceChildren(current.panel.element);

  for (const [index, button] of [...tabs.children].entries()) {
    button.setAttribute('aria-pressed', Object.keys(programs)[index] === name);
  }

  current.panel.refresh();
  announce(name, programs[name]);
}

function remember(name) {
  const definition = programs[name];

  const entry = {
    name,
    panel: controls.create({
      parameters: definition.parameters,
      groups: definition.groups,
    }),
  };

  panels.set(name, entry);

  return entry;
}
