// A set of controls built from a declaration: what each one is called, how
// to present it, and the conditions under which it applies at all -- `when`
// for the control, `choice_when` for one of its choices. A control the
// current settings rule out is taken away rather than shown dead.
// Conditions read the controls resolved before them, which is the order
// they are declared in.

export function create({ parameters, groups, initial = {}, onChange = () => {} }) {
  const element = document.createElement('div');
  const fields = new Map();

  const rendered = (groups || [{ parameters: Object.keys(parameters) }])
    .map(group => buildGroup(group, parameters, fields, initial));

  for (const group of rendered) element.append(group.element);

  const refresh = () => settle(fields, rendered);

  element.addEventListener('input', () => {
    refresh();
    onChange();
  });

  refresh();

  return {
    element,
    fields,
    refresh,
    values: () => values(fields),
  };
}

function values(fields) {
  const taken = {};

  for (const [name, field] of fields) {
    if (!field.element.hidden) {
      taken[name] = read(field);
    }
  }

  return taken;
}

function buildGroup(group, parameters, fields, initial) {
  const element = group.label ? document.createElement('details') : document.createElement('div');
  const grid = document.createElement('div');
  grid.className = 'fields';

  for (const name of group.parameters) {
    const field = buildField(name, parameters[name], initial[name]);
    fields.set(name, field);
    grid.append(field.element);
  }

  if (group.label) {
    element.className = 'group';
    element.open = group.collapsed !== true;

    const heading = document.createElement('summary');
    heading.textContent = group.label;
    element.append(heading);
  }

  element.append(grid);

  return { element, names: group.parameters };
}

function buildField(name, rule, given) {
  const element = document.createElement('label');
  element.className = rule.wide ? 'field wide' : 'field';

  const caption = document.createElement('span');
  caption.textContent = rule.label;

  const control = buildControl(rule, given);
  element.append(caption, control);

  return { element, control, caption, rule };
}

function buildControl(rule, given) {
  const chosen = given ?? rule.default;

  if (rule.type === 'code') {
    const area = document.createElement('textarea');
    area.rows = rule.rows ?? 5;
    area.spellcheck = false;
    if (rule.placeholder) area.placeholder = rule.placeholder;
    if (chosen !== null && chosen !== undefined) area.value = chosen;

    return area;
  }

  if (rule.choices) {
    const select = document.createElement('select');
    select.replaceChildren(...Object.entries(rule.choices).map(([value, label]) => new Option(label, value)));
    select.value = chosen;
    return select;
  }

  const input = document.createElement('input');
  input.type = inputType(rule);

  if (rule.scale === 'log') {
    input.min = 0;
    input.max = 1;
    input.step = 0.001;
    if (chosen !== null && chosen !== undefined) input.value = position(rule, Number(chosen));

    return input;
  }

  if (rule.low !== undefined) input.min = rule.low;
  if (rule.high !== undefined) input.max = rule.high;
  input.step = rule.step ?? (rule.type === 'integer' ? 1 : 'any');
  if (rule.placeholder) input.placeholder = rule.placeholder;
  if (chosen !== null && chosen !== undefined) input.value = chosen;

  return input;
}

// A logarithmic slider holds a position and maps it to its value, so every
// step scales the value by the same factor and the low end gets as much
// travel as the high end.
function position(rule, value) {
  return Math.log(value / rule.low) / Math.log(rule.high / rule.low);
}

function valued(rule, position) {
  return Number((rule.low * Math.pow(rule.high / rule.low, position)).toPrecision(3));
}

function inputType(rule) {
  if (rule.type === 'colour') return 'color';

  return rule.range ? 'range' : 'number';
}

function settle(fields, rendered) {
  const resolved = {};

  for (const [name, field] of fields) {
    field.element.hidden = !applies(field.rule.when, resolved);
    resolved[name] = field.element.hidden ? null : reconcile(field, resolved);
  }

  for (const group of rendered) {
    group.element.hidden = group.names.every(name => fields.get(name).element.hidden);
  }
}

function reconcile(field, resolved) {
  if (field.rule.choices) {
    offerChoices(field, resolved);
  }

  // A slider says nothing about where it sits, so its own label reads out
  // the value it is set to.
  if (field.rule.range) {
    field.caption.textContent = `${field.rule.label} ${read(field)}`;
  }

  return read(field);
}

function offerChoices(field, resolved) {
  const choices = field.rule.choices;
  const available = Object.keys(choices).filter(value => applies(condition(field.rule, value), resolved));
  const listed = [...field.control.options].map(option => option.value);
  const chosen = field.control.value;

  if (listed.join() === available.join()) return;

  field.control.replaceChildren(...available.map(value => new Option(choices[value], value)));
  field.control.value = available.includes(chosen) ? chosen : available[0];
}

function read(field) {
  const value = field.control.value;

  if (value === '') return null;
  if (field.rule.scale === 'log') return valued(field.rule, Number(value));

  return field.rule.type === 'number' || field.rule.type === 'integer' ? Number(value) : value;
}

function condition(rule, value) {
  return rule.choice_when ? rule.choice_when[value] : null;
}

function applies(condition, resolved) {
  if (!condition) return true;

  return Object.entries(condition).every(([parameter, allowed]) => allowed.includes(resolved[parameter]));
}
