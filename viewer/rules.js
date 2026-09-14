// A tiny closed language over per-node attributes. A program is an
// ordered list of assignments; every name a statement reads is either a
// source the tree computed or a name an earlier statement made, so the
// list is the linearisation of a dependency graph. Everything is a
// number: comparisons and logic yield 0 or 1. The language is closed --
// a fixed operator set and a fixed function table, nothing else -- so a
// program arriving from the address bar can compute, but never run code.
//
//   pruned = tip_links < 4          # who leaves the trunk
//   big    = cluster_mass >= 20
//   budded = big && spine_links > 2 # which clusters stand as buds

const FUNCTIONS = {
  min: { arity: 2, apply: Math.min },
  max: { arity: 2, apply: Math.max },
  abs: { arity: 1, apply: Math.abs },
  floor: { arity: 1, apply: Math.floor },
  sqrt: { arity: 1, apply: Math.sqrt },
  exp: { arity: 1, apply: Math.exp },
  log: { arity: 1, apply: Math.log },
  pow: { arity: 2, apply: Math.pow },
  clamp: { arity: 3, apply: (x, low, high) => Math.min(Math.max(x, low), high) },
  mix: { arity: 3, apply: (from, to, amount) => from + (to - from) * amount },

  smoothstep: {
    arity: 3,
    apply: (from, to, x) => {
      const t = Math.min(Math.max((x - from) / (to - from), 0), 1);
      return t * t * (3 - 2 * t);
    },
  },
};

const BINDING = {
  '?': 1,
  '||': 2,
  '&&': 3,
  '==': 4,
  '!=': 4,
  '<': 5,
  '<=': 5,
  '>': 5,
  '>=': 5,
  '+': 6,
  '-': 6,
  '*': 7,
  '/': 7,
  '%': 7,
  '**': 9,
};

const UNARY_BINDING = 8;

export function compile(text, sources) {
  const statements = parse(tokenise(text));
  check(statements, sources);

  return {
    targets: statements.map(statement => statement.target),
    run: (environment, length) => run(statements, environment, length),
  };
}

function run(statements, environment, length) {
  for (const statement of statements) {
    const value = evaluator(statement.expression, environment);
    const column = new Float64Array(length);

    for (let index = 0; index < length; index++) {
      column[index] = value(index);
    }

    environment[statement.target] = column;
  }

  return environment;
}

function check(statements, sources) {
  const known = new Set(sources);

  for (const statement of statements) {
    if (sources.includes(statement.target)) {
      throw new Error(`'${statement.target}' is a source and cannot be assigned`);
    }

    visit(statement.expression, known);
    known.add(statement.target);
  }
}

function visit(node, known) {
  if (node.kind === 'reference' && !known.has(node.name)) {
    throw new Error(`unknown name '${node.name}'`);
  }

  if (node.kind === 'call') {
    const declared = FUNCTIONS[node.name];

    if (!declared) throw new Error(`unknown function '${node.name}'`);

    if (node.args.length !== declared.arity) {
      throw new Error(`${node.name} takes ${declared.arity} argument${declared.arity === 1 ? '' : 's'}`);
    }
  }

  for (const child of children(node)) visit(child, known);
}

function children(node) {
  if (node.kind === 'unary') return [node.operand];
  if (node.kind === 'binary') return [node.left, node.right];
  if (node.kind === 'ternary') return [node.condition, node.consequent, node.alternate];
  if (node.kind === 'call') return node.args;

  return [];
}

// Each expression compiles to a closure over the columns it reads, built
// per run because a statement's inputs may be columns an earlier
// statement just made.
function evaluator(node, environment) {
  if (node.kind === 'literal') {
    const value = node.value;
    return () => value;
  }

  if (node.kind === 'reference') {
    const column = environment[node.name];
    return index => column[index];
  }

  if (node.kind === 'unary') {
    const operand = evaluator(node.operand, environment);

    return node.operator === '-'
      ? index => -operand(index)
      : index => (operand(index) === 0 ? 1 : 0);
  }

  if (node.kind === 'ternary') {
    const condition = evaluator(node.condition, environment);
    const consequent = evaluator(node.consequent, environment);
    const alternate = evaluator(node.alternate, environment);

    return index => (condition(index) !== 0 ? consequent(index) : alternate(index));
  }

  if (node.kind === 'call') {
    return called(node, environment);
  }

  return combined(node, environment);
}

function called(node, environment) {
  const apply = FUNCTIONS[node.name].apply;
  const args = node.args.map(argument => evaluator(argument, environment));

  if (args.length === 1) {
    const [a] = args;
    return index => apply(a(index));
  }

  if (args.length === 2) {
    const [a, b] = args;
    return index => apply(a(index), b(index));
  }

  const [a, b, c] = args;
  return index => apply(a(index), b(index), c(index));
}

function combined(node, environment) {
  const left = evaluator(node.left, environment);
  const right = evaluator(node.right, environment);

  switch (node.operator) {
    case '+': return index => left(index) + right(index);
    case '-': return index => left(index) - right(index);
    case '*': return index => left(index) * right(index);
    case '/': return index => left(index) / right(index);
    case '%': return index => left(index) % right(index);
    case '**': return index => Math.pow(left(index), right(index));
    case '<': return index => (left(index) < right(index) ? 1 : 0);
    case '<=': return index => (left(index) <= right(index) ? 1 : 0);
    case '>': return index => (left(index) > right(index) ? 1 : 0);
    case '>=': return index => (left(index) >= right(index) ? 1 : 0);
    case '==': return index => (left(index) === right(index) ? 1 : 0);
    case '!=': return index => (left(index) !== right(index) ? 1 : 0);
    case '&&': return index => (left(index) !== 0 && right(index) !== 0 ? 1 : 0);
    default: return index => (left(index) !== 0 || right(index) !== 0 ? 1 : 0);
  }
}

function parse(tokens) {
  const state = { tokens, at: 0 };
  const statements = [];

  skipBreaks(state);

  while (peek(state).kind !== 'end') {
    statements.push(statement(state));

    if (peek(state).kind !== 'break' && peek(state).kind !== 'end') {
      throw failure(state, 'expected the statement to end');
    }

    skipBreaks(state);
  }

  return statements;
}

function statement(state) {
  const name = take(state);

  if (name.kind !== 'name') throw failure(state, 'expected a name to assign');
  if (take(state).text !== '=') throw failure(state, `expected '=' after '${name.text}'`);

  return { target: name.text, expression: expression(state, 0) };
}

function expression(state, floor) {
  let left = prefix(state);

  while (true) {
    const token = peek(state);
    const binding = token.kind === 'operator' ? BINDING[token.text] : undefined;

    if (binding === undefined || binding <= floor) return left;

    take(state);

    if (token.text === '?') {
      const consequent = expression(state, 0);

      if (take(state).text !== ':') throw failure(state, "expected ':'");

      left = { kind: 'ternary', condition: left, consequent, alternate: expression(state, binding - 1) };
      continue;
    }

    const right = expression(state, token.text === '**' ? binding - 1 : binding);
    left = { kind: 'binary', operator: token.text, left, right };
  }
}

function prefix(state) {
  const token = take(state);

  if (token.kind === 'number') return { kind: 'literal', value: token.value };
  if (token.text === '-') return { kind: 'unary', operator: '-', operand: expression(state, UNARY_BINDING) };
  if (token.text === '!') return { kind: 'unary', operator: '!', operand: expression(state, UNARY_BINDING) };

  if (token.text === '(') {
    const inner = expression(state, 0);

    if (take(state).text !== ')') throw failure(state, "expected ')'");

    return inner;
  }

  if (token.kind === 'name') {
    if (peek(state).text !== '(') return { kind: 'reference', name: token.text };

    return { kind: 'call', name: token.text, args: args(state) };
  }

  throw failure(state, `unexpected '${token.text ?? 'end'}'`);
}

function args(state) {
  take(state);
  const gathered = [];

  if (peek(state).text === ')') {
    take(state);
    return gathered;
  }

  while (true) {
    gathered.push(expression(state, 0));
    const token = take(state);

    if (token.text === ')') return gathered;
    if (token.text !== ',') throw failure(state, "expected ',' or ')'");
  }
}

function peek(state) {
  return state.tokens[state.at];
}

function take(state) {
  return state.tokens[state.at++];
}

function skipBreaks(state) {
  while (peek(state).kind === 'break') take(state);
}

function failure(state, message) {
  return new Error(message);
}

const TWO_CHARACTER = ['**', '<=', '>=', '==', '!=', '&&', '||'];
const ONE_CHARACTER = '+-*/%<>!(),?:=';

function tokenise(text) {
  const tokens = [];
  let at = 0;

  while (at < text.length) {
    const character = text[at];

    if (character === '\n' || character === ';') {
      tokens.push({ kind: 'break' });
      at++;
      continue;
    }

    if (/\s/.test(character)) {
      at++;
      continue;
    }

    if (character === '#') {
      while (at < text.length && text[at] !== '\n') at++;
      continue;
    }

    const two = text.slice(at, at + 2);

    if (TWO_CHARACTER.includes(two)) {
      tokens.push({ kind: 'operator', text: two });
      at += 2;
      continue;
    }

    if (ONE_CHARACTER.includes(character)) {
      tokens.push({ kind: 'operator', text: character });
      at++;
      continue;
    }

    const number = /^\d+(\.\d+)?([eE][+-]?\d+)?/.exec(text.slice(at));

    if (number) {
      tokens.push({ kind: 'number', value: Number(number[0]) });
      at += number[0].length;
      continue;
    }

    const name = /^[A-Za-z_][A-Za-z0-9_]*/.exec(text.slice(at));

    if (name) {
      tokens.push({ kind: 'name', text: name[0] });
      at += name[0].length;
      continue;
    }

    throw new Error(`unexpected character '${character}'`);
  }

  tokens.push({ kind: 'end' });

  return tokens;
}
