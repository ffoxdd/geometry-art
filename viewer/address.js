// The whole view is addressable, so a particular reading of a particular
// run can be linked rather than described. Every part of the page writes its
// own parameter through here, so none of them drops another's.

export function all() {
  return Object.fromEntries(new URLSearchParams(location.search));
}

export function get(name) {
  return new URLSearchParams(location.search).get(name);
}

export function set(name, value) {
  const query = new URLSearchParams(location.search);

  if (value === null) {
    query.delete(name);
  } else {
    query.set(name, value);
  }

  history.replaceState(null, '', `?${query}`);
}
