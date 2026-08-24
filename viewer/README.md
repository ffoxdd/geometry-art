# Geometry Art studio

`index.html` renders a snapshot (`generate_globe --snapshot <name>` writes
`<name>.json`) with Three.js. Open it from any static server, drop a
snapshot onto it, or address a view directly:

```
index.html?snapshot=<file>&colour=metal|capacity|area&weight=<n>
```

`server.py` serves the same page and adds a run panel. It starts
`generate_globe` on request, one run at a time, and keeps each run under
`runs/<id>/` with its parameters, log and snapshot. The snapshot is
rewritten every second while the solver runs, so the page redraws the
tessellation as it evolves.

```
python3 viewer/server.py            # http://localhost:8731/
python3 viewer/server.py --help     # --port, --binary, --runs
```

Standard library only; needs a built `generate_globe` (defaults to
`build-release/generate_globe`). A run can be cancelled while in flight
and deleted afterwards; `?run=<id>` links to one.
