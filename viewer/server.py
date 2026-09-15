#!/usr/bin/env python3
"""Serves the studio and runs tessellate or aggregate on request.

    ./studio [--port 8731] [--binary build-release/tessellate]
             [--aggregate-binary build-release/aggregate]
             [--skeletonize-binary build-release/skeletonize] [--runs runs]

Standard library only. Runs are queued and executed one at a time, each in
its own directory under --runs with the launch parameters, the program log,
and the snapshot the viewer loads. The snapshot is rewritten while the run
is in progress, so the viewer can redraw the structure as it evolves.

Each program declares its parameters once, here, and the page builds its
form from that declaration. A parameter says which flag it becomes, how to
present it, and when it applies at all: `when` governs the parameter,
`choice_when` governs one of its choices. A condition reads parameters
resolved before it, so a parameter is declared after everything it names.

An export is a program run over a run's snapshot whose product is handed
back as a download; it declares its parameters the same way, and which
programs' runs it applies to.
"""

import argparse
import atexit
import errno
import http.server
import json
import os
import shutil
import signal
import subprocess
import threading
import time
import uuid
from collections import deque
from pathlib import Path
from urllib.parse import parse_qsl, urlparse

REPOSITORY = Path(__file__).resolve().parent.parent
VIEWER = Path(__file__).resolve().parent

SPHERE_ONLY = {"geometry": ["sphere"]}
FLAT = {"geometry": ["cylinder", "torus", "plane"]}
GRADIENT_CAPABLE = {"geometry": ["sphere", "cylinder", "plane"]}

PROGRAMS = {
    "tessellate": {
        "label": "tessellate",
        "summary": "capacity-constrained Voronoi cells",
        "headline": ["geometry", "points", "density_field"],

        "groups": [
            {
                "label": "domain",
                "parameters": ["geometry", "points", "width", "height"],
            },
            {
                "label": "density",
                "parameters": ["density_field", "contrast", "seed"],
            },
            {
                "label": "solver",
                "collapsed": True,
                "parameters": [
                    "warm_start",
                    "lloyd_passes",
                    "newton_iterations",
                    "inner_solver",
                    "newton_curvature",
                    "capacity_tolerance",
                    "max_outer_iterations",
                    "max_inner_iterations",
                ],
            },
        ],

        "parameters": {
            "geometry": {
                "flag": "--geometry",
                "label": "geometry",
                "type": "text",
                "choices": {"sphere": "sphere", "cylinder": "cylinder", "torus": "torus", "plane": "plane"},
                "default": "sphere",
            },
            "points": {
                "flag": "--points",
                "label": "points",
                "type": "integer",
                "low": 2,
                "high": 20000,
                "default": 200,
            },
            "width": {
                "flag": "--width",
                "label": "width",
                "type": "number",
                "low": 0.1,
                "high": 100.0,
                "step": 0.1,
                "default": 2.0,
                "when": FLAT,
            },
            "height": {
                "flag": "--height",
                "label": "height",
                "type": "number",
                "low": 0.1,
                "high": 100.0,
                "step": 0.1,
                "default": 1.0,
                "when": FLAT,
            },
            "density_field": {
                "flag": "--density-field",
                "label": "field",
                "type": "text",
                "choices": {
                    "noise-smooth": "noise (C1 spline)",
                    "noise": "noise (kinked)",
                    "noise-fit": "noise (global fit)",
                    "quadratic": "quadratic",
                    "quadratic-piecewise": "quadratic (tiles)",
                    "linear": "gradient",
                    "constant": "constant",
                },
                "choice_when": {
                    "noise-smooth": SPHERE_ONLY,
                    "noise-fit": SPHERE_ONLY,
                    "quadratic": SPHERE_ONLY,
                    "quadratic-piecewise": SPHERE_ONLY,
                    "linear": GRADIENT_CAPABLE,
                },
                "default": "noise-smooth",
                "wide": True,
            },
            "contrast": {
                "flag": "--contrast",
                "label": "contrast",
                "type": "number",
                "low": 1.0,
                "high": 1000.0,
                "step": 0.5,
                "default": 4.0,
                "when": {"density_field": ["linear"]},
            },
            "seed": {
                "flag": "--seed",
                "label": "seed",
                "type": "integer",
                "low": 0,
                "high": 2**31 - 1,
                "placeholder": "random",
                "default": None,
            },
            "warm_start": {
                "flag": "--warm-start",
                "label": "warm start",
                "type": "text",
                "choices": {"lloyd": "lloyd", "newton": "newton"},
                "default": "lloyd",
            },
            "lloyd_passes": {
                "flag": "--lloyd-passes",
                "label": "lloyd passes",
                "type": "integer",
                "low": 0,
                "high": 1000,
                "default": 5,
                "when": {"warm_start": ["lloyd"]},
            },
            "newton_iterations": {
                "flag": "--newton-iterations",
                "label": "newton steps",
                "type": "integer",
                "low": 0,
                "high": 10000,
                "default": 50,
                "when": {"warm_start": ["newton"]},
            },
            "inner_solver": {
                "flag": "--inner-solver",
                "label": "inner solver",
                "type": "text",
                "choices": {"lbfgs": "lbfgs", "newton": "newton"},
                "choice_when": {"lbfgs": SPHERE_ONLY},
                "default": "lbfgs",
            },
            "newton_curvature": {
                "flag": "--newton-curvature",
                "label": "curvature",
                "type": "text",
                "choices": {
                    "exact": "exact",
                    "gauss-newton": "gauss-newton",
                    "finite-difference": "finite difference",
                },
                "default": "exact",
                "when": {"inner_solver": ["newton"]},
                "wide": True,
            },
            "capacity_tolerance": {
                "flag": "--capacity-tolerance",
                "label": "tolerance",
                "type": "number",
                "low": 1e-12,
                "high": 1.0,
                "step": "any",
                "default": 1e-7,
            },
            "max_outer_iterations": {
                "flag": "--max-outer-iterations",
                "label": "max outer",
                "type": "integer",
                "low": 1,
                "high": 1000,
                "default": 30,
            },
            "max_inner_iterations": {
                "flag": "--max-inner-iterations",
                "label": "max inner",
                "type": "integer",
                "low": 1,
                "high": 100000,
                "default": 200,
            },
        },
    },

    "aggregate": {
        "label": "aggregate",
        "summary": "diffusion-limited aggregation",
        "headline": ["particles", "return_mode"],

        "groups": [
            {
                "label": "growth",
                "parameters": ["particles", "overlap", "seed"],
            },
            {
                "label": "walk",
                "collapsed": True,
                "parameters": ["spawn_margin", "return_mode", "kill_factor"],
            },
        ],

        "parameters": {
            "particles": {
                "flag": "--particles",
                "label": "particles",
                "type": "integer",
                "low": 2,
                "high": 200000,
                "default": 5000,
            },
            "overlap": {
                "flag": "--overlap",
                "label": "overlap",
                "type": "number",
                "low": 1e-4,
                "high": 0.5,
                "step": 0.001,
                "default": 0.01,
            },
            "seed": {
                "flag": "--seed",
                "label": "seed",
                "type": "integer",
                "low": 0,
                "high": 2**31 - 1,
                "placeholder": "random",
                "default": None,
            },
            "spawn_margin": {
                "flag": "--spawn-margin",
                "label": "spawn margin",
                "type": "number",
                "low": 0.1,
                "high": 100.0,
                "step": 0.1,
                "default": 2.0,
            },
            "return_mode": {
                "flag": "--return-mode",
                "label": "return mode",
                "type": "text",
                "choices": {"harmonic": "harmonic", "kill": "kill"},
                "default": "harmonic",
            },
            "kill_factor": {
                "flag": "--kill-factor",
                "label": "kill factor",
                "type": "number",
                "low": 1.0,
                "high": 1000.0,
                "step": 0.5,
                "default": 3.0,
                "when": {"return_mode": ["kill"]},
            },
        },
    },
}

EXPORTS = {
    "skeletonize": {
        "label": "model",
        "summary": "the cell edges as bars, ready to print",
        "programs": ["tessellate"],
        "parameters": {
            "format": {
                "flag": "--format",
                "label": "format",
                "type": "text",
                "choices": {"stl": "STL", "obj": "OBJ", "ply": "PLY", "off": "OFF"},
                "default": "stl",
            },
            "scale": {
                "flag": "--scale",
                "label": "scale",
                "type": "number",
                "low": 0.1,
                "high": 10000.0,
                "step": 0.5,
                "default": 50.0,
            },
            "bar_width": {
                "flag": "--bar-width",
                "label": "bar width",
                "type": "number",
                "low": 0.01,
                "high": 1000.0,
                "step": 0.1,
                "default": 1.5,
            },
            "bar_thickness": {
                "flag": "--bar-thickness",
                "label": "bar thickness",
                "type": "number",
                "low": 0.01,
                "high": 1000.0,
                "step": 0.1,
                "default": 1.5,
            },
            "resolution": {
                "flag": "--resolution",
                "label": "resolution",
                "type": "number",
                "low": 0.01,
                "high": 1000.0,
                "step": 0.1,
                "default": 1.0,
            },
        },
    },
}

CONVERTERS = {"integer": int, "number": float, "text": str}

SNAPSHOT_INTERVAL_SECONDS = 1.0
LOG_TAIL_LINES = 40


class Runs:
    def __init__(self, binaries, directory):
        self.binaries = binaries
        self.directory = directory
        self.directory.mkdir(parents=True, exist_ok=True)
        self.lock = threading.Lock()
        self.queue = deque()
        self.processes = {}
        self.worker = threading.Thread(target=self.drain, daemon=True)
        self.worker.start()
        self.adopt()
        atexit.register(self.shutdown)

    def start(self, requested):
        program, parameters = validate(requested)
        binary = self.binaries.get(program)

        if binary is None or not binary.exists():
            raise ValueError(f"{program} is not built; expected {binary}")

        run_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
        folder = self.directory / run_id
        folder.mkdir()
        record = {
            "id": run_id,
            "program": program,
            "parameters": parameters,
            "status": "queued",
            "created": time.time(),
            "started": None,
            "finished": None,
            "summary": None,
        }
        write_json(folder / "run.json", record)

        with self.lock:
            self.queue.append(run_id)

        return record

    def drain(self):
        while True:
            with self.lock:
                run_id = self.queue.popleft() if self.queue else None

            if run_id is None:
                time.sleep(0.25)
                continue

            self.execute(run_id)

    def execute(self, run_id):
        folder = self.directory / run_id
        record = read_json(folder / "run.json")

        if record is None or record["status"] != "queued":
            return

        program = record["program"]
        command = [str(self.binaries[program]), "--snapshot", str(folder / "snapshot"),
                   "--snapshot-interval", str(SNAPSHOT_INTERVAL_SECONDS)]

        for name, value in record["parameters"].items():
            if value is not None:
                command += [PROGRAMS[program]["parameters"][name]["flag"], str(value)]

        record["status"] = "running"
        record["started"] = time.time()
        record["command"] = command
        write_json(folder / "run.json", record)

        with open(folder / "log.txt", "wb") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, cwd=str(REPOSITORY))

            with self.lock:
                self.processes[run_id] = process

            code = process.wait()

            with self.lock:
                self.processes.pop(run_id, None)

        record = read_json(folder / "run.json") or record
        record["finished"] = time.time()
        record["summary"] = final_line(folder / "log.txt")

        if record["status"] == "cancelling":
            record["status"] = "cancelled"
        else:
            record["status"] = "finished" if code == 0 else "failed"
            record["exit_code"] = code

        write_json(folder / "run.json", record)

    # The solver dies with the server rather than being orphaned mid-run.
    def shutdown(self):
        with self.lock:
            self.queue.clear()
            processes = list(self.processes.values())

        for process in processes:
            process.terminate()

    def cancel(self, run_id):
        folder = self.directory / run_id
        record = read_json(folder / "run.json")

        if record is None:
            return None

        with self.lock:
            if run_id in self.queue:
                self.queue.remove(run_id)
                record["status"] = "cancelled"
                write_json(folder / "run.json", record)
                return record

            process = self.processes.get(run_id)

        if process is not None:
            record["status"] = "cancelling"
            write_json(folder / "run.json", record)
            process.send_signal(signal.SIGTERM)

        return record

    def delete(self, run_id):
        self.cancel(run_id)
        folder = self.directory / run_id

        if folder.exists():
            shutil.rmtree(folder)

    # An export runs to completion here and now, over whatever snapshot the
    # run has at this moment, and its product is named after the run.
    def export(self, run_id, name, requested):
        folder = self.directory / run_id
        record = read_json(folder / "run.json")

        if record is None:
            raise LookupError("no such run")

        declaration = EXPORTS.get(name)
        binary = self.binaries.get(name)

        if declaration is None or binary is None or not binary.exists():
            raise ValueError(f"{name} is not built; expected {binary}")

        if record["program"] not in declaration["programs"]:
            raise ValueError(f"{name} does not apply to a {record['program']} run")

        snapshot = snapshot_file(folder)

        if not snapshot.exists():
            raise ValueError("the run has no snapshot yet")

        parameters = resolve_all(declaration, requested)
        product = folder / f"{name}.{parameters['format']}"
        command = [str(binary), str(snapshot), "--output", str(product)]

        for parameter, value in parameters.items():
            if value is not None:
                command += [declaration["parameters"][parameter]["flag"], str(value)]

        completed = subprocess.run(command, capture_output=True, text=True, cwd=str(REPOSITORY))

        if completed.returncode != 0:
            raise ValueError((completed.stderr or completed.stdout).strip().splitlines()[-1])

        return product, f"{describe_run(record)}.{parameters['format']}"

    def describe(self, run_id, tail=LOG_TAIL_LINES):
        folder = self.directory / run_id
        record = read_json(folder / "run.json")

        if record is None:
            return None

        snapshot = snapshot_file(folder)
        record["snapshot"] = f"/runs/{run_id}/{snapshot.name}" if snapshot.exists() else None
        record["snapshot_version"] = snapshot.stat().st_mtime_ns if snapshot.exists() else None
        record["log"] = tail_lines(folder / "log.txt", tail)

        return record

    def list(self):
        records = []

        for folder in sorted(self.directory.iterdir(), reverse=True):
            record = read_json(folder / "run.json")

            if record is not None:
                snapshot = snapshot_file(folder)
                record["snapshot"] = f"/runs/{folder.name}/{snapshot.name}" if snapshot.exists() else None
                records.append(record)

        return records

    # The run directory outlives any one server, so what is already in it is
    # brought up to what a record holds today: a run that was in flight has
    # no process behind it any more, and a record that names no program is
    # from before there was a second one.
    def adopt(self):
        for folder in self.directory.iterdir():
            record = read_json(folder / "run.json")

            if record is None:
                continue

            record.setdefault("program", "tessellate")

            if record["status"] in ("queued", "running", "cancelling"):
                record["status"] = "interrupted"

            write_json(folder / "run.json", record)


def validate(requested):
    program = requested.get("program") or "tessellate"

    if program not in PROGRAMS:
        raise ValueError(f"program must be one of {sorted(PROGRAMS)}")

    return program, resolve_all(PROGRAMS[program], requested)


def resolve_all(declaration, requested):
    parameters = {}

    for name, rule in declaration["parameters"].items():
        parameters[name] = resolve(name, rule, requested, parameters)

    return parameters


# A file a run hands out is named the way the page names the run.
def describe_run(record):
    program = PROGRAMS.get(record["program"])

    if program is None:
        return record["program"]

    headline = [str(record["parameters"].get(name)) for name in program["headline"]
                if record["parameters"].get(name) is not None]

    return "-".join([program["label"], *headline])


# A parameter that does not apply is dropped rather than defaulted, so the
# program is never handed a flag the form did not offer.
def resolve(name, rule, requested, resolved):
    if not applies(rule.get("when"), resolved):
        return None

    given = requested.get(name)

    if given is None or given == "":
        return fallback(rule, resolved)

    return checked(name, rule, CONVERTERS[rule["type"]](given), resolved)


def fallback(rule, resolved):
    default = rule.get("default")

    if default is None or choice_applies(rule, default, resolved):
        return default

    return next((choice for choice in rule["choices"] if choice_applies(rule, choice, resolved)), None)


def checked(name, rule, value, resolved):
    if "choices" in rule and value not in rule["choices"]:
        raise ValueError(f"{name} must be one of {sorted(rule['choices'])}")

    if not choice_applies(rule, value, resolved):
        raise ValueError(f"{name} cannot be {value} alongside the other parameters")

    if "low" in rule and not (rule["low"] <= value <= rule["high"]):
        raise ValueError(f"{name} must be between {rule['low']} and {rule['high']}")

    return value


def choice_applies(rule, value, resolved):
    return applies(rule.get("choice_when", {}).get(value), resolved)


def applies(condition, resolved):
    if condition is None:
        return True

    return all(resolved.get(parameter) in allowed for parameter, allowed in condition.items())


def snapshot_file(folder):
    current = folder / "snapshot.json"

    if current.exists():
        return current

    return folder / "globe.json"


def final_line(path):
    lines = tail_lines(path, 200)

    for line in reversed(lines):
        if "Final" in line or "Grew" in line:
            return line.strip()

    return lines[-1].strip() if lines else None


def tail_lines(path, count):
    if not path.exists():
        return []

    with open(path, "rb") as handle:
        return [line.decode("utf-8", "replace").rstrip("\n") for line in deque(handle, maxlen=count)]


def read_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def write_json(path, value):
    temporary = path.with_suffix(".tmp")

    with open(temporary, "w") as handle:
        json.dump(value, handle, indent=2)

    os.replace(temporary, path)


class Handler(http.server.SimpleHTTPRequestHandler):
    runs = None
    programs = None
    exports = None

    def __init__(self, *arguments, **keywords):
        super().__init__(*arguments, directory=str(VIEWER), **keywords)

    def do_GET(self):
        path = urlparse(self.path).path

        if path == "/api/programs":
            return self.reply(200, self.programs)

        if path == "/api/exports":
            return self.reply(200, self.exports)

        if path.startswith("/api/runs/") and "/export/" in path:
            return self.send_export(path)

        if path == "/api/runs":
            return self.reply(200, self.runs.list())

        if path.startswith("/api/runs/"):
            record = self.runs.describe(path.split("/")[3])
            return self.reply(200, record) if record else self.reply(404, {"error": "no such run"})

        if path.startswith("/runs/"):
            return self.send_run_file(path)

        return super().do_GET()

    def do_POST(self):
        path = urlparse(self.path).path
        body = self.read_body()

        if path == "/api/runs":
            try:
                return self.reply(201, self.runs.start(body))
            except (ValueError, TypeError) as error:
                return self.reply(400, {"error": str(error)})

        if path.startswith("/api/runs/") and path.endswith("/cancel"):
            record = self.runs.cancel(path.split("/")[3])
            return self.reply(200, record) if record else self.reply(404, {"error": "no such run"})

        return self.reply(404, {"error": "unknown endpoint"})

    def do_DELETE(self):
        path = urlparse(self.path).path

        if path.startswith("/api/runs/"):
            self.runs.delete(path.split("/")[3])
            return self.reply(200, {"deleted": True})

        return self.reply(404, {"error": "unknown endpoint"})

    def send_export(self, path):
        parts = path.split("/")
        requested = dict(parse_qsl(urlparse(self.path).query))

        try:
            product, filename = self.runs.export(parts[3], parts[5], requested)
        except LookupError as error:
            return self.reply(404, {"error": str(error)})
        except (ValueError, TypeError, IndexError) as error:
            return self.reply(400, {"error": str(error)})

        content = product.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(content)

    def send_run_file(self, path):
        parts = path.split("/")

        if len(parts) != 4 or ".." in path:
            return self.reply(404, {"error": "not found"})

        target = self.runs.directory / parts[2] / parts[3]

        if not target.is_file():
            return self.reply(404, {"error": "not found"})

        content = target.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", "application/json" if target.suffix == ".json" else "text/plain")
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(content)

    def read_body(self):
        length = int(self.headers.get("Content-Length") or 0)

        if length == 0:
            return {}

        return json.loads(self.rfile.read(length))

    def reply(self, status, payload):
        content = json.dumps(payload).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(content)

    def end_headers(self):
        if self.path.endswith((".html", ".js", ".css")):
            self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def log_message(self, format, *arguments):
        if not self.path.startswith("/api/"):
            super().log_message(format, *arguments)


# A program with no binary behind it is dropped from the schema, so the page
# never offers a run the server would refuse.
def buildable(declared, binaries):
    return {name: declared[name] for name in declared if binaries[name].exists()}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--port", type=int, default=8731)
    parser.add_argument("--binary", type=Path, default=REPOSITORY / "build-release" / "tessellate")
    parser.add_argument("--aggregate-binary", type=Path, default=REPOSITORY / "build-release" / "aggregate")
    parser.add_argument("--skeletonize-binary", type=Path, default=REPOSITORY / "build-release" / "skeletonize")
    parser.add_argument("--runs", type=Path, default=REPOSITORY / "runs")
    arguments = parser.parse_args()

    binaries = {
        "tessellate": arguments.binary.resolve(),
        "aggregate": arguments.aggregate_binary.resolve(),
        "skeletonize": arguments.skeletonize_binary.resolve(),
    }

    programs = buildable(PROGRAMS, binaries)

    if not programs:
        parser.error("no binaries found; build one or pass --binary / --aggregate-binary")

    exports = buildable(EXPORTS, binaries)

    for name, binary in binaries.items():
        if name not in programs and name not in exports:
            print(f"note: {name} is not built ({binary}); the studio will not offer it")

    server = listen(parser, arguments.port)
    Handler.programs = programs
    Handler.exports = exports
    Handler.runs = Runs(binaries, arguments.runs.resolve())
    print(f"Geometry Art studio at http://localhost:{arguments.port}/  (runs in {arguments.runs})")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


def listen(parser, port):
    try:
        return http.server.ThreadingHTTPServer(("127.0.0.1", port), Handler)
    except OSError as error:
        if error.errno != errno.EADDRINUSE:
            raise

        parser.error(
            f"port {port} is busy; a studio is probably already running at "
            f"http://localhost:{port}/ -- reload it, or pass --port"
        )


if __name__ == "__main__":
    main()
