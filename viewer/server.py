#!/usr/bin/env python3
"""Serves the viewer and runs generate_globe on request.

    python3 viewer/server.py [--port 8731] [--binary build-release/generate_globe] [--runs runs]

Standard library only. Runs are queued and executed one at a time, each in
its own directory under --runs with the launch parameters, the solver log,
and the snapshot the viewer loads. The snapshot is rewritten while the run
is in progress, so the viewer can redraw the tessellation as it evolves.
"""

import argparse
import atexit
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
from urllib.parse import urlparse

REPOSITORY = Path(__file__).resolve().parent.parent
VIEWER = Path(__file__).resolve().parent

PARAMETERS = {
    "points": {"flag": "--points", "type": int, "low": 2, "high": 20000, "default": 200},
    "density_field": {"flag": "--density-field", "type": str, "choices": ["constant", "linear", "quadratic", "quadratic-piecewise", "noise", "noise-smooth", "noise-fit"], "default": "noise"},
    "seed": {"flag": "--seed", "type": int, "low": 0, "high": 2**31 - 1, "default": None},
    "warm_start": {"flag": "--warm-start", "type": str, "choices": ["lloyd", "newton"], "default": "lloyd"},
    "lloyd_passes": {"flag": "--lloyd-passes", "type": int, "low": 0, "high": 1000, "default": 5},
    "newton_iterations": {"flag": "--newton-iterations", "type": int, "low": 0, "high": 10000, "default": 50},
    "inner_solver": {"flag": "--inner-solver", "type": str, "choices": ["lbfgs", "newton"], "default": "lbfgs"},
    "newton_curvature": {"flag": "--newton-curvature", "type": str, "choices": ["gauss-newton", "finite-difference"], "default": "gauss-newton"},
    "capacity_tolerance": {"flag": "--capacity-tolerance", "type": float, "low": 1e-12, "high": 1.0, "default": 1e-7},
    "max_outer_iterations": {"flag": "--max-outer-iterations", "type": int, "low": 1, "high": 1000, "default": 30},
    "max_inner_iterations": {"flag": "--max-inner-iterations", "type": int, "low": 1, "high": 100000, "default": 200},
}

SNAPSHOT_INTERVAL_SECONDS = 1.0
LOG_TAIL_LINES = 40


class Runs:
    def __init__(self, binary, directory):
        self.binary = binary
        self.directory = directory
        self.directory.mkdir(parents=True, exist_ok=True)
        self.lock = threading.Lock()
        self.queue = deque()
        self.processes = {}
        self.worker = threading.Thread(target=self.drain, daemon=True)
        self.worker.start()
        self.mark_interrupted()
        atexit.register(self.shutdown)

    def start(self, requested):
        parameters = validate(requested)
        run_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
        folder = self.directory / run_id
        folder.mkdir()
        record = {
            "id": run_id,
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

        command = [str(self.binary), "--render", "false", "--snapshot", str(folder / "globe"),
                   "--snapshot-interval", str(SNAPSHOT_INTERVAL_SECONDS)]

        for name, value in record["parameters"].items():
            if value is not None:
                command += [PARAMETERS[name]["flag"], str(value)]

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

    def describe(self, run_id, tail=LOG_TAIL_LINES):
        folder = self.directory / run_id
        record = read_json(folder / "run.json")

        if record is None:
            return None

        snapshot = folder / "globe.json"
        record["snapshot"] = f"/runs/{run_id}/globe.json" if snapshot.exists() else None
        record["snapshot_version"] = snapshot.stat().st_mtime_ns if snapshot.exists() else None
        record["log"] = tail_lines(folder / "log.txt", tail)

        return record

    def list(self):
        records = []

        for folder in sorted(self.directory.iterdir(), reverse=True):
            record = read_json(folder / "run.json")

            if record is not None:
                record["snapshot"] = f"/runs/{folder.name}/globe.json" if (folder / "globe.json").exists() else None
                records.append(record)

        return records

    # A run that was in flight when the server last stopped has no process
    # behind it any more; say so rather than leaving it "running" forever.
    def mark_interrupted(self):
        for folder in self.directory.iterdir():
            record = read_json(folder / "run.json")

            if record is not None and record["status"] in ("queued", "running", "cancelling"):
                record["status"] = "interrupted"
                write_json(folder / "run.json", record)


def validate(requested):
    parameters = {}

    for name, rule in PARAMETERS.items():
        value = requested.get(name, rule["default"])

        if value is None or value == "":
            parameters[name] = None
            continue

        value = rule["type"](value)

        if "choices" in rule and value not in rule["choices"]:
            raise ValueError(f"{name} must be one of {rule['choices']}")

        if "low" in rule and not (rule["low"] <= value <= rule["high"]):
            raise ValueError(f"{name} must be between {rule['low']} and {rule['high']}")

        parameters[name] = value

    return parameters


def final_line(path):
    lines = tail_lines(path, 200)

    for line in reversed(lines):
        if "Final" in line:
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

    def __init__(self, *arguments, **keywords):
        super().__init__(*arguments, directory=str(VIEWER), **keywords)

    def do_GET(self):
        path = urlparse(self.path).path

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
        if self.path.endswith((".html", ".js")):
            self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def log_message(self, format, *arguments):
        if not self.path.startswith("/api/"):
            super().log_message(format, *arguments)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--port", type=int, default=8731)
    parser.add_argument("--binary", type=Path, default=REPOSITORY / "build-release" / "generate_globe")
    parser.add_argument("--runs", type=Path, default=REPOSITORY / "runs")
    arguments = parser.parse_args()

    if not arguments.binary.exists():
        parser.error(f"{arguments.binary} does not exist; build it or pass --binary")

    Handler.runs = Runs(arguments.binary.resolve(), arguments.runs.resolve())
    server = http.server.ThreadingHTTPServer(("127.0.0.1", arguments.port), Handler)
    print(f"Globe Art studio at http://localhost:{arguments.port}/  (runs in {arguments.runs})")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
