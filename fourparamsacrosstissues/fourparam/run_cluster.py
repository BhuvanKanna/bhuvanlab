#!/usr/bin/env python3
"""
run_cluster.py

Cluster driver: generate the same 100 fourparam tables as ``generate_all.py``,
but spread the 50 tissues across several machines instead of one.

It does **not** reimplement any of the science. Each remote node runs this
repo's own ``generate_fourparam.py`` unchanged (once per threshold), so the
tables are bit-identical to a local run.

Scheduling is a **dynamic work queue**: every node pulls the next unfinished
tissue as soon as it frees up. That self-balances across machines with
different per-core speeds and means a slow or dead node just does less work
instead of holding up the run.

Per tissue on a remote node:
    rsync the one .csv.gz over -> run generate_fourparam.py twice
    (raw + excluded) -> rsync the 2 CSVs back -> delete the remote .csv.gz

The local machine also works, with no copying.

Nodes are named by their ``~/.ssh/config`` alias. A node that fails its
preflight check is dropped and the rest carry on.

Usage
-----
    python3 run_cluster.py                    # all configured nodes
    python3 run_cluster.py --nodes local,megatron
    python3 run_cluster.py --dry-run          # show the plan, do nothing
"""
from __future__ import annotations

import argparse
import queue
import shlex
import subprocess
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA_DIR = REPO / "data"
OUT_DIR = REPO / "outputs"

ID_COL, NAME_COL = "Name", "Description"
THRESHOLDS = [None, -1.0]          # raw, then excluded <= -1
REMOTE_ROOT = "~/fourparam-work"   # per-node scratch (code + one tissue at a time)
REMOTE_PY = "~/fourparam-venv/bin/python"
MAX_ATTEMPTS = 3               # per-tissue retry cap across all nodes

# jobs = worker processes used for curve fitting on that node (≈ its core count)
NODES = [
    {"name": "local",    "jobs": 20},
    {"name": "megatron", "jobs": 16},
    {"name": "hotrod",   "jobs": 10},
    {"name": "kup",      "jobs": 10},
]

SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10"]

_print_lock = threading.Lock()


def log(msg: str) -> None:
    with _print_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run(cmd: list[str], timeout: int = 7200) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


# -- naming (mirrors generate_fourparam.output_csv_for) -----------------------

def out_names(tissue: str) -> list[str]:
    stem = f"v11_log2_{tissue}"
    return [f"{stem}_fourparam.csv", f"{stem}_fourparam_excluded_at_or_below_-1.csv"]


def tissue_done(tissue: str) -> bool:
    return all((OUT_DIR / n).exists() for n in out_names(tissue))


# -- preflight ---------------------------------------------------------------

def preflight(node: dict) -> bool:
    """Verify a node can actually fit curves, and stage the code onto it."""
    name = node["name"]
    if name == "local":
        return True

    run(SSH + [name, f"mkdir -p {REMOTE_ROOT}/data {REMOTE_ROOT}/outputs"], timeout=60)
    r = run(["rsync", "-az", "-e", "ssh -o BatchMode=yes",
             f"{HERE}/bhuvanfitter.py", f"{HERE}/generate_fourparam.py",
             f"{name}:{REMOTE_ROOT}/fourparam/"], timeout=300)
    if r.returncode != 0:
        log(f"  ! {name} code sync FAILED: {r.stderr.strip()}")
        return False

    # Import the real fit library (not just its deps) and do an actual fit, so a
    # node that would die on the first gene is dropped here instead of mid-run.
    check = (
        f"cd {REMOTE_ROOT}/fourparam && {REMOTE_PY} -c "
        "'import numpy as np; from bhuvanfitter import BhuvanFitter; "
        "d=np.random.default_rng(0).normal(3,1,200); "
        "r=BhuvanFitter(d,gene_name=\"probe\").fit(\"fourparam\",max_nfev=2000); "
        "print(\"ok\" if r[\"fit_success\"] else \"fitfail\")'"
    )
    r = run(SSH + [name, check], timeout=120)
    if r.returncode != 0 or "ok" not in r.stdout:
        tail = (r.stderr or r.stdout).strip().splitlines()[-1:] or ["(no output)"]
        log(f"  ! {name} preflight FAILED: {tail[0]}")
        return False
    log(f"  + {name} ready ({node['jobs']} jobs)")
    return True


# -- one tissue on one node --------------------------------------------------

def do_local(tissue: str, jobs: int) -> bool:
    src = DATA_DIR / f"v11_log2_{tissue}.csv.gz"
    for thr in THRESHOLDS:
        cmd = ["python3", str(HERE / "generate_fourparam.py"),
               "--input", str(src), "--id-col", ID_COL, "--name-col", NAME_COL,
               "--jobs", str(jobs)]
        if thr is not None:
            cmd += ["--threshold", str(thr)]
        r = run(cmd)
        if r.returncode != 0:
            log(f"  ! local {tissue} failed: {r.stderr.strip()[-300:]}")
            return False
    return True


def do_remote(node: str, tissue: str, jobs: int) -> bool:
    fn = f"v11_log2_{tissue}.csv.gz"
    remote_in = f"{REMOTE_ROOT}/data/{fn}"

    r = run(["rsync", "-az", "-e", "ssh -o BatchMode=yes",
             str(DATA_DIR / fn), f"{node}:{REMOTE_ROOT}/data/"], timeout=3600)
    if r.returncode != 0:
        log(f"  ! {node} send {tissue} failed: {r.stderr.strip()[-200:]}")
        return False

    for thr in THRESHOLDS:
        extra = "" if thr is None else f" --threshold {thr}"
        cmd = (f"{REMOTE_PY} {REMOTE_ROOT}/fourparam/generate_fourparam.py "
               f"--input {remote_in} --id-col {ID_COL} --name-col {NAME_COL} "
               f"--jobs {jobs}{extra}")
        r = run(SSH + [node, cmd])
        if r.returncode != 0:
            log(f"  ! {node} fit {tissue} failed: {(r.stderr or r.stdout).strip()[-300:]}")
            run(SSH + [node, f"rm -f {remote_in}"], timeout=60)
            return False

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ok = True
    for n in out_names(tissue):
        r = run(["rsync", "-az", "-e", "ssh -o BatchMode=yes",
                 f"{node}:{REMOTE_ROOT}/outputs/{shlex.quote(n)}", str(OUT_DIR) + "/"],
                timeout=1800)
        if r.returncode != 0:
            log(f"  ! {node} retrieve {n} failed: {r.stderr.strip()[-200:]}")
            ok = False
    # free the node's disk; outputs are already copied back
    run(SSH + [node, f"rm -f {remote_in} {REMOTE_ROOT}/outputs/*{tissue}*"], timeout=60)
    return ok


# -- worker thread -----------------------------------------------------------

def worker(node: dict, q: "queue.Queue[str]", stats: dict, total: int,
           attempts: dict) -> None:
    name, jobs = node["name"], node["jobs"]
    while True:
        try:
            tissue = q.get_nowait()
        except queue.Empty:
            return
        if tissue_done(tissue):
            log(f"{name:9s} skip {tissue} (already done)")
            q.task_done()
            continue

        t0 = time.time()
        log(f"{name:9s} START {tissue}")
        try:
            ok = do_local(tissue, jobs) if name == "local" else do_remote(name, tissue, jobs)
        except subprocess.TimeoutExpired:
            log(f"  ! {name} TIMEOUT on {tissue}")
            ok = False
        except Exception as e:                      # noqa: BLE001 - node must never kill the run
            log(f"  ! {name} error on {tissue}: {e}")
            ok = False

        dt = time.time() - t0
        with _print_lock:
            stats[name]["time"] += dt
            stats[name]["ok" if ok else "fail"] += 1
            done = sum(s["ok"] for s in stats.values())
        log(f"{name:9s} {'DONE ' if ok else 'FAIL '}{tissue} in {dt/60:.1f}m "
            f"({done}/{total} tissues)")
        if not ok:
            with _print_lock:
                attempts[tissue] = attempts.get(tissue, 1) + 1
                retry = attempts[tissue] <= MAX_ATTEMPTS
            if retry:
                q.put(tissue)      # hand it back for another node to retry
            else:
                log(f"  ! {tissue} giving up after {MAX_ATTEMPTS} attempts")
        q.task_done()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nodes", type=str, default=None,
                   help="Comma-separated node subset (default: all configured).")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    nodes = NODES
    if args.nodes:
        want = {n.strip() for n in args.nodes.split(",")}
        nodes = [n for n in NODES if n["name"] in want]

    tissues = sorted(t.name[len("v11_log2_"):-len(".csv.gz")]
                     for t in DATA_DIR.glob("v11_log2_*.csv.gz"))
    if not tissues:
        raise SystemExit(f"No tissue matrices in {DATA_DIR}")

    pending = [t for t in tissues if not tissue_done(t)]
    print(f"{len(tissues)} tissues -> {2*len(tissues)} tables; "
          f"{len(tissues)-len(pending)} already complete, {len(pending)} to do")

    if args.dry_run:
        print("Nodes:", ", ".join(f"{n['name']}({n['jobs']})" for n in nodes))
        print("Pending:", ", ".join(pending) or "(none)")
        return
    if not pending:
        print("Nothing to do.")
        return

    print("Preflight...")
    live = [n for n in nodes if preflight(n)]
    if not live:
        raise SystemExit("No usable nodes.")
    total_jobs = sum(n["jobs"] for n in live)
    print(f"Running on {len(live)} nodes / {total_jobs} cores: "
          f"{', '.join(n['name'] for n in live)}\n")

    q: "queue.Queue[str]" = queue.Queue()
    for t in pending:
        q.put(t)

    stats = {n["name"]: {"ok": 0, "fail": 0, "time": 0.0} for n in live}
    attempts: dict[str, int] = {}

    t0 = time.time()
    threads = [threading.Thread(target=worker, args=(n, q, stats, len(pending), attempts),
                                daemon=True)
               for n in live]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    elapsed = time.time() - t0
    made = len(list(OUT_DIR.glob("*.csv"))) if OUT_DIR.exists() else 0
    print("\n" + "=" * 62)
    print(f"Wall clock: {elapsed/60:.1f} min")
    for n in live:
        s = stats[n["name"]]
        print(f"  {n['name']:9s} {s['ok']:3d} tissues  {s['fail']} failures  "
              f"{s['time']/60:6.1f} min busy")
    print(f"Tables in outputs/: {made}/{2*len(tissues)}")
    missing = [t for t in tissues if not tissue_done(t)]
    if missing:
        print(f"INCOMPLETE ({len(missing)}): {', '.join(missing)}")
        print("Re-run this command to retry them.")
    else:
        print("All tables present.")


if __name__ == "__main__":
    main()
