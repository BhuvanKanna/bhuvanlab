# -*- coding: utf-8 -*-
"""
stamp_version.py

Stamp `index.html` with a build id and publish the same id in `version.json`.

Why this exists
---------------
GitHub Pages serves everything under `docs/` with a fixed
`Cache-Control: max-age=600` and offers no way to change it, so for ten minutes
after a deploy a returning visitor can be handed the previous build. That is
two problems wearing one coat:

* **Sub-resources** (`manifest.json`, `genes.tsv`, `genelists/`, `qc/`, `r2/`)
  are fetched by URL, so appending `?v=<build id>` is a complete fix -- a new
  build requests a URL no cache has ever seen. `localUrl()` in index.html does
  that, using the `BUILD_ID` this script writes.

* **index.html itself cannot be busted that way.** Its URL is the one people
  visit and bookmark; there is no second name for it. So instead the page asks
  whether it is current: it fetches `version.json` with `no-store` and, if the
  id there disagrees with its own compiled-in `BUILD_ID`, offers a reload.

Both halves need `BUILD_ID` and `version.json` to agree, which is exactly the
sort of thing that rots when it is maintained by hand. Hence this script: run
it after editing index.html and before committing.

    python docs/stamp_version.py            # stamp, if the content changed
    python docs/stamp_version.py --check    # exit 1 if a stamp is due (for CI)

The id is a content hash, not a timestamp, so it only moves when index.html
actually changes -- re-running this on an unchanged file rewrites nothing and a
no-op commit never appears. The hash is taken with the `BUILD_ID` line itself
blanked out, otherwise stamping the file would change the very bytes being
hashed and the value could never settle.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INDEX = HERE / "index.html"
VERSION = HERE / "version.json"

# The one line this script owns. Matched on bytes so the file's CRLF endings
# survive untouched -- reading as text and writing back would flatten all 4,900
# of them and bury the real change in the diff.
STAMP_RE = re.compile(rb'(const BUILD_ID = ")([0-9a-f]*)(";)')
DIGEST_CHARS = 12


def compute(raw: bytes) -> str:
    """Content hash of index.html with the stamp itself blanked out."""
    neutral = STAMP_RE.sub(rb'\1\3', raw)
    return hashlib.sha256(neutral).hexdigest()[:DIGEST_CHARS]


def read_index() -> bytes:
    if not INDEX.is_file():
        sys.exit(f"error: {INDEX} not found")
    raw = INDEX.read_bytes()
    if not STAMP_RE.search(raw):
        sys.exit(
            f"error: no `const BUILD_ID = \"...\";` line in {INDEX.name}.\n"
            "       The cache-busting block was removed or renamed; restore it "
            "or this script has nothing to stamp."
        )
    return raw


def current_stamp(raw: bytes) -> str:
    return STAMP_RE.search(raw).group(2).decode()


def published() -> str | None:
    if not VERSION.is_file():
        return None
    try:
        return json.loads(VERSION.read_text(encoding="utf-8")).get("build")
    except (ValueError, OSError):
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    ap.add_argument("--check", action="store_true",
                    help="report whether a stamp is due; change nothing. "
                         "Exit 1 if index.html or version.json is out of date.")
    args = ap.parse_args()

    raw = read_index()
    want = compute(raw)
    have = current_stamp(raw)
    live = published()

    if have == want and live == want:
        print(f"up to date: build {want}")
        return

    if args.check:
        if have != want:
            print(f"STALE: index.html carries {have or '(unstamped)'}, "
                  f"content hashes to {want}")
        if live != want:
            print(f"STALE: version.json carries {live or '(missing)'}, "
                  f"expected {want}")
        print("run: python docs/stamp_version.py")
        sys.exit(1)

    if have != want:
        INDEX.write_bytes(STAMP_RE.sub(
            lambda m: m.group(1) + want.encode() + m.group(3), raw))
        print(f"index.html  {have or '(unstamped)'} -> {want}")

    # Rewritten whenever it disagrees, including when only it drifted -- the
    # page compares against this file, so a wrong value here is the one that
    # actually misleads a visitor.
    if live != want:
        VERSION.write_text(json.dumps({"build": want}) + "\n", encoding="utf-8")
        print(f"version.json {live or '(missing)'} -> {want}")


if __name__ == "__main__":
    main()
