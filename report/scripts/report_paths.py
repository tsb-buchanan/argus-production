#!/usr/bin/env python3
"""Find the report's data directory, wherever the caller happens to be laid out.

WHY. These scripts resolved their inputs as <script's parent's parent>/data, which is true
in the working repo, where docs/report/all_geometry_2026-09-15/ holds data, fig, sections
and scripts side by side. It is NOT true in the delivery repo, which flattens report/ and
lifts data/ to the repo root, so `cd report && python3 scripts/make_cr_table.py` dies on a
report/data that does not exist. The documented workaround is to assemble a build directory
with the working-repo layout first, which works, but it makes every command in the guides
conditional on a setup step a reader can forget.

Project decision, 2026-09-24: make the scripts find data/ themselves.

HOW. Walk up from this file looking for a directory that contains data/rans_forces.json,
then take that directory as the data root. The SENTINEL IS A FILE, NOT THE data/ DIRECTORY,
because an empty or half-copied data/ would otherwise satisfy the search and the caller
would fail later with a confusing missing-key error instead of a clear missing-root one.

IT RAISES RATHER THAN GUESSING. If no ancestor carries the sentinel, that is a genuinely
unresolvable layout and the caller should stop with a message naming what it looked for and
where. Returning a plausible default here is how a script ends up reading the wrong repo's
data and reporting numbers from it.
"""
from pathlib import Path

SENTINEL = "data/rans_forces.json"


def data_root(start=None, sentinel=SENTINEL):
    """The directory whose `data/` holds the report's force file.

    Search order is nearest-first from `start` (default: this file), so a build directory
    that carries its own data/ wins over an outer repo that also has one. That matters when
    a build directory is assembled INSIDE the repo: the nearer copy is the one the caller
    staged and therefore the one it means.
    """
    here = Path(start or __file__).resolve()
    for d in [here, *here.parents]:
        if (d / sentinel).is_file():
            return d
    raise SystemExit(
        "HALT: cannot locate the report data directory. Looked for %r in %s and every "
        "parent. Run this from a tree that carries the report's data/, or assemble a build "
        "directory as the postprocessing guide describes." % (sentinel, here.parent))


def data(*parts, start=None):
    """Path under the resolved data directory: data('span_efficiency.json')."""
    return data_root(start).joinpath("data", *parts)
