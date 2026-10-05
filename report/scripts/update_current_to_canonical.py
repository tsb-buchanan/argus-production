#!/usr/bin/env python3
"""Bring data/rans_forces_current.json into line with data/rans_forces.json, row by row.

WHY. The report and every figure read rans_forces.json; rans_forces_current.json is the raw
harvest the table captions cite. They had drifted: four of the six Condition CR trims (C, F, H
and M) still held superseded values in the _current file (C and M their published-mesh trims,
F and H the C_L-corrected untrimmed legs), and no row outside the cruise trims carried a
mesh_generation label, so a reader could not tell which generation a row was. The README had
to warn readers off the file. This makes every row the current case instead.

WHAT CHANGES. For every case present in both files, the result fields below are set to the
rans_forces.json values; where a value actually changes, the previous row is kept verbatim
under superseded_entry (an existing superseded_entry is never overwritten), so nothing is
destroyed. Cases present in only one file are a FAILURE, not a skip (GEO-089). The file is
backed up first and written through a temporary name, re-read before it is promoted.
"""
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_paths import data_root  # noqa: E402

D = data_root() / "data"
FIELDS = ["Cd", "Cd_counts", "Cl", "alpha_deg", "L_over_D", "mesh_generation",
          "mesh_generation_basis", "source_case", "trim_status", "mesh_cells"]


def main():
    canon = json.loads((D / "rans_forces.json").read_text())["cases"]
    p = D / "rans_forces_current.json"
    cur = json.loads(p.read_text())
    cases = cur["cases"]
    only_c, only_k = sorted(set(canon) - set(cases)), sorted(set(cases) - set(canon))
    if only_c or only_k:
        sys.exit("REFUSED: case sets differ; only in canonical %s, only in current %s" % (only_c, only_k))
    changed, labelled = [], 0
    for n, c in canon.items():
        e = cases[n]
        before = {k: e.get(k) for k in FIELDS}
        after = {k: c.get(k) for k in FIELDS}
        result_diff = [k for k in ("Cd", "Cl", "alpha_deg") if before[k] != after[k]]
        if result_diff and "superseded_entry" not in e:
            e["superseded_entry"] = {k: v for k, v in e.items() if k != "superseded_entry"}
        if result_diff:
            changed.append((n, result_diff, before, after))
        for k in FIELDS:
            if after[k] is not None:
                e[k] = after[k]
        if e.get("mesh_generation"):
            labelled += 1
    if labelled != len(cases):
        sys.exit("REFUSED: %d of %d rows labelled after the update" % (labelled, len(cases)))
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    stamp = now.replace("-", "").replace(":", "")
    bak = p.with_name(p.name + ".bak_" + stamp)
    shutil.copy2(p, bak)
    cur["current_to_canonical_2026_10_05"] = dict(
        what="every row's result fields and mesh_generation set to the rans_forces.json row",
        rows_with_changed_results=[c[0] for c in changed], rows_labelled=labelled,
        fields=FIELDS, backup=bak.name, generator="scripts/update_current_to_canonical.py")
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(cur, indent=1))
    json.loads(tmp.read_text())
    os.replace(tmp, p)
    for n, diff, b, a in changed:
        print("  %-10s %s  Cd %s -> %s  Cl %s -> %s  alpha %s -> %s  gen %s -> %s"
              % (n, ",".join(diff), b["Cd"], a["Cd"], b["Cl"], a["Cl"], b["alpha_deg"],
                 a["alpha_deg"], b["mesh_generation"], a["mesh_generation"]))
    print("  %d rows, %d with changed results, %d labelled; backup %s"
          % (len(cases), len(changed), labelled, bak.name))
    return 0


if __name__ == "__main__":
    sys.exit(main())
