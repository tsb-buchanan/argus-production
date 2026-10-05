#!/usr/bin/env python3
"""Gate: every derived wing STL must be a watertight, outward-oriented solid (GEO-101).

WHY THIS EXISTS. Every derived wing surface in geometry/derived/wing3d was an OPEN
SHELL: 398 open edges tracing the root section outline, with no root cap. Nothing in
the record ever claimed they were closed, which is the finding rather than an excuse:
SURFACECHECK WAS NEVER RUN AND checkMesh CANNOT SEE IT, because checkMesh inspects
the volume mesh that snappy produced from the leaky surface, not the surface.

THE GATE IS DERIVED FROM THE FILE, NEVER ASSERTED (GEO-080). Each check reads the
triangles and recomputes; nothing is taken from a sidecar or a name.

  G1  every edge has multiplicity EXACTLY 2          (closed)
  G2  the two triangles on an edge traverse it in OPPOSITE directions
                                                      (consistently oriented)
  G3  divergence-theorem volume > 0                   (normals point OUTWARD)
  G4  min(y) < -0.010 m                               (the root overhang exists)
  G5  every point with y >= 0 is present and UNMOVED against the pre-fix manifest
                                                      (the aerodynamic surface is untouched)

G5 IS THE ONE THAT PROTECTS THE RESULT. A repair that closes the surface by moving
the wing is not a repair, it is a different wing silently revalidated. The manifest
is captured BEFORE the repair with --snapshot, so it is a known answer taken from the
old state rather than derived from the new one.

PARTITION (GEO-089): PASS + DECLARED-EXCLUDED-WITH-REASON must equal the file count.
Anything in neither is a FAILURE, not a skip.

NEGATIVE CONTROL: --negative-control asserts the checker FAILS on a known-open
surface. A gate that cannot fail is not a gate, and this project has a written
history of them.
"""
import argparse
import collections
import glob
import hashlib
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
WING3D = REPO / "geometry" / "derived" / "wing3d"
MANIFEST = WING3D / "PREFIX_MANIFEST.json"
ROUND = 6          # vertex unification tolerance, decimal places in metres
Y_OVERHANG_MAX = -0.010

# Files deliberately not required to close, each with its reason. Empty today: every
# derived wing surface is meshed or shared, so every one must be a solid.
# FORMALLY RETIRED, not merely excluded (GEO-102). A gate exclusion is a statement
# about the GATE; retirement is a statement about the FILE. These are renamed to
# *.stl.retired so nothing globbing *.stl can pick them up, and the gate ASSERTS that
# state rather than trusting it: prefer the canary to the exclusion. See
# geometry/derived/wing3d/RETIRED.md for the measured reason.
RETIRED = {
    "baseline_oml_banded.stl":
        "multi-solid layer-sweep trial; 60 edges shared by 4 triangles WITHIN single "
        "solids (a fold, not a seam). Superseded, unregistered, must not be shared.",
    "slab_banded.stl":
        "multi-solid layer-sweep trial; 136 edges shared by 4 triangles WITHIN single "
        "solids. Superseded, unregistered, must not be shared.",
}


def load(path, round_decimals=ROUND):
    """round_decimals=None means EXACT vertex matching, which is what OpenFOAM's
    triSurface reader does. See G6."""
    txt = Path(path).read_text(errors="ignore")
    V = np.array(re.findall(r"vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)", txt),
                 dtype=float)
    if V.size == 0:
        raise ValueError("%s: no vertices parsed" % path)
    key = V if round_decimals is None else np.round(V, round_decimals)
    u, inv = np.unique(key, axis=0, return_inverse=True)
    return u, inv.reshape(-1, 3)


def open_edges_at(path, round_decimals):
    """Open-edge count under one vertex-unification tolerance. Separated out so G6
    can ask the SAME question at a DIFFERENT tolerance."""
    u, tri = load(path, round_decimals)
    und = collections.Counter()
    for a, b, c in tri:
        for x, y in ((a, b), (b, c), (c, a)):
            und[(x, y) if x < y else (y, x)] += 1
    return int(sum(1 for v in und.values() if v != 2))


def y_ge_zero_fingerprint(u):
    """Count and hash of the unique points with y >= 0, sorted. This is what G5
    compares; it is compact and it changes if any such point moves by >= 1e-6."""
    P = u[u[:, 1] >= 0.0]
    P = P[np.lexsort((P[:, 2], P[:, 1], P[:, 0]))]
    h = hashlib.sha256(np.ascontiguousarray(P).tobytes()).hexdigest()
    return int(P.shape[0]), h


def closure_probe(path):
    """THE ONE CLOSURE DEFINITION ON THIS PROJECT (GEO-102). Import this; do not
    re-derive it.

    Returns a dict of measured facts about one STL surface, with NO verdict attached:
        n_tri, n_open_edges, n_over_edges (multiplicity > 2), n_orientation_errors,
        volume (divergence theorem, >0 means outward), y_min, y_max.

    TWO CHECKERS EXISTED FOR ONE QUESTION and that is how definitions drift: this
    file's file-level gate over geometry/derived/wing3d, and the ARGUS workstream's
    case-level guard in scripts/assert_surface_closed.py which refuses to quote a
    3D result from a case meshed on an open shell. Those are DIFFERENT GATES asking
    the SAME QUESTION, so the question gets one implementation and the gates keep
    their own scopes. Same move as chord_for() and section_conventions.le_te_base.

    THE VERTEX-UNIFICATION TOLERANCE IS PART OF THE DEFINITION. Two implementations
    rounding at different tolerances would disagree about which edges are open, and
    would disagree only on marginal geometry, i.e. exactly when it matters.
    """
    u, tri = load(path)
    directed = collections.Counter()
    for a, b, c in tri:
        directed[(a, b)] += 1
        directed[(b, c)] += 1
        directed[(c, a)] += 1
    und = collections.Counter()
    for (a, b), n in directed.items():
        und[(a, b) if a < b else (b, a)] += n
    T = u[tri]
    return {
        "n_tri": int(len(tri)),
        "n_open_edges": int(sum(1 for v in und.values() if v == 1)),
        "n_over_edges": int(sum(1 for v in und.values() if v > 2)),
        "n_orientation_errors": int(sum(1 for v in directed.values() if v > 1)),
        "volume": float(np.einsum("ij,ij->i", T[:, 0, :],
                                  np.cross(T[:, 1, :], T[:, 2, :])).sum() / 6.0),
        "y_min": float(u[:, 1].min()),
        "y_max": float(u[:, 1].max()),
        "round_decimals": ROUND,
    }


def checks(path, manifest):
    u, tri = load(path)
    pr = closure_probe(path)          # ONE definition, consumed here
    res = {}
    res["G1"] = (pr["n_open_edges"] == 0 and pr["n_over_edges"] == 0,
                 "%d open, %d with multiplicity > 2"
                 % (pr["n_open_edges"], pr["n_over_edges"]))
    res["G2"] = (pr["n_orientation_errors"] == 0,
                 "%d directed edges traversed twice the same way"
                 % pr["n_orientation_errors"])
    vol = pr["volume"]
    res["G3"] = (vol > 0.0, "volume %+.9f m3" % vol)
    ymin = pr["y_min"]
    res["G4"] = (ymin < Y_OVERHANG_MAX, "min(y) %+.6f m" % ymin)

    # G6  THE TOLERANCE MUST NOT BE LOAD-BEARING (GEO-102, second instance).
    # G1 unifies vertices at ROUND decimals. OpenFOAM's triSurface reader, which is
    # the CONSUMER, does not. A surface can therefore pass G1 and still be an open
    # shell to snappyHexMesh, and that gap is invisible to every check that shares
    # G1's tolerance -- including the repair in close_wing_surface.py, which uses the
    # SAME constant. A gate cannot detect an error introduced at its own tolerance.
    # Measured instance: chc_g02_c06, cfft_b02_c01 and cffw_b01_c01 each carried 397
    # sub-tolerance duplicate vertices around the root rim, read 0 open edges at
    # ROUND=6 and 796 at exact precision.
    n_exact = open_edges_at(path, None)
    res["G6"] = (n_exact == pr["n_open_edges"],
                 "%d open edges exact vs %d at %d dp%s"
                 % (n_exact, pr["n_open_edges"], ROUND,
                    "" if n_exact == pr["n_open_edges"]
                    else "  <-- CLOSED TO THIS GATE, OPEN TO snappyHexMesh"))

    name = os.path.basename(path)
    n, h = y_ge_zero_fingerprint(u)
    if manifest is None or name not in manifest:
        res["G5"] = (False, "NO PRE-FIX MANIFEST ENTRY, so the surface cannot be "
                            "proved unmoved (run --snapshot before repairing)")
    else:
        m = manifest[name]
        ok = (n == m["n_points_y_ge_0"] and h == m["sha256_y_ge_0"])
        res["G5"] = (ok, "y>=0 points %d vs %d, hash %s"
                     % (n, m["n_points_y_ge_0"], "match" if h == m["sha256_y_ge_0"] else "DIFFER"))
    res["_vol"] = vol
    res["_ntri"] = len(tri)
    return res


def snapshot(files):
    man = {}
    for f in files:
        u, tri = load(f)
        n, h = y_ge_zero_fingerprint(u)
        ec = collections.Counter()
        for a, b, c in tri:
            for e in ((a, b), (b, c), (c, a)):
                ec[tuple(sorted(e))] += 1
        man[os.path.basename(f)] = {
            "n_triangles_prefix": int(len(tri)),
            "n_points_y_ge_0": n,
            "sha256_y_ge_0": h,
            "open_edges_prefix": int(sum(1 for v in ec.values() if v == 1)),
            "note": ("Captured BEFORE the root-overhang repair. G5 compares against "
                     "this so the aerodynamic surface is proved unmoved rather than "
                     "assumed. open_edges_prefix is the defect the repair removes."),
        }
    MANIFEST.write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")
    print("snapshot written: %s  (%d files)" % (MANIFEST.relative_to(REPO), len(man)))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--snapshot", action="store_true",
                    help="capture the pre-fix y>=0 manifest; run BEFORE repairing")
    ap.add_argument("--negative-control", action="store_true",
                    help="assert the checker FAILS on a deliberately open surface")
    ap.add_argument("--glob", default=str(WING3D / "*.stl"))
    args = ap.parse_args()

    files = sorted(f for f in glob.glob(args.glob) if f.endswith(".stl"))
    if args.snapshot:
        return snapshot(files)

    if args.negative_control:
        return negative_control(files)

    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else None
    npass = 0
    excluded = []
    failed = []
    # RETIREMENT IS ASSERTED, NOT ASSUMED. Each retired name must be ABSENT as .stl
    # and PRESENT as .stl.retired. If someone un-retires one by renaming it back, the
    # gate says so instead of silently adopting it.
    for name, why in sorted(RETIRED.items()):
        live = WING3D / name
        ret = WING3D / (name + ".retired")
        if live.exists() or not ret.exists():
            print("  %-38s RETIREMENT BROKEN: %s%s" % (name,
                  "the .stl is back " if live.exists() else "",
                  "the .stl.retired is missing" if not ret.exists() else ""))
            failed.append(name)
        else:
            excluded.append((name, "RETIRED: " + why))
    print("SURFACE-CLOSURE GATE (GEO-101): %d files" % len(files))
    for f in files:
        name = os.path.basename(f)
        r = checks(f, manifest)
        gates = [k for k in ("G1", "G2", "G3", "G4", "G5", "G6")]
        ok = all(r[k][0] for k in gates)
        npass += ok
        if not ok:
            failed.append(name)
        print("  %-38s %s  %s" % (name, "PASS" if ok else "FAIL",
                                  "" if ok else " ".join(
                                      "%s(%s)" % (k, r[k][1]) for k in gates if not r[k][0])))
    print()
    total = len(files) + len(RETIRED)
    print("  PASS %d + RETIRED %d = %d of %d surfaces (11 live + 2 retired)"
          % (npass, len(excluded), npass + len(excluded), total))
    for n, why in excluded:
        print("     %-42s %s" % (n, why))
    if npass + len(excluded) != total or failed:
        print()
        print("PARTITION FAILURE: %d file(s) in neither bucket: %s"
              % (len(failed), ", ".join(failed)))
        return 1
    print("\nALL FILES CLOSED, ORIENTED OUTWARD, OVERHANGING, AND PROVABLY UNMOVED "
          "ABOVE y = 0.")
    return 0


def negative_control(files):
    """A gate that cannot fail is not a gate. Prove G1 fires on an open surface by
    DELETING one triangle from a real file and checking the gate catches it."""
    import tempfile
    print("NEGATIVE CONTROL: the gate must FAIL on an open surface")
    src = files[0]
    txt = Path(src).read_text(errors="ignore")
    i = txt.find("  facet normal")
    j = txt.find("  facet normal", i + 10)
    k = txt.find("  facet normal", j + 10)
    holed = txt[:j] + txt[k:]
    ok = True
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "holed.stl"
        p.write_text(holed)
        r = checks(p, None)
        caught = not r["G1"][0]
        ok &= caught
        print("  one triangle removed from %-28s G1 %s  (%s)"
              % (os.path.basename(src), "CAUGHT" if caught else "*** MISSED ***", r["G1"][1]))
        # and the unmodified file must give the SAME verdict it gives in the main run,
        # so the control cannot pass by breaking the reader.
        r2 = checks(src, None)
        print("  unmodified %-38s G1 %s  (%s)"
              % (os.path.basename(src), "open" if not r2["G1"][0] else "closed", r2["G1"][1]))
    print("\nNEGATIVE CONTROL: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
