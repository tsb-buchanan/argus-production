#!/usr/bin/env python3
"""weld_surface.py: merge near-coincident vertices so OpenFOAM sees ONE closed surface.

THE DEFECT. surfaceCheck on the repaired OMLs reports:
    Surface is not closed since not all edges connected to two faces
    Number of zones (connected area with consistent normal) : 2
    More than one normal orientation.
while three separate Python gates all reported it closed, because every one of them welds
vertices at 1e-6 m BEFORE counting edges. OPENFOAM DOES NOT. Measured: OpenFOAM reads
96,316 vertices and the EXACT-unique count is 96,316 too, so org-7's ascii STL reader
merges only byte-identical points. A closed genus-0 triangulation of 191,834 triangles
needs 95,919 vertices, so 397 are unmerged duplicates -- matching the 398-point root rim.

THE CAUSE, MEASURED NOT ASSUMED. 399 vertex pairs lie within 1e-5 m of each other,
separated by 68 nm to 10 um. The root skirt was placed AGAINST the rim rather than welded
INTO it, so each rim point exists twice at coordinates that differ in the last few digits
of the %.9e text. Geometrically coincident, topologically distinct.

WHY IT MATTERED ANYWAY, given the mesh came out sealed. snappyHexMesh's inside/outside test
is ray casting against triangles, which a geometrically closed surface passes whether or not
it is topologically welded -- which is why area ratio came out 1.0015 and the baffle
fraction 0.128%. THAT WAS LUCK, NOT DESIGN. surfaceFeatures extracts feature edges from
TOPOLOGY, so an unwelded seam is a feature line the mesher will faithfully refine, and the
two-zone normal inconsistency makes the inside test orientation-dependent.

WHAT THIS MOVES, AND IT IS DELIBERATELY ASYMMETRIC. Within each cluster the representative
is the member with the LARGEST y, so the SKIRT (added, at y < 0) snaps UP onto the WING
(original). The aerodynamic surface is the thing being validated and it does not move except
where two original points were already duplicated. Maximum displacement is reported, and it
is bounded by the cluster tolerance.

VALIDATE BEFORE WRITE (GEO-080 item 6). The welded surface is written to a TEMPORARY path,
surfaceCheck is run on it, and it is promoted to the real filename ONLY if OpenFOAM itself
now calls it closed. Nothing invalid ever reaches the destination name.

WHICH surfaceCheck. OpenFOAM-org 12, taken from (1) an OpenFOAM-12 environment already
loaded in the calling shell, else (2) the etc/bashrc named by $ARGUS_OF_BASHRC, else
(3) /opt/openfoam12/etc/bashrc. The version is read back from surfaceCheck's own Build line.
If surfaceCheck cannot be run under OpenFOAM-12, the script stops with exit status 2 and says
so: a check that did not run is never reported as a surface that failed it.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

DEFAULT_BASHRC = "/opt/openfoam12/etc/bashrc"
TOL = 1.0e-5          # cluster radius; measured max real separation is 9.989e-06


class SurfaceCheckUnavailable(RuntimeError):
    """surfaceCheck did not run under OpenFOAM-12, so nothing was checked."""


def openfoam_source():
    """Return (bashrc to source, or None to use the loaded environment; description).

    A set ARGUS_OF_BASHRC that does not exist is an error. It never falls through to the
    default, so a wrong override cannot be silently replaced by another installation.
    """
    override = os.environ.get("ARGUS_OF_BASHRC")
    if os.environ.get("WM_PROJECT_VERSION") == "12" and os.environ.get("WM_PROJECT_DIR"):
        how = "already loaded (WM_PROJECT_DIR=%s)" % os.environ["WM_PROJECT_DIR"]
        if override:
            how += "; ARGUS_OF_BASHRC ignored"
        return None, how
    if override:
        bashrc, origin = os.path.expanduser(override), "ARGUS_OF_BASHRC"
    else:
        bashrc, origin = DEFAULT_BASHRC, "default"
    if not os.path.isfile(bashrc):
        raise SurfaceCheckUnavailable(
            "OpenFOAM-12 is not loaded in this shell and its bashrc was not found at %s (%s)"
            % (bashrc, origin))
    return bashrc, "sourced %s (%s)" % (bashrc, origin)


def surface_check(path, bashrc):
    """Run surfaceCheck on path. Returns (verdict, OpenFOAM build).

    Raises SurfaceCheckUnavailable unless surfaceCheck ran to completion under OpenFOAM-12.
    """
    if bashrc is None:
        cmd = ["surfaceCheck", str(path)]
    else:
        # a non-login shell without rc files, so no profile can load another OpenFOAM.
        # `set --` first: OpenFOAM's bashrc reads the caller's arguments as settings (org-7
        # sources any that are files, which here recursed into itself and hung).
        cmd = ["bash", "--noprofile", "--norc", "-c",
               'rc="$1"; stl="$2"; set --; source "$rc" >/dev/null 2>&1; '
               'command -v surfaceCheck >/dev/null || exit 127; exec surfaceCheck "$stl"',
               "weld_surface", bashrc, str(path)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except FileNotFoundError:
        raise SurfaceCheckUnavailable("surfaceCheck is not on PATH although "
                                      "WM_PROJECT_VERSION is 12")
    except subprocess.TimeoutExpired:
        raise SurfaceCheckUnavailable("surfaceCheck did not finish within 900 s")
    o = r.stdout + r.stderr
    if bashrc is not None and r.returncode == 127:
        raise SurfaceCheckUnavailable("surfaceCheck not found after sourcing %s" % bashrc)
    build = re.search(r"^Build\s*:\s*(\S+)", o, re.M)
    if r.returncode != 0 or "Reading surface" not in o or not build:
        tail = " | ".join(o.strip().splitlines()[-4:]) or "no output"
        raise SurfaceCheckUnavailable("surfaceCheck on %s did not complete (exit status %d): %s"
                                      % (path, r.returncode, tail))
    if not re.match(r"12(\D|$)", build.group(1)):
        raise SurfaceCheckUnavailable("surfaceCheck ran under OpenFOAM build %s, not OpenFOAM-12"
                                      % build.group(1))
    z = re.search(r"Number of zones[^:]*:\s*(\d+)", o)
    return ({"closed": "Surface is not closed" not in o,
             "zones": int(z.group(1)) if z else None,
             "one_orientation": "More than one normal orientation" not in o,
             "illegal": "no illegal triangles" not in o}, build.group(1))


def unavailable(err):
    sys.stdout.flush()
    print("\nFAILED: surfaceCheck could not be run, so nothing was checked and nothing was "
          "promoted.\n  reason: %s\n  if OpenFOAM-12 is missing: load it (source %s) or set "
          "ARGUS_OF_BASHRC to its etc/bashrc, and rerun." % (err, DEFAULT_BASHRC),
          file=sys.stderr)
    return 2


def weld(src: Path, dst: Path, tol=TOL):
    txt = src.read_text(errors="ignore")
    name = (re.search(r"^\s*solid\s+(\S+)", txt, re.M) or [None, "wing"])[1]
    V = np.array(re.findall(r"vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)", txt),
                 dtype=float)
    tris = V.reshape(-1, 3, 3)
    u, inv = np.unique(V, axis=0, return_inverse=True)

    # cluster near-coincident vertices; representative = LARGEST y so the skirt snaps
    # up onto the wing rather than the wing down onto the skirt
    tree = cKDTree(u)
    parent = np.arange(len(u))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in tree.query_pairs(tol, output_type="ndarray"):
        ra, rb = find(a), find(b)
        if ra != rb:
            hi, lo = (ra, rb) if u[ra][1] >= u[rb][1] else (rb, ra)
            parent[lo] = hi
    roots = np.array([find(i) for i in range(len(u))])
    moved = np.linalg.norm(u - u[roots], axis=1)
    n_merged = int((roots != np.arange(len(u))).sum())

    welded = u[roots][inv].reshape(-1, 3, 3)
    # drop triangles that became degenerate (two vertices merged into one)
    keep = ~((np.all(welded[:, 0] == welded[:, 1], 1)) |
             (np.all(welded[:, 1] == welded[:, 2], 1)) |
             (np.all(welded[:, 0] == welded[:, 2], 1)))
    out = welded[keep]

    with open(dst, "w") as fh:
        fh.write(f"solid {name}\n")
        for t in out:
            n = np.cross(t[1] - t[0], t[2] - t[0])
            m = np.linalg.norm(n)
            n = n / m if m > 0 else n
            fh.write(f" facet normal {n[0]:.9e} {n[1]:.9e} {n[2]:.9e}\n  outer loop\n")
            for v in t:
                fh.write(f"   vertex {v[0]:.9e} {v[1]:.9e} {v[2]:.9e}\n")
            fh.write("  endloop\n endfacet\n")
        fh.write(f"endsolid {name}\n")
    return {"vertices_before": len(u), "merged": n_merged,
            "vertices_after": len(u) - n_merged,
            "max_displacement_m": float(moved.max()),
            "tris_before": len(tris), "tris_after": int(keep.sum()),
            "tris_dropped": int((~keep).sum())}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("stls", nargs="+")
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--in-place", action="store_true",
                    help="promote over the original ONLY if surfaceCheck then passes")
    a = ap.parse_args()
    try:
        bashrc, how = openfoam_source()
    except SurfaceCheckUnavailable as e:
        return unavailable(e)
    shown = False
    bad = 0
    for s in a.stls:
        src = Path(s)
        print(f"\n=== {src} ===")
        try:
            before, build = surface_check(src, bashrc)
        except SurfaceCheckUnavailable as e:
            return unavailable(e)
        if not shown:
            print(f"  surfaceCheck: OpenFOAM build {build}, {how}")
            shown = True
        print(f"  before: {before}")
        tmp = Path(tempfile.mkstemp(suffix=".stl", prefix="weld_")[1])
        st = weld(src, tmp, a.tol)
        print(f"  merged {st['merged']} vertices ({st['vertices_before']} -> "
              f"{st['vertices_after']}), max displacement "
              f"{1e9*st['max_displacement_m']:.1f} nm")
        print(f"  triangles {st['tris_before']} -> {st['tris_after']} "
              f"({st['tris_dropped']} degenerate dropped)")
        try:
            after, _ = surface_check(tmp, bashrc)
        except SurfaceCheckUnavailable as e:
            tmp.unlink()
            return unavailable(e)
        print(f"  after:  {after}")
        ok = after["closed"] and after["zones"] == 1 and after["one_orientation"]
        if not ok:
            bad += 1
            print(f"  NOT PROMOTED. surfaceCheck still refuses it; the temporary is at {tmp}")
            continue
        if a.in_place:
            shutil.move(str(tmp), str(src))
            print(f"  PROMOTED over {src}")
        else:
            print(f"  would promote (rerun with --in-place). temporary at {tmp}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
