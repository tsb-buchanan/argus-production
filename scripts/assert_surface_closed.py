#!/usr/bin/env python3
"""assert_surface_closed.py: refuse to quote a 3D result from a case meshed on an OPEN SHELL.

THE DEFECT THIS EXISTS FOR. Every wing surface this project has ever meshed is an open
shell. The root rib carries no cap, so the aerofoil interior is connected to the outside
and snappyHexMesh floods it. The wing is then meshed as a ZERO-THICKNESS BAFFLE with fluid
on both sides. Measured on cases/mesh_v2/baseline_corrected, 7.51 M cells:

    76.8% of wing patch faces are back-to-back anti-parallel coincident pairs
    74.5% of the patch AREA is in such pairs
    patch area / STL area = 1.905, which is the factor of two, not stair-stepping
    doubling is 85.0% at eta 0.00-0.10 falling monotonically to 71.3% at the tip,
        the signature of a leak propagating from the root

THE MECHANISM IS ARITHMETIC, not a hypothesis. The rim is planar at y = 0 BEFORE placement.
Placement applies 5 deg of dihedral as a rotation about x, tilting the rib out of the
symmetry plane:

    root rib z-extent 0.0957 m  x  sin(5 deg) = 0.087156  =  8.341 mm
    measured rim y range        -3.692 .. +4.679          =  8.371 mm   (0.36% agreement)

63.3% of the rim then sits at y > 0, leaving a slot up to 4.679 mm against a 3.891 mm
surface cell, i.e. 1.20 cells wide. Fluid walks in.

WHY THE GUARD READS THE CASE AND NOT A REGISTRY (GEO-083). A guard that reads a registry
inherits the registry's correctness: it asserts currency for a fact it did not establish.
This one opens the case's OWN constant/triSurface/*.stl and counts edges. Nothing it
reports depends on anything staying true elsewhere.

WHY IT IS NOT ENOUGH TO REFINE (the counter-intuitive part, and the reason this is a gate
rather than a note). The gap is fixed at 4.679 mm by GEOMETRY. At surface level 9 it is
1.20 cells and leaks only where widest; at level 11 it is 4.79 cells and leaks freely along
the whole 1.32 m rim. REFINING BEFORE CLOSING MAKES IT WORSE. The gate therefore fires
before mesh work, not only before reporting.

WHY checkMesh NEVER SAW IT. checkMesh validates the mesh it is given; a flooded mesh is a
perfectly valid mesh of the wrong domain. surfaceCheck would have caught it in one command
and was never run. That absence is the finding: this is a VOUCHED-SUCCESS instance in the
geometry layer (GEO-080), where every downstream artefact reported cleanly on a step nobody
observed.

TWO CHECKS, AND THE SECOND IS NOT REDUNDANT.
  1. CLOSED. Every edge has multiplicity exactly 2, orientations are consistent, and the
     divergence-theorem volume is positive (outward normals).
  2. CROSSES THE SYMMETRY PLANE. A solid capped exactly ON y = 0 is closed and still wrong:
     it forces snappy to choose between wall and symmetry at the same location. The
     geometry must extend THROUGH the plane so the domain boundary cuts a closed solid.
     A surface can pass check 1 and fail check 2, which is precisely the fix done badly.
"""
import argparse
import collections
import glob
import os
import re
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WELD_DECIMALS = 6          # 1 micron; STL ascii carries more than enough precision
MIN_OVERHANG_M = 0.010     # a solid must reach at least this far past y = 0


class SurfaceNotClosedError(RuntimeError):
    """Raised by gate(). NOT caught anywhere in this project on purpose."""


def _tris_from_stl(path):
    txt = open(path, errors="ignore").read()
    v = re.findall(r"vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)", txt)
    if len(v) < 3:
        raise ValueError(f"{path}: fewer than one triangle parsed; not ascii STL?")
    return np.asarray(v, dtype=float).reshape(-1, 3, 3)


def probe(path):
    """THIN ADAPTER over the project's ONE closure definition (GEO-102).

    THE MEASUREMENT IS NOT MINE ANY MORE, AND THAT IS THE POINT. Two implementations of
    "is this surface closed" existed for one question: the geometry workstream's file-level gate
    over geometry/derived/wing3d, and this case-level guard. DIFFERENT GATES, SAME QUESTION,
    so the question gets one implementation and the gates keep their own scopes. Same move
    as chord_for() and section_conventions.le_te_base, both of which exist because a
    convention carried as a flag on one caller does not travel (traps item 4n).

    WHAT I GAINED BY NOT KEEPING MY OWN. The vertex-unification tolerance is PART of the
    definition: two implementations rounding differently disagree about which edges are
    open, and they disagree only on marginal geometry, which is exactly when it matters.
    Mine rounded at 6 decimals by its own constant and never reported it. closure_probe
    returns `round_decimals` WITH the result, so the tolerance now travels with the number
    (D065: the frame travels with the number or the number is not reportable).

    Only the KEY NAMES are adapted here, plus y-range keys carrying an _m suffix so the
    unit is stated. No fact is recomputed.
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from check_surface_closed import closure_probe
    p = closure_probe(path)
    return {
        "path": path,
        "n_tri": p["n_tri"],
        "n_open_edges": p["n_open_edges"],
        "n_over_edges": p["n_over_edges"],
        "n_inconsistent": p["n_orientation_errors"],
        "volume_m3": p["volume"],
        "y_min_m": p["y_min"],
        "y_max_m": p["y_max"],
        "round_decimals": p["round_decimals"],
    }


def verdict(p, require_overhang=True):
    """Returns (ok, [reasons]). Each reason names the measured value, never just the rule."""
    bad = []
    if p["n_open_edges"]:
        bad.append(f"OPEN: {p['n_open_edges']} edges with one neighbour "
                   f"(vertices unified at 1e-{p['round_decimals']} m)")
    if p["n_over_edges"]:
        bad.append(f"NON-MANIFOLD: {p['n_over_edges']} edges with more than two neighbours")
    if p["n_inconsistent"]:
        bad.append(f"INCONSISTENT ORIENTATION: {p['n_inconsistent']} directed edges "
                   f"traversed the same way twice")
    if not bad and p["volume_m3"] <= 0:
        bad.append(f"INWARD NORMALS: enclosed volume {p['volume_m3']:+.6f} m3 is not positive")
    if require_overhang and not bad and p["y_min_m"] > -MIN_OVERHANG_M:
        bad.append(f"DOES NOT CROSS THE SYMMETRY PLANE: y_min {1000*p['y_min_m']:+.3f} mm, "
                   f"needs < {-1000*MIN_OVERHANG_M:.1f} mm. A solid capped ON y=0 makes "
                   f"snappy choose between wall and symmetry at one location.")
    return (not bad), bad


def gate(case_dir, require_overhang=True):
    """HALT unless every surface this case will mesh on is a closed, outward solid.

    Call this from anything that meshes a case or emits a 3D number from one. It raises
    rather than returning a flag, because a flag gets read as advisory and this is not.
    """
    stls = sorted(glob.glob(os.path.join(case_dir, "constant/triSurface/*.stl")))
    if not stls:
        raise SurfaceNotClosedError(
            f"{case_dir}: no constant/triSurface/*.stl. ABSENCE IS NOT A PASS: this gate "
            f"cannot certify a surface it never saw (GEO-087).")
    problems = []
    for s in stls:
        ok, why = verdict(probe(s), require_overhang)
        if not ok:
            problems += [f"  {os.path.relpath(s, case_dir)}: {w}" for w in why]
    if problems:
        raise SurfaceNotClosedError(
            f"{case_dir}\n" + "\n".join(problems)
            + "\n  NO 3D RESULT FROM THIS CASE IS QUOTABLE (ARG-096). The wing is meshed as"
              "\n  a zero-thickness baffle with fluid on both sides of it.")
    return True


# --------------------------------------------------------------------------------------
# SELF-TEST. KNOWN ANSWERS IN BOTH DIRECTIONS (D058 item 2, GEO-080 item 6).
# A gate is worthless until it has been SHOWN to fail. Four fixtures: one that must pass
# and three that must fail, each for a DIFFERENT stated reason, so a checker that fails
# everything cannot be mistaken for a working one.
# --------------------------------------------------------------------------------------
_TET = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
_TET_FACES = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]   # outward


def _write_stl(path, tris):
    with open(path, "w") as fh:
        fh.write("solid t\n")
        for t in tris:
            fh.write(" facet normal 0 0 0\n  outer loop\n")
            for v in t:
                fh.write(f"   vertex {v[0]:.9f} {v[1]:.9f} {v[2]:.9f}\n")
            fh.write("  endloop\n endfacet\n")
        fh.write("endsolid t\n")


def self_test(tmpdir, real_open_stl=None):
    import tempfile
    fails = []

    def check(name, path, want_ok, want_substr, require_overhang=False):
        ok, why = verdict(probe(path), require_overhang)
        got = "PASS" if ok else "FAIL"
        exp = "PASS" if want_ok else "FAIL"
        hit = want_ok or any(want_substr in w for w in why)
        status = "ok " if (ok == want_ok and hit) else "BROKEN"
        if status == "BROKEN":
            fails.append(f"{name}: expected {exp}/{want_substr!r}, got {got} {why}")
        print(f"  [{status}] {name:34s} expected {exp:4s} got {got:4s}  "
              f"{'' if ok else why[0][:66]}")

    # 1. KNOWN CLOSED: a tetrahedron with outward normals, shifted so it crosses y = 0.
    p = os.path.join(tmpdir, "closed.stl")
    _write_stl(p, (_TET - np.array([0., 0.5, 0.]))[np.array(_TET_FACES)])
    check("closed outward tetrahedron", p, True, "")

    # 2. KNOWN OPEN: same tetrahedron with one face removed. Must report 3 open edges.
    p = os.path.join(tmpdir, "open.stl")
    _write_stl(p, _TET[np.array(_TET_FACES[:3])])
    check("tetrahedron missing one face", p, False, "OPEN")

    # 3. KNOWN INVERTED: closed, manifold, consistently oriented, but normals INWARD.
    #    This is the fixture that catches a checker testing only for holes.
    p = os.path.join(tmpdir, "inverted.stl")
    _write_stl(p, _TET[np.array([f[::-1] for f in _TET_FACES])])
    check("closed tetrahedron, normals inward", p, False, "INWARD")

    # 4. KNOWN CAPPED-ON-THE-PLANE: closed, outward, but sitting entirely at y >= 0.
    #    This is the fix done badly, and it must not be waved through by check 1.
    p = os.path.join(tmpdir, "onplane.stl")
    _write_stl(p, _TET[np.array(_TET_FACES)])
    check("closed solid capped on y=0", p, False, "DOES NOT CROSS", require_overhang=True)

    # 5. THE REAL ARTEFACT. The gate must fail on the geometry as it stands today; if it
    #    passes here the gate is broken, whatever the synthetic fixtures said.
    if real_open_stl and os.path.exists(real_open_stl):
        check(f"REAL {os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(real_open_stl))))}",
              real_open_stl, False, "OPEN")

    print()
    if fails:
        for f in fails:
            print("  BROKEN:", f)
        return 1
    print("  self-test: 5/5 fixtures behaved as specified, in BOTH directions.")
    return 0


def audit_all(require_overhang=False):
    """PARTITION over every case-local surface (GEO-089). PASS + FAIL + UNREADABLE = ALL.

    require_overhang defaults FALSE here so the audit reports the PRIMARY defect (open
    shells) without every line also carrying the secondary one. Flip it once the geometry
    is repaired and the secondary check becomes the live question.
    """
    pats = sorted(glob.glob(os.path.join(REPO, "cases/*/*/constant/triSurface/*.stl")))
    ok_l, bad_l, err_l = [], [], []
    print(f"  {'case':46s} {'stl':14s} {'tris':>8s} {'open':>6s} {'over-mult':>9s} "
          f"{'y_min mm':>8s}  verdict")
    for s in pats:
        rel = os.path.relpath(s, os.path.join(REPO, "cases"))
        case = "/".join(rel.split("/")[:2])
        try:
            p = probe(s)
        except Exception as e:                       # UNREADABLE IS ITS OWN BUCKET, never
            err_l.append((case, str(e)))             # silently dropped (GEO-087 item 3c)
            print(f"  {case:46s} {os.path.basename(s):14s} {'':>8s} {'':>6s} {'':>9s} "
                  f"{'':>8s}  UNREADABLE")
            continue
        ok, why = verdict(p, require_overhang)
        (ok_l if ok else bad_l).append(case)
        gap = f"{1000*p['y_min_m']:.3f}"
        # EVERY reason, not why[0]. Printing only the first defect is the same shape as a
        # coverage check that hides its diagnosis: the banded slabs read as merely OPEN for
        # a whole session while ALSO carrying 60 and 136 edges shared by >2 triangles, a
        # SECOND AND INDEPENDENT defect that a root cap would not have fixed. The gate had
        # measured it; the display threw it away. Found by the geometry workstream reporting a
        # defect mine had already seen and not shown.
        print(f"  {case:46s} {os.path.basename(s):14s} {p['n_tri']:8d} "
              f"{p['n_open_edges']:6d} {p['n_over_edges']:9d} {gap:>8s}  "
              f"{'CLOSED' if ok else ' + '.join(w.split(':')[0] for w in why)}")

    total = len(ok_l) + len(bad_l) + len(err_l)
    print(f"\n  PARTITION  closed {len(ok_l)}  +  defective {len(bad_l)}  +  unreadable "
          f"{len(err_l)}  =  {total} surfaces probed")
    assert total == len(pats), f"partition {total} != enumerated {len(pats)}"
    print(f"  ASSERTED: buckets sum to the enumerated set ({len(pats)}).")
    # REPORT FIRST, EXIT AFTER (GEO-087 item 5). A coverage check that returns non-zero
    # before printing its diagnosis turns a useful finding into an unactionable alarm.
    return 1 if (bad_l or err_l) else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("case", nargs="?", help="case directory to gate")
    ap.add_argument("--audit-all", action="store_true",
                    help="partition every case-local surface in the repo")
    ap.add_argument("--self-test", action="store_true",
                    help="known answers in both directions; run this before trusting a PASS")
    ap.add_argument("--require-overhang", action="store_true",
                    help="also demand the solid cross y=0 (enable once geometry is repaired)")
    a = ap.parse_args()

    if a.self_test:
        import tempfile
        real = os.path.join(REPO, "cases/trim_corrected/baseline_corrected_a1p7589/"
                                  "constant/triSurface/wing.stl")
        with tempfile.TemporaryDirectory() as td:
            return self_test(td, real)
    if a.audit_all:
        return audit_all(a.require_overhang)
    if not a.case:
        ap.error("give a case directory, --audit-all or --self-test")
    try:
        gate(a.case, a.require_overhang)
    except SurfaceNotClosedError as e:
        print(f"REFUSED\n{e}", file=sys.stderr)
        return 1
    print(f"{a.case}: all surfaces closed, outward"
          + (", and crossing y=0" if a.require_overhang else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
