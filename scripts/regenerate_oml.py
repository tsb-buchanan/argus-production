#!/usr/bin/env python3
"""Regeneration gate for the production 3D OMLs.

RULING (D050): the derived STLs are git-ignored (.gitignore line 13, *.stl), and
the fix is NOT to track binaries. If the promoted surface regenerates
BIT-FOR-BIT from its recorded recipe, the binary is a CACHE and the recipe plus
checksum is the artifact. That is stronger than storing a blob, because it
proves the chain rather than its output.

Chain, per candidate:
    delivered .vsp3
      -> scripts/vsp3_to_stl.py --nchord 200 --nspan 240 --blunt-te 0.0063
      -> scripts/apply_placement_transform.py --order xy   (Rx(5.000) @ Ry(1.873))

--check  regenerate into a temporary directory and compare sha256 against
         registry/candidates.yaml. Non-zero exit and a loud message on any
         mismatch. This is what regenerate_results.py calls.
--write  regenerate in place (used only when a recipe change is intended).

DETERMINISM IS TESTED, NOT ASSUMED. Floating-point operation ordering in the
rotation can break bit-for-bit reproducibility; if it does, this gate says so
rather than being relaxed to a tolerance.
"""
import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]

CHAIN = {
    "baseline": {
        "vsp3": "dso_reference/baseline_wing_only_refined.vsp3",
        "out": "geometry/derived/wing3d/baseline_oml_placed.stl",
    },
    "mbr_c01_l011": {
        "vsp3": "dso_reference/mbr_c01_l011.vsp3",
        "out": "geometry/derived/wing3d/mbr_c01_l011_oml_placed.stl",
    },
    # Corrected pair, 2026-07-29 delivery. Same route, same arguments; only the
    # source .vsp3 differs, which is the point: the correction touched the
    # inserted aerofoils and nothing else (planform identical to 0.00e+00).
    "baseline_corrected": {
        "vsp3": "dso_reference/2026-07-29_correction_pair/baseline_corrected/baseline_corrected.vsp3",
        "out": "geometry/derived/wing3d/baseline_corrected_oml_placed.stl",
    },
    "mcv2_i002_c01": {
        "vsp3": "dso_reference/2026-07-29_correction_pair/mcv2_i002_c01/mcv2_i002_c01.vsp3",
        "out": "geometry/derived/wing3d/mcv2_i002_c01_oml_placed.stl",
    },
    "cte_i002_c04": {
        "vsp3": "dso_reference/2026-07-29_correction_pair/cte_i002_c04/cte_i002_c04.vsp3",
        "out": "geometry/derived/wing3d/cte_i002_c04_oml_placed.stl",
    },
    # Phase 0.5 SYNTHETIC wings (GEO-082). DIFFERENT ROUTE, and that is the point:
    # they are analytic, so there is no .vsp3 and NO PLACEMENT TRANSFORM. Handled by
    # a separate builder rather than bent into the ARGUS chain.
    "syn_elliptical": {
        "synthetic": "elliptical",
        "out": "geometry/derived/wing3d/syn_elliptical_oml.stl",
    },
    "syn_rectangular": {
        "synthetic": "rectangular",
        "out": "geometry/derived/wing3d/syn_rectangular_oml.stl",
    },
}
# DECLARED EXCLUSIONS, NOT SILENT ONES (GEO-087). These ARE registered with
# checksums, so a gate that iterates CHAIN alone covers 7 of 10 and reports "all
# registered OMLs" over an unstated subset. That is the defect gate 0c names, and
# it was in this file's own success message. Naming them here makes the partition
# CHAIN + DECLARED = registry, which is then ASSERTED below, so a future registered
# surface that lands in neither set FAILS the gate instead of vanishing from it.
DECLARED_NOT_REGENERABLE = {
    "mbr_c01_l011_a010": "amplitude-ladder blend, superseded lineage; built by a "
                         "different route (vsp3_to_stl --blend-with) that this gate "
                         "does not carry. Valid as calibration, not for new work.",
    "mbr_c01_l011_a025": "as mbr_c01_l011_a010",
    "mbr_c01_l011_a050": "as mbr_c01_l011_a010",
    # TIP-CAP-REPAIRED MESHING SURFACES (ARG-098). Excluded from THIS gate by design,
    # not by oversight: the chain above runs --cap fan so the anchors keep reproducing
    # bit-for-bit, and these are --cap delaunay --cap-edge 0.0030 plus a weld_surface.py
    # step the chain does not carry. Their reproduction status is PINNED in the registry
    # rather than merely named here, which is what an exclusion owes (GEO-092).
    "mcv2_i002_c01::geometry_3d_laddercap":
        "delaunay cap + weld_surface step, outside this fan-cap chain. Route VERIFIED "
        "bit-for-bit 2026-08-03 and recorded in registry route_status.",
    "baseline_corrected::geometry_3d_laddercap":
        "delaunay cap + weld_surface step, outside this fan-cap chain. Route VERIFIED "
        "bit-for-bit 2026-08-03 from baseline_corrected.vsp3. NOT YET MESHED BY ANY CASE.",
    "baseline::geometry_3d_laddercap":
        "delaunay cap + weld_surface step, outside this fan-cap chain. Route VERIFIED "
        "bit-for-bit from baseline_wing_only_refined.vsp3. THIS is what the whole Phase A "
        "baseline leg was meshed on, in error; see ARG-099 and registry MUST_NOT.",
    # FOUND BY THE ENUMERATOR FIX ABOVE, not by being bitten (ARG-098). It was registered
    # with a checksum and invisible to this gate for its whole life, so the gate has been
    # reporting over 10 of 11 pre-existing surfaces without saying so. It is PRE-DATING
    # my two entries: exposing it is the fix working, not a regression from it.
    "baseline::geometry_3d_superseded":
        "the UNPLACED loft, superseded by baseline_oml_placed.stl which IS in the chain. "
        "NOT regenerable by this gate BY CONSTRUCTION: current code always applies the "
        "placement transform, so the chain cannot emit a pre-transform surface. Its own "
        "route is pinned in the registry at commit 25db7509 (D050), and it is retained "
        "because the A1 noise-floor test was meshed from it.",
}
LOFT_ARGS = ["--nchord", "200", "--nspan", "240", "--blunt-te", "0.0063"]
XFORM_ARGS = ["--order", "xy"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def close_surface(path):
    """Root-overhang closure (GEO-101). PART OF THE RECIPE, not a post-hoc edit.
    The lofts are open shells at the root; the promoted files are closed solids, so
    the chain that reproduces them must include this step or the regeneration gate
    stops proving anything."""
    subprocess.run(
        [sys.executable, str(REPO / "scripts" / "close_wing_surface.py"),
         str(path), "--write"],
        cwd=REPO, check=True, capture_output=True)
    return path


def regenerate(cid, spec, workdir, cap="fan", cap_edge=None):
    if "synthetic" in spec:
        out = workdir / ("%s.stl" % cid)
        subprocess.run(
            [sys.executable, str(REPO / "geometry" / "scripts" / "make_synthetic_wing.py"),
             "--kind", spec["synthetic"], "--out", str(out),
             "--nchord", "200", "--nspan", "240"],
            cwd=REPO, check=True, capture_output=True)
        return close_surface(out)
    loft = workdir / ("%s_loft.stl" % cid)
    placed = workdir / ("%s_placed.stl" % cid)
    subprocess.run(
        [sys.executable, str(REPO / "scripts" / "vsp3_to_stl.py"),
         str(REPO / spec["vsp3"]), "--out", str(loft)] + LOFT_ARGS
        + ["--cap", cap] + (["--cap-edge", str(cap_edge)] if cap_edge else []),
        cwd=REPO, check=True, capture_output=True)
    subprocess.run(
        [sys.executable, str(REPO / "scripts" / "apply_placement_transform.py"),
         "--loft", str(loft), "--out", str(placed)] + XFORM_ARGS,
        cwd=REPO, check=True, capture_output=True)
    return close_surface(placed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--write", action="store_true")
    # PRODUCE MODE. Default is unchanged so the bit-for-bit gate keeps working:
    # --cap defaults to "fan", which is what every registered checksum was built
    # with. --emit switches this from a GATE (does the chain reproduce the
    # registry?) into a BUILDER (run the same chain, write somewhere new), so the
    # repaired geometry comes off the SAME chain definition rather than a second
    # copy of it that will not learn when this one changes.
    ap.add_argument("--emit", type=Path, default=None,
                    help="produce mode: write the chain's output here instead of "
                         "checking it against the registry")
    ap.add_argument("--cap", choices=["fan", "delaunay"], default="fan")
    ap.add_argument("--cap-edge", type=float, default=None)
    args = ap.parse_args()
    if not (args.check or args.write or args.emit):
        ap.error("pass --check, --write or --emit")

    if args.emit:
        args.emit.mkdir(parents=True, exist_ok=True)
        print("PRODUCE MODE: cap=%s  ->  %s" % (args.cap, args.emit))
        print("NOT a bit-for-bit gate: --cap %s changes the tip by design.\n" % args.cap)
        rows = []
        with tempfile.TemporaryDirectory(prefix="argus_oml_") as td:
            wd = Path(td)
            for cid, spec in CHAIN.items():
                if "synthetic" in spec:
                    print("   %-22s SKIPPED: analytic route, its cap lives in "
                          "make_synthetic_wing.py (also a fan; not fixed here)" % cid)
                    continue
                out = regenerate(cid, spec, wd, cap=args.cap, cap_edge=args.cap_edge)
                dst = args.emit / ("%s_oml_placed.stl" % cid)
                shutil.copy(out, dst)
                h = sha256(dst)
                rows.append((cid, dst.name, h))
                print("   %-22s %s  %s" % (cid, h[:16], dst.name))
        print("\n  registry sha256 values, for candidates.yaml:")
        for cid, nm, h in rows:
            print("    %-22s sha256: %s" % (cid, h))
        print("\n  regenerate_with: scripts/regenerate_oml.py --emit <dir> --cap %s"
              % args.cap)
        return 0

    reg = yaml.safe_load(open(REPO / "registry" / "candidates.yaml"))

    # COVERAGE ASSERTION (GEO-087): a gate must assert it saw everything it claims
    # to have checked. A verdict over an unstated subset is not a verdict.
    def registered(node, path=""):
        if isinstance(node, dict):
            for k, v in node.items():
                # ANY geometry_3d* KEY IS A REGISTERED SURFACE, not just the bare one.
                # Matching k == "geometry_3d" EXACTLY meant geometry_3d_laddercap -- the
                # surface every 3D case actually meshes -- was enumerated by NOTHING, so
                # it could not even be reported as unaccounted. The DECLARED_NOT_
                # REGENERABLE comment promises a new registered surface FAILS rather than
                # vanishes; that held for a new CANDIDATE and silently did not for a new
                # KEY on an existing candidate (ARG-098, GEO-087 species: a gate that
                # narrows its own input set and reports as though over the whole).
                if k.startswith("geometry_3d") and isinstance(v, dict) and v.get("sha256"):
                    yield (path or "?") if k == "geometry_3d" \
                        else "%s::%s" % (path or "?", k)
                else:
                    yield from registered(v, path or k)
        elif isinstance(node, list):
            for it in node:
                if isinstance(it, dict) and it.get("sha256") and "candidate_id" in it:
                    yield it["candidate_id"]
                else:
                    yield from registered(it, path)

    reg_surfaces = sorted(set(registered(reg)))
    covered = set(CHAIN) | set(DECLARED_NOT_REGENERABLE)
    unaccounted = [s for s in reg_surfaces if s not in covered]
    print("COVERAGE: %d registered surfaces = %d regenerated + %d declared "
          "not-regenerable + %d unaccounted"
          % (len(reg_surfaces), len(CHAIN), len(DECLARED_NOT_REGENERABLE),
             len(unaccounted)))
    for cid, why in sorted(DECLARED_NOT_REGENERABLE.items()):
        print("   declared %-20s %s" % (cid, why[:88]))
    if unaccounted:
        print()
        print("COVERAGE FAILURE: these are registered with a checksum but appear in")
        print("NEITHER the regeneration chain NOR the declared exclusions, so this")
        print("gate would have reported a verdict over a subset without saying so:")
        for cid in unaccounted:
            print("   %s" % cid)
        return 1
    print()

    failures = []

    with tempfile.TemporaryDirectory(prefix="argus_oml_") as td:
        workdir = Path(td)
        for cid, spec in CHAIN.items():
            expected = reg[cid]["geometry_3d"]["sha256"]
            placed = regenerate(cid, spec, workdir)
            got = sha256(placed)
            ok = got == expected
            print("%-14s regenerated %s  registry %s  %s"
                  % (cid, got[:16], expected[:16], "MATCH" if ok else "MISMATCH"))
            if not ok:
                failures.append((cid, expected, got))
            if args.write:
                shutil.copy2(placed, REPO / spec["out"])

    if failures:
        print()
        print("REGENERATION GATE FAILED. The promoted binary is NOT reproducible from")
        print("its recorded recipe, so the recipe is not a substitute for the file.")
        for cid, exp, got in failures:
            print("   %s: registry %s but regenerated %s" % (cid, exp, got))
        return 1

    print()
    print("REGENERATION GATE PASSED: all %d surfaces in the regeneration chain "
          "reproduce bit-for-bit, and the chain plus the %d declared exclusions "
          "account for every one of the %d registered surfaces."
          % (len(CHAIN), len(DECLARED_NOT_REGENERABLE), len(reg_surfaces)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
