#!/usr/bin/env python3
"""assert_case_postpro.py: will this case produce what the post-processing needs?

AN INPUT GATE, RUN BEFORE SUBMISSION. Every gate this project already has looks at
what a case PRODUCED. checkMesh passed every flooded 3D mesh ever built here, because
a flooded mesh is a perfectly valid mesh of the wrong domain (GEO-102). The same shape
applies to post-processing: a chain that runs cleanly on whatever it was handed cannot
tell you the case never sampled the field you needed.

WHAT IT COSTS TO FIND OUT LATE. The wall-shear VECTOR is not in any existing 3D export
-- results/vtk_export/*.vtk carries cf as a MAGNITUDE, because the export was written
for making pictures. The hinge-moment viscous term needs the vector. So that term is
`null` on every case this project has run, and on the 36 cases whose processor
directories were purged on 2026-08-21 it is unrecoverable at any price. A field costs
nothing to write during a solve and cannot be bought afterwards.

TWO MODES, because the two questions are different:
    --before   read system/controlDict and say what WILL be written   (input gate)
    --after    read postProcessing/ and say what WAS written          (output check)
A case can pass the first and fail the second: a function object can be configured and
still produce nothing if the patch name is wrong. Neither substitutes for the other.

REPORT FIRST, EXIT AFTER (GEO-089). An earlier check in this repo returned non-zero
before printing which items failed, announcing a count while suppressing the diagnosis.

Usage:
    python3 scripts/assert_case_postpro.py --before <case> [<case> ...]
    python3 scripts/assert_case_postpro.py --after  <case> [--time 5000]
    python3 scripts/assert_case_postpro.py --self-test
"""
import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Each requirement names the CONSUMER that needs it. A requirement with no consumer
# should not exist; a consumer with no requirement here is an unguarded dependency.
REQUIREMENTS = [
    dict(key="forces_split",
         consumer="hpc/make_run_card.py, scripts/make_validation_run_card.py",
         why="results.CD_pressure and CD_viscous. forceCoeffs CANNOT supply them: its "
             "Cd(f)/Cd(r) are FRONT and REAR, not pressure and viscous.",
         before=lambda t: re.search(r"type\s+forces\s*;", t) is not None,
         after=lambda d: any(d.glob("*orce*/*/forces*.dat"))),

    dict(key="wing_surface_p",
         consumer="scripts/hinge_moment_integrate.py, scripts/export_wing_vtk.py",
         why="pressure on the wing, for every surface integral",
         before=lambda t: _surface_has(t, "p"),
         after=lambda d: _after_surface(d, "p")),

    dict(key="wing_surface_tau_VECTOR",
         consumer="scripts/hinge_moment_integrate.py (viscous term), "
                  "scripts/wing_surface_maps.py (separation)",
         why="THE VECTOR, not the magnitude. Reverse flow is the SIGN of tau, which |cf| "
             "cannot express, and the viscous moment needs the direction. This is the "
             "requirement every existing 3D case fails.",
         before=lambda t: _surface_has(t, "wallShearStress"),
         after=lambda d: _after_surface(d, "wallShearStress")),

    dict(key="wing_surface_yplus",
         consumer="run cards (wall_treatment.yplus_*)",
         why="a wall-treatment claim is only checkable against a measured y+",
         before=lambda t: _surface_has(t, "yPlus"),
         after=lambda d: _after_surface(d, "yPlus")),

    dict(key="surface_format_carries_topology",
         consumer="scripts/hinge_moment_integrate.py",
         why="face AREAS and NORMALS come from the polygons. `raw` writes face-centre "
             "POINTS ONLY, so no surface integral can be formed from it.",
         before=lambda t: re.search(r"surfaceFormat\s+(vtk|vtp|ensight)\s*;", t) is not None,
         after=lambda d: any(d.glob("*urface*/**/*.vtk")) or any(d.glob("*urface*/**/*.vtp"))),

    dict(key="compressible_cp_frame",
         consumer="every cp figure and every surface integral formed from cp",
         why="FOR A COMPRESSIBLE SOLVER p IS ABSOLUTE PASCALS. `pInf 0` leaves the whole "
             "freestream in the numerator and `rhoInf 1` drops the freestream density "
             "out of the dynamic pressure, 16.6% on the ONERA M6 case. Both are correct "
             "for simpleFoam, where p is kinematic gauge, and the template carried them "
             "unconditionally until 2026-08-22. A cp wrong by the freestream still plots, "
             "which is why this is a gate and not a comment.",
         before=lambda t: _cp_frame_ok(t),
         after=lambda d: True),   # nothing WRITTEN can show this; it is a setup property

    dict(key="trefftz_planes",
         consumer="scripts/trefftz.py",
         why="induced drag as a VLM defines it. Sampled IN-SOLVER because the volume "
             "field is the first thing deleted when scratch fills, and CDi decays with "
             "distance so SEVERAL planes are needed to show the placement sensitivity.",
         before=lambda t: t.count("cuttingPlane") >= 2,
         after=lambda d: len(list(d.glob("*refftz*/*/*"))) >= 2),

    dict(key="residuals",
         consumer="run cards (convergence.final_residuals), D023 tail drift",
         why="org-7 calls this `residuals`; ESI calls it `solverInfo`, and porting "
             "between them has already cost this project a job. AND THE FIELD LIST IS "
             "PART OF THE REQUIREMENT, not a detail: a bare `#includeFunc residuals` "
             "writes a residuals.dat containing TIMESTAMPS AND NOTHING ELSE, with no "
             "warning. Measured on OF-12 cavity, 2026-08-26: bare -> 1 column, "
             "`residuals(U, p)` -> 4 columns. CMPB_a0p87 shipped with the bare form.",
         before=lambda t: _residuals_ok(t),
         after=lambda d: _residuals_written(d)),
]


def _residuals_written(d):
    """True only if a residuals file under postProcessing/ carries RESIDUAL COLUMNS.

    Looks in residuals*/, argusResiduals*/ (the argusPostPro block's own directory)
    and solverInfo*/. Existence is not enough: a bare `#includeFunc residuals` writes
    a residuals.dat whose rows hold the time and nothing else, which is exactly the
    silent failure the `before` check exists to prevent, so the after check must not
    pass on it either. A data row with more than one column is required.
    """
    for pat in ("residuals*/*/*.dat", "argusResiduals*/*/*.dat", "solverInfo*/*/*.dat"):
        for f in d.glob(pat):
            try:
                with open(f) as fh:
                    for line in fh:
                        if not line.strip() or line.lstrip().startswith("#"):
                            continue
                        if len(line.split()) > 1:
                            return True
            except OSError:
                continue
    return False


def _residuals_ok(t):
    """True only if residuals are configured WITH A FIELD LIST.

    THREE FORMS EXIST AND ONLY TWO OF THEM WORK:

        type residuals; ... fields (U p k omega);   WORKS
        #includeFunc residuals(U, p, k, omega)      WORKS
        #includeFunc residuals                      WRITES AN EMPTY FILE, SILENTLY

    The third is the dangerous one. OF-12's caseDicts template carries
    `fields  (<fieldNames>);` and supplies NO default, so the bare form leaves the
    placeholder unfilled. The run does not fail. The function object does not warn.
    A file appears at the right path with the right name and the right number of
    ROWS, and it has a single column. Anything reading it for a convergence claim
    gets nothing and cannot tell that it got nothing.

    The old predicate matched only `type residuals;`, which meant it ALSO reported
    MISSING for the perfectly good `#includeFunc residuals(U, p)` form. Right verdict
    on CMPB, wrong route, and a false negative on a correct case: exactly the pair
    GEO-080 item 5 warns about.
    """
    # Form 1: an explicit block. Require a fields entry somewhere in the same text.
    if re.search(r"type\s+(residuals|solverInfo)\s*;", t):
        if re.search(r"\bfields\s*\(([^)]*\S[^)]*)\)", t):
            return True
        return False
    # Form 2: includeFunc WITH a non-empty argument list.
    m = re.search(r"#includeFunc\s+(residuals|solverInfo)\s*\(([^)]*)\)", t)
    if m and m.group(2).strip():
        return True
    # Form 3, or absent: not acceptable.
    return False


def _cp_frame_ok(t):
    """Incompressible: pInf 0 and rhoInf 1 are right. Compressible: both must be real.

    THE CHECK IS CONDITIONAL ON THE SOLVER, so it cannot fire on the cases it does not
    apply to, and it cannot pass by not applying either: a compressible case with either
    value left at its incompressible default FAILS. That is the difference between an
    exemption and a bypass (GEO-092).
    """
    if not re.search(r"^\s*application\s+(rho|sonic)\w*\s*;", t, re.M):
        return True
    m = re.search(r"pInf\s+([-\d.eE+]+)\s*;", t)
    r = re.search(r"rhoInf\s+([-\d.eE+]+)\s*;", t)
    if not m or not r:
        return False
    return abs(float(m.group(1))) > 1.0 and abs(float(r.group(1)) - 1.0) > 1e-6


def _surface_has(txt, field):
    """Is `field` in the fields list of a `surfaces` function object?"""
    for m in re.finditer(r"type\s+surfaces\s*;(.*?)(?=\n\s*\w+\s*\n\s*\{|\Z)", txt, re.S):
        blk = m.group(1)
        f = re.search(r"fields\s*\(([^)]*)\)", blk)
        if f and re.search(r"\b%s\b" % re.escape(field), f.group(1)):
            return True
    return False


def _after_surface(d, field):
    for p in d.glob("*urface*/**/*"):
        if p.is_file() and field.lower() in p.name.lower():
            return True
    # a single vtk can carry several fields
    for p in d.glob("*urface*/**/*.vtk"):
        try:
            if re.search(r"\b%s\b" % re.escape(field), p.read_text(errors="replace")[:20000], re.I):
                return True
        except OSError:
            pass
    return False


def resolve_controldict(case):
    """controlDict with any #include pulled in, because the requirements may live in
    an included block and a check that reads only the top file would miss them."""
    cd = Path(case) / "system/controlDict"
    if not cd.exists():
        raise FileNotFoundError(cd)
    txt = cd.read_text(errors="replace")
    for m in re.finditer(r'#include(?:Func|Etc)?\s+"([^"]+)"', txt):
        inc = Path(case) / "system" / m.group(1)
        if inc.exists():
            txt += "\n" + inc.read_text(errors="replace")
    return txt


def check(case, mode, time=None):
    case = Path(case)
    ok, bad = [], []
    if mode == "before":
        txt = resolve_controldict(case)
        for r in REQUIREMENTS:
            (ok if r["before"](txt) else bad).append(r)
    else:
        d = case / "postProcessing"
        if not d.is_dir():
            print("  %s: no postProcessing/ at all" % case.name)
            return [], REQUIREMENTS
        for r in REQUIREMENTS:
            (ok if r["after"](d) else bad).append(r)
    return ok, bad


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="*")
    ap.add_argument("--before", action="store_true")
    ap.add_argument("--after", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        return self_test()
    if not a.cases or (a.before == a.after):
        ap.error("give --before or --after (not both) and at least one case")

    mode = "before" if a.before else "after"
    rc = 0
    for c in a.cases:
        print("\n=== %s (%s) ===" % (c, mode))
        try:
            ok, bad = check(c, mode)
        except FileNotFoundError as e:
            print("  cannot read %s" % e); rc = 1; continue
        for r in ok:
            print("  ok      %s" % r["key"])
        for r in bad:
            print("  MISSING %s" % r["key"])
            print("            consumer: %s" % r["consumer"])
            print("            why:      %s" % r["why"])
        # the buckets must partition the requirement list (GEO-089)
        assert len(ok) + len(bad) == len(REQUIREMENTS), "buckets do not partition"
        print("  ---- %d of %d requirements satisfied ----" % (len(ok), len(REQUIREMENTS)))
        if bad:
            rc = 1
    return rc


def self_test():
    """Known answers in BOTH directions: a case that satisfies a requirement and one
    that does not. A gate only ever exercised on passing input is not a gate."""
    import tempfile
    ok = True
    good = """
functions {
  f1 { type forces; patches (wing); }
  s1 { type surfaces; surfaceFormat vtk; fields (p wallShearStress yPlus);
       surfaces { wingSurface { type patch; patches (wing); } } }
  s2 { type surfaces; surfaceFormat vtk; fields (U);
       surfaces { a { type cuttingPlane; } b { type cuttingPlane; } } }
  r1 { type residuals; fields (U p); }
}"""
    bad = """
functions {
  f1 { type forceCoeffs; patches (wing); }
  s1 { type surfaces; surfaceFormat raw; fields (p);
       surfaces { wingSurface { type patch; patches (wing); } } }
}"""
    # A COMPRESSIBLE CASE CARRYING THE INCOMPRESSIBLE DEFAULTS. The compressible_cp_frame
    # requirement would otherwise be exercised only where it does not apply, which is a
    # gate that passes by never firing (GEO-092). This fixture makes it fire.
    comp_bad = good.replace("functions {",
                            "application rhoSimpleFoam;\nfunctions {\n"
                            "  cp1 { type pressure; calcTotal no; calcCoeff yes;\n"
                            "        rho rho; rhoInf 1; pInf 0; }")
    comp_ok = good.replace("functions {",
                           "application rhoSimpleFoam;\nfunctions {\n"
                           "  cp1 { type pressure; calcTotal no; calcCoeff yes;\n"
                           "        rho rho; rhoInf 1.166399; pInf 100000; }")
    for name, txt, want in (("complete case", good, 8), ("bare case", bad, 2),
                            ("compressible, bad frame", comp_bad, 7),
                            ("compressible, good frame", comp_ok, 8)):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td) / "system"; d.mkdir(parents=True)
            (d / "controlDict").write_text(txt)
            got, miss = check(Path(td), "before")
            good_res = len(got) == want
            ok = ok and good_res
            print("  [%s] %-16s %d of %d satisfied (want %d)"
                  % ("ok " if good_res else "FAIL", name, len(got), len(REQUIREMENTS), want))
            if name == "compressible, bad frame":
                fired = "compressible_cp_frame" in {r["key"] for r in miss}
                ok = ok and fired
                print("  [%s] %-16s the cp-frame gate FIRES on pInf 0 with rhoInf 1"
                      % ("ok " if fired else "FAIL", ""))
            if name == "bare case":
                keys = {r["key"] for r in miss}
                need = "wing_surface_tau_VECTOR" in keys and "surface_format_carries_topology" in keys
                ok = ok and need
                print("  [%s] %-16s correctly flags the tau vector and the raw format"
                      % ("ok " if need else "FAIL", ""))
    print("\n  self-test: %s" % ("PASS, both directions" if ok else "*** FAILED ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
