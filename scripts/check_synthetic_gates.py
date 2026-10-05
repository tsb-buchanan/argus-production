#!/usr/bin/env python3
"""Hard gates on the Phase 0.5 SYNTHETIC wings, measured FROM THE SURFACE (GEO-082).

WHY THIS IS A SEPARATE SUITE AND NOT check_oml_gates.py. That suite's six gates
target the ARGUS wing's placement transform: dihedral 0.139 m rise, incidence
1.873 deg, LE sweep 28.86199 deg, tip twist -4.11495 deg. The synthetic wings
have NONE of those by construction, so every one of those gates fails on a
CORRECT synthetic surface. Inheriting them would be D077 in its plainest form:
a test valid in one configuration reused where the thing that made it true no
longer holds. --demo-argus-suite runs it anyway and shows it failing, so the
claim is demonstrated rather than asserted.

AND THE RE-DERIVED SUITE IS STRICTLY STRONGER HERE, which is the compensation.
Every target below is ANALYTIC, so there is no measurement uncertainty anywhere
in the gate: the wing either is the shape it was defined to be or it is not.

GATES:
  1. SEMI-SPAN         tip station at exactly b/2 = 1.8288 m.
  2. PLANFORM AREA     half-wing b^2/(2 AR) = 0.557418240 m^2. A SIGN TEST, not
                       a tolerance (D060): chordwise cosine sampling hits the
                       exact LE and TE, and spanwise trapezoid across a CONCAVE
                       elliptical chord law can only UNDER-read. A measured area
                       ABOVE the analytic value is therefore a defect, not
                       scatter, and no tolerance can rescue it.
  3. CHORD LAW         measured chord at every ring against c(y), max abs error.
  4. ZERO SWEEP        every ring's LE at x = 0.
  5. ZERO DIHEDRAL     every ring's mid-thickness line at z = 0. THIS IS THE
                       INVERTED ARGUS GATE 1, and it is deliberately marked
                       WEAK: zero is also what a silently-dropped transform
                       produces, so it cannot distinguish "correctly planar"
                       from "transform lost". It is gate 3 and gate 8 that carry
                       the load here, because they have non-trivial targets.
  6. ZERO TWIST        section chord line parallel to x at every ring.
  7. SECTION t/c       10.00% at every ring, and CONSTANT along the span.
  8. NACA 0010 SHAPE   RMS departure of the measured half-thickness from the
                       analytic polynomial, at every ring, in unit-chord
                       coordinates. Non-trivial target, so this is the gate that
                       actually proves the aerofoil is what it claims.
  9. ROOT OPEN         no facets at y = 0 (symmetry plane, taken from the
                       background mesh), matching the ARGUS convention.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "geometry" / "scripts"))
from check_oml_gates import read_stl_vertices                      # noqa: E402
from vsp3_to_stl import cosine                                     # noqa: E402
from section_conventions import le_te_base                         # noqa: E402
from make_synthetic_wing import (SPAN_M, ASPECT_RATIO, ETA_TIP,    # noqa: E402
                                 naca_thickness, cut_ellipse_fn, elliptical_c0)

# Supremum of 2*naca_thickness, in per cent, found on a dense grid. NOT 10.00.
NACA_TC_MAX_PCT = float(200.0 * naca_thickness(np.linspace(0.0, 1.0, 4000001)).max())

# Resolution floor of the ASCII STL itself. write_stl formats vertices as %.9e, so
# each coordinate carries ~1e-9 relative error; t/c = 2 z / c divides two such
# numbers and is quoted in per cent, so the floor is ~1e-8 %pt. Gating tighter than
# the file format is gating on the writer, not the geometry: the elliptical wing,
# whose chord differs at every station, read 5.2e-9 against a guessed 1e-9.
STL_TC_FLOOR_PCT = 2e-8

TOL = {"semi_span_m": 1e-9, "chord_m": 5e-7, "le_x_m": 1e-9, "dihedral_z_m": 1e-9,
       "twist_deg": 1e-9, "tc_pct": 1e-4, "naca_rms_c": 1e-9, "area_rel": 2e-4}


def rings_from_vertices(v, tol=1e-9):
    """Group vertices into spanwise rings by their y value."""
    ys = np.unique(np.round(v[:, 1] / tol) * tol)
    return [(y, v[np.abs(v[:, 1] - y) < 1e-9]) for y in ys]


def check(path, kind, nchord=200):
    # The t/c gate compares against the polynomial's maximum ON THE ACTUAL SAMPLING
    # GRID, not against its supremum and not against a tolerance. A cosine grid does
    # not land on x = 0.3, so a correct surface reads 3.3e-4 %pt low; that is a
    # COMPUTABLE quantity, so it is computed and the gate becomes exact.
    tc_on_grid = float(200.0 * naca_thickness(cosine(nchord)).max())
    v = read_stl_vertices(path)
    # THE ROOT OVERHANG IS NOT PART OF THE AERODYNAMIC SURFACE (GEO-101 I4, GEO-102).
    # Since the closure repair, every file carries a skirt and a cap below y = 0 that
    # exist only so snappy sees a solid. They lie OUTSIDE the half-model domain and
    # contribute no wetted area, so every gate below must measure the y >= 0 surface
    # and nothing else. Without this filter the cap projects onto the planform and
    # gate 2 over-reads by 1.38%, which is the overhang being counted as wing.
    n_all = len(v)
    tri = v.reshape(-1, 3, 3)
    keep = (tri[:, :, 1] >= -1e-9).all(axis=1)
    v = tri[keep].reshape(-1, 3)
    n_dropped = n_all - len(v)
    rings = rings_from_vertices(v)
    span, ar = SPAN_M, ASPECT_RATIO
    if kind == "elliptical":
        c0 = elliptical_c0(span, ar, ETA_TIP)
        chord_of = cut_ellipse_fn(c0, span, ETA_TIP)
    else:
        c0 = (span ** 2 / ar) / span
        chord_of = lambda y: c0                                    # noqa: E731

    res = []
    ys = np.array([y for y, _ in rings])

    # 1. semi-span
    got = float(ys.max())
    res.append(("1 semi-span (m)", got, span / 2.0, abs(got - span / 2.0),
                TOL["semi_span_m"], abs(got - span / 2.0) <= TOL["semi_span_m"]))

    chord_err, le_err, z_err, twist_err, tc_err, naca_rms = [], [], [], [], [], []
    widths, centres, n_interior = [], [], []
    # CANARY, NOT AN EXCLUSION (GEO-086). These two conditions used to `continue`
    # silently, so a degenerate ring dropped out of gates 3 to 8 and the suite still
    # reported PASS over the rings that remained. A known-bad input CARVED OUT of an
    # audit is a hole in it; the same input RETAINED AS A REQUIRED ABSENCE tests the
    # detector on every run. Neither can fire on a correct wing, so the count is
    # gated at zero rather than the rings being quietly dropped.
    skipped_thin, skipped_zero_chord = [], []
    for y, pts in rings:
        pts = np.unique(np.round(pts, 12), axis=0)
        if pts.shape[0] < 4:
            skipped_thin.append((float(y), int(pts.shape[0])))
            continue
        x, z = pts[:, 0], pts[:, 2]
        c_meas = float(x.max() - x.min())
        chord_err.append(abs(c_meas - chord_of(y)))
        le_err.append(abs(float(x.min())))
        widths.append(c_meas)
        centres.append(y)
        if c_meas <= 0:
            skipped_zero_chord.append(float(y))
            continue
        xs, zt = (x - x.min()) / c_meas, z / c_meas

        # THE FLAT TIP CAP CONTRIBUTES ONE VERTEX THAT IS NOT ON THE AEROFOIL: the
        # ring centroid, at x/c 0.5 and z/c ~ 3e-6. Comparing it against the NACA
        # polynomial gave a 2.36e-02 c "shape error" that was entirely the gate's.
        # Identify interior vertices structurally and ASSERT THE COUNT rather than
        # filtering silently: exactly one on a capped ring, none anywhere else.
        body = (xs > 0.02) & (xs < 0.98)
        interior = body & (np.abs(zt) < 0.5 * naca_thickness(np.clip(xs, 0, 1)))
        n_interior.append((y, int(interior.sum())))
        keep = ~interior

        # LE vertex and TE-face midpoint. THE TE REFERENCE MUST BE THE MIDPOINT:
        # this aerofoil has a 0.21%c blunt TE, so taking either TE corner tilts the
        # chord line by atan(0.00105) = 0.0602 deg and reads as twist. That is the
        # SAME trap check_oml_gates.py already solved with --te-ref midpoint (D048);
        # it did not travel, which is D058 exactly.
        # SHARED CONVENTION, not a local re-derivation (GEO-084). This is where the
        # 0.0602 deg phantom twist came from the first time.
        outline = pts[keep]
        le_p, te_p, _base = le_te_base(outline)
        z_err.append(max(abs(float(le_p[2])), abs(float(te_p[2]))))
        twist_err.append(abs(math.degrees(math.atan2(te_p[2] - le_p[2],
                                                     te_p[0] - le_p[0]))))

        # t/c: the target is NOT 10.00%. The NACA 4-digit polynomial's actual
        # maximum is 10.00246%, because the "10" names the nominal thickness the
        # polynomial approximates, not its supremum. Gating on 10.00 flagged a
        # correct aerofoil at 2.5e-3 %pt. Target = the polynomial's own maximum,
        # and sampling a smooth maximum on a discrete grid can only UNDER-read it,
        # so the sign test runs against THAT (D060).
        tc_err.append(2.0 * float(zt[keep].max()) * 100.0 - tc_on_grid)

        up = keep & (zt > 0)
        naca_rms.append(float(np.sqrt(np.mean(
            (zt[up] - naca_thickness(np.clip(xs[up], 0, 1))) ** 2))))

    def bundle(n, label, arr, tol):
        got = float(np.max(arr))
        res.append((("%d %s" % (n, label)), got, 0.0, got, tol, got <= tol))

    bundle(3, "chord law, max |err| (m)", chord_err, TOL["chord_m"])
    bundle(4, "zero sweep, max |x_LE| (m)", le_err, TOL["le_x_m"])
    bundle(5, "zero dihedral, max |z| LE/TE-mid (m) [WEAK]", z_err, TOL["dihedral_z_m"])
    bundle(6, "zero twist, max |deg| (TE MIDPOINT ref)", twist_err, TOL["twist_deg"])

    tc_err = np.asarray(tc_err)
    res.append(("7 section t/c (%%) vs grid max %.6f" % tc_on_grid, float(tc_err.min()),
                0.0, float(np.abs(tc_err).max()), STL_TC_FLOOR_PCT,
                float(np.abs(tc_err).max()) <= STL_TC_FLOOR_PCT))
    res.append(("7b t/c vs polynomial supremum %.5f, SIGN (under-read only)"
                % NACA_TC_MAX_PCT, tc_on_grid - NACA_TC_MAX_PCT, 0.0,
                tc_on_grid - NACA_TC_MAX_PCT, 0.0,
                tc_on_grid - NACA_TC_MAX_PCT <= 0.0))
    res.append(("7c t/c spanwise spread (%)", float(tc_err.max() - tc_err.min()), 0.0,
                float(tc_err.max() - tc_err.min()), TOL["tc_pct"],
                float(tc_err.max() - tc_err.min()) <= TOL["tc_pct"]))
    bundle(8, "NACA 0010 shape, max RMS (c)", naca_rms, TOL["naca_rms_c"])

    capped = [(y, n) for y, n in n_interior if n]
    ok_cap = (len(capped) == 1 and abs(capped[0][0] - span / 2.0) < 1e-9
              and capped[0][1] == 1)
    res.append(("8b cap vertices off-aerofoil (structural: exactly 1, at the tip)",
                sum(n for _y, n in n_interior), 1, len(capped), 0, ok_cap))

    # CANARY GATES. Both must be zero on any wing this suite can adjudicate; a
    # non-zero count means gates 3 to 8 silently covered fewer rings than the file
    # contains, so the suite's own coverage is the thing being asserted here.
    res.append(("0a rings dropped as degenerate (<4 vertices) [CANARY]",
                len(skipped_thin), 0, len(skipped_thin), 0, not skipped_thin))
    res.append(("0b rings dropped as zero-chord [CANARY]",
                len(skipped_zero_chord), 0, len(skipped_zero_chord), 0,
                not skipped_zero_chord))
    res.append(("0c rings actually adjudicated by gates 3-8",
                len(naca_rms), len(rings), len(rings) - len(naca_rms), 0,
                len(naca_rms) == len(rings)))

    # 2. planform area: SIGN TEST
    order = np.argsort(centres)
    yy, ww = np.asarray(centres)[order], np.asarray(widths)[order]
    area_half = float(np.trapz(ww, yy))
    area_ref = span ** 2 / (2.0 * ar)
    rel = (area_half - area_ref) / area_ref
    # SIGN FLOOR. The under-read argument applies to a CONCAVE chord law, where the
    # trapezoid rule cuts corners. For the CONSTANT-chord wing the rule is EXACT, so
    # only round-off remains and the residual can fall either way at the 1e-16 level.
    # A sign test with no floor fails on one ulp, which is the test being wrong, not
    # the geometry.
    sign_floor = 1e-12
    res.insert(1, ("2 planform area, half (m2)", area_half, area_ref, rel,
                   TOL["area_rel"], rel <= sign_floor and abs(rel) <= TOL["area_rel"]))
    res.append(("2b area SIGN (<= 0 for a concave chord law; round-off floor 1e-12)",
                rel, 0.0, rel, sign_floor, rel <= sign_floor))

    # 9. root open
    # GATE 9 IS INVERTED SINCE GEO-101. The root used to be required OPEN, because
    # write_stl left it open for the symmetry plane. That premise was the defect: an
    # open root is what let snappy flood the aerofoil interior. The root is now
    # required CLOSED BY AN OVERHANG, and the evidence is that triangles exist below
    # y = 0 which this gate deliberately filtered out above.
    res.append(("9 root CLOSED: overhang triangles below y=0 (was: root open)",
                n_dropped, 0, n_dropped, 0, n_dropped > 0))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stl", required=True)
    ap.add_argument("--kind", choices=["elliptical", "rectangular"], required=True)
    ap.add_argument("--nchord", type=int, default=200,
                    help="chordwise points per surface used at build time")
    ap.add_argument("--demo-argus-suite", action="store_true",
                    help="also run the ARGUS gates, which MUST fail here (D077 demo)")
    args = ap.parse_args()

    res = check(args.stl, args.kind, args.nchord)
    print("SYNTHETIC WING GATES (GEO-082): %s [%s]" % (Path(args.stl).name, args.kind))
    print("  %-52s %14s %14s %10s  %s" % ("gate", "measured", "target", "err", ""))
    npass = 0
    for name, got, target, err, tol, ok in res:
        npass += bool(ok)
        print("  %-52s %14.9g %14.9g %10.2e  %s"
              % (name, got, target, err, "PASS" if ok else "*** FAIL ***"))
    print("\n  %d/%d gates PASS" % (npass, len(res)))

    if args.demo_argus_suite:
        import subprocess
        print("\nD077 DEMONSTRATION: the ARGUS suite on a CORRECT synthetic wing.")
        print("It targets dihedral 0.139 m, incidence 1.873 deg, sweep 28.86199 deg")
        print("and tip twist -4.11495 deg, none of which this wing has BY DESIGN.")
        r = subprocess.run([sys.executable, str(REPO / "scripts" / "check_oml_gates.py"),
                            args.stl], capture_output=True, text=True)
        tail = (r.stdout + r.stderr).strip().splitlines()
        for line in tail[-14:]:
            print("    " + line)
        overall = [l for l in (r.stdout + r.stderr).splitlines() if "OVERALL" in l]
        print("    -> %s ; exit status %d. FAIL IS THE EXPECTED RESULT HERE."
              % (overall[0].strip() if overall else "no OVERALL line", r.returncode))
    return 0 if npass == len(res) else 1


if __name__ == "__main__":
    sys.exit(main())
