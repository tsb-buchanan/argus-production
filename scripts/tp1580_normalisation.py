#!/usr/bin/env python3
"""tp1580_normalisation.py: TP-1580 <-> DSO coefficient renormalisation, derived not inherited.

WHY THIS IS A TESTED FUNCTION AND NOT AN INLINE FACTOR. TP-1580 and the DSO package
normalise on DIFFERENT reference areas for the SAME wing, and the ratio is 1.1131, so a
raw comparison reads 11.31% high in CL. That is large enough to look like a physical
disagreement and small enough to be believed. The factor is therefore derived here from
first principles, asserted against the report's published value, and applied only through
these functions.

THE DERIVATION, so the rescaling is defensible rather than asserted.
TP-1580 page 2 states the AR-12 reference area as 1.11 m^2 (12 ft^2) and the reference
span as 3.66 m (12.0 ft), and states that "These reference wing areas are based on the
trapezoidal planform which extends from the model center line to the wing tip".
Reconstructing that trapezoid from Zheng's delivered section table, by extending the
outboard panel (chord 1.116670 ft at y 2.297998 ft, the eta 0.383 crank, to chord
0.500000 ft at y 5.999996 ft, the tip) linearly inboard to the centre line:
    taper slope      (0.500000 - 1.116670) / (5.999996 - 2.297998) = -0.1665777 ft/ft
    root chord       1.116670 + 0.1665777 * 2.297998                =  1.499465 ft
    trapezoid area   (1.499465 + 0.500000)/2 * 12.0                 = 11.996793 ft^2
    trapezoid AR     12.0^2 / 11.996793                             = 12.0032
11.996793 against the published 12.0 ft^2 is 0.027%, which CONFIRMS the construction and
explains the "AR 12" designation. The actual WETTED planform is 13.3572 ft^2 at AR 10.781.
The two are different references for the same wing, not a conflict.

A D036 TRAP LIVES EXACTLY HERE. The area ratio 13.3572/12.0 = 1.113100 (DIMENSIONLESS)
and Zheng's mean geometric chord S/b = 13.3572/12.0 = 1.113100 ft (a LENGTH) are the SAME
NUMERAL, because bref = 12.0 ft and Sref_TP = 12.0 ft^2 share the numeral 12. Any code or
prose landing on 1.1131 must state which one it means.

REFERENCE CHORDS ARE A SEPARATE TRAP, and the report is explicit about it: the AR-12 wing
has a mean geometric chord of 33.02 cm (13.00 in.) and the AR-10 wing 34.90 cm, but "an
approximate average value of 34.04 cm (13.40 in.) was used as the reference mean geometric
chord for BOTH". So TP-1580's cref is NOT the AR-12 mean geometric chord. Moment
coefficients carry both the area and the chord, so Cm needs both factors.

Usage:
  python3 scripts/tp1580_normalisation.py --selftest
"""

import argparse

FT2_M2 = 0.09290304
FT_M = 0.3048

# --- as PUBLISHED, not converted (item 5 rule: store as published) ----------------------
TP1580 = {
    "AR12": {"Sref_ft2": 12.0, "bref_ft": 12.0, "Sref_m2": 1.11, "bref_m": 3.66},
    "AR10": {"Sref_ft2": 11.21, "bref_ft": 10.59, "Sref_m2": 1.04, "bref_m": 3.23},
    "cref_in": 13.40, "cref_cm": 34.04,
    "cref_note": "approximate average used for BOTH AR-10 and AR-12; equals the local "
                 "wing chord at the trailing-edge break station. NOT the AR-12 mean "
                 "geometric chord, which the report gives as 33.02 cm (13.00 in.).",
    "mgc_AR12_cm": 33.02, "mgc_AR12_in": 13.00, "mgc_AR10_cm": 34.90,
    "moment_centre": "model plane of symmetry, 50.14 cm (19.74 in.) longitudinally aft of "
                     "the wing leading edge and 6.60 cm (2.60 in.) vertically below the "
                     "wing reference plane; corresponds to the quarter-chord of the mean "
                     "geometric chord of the AR-12 wing",
}
DSO = {"Sref_ft2": 13.3572, "Cref_ft": 1.29251, "bref_ft": 12.0}

# outboard panel as delivered (refined_section_table.csv), for the trapezoid derivation
PANEL = {"y_crank_ft": 2.297998, "c_crank_ft": 1.116670,
         "y_tip_ft": 5.999996, "c_tip_ft": 0.500000}


def trapezoid_from_outboard_panel(p=PANEL, bref_ft=12.0):
    """Extend the outboard panel linearly to the centre line. Returns root chord, area, AR."""
    slope = (p["c_tip_ft"] - p["c_crank_ft"]) / (p["y_tip_ft"] - p["y_crank_ft"])
    c_root = p["c_crank_ft"] - slope * p["y_crank_ft"]
    area = 0.5 * (c_root + p["c_tip_ft"]) * bref_ft
    return c_root, area, bref_ft ** 2 / area


def area_ratio():
    """S_wetted / S_TP1580(AR12). DIMENSIONLESS. See the D036 note in the docstring."""
    return DSO["Sref_ft2"] / TP1580["AR12"]["Sref_ft2"]


def cl_tp1580_to_dso(cl):
    """A CL published on TP-1580's 12 ft^2 expressed on the DSO 13.3572 ft^2 basis.
    DIVIDES by the area ratio, so the number gets SMALLER by 11.31%."""
    return cl / area_ratio()


def cl_dso_to_tp1580(cl):
    return cl * area_ratio()


def cm_tp1580_to_dso(cm):
    """Cm carries q*S*c, so BOTH the area and the reference chord change."""
    f = (TP1580["AR12"]["Sref_ft2"] * TP1580["cref_in"] / 12.0) / \
        (DSO["Sref_ft2"] * DSO["Cref_ft"])
    return cm * f


def cm_dso_to_tp1580(cm):
    f = (TP1580["AR12"]["Sref_ft2"] * TP1580["cref_in"] / 12.0) / \
        (DSO["Sref_ft2"] * DSO["Cref_ft"])
    return cm / f


# --- Condition WT Reynolds definition (D042 item 3) -------------------------------------
WT = {"U_m_s": 40.8, "nu_m2_s": 1.46e-5, "Re_stated": 0.95e6,
      "chord_named_in_definition_m": 0.3404,
      "definition": "ARGUS_reference_data.md section 9 item 2, VERBATIM: 'Condition WT "
                    "(tunnel-matched validation): M 0.12, U = 40.8 m/s, Re_cref = 0.95e6 "
                    "on cref = 0.3404 m.' The chord is NAMED at the point of definition, "
                    "so it is resolved, not inferred from the number."}
DSO_CREF_M = 0.393957


def re_on(chord_m, U=WT["U_m_s"], nu=WT["nu_m2_s"]):
    return U * chord_m / nu


def selftest():
    ok = True

    def chk(name, got, exp, tol):
        nonlocal ok
        good = abs(got - exp) <= tol
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] {name:46s} {got:.6f} vs {exp:.6f} "
              f"(tol {tol:g})")

    c_root, area, ar = trapezoid_from_outboard_panel()
    chk("trapezoidal root chord, ft", c_root, 1.499465, 1e-5)
    chk("trapezoidal area, ft^2", area, 11.996793, 1e-5)
    chk("trapezoidal aspect ratio", ar, 12.0032, 1e-3)
    # THE ASSERTION THAT MATTERS: the reconstruction must land on the PUBLISHED 12.0 ft^2
    rel = abs(area - TP1580["AR12"]["Sref_ft2"]) / TP1580["AR12"]["Sref_ft2"]
    chk("reconstruction vs published Sref, rel err", rel, 0.0, 5e-4)
    chk("area ratio (dimensionless)", area_ratio(), 1.113100, 1e-6)
    chk("CL overstatement if used raw, percent", (area_ratio() - 1) * 100, 11.31, 5e-3)
    chk("round trip CL", cl_dso_to_tp1580(cl_tp1580_to_dso(1.23)), 1.23, 1e-12)
    chk("round trip Cm", cm_dso_to_tp1580(cm_tp1580_to_dso(-0.0965)), -0.0965, 1e-12)
    # published unit conversions, as a check on the stored PUBLISHED values
    chk("AR12 Sref ft^2 -> m^2 vs published 1.11",
        TP1580["AR12"]["Sref_ft2"] * FT2_M2, 1.11, 5e-3)
    chk("cref 13.40 in -> m", TP1580["cref_in"] * 0.0254, 0.340360, 1e-6)
    chk("AR12 mean geometric chord 13.00 in -> m",
        TP1580["mgc_AR12_in"] * 0.0254, 0.330200, 1e-6)
    # the D036 numeral collision, asserted so it cannot be forgotten
    chk("D036 collision: area ratio == S/b in ft, same numeral",
        area_ratio(), DSO["Sref_ft2"] / DSO["bref_ft"], 1e-12)
    # --- the Re_cref naming collision, asserted so it cannot drift back ----------------
    chk("Condition WT Re on the NAMED chord 0.3404 m",
        re_on(WT["chord_named_in_definition_m"]) / 1e6, 0.9512, 1e-3)
    chk("Condition WT Re on the DSO Cref 0.393957 m instead",
        re_on(DSO_CREF_M) / 1e6, 1.1009, 1e-3)
    gap = (re_on(DSO_CREF_M) / re_on(WT["chord_named_in_definition_m"]) - 1) * 100
    chk("the two differ by (percent) - the field is NAMED Re_cref and the "
        "number is NOT on Cref", gap, 15.74, 0.05)
    # the two CANDIDATE mean-geometric chords are NOT separable at this precision
    chk("TP-1580 blended cref vs DSO S/b, percent apart",
        abs(0.3404 - 0.339273) / 0.339273 * 100, 0.332, 5e-3)
    # D042 item 3c: never re-derive Re from Mach with OUR nu
    a_ms = 340.0
    re_derived = re_on(WT["chord_named_in_definition_m"], U=0.118 * a_ms) / 1e6
    chk("Re re-derived from the tabulated M 0.118 with OUR nu (must NOT be used)",
        re_derived, 0.9354, 5e-3)
    chk("...against the report's stated 0.97e6, percent LOW",
        (1 - re_derived / 0.97) * 100, 3.57, 0.05)
    print(f"\nSELFTEST: {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        raise SystemExit(0 if selftest() else 1)
    c_root, area, ar = trapezoid_from_outboard_panel()
    print(f"trapezoid: root {c_root:.6f} ft, area {area:.6f} ft^2, AR {ar:.4f}")
    print(f"area ratio {area_ratio():.6f}; a raw TP-1580 CL reads "
          f"{(area_ratio()-1)*100:.2f}% high on the DSO basis")


if __name__ == "__main__":
    main()
