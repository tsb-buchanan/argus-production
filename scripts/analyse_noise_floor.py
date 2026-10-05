#!/usr/bin/env python3
"""analyse_noise_floor.py: A1 repeatability verdict against the D047 pre-committed bands.

REPORTS THE NUMBER AND THE BAND, NOTHING ELSE. The bands were written to the decision log
BEFORE the test ran (D047 item 6) and are reproduced here verbatim so the comparison
cannot drift:

    UNDER 1 COUNT   -> snappy CARRIES the candidate matrix; the architecture stands.
    1 TO 2 COUNTS   -> MARGINAL. snappy usable ONLY with the measured noise carried as an
                       explicit uncertainty band on EVERY reported delta.
    OVER 2 COUNTS   -> re-meshing candidates independently is UNSAFE; the whole matrix
                       goes structured.
    A result of 1.9 counts is MARGINAL and not "essentially under 2".

CONVERGENCE IS D023 FORCE STATIONARITY, NOT RESIDUALS. The window mean drift must be
under 0.01 counts over the last N iterations; residuals are meaningless here for the
reason measured in mesh_requirements.md section 5. Both cases must pass INDEPENDENTLY
before their difference means anything: differencing two unconverged runs measures the
transient, not the mesh.

Usage:
  python3 scripts/analyse_noise_floor.py [--window 1000]
"""

import argparse
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
CASES = {"A1_ref": "cases/noise_floor/A1_ref", "A1_pert": "cases/noise_floor/A1_pert"}
DRIFT_GATE = 0.01          # counts, D023

# --- REFERENCE QUANTITIES, taken from the SAME values forceCoeffs uses ------------------
# D052: an earlier version hard-coded AR = 12.0, which is the TRAPEZOIDAL-reference
# aspect ratio belonging to Sref 11.9968 ft^2. Our coefficients are normalised on the
# WETTED area (Aref 0.620462 m2 = half of 1.240924), whose AR is 10.7807. The two differ
# by 1.1134, which is the 1.1131 area ratio yet again (D036, fourth occurrence). Mixing a
# CL from one reference area with an AR from the other silently inflates the induced term
# by 11%. The assertion below makes that impossible rather than merely documented.
AREF_HALF = 0.620462       # m2, EXACTLY the forceCoeffs Aref
BREF = 3.6576              # m, full span
LREF = 0.39396             # m, DSO Cref, for cross-checking the dictionary


def aspect_ratio(aref_half=None):
    """AR = Bref^2 / Sref on the SAME Sref the coefficients use. Never hard-code this.

    aref_half is READ FROM THE CASE when available, so the AR and the coefficients cannot
    come from different reference areas even if the dictionary is edited."""
    a = AREF_HALF if aref_half is None else aref_half
    return BREF ** 2 / (2.0 * a)


def assert_consistent_reference(case: Path):
    """RAISE if the case's forceCoeffs Aref/lRef differ from what this script assumes.

    The failure this guards against is not a typo, it is a coefficient set and an aspect
    ratio that belong to DIFFERENT reference areas while both looking entirely plausible.
    e > 1 on a planar wing is the physical tell, and it only shows up if someone thinks
    to check; this raises instead."""
    f = sorted((case / "postProcessing/forceCoeffs").glob("*/forceCoeffs.dat"))
    if not f:
        return None
    aref_seen = []
    for ln in f[0].read_text().splitlines():
        if ln.startswith("# Aref"):
            v = float(ln.split(":")[1])
            # tolerance is 0.1% RELATIVE: it must catch reference-area MIXING, which is
            # 11.3%, and must not trip on dictionary rounding, which is 3e-6 here.
            if abs(v - AREF_HALF) / AREF_HALF > 1e-3:
                raise AssertionError(
                    f"{case.name}: forceCoeffs Aref {v} differs from {AREF_HALF} by "
                    f"{abs(v-AREF_HALF)/AREF_HALF*100:.2f}%. AR and the coefficients "
                    f"would belong to DIFFERENT reference areas (D052).")
            aref_seen.append(v)
        if ln.startswith("# lRef"):
            v = float(ln.split(":")[1])
            if abs(v - LREF) / LREF > 1e-3:
                raise AssertionError(f"{case.name}: forceCoeffs lRef {v} != {LREF}")
    return aref_seen[0] if aref_seen else None


def read_forces(case: Path):
    f = sorted((case / "postProcessing/forceCoeffs").glob("*/forceCoeffs.dat"))
    t, cd, cl, cm = [], [], [], []
    for p in f:
        for ln in p.read_text().splitlines():
            if ln.startswith("#"):
                continue
            v = ln.split()
            if len(v) < 4:
                continue
            t.append(float(v[0])); cm.append(float(v[1]))
            cd.append(float(v[2])); cl.append(float(v[3]))
    o = np.argsort(t)
    return (np.array(t)[o], np.array(cd)[o], np.array(cl)[o], np.array(cm)[o])


def stationarity(t, y, window):
    """D023: window-MEAN drift between the two halves of the last `window` iterations."""
    m = t >= t.max() - window
    yy = y[m]
    if len(yy) < 20:
        return None, None, None
    h = len(yy) // 2
    d = (yy[h:].mean() - yy[:h].mean()) * 1e4          # counts
    return d, yy.mean(), yy.std() * 1e4


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--window", type=float, default=1000)
    a = ap.parse_args()

    res = {}
    print("A1 NOISE FLOOR, against the D047 pre-committed bands\n")
    ok = True
    for nm, rel in CASES.items():
        aref_case = assert_consistent_reference(REPO / rel)
        t, cd, cl, cm = read_forces(REPO / rel)
        d, mean_cd, sd = stationarity(t, cd, a.window)
        _, mean_cl, _ = stationarity(t, cl, a.window)
        _, mean_cm, _ = stationarity(t, cm, a.window)
        conv = d is not None and abs(d) <= DRIFT_GATE
        ok &= bool(conv)
        res[nm] = (mean_cd, mean_cl, mean_cm)
        print(f"{nm}: {len(t)} samples, t {t.min():.0f}-{t.max():.0f}")
        print(f"   cd window mean {mean_cd:.8f} ({mean_cd*1e4:.3f} counts), "
              f"sd {sd:.3f} counts")
        print(f"   D023 drift {d:+.4f} counts vs gate {DRIFT_GATE}  -> "
              f"{'CONVERGED' if conv else 'NOT CONVERGED'}")
        print(f"   cl {mean_cl:.6f}   cm {mean_cm:.6f}")
    if not ok:
        print("\nAT LEAST ONE CASE IS NOT STATIONARY. No verdict is issued: differencing "
              "two unconverged runs measures the transient, not the mesh.")
        return 1

    dcd = (res["A1_pert"][0] - res["A1_ref"][0]) * 1e4
    dcl = res["A1_pert"][1] - res["A1_ref"][1]
    dcm = res["A1_pert"][2] - res["A1_ref"][2]
    print(f"\n{'='*62}\nNOISE FLOOR  |dCD| = {abs(dcd):.3f} COUNTS")
    print(f"   (dCD {dcd:+.3f} counts, dCL {dcl:+.6f}, dCm {dcm:+.6f})")
    band = ("UNDER 1 COUNT: snappy CARRIES the candidate matrix" if abs(dcd) < 1.0
            else "1 TO 2 COUNTS: MARGINAL, noise must be carried as an explicit band "
                 "on every reported delta" if abs(dcd) <= 2.0
            else "OVER 2 COUNTS: re-meshing candidates independently is UNSAFE; "
                 "the matrix goes structured")
    print(f"   BAND: {band}\n{'='*62}")
    # THE FIXED-ALPHA NUMBER IS NOT THE WHOLE STORY, and this is reported alongside
    # rather than instead: the two meshes give slightly different CL at the same alpha,
    # and at fixed alpha a CL difference produces a CD difference through INDUCED drag.
    # The production candidate comparison is at MATCHED CL, so the operative noise is
    # the part that survives after that term is removed.
    cl = res["A1_ref"][1]
    AR = aspect_ratio(aref_case)
    print(f"   AR = Bref^2/Sref = {BREF}^2/{2*AREF_HALF:.6f} = {AR:.4f} "
          f"(NOT 12.0032, which is the trapezoidal reference)")
    for e in (0.85, 0.90, 0.95):
        dcdi_dcl = 2 * cl / (np.pi * AR * e)
        corr = dcd - dcdi_dcl * dcl * 1e4
        print(f"   at MATCHED CL (e={e:.2f}, dCDi/dCL={dcdi_dcl:.5f}): "
              f"residual {corr:+.3f} counts")
    print("   Report the fixed-alpha number against the bands; the matched-CL residual")
    print("   is context for what the candidate comparison will actually see.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
