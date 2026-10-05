#!/usr/bin/env python3
"""SHARED section conventions. One definition, called by every gate (GEO-084).

WHY THIS MODULE EXISTS, and it is the general finding of the Phase 0.5 delivery:
A CONVENTION THAT LIVES AS A FLAG ON ONE CALLER DOES NOT TRAVEL; A CONVENTION THAT
LIVES IN A SHARED MODULE CANNOT FAIL TO.

The blunt-TE chord reference was solved in D048 as `--te-ref midpoint` on
check_oml_gates.py. It was then RE-DERIVED WRONGLY in check_synthetic_gates.py
weeks later, where taking a TE corner instead of the base midpoint tilted the
chord line by atan(0.00105) = 0.0602 deg and read as twist on a wing built with
exactly zero twist. The knowledge existed, in the same repository, with a docstring
explaining it, and it could not travel because it was a flag rather than a
function. That is D058 with a second worked example, and chord_for(purpose) is the
precedent for the fix: put the rule where it cannot be bypassed.

THE RULE ITSELF. On a BLUNT trailing edge the two base corners sit at essentially
the same x, so `argmax(x)` picks one arbitrarily and tilts the chord line by up to
atan(t_base/c) about the true chord. On the ARGUS 0.63%c base that is +/- 0.180
deg against a 0.25 deg gate tolerance; on the synthetic 0.21%c base it is 0.0602
deg against a gate that wanted exact zero. THE MIDPOINT REFERENCE IS TE-CONVENTION
INSENSITIVE: on a blunt base it is the base midpoint, the standard chord
definition, and on a sharp edge it degenerates to the edge itself. ONE DEFINITION
THEREFORE SERVES BOTH SURFACES, which is required, because this programme holds a
blunt-TE surface and a sharp-TE reference simultaneously.
"""
import numpy as np

DEFAULT_TE_REF = "midpoint"
DEFAULT_BAND_FRAC = 0.005


def le_te_base(sec, te_ref=DEFAULT_TE_REF, band_frac=DEFAULT_BAND_FRAC):
    """Leading edge, trailing-edge reference point, and blunt-base thickness.

    sec: (n, 3) vertices of ONE section, x chordwise.

    te_ref='midpoint' (default) takes the centroid of the vertices within
    band_frac of the most-aft station. te_ref='maxx' takes the single most-aft
    vertex and is retained ONLY so that the two can be compared; it is not a
    defensible default on any surface in this programme.
    """
    sec = np.asarray(sec, dtype=float)
    x = sec[:, 0]
    le = sec[np.argmin(x)]
    chord_est = float(x.max() - x.min())
    band = sec[x >= x.max() - band_frac * chord_est]
    if te_ref == "maxx":
        te = sec[np.argmax(x)]
    elif te_ref == "midpoint":
        te = band.mean(axis=0)
    else:
        raise ValueError("te_ref must be 'maxx' or 'midpoint', got %r" % (te_ref,))
    if len(band) > 1:
        d = band[:, None, :] - band[None, :, :]
        base_thk = float(np.sqrt((d ** 2).sum(axis=2)).max())
    else:
        base_thk = 0.0
    return le, te, base_thk


def self_test():
    """Known answer: a symmetric blunt-TE section has ZERO twist by construction,
    and only the midpoint reference recovers it."""
    import math
    n = 200
    s = 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, n)))
    half_t = 0.00105                                   # 0.21%c blunt base, NACA 0010
    t = 5 * 0.10 * (0.2969 * np.sqrt(s) - 0.1260 * s - 0.3516 * s ** 2
                    + 0.2843 * s ** 3 - 0.1015 * s ** 4)
    xs = np.concatenate([s, s[::-1][:-1]])
    zs = np.concatenate([t, -t[::-1][:-1]])
    sec = np.column_stack([xs, np.zeros_like(xs), zs])

    ok = True
    for ref, want_zero in (("midpoint", True), ("maxx", False)):
        le, te, base = le_te_base(sec, te_ref=ref)
        twist = abs(math.degrees(math.atan2(te[2] - le[2], te[0] - le[0])))
        good = (twist < 1e-12) if want_zero else (twist > 1e-3)
        ok &= good
        print("  te_ref=%-9s twist %.6f deg  base %.6f c  %s"
              % (ref, twist, base, "PASS" if good else "FAIL"))
    print("  (maxx is EXPECTED to be non-zero: that is the defect being guarded, and")
    print("   atan(%.5f) = %.4f deg is exactly what it reads.)"
          % (half_t, math.degrees(math.atan(half_t))))
    print("\nSELF TEST: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(self_test())
