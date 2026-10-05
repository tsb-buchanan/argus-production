#!/usr/bin/env python3
"""audit_span_efficiency.py: re-derive the VLM-vs-RANS span-efficiency claim from source.

MEETING, 2026-08-13: the room doubted e = 0.55 for the VLM, because a span LOADING that agrees
with ours to +/-0.025 cannot also carry 65% more induced drag than ours. They were right.

WHAT WAS WRONG. VSPAERO's .polar carries TWO induced-drag sets and we quoted the wrong one:

  SURFACE INTEGRATION (near field)   CDi,  E   -- pressure integrated on the panels
  WAKE INDUCED        (far field)    CDiw, Ew  -- Trefftz plane, i.e. the wake

They are both in the same row of the same file, 22 columns apart, and the header names both
blocks. Near-field induced drag in a vortex-lattice method is the unreliable one: it depends on
resolving leading-edge suction, which a panel method does poorly. The far-field value is the
trustworthy one and is what any RANS Trefftz-plane number must be compared against.

OUR RANS NUMBER IS A TREFFTZ-PLANE INTEGRAL. So the like-for-like partner is CDiw / Ew, and we
used CDi / E. That is the frame rule (D068) with a new instance: two quantities with the same
NAME and the same UNITS, computed by the same code, in the same file, differing by a factor of
1.74 because one is near-field and the other far-field.

THE THIRD ROUTE IS THE ARBITER. Neither of VSPAERO's own columns is allowed to adjudicate
between them, so this script also computes the induced drag INDEPENDENTLY from the .lod
spanwise loading by classical lifting-line (Fourier series in the Trefftz plane). That number
knows nothing about either column and says which one the loading actually supports.
"""
import argparse
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from parse_vspaero_lod import parse_lod  # noqa: E402

DSO = REPO / "dso_reference/2026-07-29_correction_pair"
CASES = ("baseline_corrected", "mcv2_i002_c01", "cte_i002_c04")
CL_TRIM = 0.428277635108

# RANS Trefftz-plane induced drag, ARG-103 item 1. Counts, at each case's OWN CL.
RANS = {"baseline_corrected": (67.583, 0.455051),
        "mcv2_i002_c01": (78.717, 0.492932)}


def polar(tag):
    L = [l for l in (DSO / tag / ("%s.polar" % tag)).read_text().splitlines() if l.strip()]
    return dict(zip(L[2].split(), (float(x) for x in L[3].split())))


def vspaero_ref(tag):
    d = {}
    for ln in (DSO / tag / ("%s.vspaero" % tag)).read_text().splitlines():
        if "=" in ln:
            k, _, v = ln.partition("=")
            try:
                d[k.strip()] = float(v)
            except ValueError:
                pass
    return d


def lifting_line_e(tag):
    """e from the SPANWISE LOADING alone, by classical lifting line. Independent of the polar.

    y = -(b/2) cos(theta); the circulation is expanded as Gamma = 2 b V sum A_n sin(n theta),
    which gives CL = pi AR A_1 and CDi = pi AR sum n A_n^2, hence e = A_1^2 / sum n A_n^2.

    THE HALF-WING FACTOR OF TWO (D068 instance 1) is handled by MIRRORING the strips rather
    than by scaling: the .lod carries the starboard half only, and a Fourier fit needs the
    whole span. The CL closure below is what proves the mirroring is right.
    """
    h, c, A = parse_lod(DSO / tag / ("%s.lod" % tag))
    iy, ic, icl, ida = c.index("Yavg"), c.index("Chord"), c.index("Cl"), c.index("dArea")
    b, Sref = h["Bref"], h["Sref"]
    y, ch, cl = A[:, iy], A[:, ic], A[:, icl]
    k = np.argsort(y)
    y, ch, cl = y[k], ch[k], cl[k]
    CL_closure = 2.0 * np.sum(A[:, icl] * A[:, ida]) / Sref

    # mirror to the full span, then fit on theta
    yf = np.concatenate([-y[::-1], y])
    gf = np.concatenate([(ch * cl)[::-1], ch * cl])          # loading  c * cl
    th = np.arccos(np.clip(-2.0 * yf / b, -1.0, 1.0))
    o = np.argsort(th)
    th, gf = th[o], gf[o]

    # Gamma / (2 b V) = c*cl / (4 b)  =  sum A_n sin(n theta)
    g = gf / (4.0 * b)
    N = 24
    M = np.stack([np.sin(n * th) for n in range(1, N + 1)], axis=1)
    An, *_ = np.linalg.lstsq(M, g, rcond=None)
    AR = b * b / Sref
    CL_ll = np.pi * AR * An[0]
    n = np.arange(1, N + 1)
    CDi_ll = np.pi * AR * float(np.sum(n * An ** 2))
    e_ll = An[0] ** 2 / float(np.sum(n * An ** 2))
    return CL_closure, CL_ll, CDi_ll, e_ll, AR


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.parse_args()

    print("=" * 96)
    print("1. THE TWO INDUCED-DRAG SETS IN VSPAERO'S OWN .polar, side by side")
    print("=" * 96)
    print("  %-20s %10s %10s %10s %10s %10s" %
          ("case", "CDi ct", "E", "CDiw ct", "Ew", "AR"))
    P, R = {}, {}
    for t in CASES:
        p, ref = polar(t), vspaero_ref(t)
        AR = ref["Bref"] ** 2 / ref["Sref"]
        P[t], R[t] = p, AR
        print("  %-20s %10.3f %10.4f %10.3f %10.4f %10.4f"
              % (t, p["CDi"] * 1e4, p["E"], p["CDiw"] * 1e4, p["Ew"], AR))

    print("\n  KNOWN-ANSWER CHECK: reproduce VSPAERO's own E and Ew from CLtot, AR and each CDi.")
    bad = 0
    for t in CASES:
        p, AR = P[t], R[t]
        for col, cd in (("E", "CDi"), ("Ew", "CDiw")):
            got = p["CLtot"] ** 2 / (np.pi * AR * p[cd])
            ok = abs(got - p[col]) < 5e-4
            bad += not ok
            print("    %-20s %-3s recomputed %.4f against %.4f  %s"
                  % (t, col, got, p[col], "OK" if ok else "*** MISMATCH ***"))
    if bad:
        print("    the recomputation disagrees with the file: STOP, the reading is wrong")
        return 1
    print("    Both reproduce exactly, so E and Ew differ ONLY in which CDi they use.")

    print("\n" + "=" * 96)
    print("2. THE INDEPENDENT ARBITER: induced drag from the LOADING alone (lifting line)")
    print("=" * 96)
    for t in CASES:
        clc, cll, cdi, e, AR = lifting_line_e(t)
        p = P[t]
        near, far = p["CDi"] * 1e4, p["CDiw"] * 1e4
        pick = "FAR FIELD (CDiw)" if abs(cdi * 1e4 - far) < abs(cdi * 1e4 - near) else \
               "NEAR FIELD (CDi)"
        print("  %-20s CL closure %.6f (target %.6f)" % (t, clc, CL_TRIM))
        print("      lifting-line CDi %8.3f ct,  e = %.4f" % (cdi * 1e4, e))
        print("      near field %8.3f ct   far field %8.3f ct   -> the loading supports %s"
              % (near, far, pick))

    print("\n" + "=" * 96)
    print("3. WHAT THE DELTA BECOMES. The DSO's predicted benefit, both ways.")
    print("=" * 96)
    b = P["baseline_corrected"]
    print("  %-20s %14s %14s" % ("candidate", "near dCDi ct", "far dCDi ct"))
    for t in ("mcv2_i002_c01", "cte_i002_c04"):
        print("  %-20s %+14.3f %+14.3f"
              % (t, (P[t]["CDi"] - b["CDi"]) * 1e4, (P[t]["CDiw"] - b["CDiw"]) * 1e4))
    print("\n  The -6.171 (mcv2) and -6.445 (cte) this project has quoted for two weeks are the")
    print("  NEAR-FIELD column. The far-field numbers from the SAME RUNS are ~25x smaller.")

    print("\n" + "=" * 96)
    print("4. THE LIKE-FOR-LIKE COMPARISON. Our RANS number is a Trefftz-plane integral.")
    print("=" * 96)
    kb = RANS["baseline_corrected"][0] / RANS["baseline_corrected"][1] ** 2
    km = RANS["mcv2_i002_c01"][0] / RANS["mcv2_i002_c01"][1] ** 2
    AR = R["baseline_corrected"]
    for nm, k in (("baseline", kb), ("morphed", km)):
        print("  RANS %-9s CDi/CL^2 %8.3f ct -> at CL %.6f: %7.3f ct, e = %.4f"
              % (nm, k, CL_TRIM, k * CL_TRIM ** 2, 1.0 / (np.pi * AR * k * 1e-4)))
    print("  VLM  baseline  far field                          %7.3f ct, e = %.4f"
          % (b["CDiw"] * 1e4, b["Ew"]))
    print("  VLM  morphed   far field                          %7.3f ct, e = %.4f"
          % (P["mcv2_i002_c01"]["CDiw"] * 1e4, P["mcv2_i002_c01"]["Ew"]))
    print("\n  RANS dCDi (ARG-103, Trefftz)                     %+7.3f ct"
          % (km * CL_TRIM ** 2 - kb * CL_TRIM ** 2))
    print("  VLM  dCDi far field, like for like               %+7.3f ct"
          % ((P["mcv2_i002_c01"]["CDiw"] - b["CDiw"]) * 1e4))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
