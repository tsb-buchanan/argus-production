#!/usr/bin/env python3
"""Prandtl lifting-line (Glauert Fourier series) for the Phase 0.5 synthetic wings.

GEO-082. THIS IS THE REFERENCE, NOT A POST-HOC EXPLANATION. The two synthetic
wings exist to give the 3D RANS route an answer that is known WITHOUT
MEASUREMENT, so the reference has to be computed here and recorded BEFORE any
solve, not fitted to one afterwards.

    ELLIPTICAL PLANFORM  span efficiency e = 1.000000 EXACTLY. For c(theta) =
                         c0 sin(theta) the monoplane equation is satisfied by A_1
                         alone, every higher A_n is identically zero, so
                         delta = 0 and e = 1. That is an ANALYTIC identity, not a
                         converged number, and it is this module's known-answer
                         gate (D056: a null whose answer is known independently).
    RECTANGULAR PLANFORM e < 1 STRICTLY. Physics fixes the SIGN, so the gate is a
                         sign test, not a tolerance (D060): a rectangular wing
                         measuring e at or above the elliptical wing's is a
                         defect, and no tolerance argument can rescue it.

METHOD. With mu(theta) = a0 c(theta) / (4 b) the monoplane equation is
    sum_n A_n sin(n theta) [ n mu(theta) + sin(theta) ] = mu(theta) (alpha - alpha_L0) sin(theta)
solved by collocation at interior theta. Symmetric loading uses ODD n only.
Then CL = pi AR A_1, CDi = pi AR sum_n n A_n^2, delta = sum_{n>=2} n (A_n/A_1)^2,
and e = 1 / (1 + delta).

SCOPE, stated so no one over-reads it. This is INVISCID, INCOMPRESSIBLE, and
assumes a straight lifting line with small-angle wake. It predicts INDUCED drag
only; it says nothing about profile drag, and the RANS CD must have its profile
part removed before comparison. It is the right reference for e and for nothing
else here.
"""
import argparse
import json
import math
import sys

import numpy as np

A0_THIN = 2.0 * math.pi          # thin-aerofoil lift-curve slope, per radian


def planform_area(chord_fn, span, n=20001):
    """Half-span-doubled planform area, integrated WITHOUT the tip singularity.

    A plain trapezoid rule in y is wrong for an ellipse: c(y) ~ sqrt(1 - (2y/b)^2)
    has an INFINITE DERIVATIVE at the tip, so the rule under-reads and the recovered
    AR came out 12.000001588 instead of 12. Substituting y = (b/2) cos(t) gives
        area = 2 * int_0^{pi/2} c((b/2) cos t) * (b/2) * sin t dt
    whose integrand is smooth for the ellipse (it becomes c0 sin^2 t), so the same
    rule is then exact to round-off. This is the D060 shape: the error had a KNOWN
    SIGN (a concave chord distribution can only be under-read), so the disagreement
    was evidence of a defect in the quadrature rather than scatter to be tolerated.

    Any y where chord_fn JUMPS must be passed as chord_fn.breaks, and the integral
    is split there: a trapezoid rule across a discontinuity smears it (D065, the
    edges-versus-centres rule, in its integration form).
    """
    breaks = sorted(set([0.0] + list(getattr(chord_fn, "breaks", [])) + [span / 2.0]))
    total = 0.0
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        if hi <= lo:
            continue
        t_hi = math.acos(min(1.0, max(-1.0, 2.0 * lo / span)))
        t_lo = math.acos(min(1.0, max(-1.0, 2.0 * hi / span)))
        t = np.linspace(t_lo, t_hi, n)
        yv = (span / 2.0) * np.cos(t)
        # Evaluate strictly inside the sub-interval so a jump at an endpoint cannot
        # contribute the wrong branch.
        cv = np.asarray([chord_fn(min(max(v, lo + 1e-12), hi - 1e-12)) for v in yv])
        total += 2.0 * float(np.trapz(cv * (span / 2.0) * np.sin(t), t))
    return total


def solve(chord_fn, span, n_terms=40, alpha_rad=math.radians(5.0),
          alpha_l0_rad=0.0, a0=A0_THIN):
    """Return dict with A_n, CL, CDi, delta, e for one planform.

    chord_fn: callable y (m, 0 at root, +span/2 at tip) -> chord (m).
    Uses ODD harmonics only, which is exact for a symmetric untwisted wing and
    avoids solving for coefficients that are identically zero.
    """
    ns = np.arange(1, 2 * n_terms, 2)                      # 1, 3, 5, ...
    # Collocation at interior points; theta = pi/2 is the root, theta -> 0 the tip.
    theta = np.linspace(0.0, math.pi / 2.0, n_terms + 2)[1:-1]
    y = (span / 2.0) * np.cos(theta)                       # theta=pi/2 -> y=0 (root)
    c = np.asarray([chord_fn(abs(yy)) for yy in y], dtype=float)
    if np.any(c <= 0):
        raise ValueError("chord must be strictly positive at every collocation point")
    mu = a0 * c / (4.0 * span)

    # alpha_l0_rad may be a SCALAR or a CALLABLE of |y|, which is what lets a
    # spanwise camber schedule (the aft-camber morph) be carried directly rather
    # than smeared into one equivalent angle. A scalar behaves exactly as before.
    if callable(alpha_l0_rad):
        al0 = np.asarray([alpha_l0_rad(abs(yy)) for yy in y], dtype=float)
    else:
        al0 = np.full(theta.shape, float(alpha_l0_rad))

    M = np.sin(np.outer(theta, ns)) * (ns[None, :] * mu[:, None] + np.sin(theta)[:, None])
    rhs = mu * (alpha_rad - al0) * np.sin(theta)
    A, *_ = np.linalg.lstsq(M, rhs, rcond=None)

    # AR from the SAME chord function that produced the solution, never passed in
    # separately: an AR from one planform against coefficients from another is the
    # frame-rule failure this project keeps hitting (project frame rule, instance 2).
    area = planform_area(chord_fn, span)
    ar = span ** 2 / area

    cl = math.pi * ar * A[0]
    cdi = math.pi * ar * float(np.sum(ns * A ** 2))
    delta = float(np.sum(ns[1:] * (A[1:] / A[0]) ** 2))
    return {"A": A.tolist(), "CL": cl, "CDi": cdi, "delta": delta,
            "e": 1.0 / (1.0 + delta), "AR": ar, "area_m2": area,
            "n_terms": n_terms, "alpha_deg": math.degrees(alpha_rad)}


def elliptical_chord(c0, span):
    def f(y):
        r = 2.0 * y / span
        return c0 * math.sqrt(max(0.0, 1.0 - r * r))
    return f


def truncated_elliptical_chord(c0, span, eta_tip):
    """Ellipse of root chord c0 cut off at |2y/b| = eta_tip.

    A TRUE ellipse has ZERO chord at the tip, which cannot be meshed and cannot
    even be written as a surface: the tip ring collapses to a point and its cap
    triangles are degenerate. So the as-built wing is truncated, and the cost of
    truncating is COMPUTED here rather than assumed negligible.
    """
    base = elliptical_chord(c0, span)

    def f(y):
        r = 2.0 * y / span
        return base(y) if r <= eta_tip else 0.0
    f.breaks = [eta_tip * span / 2.0]      # integrate up to the cut, not across it
    return f


def rectangular_chord(c):
    return lambda y: c


def self_test():
    """Known answers, in both directions."""
    span, ar = 3.6576, 12.0
    area = span ** 2 / ar
    c0 = 4.0 * area / (math.pi * span)
    ok = True

    ell = solve(elliptical_chord(c0, span), span)
    # e must be 1 to floating point, and the higher harmonics identically zero.
    e_err = abs(ell["e"] - 1.0)
    a_hi = max(abs(a) for a in ell["A"][1:])
    print("KNOWN-ANSWER GATE (analytic identity, not a converged number)")
    print("  elliptical e                 : %.12f  (|e-1| = %.2e)" % (ell["e"], e_err))
    print("  largest higher harmonic |A_n|: %.3e" % a_hi)
    print("  AR recovered from chord fn   : %.9f (want 12)" % ell["AR"])
    good = e_err < 1e-9 and a_hi < 1e-12 and abs(ell["AR"] - 12.0) < 1e-9
    print("  verdict                      : %s" % ("PASS" if good else "FAIL"))
    ok &= good

    rec = solve(rectangular_chord(area / span), span)
    print()
    print("SIGN TEST (physics fixes the direction, so no tolerance is argued)")
    print("  rectangular e                : %.6f" % rec["e"])
    print("  rectangular delta            : %.6f" % rec["delta"])
    print("  e_rect < e_ell               : %s"
          % ("PASS" if rec["e"] < ell["e"] else "*** FAIL ***"))
    ok &= rec["e"] < ell["e"]

    # Convergence: e must be stable in the number of terms, or the "known" answer
    # is a discretisation artefact.
    es = [solve(rectangular_chord(area / span), span, n_terms=n)["e"]
          for n in (10, 20, 40, 80)]
    spread = max(es) - min(es)
    print("  e_rect vs n_terms 10/20/40/80: %s" % " ".join("%.6f" % v for v in es))
    print("  spread                       : %.2e  %s"
          % (spread, "PASS" if spread < 1e-4 else "*** FAIL ***"))
    ok &= spread < 1e-4

    print("\nSELF TEST: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--planform", choices=["elliptical", "rectangular",
                                           "truncated_elliptical"])
    ap.add_argument("--span", type=float, default=3.6576)
    ap.add_argument("--ar", type=float, default=12.0)
    ap.add_argument("--eta-tip", type=float, default=0.98)
    args = ap.parse_args()
    if args.self_test or not args.planform:
        return self_test()

    area = args.span ** 2 / args.ar
    if args.planform == "elliptical":
        fn = elliptical_chord(4.0 * area / (math.pi * args.span), args.span)
    elif args.planform == "truncated_elliptical":
        fn = truncated_elliptical_chord(4.0 * area / (math.pi * args.span),
                                        args.span, args.eta_tip)
    else:
        fn = rectangular_chord(area / args.span)
    print(json.dumps(solve(fn, args.span), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
