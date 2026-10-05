#!/usr/bin/env python3
"""THE splitter. Separating an aerofoil section into upper and lower surfaces.

Import this. Do not write a fourth implementation.

    from airfoil_surfaces import split_branches, surface_at

WHY THIS EXISTS AS A SHARED UTILITY (D056, D058)
------------------------------------------------
On an AFT-CAMBERED SUPERCRITICAL SECTION the obvious test is wrong. The EET
sections carry enough aft camber that the LOWER surface rises ABOVE the chord
line and stays above to the trailing edge: at defining station 4, 32 of its 87
stored points have z > 0, from x/c 0.78 onward. Anything that separates the
surfaces by the SIGN of z therefore mislabels a third of the lower surface as
upper, and corrupts the comparison over exactly x/c 0.78 to 1.00.

That band is where aft-camber morphing does its work, so THE ARTEFACT AND THE
SIGNAL OCCUPY THE SAME PLACE. The failure does not look like a bug; it looks
like a result.

THREE IMPLEMENTATIONS WERE TRIED AND REJECTED, each failing in a way that
produced a plausible number rather than an error:

  SIGN OF z
      Wrong for this section family, as above. This is the one that shipped
      briefly and inflated a reported surface residual roughly threefold.

  MAX/MIN z PER CHORDWISE BIN
      Density-dependent. Correct on a dense slice (1592 points) and wrong on
      the sparse stored ordinates (151 points), where most bins hold a single
      point and the max is simply whichever branch that point came from.
      Rejected by a null case at 5.6e-02.

  ANGULAR SORT ABOUT THE CENTROID
      An aerofoil is too thin near the trailing edge for a star-shaped sort:
      points on opposite surfaces get near-identical angles. Misassigned 30 of
      151 points and produced a non-monotonic branch. Rejected by the same null
      case at 1.6e-02.

  ADOPTED: walk the points in x order and assign each to whichever branch its z
      is nearer. This tracks CONTINUITY and assumes nothing about camber,
      sign, or sampling density.
      VALIDATED against stored branch membership at all 19 defining stations of
      the delivered baseline: ZERO misassignments, every branch monotonic in x.

RELATED PRIOR ART IN THIS REPO, and read this before assuming the problem is
new: scripts/collect_trim_deltas.py solved the same trap independently for the
2D programme, assigning wall faces to a surface by PROXIMITY TO THE MESHED .dat
SURFACES and noting in its docstring that "a z-sign test fails on drooped-TE
morphed sections". That script is deliberately NOT refactored to use this
module: it underpins D019, D021 and D022, and changing validated code to adopt
a new utility would mean re-validating signed results for no gain. The two
approaches are equivalent in intent; use that one for wall-face data that has a
reference surface, and this one for section point sets that do not.
"""
import numpy as np


def split_branches(xr, zr):
    """Split a section point set into (upper, lower) branches.

    xr, zr are chord-frame coordinates (any scaling; only ordering matters).
    Returns ((x_up, z_up), (x_lo, z_lo)), or None if there are too few points.

    Assigns each point, in x order, to whichever branch its z is nearer to.
    Near the leading edge the two branches share points; that is harmless
    because the branches coincide there.
    """
    xr = np.asarray(xr, float)
    zr = np.asarray(zr, float)
    if len(xr) < 10:
        return None
    o = np.argsort(xr, kind="stable")
    xs, zs = xr[o], zr[o]
    up_i, lo_i = [0], [0]
    for k in range(1, len(xs)):
        du = abs(zs[k] - zs[up_i[-1]])
        dl = abs(zs[k] - zs[lo_i[-1]])
        (up_i if du <= dl else lo_i).append(k)
    upper = (xs[up_i], zs[up_i])
    lower = (xs[lo_i], zs[lo_i])
    if upper[1].mean() < lower[1].mean():
        upper, lower = lower, upper
    return upper, lower


def surface_at(xr, zr, xq, upper):
    """Evaluate one surface at the abscissae xq. NaN outside the branch range.

    Returns NaN rather than extrapolating: a point outside the branch is
    UNMEASURED, not zero, and the caller must be able to tell the difference.
    """
    br = split_branches(xr, zr)
    if br is None:
        return None
    xs, zs = br[0] if upper else br[1]
    o = np.argsort(xs)
    xs, zs = xs[o], zs[o]
    xu, idx = np.unique(xs, return_index=True)
    zu = zs[idx]
    if len(xu) < 5:
        return None
    xq = np.asarray(xq, float)
    inb = (xq >= xu.min()) & (xq <= xu.max())
    out = np.full(len(xq), np.nan)
    out[inb] = np.interp(xq[inb], xu, zu)
    return out


def self_test(vsp3="dso_reference/baseline_wing_only_refined.vsp3"):
    """Re-run the validation: zero misassignments at all 19 defining stations.

    Kept in the module so the claim in the docstring is checkable rather than
    asserted. Run: python3 scripts/airfoil_surfaces.py
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from vsp3_sections import extract_wing

    secs = extract_wing(vsp3)
    secs = secs["sections"] if isinstance(secs, dict) else secs
    total_mis, n_stations, non_monotonic = 0, 0, 0
    for sec in secs:
        if sec.get("upper") is None:
            continue
        up = np.asarray(sec["upper"], float)
        lo = np.asarray(sec["lower"], float)
        xr = np.concatenate([up[:, 0], lo[:, 0]])
        zr = np.concatenate([up[:, 1], lo[:, 1]])
        U, _ = split_branches(xr, zr)
        stored_upper = set(zip(np.round(up[:, 0], 9), np.round(up[:, 1], 9)))
        recovered = set(zip(np.round(U[0], 9), np.round(U[1], 9)))
        total_mis += len(recovered - stored_upper)
        non_monotonic += 0 if np.all(np.diff(U[0]) >= -1e-12) else 1
        n_stations += 1
    print("stations checked      : %d" % n_stations)
    print("misassigned points    : %d" % total_mis)
    print("non-monotonic branches: %d" % non_monotonic)
    ok = total_mis == 0 and non_monotonic == 0 and n_stations >= 19
    print("SELF TEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(self_test())
