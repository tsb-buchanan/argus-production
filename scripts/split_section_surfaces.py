#!/usr/bin/env python3
"""Add the upper/lower surface flag to the raw band CSVs written on HPC12.

section_slices_from_surface.py deliberately stops at the RAW BAND and says why: HPC12 runs
numpy 1.12, which has no `argsort(kind="stable")`, so the shared splitter cannot run there,
and D056 forbids writing a fourth implementation to work around that. The split therefore
happens here, on the laptop, with airfoil_surfaces.split_branches unmodified.

RECOVERING INDICES WITHOUT REIMPLEMENTING THE ASSIGNMENT RULE
-------------------------------------------------------------
split_branches returns COORDINATES, not indices, because it was written for geometry point
sets. Attaching a flag to a DATA row needs indices. Replaying its documented rule here to
get them would be exactly the fourth implementation D056 forbids, and would be free to drift
from the validated one.

So the membership is RECOVERED BY A MERGE WALK over the returned branches instead. Both
branches are subsequences of the x-sorted point list, in increasing order, so stepping
through the sorted points and consuming whichever branch matches at (x, z) reconstructs the
partition without knowing how it was decided. The recovery is then ASSERTED against the
arrays split_branches actually returned: if the reconstructed upper does not equal the
returned upper element for element, this HALTS. That makes the recovery a within-file proof
(D052 item 4) rather than an assumption about the walk's internals.

THE CLEANLINESS TEST, AND WHY ITS FIRST VERSION WAS WRONG (D068, caught 2026-09-15)
-----------------------------------------------------------------------------------
My first criterion was min(z_upper) > max(z_lower) bin by bin. It fired on 16 of 55 files,
always in the SINGLE bin x/c 0.875-0.900, and always at the outboard stations. THE GATE WAS
WRONG, NOT THE DATA, and the reason is this project's own standing failure mode.

A cut is a BAND, not a plane: it pools every face within +-0.004 semispan in y, which at
eta 0.95 is 0.084 of the LOCAL chord. So min(z_upper) and max(z_lower) are generally
measured at OPPOSITE EDGES OF THE BAND, i.e. at different spanwise stations, and differencing
them compares two quantities in different frames. Toward the trailing edge, where the section
thins, the spanwise z variation across the band exceeds the local half-thickness and the
TAILS overlap while the surfaces themselves never touch.

MEASURED, because the claim needs a number rather than an argument. At a FIXED x/c +-0.0015,
where a 2D section has exactly one point per surface and therefore ZERO z spread, the
within-branch spread is 0.002 to 0.012 c. That spread can only be spanwise. Against it the
branch MEDIAN separation is 0.096 to 0.133 c at x/c 0.30 and 0.016 to 0.021 c at x/c 0.88:
never negative, and never smaller than the spread.

SO THE TEST COMPARES LIKE WITH LIKE: median(z_upper) > median(z_lower) per bin, which is
insensitive to which band edge a tail came from. The zero margin is kept, because the
surfaces crossing is a one-directional defect and a sign test beats a tolerance (D060).

A SECOND, MORE USEFUL NUMBER IS REPORTED BESIDE IT: whether the branch separation exceeds the
within-branch spanwise spread. Where it does not, the band pools more span than the local
thickness can resolve, and THAT is a real statement about the data rather than about the
splitter. It is reported per file rather than used to fail one, because the cp and cf values
are still each attached to their own correct surface.

DO NOT "FIX" A FIRING GATE BY WIDENING ITS TOLERANCE. The first version was replaced because
its ESTIMATOR pooled frames, not because its threshold was tight.
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from airfoil_surfaces import split_branches

XMIN, XMAX, NBIN = 0.10, 0.90, 32


def recover_membership(x, z):
    """-> boolean array, True = upper. Halts if the recovery disagrees with split_branches."""
    br = split_branches(x, z)
    if br is None:
        return None, "too-few-points"
    (xu, zu), (xl, zl) = br

    o = np.argsort(np.asarray(x, float), kind="stable")
    xs, zs = np.asarray(x, float)[o], np.asarray(z, float)[o]

    is_up = np.zeros(len(xs), dtype=bool)
    iu = il = 0
    got_u, got_l = [], []
    for k in range(len(xs)):
        mu = iu < len(xu) and xu[iu] == xs[k] and zu[iu] == zs[k]
        ml = il < len(xl) and xl[il] == xs[k] and zl[il] == zs[k]
        if mu:
            is_up[k] = True
            got_u.append(k)
            iu += 1
        if ml:
            got_l.append(k)
            il += 1
        if not mu and not ml:
            return None, "merge-walk-desync at sorted index %d" % k
    if iu != len(xu) or il != len(xl):
        return None, "merge-walk consumed %d/%d upper, %d/%d lower" % (iu, len(xu), il, len(xl))

    # The recovery is PROVED against what the splitter returned, not trusted.
    if not (np.array_equal(xs[got_u], xu) and np.array_equal(zs[got_u], zu)):
        return None, "recovered upper does not match split_branches"
    if not (np.array_equal(xs[got_l], xl) and np.array_equal(zs[got_l], zl)):
        return None, "recovered lower does not match split_branches"

    # Points at the LE belong to both branches; the flag has to pick one, and `upper`
    # wins. Report how many so the ambiguity is counted rather than hidden.
    shared = len(set(got_u) & set(got_l))
    out = np.zeros(len(x), dtype=bool)
    out[o] = is_up
    return out, "shared=%d" % shared


def cleanliness(x, z, is_up):
    """Bin-by-bin branch ordering over the mid-chord.

    -> (nbad, nbins_tested, worst_median_gap, nbins_unresolved)

    `nbad` counts bins where the branch MEDIANS are out of order, which is a genuine split
    failure. `nbins_unresolved` counts bins where the median separation is smaller than the
    larger within-branch spanwise spread: there the branches are correctly labelled but the
    band pools more span than the local thickness resolves. The two are counted separately
    because they mean different things and only the first is a defect.
    """
    edges = np.linspace(XMIN, XMAX, NBIN + 1)
    nbad = ntest = nunres = 0
    worst = None
    for i in range(NBIN):
        m = (x >= edges[i]) & (x < edges[i + 1])
        zu, zl = z[m & is_up], z[m & ~is_up]
        if len(zu) < 3 or len(zl) < 3:
            continue
        ntest += 1
        gap = float(np.median(zu) - np.median(zl))
        if worst is None or gap < worst:
            worst = gap
        if gap <= 0:
            nbad += 1
        spread = max(zu.max() - zu.min(), zl.max() - zl.min())
        if gap <= spread:
            nunres += 1
    return nbad, ntest, worst, nunres


def process(path):
    head, rows = [], []
    with open(path) as fh:
        for r in csv.reader(fh):
            if r and r[0].startswith("#"):
                head.append(r)
            elif r and r[0] == "x_over_c":
                cols = r
            elif r:
                rows.append(r)
    # REFUSE A SECOND PASS rather than appending a second flag column and a second
    # provenance header. Re-running is the obvious thing to do after editing this file,
    # and silently doubling the columns is the obvious way for that to go wrong.
    if "surface" in cols:
        return {"file": path.name, "status": "ALREADY-SPLIT",
                "detail": "carries a 'surface' column; re-transfer the raw band to redo"}

    a = np.array(rows, dtype=float)
    x, z = a[:, 0], a[:, 1]

    is_up, note = recover_membership(x, z)
    if is_up is None:
        return {"file": path.name, "status": "SPLIT-FAILED", "detail": note}

    nbad, ntest, worst, nunres = cleanliness(x, z, is_up)
    status = "ok" if nbad == 0 and ntest >= 20 else (
        "BRANCHES-OUT-OF-ORDER" if nbad else "too-few-testable-bins")

    # THE FLAG IS WRITTEN EVEN WHEN THE VERDICT IS NOT "ok", AND THE VERDICT GOES IN THE
    # HEADER. An earlier version of this file carried the comment "promote only after the
    # checks pass" above an UNCONDITIONAL replace(), which is the vouched-success shape
    # (GEO-080): a claim authored beside the action rather than derived from it. Either
    # gate, or say plainly that you do not. This says plainly.
    tmp = path.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh)
        for h in head:
            if not h[0].startswith("# surface_split"):
                w.writerow(h)
        w.writerow(["# surface_split",
                    "airfoil_surfaces.split_branches (index recovery proved against its "
                    "returned branches); LE %s; verdict %s; branch-median ordering clean in "
                    "%d/%d bins over x/c %.2f-%.2f; %d of those bins have median separation "
                    "below the within-branch spanwise spread of the band"
                    % (note, status, ntest - nbad, ntest, XMIN, XMAX, nunres)])
        w.writerow(cols + ["surface"])
        for i in range(len(a)):
            w.writerow(rows[i] + ["upper" if is_up[i] else "lower"])
    tmp.replace(path)

    return {"file": path.name, "status": status, "n": len(a),
            "n_upper": int(is_up.sum()), "n_lower": int((~is_up).sum()),
            "bins_bad": nbad, "bins_tested": ntest, "worst_gap": worst,
            "bins_unresolved": nunres, "detail": note}


def main():
    d = Path(sys.argv[1])
    files = sorted(d.glob("*_eta*.csv"))
    if not files:
        sys.exit("no section CSVs under %s" % d)
    res = [process(f) for f in files]
    bad = [r for r in res if r["status"] != "ok"]
    for r in res:
        if r["status"] == "ok":
            print("ok    %-26s n=%6d  up=%6d lo=%6d  bins %2d/%2d  min median gap %+.5f  "
                  "band-unresolved bins %d"
                  % (r["file"], r["n"], r["n_upper"], r["n_lower"],
                     r["bins_tested"] - r["bins_bad"], r["bins_tested"], r["worst_gap"],
                     r["bins_unresolved"]))
        else:
            print("BAD   %-26s %s  %s" % (r["file"], r["status"], r.get("detail", "")))
    nun = sum(r.get("bins_unresolved", 0) for r in res)
    print("\nPARTITION: %d files = %d clean + %d not-clean" % (len(res), len(res) - len(bad), len(bad)))
    print("band-unresolved bins (correctly labelled, but band pools more span than the "
          "local thickness resolves): %d across all files" % nun)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
