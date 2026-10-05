#!/usr/bin/env python3
"""wing_outlines_of12.py: smooth aerofoil outlines cut from the PLACED STL, one per station.

    python3 scripts/wing_outlines_of12.py --all
    python3 scripts/wing_outlines_of12.py --self-test

WHY THIS EXISTS. Review, 2026-09-16: the appendix aerofoil plots made the leading and
trailing edges wiggle on a smooth geometry, a false impression.

He is right, and the cause is that the outline was never geometry. plot_section_slices.py drew
it from the CFD BAND SAMPLE: the mean z/c of the faces inside a +/-0.004 semispan cut. That
band is 7.32 mm wide regardless of station, so on a swept, tapered wing it pools faces whose
true z genuinely differs, and the pooled mean wiggles wherever the surface curves fastest.
MEASURED, same treatment both sides (60 bins, median, second difference of z/c):

    eta    source            LE 0-5%   mid 20-80%   TE 95-100%
    0.20   CFD band         1.47e-03     3.89e-05     2.81e-04
    0.60   CFD band         8.90e-03     4.35e-05     4.60e-03
    0.95   CFD band         1.08e-02     5.33e-05     5.01e-03
    0.20   STL planar cut   3.50e-04     7.93e-05     1.75e-04
    0.60   STL planar cut   4.72e-04     6.76e-05     1.98e-04
    0.95   STL planar cut   1.82e-04     8.74e-05     8.47e-05

The band is 4x rougher at the root leading edge, 19x at mid-span and 59x at eta 0.95. The
DIAGNOSTIC detail is not the ratio but the trend: the band degrades SEVENFOLD outboard while
the cut is flat across span, which is what a fixed-width band sampling a shrinking chord must
do and what a zero-width plane cannot do.

WHERE THE CUT IS WORSE, STATED RATHER THAN OMITTED: mid-chord, 7.9e-05 against the band's
3.9e-05, because a plane crosses ~1,600 triangles where the band pools 60,114 faces, so the
per-bin median is noisier. Mid-chord is not where the wiggle was visible and 7.9e-05 of chord
is 0.02 mm at the tip, but a fix that is better in three places and worse in one should say so.

THE CUT IS A PLANE, NOT A SLAB. Every returned point lies exactly on y = eta * semi, obtained
by linear interpolation along each crossed triangle edge. There is no width, so there is no
pooling, so there is nothing to average.

THE UPPER/LOWER SPLIT IS THE SHARED airfoil_surfaces.split_branches, IMPORTED AND NOT
REIMPLEMENTED (D056 records three implementations tried and rejected; the sign-of-z one
shipped briefly and was wrong on exactly these aft-cambered sections). My first version of the
smoothness test DID write a throwaway splitter, and it reported the planar cut as 3000x rougher
than the band, which is how I know it matters: the measurement convicted my own shortcut.

FRAME. The STL is the PLACED surface in metres, the same object snappyHexMesh was given, so
its x/c is the CFD's x/c with sweep, twist and dihedral already in it. The unit-chord
ordinates in the .vsp3 are NOT usable here: they are untwisted and unplaced, and a smooth
section in the wrong frame would sit under the Cp curves looking authoritative and aligned
with nothing.
"""
import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from airfoil_surfaces import split_branches  # noqa: E402

STL = REPO / "geometry/derived/wing3d"
OUT = REPO / "results/wing_outlines_of12"
SEMI = 1.8288                      # bref/2, the registered value, common to every case
ETAS = (0.20, 0.40, 0.60, 0.80, 0.95)

# Verified against each case's own CASE_PROVENANCE.json and the registry sha256, not from
# any prose document: two of those disagree (see scripts/render_concept_3d.py).
GEOM = {"B": "baseline", "C": "cte_i002_c04", "F": "cfft_b02_c01",
        "H": "chc_g02_c06", "M": "mcv2_i002_c01", "W": "cffw_b01_c01"}


def load_stl(name):
    p = STL / ("%s_oml_placed_laddercap.stl" % name)
    if not p.exists():
        raise SystemExit("%s: no placed STL at %s" % (name, p))
    txt = p.read_text(errors="replace")
    V = np.array(re.findall(r"vertex\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)", txt),
                 dtype=float).reshape(-1, 3, 3)
    return V, p


def planar_cut(V, y0):
    """Exact triangle-plane intersection at y = y0. Returns the intersection points."""
    s = V[:, :, 1] - y0
    hit = ~((s > 0).all(axis=1) | (s < 0).all(axis=1))
    pts = []
    for tri, sg in zip(V[hit], s[hit]):
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if (sg[a] < 0) != (sg[b] < 0):
                t = sg[a] / (sg[a] - sg[b])
                pts.append(tri[a] + t * (tri[b] - tri[a]))
    return np.array(pts)


def outline(V, eta):
    """-> dict with upper/lower x/c and z/c, in the SAME frame the section CSVs use."""
    P = planar_cut(V, eta * SEMI)
    if len(P) < 100:
        raise SystemExit("eta %.2f: only %d intersection points; refusing" % (eta, len(P)))
    x, z = P[:, 0], P[:, 2]
    chord = float(x.max() - x.min())
    if chord <= 0:
        raise SystemExit("eta %.2f: degenerate chord" % eta)
    # x/c from the section's own extent and z/c about its own mean z: EXACTLY the convention
    # section_slices_from_surface.py uses, so the outline and the Cp curves share one frame.
    #
    # BUT THAT FRAME IS PER-GEOMETRY, AND SIX OF THEM GO ON ONE AXIS, so the frame must
    # travel with the numbers (D068). x_min_m and z_mean_m are recorded for exactly that.
    #
    # WHY IT IS NOT OPTIONAL, measured. The six wings are identical forward of the hinge, so
    # their outlines MUST coincide there; that is a null with an independently known answer
    # (D058). Taking z about each section's OWN mean fails it, because aft camber moves that
    # mean and a purely aft change then leaks forward as a rigid offset: 0.000e+00 at eta
    # 0.20 and 0.40 where the morph is zero, but 0.43, 1.01 and 0.44 per cent of chord at
    # eta 0.60, 0.80 and 0.95 where it is not. At eta 0.80 that spurious forward offset is
    # 1.01 %c against a REAL aft difference of 2.36 %c: 43 per cent of the signal the figure
    # exists to show, and the worst pair is M against W, the largest camber command against
    # no aft camber at all. A consumer re-datums onto ONE geometry's frame using these two
    # fields; z_abs = z_over_c * chord_m + z_mean_m recovers the placed coordinate exactly.
    xc, zc = (x - x.min()) / chord, (z - z.mean()) / chord
    br = split_branches(xc, zc)
    if br is None:
        raise SystemExit("eta %.2f: split_branches returned None; refusing to guess a split" % eta)
    (ux, uz), (lx, lz) = br
    ou, ol = np.argsort(ux), np.argsort(lx)
    return dict(eta=eta, chord_m=chord, x_min_m=float(x.min()), z_mean_m=float(z.mean()),
                n_points=int(len(P)),
                upper_xc=np.asarray(ux)[ou].tolist(), upper_zc=np.asarray(uz)[ou].tolist(),
                lower_xc=np.asarray(lx)[ol].tolist(), lower_zc=np.asarray(lz)[ol].tolist())


def binned_rough(xc, zc, lo, hi, nb=60):
    e = np.linspace(lo, hi, nb + 1)
    zb = np.array([np.median(zc[(xc >= e[k]) & (xc < e[k + 1])])
                   if ((xc >= e[k]) & (xc < e[k + 1])).sum() else np.nan for k in range(nb)])
    g = ~np.isnan(zb)
    return float(np.nanmean(np.abs(np.diff(zb[g], 2)))) if g.sum() > 3 else float("nan")


def self_test():
    """KNOWN ANSWER (D058): the cut must be smoother than the band at LE and TE, and its
    roughness must NOT grow outboard. Both are checked, because the second is the one that
    identifies the mechanism rather than merely the symptom."""
    V, _ = load_stl(GEOM["B"])
    bad = 0
    le = []
    for eta in (0.20, 0.60, 0.95):
        o = outline(V, eta)
        ux, uz = np.array(o["upper_xc"]), np.array(o["upper_zc"])
        r_le = binned_rough(ux, uz, 0.0, 0.05)
        r_te = binned_rough(ux, uz, 0.95, 1.0)
        le.append(r_le)
        ok = r_le < 1.0e-03 and r_te < 1.0e-03
        print("  %s eta %.2f  LE %.2e  TE %.2e  (band was up to 1.08e-02 / 5.01e-03)"
              % ("PASS" if ok else "FAIL", eta, r_le, r_te))
        bad += 0 if ok else 1
    grow = le[-1] / le[0] if le[0] else float("inf")
    ok = grow < 3.0
    print("  %s LE roughness outboard/root ratio %.2f (band's was 7.3; a plane has no width "
          "so it must not grow)" % ("PASS" if ok else "FAIL", grow))
    bad += 0 if ok else 1
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not a.all:
        return ap.error("give --all or --self-test")
    a.out.mkdir(parents=True, exist_ok=True)
    ok, bad = 0, []
    for letter, name in sorted(GEOM.items()):
        t0 = time.time()
        try:
            V, p = load_stl(name)
            secs = [outline(V, e) for e in ETAS]
        except SystemExit as e:
            print("  %-2s REFUSED: %s" % (letter, e)); bad.append(letter); continue
        rec = dict(letter=letter, geometry=name, stl=str(p.relative_to(REPO)),
                   stl_bytes=p.stat().st_size, semispan_m=SEMI,
                   method="exact triangle-plane intersection; no band, no smoothing",
                   split="airfoil_surfaces.split_branches, imported unmodified",
                   sections=secs, seconds=round(time.time() - t0, 1))
        (a.out / ("outline_%s.json" % letter)).write_text(json.dumps(rec) + "\n")
        print("  %-2s %-16s %d stations, %d..%d pts, %.0f s"
              % (letter, name, len(secs), min(s["n_points"] for s in secs),
                 max(s["n_points"] for s in secs), rec["seconds"]))
        ok += 1
    print("  PARTITION: %d requested = %d written + %d refused%s"
          % (len(GEOM), ok, len(bad), (": " + ", ".join(bad)) if bad else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
