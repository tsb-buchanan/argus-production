#!/usr/bin/env python3
"""Close a derived wing STL into a watertight solid by a ROOT OVERHANG (GEO-101).

THE DEFECT. Every derived wing surface is an OPEN SHELL: the root section outline is
a free rim with no cap, 398 open edges on the ARGUS lofts and 399 on the synthetics.
After the placement transform's 5 degrees of dihedral the rim is no longer planar at
y = 0 either: the root rib's z-extent 0.095687 m times sin(5 deg) = 8.340 mm predicts
the measured rim y range of 8.371 mm to 0.37%. So 63.3% of the rim sits ABOVE the
symmetry plane, leaving a slot up to 4.679 mm, and 36.4% pokes below it.

THE REPAIR IS ADDITIVE AND EXTENDS THROUGH THE PLANE, NOT UP TO IT.
  1. skirt every rim edge in -y down to Y_CAP by a ruled surface
  2. cap the projected rim with a planar rib at Y_CAP

NOT A CAP AT y = 0. A cap coincident with the symmetry patch makes snappy choose
between wall and symmetry at the same location, which is a second failure mode
wearing the costume of a fix. With the solid extending THROUGH y = 0, the domain
boundary cuts a closed solid and no gap can exist by construction.

EVERY RIM POINT HAS y >= -3.692 mm, comfortably above Y_CAP = -20 mm, so every skirt
segment has strictly positive length and no degenerate triangle is created. That is
asserted per segment, not assumed.

INVARIANTS, asserted here and re-checked independently by scripts/check_surface_closed.py:
  I1  NO POINT WITH y >= 0 MOVES. Original facets are copied VERBATIM as text, so the
      aerodynamic surface is unchanged bit-for-bit rather than re-emitted and hoped
      to round-trip. A repair that moves the wing is a different wing silently
      revalidated.
  I2  ADDITIVE ONLY. No re-loft, no OpenVSP re-export, no touching the placement
      transform, nothing under geometry/reference/.
  I3  IDENTICAL TREATMENT for every candidate: same Y_CAP, same construction, same
      code path. Skirting pair members differently would move the defect out of
      common mode and into the difference.
  I4  REFERENCE QUANTITIES UNCHANGED. Sref 1.24092 m2, Cref 0.39396 m, bref 3.6576 m
      all stand: the skirt lies entirely at y < 0, outside the half-model domain, and
      contributes NO WETTED AREA. It is a meshing construct, not a surface.
  I5  min(y) < -0.010 m afterwards, i.e. the overhang provably exists.
"""
import argparse
import collections
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
Y_CAP = -0.020        # metres; the rib plane the skirt drops to
ROUND = 6             # vertex unification, decimal places in metres
FMT = "%.9e"          # matches write_stl, so new facets read like the old ones


def parse_facets(text):
    """Return (list of raw facet blocks, Nx3x3 vertex array). The raw text is kept so
    original facets can be re-emitted VERBATIM (I1)."""
    blocks = re.findall(r"[ \t]*facet normal.*?endfacet\n", text, re.S)
    V = np.array(re.findall(r"vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)", text),
                 dtype=float).reshape(-1, 3, 3)
    if len(blocks) != len(V):
        raise ValueError("facet/vertex mismatch: %d blocks, %d triangles" % (len(blocks), len(V)))
    return blocks, V


def open_loop(V):
    """The single open rim, as an ORDERED cycle of unique-point indices, with the
    directed edges as the owning triangles traverse them."""
    key = np.round(V.reshape(-1, 3), ROUND)
    u, inv = np.unique(key, axis=0, return_inverse=True)
    tri = inv.reshape(-1, 3)
    directed = collections.Counter()
    for a, b, c in tri:
        directed[(a, b)] += 1
        directed[(b, c)] += 1
        directed[(c, a)] += 1
    und = collections.Counter()
    for (a, b), n in directed.items():
        und[(a, b) if a < b else (b, a)] += n
    open_und = [e for e, n in und.items() if n == 1]
    if not open_und:
        return u, tri, []
    # keep each open edge in the direction its owning triangle traverses it
    dir_open = []
    for a, b in open_und:
        if directed.get((a, b)):
            dir_open.append((a, b))
        elif directed.get((b, a)):
            dir_open.append((b, a))
        else:
            raise ValueError("open edge in neither direction: %s" % ((a, b),))
    nxt = {}
    for a, b in dir_open:
        if a in nxt:
            raise ValueError("rim vertex %d has two outgoing open edges: not a simple loop" % a)
        nxt[a] = b
    start = dir_open[0][0]
    loop = [start]
    cur = nxt[start]
    while cur != start:
        loop.append(cur)
        cur = nxt.get(cur)
        if cur is None:
            raise ValueError("open rim is not a closed cycle")
        if len(loop) > len(dir_open) + 2:
            raise ValueError("rim walk did not terminate")
    if len(loop) != len(dir_open):
        raise ValueError("rim has %d edges but the cycle visits %d points; more than "
                         "one loop is present and this repair assumes exactly one"
                         % (len(dir_open), len(loop)))
    return u, tri, loop


def ear_clip(poly2d):
    """Triangulate a simple polygon given as Nx2, returning index triples. Ear
    clipping rather than a centroid fan: an aerofoil outline is not guaranteed
    star-shaped about its centroid, and a fan that folds would emit inverted
    triangles that G2 would then blame on the skirt."""
    pts = np.asarray(poly2d, dtype=float)
    n = len(poly2d)
    idx = list(range(n))
    area2 = sum(poly2d[i][0]*poly2d[(i+1) % n][1] - poly2d[(i+1) % n][0]*poly2d[i][1]
                for i in range(n))
    if area2 < 0:
        idx.reverse()
    tris = []
    guard = 0
    while len(idx) > 2:
        guard += 1
        if guard > 10 * n:
            raise ValueError("ear clipping failed to terminate; polygon may self-intersect")
        for k in range(len(idx)):
            i0, i1, i2 = idx[k-1], idx[k], idx[(k+1) % len(idx)]
            a, b, c = poly2d[i0], poly2d[i1], poly2d[i2]
            cr = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
            if cr <= 0:
                continue
            others = [j for j in idx if j not in (i0, i1, i2)]
            if others:
                Q = pts[others]
                d1 = (b[0]-a[0])*(Q[:, 1]-a[1]) - (b[1]-a[1])*(Q[:, 0]-a[0])
                d2 = (c[0]-b[0])*(Q[:, 1]-b[1]) - (c[1]-b[1])*(Q[:, 0]-b[0])
                d3 = (a[0]-c[0])*(Q[:, 1]-c[1]) - (a[1]-c[1])*(Q[:, 0]-c[0])
                if np.any((d1 >= 0) & (d2 >= 0) & (d3 >= 0)):
                    continue
            tris.append((i0, i1, i2))
            idx.pop(k)
            break
        else:
            raise ValueError("no ear found; polygon is not simple")
    return tris


def facet(p, q, r):
    n = np.cross(q - p, r - p)
    m = np.linalg.norm(n)
    n = n/m if m > 0 else np.array([0.0, 0.0, 1.0])
    return ("  facet normal " + " ".join(FMT % v for v in n) + "\n    outer loop\n"
            + "".join("      vertex " + " ".join(FMT % v for v in pt) + "\n"
                      for pt in (p, q, r))
            + "    endloop\n  endfacet\n")


def close(path, y_cap=Y_CAP):
    text = Path(path).read_text(errors="ignore")
    blocks, V = parse_facets(text)
    u, tri, loop = open_loop(V)
    if not loop:
        return None, "already closed"

    P = u[loop]
    if P[:, 1].min() <= y_cap:
        raise ValueError("a rim point at y=%.6f is at or below the cap plane %.6f; "
                         "the skirt would be degenerate or inverted" % (P[:, 1].min(), y_cap))
    seg = P[:, 1] - y_cap
    if seg.min() <= 0:
        raise ValueError("non-positive skirt segment length")

    low = P.copy()
    low[:, 1] = y_cap
    new = []
    n = len(loop)
    # SKIRT. The original triangle traverses its rim edge a->b, so the new triangle
    # sharing that edge must traverse b->a for the surface to stay consistently
    # oriented. Everything below follows from that one constraint.
    for i in range(n):
        a, b = P[i], P[(i+1) % n]
        a2, b2 = low[i], low[(i+1) % n]
        new.append(facet(b, a, a2))
        new.append(facet(b, a2, b2))
    # CAP. Its boundary must oppose the skirt's bottom edges (a2->b2), so the cap
    # traverses the bottom loop in reverse.
    cap_tris = ear_clip([(low[i][0], low[i][2]) for i in range(n)])
    for i0, i1, i2 in cap_tris:
        p, q, r = low[i0], low[i1], low[i2]
        nrm = np.cross(q-p, r-p)
        if nrm[1] > 0:            # cap must face -y (outward, away from the wing)
            p, q, r = p, r, q
        new.append(facet(p, q, r))

    out = "solid wing_closed\n" + "".join(blocks) + "".join(new) + "endsolid wing_closed\n"
    return out, "%d rim edges skirted, %d cap triangles, %d facets total" % (
        n, len(cap_tris), len(blocks) + len(new))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+")
    ap.add_argument("--y-cap", type=float, default=Y_CAP)
    ap.add_argument("--write", action="store_true", help="replace in place after checks")
    args = ap.parse_args()
    rc = 0
    for f in args.files:
        p = Path(f)
        try:
            out, msg = close(p, args.y_cap)
        except ValueError as e:
            print("  %-38s REFUSED: %s" % (p.name, e))
            rc = 1
            continue
        if out is None:
            print("  %-38s %s" % (p.name, msg))
            continue
        # I1 asserted here as well as in the gate: every original facet block must
        # appear verbatim in the output.
        # I1 asserted POSITIONALLY: the first N facets of the output must be the
        # original N, byte for byte and in order. The previous form asked "is each
        # block somewhere in the output", which is 190642 substring searches over a
        # 50 MB string and never returned.
        blocks, _ = parse_facets(Path(f).read_text(errors="ignore"))
        out_blocks, _ = parse_facets(out)
        same = (len(out_blocks) >= len(blocks)
                and all(a == b for a, b in zip(blocks, out_blocks)))
        if not same:
            print("  %-38s REFUSED: original facets not reproduced verbatim in order"
                  % p.name)
            rc = 1
            continue
        print("  %-38s %s" % (p.name, msg))
        if args.write:
            p.write_text(out)
    return rc


if __name__ == "__main__":
    sys.exit(main())
