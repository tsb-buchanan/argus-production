#!/usr/bin/env python3
"""wing_sections.py: spanwise Cp sections from a 3D wing surface sample.

Reads the `wingSurface` sampled-surface VTK that the solve writes, cuts it at constant
eta = y/(b/2), and produces Cp(x/c) at each station plus the spanwise load distribution.

WHY THIS EXISTS SEPARATELY FROM make_report_figures.py: the inputs are ~400 MB ASCII VTK
files that live on the cluster, so this script consumes them once and writes small CSVs.
The figure script then plots the CSVs and stays runnable with no cluster access.

FRAME (D068). Cp = p / (0.5 U^2) with OpenFOAM's KINEMATIC p (incompressible, so p is
already divided by rho). x/c is on the LOCAL chord at that station, taken from the
surface itself, because this wing is tapered and swept -- a global chord would misplace
every station.
"""
import argparse
import csv
import gzip
import re
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]


def read_vtk_cell_centres(path, as_vector=False, with_area=False):
    """Legacy-VTK POLYDATA -> (cell centroid xyz, cell value).

    as_vector=True returns the FULL component array (N,ncomp) instead of the magnitude.
    The default collapses a vector to its magnitude, which is right for Cp-style scalars
    and WRONG for anything whose SIGN carries the answer: reverse flow is wallShearStress
    with a NEGATIVE streamwise component, and a magnitude cannot express it (D060 -- where
    an error is one-directional, test the sign).

    THE DATA IS CELL_DATA, NOT POINT_DATA, and that is not a detail. OpenFOAM's
    sampledSurface writes one value per FACE:

        POINTS      5679704
        POLYGONS    5276990          <- note: FEWER polygons than points
        CELL_DATA   5276990
        FIELD attributes 1
        p 1 5276990 float

    Pairing the value array against the POINT array would silently misalign every value
    by a growing offset and still produce a plausible-looking Cp curve. So the polygon
    connectivity is read and each value is attached to its own face CENTROID.

    Polygons are mixed 3/4/5-sided (snappyHexMesh surfaces are not all triangles), so the
    connectivity is read as a flat stream with a leading count per face rather than
    assumed to be a fixed-width table.
    """
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", errors="ignore") as fh:
        pts = None
        npts = ncell = 0
        ncomp = 1
        mode = None
        buf = []
        conn = []
        vals = []
        for line in fh:
            s = line.strip()
            if not s:
                continue
            head = s.split()[0]
            if head == "POINTS":
                npts = int(s.split()[1]); mode = "pts"; buf = []; continue
            if head == "POLYGONS":
                pts = np.fromiter(buf, dtype=float, count=3 * npts).reshape(npts, 3)
                ncell = int(s.split()[1]); mode = "conn"; buf = []; continue
            if head in ("CELL_DATA", "POINT_DATA"):
                conn = buf; mode = "await"; buf = []; continue
            if mode == "await":
                if head in ("FIELD", "SCALARS", "LOOKUP_TABLE"):
                    if head == "FIELD":
                        mode = "fieldhdr"
                    elif head == "LOOKUP_TABLE":
                        mode = "val"
                    continue
            if mode == "fieldhdr":
                # "<name> <ncomp> <ntuples> <type>". READ ncomp: wallShearStress is a
                # VECTOR (3 components), so assuming 1 would take every third value and
                # silently produce a curve that looks plausible and means nothing.
                f = s.split()
                ncomp = int(f[1]) if len(f) >= 3 and f[1].isdigit() else 1
                mode = "val"; continue
            if mode in ("pts", "conn", "val"):
                buf.extend(s.split())
                continue
        if mode == "val":
            vals = buf

    # walk the flat connectivity: [n, i0..i(n-1)] repeated
    conn = np.asarray(conn, dtype=np.int64)
    V = np.asarray(vals[:ncell * ncomp], dtype=float)
    if ncomp > 1:
        V = V.reshape(-1, ncomp)
        if not as_vector:               # vector field -> magnitude
            V = np.linalg.norm(V, axis=1)
    cx = np.empty((ncell, 3), dtype=float)
    ar = np.empty(ncell, dtype=float) if with_area else None
    k = 0
    for c in range(ncell):
        n = conn[k]; idx = conn[k + 1: k + 1 + n]; k += 1 + n
        q = pts[idx]
        cen = q.mean(axis=0)
        cx[c] = cen
        if with_area:
            # FAN TRIANGULATION ABOUT THE CENTROID, valid for the mixed 3/4/5-sided
            # polygons snappyHexMesh produces and for non-convex ones, where the shoelace
            # formula about an arbitrary vertex is not.
            v = q - cen
            ar[c] = 0.5 * np.linalg.norm(np.cross(v, np.roll(v, -1, axis=0)).sum(axis=0))
    m = min(len(cx), len(V))
    if with_area:
        return cx[:m], V[:m], ar[:m]
    return cx[:m], V[:m]


def section(P, V, eta, semi, band, uinf):
    """Cp(x/c) on a thin spanwise band. Returns (xc, cp, z_over_c) sorted along the chord."""
    y0 = eta * semi
    m = np.abs(P[:, 1] - y0) <= band * semi
    if m.sum() < 50:
        return None
    x, z, p = P[m, 0], P[m, 2], V[m]
    c = x.max() - x.min()
    if c <= 0:
        return None
    xc = (x - x.min()) / c
    cp = p / (0.5 * uinf ** 2)
    return xc, cp, (z - z.mean()) / c


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--base", required=True, help="baseline p_wing.vtk[.gz]")
    ap.add_argument("--morph", help="morphed p_wing.vtk[.gz]")
    ap.add_argument("--eta", type=float, nargs="+", default=[0.30, 0.50, 0.70, 0.80, 0.90])
    ap.add_argument("--band", type=float, default=0.004, help="half-width as a fraction of semispan")
    ap.add_argument("--uinf", type=float, default=40.8)
    ap.add_argument("--out", type=Path, default=REPO / "results/wing_sections")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    for tag, f in (("baseline", a.base), ("morphed", a.morph)):
        if not f:
            continue
        print(f"reading {tag}: {f}")
        P, V = read_vtk_cell_centres(f)
        semi = P[:, 1].max()
        print(f"  {len(P):,} points, semispan {semi:.4f} m, "
              f"p range {V.min():.1f} to {V.max():.1f} (kinematic)")
        for e in a.eta:
            r = section(P, V, e, semi, a.band, a.uinf)
            if r is None:
                print(f"    eta {e:.2f}: too few points in band, skipped")
                continue
            xc, cp, zc = r
            o = np.argsort(xc)
            p = a.out / (f"{tag}_eta{e:.2f}".replace(".", "p") + ".csv")
            with open(p, "w", newline="") as fh:
                w = csv.writer(fh); w.writerow(["x_over_c", "cp", "z_over_c"])
                w.writerows(zip(xc[o], cp[o], zc[o]))
            print(f"    eta {e:.2f}: {len(xc):,} pts, Cp {cp.min():+.3f} to {cp.max():+.3f} -> {p.name}")


if __name__ == "__main__":
    raise SystemExit(main())
