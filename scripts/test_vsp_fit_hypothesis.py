#!/usr/bin/env python3
"""D049 amendment 1 item 3: does OpenVSP's surface pass through the STORED FILE
AIRFOIL POINTS?

Pre-committed branches are in docs/decisions.md D049 amendment 1 item 3 and are
deliberately not restated here.

Method: map the OpenVSP export BACK into the wing frame by inverting the parent
placement transform, slice it at a DEFINING station's exact y, normalise the
section to unit chord, and compare against the .vsp3's stored FILE airfoil
ordinates for that station.

The comparison is clean of the TE-convention problem that dogged the loft
comparison: the stored file points are ALREADY CLOSED (vsp3_to_stl.py's
blunt_te term exists precisely to re-open them), and the VSP export is closed
too, so both sides are sharp and no blunt/sharp reconciliation is involved.

!!! WARNING, D056 2026-07-28: THIS SCRIPT'S UPPER/LOWER SPLIT IS DEFECTIVE !!!
It splits surfaces by the SIGN OF z in the chord frame. That is wrong for the
EET aft-cambered supercritical sections, whose LOWER surface rises ABOVE the
chord line from about x/c 0.78 to the trailing edge (32 of 87 stored points at
station 4). Those points are labelled UPPER and corrupt the comparison over
exactly x/c 0.78-1.00. Every number this script produced in the x/c 0.7-1.0
band is SUPERSEDED; see docs/decisions.md D056 and use scripts/close_d031.py,
whose split_branches() is validated against stored branch membership at all 19
defining stations with zero misassignments.
"""
import argparse
import json
import math
import struct
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from vsp3_sections import extract_wing  # noqa: E402

FT = 0.3048
DIHEDRAL_DEG, INCIDENCE_DEG = 5.000, 1.873
TRANSLATE_FT = np.array([3.860, 0.0, -0.400])


def read_stl_triangles(path):
    with open(path, "rb") as fh:
        if fh.read(5) == b"solid":
            fh.seek(0)
            v = []
            for line in fh:
                s = line.split()
                if s and s[0] == b"vertex":
                    v.append((float(s[1]), float(s[2]), float(s[3])))
            return np.asarray(v, float).reshape(-1, 3, 3)
        fh.seek(80)
        n = struct.unpack("<I", fh.read(4))[0]
        buf = fh.read(n * 50)
        return np.frombuffer(
            np.frombuffer(buf, np.uint8).reshape(n, 50)[:, 12:48].tobytes(),
            "<f4").reshape(n, 3, 3).astype(float)


def rot_x(d):
    c, s = math.cos(math.radians(d)), math.sin(math.radians(d))
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def rot_y(d):
    c, s = math.cos(math.radians(d)), math.sin(math.radians(d))
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def slice_y(tris, y0):
    d = tris[:, :, 1] - y0
    keep = ~((d > 0).all(axis=1) | (d < 0).all(axis=1))
    t, dd = tris[keep], d[keep]
    pts = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        m = (dd[:, a] * dd[:, b]) < 0
        if m.any():
            w = (dd[m][:, a] / (dd[m][:, a] - dd[m][:, b]))[:, None]
            pts.append(t[m][:, a, :] + w * (t[m][:, b, :] - t[m][:, a, :]))
    return np.vstack(pts) if pts else np.empty((0, 3))


def to_unit_chord(pts, band_frac=0.005):
    x, z = pts[:, 0], pts[:, 2]
    le = pts[np.argmin(x)]
    band = pts[x >= x.max() - band_frac * (x.max() - x.min())]
    te = band.mean(axis=0)
    ch = math.hypot(te[0] - le[0], te[2] - le[2])
    ca, sa = (te[0] - le[0]) / ch, (te[2] - le[2]) / ch
    xr = ((x - le[0]) * ca + (z - le[2]) * sa) / ch
    zr = (-(x - le[0]) * sa + (z - le[2]) * ca) / ch
    return xr, zr, ch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vsp-stl", required=True)
    ap.add_argument("--vsp3", default="dso_reference/baseline_wing_only_refined.vsp3")
    ap.add_argument("--stations", default="4,9,14")
    ap.add_argument("--json-out")
    args = ap.parse_args()

    tris = read_stl_triangles(args.vsp_stl)
    tris = tris[tris[:, :, 1].mean(axis=1) >= 0]          # right half

    # Invert the parent placement transform: p_wing = R^T (p_export - t)
    R = rot_x(DIHEDRAL_DEG) @ rot_y(INCIDENCE_DEG)
    t = TRANSLATE_FT * FT
    wing = (tris - t) @ R

    data = extract_wing(args.vsp3)
    secs = data["sections"] if isinstance(data, dict) else data

    out = []
    for i in [int(v) for v in args.stations.split(",")]:
        sec = secs[i]
        if sec.get("upper") is None:
            continue
        y0 = float(sec["y"]) * FT
        pts = slice_y(wing, y0)
        if len(pts) < 50:
            print("station %d: slice returned %d points, SKIPPED" % (i, len(pts)))
            continue
        xr, zr, ch = to_unit_chord(pts)

        res = {"station": i, "eta": sec.get("eta"), "n_slice_pts": int(len(pts))}
        for surf in ("upper", "lower"):
            fp = np.asarray(sec[surf], float)
            m = zr >= 0 if surf == "upper" else zr < 0
            xs, zs = xr[m], zr[m]
            o = np.argsort(xs)
            xs, zs = xs[o], zs[o]
            xu, idx = np.unique(xs, return_index=True)
            zu = zs[idx]
            # Evaluate the VSP surface AT the stored file abscissae.
            inb = (fp[:, 0] >= xu.min()) & (fp[:, 0] <= xu.max())
            dz = np.interp(fp[inb, 0], xu, zu) - fp[inb, 1]
            res[surf] = {
                "n_file_points_compared": int(inb.sum()),
                "rms_pct_c": float(100 * np.sqrt(np.mean(dz ** 2))),
                "max_pct_c": float(100 * np.abs(dz).max()),
                "x_of_max": float(fp[inb][np.argmax(np.abs(dz)), 0]),
            }
        out.append(res)

    print("%5s %7s %10s %10s %9s %10s %10s %9s"
          % ("stn", "eta", "up_rms%c", "up_max%c", "x@up_max", "lo_rms%c", "lo_max%c", "x@lo_max"))
    for r in out:
        print("%5d %7.3f %10.4f %10.4f %9.3f %10.4f %10.4f %9.3f"
              % (r["station"], r["eta"], r["upper"]["rms_pct_c"], r["upper"]["max_pct_c"],
                 r["upper"]["x_of_max"], r["lower"]["rms_pct_c"], r["lower"]["max_pct_c"],
                 r["lower"]["x_of_max"]))
    print()
    allmax = max(max(r["upper"]["max_pct_c"], r["lower"]["max_pct_c"]) for r in out)
    allrms = max(max(r["upper"]["rms_pct_c"], r["lower"]["rms_pct_c"]) for r in out)
    print("worst RMS %.4f %%c, worst max %.4f %%c, against the 0.277 %%c camber residual"
          % (allrms, allmax))

    if args.json_out:
        json.dump(out, open(args.json_out, "w"), indent=2)


if __name__ == "__main__":
    main()
