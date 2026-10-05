#!/usr/bin/env python3
"""Apply the parent-Blank placement transform to an in-house lofted wing STL.

The defect in geometry/derived/wing3d/baseline_oml.stl is a MISSING RIGID
TRANSFORM, not a wrong surface: the 5 deg dihedral and 1.873 deg incidence are
carried by the parent Blank geom 'wing_comps', which the loft route never walks
(D048 item 6). This applies that transform and nothing else, so the blunt TE,
the D004 fairing and the planform validated to 5e-7 all survive untouched.

ROTATION ORDER IS DETERMINED EMPIRICALLY, not assumed. VSP's composition order
is a fact about VSP, so both orderings are built and scored against the OpenVSP
export; --order auto picks the winner and reports the margin.

Rotation is about the wing root leading edge, which in the loft frame is the
origin: the loft places the root LE at (0, 0, 0) and the parent Blank carries
the translation to model coordinates.
"""
import argparse
import hashlib
import json
import math
import struct

import numpy as np

FT_TO_M = 0.3048
# Parent Blank 'wing_comps', read from the .vsp3 by scripts/vsp_probe.py.
DIHEDRAL_DEG = 5.000
INCIDENCE_DEG = 1.873
TRANSLATE_FT = (3.860, 0.0, -0.400)


def read_stl_triangles(path):
    with open(path, "rb") as fh:
        if fh.read(5) == b"solid":
            fh.seek(0)
            verts = []
            for line in fh:
                s = line.split()
                if s and s[0] == b"vertex":
                    verts.append((float(s[1]), float(s[2]), float(s[3])))
            return np.asarray(verts, dtype=float).reshape(-1, 3, 3)
        fh.seek(80)
        n = struct.unpack("<I", fh.read(4))[0]
        buf = fh.read(n * 50)
        return np.frombuffer(
            np.frombuffer(buf, dtype=np.uint8).reshape(n, 50)[:, 12:48].tobytes(),
            dtype="<f4").reshape(n, 3, 3).astype(float)


def write_stl_ascii(path, tris, name="wing"):
    with open(path, "w") as fh:
        fh.write("solid %s\n" % name)
        for t in tris:
            fh.write("  facet normal 0 0 0\n    outer loop\n")
            for v in t:
                fh.write("      vertex %.9e %.9e %.9e\n" % tuple(v))
            fh.write("    endloop\n  endfacet\n")
        fh.write("endsolid %s\n" % name)


def rot_x(deg):
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)


def rot_y(deg):
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=float)


def build_R(order, dihedral, incidence):
    """order 'xy' means Rx applied LAST, i.e. R = Rx @ Ry."""
    if order == "xy":
        return rot_x(dihedral) @ rot_y(incidence)
    if order == "yx":
        return rot_y(incidence) @ rot_x(dihedral)
    raise ValueError("order must be xy or yx")


def landmarks(tris):
    """Root and tip LE/TE, used to score a candidate transform."""
    v = tris.reshape(-1, 3)
    v = v[v[:, 1] >= -1e-9]
    y = v[:, 1]
    out = {}
    for tag, sel in (("root", y < y.min() + 0.006), ("tip", y > y.max() - 0.006)):
        sec = v[sel]
        out[tag + "_le"] = sec[np.argmin(sec[:, 0])]
        out[tag + "_te"] = sec[np.argmax(sec[:, 0])]
    return out


def score(tris, ref):
    """RMS landmark distance to the reference surface's landmarks."""
    a, b = landmarks(tris), landmarks(ref)
    d = [np.linalg.norm(a[k] - b[k]) for k in a]
    return float(np.sqrt(np.mean(np.square(d)))), {k: float(np.linalg.norm(a[k] - b[k])) for k in a}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loft", required=True, help="in-house loft STL, metres")
    ap.add_argument("--out", required=True)
    ap.add_argument("--reference", help="OpenVSP export STL (metres) used to pick the order")
    ap.add_argument("--order", default="auto", choices=("auto", "xy", "yx"))
    ap.add_argument("--dihedral-deg", type=float, default=DIHEDRAL_DEG)
    ap.add_argument("--incidence-deg", type=float, default=INCIDENCE_DEG)
    ap.add_argument("--json-out")
    args = ap.parse_args()

    tris = read_stl_triangles(args.loft)
    t = np.array(TRANSLATE_FT, dtype=float) * FT_TO_M

    # The loft must place the root LE at the origin for a rotation about the
    # origin to be a rotation about the root LE. Verify rather than assume.
    lm = landmarks(tris)
    root_le = lm["root_le"]
    if np.linalg.norm(root_le) > 2e-3:
        raise RuntimeError(
            "loft root LE is at %s, not the origin; the rotation centre "
            "assumption does not hold for this file" % root_le.tolist())

    results = {}
    if args.order == "auto":
        if not args.reference:
            raise SystemExit("--order auto requires --reference")
        ref = read_stl_triangles(args.reference)
        for o in ("xy", "yx"):
            R = build_R(o, args.dihedral_deg, args.incidence_deg)
            cand = (tris @ R.T) + t
            rms, per = score(cand, ref)
            results[o] = {"landmark_rms_m": rms, "per_landmark_m": per}
        order = min(results, key=lambda k: results[k]["landmark_rms_m"])
        margin = abs(results["xy"]["landmark_rms_m"] - results["yx"]["landmark_rms_m"])
    else:
        order = args.order
        margin = None

    R = build_R(order, args.dihedral_deg, args.incidence_deg)
    out_tris = (tris @ R.T) + t
    write_stl_ascii(args.out, out_tris, name="wing_transformed")

    meta = {
        "source_loft": args.loft,
        "source_loft_sha256": sha256(args.loft),
        "transform": {
            "dihedral_deg_about_x": args.dihedral_deg,
            "incidence_deg_about_y": args.incidence_deg,
            "rotation_centre": "wing root LE = loft origin (verified, |root_LE| = %.2e m)"
                               % float(np.linalg.norm(root_le)),
            "order_applied": order,
            "order_meaning": "xy => R = Rx @ Ry (Rx applied last)",
            "order_selection": "empirical against the OpenVSP export"
                               if args.order == "auto" else "forced by --order",
            "translation_ft": list(TRANSLATE_FT),
            "translation_m": t.tolist(),
            "rotation_matrix": R.tolist(),
        },
        "order_scores_landmark_rms_m": {k: v["landmark_rms_m"] for k, v in results.items()},
        "order_margin_m": margin,
        "order_detail": results,
        "triangles": int(len(out_tris)),
        "out": args.out,
        "out_sha256": sha256(args.out),
    }
    print(json.dumps(meta, indent=2))
    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(meta, fh, indent=2)


if __name__ == "__main__":
    main()
