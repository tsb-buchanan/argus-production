#!/usr/bin/env python3
"""Section-wise comparison of the OpenVSP export against the in-house loft.

WHY NOT A POINT CLOUD METRIC: nearest-neighbour distance between two
independently tessellated surfaces has a floor of roughly half the vertex
spacing. At ~92k vertices over ~2.5 m2 that floor is ~2.6 mm, which is the same
order as the deviation being looked for, so a cloud metric cannot resolve the
question. It reports its own discretisation and looks like a real result.

This script instead slices BOTH surfaces with the same plane, builds the two
section curves, and compares them as curves on a common x/c grid. That has no
sampling floor: an exactly coincident surface returns zero regardless of how
either side was tessellated.

The comparison is done in the LOFT frame: the VSP export is mapped back through
the fitted rigid transform, so sections are at constant y for both and the
placement difference (5 deg dihedral, 1.873 deg incidence) is removed. What is
left is lofting error. Run compare_oml_surfaces.py for the placement term.

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
import struct

import numpy as np


def read_stl_triangles(path):
    with open(path, "rb") as fh:
        if fh.read(5) == b"solid":
            fh.seek(0)
            verts = []
            for line in fh:
                s = line.split()
                if s and s[0] == b"vertex":
                    verts.append((float(s[1]), float(s[2]), float(s[3])))
            a = np.asarray(verts, dtype=float)
            return a.reshape(-1, 3, 3)
        fh.seek(80)
        n = struct.unpack("<I", fh.read(4))[0]
        buf = fh.read(n * 50)
        arr = np.frombuffer(
            np.frombuffer(buf, dtype=np.uint8).reshape(n, 50)[:, 12:48].tobytes(),
            dtype="<f4").reshape(n, 3, 3)
        return arr.astype(float)


def slice_plane_y(tris, y0):
    """Intersect triangles with the plane y = y0, return segment endpoints."""
    v0, v1, v2 = tris[:, 0], tris[:, 1], tris[:, 2]
    d = tris[:, :, 1] - y0
    pos = d > 0
    ncross = pos.sum(axis=1)
    keep = (ncross == 1) | (ncross == 2)
    t = tris[keep]
    dd = d[keep]

    pts = []
    for edge in ((0, 1), (1, 2), (2, 0)):
        a, b = edge
        da, db = dd[:, a], dd[:, b]
        m = (da * db) < 0
        if not m.any():
            continue
        w = (da[m] / (da[m] - db[m]))[:, None]
        p = t[m][:, a, :] + w * (t[m][:, b, :] - t[m][:, a, :])
        pts.append(p)
    if not pts:
        return np.empty((0, 3))
    return np.vstack(pts)


def section_curves(pts, n=400, te_ref="midpoint", band_frac=0.005):
    """Split a section point set into upper and lower surfaces on a common x/c.

    te_ref='midpoint' references the chord to the CENTROID of the aft band
    rather than the single most-aft point. On the in-house blunt TE (0.63%c)
    the two base corners sit at the same x, so a max-x reference picks one
    arbitrarily and rotates the whole normalisation by up to +/- 0.180 deg,
    which contaminates every x/c on the section and not just the TE. That is
    the D031 contaminant; the midpoint reference removes it and degenerates to
    the same point on VSP's sharp closure.
    """
    if len(pts) < 20:
        return None
    x, z = pts[:, 0], pts[:, 2]
    ile = np.argmin(x)
    xle, zle = x[ile], z[ile]
    if te_ref == "midpoint":
        band = pts[x >= x.max() - band_frac * (x.max() - x.min())]
        xte, zte = band[:, 0].mean(), band[:, 2].mean()
    else:
        ite = np.argmax(x)
        xte, zte = x[ite], z[ite]
    chord = np.hypot(xte - xle, zte - zle)
    if chord <= 0:
        return None

    # Rotate into the chord frame so upper/lower split is unambiguous under
    # twist, then normalise by chord.
    ca = (xte - xle) / chord
    sa = (zte - zle) / chord
    xr = ((x - xle) * ca + (z - zle) * sa) / chord
    zr = (-(x - xle) * sa + (z - zle) * ca) / chord

    grid = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))  # cosine, clustered LE/TE
    out = {}
    for name, m in (("upper", zr >= 0), ("lower", zr < 0)):
        xs, zs = xr[m], zr[m]
        if len(xs) < 5:
            return None
        o = np.argsort(xs)
        xs, zs = xs[o], zs[o]
        xs, idx = np.unique(xs, return_index=True)
        zs = zs[idx]
        out[name] = np.interp(grid, xs, zs, left=np.nan, right=np.nan)
    out["grid"] = grid
    out["chord"] = chord
    out["xle"], out["zle"] = xle, zle
    return out


def _probe(grid, arr, xq):
    m = np.isfinite(arr)
    if m.sum() < 5:
        return float("nan")
    return float(100.0 * np.interp(xq, grid[m], arr[m]))


def _slope(grid, arr):
    """Least-squares slope of arr vs x/c through the origin, in %c per x/c."""
    m = np.isfinite(arr)
    if m.sum() < 5:
        return float("nan")
    x, y = grid[m], arr[m]
    return float(100.0 * (x @ y) / (x @ x))


def et_of(rows, key):
    i = int(np.argmax([r[key] for r in rows]))
    return rows[i]["eta"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vsp", required=True)
    ap.add_argument("--loft", required=True)
    ap.add_argument("--transform-json", required=True,
                    help="loft_vs_vsp.json from compare_oml_surfaces.py")
    ap.add_argument("--n-stations", type=int, default=41)
    ap.add_argument("--json-out")
    ap.add_argument("--te-ref", default="midpoint", choices=("midpoint", "maxx"))
    ap.add_argument("--etas", help="explicit comma-separated eta list, overrides --n-stations")
    args = ap.parse_args()

    tf = json.load(open(args.transform_json))["fitted_rigid_transform"]
    # Rebuild R from the reported angles is lossy; recompute from scratch by
    # re-fitting is overkill, so instead map the LOFT forward using the same
    # convention used to derive it: here we only need y-planes to line up, and
    # the fit's translation/rotation is applied to the VSP export inversely.
    # Simplest robust route: fit is unnecessary for section comparison if we
    # slice each surface in ITS OWN frame at the same normalised station.
    vsp = read_stl_triangles(args.vsp)
    loft = read_stl_triangles(args.loft)

    # Keep the y >= 0 half of whichever input is a symmetric full-span export.
    # This MUST be applied to both inputs: applying it to one only maps eta onto
    # opposite halves and the metric silently reports garbage. Caught by the
    # self-comparison null test, which is why that test exists.
    def right_half(t):
        if t[:, :, 1].min() < -1e-6 and t[:, :, 1].max() > 1e-6:
            return t[t[:, :, 1].mean(axis=1) >= 0]
        return t

    vsp = right_half(vsp)
    loft = right_half(loft)

    def yrange(t):
        return t[:, :, 1].min(), t[:, :, 1].max()

    vy0, vy1 = yrange(vsp)
    ly0, ly1 = yrange(loft)

    if args.etas:
        etas = np.array([float(v) for v in args.etas.split(",")])
    else:
        etas = np.linspace(0.03, 0.97, args.n_stations)
    rows, skipped = [], []
    for eta in etas:
        pv = slice_plane_y(vsp, vy0 + eta * (vy1 - vy0))
        pl = slice_plane_y(loft, ly0 + eta * (ly1 - ly0))
        cv = section_curves(pv, te_ref=args.te_ref)
        cl = section_curves(pl, te_ref=args.te_ref)
        if cv is None or cl is None:
            skipped.append((float(eta), "section extraction returned too few points"))
            continue
        du = cv["upper"] - cl["upper"]
        dl = cv["lower"] - cl["lower"]
        both = np.concatenate([du, dl])
        nan_frac = float(np.mean(np.isnan(both)))
        both = both[~np.isnan(both)]
        if not len(both):
            skipped.append((float(eta), "all-NaN after interpolation"))
            continue

        # DECOMPOSITION, the point of this comparison.
        #   camber    = (du + dl)/2  a mean-line difference: REAL loft error.
        #   thickness = (du - dl)    a symmetric difference: this is where the
        #                            blunt-vs-sharp TE convention MUST land,
        #                            because vsp3_to_stl.py's blunt_te term adds
        #                            +b*(x/c)/2 to the upper and -b*(x/c)/2 to
        #                            the lower, which is pure thickness and
        #                            exactly zero camber.
        camber = 0.5 * (du + dl)
        thick = du - dl
        g = cv["grid"]
        iu = int(np.nanargmax(np.abs(du)))
        ic = int(np.nanargmax(np.abs(camber)))
        rows.append({
            "eta": float(eta),
            "chord_vsp_m": cv["chord"],
            "chord_loft_m": cl["chord"],
            "chord_delta_pct": 100.0 * (cl["chord"] - cv["chord"]) / cv["chord"],
            "dz_rms_pct_c": float(100.0 * np.sqrt(np.mean(both ** 2))),
            "dz_max_pct_c": float(100.0 * np.max(np.abs(both))),
            "x_of_max_upper": float(g[iu]),
            "dz_max_upper_pct_c": float(100.0 * du[iu]),
            "camber_rms_pct_c": float(100.0 * np.sqrt(np.nanmean(camber ** 2))),
            "camber_max_pct_c": float(100.0 * np.nanmax(np.abs(camber))),
            "x_of_camber_max": float(g[ic]),
            "thick_rms_pct_c": float(100.0 * np.sqrt(np.nanmean(thick ** 2))),
            # Probe the thickness channel at fixed x/c and fit its slope.
            # vsp3_to_stl.py's blunt_te term is +b*(x/c)/2 on the upper surface
            # and -b*(x/c)/2 on the lower, so (upper-lower) MUST come out as a
            # straight line through the origin with slope -0.63%c. Anything
            # else in this channel is not the TE convention.
            "thick_probe_pct_c": {("%.2f" % xq): _probe(g, thick, xq)
                                  for xq in (0.25, 0.50, 0.75, 0.95)},
            "thick_slope_pct_c_per_xc": _slope(g, thick),
            "nan_fraction": nan_frac,
        })

    rms = np.array([r["dz_rms_pct_c"] for r in rows])
    mx = np.array([r["dz_max_pct_c"] for r in rows])
    et = np.array([r["eta"] for r in rows])

    cam = np.array([r["camber_rms_pct_c"] for r in rows])
    thk = np.array([r["thick_slope_pct_c_per_xc"] for r in rows])
    summary = {
        "te_reference": args.te_ref,
        "n_stations_requested": int(len(etas)),
        "n_stations_compared": len(rows),
        "n_stations_SKIPPED": len(skipped),
        "skipped_detail": skipped,
        "max_nan_fraction_any_station": float(max(r["nan_fraction"] for r in rows)),
        "camber_rms_pct_chord": {"mean": float(cam.mean()), "max": float(cam.max()),
                                 "eta_of_max": float(et_of(rows, "camber_rms_pct_c"))},
        "camber_max_pct_chord": float(max(r["camber_max_pct_c"] for r in rows)),
        "thickness_slope_pct_c_per_xc": {"mean": float(np.nanmean(thk)),
                                         "min": float(np.nanmin(thk)),
                                         "max": float(np.nanmax(thk)),
                                         "EXPECTED_from_blunt_te_convention": -0.63},
        "dz_rms_pct_chord": {"mean": float(rms.mean()), "max": float(rms.max()),
                             "eta_of_max": float(et[np.argmax(rms)])},
        "dz_max_pct_chord": {"max": float(mx.max()), "eta_of_max": float(et[np.argmax(mx)])},
        "chord_delta_pct": {"max_abs": float(max(abs(r["chord_delta_pct"]) for r in rows))},
    }

    print(json.dumps(summary, indent=2))
    print()
    if skipped:
        print("WARNING: %d of %d stations were SKIPPED, listed above. A skipped"
              % (len(skipped), args.n_stations))
        print("station is not a zero deviation; it is an unmeasured one.")
        print()
    print("%6s %10s %10s %10s %10s %10s %10s"
          % ("eta", "chord_d_%", "dz_rms_%c", "camb_rms", "camb_max", "x@cmax", "thk_slope"))
    for r in rows:
        print("%6.3f %10.4f %10.4f %10.4f %10.4f %10.3f %10.4f"
              % (r["eta"], r["chord_delta_pct"], r["dz_rms_pct_c"],
                 r["camber_rms_pct_c"], r["camber_max_pct_c"],
                 r["x_of_camber_max"], r["thick_slope_pct_c_per_xc"]))

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump({"summary": summary, "stations": rows}, fh, indent=2)


if __name__ == "__main__":
    main()
