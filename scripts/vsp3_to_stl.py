#!/usr/bin/env python3
"""vsp3_to_stl.py: build a meshable half-wing STL from a delivered .vsp3 OML.

WHY THIS EXISTS, and its one limitation, stated up front: OpenVSP is not
installed and could not be obtained in this environment, so the OML cannot be
exported by VSP's own tessellator. This script lofts the wing from the section
definitions parsed out of the .vsp3 (scripts/vsp3_sections.py). At the 19
DEFINING stations the surface reproduces the delivered section coordinates
exactly, by construction. BETWEEN stations it interpolates LINEARLY in the
chord-normalized ordinates, which is what an OpenVSP WingGeom does for a wing
section, but that equivalence is asserted from the VSP wing-skinning convention
and is NOT verified here against VSP itself. That residual is the reason this
route needs a ruling before it becomes the Phase 1 geometry of record.

Conventions read from the file and applied here:
  sweep at sweep_location 0.0 (leading edge), constant 28.86199 deg
  dihedral 0, twist about twist_location 0.25 (quarter chord)
  native units FEET; --scale 0.3048 converts to metres (project rules, hard fact 1)

Half-span only, root at y = 0, for a symmetry-plane 3D case.

Usage:
  python3 scripts/vsp3_to_stl.py dso_reference/baseline_wing_only_refined.vsp3 \
      --out geometry/derived/wing3d/baseline_oml.stl --nchord 200 --nspan 240
"""

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from vsp3_sections import extract_wing as extract  # noqa: E402


def cosine(n):
    """Chordwise parameter, clustered at LE and TE."""
    return 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, n)))


def resample(pts, s):
    """Airfoil surface (LE->TE, unit chord) resampled onto chordwise fractions s."""
    a = np.asarray(pts, dtype=float)
    x, z = a[:, 0], a[:, 1]
    o = np.argsort(x, kind="stable")
    return np.interp(s, x[o], z[o])


def blend_sections(base, other, f):
    """Return sections at FRACTION f of the morph carried by `other`.

    The morph law is a pure z-displacement applied equally to both surfaces
    (delta_zc = -A S(xi), thickness-preserving), so it is ADDITIVE in the
    ordinates: base + f*(other - base) is exactly f times the amplitude
    schedule. This is not an approximation, it is the definition of the law.

    Everything that is NOT camber must be identical between the two files, so
    it is asserted rather than blended. If a planform quantity ever differs the
    blend is meaningless and this raises instead of silently averaging it.
    """
    if len(base) != len(other):
        raise ValueError("section count differs: %d vs %d" % (len(base), len(other)))
    out = []
    for i, (a, b) in enumerate(zip(base, other)):
        for k in ("eta", "y", "chord", "twist", "sweep", "span"):
            va, vb = a.get(k), b.get(k)
            if va is None or vb is None:
                continue
            if abs(float(va) - float(vb)) > 1e-9 * max(1.0, abs(float(va))):
                raise ValueError("section %d: %s differs between files (%r vs %r); "
                                 "the two geometries are not a pure camber pair"
                                 % (i, k, va, vb))
        c = dict(a)
        for surf in ("upper", "lower"):
            if a.get(surf) is None:
                continue
            pa, pb = np.asarray(a[surf], float), np.asarray(b[surf], float)
            if pa.shape != pb.shape:
                raise ValueError("section %d: %s point counts differ (%s vs %s)"
                                 % (i, surf, pa.shape, pb.shape))
            dx = np.abs(pa[:, 0] - pb[:, 0]).max()
            if dx > 1e-9:
                raise ValueError("section %d: %s abscissae differ by %.3e; cannot blend "
                                 "ordinates on mismatched x" % (i, surf, dx))
            z = pa[:, 1] + f * (pb[:, 1] - pa[:, 1])
            c[surf] = np.column_stack([pa[:, 0], z])
        out.append(c)
    return out


def build(sections, nchord, blunt_te=0.0):
    """Per defining station: chord-normalized upper/lower on a common distribution,
    plus the station's chord, twist, and leading-edge position.

    blunt_te INVERTS the full-chord linear closure that OpenVSP applied to make the
    surface close (D031). VSP removed blunt_te*(x/c) of thickness; this adds back
    blunt_te*(x/c)/2 on each surface, normal to the chord, in unit-chord coordinates
    so it is automatically 0.63% of the LOCAL chord at every station. It recovers the
    geometry D004 validated rather than inventing a second blunting rule: a fresh rule
    would give a different blunt geometry and break the link to the 2D programme.
    LE closure is untouched (the offset vanishes at x/c = 0)."""
    s = cosine(nchord)
    out = []
    x_le = 0.0
    for i, sec in enumerate(sections):
        if i > 0:
            sw = sec["sweep"] if sec["sweep"] is not None else 0.0
            x_le += sec["span"] * math.tan(math.radians(sw))
        up, lo = resample(sec["upper"], s), resample(sec["lower"], s)
        if blunt_te:
            up = up + blunt_te * s / 2.0
            lo = lo - blunt_te * s / 2.0
        out.append({"eta": sec["eta"], "y": sec["y"], "chord": sec["chord"],
                    "twist": sec["twist"] or 0.0, "x_le": x_le, "up": up, "lo": lo})
    return s, out


def station_points(st, s, nchord):
    """3D points of one station: upper LE->TE then lower TE->LE, twisted and placed."""
    c, t = st["chord"], math.radians(st["twist"])
    xs = np.concatenate([s, s[::-1][1:-1]])
    zs = np.concatenate([st["up"], st["lo"][::-1][1:-1]])
    # twist about the quarter-chord of the unit airfoil, then scale and place
    xr = (xs - 0.25) * math.cos(t) + zs * math.sin(t) + 0.25
    zr = -(xs - 0.25) * math.sin(t) + zs * math.cos(t)
    return np.column_stack([st["x_le"] + c * xr,
                            np.full(xs.size, st["y"]),
                            c * zr])


def interp_station(a, b, f, s):
    return {"eta": a["eta"] + f * (b["eta"] - a["eta"]),
            "y": a["y"] + f * (b["y"] - a["y"]),
            "chord": a["chord"] + f * (b["chord"] - a["chord"]),
            "twist": a["twist"] + f * (b["twist"] - a["twist"]),
            "x_le": a["x_le"] + f * (b["x_le"] - a["x_le"]),
            "up": a["up"] + f * (b["up"] - a["up"]),
            "lo": a["lo"] + f * (b["lo"] - a["lo"])}


def spanwise(stations, nspan, s):
    """Defining stations ALWAYS retained; extra stations inserted per panel in
    proportion to panel span so tessellation is even, not clustered."""
    spans = [stations[i + 1]["y"] - stations[i]["y"] for i in range(len(stations) - 1)]
    total = sum(spans)
    extra = [max(1, int(round((nspan - len(stations)) * sp / total))) for sp in spans]
    out = [stations[0]]
    for i, n in enumerate(extra):
        for k in range(1, n + 1):
            out.append(interp_station(stations[i], stations[i + 1], k / (n + 1), s))
        out.append(stations[i + 1])
    return out


def cap_triangles(ring, prev_ring, target_edge):
    """A WELL-SHAPED tip cap, replacing the centroid fan.

    WHAT THE FAN DID, AND WHY IT HAD TO GO. The original cap was
        ctr = tip.mean(axis=0)
        for j in range(npts): tris.append((tip[j], tip[j+1], ctr))
    Every triangle therefore had two spokes running from the section centroid out to
    the outline (up to ~100 mm on a 199 mm tip chord) and a base equal to the local
    outline spacing (~1 mm, and far tighter where cosine clustering packs points at the
    LE and TE). MEASURED CONSEQUENCE on the delivered STL: one vertex of valence 398,
    tip-band aspect ratios to 35,212, 4.9% of tip triangles above AR 100, and 219
    triangles carrying edges over 50 mm. snappyHexMesh reproduced that surface
    faithfully, which is why checkMesh never complained and why the terracing survived
    every mesh setting we tried (GEO-102: a gate on the output cannot see a wrong input).

    WHY NOT DELAUNAY, MEASURED RATHER THAN ARGUED. The first attempt seeded interior
    points, ran an unconstrained Delaunay over boundary+seeds and culled triangles whose
    centroid fell outside the outline. Cap quality was excellent -- median aspect ratio
    2.94 against the fan's 315.95, longest edge 7.8 mm against 76.2 mm -- and
    surfaceCheck then reported "Surface is not closed", 4 unconnected parts, 48 zones.
    UNCONSTRAINED DELAUNAY DOES NOT GUARANTEE THE POLYGON BOUNDARY APPEARS IN THE
    TRIANGULATION, so culling by centroid leaves gaps along the rim. Closure is not
    negotiable here: it is the defect that suspended every 3D result in this project.

    SO THE CAP IS BUILT STRUCTURALLY. The tip ring is an aerofoil outline: upper LE->TE
    followed by lower TE->LE. Pairing each upper point with the lower point at the same
    chordwise station gives a ladder of ribs, and stitching adjacent ribs consumes every
    boundary edge EXACTLY ONCE by construction. Closure cannot break, because the rim is
    an output of the algorithm rather than something it has to rediscover.
    """
    P = np.asarray(ring, float)
    n = len(P)
    nc = (n + 2) // 2          # points per surface; ring = upper(nc) + lower(nc-2)

    def lower_of(i):
        """Ring index of the lower-surface point at the same chordwise station as
        upper point i. LE (i=0) and TE (i=nc-1) are SHARED single points."""
        if i == 0 or i == nc - 1:
            return i
        return nc + (nc - 2 - i)

    def rib(i, m):
        """m+1 points from the upper surface down to the lower at station i."""
        a, b = P[i], P[lower_of(i)]
        return [a + (b - a) * (k / m) for k in range(m + 1)]

    # Rib subdivision follows LOCAL thickness, so a 24 mm rib at mid-chord is cut into
    # target-sized pieces while a vanishing rib at the LE stays a single point.
    ms = []
    for i in range(nc):
        h = float(np.linalg.norm(P[i] - P[lower_of(i)]))
        ms.append(max(1, int(round(h / target_edge))))

    tris = []
    for i in range(nc - 1):
        A, B = rib(i, ms[i]), rib(i + 1, ms[i + 1])
        # Two-pointer stitch between polylines of unequal length: advance whichever
        # side is behind in normalised arc position. Every segment of both ribs is
        # consumed exactly once, so no edge is left unpaired.
        ia = ib = 0
        while ia < ms[i] or ib < ms[i + 1]:
            fa, fb = ia / ms[i], ib / ms[i + 1]
            if ib >= ms[i + 1] or (ia < ms[i] and fa <= fb):
                tris.append((A[ia], A[ia + 1], B[ib])); ia += 1
            else:
                tris.append((A[ia], B[ib + 1], B[ib])); ib += 1

    # OUTWARD is away from the wing: along (tip centroid - previous station centroid).
    out = P.mean(axis=0) - np.asarray(prev_ring, float).mean(axis=0)
    fixed = []
    for p0, p1, p2 in tris:
        nrm = np.cross(p1 - p0, p2 - p0)
        if np.dot(nrm, out) < 0:
            p1, p2 = p2, p1
        if np.linalg.norm(np.cross(p1 - p0, p2 - p0)) > 0.0:
            fixed.append((p0, p1, p2))
    return fixed


def write_stl(path, rings, scale, name="wing", bands=None, etas=None,
              cap="fan", cap_edge=None):
    """Triangulated skin plus a flat tip cap. Root is left OPEN: it is the
    symmetry plane and snappy takes it from the background mesh.

    bands=None writes ONE solid (the registered production geometry, byte-for-byte
    unchanged by this option existing). bands=(n, lo, hi) instead writes n named
    spanwise solids between eta lo and hi plus root/tip guards, so snappy sees them
    as separate patches and each can carry its own layer specification. That is what
    lets a single mesh run sweep the requested layer thickness."""
    groups = {}

    def quad(g, p0, p1, p2, p3):
        groups.setdefault(g, []).append((p0, p1, p2))
        groups[g].append((p0, p2, p3))

    R = [r * scale for r in rings]
    npts = R[0].shape[0]

    def band_of(e):
        if bands is None:
            return name
        n, lo, hi = bands
        if e < lo:
            return "root_guard"
        if e >= hi:
            return "tip_guard"
        return f"band{int((e - lo) / (hi - lo) * n) + 1:02d}"

    for i in range(len(R) - 1):
        A, B = R[i], R[i + 1]
        e = 0.5 * (etas[i] + etas[i + 1]) if etas is not None else 0.0
        g = band_of(e)
        for j in range(npts):
            k = (j + 1) % npts
            quad(g, A[j], A[k], B[k], B[j])
    tip = R[-1]
    gtip = band_of(1.0)
    if cap == "fan":
        # THE ORIGINAL, KEPT ON PURPOSE AND NOT AS A COURTESY. registry/candidates.yaml
        # records a `regenerate_with` command for every committed geometry, and those
        # commands must keep reproducing their recorded sha256. Changing the DEFAULT
        # would silently invalidate the provenance of all eleven at once. New geometry
        # is generated with --cap delaunay; the old files stay regenerable.
        ctr = tip.mean(axis=0)
        for j in range(npts):
            k = (j + 1) % npts
            groups.setdefault(gtip, []).append((tip[j], tip[k], ctr))
    elif cap == "delaunay":
        # Default target edge is the median outline spacing, so the cap is as fine as
        # the rim it stitches to rather than an arbitrary length.
        te = cap_edge
        if te is None:
            seg = np.linalg.norm(tip - np.roll(tip, -1, axis=0), axis=1)
            te = float(np.median(seg))
        groups.setdefault(gtip, []).extend(cap_triangles(tip, R[-2], te))
    else:
        raise ValueError(f"unknown cap mode {cap!r}")

    with open(path, "w") as f:
        for g, tris in groups.items():
            f.write(f"solid {g}\n")
            _facets(f, tris)
            f.write(f"endsolid {g}\n")
    return sum(len(v) for v in groups.values())


def _facets(f, tris):
        for a, b, c in tris:
            n = np.cross(b - a, c - a)
            m = np.linalg.norm(n)
            n = n / m if m > 0 else np.array([0.0, 0.0, 1.0])
            f.write(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n    outer loop\n")
            for p in (a, b, c):
                f.write(f"      vertex {p[0]:.9e} {p[1]:.9e} {p[2]:.9e}\n")
            f.write("    endloop\n  endfacet\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("vsp3")
    ap.add_argument("--geom", default=None)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nchord", type=int, default=200, help="points per surface")
    ap.add_argument("--nspan", type=int, default=240, help="total spanwise stations")
    ap.add_argument("--scale", type=float, default=0.3048, help="feet -> metres")
    ap.add_argument("--bands", type=int, default=None,
                    help="split the span into N named solids (layer-floor trial only)")
    ap.add_argument("--band-eta", default="0.12,0.88", help="lo,hi eta for the banded region")
    ap.add_argument("--blend-with", default=None,
                    help="second .vsp3; loft at --blend-fraction of the morph it carries")
    ap.add_argument("--blend-fraction", type=float, default=1.0,
                    help="fraction of the --blend-with amplitude schedule (0 = this file)")
    # DEFAULT IS "fan" SO EVERY registry regenerate_with COMMAND STILL REPRODUCES ITS
    # RECORDED sha256. The fan is a measured defect (valence 398, tip AR to 35,212);
    # it stays the default only until the registry is migrated to delaunay entries.
    ap.add_argument("--cap", choices=["fan", "delaunay"], default="fan",
                    help="tip cap: 'fan' reproduces the committed geometry byte-for-byte; "
                         "'delaunay' is the repaired cap (max AR ~4)")
    ap.add_argument("--cap-edge", type=float, default=None,
                    help="target cap edge length in METRES (rings are scaled before capping). "
                         "Default: median outline spacing of the tip ring.")
    ap.add_argument("--blunt-te", type=float, default=0.0,
                    help="invert VSP's full-chord linear closure; 0.0063 restores the "
                         "measured 0.63%%c blunt TE (D031)")
    args = ap.parse_args()

    data = extract(args.vsp3, args.geom)
    secs = data["sections"]
    if args.blend_with:
        secs = blend_sections(secs, extract(args.blend_with, args.geom)["sections"],
                              args.blend_fraction)
    s, defining = build(secs, args.nchord, args.blunt_te)
    sts = spanwise(defining, args.nspan, s)
    rings = [station_points(st, s, args.nchord) for st in sts]
    etas = [st["eta"] for st in sts]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    bands = None
    if args.bands:
        lo, hi = (float(v) for v in args.band_eta.split(","))
        bands = (args.bands, lo, hi)
    ntri = write_stl(args.out, rings, args.scale, bands=bands, etas=etas,
                     cap=args.cap, cap_edge=args.cap_edge)
    sha = hashlib.sha256(args.out.read_bytes()).hexdigest()

    # geometry of record: planform in METRES, recomputed from the lofted stations.
    # rings are still in NATIVE units here (write_stl scales a copy), so scale both
    # coordinates before integrating: mixing them silently costs one factor of 0.3048.
    ys = np.array([r[0, 1] for r in rings]) * args.scale
    chords = np.array([st["chord"] for st in spanwise(defining, args.nspan, s)]) * args.scale
    half_area = float(np.trapz(chords, ys))
    meta = {
        "source_vsp3": args.vsp3, "geom": data["geom"],
        "tessellation": {"n_chordwise_points_per_surface": args.nchord,
                         "n_spanwise_stations": len(rings),
                         "chordwise_distribution": "cosine (clustered LE and TE)",
                         "spanwise_rule": "all 19 defining stations retained; extras "
                                          "inserted per panel in proportion to panel span",
                         "triangles": ntri},
        "units": {"native": "feet", "scale_applied": args.scale, "output": "metres"},
        "half_span_m": float(ys.max()), "half_area_m2": half_area,
        "full_span_m": float(2 * ys.max()), "full_area_m2": 2 * half_area,
        "root_chord_m": float(chords[0]), "tip_chord_m": float(chords[-1]),
        "stl_sha256": sha,
        "blunt_te_c": args.blunt_te,
        "blend": ({"with": args.blend_with, "fraction": args.blend_fraction,
                   "law": "base + f*(other - base) on the ordinates; exact because the "
                          "morph is a pure additive z-displacement"}
                  if args.blend_with else None),
        "te_convention": ("blunt 0.63%c, VSP's full-chord linear closure INVERTED (D031): "
                          "matches the 2D programme and the measured hardware"
                          if args.blunt_te else "sharp-closed as delivered by VSP"),
        "limitation": "lofted in-house; exact at the 19 defining stations, LINEAR between "
                      "them. Not verified against OpenVSP's own tessellator (not installed).",
    }
    args.out.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
