#!/usr/bin/env python3
"""vsp3_sections.py: extract wing cross-sections from an OpenVSP .vsp3 file.

Reads the .vsp3 XML directly (no OpenVSP install needed) and returns, per
wing XSec: station eta, y, chord, twist, sweep, dihedral, and the FILE-type
airfoil point sets (unit-chord upper/lower surfaces). Also computes planform
totals (span, area, MAC, aspect ratio) from the section parameters.

Units: values are returned in the NATIVE model unit of the file. OpenVSP
exports for ARGUS are in FEET (project rules, hard fact 1); this script does NOT
scale. Airfoil point sets are unit-chord and therefore unit-free.

Conventions verified against dso_reference/baseline_wing_only_refined.vsp3
(2026-07-24):
1. WingGeom/XSecSurf holds N XSec nodes. XSec[0] is the root station; its
   chord is its Tip_Chord parameter. XSec[i>=1] describes panel i, whose
   outboard (tip) station carries chord = Tip_Chord, twist = Twist, and the
   panel Span / Sweep / Dihedral parameters.
2. eta_i = cumsum(Span[1:]) / halfspan; y_i = cumsum(Span[1:]).
3. FILE airfoils store UpperPnts / LowerPnts as comma-separated x,z,0
   triples, LE (0,0) to TE (~1,~0), unit chord.
4. Airfoil names are NOT reliable station identifiers (delivered refined
   files repeat names across inserted stations); match stations by eta.

CLI:
    python3 scripts/vsp3_sections.py FILE.vsp3 [--geom NAME] [--list]
        [--csv OUT.csv] [--json OUT.json] [--dump-dat DIR]

API:
    list_wing_geoms(path) -> [names]
    extract_wing(path, geom=None) -> {"file", "geom", "planform", "sections"}
    section_properties(sec, n=400) -> t/c max, x_tmax, camber max, x_cmax, ...
"""

import argparse
import csv
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


def _parm(node, name):
    p = node.find(name)
    return float(p.attrib["Value"]) if p is not None else None


def _points(text):
    vals = [float(v) for v in text.strip().rstrip(",").split(",")]
    pts = np.asarray(vals).reshape(-1, 3)
    return pts[:, :2].copy()  # (x, z), third column is 0


def _wing_geoms(root):
    out = []
    for g in root.iter("Geom"):
        pc, wg = g.find("ParmContainer"), g.find("WingGeom")
        if pc is None or wg is None:
            continue
        nm = pc.find("Name")
        out.append((nm.text if nm is not None else "?", g, wg))
    return out


def list_wing_geoms(path):
    root = ET.parse(str(path)).getroot()
    return [nm for nm, _, _ in _wing_geoms(root)]


def extract_wing(path, geom=None):
    """Extract sections of one wing geom. geom=None: sole wing geom or error."""
    root = ET.parse(str(path)).getroot()
    wings = _wing_geoms(root)
    if geom is not None:
        wings = [w for w in wings if w[0] == geom]
    if len(wings) != 1:
        names = [nm for nm, _, _ in _wing_geoms(root)]
        raise ValueError(f"{path}: need exactly one wing geom, found "
                         f"{[w[0] for w in wings]} (available: {names})")
    name, _, wg = wings[0]

    xsecs = wg.find("XSecSurf").findall("XSec")
    spans = []
    secs = []
    for i, xs in enumerate(xsecs):
        wp = xs.find("ParmContainer/XSec")
        fa = xs.find("XSec/XSecCurve/FileAirfoil")
        sec = {
            "index": i,
            "airfoil_name": (fa.find("AirfoilName").text
                             if fa is not None else None),
            "chord": _parm(wp, "Tip_Chord"),
            "twist": _parm(wp, "Twist") if i else (_parm(wp, "Twist") or 0.0),
            "twist_location": _parm(wp, "Twist_Location"),
            "span": _parm(wp, "Span") if i else 0.0,
            "sweep": _parm(wp, "Sweep") if i else None,
            "sweep_location": _parm(wp, "Sweep_Location") if i else None,
            "dihedral": _parm(wp, "Dihedral") if i else None,
            "panel_area": _parm(wp, "Area") if i else 0.0,
            "upper": _points(fa.find("UpperPnts").text) if fa is not None else None,
            "lower": _points(fa.find("LowerPnts").text) if fa is not None else None,
        }
        spans.append(sec["span"])
        secs.append(sec)

    y = np.cumsum(spans)
    half = y[-1]
    for sec, yi in zip(secs, y):
        sec["y"] = float(yi)
        sec["eta"] = float(yi / half)

    area_half = sum(s["panel_area"] for s in secs[1:])
    span_full = 2.0 * half
    area_full = 2.0 * area_half
    planform = {
        "halfspan": float(half),
        "span": float(span_full),
        "area": float(area_full),
        "aspect_ratio": float(span_full**2 / area_full),
        "root_chord": secs[0]["chord"],
        "tip_chord": secs[-1]["chord"],
        "n_sections": len(secs),
    }
    return {"file": str(path), "geom": name,
            "planform": planform, "sections": secs}


def resample_surfaces(sec, x=None, n=400):
    """Cubic-spline upper/lower z(x) on a common unit-chord grid (LE-clustered)."""
    from scipy.interpolate import CubicSpline
    if x is None:
        x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, n)))
        x = x * (sec["upper"][-1, 0] - sec["upper"][0, 0]) + sec["upper"][0, 0]
    zu = CubicSpline(sec["upper"][:, 0], sec["upper"][:, 1])(x)
    zl = CubicSpline(sec["lower"][:, 0], sec["lower"][:, 1])(x)
    return x, zu, zl


def section_properties(sec, n=2000):
    """t/c max, x(t max), max camber, x(max camber) from the unit-chord section."""
    x, zu, zl = resample_surfaces(sec, n=n)
    t = zu - zl
    cam = 0.5 * (zu + zl)
    it, ic = int(np.argmax(t)), int(np.argmax(np.abs(cam)))
    return {"t_over_c_max": float(t[it]), "x_tmax": float(x[it]),
            "camber_max": float(cam[ic]), "x_cmax": float(x[ic]),
            "te_thickness": float(t[-1])}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("vsp3", type=Path)
    ap.add_argument("--geom", default=None, help="wing geom name (e.g. cruise_wing)")
    ap.add_argument("--list", action="store_true", help="list wing geoms and exit")
    ap.add_argument("--csv", type=Path, help="write section table CSV")
    ap.add_argument("--json", type=Path, help="write full extraction JSON")
    ap.add_argument("--dump-dat", type=Path, metavar="DIR",
                    help="write per-station Selig .dat files (unit chord)")
    args = ap.parse_args(argv)

    if args.list:
        for nm in list_wing_geoms(args.vsp3):
            print(nm)
        return 0

    wing = extract_wing(args.vsp3, args.geom)
    rows = []
    for s in wing["sections"]:
        props = section_properties(s) if s["upper"] is not None else {}
        rows.append({"index": s["index"], "eta": s["eta"], "y": s["y"],
                     "chord": s["chord"], "twist": s["twist"],
                     "sweep": s["sweep"], "dihedral": s["dihedral"],
                     "airfoil_name": s["airfoil_name"], **props})

    hdr = ["index", "eta", "y", "chord", "twist", "sweep", "dihedral",
           "t_over_c_max", "x_tmax", "camber_max", "x_cmax", "te_thickness",
           "airfoil_name"]
    print(f"# {wing['file']} geom={wing['geom']} (native model units)")
    print("# planform:", json.dumps(wing["planform"]))
    print(",".join(hdr))
    for r in rows:
        print(",".join("" if r.get(k) is None else
                       (r[k] if isinstance(r[k], str) else f"{r[k]:.10g}")
                       for k in hdr))

    if args.csv:
        with args.csv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=hdr, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    if args.json:
        blob = {"file": wing["file"], "geom": wing["geom"],
                "planform": wing["planform"],
                "sections": [{**{k: v for k, v in s.items()
                                 if k not in ("upper", "lower")},
                              "upper": s["upper"].tolist() if s["upper"] is not None else None,
                              "lower": s["lower"].tolist() if s["lower"] is not None else None}
                             for s in wing["sections"]]}
        args.json.write_text(json.dumps(blob, indent=1) + "\n")
    if args.dump_dat:
        args.dump_dat.mkdir(parents=True, exist_ok=True)
        for s in wing["sections"]:
            if s["upper"] is None:
                continue
            contour = np.vstack([s["upper"][::-1], s["lower"][1:]])
            tag = f"{Path(wing['file']).stem}_{wing['geom']}_eta{s['eta']:.4f}"
            out = args.dump_dat / (tag.replace(".", "p") + ".dat")
            with out.open("w") as f:
                f.write(f"{wing['geom']} eta {s['eta']:.6f} "
                        f"chord {s['chord']:.6f} twist {s['twist']:.5f}\n")
                for px, pz in contour:
                    f.write(f" {px: .9f}  {pz: .9f}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
