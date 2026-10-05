#!/usr/bin/env python3
"""parse_vspaero_lod.py: VSPAERO .lod -> machine-readable strip loads.

DSO reference package (dso_reference/, Zheng 2026-07-24). Extracts the header
reference block and the final-iteration strip table into:
1. <out>/<case>_strips.csv: eta (strip center, Yavg/(Bref/2)), dSpan_ft,
   chord_ft, cl, cd, cl_induced, cd_induced, plus raw columns.
2. <out>/<case>_strips_binned20.csv: re-binned to the wing3d forceCoeffs
   binData 20-strip convention (equal spanwise bins over the half-span,
   area-weighted) for the eventual RANS overlay.
3. <out>/<case>_header.json: reference quantities as-run (NOTE: Sref/Cref are
   the geometry-derived values 13.3572 / 1.29251 ft, NOT the TP-1580
   reference trapezoid 12.0 / 1.11667 ft; see the D017 gate entry).

Dimensionless columns are the agreed comparison basis (Zheng README);
dimensional columns in the package predate his unit audit and are not used.

Usage: python3 scripts/parse_vspaero_lod.py --lod <file.lod> --out results/dso_reference
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]


def parse_lod(path):
    header = {}
    strips = []
    cols = None
    for ln in Path(path).read_text().splitlines():
        t = ln.split()
        if not t:
            continue
        if len(t) == 3 and t[0].endswith("_") and t[2] != "":
            try:
                header[t[0].rstrip("_")] = float(t[1])
                continue
            except ValueError:
                pass
        if t[0] == "Iter":
            cols = t
            continue
        if cols and t[0].isdigit():
            try:
                strips.append([float(v) for v in t[: len(cols)]])
            except ValueError:
                pass
    arr = np.array(strips)
    # keep the final iteration only
    last_iter = arr[:, 0].max()
    arr = arr[arr[:, 0] == last_iter]
    return header, cols, arr


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lod", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nbins", type=int, default=20)
    args = ap.parse_args()

    header, cols, arr = parse_lod(args.lod)
    out = REPO / args.out
    out.mkdir(parents=True, exist_ok=True)
    case = args.lod.stem

    ci = {c: i for i, c in enumerate(cols)}
    half_span = header["Bref"] / 2.0
    eta = arr[:, ci["Yavg"]] / half_span
    rows = []
    for k in range(len(arr)):
        rows.append({
            "strip": int(arr[k, ci["TrailVort"]]),
            "eta": round(float(eta[k]), 6),
            "dSpan_ft": float(arr[k, ci["dSpan"]]),
            "chord_ft": float(arr[k, ci["Chord"]]),
            "dArea_ft2": float(arr[k, ci["dArea"]]),
            "cl": float(arr[k, ci["Cl"]]),
            "cd": float(arr[k, ci["Cd"]]),
            "cl_induced": float(arr[k, ci["Cli"]]),
            "cd_induced": float(arr[k, ci["Cdi"]]),
        })
    with (out / f"{case}_strips.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # area-weighted rebin to n equal spanwise bins (wing3d binData convention)
    edges = np.linspace(0.0, 1.0, args.nbins + 1)
    binned = []
    for b in range(args.nbins):
        m = (eta >= edges[b]) & (eta < edges[b + 1])
        if not m.any():
            binned.append({"bin": b + 1, "eta_lo": edges[b], "eta_hi": edges[b + 1],
                           "cl": "", "cd": "", "cl_induced": "", "cd_induced": ""})
            continue
        wgt = arr[m, ci["dArea"]]
        binned.append({
            "bin": b + 1, "eta_lo": round(float(edges[b]), 3),
            "eta_hi": round(float(edges[b + 1]), 3),
            "cl": float(np.average(arr[m, ci["Cl"]], weights=wgt)),
            "cd": float(np.average(arr[m, ci["Cd"]], weights=wgt)),
            "cl_induced": float(np.average(arr[m, ci["Cli"]], weights=wgt)),
            "cd_induced": float(np.average(arr[m, ci["Cdi"]], weights=wgt)),
        })
    with (out / f"{case}_strips_binned20.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(binned[0].keys()))
        w.writeheader(); w.writerows(binned)

    (out / f"{case}_header.json").write_text(json.dumps(header, indent=2) + "\n")
    print(f"{case}: {len(rows)} strips (iter {int(arr[0,0])}), Sref {header['Sref']}, "
          f"Cref {header['Cref']}, AoA {header['AoA']}; wrote strips + binned20 + header")


if __name__ == "__main__":
    main()
