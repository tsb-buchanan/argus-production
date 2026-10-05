#!/usr/bin/env python3
"""measure_lineage_pair.py: the measured difference between the UNCORRECTED and CORRECTED
baseline lineages, station by station.

    python3 scripts/measure_lineage_pair.py > docs/report/.../data/lineage_measured.json

WHY THIS EXISTS. The report asserted "the two surfaces differ by at most 0.18 per cent of
chord, by zero inboard of eta = 0.4". Neither number came from a file: 0.18 appears nowhere
in any measurement, and the difference is not zero inboard of 0.4, it is zero inboard of
0.450 because 0.450 is a DEFINING STATION. Replacing one unsourced number with another from
a working note would repeat the defect, so this measures the pair and the report cites it.

THREE NUMBERS WERE IN CIRCULATION AND THEY ARE NOT THE SAME QUANTITY (D065, the frame rule
on a geometric quantity).
  1. 0.18 per cent of chord -- the report's, unsourced, and wrong for any of the three.
  2. 0.407 percentage points of t/c -- correct, and it is THIS measurement, at eta 0.515.
  3. 0.495 percentage points of t/c at eta 0.7100 -- the registry's MANDATORY RESULT LABEL,
     and a DIFFERENT QUANTITY: it is the as-delivered geometry's deviation from its own
     intended loft, caused by the 9-to-19 aerofoil refinement. This script returns EXACTLY
     0.0000 at eta 0.7100, because 0.7100 is a station both files define, so the two
     deliveries agree there precisely while both differ from the intended loft. A t/c
     deviation is meaningless until you say deviation FROM WHAT.

THE STAIRCASE IS THE RESULT, not an artefact: the uncorrected file has 9 distinct aerofoils
across 19 stations, so interior stations INHERIT a neighbour's section. The difference is
therefore zero at every shared station and grows between them.

THE 400-POINT GRID IS CHECKED, NOT ASSUMED, and that check is why this note exists. A uniform
grid in x/c puts only about four points inboard of x/c 0.01, so it CAN step straight over a
leading-edge difference: measuring the uncorrected file's OWN eta 0.8520 -> 0.8820 station step
gives 6.9e-04 c on uniform-400 and 2.3e-03 c on a 20k or cosine grid, where the maximum sits at
x/c 0.001. That is a factor of 3.3 hidden by grid choice alone.
It does NOT happen here: the lineage difference reads 0.2089 %c at eta 0.5150 on uniform-400,
uniform-20k and cosine-4000 alike, to four decimals, because its maximum is at x/c ~ 0.345 and
is well resolved. Only the reported x/c location moves, 0.3559 to 0.3453. d(t/c) carries no
grid at all, being a property of each section. Anyone changing the grid should re-run that
comparison rather than assume the insensitivity carries over to a different pair of surfaces.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vsp3_sections import extract_wing, section_properties

REPO = Path(__file__).resolve().parents[1]
UNC = REPO / "geometry/source/baseline_wing_only_refined.vsp3"
COR = REPO / "dso_reference/2026-07-29_correction_pair/baseline_corrected/baseline_corrected.vsp3"


def distinct(w):
    return len({(s["upper"].tobytes(), s["lower"].tobytes()) for s in w["sections"]})


def main():
    a, b = extract_wing(str(UNC)), extract_wing(str(COR))
    ea = [s["eta"] for s in a["sections"]]
    eb = [s["eta"] for s in b["sections"]]
    if not np.allclose(ea, eb, atol=1e-9):
        sys.exit("FATAL: the two files do not share station eta; a per-station difference "
                 "would be comparing different spanwise locations.")
    rows = []
    for sa, sb in zip(a["sections"], b["sections"]):
        pa, pb = section_properties(sa), section_properties(sb)
        # max |dz| on a COMMON chordwise grid. The two files carry different abscissae, and
        # differencing raw ordinates would compare unrelated points (the .vsp3 resample trap).
        g = np.linspace(0.0, 1.0, 400)
        dz = max(float(np.max(np.abs(np.interp(g, sa[k][:, 0], sa[k][:, 1])
                                     - np.interp(g, sb[k][:, 0], sb[k][:, 1]))))
                 for k in ("upper", "lower"))
        rows.append(dict(eta=round(sa["eta"], 4),
                         dtc_pt=100.0 * (pb["t_over_c_max"] - pa["t_over_c_max"]),
                         dcam_pc=100.0 * (pb["camber_max"] - pa["camber_max"]),
                         dz_pc=100.0 * dz))
    nz = [r["eta"] for r in rows if r["dz_pc"] > 1e-9]
    out = dict(
        uncorrected=str(UNC.relative_to(REPO)), corrected=str(COR.relative_to(REPO)),
        n_stations=len(rows),
        distinct_aerofoils=dict(uncorrected=distinct(a), corrected=distinct(b)),
        worst_dtc_pt=max(rows, key=lambda r: abs(r["dtc_pt"])),
        worst_dcam_pc=max(rows, key=lambda r: abs(r["dcam_pc"])),
        worst_dz_pc=max(rows, key=lambda r: r["dz_pc"]),
        first_eta_that_differs=(min(nz) if nz else None),
        last_eta_identical_inboard=max([r["eta"] for r in rows
                                        if r["eta"] < (min(nz) if nz else 9)], default=None),
        stations=rows)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
