#!/usr/bin/env python3
"""Trefftz-plane induced drag for every trimmed case, using the EXISTING instrument.

WHY THIS EXISTS. scripts/case_derived_quantities.py already computes C_Di and span
efficiency on the sampled Trefftz planes, and the instrument behind it was validated
against closed-form elliptical and rectangular wings in ARG-174. But NOT ONE of the
eighteen trimmed cases carries the result: every trim's derived record reads

    "available": false, "note": "no C_L available, so e cannot be formed"

That is a gate discarding a quantity it could have computed. C_Di does NOT depend on
C_L; only the span efficiency e = C_L^2 / (pi AR C_Di) does. So a missing C_L
suppressed BOTH, and the trims are the only quotable cases in the campaign.

C_L is now known for all eighteen to seven figures, from each case's own convergence
marker via data/rans_forces.json. This script supplies it and runs the existing
span_efficiency() unchanged. IT REIMPLEMENTS NOTHING: duplicating a validated
instrument is how you end up with two answers and no way to choose between them.

WHAT THE NUMBER IS FOR. Our surface integral gives TOTAL drag: induced + form +
viscous + wave. The VSPAERO study reports C_Diw, induced ALONE, because a
vortex-lattice method cannot compute the rest. Setting one against the other is a
category error, and it once produced an outright SIGN VIOLATION (ARG-103): VLM
induced 98.79 ct against a RANS induced-PLUS-form 88.29 ct, where ours contains
strictly more terms and cannot be smaller. A Trefftz plane isolates OUR induced drag
under the VLM's own definition so the comparison becomes like-for-like.

WHAT IT IS NOT FOR. It does not replace, check or improve the measured total drag.
That comes from the wing-surface integral and is not in question here.

TWO CAVEATS THAT TRAVEL WITH EVERY NUMBER (ARG-103 item 6, and the later withdrawal):
1. THE STATION SPREAD IS THE UNCERTAINTY, NOT SCATTER TO AVERAGE AWAY. Induced drag
   on a wake plane has no plateau: diffusion biases the far planes DOWN while the
   near plane has not finished rolling up and is biased the other way. The truth is
   BRACKETED, not approached from one side, so the nearest station is the lowest e
   rather than the most correct one.
2. e IS NOT USABLE AS AN ABSOLUTE, only as a same-station difference. The measured
   station spread is of order 11% against a morphing effect of order 0.4%.
And it is the INDUCED-DRAG e, not the Oswald e0; quoting it against the 0.75-0.85
airliner band is the comparison that has already been made in error once.

Usage:
    python3 scripts/trefftz_trims.py [--json out.json]
"""

import argparse
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from case_derived_quantities import span_efficiency, AR_FULL  # noqa: E402

RANS = REPO / "docs/report/all_geometry_2026-09-15/data/rans_forces.json"
CONDS = ["condition_CR", "early_cruise", "late_cruise"]
FAMS = ["B", "C", "M", "W", "H", "F"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    d = json.loads(RANS.read_text())
    rows, missing = [], []

    for cond in CONDS:
        tr = d["trims"][cond]
        frame = tr["coefficient_frame"]
        for fam in FAMS:
            g = tr["geometries"][fam]
            case = g["trim_case"]
            # Every input is read from the harvest, not typed here.
            meta = {
                "U_inf": frame["magUInf_m_s"],
                "alpha_deg": g["alpha_deg"],
                "Aref": frame["Aref_m2"],
            }
            se = span_efficiency(case, meta, g["Cl"])
            if not se.get("available"):
                missing.append((cond, fam, case, se.get("note", "?")))
                continue
            rows.append(dict(condition=cond, family=fam, case=case,
                             CL=g["Cl"], alpha_deg=g["alpha_deg"],
                             CD_total_counts=g["Cd_counts"],
                             dCD_total_counts=g["delta_Cd_counts_vs_baseline"],
                             span_efficiency=se))

    # ---- report -----------------------------------------------------------
    print("Trefftz-plane induced drag at every trim, four stations each.")
    print("AR (full wing, DSO basis) = %.4f\n" % AR_FULL)
    hdr = "%-14s %-3s %8s %10s %10s %10s %10s %8s %8s"
    print(hdr % ("case", "fam", "CL", "CDi_x1", "CDi_x2", "CDi_x3", "CDi_x4", "e_x1", "e_x4"))
    print(hdr % ("", "", "", "[ct]", "[ct]", "[ct]", "[ct]", "", ""))
    for cond in CONDS:
        sub = [r for r in rows if r["condition"] == cond]
        if not sub:
            continue
        print("-- %s --" % cond)
        for r in sub:
            p = r["span_efficiency"]["planes"]
            cd = [q["CDi_counts"] for q in p] + [float("nan")] * (4 - len(p))
            es = [q["span_efficiency"] for q in p] + [None] * (4 - len(p))
            print(hdr % (r["case"], r["family"], "%.6f" % r["CL"],
                         "%.3f" % cd[0], "%.3f" % cd[1], "%.3f" % cd[2], "%.3f" % cd[3],
                         "%.4f" % es[0] if es[0] else "-",
                         "%.4f" % es[3] if es[3] else "-"))

    # ---- induced-drag deltas, formed WITHIN one condition and ONE station --
    # D076: never across operating points. And never across stations, because the
    # station spread is the uncertainty (see the caveats above), so a delta taken
    # between two different stations would be reading that uncertainty as signal.
    print("\nInduced-drag increment vs that condition's own B trim, PER STATION.")
    print("%-14s %-3s %10s %10s %10s %10s   %12s" %
          ("case", "fam", "dCDi_x1", "dCDi_x2", "dCDi_x3", "dCDi_x4", "dCD_total"))
    deltas = {}
    for cond in CONDS:
        sub = {r["family"]: r for r in rows if r["condition"] == cond}
        if "B" not in sub:
            print("-- %s -- no baseline trim with Trefftz data, no delta formed" % cond)
            continue
        base = [q["CDi_counts"] for q in sub["B"]["span_efficiency"]["planes"]]
        print("-- %s --" % cond)
        for fam in FAMS:
            if fam not in sub:
                continue
            p = [q["CDi_counts"] for q in sub[fam]["span_efficiency"]["planes"]]
            dd = [p[i] - base[i] for i in range(min(len(p), len(base)))]
            deltas[(cond, fam)] = dd
            cells = ["%+10.3f" % v for v in dd] + ["%10s" % "-"] * (4 - len(dd))
            print("%-14s %-3s %s   %+12.3f"
                  % (sub[fam]["case"], fam, " ".join(cells),
                     sub[fam]["dCD_total_counts"]))

    if missing:
        print("\nNO TREFFTZ RESULT (%d of %d):" % (len(missing), len(CONDS) * len(FAMS)))
        for cond, fam, case, note in missing:
            print("   %-14s %-3s %-14s %s" % (case, fam, cond, note))
    print("\nCOVERAGE: %d computed + %d missing = %d of %d trims"
          % (len(rows), len(missing), len(rows) + len(missing), len(CONDS) * len(FAMS)))

    if a.json:
        out = {
            "_what": "Trefftz-plane induced drag at every trimmed ARGUS case.",
            "_generator": "scripts/trefftz_trims.py, calling span_efficiency() from "
                          "scripts/case_derived_quantities.py unchanged.",
            "_instrument_validated": "ARG-174, against closed-form elliptical and "
                                     "rectangular wings; best agreement 0.84%.",
            "_caveats": [
                "The station SPREAD is the uncertainty, not scatter to average away. "
                "Diffusion biases far planes down; the near plane has not finished "
                "rolling up and is biased the other way. The truth is BRACKETED.",
                "e is the INDUCED-DRAG e, not the Oswald e0, and is usable only as a "
                "same-station difference, not as an absolute.",
                "Deltas are formed within one condition AND one station only.",
            ],
            "AR_full_DSO": AR_FULL,
            "coverage": {"computed": len(rows), "missing": len(missing),
                         "total": len(CONDS) * len(FAMS)},
            "missing": [{"condition": c, "family": f, "case": k, "note": n}
                        for c, f, k, n in missing],
            "trims": rows,
            "dCDi_counts_by_station": {"%s/%s" % k: v for k, v in deltas.items()},
        }
        Path(a.json).write_text(json.dumps(out, indent=1, sort_keys=True))
        print("wrote %s" % a.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
