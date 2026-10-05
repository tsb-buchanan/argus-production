#!/usr/bin/env python3
"""Update the report's Condition CR entries to the current meshes, and add induced drag.

WHAT THIS CHANGES. The four Condition CR trim cases (SWB, SWW, SWH, SWF) are replaced with
the values measured on the September 2026 meshes. Every other case in the file, including
all of C and M and all twelve cruise cases, is left exactly as it was.

WHY C AND M ARE NOT TOUCHED. Their meshes match the baseline's to 0.10% on level-7 cells
and 0.15% on level-8, so whatever grid error they carry is common to both sides of C-B and
M-B and cancels in the difference (D071). Re-running them would move both sides together
and leave the deltas where they are.

WHAT IS ADDED, AND IT IS NEW RATHER THAN REVISED. No induced-drag or span-efficiency field
existed in this file before. The Trefftz planes now supply C_Di at four downstream stations
and the span efficiency formed on it. That is what lets Section 5 compare against the
vortex-lattice study on the quantity the study actually ranks on, which the report
previously recorded as not possible.

PROVENANCE IS PRESERVED, NOT OVERWRITTEN. The previous entry for each replaced case is kept
verbatim under `superseded_entry`, and the file is backed up before any write. Nothing is
destroyed; a reader can always recover what the number used to be.

C_Di IS SCALED TO THE TARGET C_L by (C_L,target/C_L,achieved)^2 wherever the case is not
directly trimmed, and the field `cdi_basis` records which of the two it is. Span efficiency
is INVARIANT under that scaling, since e = C_L^2/(pi AR C_Di) and C_Di goes as C_L^2, so
the e values stand either way.
"""
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

# THE DATA ROOT IS FOUND, NOT ASSUMED. This was parent.parent, which is the report
# root only when data/, fig/ and scripts/ are siblings. The delivery repo lifts data/
# to its own root, so the assumption held in one layout and failed silently in the
# other. data_root() walks up for the directory whose data/ carries rans_forces.json
# and RAISES if there is none, rather than returning a plausible wrong tree.
# In the working repo and in an assembled build directory it returns exactly what
# parent.parent returned, so this changes no behaviour there.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_paths import data_root  # noqa: E402
HERE = data_root()
# BOTH FILES, AND THAT IS THE POINT. Every figure generator reads rans_forces.json while
# the report's own table captions cite rans_forces_current.json. Updating only the one
# named "current" left every figure silently plotting the previous numbers, which is the
# ARG-192 shape exactly: a file carrying the word "current" in its name that nothing
# current actually reads. They are written together so they cannot disagree.
DATA_FILES = [HERE / "data/rans_forces.json", HERE / "data/rans_forces_current.json"]
TARGET_CL = 0.428277635108
STATIONS = [2.394, 2.626, 2.916, 3.206]

# Measured 2026-09-21. Sources: each case's own forceCoeffs convergence window and its
# Trefftz planes through span_efficiency() in scripts/case_derived_quantities.py.
#   status "directly_trimmed": the case was re-trimmed and converged at the target.
#   status "CL_corrected":     the untrimmed leg, with C_D on the geometry's own published
#                              alpha-bracket slope and C_Di scaled by (CL_t/CL_a)^2.
R2 = {
    "SWB_trim": dict(status="directly_trimmed", Cl=0.4282682, Cd_counts=188.031,
                     alpha_deg=1.747608659, cdi=[61.834, 58.673, 56.936, 55.984],
                     e_x4=0.9673, decay_pct=9.46, drift_ct=0.0009, cells=99111506),
    "SWW_trim": dict(status="directly_trimmed", Cl=0.4282739, Cd_counts=187.163,
                     alpha_deg=1.520544715, cdi=[60.956, 57.738, 55.995, 55.055],
                     e_x4=0.9837, decay_pct=9.68, drift_ct=0.0055, cells=99203781),
    "SWH_trim": dict(status="CL_corrected", Cl=0.425731210, Cd_counts=188.065,
                     alpha_deg=1.522701, cdi=[61.581, 58.200, 56.400, 55.425],
                     e_x4=0.9771, decay_pct=10.00, drift_ct=None, cells=99121164),
    "SWF_trim": dict(status="CL_corrected", Cl=0.425774800, Cd_counts=187.257,
                     alpha_deg=1.446357, cdi=[61.015, 57.677, 55.895, 54.937],
                     e_x4=0.9858, decay_pct=9.96, drift_ct=None, cells=99150000),
}


def patch(DATA):
    d = json.loads(DATA.read_text())
    cases = d["cases"]

    missing = [k for k in R2 if k not in cases]
    if missing:
        sys.exit("cases absent from the data file, refusing to invent them: %s" % missing)

    bak = DATA.with_suffix(".json.bak_%s" % datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    shutil.copy2(DATA, bak)

    for name, new in R2.items():
        old = dict(cases[name])
        e = cases[name]
        e["superseded_entry"] = old            # provenance: never destroy the prior value
        e["Cl"] = new["Cl"]
        e["Cd_counts"] = new["Cd_counts"]
        e["Cd"] = new["Cd_counts"] * 1e-4
        e["alpha_deg"] = new["alpha_deg"]
        e["L_over_D"] = new["Cl"] / (new["Cd_counts"] * 1e-4)
        e["mesh_cells"] = new["cells"]
        e["trim_status"] = new["status"]
        e["CDi_counts_by_station"] = dict(zip(
            ["x1", "x2", "x3", "x4"], new["cdi"]))
        e["CDi_station_x_m"] = STATIONS
        e["CDi_counts_x4"] = new["cdi"][3]
        e["span_efficiency_x4"] = new["e_x4"]
        e["CDi_decay_x1_to_x4_pct"] = new["decay_pct"]
        e["CDi_drift_last_two_writes_ct"] = new["drift_ct"]
        e["cdi_basis"] = ("C_Di as measured at the trimmed condition"
                          if new["status"] == "directly_trimmed" else
                          "C_Di scaled to target C_L by (CL_t/CL_a)^2; e is invariant "
                          "under that scaling")
        e["measured"] = "2026-09-21"
    d["condition_CR_update_2026_09_21"] = {
        "what": "Condition CR trim cases for B, W, H and F replaced with values measured "
                "on the September 2026 meshes; induced drag and span efficiency added.",
        "untouched": "C and M, and every cruise case. C and M meshes match the baseline's "
                     "to 0.10% on level-7 cells, so their deltas are unaffected.",
        "target_CL": TARGET_CL,
        "stations_m": STATIONS,
        "backup": bak.name,
    }
    DATA.write_text(json.dumps(d, indent=1))
    print("  backed up -> %s" % bak.name)
    print("  %-10s %-18s %9s %9s %9s %8s" %
          ("case", "status", "Cd ct", "CDi x4", "e_x4", "cells"))
    for name, new in R2.items():
        e = cases[name]
        print("  %-10s %-18s %9.3f %9.3f %9.4f %8.1fM" %
              (name, e["trim_status"], e["Cd_counts"], e["CDi_counts_x4"],
               e["span_efficiency_x4"], e["mesh_cells"] / 1e6))
    print("  untouched: C, M and all %d cruise cases"
          % sum(1 for k in cases if k.startswith(("CMP", "LC"))))


# The trims block, which is what every figure generator actually plots. Values measured
# 2026-09-21; each geometry stores its OWN delta and the baseline that delta was formed
# against, so the file is self-describing rather than relying on one top-level baseline.
TRIMS_CR = {
    "B": dict(Cd=188.031, Cl=0.4282682, alpha=1.747608659, delta=0.0,     base=188.031),
    "W": dict(Cd=187.163, Cl=0.4282739, alpha=1.520544715, delta=-0.868,  base=188.031),
    "H": dict(Cd=188.065, Cl=0.425731210, alpha=1.522701,  delta=+0.034,  base=188.031),
    "F": dict(Cd=187.257, Cl=0.425774800, alpha=1.446357,  delta=-0.774,  base=188.031),
    # C and M are NOT re-run and their deltas stand unchanged. Their meshes match the
    # baseline's to 0.10% on level-7 cells, so the grid error is common to both sides of
    # C-B and M-B and cancels (D071). Each records the baseline it was formed against,
    # which is the PUBLISHED baseline, because that is the arithmetic that is true.
    "C": dict(Cd=190.857, Cl=0.4282704, alpha=1.2633719979465325, delta=-0.149, base=191.00600092),
    "M": dict(Cd=191.246, Cl=0.4282657, alpha=1.252018001085504,  delta=+0.240, base=191.00600092),
}


def patch_trims(DATA):
    d = json.loads(DATA.read_text())
    cr = d["trims"]["condition_CR"]
    for g, v in TRIMS_CR.items():
        e = cr["geometries"][g]
        e.setdefault("superseded_trim_entry", {k: e.get(k) for k in
                     ("Cd_counts", "Cl", "alpha_deg", "L_over_D",
                      "delta_Cd_counts_vs_baseline")})
        e["Cd_counts"] = v["Cd"]
        e["Cl"] = v["Cl"]
        e["alpha_deg"] = v["alpha"]
        e["L_over_D"] = v["Cl"] / (v["Cd"] * 1e-4)
        e["delta_Cd_counts_vs_baseline"] = v["delta"]
        e["baseline_Cd_counts_used"] = v["base"]
        if v["base"] == 188.031:
            e["source_Cd"] = ("/home/scratch/$USER/argus/solve_r2/SW%s_trim/.converged"
                              % g)
    cr["baseline_Cd_counts"] = 188.031
    cr["delta_convention"] = (
        "delta_Cd_counts = Cd(geometry) - Cd(baseline), same condition, same Aref, same "
        "target CL. Negative is a drag saving. Each geometry records "
        "baseline_Cd_counts_used, the baseline its delta was actually formed against.")
    DATA.write_text(json.dumps(d, indent=1))
    print("  trims/condition_CR patched:")
    for g in ("B", "C", "M", "W", "H", "F"):
        e = cr["geometries"][g]
        print("    %s Cd %9.3f  delta %+7.3f  (baseline %9.3f)"
              % (g, e["Cd_counts"], e["delta_Cd_counts_vs_baseline"],
                 e["baseline_Cd_counts_used"]))


def main():
    for f in DATA_FILES:
        if not f.exists():
            sys.exit("missing %s" % f)
        print("== %s" % f.name)
        patch(f)
        patch_trims(f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
