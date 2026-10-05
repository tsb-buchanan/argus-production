#!/usr/bin/env python3
"""collect_2d_grid.py: harvest an ARGUS 2D O-grid family (grid-convergence study).

Generic over build_2d_case.py cases: run cards, Richardson/GCI for cl and cd,
convergence CSV + JSON + figure. No external reference band (that is the
TMR Layer-1 job); this is discretization-error quantification (D012 ladder:
grid convergence on the EET section).

Usage:
  python3 scripts/collect_2d_grid.py --glob 'cases/baseline/CR2d_grid_*' \
      --out results/eet2d_cr
"""

import argparse
import csv
import datetime
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "validation_data/naca0012_tmr"))
from collect_results import history  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
from derive_run_card_fields import derive as derive_state  # noqa: E402
from atomic_artifact import (write_if_valid, run_card_validator,  # noqa: E402
                             preserve_run_time_fields)


def parse_case(case):
    meta = json.loads((case / "build_meta.json").read_text())
    a = math.radians(meta["alpha_deg"])
    drag = np.array([math.cos(a), 0.0, math.sin(a)])
    qref = 0.5 * 1.225 * meta["U_m_s"] ** 2 * meta["chord_m"]

    fc = np.loadtxt(history(case, "forces", "forceCoeffs*.dat"))
    t, cm, cd, cl = fc[:, 0], fc[:, 1], fc[:, 2], fc[:, 3]
    n10 = max(5, len(t) // 10)
    flat = float((cd[-n10:].max() - cd[-n10:].min()) * 1e4)

    last = history(case, "forcesAll", "forces*.dat")[-1]
    vecs = [np.array([float(v) for v in m.split()])
            for m in re.findall(r"\(([-0-9.eE+ ]+)\)", last)]
    cdp, cdv = float(vecs[0] @ drag / qref), float(vecs[1] @ drag / qref)
    yp = [float(v) for v in history(case, "yPlus", "yPlus*.dat")[-1].split()[2:5]]

    chk = (case / "log.checkMesh").read_text()
    # ACHIEVED residuals from the final solver block (the card used to record the
    # residualControl TARGETS here, which is not a record of what the run did).
    log = (case / "log.simpleFoam").read_text()
    res = {}
    for f in ("p", "Ux", "Uz", "k", "omega"):
        hits = re.findall(rf"Solving for {f}, Initial residual = ([0-9.eE+-]+)", log)
        if hits:
            res[f] = float(hits[-1])
    if "Ux" in res or "Uz" in res:
        res["U"] = max(res.get("Ux", 0.0), res.get("Uz", 0.0))
    # D023 stationarity gate: drift of the MEAN across the last 1000 iterations,
    # measured first half vs second half. Three rejected alternatives, each for a
    # concrete reason: an endpoint difference is blind to oscillation; a
    # peak-to-peak gate is unachievable for a case in a stable numerical limit
    # cycle (CR2d_trimm_eta0p80_fine_am1 holds p2p ~0.011 counts indefinitely while
    # its mean is steady to 5e-4); and comparing against the PRECEDING 1000
    # iterations reaches back into the initial transient on short runs. p2p is
    # recorded separately as the unsteadiness amplitude on the reported values.
    w = t >= t[-1] - 1000
    tw, cw = t[w], cd[w]
    mid = tw[0] + (tw[-1] - tw[0]) / 2
    mean_last = float(cw.mean()) * 1e4
    mean_drift = float(abs(cw[tw >= mid].mean() - cw[tw < mid].mean())) * 1e4
    p2p = float(cw.max() - cw.min()) * 1e4
    return {
        "meta": meta, "cl": float(cl[-1]), "cd": float(cd[-1]), "cm": float(cm[-1]),
        "cdp": cdp, "cdv": cdv, "iterations": int(t[-1]), "flat_counts": flat,
        "mean_drift_counts": mean_drift, "p2p_counts": p2p,
        "mean_last1000_counts": mean_last, "residuals": res, "yplus": yp,
        "cells": int(re.search(r"cells:\s+(\d+)", chk).group(1)),
        "nonorth": float(re.search(r"non-orthogonality Max: ([0-9.]+)", chk).group(1)),
        "skew": float(re.search(r"Max skewness = ([0-9.]+)", chk).group(1)),
        "chk_ok": "Mesh OK" in chk,
        "wall_s": (lambda m: float(m[-1]) if m else 0.0)(
            re.findall(r"ClockTime = ([0-9.]+) s", (case / "log.simpleFoam").read_text())),
    }


def run_card(case, r):
    m = r["meta"]
    st = derive_state(case)
    # GEO-081. Regeneration may recompute anything the schema classes as DERIVED
    # and must preserve everything it classes as RECORDED-AT-RUN-TIME. The class
    # lives in run_card.schema.json, so a field added later inherits the protection
    # without this collector being touched. Patching git_commit alone would have
    # left machine, wall_clock, core_hours and the build waiting to repeat it.
    prior = {}
    if (case / "run_card.json").exists():
        try:
            prior = json.loads((case / "run_card.json").read_text())
        except (ValueError, OSError):
            prior = {}
    git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                         text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=REPO,
                                capture_output=True, text=True).stdout.strip())
    np_ranks = int(re.search(r"numberOfSubdomains (\d+)",
                             (case / "system/decomposeParDict").read_text()).group(1))
    card = {
        "schema_version": 1,
        "candidate": m["candidate"],
        "case": {"path": str(case.relative_to(REPO)), "dimensionality": "2D",
                 "condition": m["condition"], "eta": None,
                 "description": f"2D grid-convergence family, {m['level']} level, "
                                f"alpha {m['alpha_deg']} deg, Condition {m['condition']}"},
        "geometry": {"file": m["geometry"], "sha256": m["geometry_sha256"]},
        "provenance": {"git_commit": git, "git_dirty": dirty,
                       # GEO-080: read from the log's own Build line, not asserted.
                       # The full build hash (st[...]["openfoam_build"]) is checked
                       # by the derivation but not stored: schema v1 forbids the key.
                       "openfoam_version": st["provenance"]["openfoam_version"],
                       # RUN date (mtime of the final solver log), not the card-write
                       # date: regenerating a card must not re-date the run it records
                       "date": datetime.date.fromtimestamp(
                           (case / "log.simpleFoam").stat().st_mtime).isoformat()},
        "mesh": {"generator": m["mesh"] + " (scripts/ogrid2d.py)", "cells": r["cells"],
                 "max_skewness": r["skew"], "max_non_orthogonality": r["nonorth"],
                 "yplus": {"min": r["yplus"][0], "mean": r["yplus"][2], "max": r["yplus"][1]},
                 "checkMesh_pass": r["chk_ok"],
                 # DERIVED FROM THE CASE, NOT INHERITED FROM A PRIOR CARD (GEO-080).
                 # This block and `hpc` below used to exist only in cards written by an
                 # older generator and were carried forward only because nothing
                 # regenerated them. Neither is in the recorded-at-run-time set, so
                 # nothing preserved them either: A CASE WITH NO PRIOR CARD COULD NEVER
                 # PRODUCE A VALID ONE. The generator was correct on every steady state
                 # (regenerating a card that already existed) and wrong on the one
                 # transition that creates the thing it protects, which is the
                 # GEO-089 transition rule in the run-card layer.
                 "mesh_generation": mesh_generation(case)},
        # GEO-080: DERIVED from controlDict + the log's Exec line, from
        # constant/turbulenceProperties, and from the nut wall BC. NOT authored
        # here. Off-recipe values raise rather than being relabelled.
        "model": st["model"],
        "bc": {"farfield": {"U": "freestream (circular, ~500c)", "p": "freestreamPressure",
                            "k": "inletOutlet (condition value)", "omega": "inletOutlet"},
               "airfoil": {"U": "noSlip", "k": "1e-12", "omega": "omegaWallFunction",
                           "nut": "nutLowReWallFunction"}},
        "numerics": {"schemes_summary": "template_case/airfoil2d VERBATIM; turbulent-branch "
                                        "k init (D012 Layer-1 lesson)",
                     # GEO-080: read from system/fvSolution residualControl
                     "convergence_criteria": st["numerics"]["convergence_criteria"]},
        "convergence": {"iterations": r["iterations"],
                        # ACHIEVED initial residuals at the last solver iteration
                        "final_residuals": r["residuals"],
                        # D023: basis is force stationarity, NOT residualControl
                        # (no 2D case in this project meets U 1e-8)
                        "converged": r["mean_drift_counts"] <= 0.01,
                        "cd_mean_drift_1000iter_counts": r["mean_drift_counts"],
                        "cd_p2p_last1000_counts": r["p2p_counts"],
                        "cd_mean_last1000_counts": r["mean_last1000_counts"],
                        "criterion": "D023 (2D): |mean(cd, 2nd half) - mean(cd, 1st "
                                     "half)| over the last 1000 iterations <= 0.01 counts. "
                                     "Reported coefficients are the instantaneous "
                                     "final-write values (keeps CD = CD_pressure + "
                                     "CD_viscous exact); cd_p2p_last1000_counts is the "
                                     "unsteadiness amplitude on them. residualControl "
                                     "targets (p 1e-7, U 1e-8) are NOT met; "
                                     "final_residuals holds ACHIEVED values"},
        "conditions": {"U_mag_m_s": m["U_m_s"], "alpha_deg": m["alpha_deg"],
                       "nu_m2_s": m["nu_m2_s"], "Re_ref": m["Re"],
                       "mach_nominal": m["mach_nominal"], "rho_ref_kg_m3": 1.225,
                       "reference": {"Aref_m2": m["chord_m"], "lRef_m": m["chord_m"],
                                     "CofR_m": [m["chord_m"] * 0.25, 0.0, 0.0],
                                     "convention": "per unit span: Aref = chord x 1 m span (2D)"}},
        "results": {"CL": r["cl"], "CD": r["cd"], "CD_pressure": r["cdp"],
                    "CD_viscous": r["cdv"], "Cm": r["cm"],
                    "cd_range_last10pct_counts": r["flat_counts"]},
        # MANDATORY since schema commit d1abb50, and absent from every card ever
        # written because that commit changed the schema and asserted "every card
        # now names nut, k and omega" without changing any generator. Derived.
        "wall_treatment": dict(st["wall_treatment_block"],
                               yplus_min=r["yplus"][0], yplus_mean=r["yplus"][2],
                               yplus_max=r["yplus"][1],
                               yplus_convention="cell-centre"),
        "trim": None,
        # DERIVED: an HPC block describes a SCHEDULER-DISPATCHED run. These 2D cases are
        # driven by the case's own Allrun on the laptop, so there is no job to describe
        # and the field is null. Presence of a slurm script in the case is what would
        # make it non-null, so this is read from the case rather than assumed.
        "hpc": hpc_block(case),
        # machine is a property of the executing HOST, not of the case, so it is
        # the one field here that stays an input. The log's own Host line is
        # recorded beside it so a reader can compare the two (GEO-080).
        "cost": {"wall_clock_s": r["wall_s"], "n_cores": np_ranks,
                 "core_hours": r["wall_s"] * np_ranks / 3600,
                 "machine": "laptop-WSL2"},
        "notes": "Grid study (D012 ladder); fixed alpha near the operating cl; "
                 "trim study follows separately.",
    }
    kept = preserve_run_time_fields(card, prior)
    # GEO-081: VALIDATE BEFORE WRITE. The card only reaches run_card.json if it
    # passes; a rejected card never exists under its real name, so a validation
    # failure can no longer present as a crash that leaves an invalid artifact.
    write_if_valid(case / "run_card.json", json.dumps(card, indent=2) + "\n",
                   run_card_validator)
    return kept



def mesh_generation(case):
    """The mesh-generation block, DERIVED from what the case actually contains.

    A case carrying no system/snappyHexMeshDict was built by the in-house serial
    structured generator, so its mesh rank count is not merely unknown, it is
    inapplicable. The processor* directories present in a solved case are the SOLVE
    decomposition; inferring mesh ranks from them would give a confidently wrong number,
    which is why the reason is recorded inside the block rather than left to a reader.
    """
    snappy = (case / "system/snappyHexMeshDict").exists()
    if snappy:
        raise SystemExit(
            f"{case}: carries system/snappyHexMeshDict, so it is NOT an in-house 2D "
            "O-grid case and this collector cannot state its mesh rank count. Rank count "
            "changes the mesh parallel snappyHexMesh produces, so guessing it here would "
            "put a wrong number into the recipe half of the card.")
    return {"ranks": None,
            "method": "none (serial in-house structured generator, not snappyHexMesh)",
            "inferred": True,
            # `machine` is NOT recoverable from the case: nothing in a 2D O-grid case
            # records where ogrid2d.py ran. It is marked inferred and the basis is stated
            # here rather than presented as a recorded fact, which is the distinction the
            # schema's own `inferred` flag exists to carry.
            "machine": "laptop-WSL2",
            "inferred_from": "case carries NO system/snappyHexMeshDict, so it was built "
                             "by the in-house serial 2D/O-grid generator (scripts/"
                             "ogrid2d.py). NOTE the processor* directories present in "
                             "this case are the SOLVE decomposition and are NOT evidence "
                             "of mesh rank count; inferring mesh ranks from them would "
                             "give a confidently wrong number. `machine` is likewise "
                             "inferred: the 2D generator records no build host, and every "
                             "2D O-grid in this repo was built on the laptop."}


def hpc_block(case):
    """Null unless the case was dispatched through a scheduler, read from the case."""
    for pat in ("*.slurm", "*.sbatch", "slurm-*.out"):
        hit = sorted(case.glob(pat))
        if hit:
            return {"scheduler": "slurm", "script": hit[0].name}
    return None


def richardson(vals, cells):
    fc_, fm, ff = vals
    r = 0.5 * (math.sqrt(cells[2] / cells[1]) + math.sqrt(cells[1] / cells[0]))
    e21, e32 = fc_ - fm, fm - ff
    if e32 == 0 or e21 / e32 <= 0:
        return {"monotone": False, "value_fine": ff, "note": "non-monotone"}
    p = math.log(e21 / e32) / math.log(r)
    if not (0.5 < p < 4.0):
        return {"monotone": True, "value_fine": ff, "p_observed": p,
                "note": "outside asymptotic range; fine value used"}
    ext = ff + (ff - fm) / (r ** p - 1.0)
    return {"monotone": True, "p_observed": p, "extrapolated": ext,
            "gci_fine_pct": 100 * 1.25 * abs(e32) / (r ** p - 1) / abs(ff),
            "value_fine": ff}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--glob", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    out = REPO / args.out
    out.mkdir(parents=True, exist_ok=True)

    data = {}
    rows = []
    for case in sorted(REPO.glob(args.glob)):
        if not (case / "postProcessing/yPlus").exists():
            print(f"skip {case.name} (incomplete)")
            continue
        r = parse_case(case)
        run_card(case, r)
        data[r["meta"]["level"]] = r
        rows.append({"level": r["meta"]["level"], "cells": r["cells"],
                     "cl": r["cl"], "cd": r["cd"], "cd_pressure": r["cdp"],
                     "cd_viscous": r["cdv"], "cm": r["cm"],
                     "yplus_mean": r["yplus"][2], "yplus_max": r["yplus"][1],
                     "iterations": r["iterations"], "flat_counts": r["flat_counts"],
                     "checkMesh_pass": r["chk_ok"], "wall_clock_s": r["wall_s"]})
        print(f"{case.name}: cells {r['cells']} cl {r['cl']:+.6f} cd {r['cd']:.6f} "
              f"(p {r['cdp']:.6f} v {r['cdv']:.6f}) y+max {r['yplus'][1]:.2f} "
              f"flat {r['flat_counts']:.3f}cts chk {'OK' if r['chk_ok'] else 'FLAGGED'}")

    rows.sort(key=lambda r: r["cells"])
    with (out / "grid_family.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    summary = {}
    if all(l in data for l in ("coarse", "medium", "fine")):
        cells = [data[l]["cells"] for l in ("coarse", "medium", "fine")]
        for q in ("cl", "cd", "cdp", "cdv"):
            summary[q] = richardson([data[l][q] for l in ("coarse", "medium", "fine")], cells)
    (out / "grid_richardson.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=1))



if __name__ == "__main__":
    main()
