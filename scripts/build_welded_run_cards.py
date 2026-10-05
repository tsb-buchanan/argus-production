#!/usr/bin/env python3
"""Assemble a schema-valid run card for each of the twelve welded cruise trims (ARG-196).

usage: build_welded_run_cards.py [<C> ...]    (C in CMPB ... LCW; default all twelve)

The sibling of build_r2_run_cards.py, and the same split: facts are READ from the case, the
card is assembled and VALIDATED here, where the schema lives. Two sources:
  1. LOCAL, the hash-checked pull of each case (scripts/pull_welded_postpro.sh) under
     results/postpro_latest/<C>_welded/: convergence marker, trim history, controlDict, the
     argusForces / argusYPlus / argusResiduals / forceCoeffs histories, MESH_INFO.txt.
  2. REMOTE, one read-only ssh per case (argv form, no local shell, so nothing expands here):
     derive_run_card_fields.py on the case itself (solver, turbulence model, wall treatment,
     convergence criteria, OpenFOAM version), CONDITION.json (rho, model-scale mu), both
     decomposeParDicts, the solve case's own wing.stl checksum, the run date (mtime of the
     newest leg's foamRun.log) and the wall clock of every leg from its own PBS stdout.
data/welded_mesh_facts.tsv supplies the checkMesh and layer facts, read from the mesh logs.

RECORDED-AT-RUN-TIME fields are measured from the run, never from the card-writing moment:
date from the solver log, wall clock from each leg's own stdout (a PBS walltime kill is read
from its "walltime N exceeded" line, since the leg's own time stamps stop when it is killed).
git_commit is the working-repo HEAD at card writing, as on the r2 cards, and the notes say so.

THE FORCE SPLIT IS DERIVED AND CHECKED. argusForces is dimensional on the compressible path;
pressure and viscous force vectors are averaged over the same 200-sample window as the gate,
projected on the case's own dragDir and divided by q Aref with q = rho U^2 / 2 from the case.
The split must reconstruct the window-mean C_D to 0.5 per cent or the card is refused.
VALIDATE BEFORE WRITE: each card goes to a .tmp, is validated by scripts/validate_run_card.py,
and is renamed only on PASS. Report first, exit after.
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PULL = REPO / "results/postpro_latest"
OUT = REPO / "results/run_cards_welded"
FACTS = REPO / "docs/report/all_geometry_2026-09-15/data/welded_mesh_facts.tsv"
HOST = "$USER@hpc12.tudelft.net"
WF = "/home/scratch/$USER/argus_wfix"
TOOLS = WF + "/runcard_tools"
ALL = ["CMPB", "CMPC", "CMPF", "CMPH", "CMPM", "CMPW", "LCB", "LCC", "LCF", "LCH", "LCM", "LCW"]
COND = {"CMP": dict(cond="EC", name="early cruise", target=0.52929708745, layer=4.37e-5),
        "LC": dict(cond="LC", name="late cruise", target=0.54411480021, layer=5.61e-5)}
CANDIDATE = {"B": "baseline", "C": "cte_i002_c04", "F": "cfft_b02_c01", "H": "chc_g02_c06",
             "M": "mcv2_i002_c01", "W": "cffw_b01_c01"}
AREF, LREF = 0.620462, 0.393957

REMOTE = r'''
set -u
C="$1"; D=__WF__/solve/${C}_welded; M=__WF__/mesh/${C}_welded
cd "$D" || { echo '{"error":"no case"}'; exit 2; }
echo "@@DERIVED"; python3 __TOOLS__/derive_run_card_fields.py --case "$D" 2>&1
echo "@@CONDITION"; cat "$M/CONDITION.json"
echo "@@MESHDECOMP"; grep -E '^\s*(numberOfSubdomains|method)\s' "$M/system/decomposeParDict"
echo "@@SOLVEDECOMP"; grep -E '^\s*(numberOfSubdomains|method)\s' "$D/system/decomposeParDict"
echo "@@STLSHA"; sha256sum "$D/constant/triSurface/wing.stl" | cut -c1-64
echo "@@RUNDATE"; date -u -r "$(ls -t log/leg*/foamRun.log | head -1)" +%Y-%m-%d
echo "@@LEGS"
for o in $(ls argus-${C}_welded*.o* _killed_walltime_*/argus-${C}_welded*.o* 2>/dev/null); do
  s=$(grep -aoE '^\[[0-9:]+\]' "$o" | head -1 | tr -d '[]')
  e=$(grep -aoE '^\[[0-9:]+\]' "$o" | tail -1 | tr -d '[]')
  k=$(grep -aoE 'walltime [0-9]+ exceeded' "$o" | grep -oE '[0-9]+' | head -1)
  h=$(grep -aoE ' on n12-[0-9]+ ' "$o" | head -1 | tr -d ' ' | sed 's/^on//')
  echo "$o|${s:-}|${e:-}|${k:-}|${h:-}"
done
'''.replace("__WF__", WF).replace("__TOOLS__", TOOLS)


def fail(msg):
    raise SystemExit("REFUSED: " + msg)


def remote(case):
    p = subprocess.run(["ssh", "-o", "BatchMode=yes", HOST, "bash", "-s", "--", case],
                       input=REMOTE, capture_output=True, text=True, timeout=600)
    blocks, cur = {}, None
    for line in p.stdout.splitlines():
        if line.startswith("@@"):
            cur = line[2:]; blocks[cur] = []
        elif cur:
            blocks[cur].append(line)
    need = ["DERIVED", "CONDITION", "MESHDECOMP", "SOLVEDECOMP", "STLSHA", "RUNDATE", "LEGS"]
    miss = [b for b in need if b not in blocks]
    if miss:
        fail("%s: remote blocks missing %s; stderr %s" % (case, miss, p.stderr[-300:]))
    return {k: "\n".join(v) for k, v in blocks.items()}


def decomp(txt):
    n = re.search(r"numberOfSubdomains\s+(\d+)", txt); m = re.search(r"method\s+(\w+)", txt)
    if not n or not m:
        fail("cannot read decomposition from %r" % txt)
    return int(n.group(1)), m.group(1)


def hms(s):
    h, m, x = (int(v) for v in s.split(":"))
    return 3600 * h + 60 * m + x


def legs(txt):
    out = []
    for line in txt.splitlines():
        if not line.strip():
            continue
        f, s, e, killed, host = line.split("|")
        if killed:
            sec, how = int(killed), "PBS walltime kill (from the job's own line)"
        elif s and e:
            sec = hms(e) - hms(s)
            sec, how = (sec + 86400 if sec < 0 else sec), "first to last time stamp"
        else:
            sec, how = 0, "no time stamps: the job stopped before its first stamped step"
        out.append(dict(stdout=f.split("/")[-1], seconds=sec, how=how, host=host or None,
                        killed_run=f.startswith("_killed")))
    if not out:
        fail("no PBS stdout found for the case")
    return out


def table(path, col_need):
    rows, cols = [], None
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            tok = line.lstrip("#").split()
            if col_need in tok:
                cols = tok
            continue
        if line.strip():
            rows.append(line)
    return cols, rows


def latest(case_dir, fo, stem):
    fs = sorted((case_dir / "postProcessing" / fo).glob("*/" + stem), key=lambda p: float(p.parent.name))
    if not fs:
        fail("%s: no %s/%s" % (case_dir.name, fo, stem))
    return fs[-1]


def window(case_dir, n=200):
    cols, rows = table(latest(case_dir, "forceCoeffs1", "forceCoeffs.dat"), "Cd")
    vals = [[float(x) for x in r.split()] for r in rows[-n:]]
    mean = lambda k: sum(v[cols.index(k)] for v in vals) / len(vals)
    return dict(Cd=mean("Cd"), Cl=mean("Cl"), Cm=mean("Cm"), n=len(vals))


def force_split(case_dir, drag, q, cd_total, n=200):
    _, rows = table(latest(case_dir, "argusForces", "forces.dat"), "forces(pressure")
    fp, fv = [0.0] * 3, [0.0] * 3
    use = rows[-n:]
    for r in use:
        nums = [float(x) for x in r.split(None, 1)[1].replace("(", " ").replace(")", " ").split()[:6]]
        fp = [a + b for a, b in zip(fp, nums[0:3])]
        fv = [a + b for a, b in zip(fv, nums[3:6])]
    fp = [x / len(use) for x in fp]; fv = [x / len(use) for x in fv]
    cdp = sum(a * b for a, b in zip(fp, drag)) / (q * AREF)
    cdv = sum(a * b for a, b in zip(fv, drag)) / (q * AREF)
    if abs(cdp + cdv - cd_total) > 0.005 * abs(cd_total):
        fail("%s: split %.7f + %.7f = %.7f does not reconstruct Cd %.7f"
             % (case_dir.name, cdp, cdv, cdp + cdv, cd_total))
    return cdp, cdv


def kv(path):
    return dict(l.split(" ", 1) for l in Path(path).read_text().splitlines() if " " in l)


def build(C):
    pre, L = ("CMP", C[3]) if C.startswith("CMP") else ("LC", C[2])
    c = COND[pre]
    d = PULL / ("%s_welded" % C)
    rem = remote(C)
    derived = json.loads(rem["DERIVED"])
    condj = json.loads(rem["CONDITION"])
    ctl = (d / "controlDict.pulled").read_text()
    g = lambda k: float(re.search(r"^\s*%s\s+([0-9.eE+-]+);" % k, ctl, re.M).group(1))
    rho, U = g("rhoInf"), g("magUInf")
    dd = [float(x) for x in re.search(r"dragDir\s+\(([^)]*)\)", ctl).group(1).split()]
    alpha = math.degrees(math.atan2(dd[2], dd[0]))
    mu = condj["scale_resolution"]["mu_model_scale"]
    if abs(condj["pinned"]["rho"] - rho) > 1e-9 or abs(condj["derived"]["U"] - U) > 1e-5:
        fail("%s: CONDITION.json rho/U %s/%s disagree with controlDict %s/%s"
             % (C, condj["pinned"]["rho"], condj["derived"]["U"], rho, U))
    nu = mu / rho
    marker = (d / "CONVERGED_MARKER.txt").read_text().strip()
    m = re.search(r"REPORT THESE: Cd ([\d.]+)\s+Cl ([\d.]+)", marker)
    cd, cl = float(m.group(1)), float(m.group(2))
    w = window(d)
    if abs(w["Cd"] - cd) > 5e-8 or abs(w["Cl"] - cl) > 5e-8:
        fail("%s: marker not reproduced by forceCoeffs" % C)
    q = 0.5 * rho * U * U
    cdp, cdv = force_split(d, dd, q, cd)
    ycols, yrows = table(latest(d, "argusYPlus", "yPlus.dat"), "patch")
    yw = [r.split() for r in yrows if r.split()[1] == "wing"][-1]
    ymin, ymax, ymean = (float(yw[ycols.index(k)]) for k in ("min", "max", "average"))
    rcols, rrows = table(latest(d, "argusResiduals", "residuals.dat"), "Ux")
    rv = dict(zip(rcols[1:], (float(x) for x in rrows[-1].split()[1:])))
    resid = {"p": rv["p"], "U": max(rv["Ux"], rv["Uy"], rv["Uz"])}
    resid.update({k: rv[k] for k in ("e", "nuTilda") if k in rv})
    th = [l.split("\t") for l in (d / "trim_history.tsv").read_text().splitlines()[1:] if l.strip()]
    hist = [{"alpha_deg": float(r[1]), "CL": float(r[2]), "CD": float(r[3])} for r in th]
    hist.append({"alpha_deg": alpha, "CL": cl, "CD": cd, "Cm": w["Cm"]})
    mi = kv(d / "MESH_INFO.txt")
    if rem["STLSHA"].strip() != mi["wing_stl_sha256"].strip():
        fail("%s: solve-case wing.stl %s differs from the mesh case's %s"
             % (C, rem["STLSHA"].strip()[:16], mi["wing_stl_sha256"][:16]))
    with open(FACTS) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        mf = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
    mf = {r["source_case"]: r for r in mf}["%s_welded" % C]
    mranks, mmethod = decomp(rem["MESHDECOMP"])
    sranks, smethod = decomp(rem["SOLVEDECOMP"])
    lg = legs(rem["LEGS"])
    wall = sum(x["seconds"] for x in lg)
    final_iter = int(float(kv(d / "PULL_INFO.txt")["final_time"]))
    wt = dict(derived["wall_treatment_block"], yplus_min=ymin, yplus_mean=ymean, yplus_max=ymax,
              yplus_convention="cell-centre")
    card = {
        "schema_version": 1,
        "candidate": CANDIDATE[L],
        "case": {"path": "%s/solve/%s_welded" % (WF.replace("$USER", "$USER"), C),
                 "dimensionality": "3D", "condition": c["cond"], "eta": None,
                 "description": "%s trim of geometry %s at M 0.78, trimmed to target C_L by "
                                "rotating the freestream and the force axes together at "
                                "fixed |U|." % (c["name"].capitalize(), L)},
        "geometry": {"file": mi["wing_stl_registered_as"], "sha256": mi["wing_stl_sha256"]},
        "provenance": {
            "git_commit": subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                                         capture_output=True, text=True).stdout.strip(),
            "git_dirty": subprocess.run(["git", "-C", str(REPO), "status", "--porcelain"],
                                        capture_output=True, text=True).stdout.strip() != "",
            "openfoam_version": derived["provenance"]["openfoam_version"],
            "date": rem["RUNDATE"].strip()},
        "mesh": {
            "generator": "snappyHexMesh",
            "cells": int(mf["cells"]),
            "max_skewness": float(mf["skewness_max"]),
            "max_non_orthogonality": float(mf["nonortho_max_deg"]),
            "avg_non_orthogonality": float(mf["nonortho_avg_deg"]),
            "yplus": {"min": ymin, "mean": ymean, "max": ymax},
            "layers": {"n_surface_layers": 5, "first_layer_thickness_m": c["layer"],
                       "expansion_ratio": 1.3},
            "checkMesh_pass": int(mf["checkmesh_failed"]) == 0,
            "mesh_generation": {"ranks": mranks, "method": mmethod, "inferred": False,
                                "inferred_from": "mesh/%s_welded/system/decomposeParDict, read "
                                                 "directly" % C, "machine": "hpc12"}},
        "model": {"solver": derived["model"]["solver"],
                  "solver_module": derived["model"]["solver_module"],
                  "turbulence_model": derived["model"]["turbulence_model"],
                  "wall_treatment": derived["model"]["wall_treatment"],
                  "decision_ref": derived["wall_functions_sanctioned_by"]},
        "bc": {"farfield": "freestreamVelocity / freestreamPressure on inlet, outlet, topBottom "
                           "and outboard, freestream direction rotated to the trim alpha",
               "airfoil": "noSlip on wing, nutUSpaldingWallFunction"},
        "numerics": {"schemes_summary": "steady compressible foamRun fluid, Spalart-Allmaras, "
                                        "recipe numerics of M6_wallModelled_CRUISE",
                     "convergence_criteria": derived["numerics"]["convergence_criteria"]},
        "convergence": {"iterations": final_iter, "final_residuals": resid,
                        "converged": "-> CONVERGED" in marker,
                        "criterion": "200-sample window: mean drift Cd <= 0.05 ct, Cl <= 1.0 ct, "
                                     "Cl span <= 1.0 ct"},
        "conditions": {"U_mag_m_s": U, "alpha_deg": alpha, "nu_m2_s": nu,
                       "Re_ref": U * LREF / nu, "mach_nominal": 0.78, "rho_ref_kg_m3": rho,
                       "reference": {"Aref_m2": AREF, "lRef_m": LREF, "CofR_m": [1.34, 0.0, 0.0],
                                     "convention": "DSO basis, S_ref 1.24092 m2 and c_ref "
                                                   "0.39396 m, b_ref 3.6576 m. Half-wing A_ref "
                                                   "0.620462 m2 = S_ref/2. Model scale with mu "
                                                   "set so Re on c_ref matches aircraft cruise "
                                                   "(CONDITION.json scale_resolution)."}},
        "results": {"CL": cl, "CD": cd, "CD_pressure": cdp, "CD_viscous": cdv, "Cm": w["Cm"],
                    "CD_counts": 1e4 * cd},
        "trim": {"target_CL": c["target"], "tolerance": 1e-4, "alpha_trim_deg": alpha,
                 "history": hist, "converged": abs(cl - c["target"]) < 1e-4},
        "cost": {"wall_clock_s": wall, "n_cores": sranks,
                 "core_hours": round(wall * sranks / 3600.0, 1),
                 "machine": lg[-1]["host"] or "hpc12"},
        "decomposition": {"method": smethod, "n_subdomains": sranks},
        "hpc": {"machine": "hpc12.tudelft.net", "solve_ranks": sranks,
                "solve_decomposition_method": smethod,
                "walltime_requested": "12:00:00 per leg on fpt-medium, 16:00:00 on fpt-large"},
        "wall_treatment": wt,
        "notes": ("Welded cruise trim (ARG-196). checkMesh fails %s checks, of the classes the "
                  "pipeline's mesh check accepts (%s). Trim reached by "
                  "continuation over %d leg(s), history in trim. Wall clock summed over every "
                  "PBS leg of this case, read from each leg's own stdout: %s. git_commit is the "
                  "working-repo HEAD at card writing; the cluster ran the argus-production "
                  "clone. Convergence marker: %s"
                  % (mf["checkmesh_failed"], "listed per mesh in ARG-196", len(hist),
                     "; ".join("%s %ds (%s)%s" % (x["stdout"], x["seconds"], x["how"],
                                                  ", killed run kept in _killed_walltime_*" if x["killed_run"] else "")
                               for x in lg), marker.replace("\n", " | "))),
    }
    return card


def main():
    cases = sys.argv[1:] or ALL
    bad = [c for c in cases if c not in ALL]
    if bad:
        fail("unknown cases %s" % bad)
    OUT.mkdir(parents=True, exist_ok=True)
    written, refused = [], []
    for C in cases:
        try:
            card = build(C)
        except SystemExit as e:
            refused.append((C, str(e))); print("  %-5s REFUSED %s" % (C, e)); continue
        p = OUT / ("%s_welded_trim.json" % C)
        tmp = p.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(card, indent=2) + "\n")
        rc = subprocess.run([sys.executable, str(REPO / "scripts/validate_run_card.py"), str(tmp)],
                            capture_output=True, text=True)
        if rc.returncode == 0:
            tmp.rename(p); written.append(C)
            print("  %-5s VALID   %s  CL %.7f  CD %.7f (p %.7f + v %.7f)  %.1f core-h  %s"
                  % (C, p.name, card["results"]["CL"], card["results"]["CD"],
                     card["results"]["CD_pressure"], card["results"]["CD_viscous"],
                     card["cost"]["core_hours"], card["provenance"]["date"]))
        else:
            tmp.unlink(missing_ok=True); refused.append((C, (rc.stdout + rc.stderr)[-400:]))
            print("  %-5s INVALID, not promoted:\n     %s" % (C, (rc.stdout + rc.stderr)[-400:]))
    print("  ---- %d written, %d refused, of %d ----" % (len(written), len(refused), len(cases)))
    return 1 if refused else 0


if __name__ == "__main__":
    sys.exit(main())
