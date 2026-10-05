#!/usr/bin/env python3
"""Move the report's early- and late-cruise trims onto the welded meshes (ARG-196).

WHAT THIS CHANGES. The twelve cruise trim records, CMP<L>_trim and LC<L>_trim for L in
B C F H M W, are replaced with the trimmed states of the twelve welded cases
(<prefix><L>_welded, HPC12, September 2026). The logical names are kept, the precedent set
when Condition CR moved to the wake-refined meshes (update_cr_from_r2.py): every figure and
table generator builds "<prefix><L>_trim", and the generation now travels in the record
(mesh_generation "welded") instead of in the name.

WHY ALL TWELVE AND NOT ONLY F, H, W. The published F, H and W cruise meshes were built
from surfaces with an unwelded root seam and were starved of cells; B, C and M were not.
But every welded case is on the 130M-cap recipe, so re-basing only three would set a
welded F beside a published B, and the delta would cross a mesh generation (D071: the
generation then differs between the two sides and does not cancel). All six geometries move
together, so every within-condition delta is formed inside one generation.

WHAT IS ADDED. The welded campaign trimmed by nudging alpha after each converged 2000-
iteration leg. Every pre-nudge leg converged on the same gate as the trim, so each is an
off-trim welded state at the same mesh, condition and recipe. They are added as
<prefix><L>_welded_leg<n> (is_trim_case false) and are the ONLY points the welded
dCd/dCl fits use. The published alpha-bracket cases are left in the file untouched,
labelled mesh_generation "published", and are no longer in any cruise fit, because a fit
through bracket and trim would be a mixed-generation object (the defect
_provenance_notes.mesh_generation_field records for Condition CR).

EVERY NUMBER IS READ FROM A FILE (GEO-080). Inputs are the hash-checked pulls under
results/postpro_latest/<prefix><L>_welded/ (scripts/pull_welded_postpro.sh):
CONVERGED_MARKER.txt, leg_markers/converged.leg*, trim_history.tsv, controlDict.pulled,
MESH_INFO.txt, SOLVER_INFO.txt, postProcessing/forceCoeffs1/*/forceCoeffs.dat.

GATES, all asserted BEFORE anything is written (validate before write):
  1. The marker's window means are reproduced from the raw forceCoeffs history by an
     independent route (last 200 samples of the numerically latest leg) to 5e-8.
  2. alpha from the controlDict dragDir equals alpha from liftDir and equals the
     trim_history final alpha plus its last nudge, to 1e-6 deg.
  3. Aref 0.620462 and lRef 0.393957 (DSO basis, half model) and magUInf per condition.
  4. |C_L - target| < 1e-4 (the project trim tolerance) on every trim.
  5. Partition: 12 trims replaced + every leg record added = every welded marker read.

PROVENANCE IS PRESERVED. Each replaced record keeps its previous content verbatim under
superseded_entry; each replaced trims-block geometry under superseded_trim_entry; each
replaced fit under superseded_fit. Both data files are backed up before the write, and the
write goes to a temporary file that is renamed only after it re-reads cleanly.
"""
import json
import math
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_paths import data_root  # noqa: E402

HERE = data_root()
DATA_FILES = [HERE / "data/rans_forces.json", HERE / "data/rans_forces_current.json"]
REPO = Path(__file__).resolve().parents[4]
PULL = REPO / "results/postpro_latest"

LETTERS = "BCFHMW"
COND = {
    "early_cruise": dict(prefix="CMP", target=0.52929708745, magU=233.934028),
    "late_cruise": dict(prefix="LC", target=0.54411480021, magU=230.501789),
}
AREF, LREF = 0.620462, 0.393957
TOL_CL = 1e-4
# A fit over two points closer than this in C_L is noise-dominated (measured 2026-10-05:
# C at early cruise, span 3.2e-4, gave 78 ct/CL against 300 to 400 for the wide fits).
MIN_FIT_SPAN = 2e-3
CLUSTER = "/home/scratch/$USER/argus_wfix"
MARKER_RE = re.compile(
    r"finalised (\S+) on (\S+), leg (\d+).*?drift Cd ([\d.]+) ct.*?Cl ([\d.]+) ct.*?"
    r"Cl SPAN ([\d.]+) ct.*?-> (\w+).*?REPORT THESE: Cd ([\d.]+)\s+Cl ([\d.]+)", re.S)


def fail(msg):
    sys.exit("REFUSED: " + msg)


def kv(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        if " " in line:
            k, v = line.split(" ", 1)
            out[k] = v.strip()
    return out


def parse_marker(text, where):
    m = MARKER_RE.search(text)
    if not m:
        fail("cannot parse convergence marker %s" % where)
    utc, node, leg, dcd, dcl, span, verdict, cd, cl = m.groups()
    if verdict != "CONVERGED":
        fail("%s verdict is %s" % (where, verdict))
    return dict(finalised_utc=utc, finalised_node=node, leg=int(leg), drift_Cd_ct=float(dcd),
                drift_Cl_ct=float(dcl), Cl_span_ct=float(span), gate_verdict=verdict,
                Cd=float(cd), Cl=float(cl), marker_text=text.strip())


def control(case_dir):
    t = (case_dir / "controlDict.pulled").read_text()
    g = lambda k: re.search(r"^\s*%s\s+([^;]+);" % k, t, re.M).group(1).strip()
    vec = lambda s: [float(v) for v in s.strip("()").split()]
    dd, ld = vec(g("dragDir")), vec(g("liftDir"))
    a_d = math.degrees(math.atan2(dd[2], dd[0]))
    a_l = math.degrees(math.atan2(-ld[0], ld[2]))
    return dict(magUInf=float(g("magUInf")), Aref=float(g("Aref")), lRef=float(g("lRef")),
                endTime=float(g("endTime")), dragDir=dd, alpha_dragDir=a_d, alpha_liftDir=a_l)


def window_means(case_dir, n=200):
    """Independent route to the marker: last n samples of the numerically latest leg."""
    files = sorted((case_dir / "postProcessing/forceCoeffs1").glob("*/forceCoeffs.dat"),
                   key=lambda p: float(p.parent.name))
    rows, cols = [], None
    for line in files[-1].read_text().splitlines():
        if line.startswith("#"):
            tok = line.lstrip("#").split()
            if "Cd" in tok:
                cols = tok
            continue
        p = line.split()
        if cols and len(p) == len(cols):
            rows.append([float(v) for v in p])
    a = np.array(rows[-n:])
    return dict(Cd=a[:, cols.index("Cd")].mean(), Cl=a[:, cols.index("Cl")].mean(),
                t_first=a[0, 0], t_last=a[-1, 0], n=len(a), file=str(files[-1]))


def trim_history(case_dir):
    lines = (case_dir / "trim_history.tsv").read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:] if l.strip()]


def gather():
    """Read and gate every welded case. Returns {(cond, L): {...}}; refuses on any failure."""
    out = {}
    for cond, c in COND.items():
        for L in LETTERS:
            name = "%s%s_welded" % (c["prefix"], L)
            d = PULL / name
            if not d.is_dir():
                fail("%s not pulled (%s)" % (name, d))
            mk = parse_marker((d / "CONVERGED_MARKER.txt").read_text(), name)
            ctl = control(d)
            wm = window_means(d)
            th = trim_history(d)
            mesh, solver = kv(d / "MESH_INFO.txt"), kv(d / "SOLVER_INFO.txt")
            info = kv(d / "PULL_INFO.txt")
            # gate 1: marker reproduced by an independent route
            if abs(wm["Cd"] - mk["Cd"]) > 5e-8 or abs(wm["Cl"] - mk["Cl"]) > 5e-8:
                fail("%s marker Cd/Cl %s/%s not reproduced by forceCoeffs %s/%s"
                     % (name, mk["Cd"], mk["Cl"], wm["Cd"], wm["Cl"]))
            # gate 2: alpha three ways
            a_hist = float(th[-1]["alpha_deg"]) + float(th[-1]["dalpha_next"])
            if not (abs(ctl["alpha_dragDir"] - ctl["alpha_liftDir"]) < 1e-6
                    and abs(ctl["alpha_dragDir"] - a_hist) < 1e-6):
                fail("%s alpha disagrees: dragDir %.9f liftDir %.9f history %.9f"
                     % (name, ctl["alpha_dragDir"], ctl["alpha_liftDir"], a_hist))
            # gate 3: frame
            if (abs(ctl["Aref"] - AREF) > 1e-9 or abs(ctl["lRef"] - LREF) > 1e-9
                    or abs(ctl["magUInf"] - c["magU"]) > 1e-6):
                fail("%s frame Aref %s lRef %s magUInf %s" % (name, ctl["Aref"], ctl["lRef"], ctl["magUInf"]))
            # gate 4: trim tolerance
            if abs(mk["Cl"] - c["target"]) >= TOL_CL:
                fail("%s Cl %s outside 1e-4 of target %s" % (name, mk["Cl"], c["target"]))
            cells = int(re.search(r"nCells:(\d+)", mesh["owner_note"]).group(1))
            legs = []
            for f in sorted((d / "leg_markers").glob("converged.leg*"),
                            key=lambda p: int(re.search(r"leg(\d+)", p.name).group(1))):
                lm = parse_marker(f.read_text(), str(f))
                n = int(re.search(r"leg(\d+)_alpha", f.name).group(1))
                a_name = float(re.search(r"_alpha([\d.]+)$", f.name).group(1))
                row = [r for r in th if int(r["leg"]) == n]
                if len(row) != 1:
                    fail("%s leg %d has %d trim_history rows" % (name, n, len(row)))
                a_leg = float(row[0]["alpha_deg"])
                if abs(a_leg - a_name) > 5e-7 or abs(float(row[0]["Cl"]) - lm["Cl"]) > 5e-8 \
                        or abs(float(row[0]["Cd"]) - lm["Cd"]) > 5e-8:
                    fail("%s leg %d marker and trim_history disagree" % (name, n))
                legs.append(dict(n=n, alpha=a_leg, marker=lm, marker_file=f.name,
                                 latest_time=float(row[0]["latest_time"])))
            if len(legs) != len(th):
                fail("%s: %d leg markers but %d trim_history rows" % (name, len(legs), len(th)))
            out[(cond, L)] = dict(name=name, d=d, marker=mk, ctl=ctl, wm=wm, legs=legs, cells=cells,
                                  mesh=mesh, solver=solver, final_time=float(info["final_time"]))
    return out


def fit_slope(points):
    """Least-squares dCd/dCl in counts per unit CL over (Cl, Cd) points."""
    cl = np.array([p[0] for p in points]); cd = np.array([p[1] for p in points]) * 1e4
    return float(np.polyfit(cl, cd, 1)[0])


def case_record(cond, L, w, old, now):
    c, mk, ctl = COND[cond], w["marker"], w["ctl"]
    src = "%s/solve/%s" % (CLUSTER, w["name"])
    rec = dict(old) if old else {}
    for k in ("superseded_entry", "dragDir", "header_mismatches_across_legs",
              "controlDict_comment_alpha_STALE", "alpha_deg_resid_vs_provenance",
              "provenance_alpha_deg", "provenance_tag_classification",
              "provenance_tag_not_a_case_prefix", "recipe_manifest_sha256",
              "iterations_marker_stated", "endTime_marker_stated"):
        rec.pop(k, None)
    rec.update(
        case=rec.get("case", "%s%s_trim" % (c["prefix"], L)),
        source_case=w["name"], condition=cond, family_letter=L, converged=True,
        status="converged", is_trim_case=True,
        Cd=mk["Cd"], Cd_counts=round(mk["Cd"] * 1e4, 6), Cl=mk["Cl"], L_over_D=mk["Cl"] / mk["Cd"],
        alpha_deg=w["ctl"]["alpha_dragDir"], Aref_m2=ctl["Aref"], lRef_m=ctl["lRef"],
        magUInf_m_s=ctl["magUInf"], dragDir=ctl["dragDir"], endTime=w["final_time"],
        iterations=w["final_time"], leg=mk["leg"], restart_legs=mk["leg"],
        finalised_utc=mk["finalised_utc"], finalised_node=mk["finalised_node"],
        drift_Cd_ct=mk["drift_Cd_ct"], drift_Cl_ct=mk["drift_Cl_ct"], Cl_span_ct=mk["Cl_span_ct"],
        gate_samples=200, gate_verdict=mk["gate_verdict"], marker_format="window_means",
        marker_source="argus_wfix_solve", marker_text=mk["marker_text"],
        window_t_first=w["wm"]["t_first"], window_t_last=w["wm"]["t_last"],
        window_samples_used=w["wm"]["n"],
        Cl_minus_target_counts=(mk["Cl"] - c["target"]) * 1e4,
        mesh_generation="welded",
        mesh_generation_basis=("welded registered surface, 130M-cell-cap recipe, HPC12 "
                               "September 2026 (ARG-196); supersedes the published-mesh trim"),
        mesh_cells=w["cells"], mesh_owner_sha256=w["mesh"]["owner_sha256"],
        geometry_sha256=w["mesh"]["wing_stl_sha256"],
        geometry_source=str(REPO / w["mesh"]["wing_stl_registered_as"]),
        turbulence_model=w["solver"]["turbulence_model"], wall_treatment_nut=w["solver"]["nut_wing"],
        trim_history=[dict(leg=l["n"], alpha_deg=l["alpha"], Cl=l["marker"]["Cl"],
                           Cd=l["marker"]["Cd"], latest_time=l["latest_time"]) for l in w["legs"]],
        sources=dict(
            Cd=src + "/.converged (gate window mean, REPORT THESE)",
            Cl=src + "/.converged (gate window mean, REPORT THESE)",
            Cd_Cl_independent=src + "/postProcessing/forceCoeffs1/<latest leg>/forceCoeffs.dat "
                                    "(last 200 samples, reproduced to 5e-8)",
            alpha_deg=src + "/system/controlDict (dragDir, equal to liftDir and to trim_history)",
            Aref_m2=src + "/system/controlDict", lRef_m=src + "/system/controlDict",
            magUInf_m_s=src + "/system/controlDict", endTime=src + "/system/controlDict",
            marker=src + "/.converged", trim_history=src + "/trim_history.tsv",
            mesh_cells="%s/mesh/%s/constant/polyMesh/owner (header)" % (CLUSTER, w["name"]),
            geometry_sha256="%s/mesh/%s/constant/triSurface/wing.stl" % (CLUSTER, w["name"]),
            turbulence_model=src + "/constant/momentumTransport",
            local_pull=str(w["d"])),
        measured=now)
    if old:
        rec["superseded_entry"] = old
    return rec


def leg_record(cond, L, w, leg, now):
    c, lm = COND[cond], leg["marker"]
    src = "%s/solve/%s" % (CLUSTER, w["name"])
    return dict(
        case="%s_leg%d" % (w["name"], leg["n"]), source_case=w["name"], condition=cond,
        family_letter=L, converged=True, status="converged", is_trim_case=False,
        Cd=lm["Cd"], Cd_counts=round(lm["Cd"] * 1e4, 6), Cl=lm["Cl"], L_over_D=lm["Cl"] / lm["Cd"],
        alpha_deg=leg["alpha"], Aref_m2=w["ctl"]["Aref"], lRef_m=w["ctl"]["lRef"],
        magUInf_m_s=w["ctl"]["magUInf"], iterations=leg["latest_time"], leg=leg["n"],
        finalised_utc=lm["finalised_utc"], finalised_node=lm["finalised_node"],
        drift_Cd_ct=lm["drift_Cd_ct"], drift_Cl_ct=lm["drift_Cl_ct"], Cl_span_ct=lm["Cl_span_ct"],
        gate_samples=200, gate_verdict=lm["gate_verdict"], marker_format="window_means",
        marker_source="argus_wfix_solve", marker_text=lm["marker_text"],
        mesh_generation="welded",
        mesh_generation_basis="pre-nudge converged leg of the welded trim case %s" % w["name"],
        mesh_cells=w["cells"], geometry_sha256=w["mesh"]["wing_stl_sha256"],
        turbulence_model=w["solver"]["turbulence_model"],
        role=("off-trim welded state: converged at the pre-nudge alpha on the same mesh, "
              "condition and recipe as the trim; used for the welded dCd/dCl fit and the "
              "welded C_L-alpha points"),
        sources=dict(marker="%s/%s" % (src, "." + leg["marker_file"]),
                     alpha_deg="%s/trim_history.tsv and the marker file name" % src),
        measured=now)


def patch(path, W, now, bak_stamp):
    d = json.loads(path.read_text())
    cases = d["cases"]
    is_full = "trims" in d
    bak = path.with_name(path.name + ".bak_" + bak_stamp)
    shutil.copy2(path, bak)
    replaced, added = [], []
    for (cond, L), w in W.items():
        name = "%s%s_trim" % (COND[cond]["prefix"], L)
        if name not in cases:
            fail("%s absent from %s; refusing to invent a logical trim" % (name, path.name))
        cases[name] = case_record(cond, L, w, cases[name], now)
        replaced.append(name)
        for leg in w["legs"]:
            r = leg_record(cond, L, w, leg, now)
            cases[r["case"]] = r
            added.append(r["case"])
    n_markers = sum(1 + len(w["legs"]) for w in W.values())
    if len(replaced) + len(added) != n_markers:
        fail("partition: %d replaced + %d added != %d markers read"
             % (len(replaced), len(added), n_markers))

    if is_full:
        fits = d["_dCd_dCl_by_family_condition"]["fits"]
        for (cond, L), w in W.items():
            pts = [(l["marker"]["Cl"], l["marker"]["Cd"]) for l in w["legs"]]
            pts.append((w["marker"]["Cl"], w["marker"]["Cd"]))
            key = "%s/%s" % (cond, L)
            old = fits.get(key)
            span = max(p[0] for p in pts) - min(p[0] for p in pts)
            fits[key] = dict(
                cases_in_fit=["%s_leg%d" % (w["name"], l["n"]) for l in w["legs"]]
                + ["%s%s_trim" % (COND[cond]["prefix"], L)],
                dCd_dCl_counts=fit_slope(pts), n_points=len(pts), mesh_generation="welded",
                Cl_span_of_fit=span,
                physical_estimate=span >= MIN_FIT_SPAN,
                note=("welded points only: pre-nudge legs plus the trim. The published "
                      "alpha-bracket cases are a different mesh generation and are excluded. "
                      + ("" if span >= MIN_FIT_SPAN else
                         "THE C_L SPAN OF THIS FIT IS %.2e, below %.0e: the nudge was tiny, so "
                         "the slope is noise-dominated and is NOT a drag-due-to-lift estimate. "
                         "It is used only to convert the trim residual into a drag bias, which "
                         "it bounds well below 0.1 count either way." % (span, MIN_FIT_SPAN))),
                superseded_fit=old)
        for cond, c in COND.items():
            t = d["trims"][cond]
            base = W[(cond, "B")]["marker"]["Cd"] * 1e4
            for L in LETTERS:
                w, e = W[(cond, L)], t["geometries"][L]
                prev = {k: v for k, v in e.items() if k != "superseded_trim_entry"}
                fit = fits["%s/%s" % (cond, L)]
                cd = w["marker"]["Cd"] * 1e4
                e.clear()
                e.update(prev)
                e.update(
                    Cd_counts=round(cd, 6), Cl=w["marker"]["Cl"], alpha_deg=w["ctl"]["alpha_dragDir"],
                    L_over_D=w["marker"]["Cl"] / w["marker"]["Cd"],
                    Cl_minus_target_counts=(w["marker"]["Cl"] - c["target"]) * 1e4,
                    delta_Cd_counts_vs_baseline=round(cd - base, 6),
                    baseline_Cd_counts_used=round(base, 6),
                    dCd_dCl_counts_per_unit_CL=fit["dCd_dCl_counts"],
                    dCd_dCl_physical_estimate=fit["physical_estimate"],
                    dCd_dCl_fit_cases=fit["cases_in_fit"],
                    Cd_bias_from_trim_residual_counts=fit["dCd_dCl_counts"] * (c["target"] - w["marker"]["Cl"]),
                    trim_case="%s%s_trim" % (c["prefix"], L), source_case=w["name"],
                    trim_tolerance_counts=1.0,
                    trim_within_tolerance=abs(w["marker"]["Cl"] - c["target"]) < TOL_CL,
                    turbulence_model=w["solver"]["turbulence_model"], mesh_generation="welded",
                    mesh_cells=w["cells"],
                    source_Cd="%s/solve/%s/.converged (gate window mean, REPORT THESE)" % (CLUSTER, w["name"]),
                    source_marker="%s/solve/%s/.converged" % (CLUSTER, w["name"]),
                    superseded_trim_entry=prev)
            t["baseline_Cd_counts"] = round(base, 6)
            t["baseline_case"] = "%sB_trim" % c["prefix"]
            t["baseline_source_case"] = W[(cond, "B")]["name"]
            t["mesh_generation"] = "welded"
            welded = [k for k, v in cases.items()
                      if v.get("condition") == cond and v.get("mesh_generation") == "welded"]
            t["coefficient_frame"]["n_converged"] = len(welded)
            t["coefficient_frame"]["n_converged_basis"] = "welded trims plus their pre-nudge legs"
            d["_frame_check_per_condition"][cond]["n_converged"] = len(welded)
            d["_frame_check_per_condition"][cond]["n_converged_basis"] = t["coefficient_frame"]["n_converged_basis"]
        # trim quality, recomputed over every trims block, not patched
        tq, out_tol, worst, n = d["_trim_quality"], [], 0.0, 0
        for cond, t in d["trims"].items():
            tgt = t["target_CL"]
            for L, e in t["geometries"].items():
                n += 1
                r = (e["Cl"] - tgt) * 1e4
                worst = max(worst, abs(r))
                if abs(r) >= 1.0:
                    out_tol.append("%s (%s, %s): Cl %.7f against target %s is %+.3f counts"
                                   % (e.get("trim_case"), cond, L, e["Cl"], tgt, r))
        tq["superseded_OUT_OF_TOLERANCE"] = tq.get("OUT_OF_TOLERANCE")
        tq.update(OUT_OF_TOLERANCE=out_tol, n_trims_adjudicated=n, n_in_tolerance=n - len(out_tol),
                  worst_abs_residual_counts=worst, recomputed=now)
        if tq["n_in_tolerance"] + len(out_tol) != n:
            fail("trim-quality partition broken")
    d["cruise_welded_update_2026_10_05"] = dict(
        what=("Early- and late-cruise trims for all six geometries replaced with the welded-mesh "
              "trims; pre-nudge welded legs added; cruise dCd/dCl fits re-formed on welded points only."),
        decision="ARG-196",
        replaced=sorted(replaced), added=sorted(added),
        untouched=("every Condition CR record, and the published cruise alpha-bracket cases, "
                   "which keep mesh_generation 'published' and are in no cruise fit"),
        gates=["marker reproduced from forceCoeffs to 5e-8", "alpha from dragDir = liftDir = "
               "trim_history to 1e-6 deg", "Aref/lRef DSO half-model, magUInf per condition",
               "|C_L - target| < 1e-4 on all twelve", "replaced + added = markers read"],
        backup=bak.name, generator="docs/report/all_geometry_2026-09-15/scripts/update_cruise_from_welded.py")
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(d, indent=1))
    json.loads(tmp.read_text())
    os.replace(tmp, path)
    return d, replaced, added


def main():
    W = gather()
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    stamp = now.replace("-", "").replace(":", "")
    for f in DATA_FILES:
        if not f.exists():
            fail("missing %s" % f)
    for f in DATA_FILES:
        d, rep, add = patch(f, W, now, stamp)
        print("== %s: %d trims replaced, %d leg records added" % (f.name, len(rep), len(add)))
        if "trims" in d:
            for cond in COND:
                t = d["trims"][cond]
                print("   %s  baseline %s = %.3f ct" % (cond, t["baseline_source_case"], t["baseline_Cd_counts"]))
                for L in LETTERS:
                    e = t["geometries"][L]
                    print("     %s  alpha %.6f  Cl %.7f (%+.3f ct)  Cd %8.3f  dCd %+7.3f  "
                          "slope %6.1f  bias %+.4f  cells %d"
                          % (L, e["alpha_deg"], e["Cl"], e["Cl_minus_target_counts"], e["Cd_counts"],
                             e["delta_Cd_counts_vs_baseline"], e["dCd_dCl_counts_per_unit_CL"],
                             e["Cd_bias_from_trim_residual_counts"], e["mesh_cells"]))
            print("   trim quality: %d of %d in tolerance, worst %.3f ct"
                  % (d["_trim_quality"]["n_in_tolerance"], d["_trim_quality"]["n_trims_adjudicated"],
                     d["_trim_quality"]["worst_abs_residual_counts"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
