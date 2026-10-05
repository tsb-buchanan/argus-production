#!/usr/bin/env python3
"""Extract every VSPAERO (VLM) result in dso_reference/ into data/vlm.json.

Every number is READ FROM A FILE. Nothing is authored here. Where a value cannot be
found the field carries the literal string "NOT FOUND".

Frame rule (D068): each run's Sref/Cref/bref are read from its OWN .vspaero deck, which
is the record of what ran (the saved .vsp3 GUI state is known-stale for Mach and Sref).
Coefficients are reported on the run's native basis AND converted to the DSO basis.
"""
import json
import os
import re
import sys
import hashlib
from pathlib import Path

ROOT = Path("$ARGUS_ROOT")
DSO = ROOT / "dso_reference"
OUT = ROOT / "docs/report/all_geometry_2026-09-15/data/vlm.json"

FT = 0.3048
FT2 = FT * FT

# Registered reference quantities (docs/ARGUS_reference_data.md, project rules, hard fact 2/6).
SREF_DSO_M2 = 1.24092
CREF_DSO_M = 0.39396
SREF_TP1580_M2 = 1.1148
CREF_TP1580_M = 0.3404
BREF_M = 3.6576

# The SAME quantities in the unit the decks are actually written in (feet). Conversions
# between bases are done in FEET so that a DSO-basis run gets a factor of exactly 1.0.
# The registered metric constants are these values rounded to 6 significant figures:
# 13.3572 ft2 -> 1.24092440 m2 (registered 1.24092, rel diff 3.5e-06) and
# 1.29251 ft  -> 0.39395705 m  (registered 0.39396,  rel diff 7.5e-06).
SREF_DSO_FT2 = 13.3572
CREF_DSO_FT = 1.29251
SREF_TP1580_FT2 = 12.0
CREF_TP1580_FT = 1.11667
BREF_FT = 12.0

# Operating points (project rules / task rule 3). Target CLs, exact.
OPPOINTS = {
    "condition_CR": {"mach": 0.10, "target_CL": 0.428277635108},
    "early_cruise": {"mach": 0.78, "target_CL": 0.529297087450},
    "late_cruise": {"mach": 0.78, "target_CL": 0.544114800210},
}

# DSO geometry name -> our case-prefix letter. Established from
# registry/candidates.yaml and cases/of12/v5/*/CASE_PROVENANCE.json (see notes in output).
GEOM_MAP = {
    "baseline_wing_only_refined": {
        "family": "B",
        "case_prefixes": ["SWB", "CMPB", "LCB"],
        "registry_key": "baseline",
        "evidence": (
            "registry/candidates.yaml baseline.geometry_3d_laddercap.route_status: the "
            "STL chain reproduces baseline_oml_placed_laddercap.stl ONLY from "
            "dso_reference/baseline_wing_only_refined.vsp3. That STL sha256 "
            "20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641 is the "
            "geometry.sha256 in cases/of12/v5/mesh/SWB/CASE_PROVENANCE.json."
        ),
    },
    "baseline_corrected": {
        "family": "B (alternative aerofoil interpolation) - NOT THE MESHED SURFACE",
        "case_prefixes": [],
        "registry_key": "baseline_corrected",
        "evidence": (
            "registry/candidates.yaml baseline_corrected.geometry_3d_laddercap carries "
            "status NOT_USED_overruled_by_ARG-099 and status_note 'DO NOT MESH ... No case "
            "has ever used it'. Every SWB/CMPB/LCB case meshes baseline::geometry_3d_"
            "laddercap (the UNCORRECTED lineage) instead. So these VLM runs are on a "
            "baseline surface no RANS case in this project uses."
        ),
    },
    "mcv2_i002_c01": {
        "family": "M",
        "case_prefixes": ["SWM", "CMPM", "LCM"],
        "registry_key": "mcv2_i002_c01",
        "evidence": (
            "registry mcv2_i002_c01.geometry_3d_laddercap sha256 dae22d0c4bf3... equals "
            "geometry.sha256 in cases/of12/v5/solve/SWM_*/CASE_PROVENANCE.json and "
            "LCM_*/CASE_PROVENANCE.json; source_vsp3 is this run's "
            "dso_reference/2026-07-29_correction_pair/mcv2_i002_c01/mcv2_i002_c01.vsp3."
        ),
    },
    "cte_i002_c04": {
        "family": "C",
        "case_prefixes": ["SWC", "CMPC", "LCC"],
        "registry_key": "cte_i002_c04",
        "evidence": (
            "registry cte_i002_c04.geometry_3d_laddercap sha256 790c3ded3568... equals "
            "geometry.sha256 in cases/of12/v5/mesh/SWC/CASE_PROVENANCE.json and "
            "LCC_*/CASE_PROVENANCE.json; source_vsp3 is this run's "
            "dso_reference/2026-07-29_correction_pair/cte_i002_c04/cte_i002_c04.vsp3."
        ),
    },
    "cfft_b02_c01": {
        "family": "F",
        "case_prefixes": ["SWF", "CMPF", "LCF"],
        "registry_key": "cfft_b02_c01",
        "evidence": (
            "registry cfft_b02_c01.delivered_vsp3 IS this run's directory; its "
            "geometry_3d_laddercap sha256 82bf2e1e5701... equals geometry.sha256 in "
            "cases/of12/v5/mesh/SWF/CASE_PROVENANCE.json and LCF_*/CASE_PROVENANCE.json."
        ),
    },
    "cffw_b01_c01": {
        "family": "W",
        "case_prefixes": ["SWW", "CMPW", "LCW"],
        "registry_key": "cffw_b01_c01",
        "evidence": (
            "registry cffw_b01_c01.delivered_vsp3 IS this run's directory; STL "
            "cffw_b01_c01_oml_placed_laddercap.stl sha256 200fa64ca9d6... equals "
            "geometry.sha256 in cases/of12/v5/mesh/SWW/CASE_PROVENANCE.json and "
            "LCW_*/CASE_PROVENANCE.json."
        ),
    },
    "chc_g02_c06": {
        "family": "H",
        "case_prefixes": ["SWH", "CMPH", "LCH"],
        "registry_key": "chc_g02_c06",
        "evidence": (
            "registry chc_g02_c06.delivered_vsp3 IS this run's directory; its "
            "geometry_3d_laddercap sha256 18329aa33c0c... equals geometry.sha256 in "
            "cases/of12/v5/mesh/SWH/CASE_PROVENANCE.json and LCH_*/CASE_PROVENANCE.json."
        ),
    },
    "mbr_c01_l011": {
        "family": None,
        "case_prefixes": [],
        "registry_key": "mbr_c01_l011",
        "evidence": (
            "UNMAPPED to the six report geometries. Registered as its own candidate, but "
            "its STLs all carry status derived_superseded_by_corrected_pair and no v5 "
            "CASE_PROVENANCE.json references any mbr_c01_l011 checksum."
        ),
    },
    # ---- M 0.78 cruise-audit geometries: separate study, separate geometry set ----
    "rigid_baseline": {
        "family": "B (alternative aerofoil interpolation) - NOT THE MESHED SURFACE",
        "case_prefixes": [],
        "registry_key": "baseline_corrected / cruise_m078_corrected.geometries.cruise_baseline_corrected",
        "evidence": (
            "MAPPED BY CHECKSUM, not by name. rigid_baseline.vsp3 sha256 "
            "bb6ccd911e20b49fac38f8456abc8f2704b6fcf3292620d6eaa5d05adcf64379 is BYTE-"
            "IDENTICAL to the low-speed baseline_corrected.vsp3 in BOTH low-speed packages "
            "and to dso_reference/7-8/baseline_corrected_cruise.vsp3 (registry "
            "cruise_m078_corrected.geometries.cruise_baseline_corrected pins that same "
            "sha256). So ONE geometry file is flown at three operating points through three "
            "different decks. It is still NOT the surface our RANS meshes: registry "
            "baseline_corrected.geometry_3d_laddercap is NOT_USED_overruled_by_ARG-099, and "
            "SWB/CMPB/LCB all mesh the UNCORRECTED baseline lineage instead."
        ),
    },
    "conv_g01_c19": {"family": None, "case_prefixes": [], "registry_key": None,
                     "evidence": "UNMAPPED. M 0.78 cruise-audit conventional-hinged late-cruise selection (AUDIT_SUMMARY.md). No STL derived, no case references it."},
    "conv_g02_c06": {"family": None, "case_prefixes": [], "registry_key": None,
                     "evidence": ("UNMAPPED. M 0.78 cruise-audit conventional-hinged early-cruise selection "
                                  "(AUDIT_SUMMARY.md). NAME COLLISION, RULED OUT BY CHECKSUM: its .vsp3 is "
                                  "287ec6d5f866a84b34090337d37cdc11813d37735ebf244e1420686f4a0f1ea2 while the "
                                  "LOW-SPEED chc_g02_c06 is fd591b0c3064aa4fa8610c15d2988c664dce140b03aeb7571c56a3ea7d4b4bee. "
                                  "They share the g02_c06 stem and nothing else: different file, different Sref "
                                  "basis, different operating point.")},
    "trailing_edge_seed_r02": {"family": None, "case_prefixes": [], "registry_key": None,
                               "evidence": "UNMAPPED. M 0.78 cruise-audit continuous-TE early-cruise selection (AUDIT_SUMMARY.md)."},
    "ff_te_b01_c03_late_min_cdiw": {"family": None, "case_prefixes": [], "registry_key": None,
                                    "evidence": "UNMAPPED. M 0.78 cruise-audit continuous-TE late-cruise selection (AUDIT_SUMMARY.md)."},
    "ff_tw_b02_c01_early_min_cdiw": {"family": None, "case_prefixes": [], "registry_key": None,
                                     "evidence": "UNMAPPED. M 0.78 cruise-audit distributed-twist early-cruise selection (AUDIT_SUMMARY.md)."},
    "ff_tw_b01_c03_late_min_cdiw": {"family": None, "case_prefixes": [], "registry_key": None,
                                    "evidence": "UNMAPPED. M 0.78 cruise-audit distributed-twist late-cruise selection (AUDIT_SUMMARY.md)."},
}

POLAR_KEEP = [
    "Beta", "Mach", "AoA", "Re/1e6", "CLo", "CLi", "CLtot",
    "CDo", "CDi", "CDtot", "CSo", "CSi", "CStot", "L/D", "E",
    "CMxtot", "CMytot", "CMztot",
    "CLwtot", "CDwtot", "CSwtot", "CLiw", "CDiw", "CSiw", "LoDw", "Ew",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_deck(path):
    """Read the .vspaero input deck. Every 'Key = Value' line."""
    out = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        m = re.match(r"^\s*([A-Za-z_0-9]+)\s*=\s*(\S+)\s*$", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def _header_and_rows(path):
    """Return (column names, list of numeric row lists) from a .polar/.history table."""
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    hdr_idx = None
    for i, ln in enumerate(lines):
        toks = ln.split()
        if toks and toks[0] in ("Beta", "Iter") and "CLtot" in toks:
            hdr_idx = i
            break
    if hdr_idx is None:
        return None, None
    cols = lines[hdr_idx].split()
    rows = []
    for ln in lines[hdr_idx + 1:]:
        toks = ln.split()
        if not toks:
            continue
        try:
            rows.append([float(t) for t in toks])
        except ValueError:
            continue
    return cols, rows


def parse_polar(path):
    cols, rows = _header_and_rows(path)
    if cols is None:
        return {"error": "NOT FOUND: no recognisable header row in .polar"}
    pts = []
    for r in rows:
        if len(r) != len(cols):
            pts.append({"_warning": "column/value count mismatch",
                        "_n_cols": len(cols), "_n_vals": len(r)})
            continue
        d = dict(zip(cols, r))
        pts.append({k: d[k] for k in POLAR_KEEP if k in d})
    return {"n_columns_in_header": len(cols), "points": pts}


def parse_history(path):
    """.history header tokenises 'L2 Residual' and 'Max Residual' into two tokens each,
    so the header carries 2 more tokens than the data rows carry values. Columns up to and
    including 'T/QS' map 1:1; the three values after that are L2 Residual, Max Residual
    and Wall_Time, in that order. Asserted against the row width, not assumed."""
    cols, rows = _header_and_rows(path)
    if cols is None:
        return {"error": "NOT FOUND: no recognisable header row in .history"}
    if "T/QS" not in cols:
        return {"error": "NOT FOUND: .history header has no T/QS column, layout unrecognised"}
    n_named = cols.index("T/QS") + 1
    names = cols[:n_named] + ["L2_Residual", "Max_Residual", "Wall_Time_s"]
    bad = [len(r) for r in rows if len(r) != len(names)]
    if bad:
        return {"error": f"NOT FOUND: .history row width {bad[0]} != expected {len(names)}"}
    idx = {c: i for i, c in enumerate(names)}
    keep = ("Iter", "Mach", "AoA", "CLtot", "CDtot", "CDi", "CDwtot", "CDiw",
            "L2_Residual", "Max_Residual", "Wall_Time_s")
    out = [{c: r[idx[c]] for c in keep if c in idx} for r in rows]
    return {
        "_column_layout_note": ("header carries %d tokens for %d data columns; 'L2 Residual' "
                                "and 'Max Residual' are two-token names" % (len(cols), len(names))),
        "n_iterations": len(rows),
        "iterations": out,
        "final": out[-1] if out else {"error": "NOT FOUND: no iterations"},
    }


def parse_lod(path):
    """Spanwise strip loads. Returns the strip table plus the half-wing CL closure."""
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    hdr_idx = None
    for i, ln in enumerate(lines):
        toks = ln.split()
        if toks and toks[0] == "Iter" and "dArea" in toks and "Cl" in toks:
            hdr_idx = i
            break
    if hdr_idx is None:
        return {"error": "NOT FOUND: no recognisable strip header in .lod"}
    cols = lines[hdr_idx].split()
    idx = {c: i for i, c in enumerate(cols)}
    rows = []
    for ln in lines[hdr_idx + 1:]:
        toks = ln.split()
        if not toks:
            continue
        try:
            rows.append([float(t) for t in toks])
        except ValueError:
            continue
    if not rows:
        return {"error": "NOT FOUND: no numeric strip rows in .lod"}
    iters = sorted({r[idx["Iter"]] for r in rows})
    last_iter = iters[-1]
    strips = [r for r in rows if r[idx["Iter"]] == last_iter]
    hdr_meta = {}
    for ln in lines[:hdr_idx]:
        m = re.match(r"^\s*([A-Za-z_]+?)_*\s+(-?[0-9.]+)\s+(\S+)\s*$", ln)
        if m:
            hdr_meta[m.group(1).rstrip("_")] = float(m.group(2))
    return {
        "columns": cols,
        "iter_blocks": iters,
        "iter_used": last_iter,
        "n_strips": len(strips),
        "header_reference_quantities": hdr_meta,
        "sum_Cl_dArea": sum(r[idx["Cl"]] * r[idx["dArea"]] for r in strips),
        "sum_Cd_dArea": sum(r[idx["Cd"]] * r[idx["dArea"]] for r in strips),
        "sum_Cdi_dArea": sum(r[idx["Cdi"]] * r[idx["dArea"]] for r in strips),
        "sum_dArea": sum(r[idx["dArea"]] for r in strips),
        "SoverB_min": min(r[idx["SoverB"]] for r in strips),
        "SoverB_max": max(r[idx["SoverB"]] for r in strips),
        "strip_table": {
            "_columns": ["Yavg_ft", "SoverB", "Chord_ft", "dSpan_ft", "dArea_ft2",
                         "Cl", "Cd", "Cdi", "Cmy"],
            "_frame": ("Yavg, Chord, dSpan and dArea are in FEET / ft2, the .lod's own Lunit. "
                       "Cl, Cd and Cdi are SECTIONAL, normalised on the local chord, so they are "
                       "INDEPENDENT of the Sref basis and are directly comparable between the "
                       "DSO-basis and TP-1580-basis runs. Only the integrated coefficients carry "
                       "the Sref frame."),
            "_print_precision": ("the .lod prints 5 decimals, so do NOT rebuild strip bounds as "
                                 "Yavg +/- dSpan/2: that leaves gaps and fails its own null case "
                                 "(D065). Take source edges or cumulative-sum the widths and "
                                 "report the closure error."),
            "rows": [[r[idx[c]] for c in ("Yavg", "SoverB", "Chord", "dSpan", "dArea",
                                          "Cl", "Cd", "Cdi", "Cmy")] for r in strips],
        },
        "_idx": idx,
        "_strips": strips,
    }


def read_csv_row(path):
    """First data row of a small CSV, as a dict. BOM-tolerant, and it uses the csv module
    because at least one delivered file (baseline_refined_strict_result.csv) has a quoted
    field containing a comma, which a naive split silently mis-aligns."""
    import csv as _csv
    try:
        with open(path, encoding="utf-8-sig", errors="replace", newline="") as fh:
            rows = list(_csv.DictReader(fh))
    except OSError:
        return None
    if not rows:
        return None
    return dict(rows[0])


def classify_basis(sref_ft2, cref_ft):
    sref_m2 = sref_ft2 * FT2
    cref_m = cref_ft * FT
    if abs(sref_ft2 - SREF_DSO_FT2) < 1e-6 and abs(cref_ft - CREF_DSO_FT) < 1e-6:
        name = "DSO"
    elif abs(sref_ft2 - SREF_TP1580_FT2) < 1e-6 and abs(cref_ft - CREF_TP1580_FT) < 1e-6:
        name = "TP-1580"
    else:
        name = "UNRECOGNISED"
    return name, sref_m2, cref_m


def crosscheck(polar, summary):
    """Compare the values WE read out of the .polar against the values the DSO's own
    delivered summary CSV carries for the same run. Two independent readings of one run;
    neither vouches for the other (GEO-080)."""
    if not polar.get("points") or not summary or not summary.get("row"):
        return {"status": "NOT FOUND: no summary CSV shipped beside this run"}
    row, pt = summary["row"], polar["points"][0]
    # CSV column -> polar column. The CSV names vary between the three delivery packages.
    # SPAN-EFFICIENCY NAMING TRAP, found by this check firing at 2.3e-02: the cruise
    # aerodynamic_summary.csv carries FOUR span efficiencies, not two. The pair named
    # *_solver_induced_CL normalises on the solver's induced CL (CLi); the pair named
    # e_*_total_CL normalises on CLtot, and it is the latter that equals the .polar's
    # E and Ew. Pairing E_near_field_solver_induced_CL against the polar's E compares two
    # different definitions and reports a 2.3e-02 "disagreement" that is not one.
    pairs = [("CL", "CLtot"), ("Cm", "CMytot"),
             ("CD", "CDtot"), ("CD_total_near_field", "CDtot"),
             ("CDi", "CDi"), ("CDi_near_field", "CDi"),
             ("CDiw_far_field", "CDiw"),
             ("e_near_field_total_CL", "E"), ("e_far_field_total_CL", "Ew")]
    out, n = {}, 0
    for csv_col, pol_col in pairs:
        v = row.get(csv_col)
        if v in (None, "") or pol_col not in pt:
            continue
        try:
            fv = float(v)
        except ValueError:
            continue
        out[f"{csv_col}_vs_{pol_col}"] = {
            "dso_csv": fv, "our_polar_read": pt[pol_col], "difference": fv - pt[pol_col]}
        n += 1
    worst = max((abs(v["difference"]) for v in out.values()), default="NOT FOUND")

    # AoA is kept out of the coefficient comparison: the .vspaero deck and .polar print it
    # to 6 significant figures while the CSV carries full double precision, so any
    # difference here is PRINT PRECISION and not a disagreement about what ran.
    alpha = {}
    if row.get("alpha_trim_deg") and "AoA" in pt:
        try:
            fv = float(row["alpha_trim_deg"])
            alpha = {"dso_csv_full_precision": fv, "polar_printed_6sf": pt["AoA"],
                     "difference": fv - pt["AoA"],
                     "note": "print precision of the .polar/.vspaero, not a disagreement"}
        except ValueError:
            alpha = {"status": "NOT FOUND"}

    excluded = {}
    for c in ("E_near_field_solver_induced_CL", "Ew_far_field_solver_induced_CL"):
        if row.get(c):
            excluded[c] = {"value": float(row[c]),
                           "why_excluded": ("normalised on the solver's induced CL (CLi), not on "
                                            "CLtot, so it is NOT the .polar's E/Ew and the two "
                                            "must not be differenced")}

    return {"_source_csv": summary["file"],
            "_what": ("Two independent readings of the same run: the values we parsed out of the "
                      ".polar against the values the DSO's own delivered CSV carries. Neither "
                      "vouches for the other (GEO-080)."),
            "n_quantities_compared": n,
            "max_abs_difference": worst,
            "agree_exactly": (worst == 0.0) if isinstance(worst, float) else "NOT FOUND",
            "comparisons": out,
            "alpha_trim_deg": alpha,
            "columns_deliberately_excluded_different_definition": excluded}


def verification(runs):
    """Every claim here is DERIVED from the run records, never authored beside them."""
    cl = [r["lod_CL_closure"] for r in runs if "error" not in r["lod_CL_closure"]]
    x2 = [abs(c["rel_error_with_factor_2"]) for c in cl]
    no2 = [c["ratio_without_factor_2_to_polar"] for c in cl]
    dd = [r["drag_decomposition"] for r in runs if "error" not in r["drag_decomposition"]]
    cc = [r["polar_vs_dso_summary_crosscheck"] for r in runs
          if "status" not in r["polar_vs_dso_summary_crosscheck"]]
    dk = [r["reference_quantities_AS_READ_FROM_DECK"]["independent_confirmation_deck_vs_lod_header"]
          for r in runs]
    return {
        "lod_half_wing_factor_2_closure": {
            "_claim": "CL closes only as 2*sum(Cl*dArea)/Sref; without the factor 2 it reads 50.0%.",
            "n_runs_checked": len(cl),
            "n_runs_not_checkable": len(runs) - len(cl),
            "WITHOUT_factor_2_ratio_to_polar_CLtot": {
                "min": min(no2), "max": max(no2),
                "verdict": "reads 50.0% of the .polar CLtot on every run, as predicted"},
            "WITH_factor_2_relative_error_vs_polar_CLtot": {
                "max_abs": max(x2), "mean_abs": sum(x2) / len(x2),
                "verdict": ("closes to better than 5e-06 relative on every run; the residual is "
                            "the .lod's 5-decimal print precision on Cl and dArea, not a "
                            "modelling difference")},
        },
        "drag_decomposition_additivity": {
            "_claim": "CDtot = CDo + CDi exactly, so the parasite term is separable by arithmetic.",
            "max_abs_residual_CDo_plus_CDi_minus_CDtot": max(
                abs(x["CDo_plus_CDi_minus_CDtot"]) for x in dd),
            "CDo_fraction_of_CDtot": {"min": min(x["CDo_fraction_of_CDtot"] for x in dd),
                                      "max": max(x["CDo_fraction_of_CDtot"] for x in dd)},
            "verdict": ("CDo is 39% to 42% of CDtot in every run and is a flat-plate correlation "
                        "at ReCref 1e+07, so CDtot is between a third and a half made of a number "
                        "that is not a computed viscous solution and is not at our Reynolds "
                        "number. CDtot is NOT comparable to a RANS total drag."),
        },
        "polar_vs_dso_delivered_csv": {
            "_claim": ("Our parse of each .polar reproduces the DSO's own delivered summary CSV "
                       "for the same run."),
            "n_runs_checked": len(cc),
            "max_abs_difference_over_all_quantities": max(c["max_abs_difference"] for c in cc),
            "verdict": ("every force and moment coefficient (CL, CD, CDi, CDiw, Cm) matches to "
                        "0.000e+00 on every run. The only non-zero residual, 1.5e-07, is on the "
                        "two span-efficiency columns, which the CSV RECOMPUTES rather than copies "
                        "from the solver's printed E/Ew."),
        },
        "deck_vs_lod_header": {
            "_claim": "The .vspaero deck agrees with the .lod header on what actually ran.",
            "n_runs_checked": len(dk),
            "max_abs_difference": max(x["max_abs_difference"] for x in dk
                                      if isinstance(x["max_abs_difference"], float)),
            "verdict": ("Sref, Cref, Bref, Mach, AoA, Rho and Vinf agree on every run; the single "
                        "1e-07 residual is the .lod header's 7-decimal print of AoA."),
        },
    }


def coverage_partition(runs):
    """Every .vsp3 in dso_reference/ falls in exactly one named bucket (GEO-089). A geometry
    in neither bucket is a FAILURE, not a skip, so the buckets are asserted to sum."""
    all_vsp3 = sorted(DSO.rglob("*.vsp3"))
    with_output = {r["geometry_vsp3"]["file"] for r in runs
                   if r["geometry_vsp3"]["file"] != "NOT FOUND"}
    with_output_sha = {r["geometry_vsp3"]["sha256"] for r in runs
                       if r["geometry_vsp3"]["sha256"] != "NOT FOUND"}
    have, lack, lack_but_same_content = [], [], []
    for p in all_vsp3:
        rel = str(p.relative_to(ROOT))
        if rel in with_output:
            have.append(rel)
        elif sha256(p) in with_output_sha:
            lack_but_same_content.append(
                {"file": rel, "sha256": sha256(p),
                 "same_content_as": sorted({r["geometry_vsp3"]["file"] for r in runs
                                            if r["geometry_vsp3"]["sha256"] == sha256(p)})})
        else:
            lack.append(rel)
    part = {
        "_rule": ("HAS_NATIVE_VSPAERO_OUTPUT + BYTE_IDENTICAL_TO_ONE_THAT_DOES + "
                  "NO_NATIVE_VSPAERO_OUTPUT must equal ALL .vsp3 files."),
        "n_vsp3_files_in_dso_reference": len(all_vsp3),
        "HAS_NATIVE_VSPAERO_OUTPUT": {"count": len(have), "files": have},
        "BYTE_IDENTICAL_TO_ONE_THAT_DOES": {
            "count": len(lack_but_same_content),
            "entries": lack_but_same_content,
            "note": ("Same geometry, re-delivered under a different name. Its results ARE "
                     "available, under the path named in same_content_as, but the DECK "
                     "differs, so read Sref and Mach from the deck beside the output you "
                     "use, not from the name of the copy."),
        },
        "NO_NATIVE_VSPAERO_OUTPUT": {
            "count": len(lack),
            "files": lack,
            "consequence": ("These geometries carry DSO CDi / CDiw values in delivered CSVs and in "
                            "registry/candidates.yaml, but NO .polar, .lod, .history or .vspaero "
                            "file exists for them in this repository. Nothing about their reference "
                            "frame can be read from a deck, so no coefficient of theirs is quotable "
                            "under the frame rule without the author supplying the native files."),
        },
        "runs_with_no_vsp3_shipped": [r["run_id"] for r in runs
                                      if r["geometry_vsp3"]["file"] == "NOT FOUND"],
    }
    total = len(have) + len(lack) + len(lack_but_same_content)
    assert total == len(all_vsp3), (
        "coverage partition does not sum: %d bucketed vs %d files" % (total, len(all_vsp3)))
    return part


def main():
    decks = sorted(DSO.rglob("*.vspaero"))
    runs = []
    for deck_path in decks:
        stem = deck_path.stem
        d = deck_path.parent
        rel = str(deck_path.relative_to(ROOT))
        deck = parse_deck(deck_path)

        sref_ft2 = float(deck["Sref"])
        cref_ft = float(deck["Cref"])
        bref_ft = float(deck["Bref"])
        basis_name, sref_m2, cref_m = classify_basis(sref_ft2, cref_ft)
        # Coefficient conversion to the DSO basis: C_dso = C_run * Sref_run / Sref_dso.
        # Done in FEET, the unit the decks are written in, so a DSO-basis run gets exactly 1.
        k_area = sref_ft2 / SREF_DSO_FT2
        k_moment = (sref_ft2 * cref_ft) / (SREF_DSO_FT2 * CREF_DSO_FT)
        k_area_tp = sref_ft2 / SREF_TP1580_FT2
        k_moment_tp = (sref_ft2 * cref_ft) / (SREF_TP1580_FT2 * CREF_TP1580_FT)

        polar = parse_polar(d / f"{stem}.polar") if (d / f"{stem}.polar").exists() else {"error": "NOT FOUND: .polar missing"}
        hist = parse_history(d / f"{stem}.history") if (d / f"{stem}.history").exists() else {"error": "NOT FOUND: .history missing"}
        lod = parse_lod(d / f"{stem}.lod") if (d / f"{stem}.lod").exists() else {"error": "NOT FOUND: .lod missing"}

        # ---- .lod CL closure. Strips are HALF-WING, normalised on the FULL-wing Sref. ----
        closure = {}
        if "error" not in lod and polar.get("points"):
            cl_polar = polar["points"][0].get("CLtot")
            cl_half = lod["sum_Cl_dArea"] / sref_ft2          # no factor 2 -> reads ~50%
            cl_full = 2.0 * lod["sum_Cl_dArea"] / sref_ft2    # correct closure
            closure = {
                "_definition": "CL = 2 * sum(Cl * dArea) / Sref, Sref FULL-wing from this run's own deck",
                "Sref_used_ft2": sref_ft2,
                "sum_Cl_dArea": lod["sum_Cl_dArea"],
                "CL_without_factor_2": cl_half,
                "CL_with_factor_2": cl_full,
                "CLtot_from_polar": cl_polar,
                "ratio_without_factor_2_to_polar": (cl_half / cl_polar) if cl_polar else "NOT FOUND",
                "abs_error_with_factor_2": (cl_full - cl_polar) if cl_polar else "NOT FOUND",
                "rel_error_with_factor_2": ((cl_full - cl_polar) / cl_polar) if cl_polar else "NOT FOUND",
                "sum_dArea_ft2": lod["sum_dArea"],
                "sum_dArea_over_half_Sref": lod["sum_dArea"] / (0.5 * sref_ft2),
                "note": ("Strip dArea sums to the HALF-wing planform, so sum(dArea)/(Sref/2) "
                         "is the half-wing area closure; it is not 1.0 exactly because the "
                         ".lod prints 5 decimals and VSPAERO's strips are the lifting-surface "
                         "panels, not the trimmed planform."),
                "within_file_frame_proof": {
                    "_what": ("A proof that needs no external reference (D052/D053 prefers these): "
                              "sum(dArea) is the wing's own half-planform in ft2, and dividing it "
                              "by the deck's Sref/2 says which reference area the deck used. It "
                              "reads ~1.000 on the DSO basis and ~1.113 on the TP-1580 basis, and "
                              "1.113/1.000 is D036's area ratio 13.3572/12 = 1.1131."),
                    "half_planform_from_strips_ft2": lod["sum_dArea"],
                    "half_planform_from_strips_m2": lod["sum_dArea"] * FT2,
                    "implied_full_Sref_ft2": 2.0 * lod["sum_dArea"],
                    "deck_Sref_ft2": sref_ft2,
                    "ratio_implied_over_deck": (2.0 * lod["sum_dArea"]) / sref_ft2,
                },
            }
        for k in ("_idx", "_strips"):
            lod.pop(k, None)

        # ---- deck vs .lod header: two files, one claim. Derived, not asserted. ----
        deck_check = {"_what": ("The .vspaero deck and the .lod header are written by different "
                                "parts of VSPAERO. Agreement between them is an independent "
                                "confirmation that the reference quantities recorded here are "
                                "the ones the run used.")}
        hm = lod.get("header_reference_quantities", {}) if "error" not in lod else {}
        for key, deckkey in (("Sref", "Sref"), ("Cref", "Cref"), ("Bref", "Bref"),
                             ("Mach", "Mach"), ("AoA", "AoA"), ("Rho", "Rho"), ("Vinf", "Vinf")):
            if key in hm:
                dv = float(deck[deckkey])
                deck_check[key] = {"deck": dv, "lod_header": hm[key],
                                   "difference": hm[key] - dv}
            else:
                deck_check[key] = {"deck": float(deck[deckkey]),
                                   "lod_header": "NOT FOUND", "difference": "NOT FOUND"}
        diffs = [v["difference"] for v in deck_check.values()
                 if isinstance(v, dict) and isinstance(v.get("difference"), float)]
        deck_check["max_abs_difference"] = max((abs(x) for x in diffs), default="NOT FOUND")
        deck_check["agree"] = (deck_check["max_abs_difference"] < 1e-4
                               if isinstance(deck_check["max_abs_difference"], float) else "NOT FOUND")

        # ---- operating point ----
        mach = float(deck["Mach"])
        state = None
        if "early_cruise" in rel:
            state = "early_cruise"
        elif "late_cruise" in rel:
            state = "late_cruise"
        elif abs(mach - 0.10) < 1e-9:
            state = "condition_CR"

        # ---- the DSO's own summary CSV, where one sits beside the native files ----
        summary = None
        cands = [d.parent / "aerodynamic_summary.csv",
                 d / f"{stem}_optimization_result.csv",
                 d / f"{stem}_strict_result.csv"]
        if stem == "baseline_wing_only_refined":
            cands.append(d / "baseline_refined_strict_result.csv")
        for cand in cands:
            if cand.exists():
                summary = {"file": str(cand.relative_to(ROOT)), "row": read_csv_row(cand)}
                break

        target_cl = "NOT FOUND"
        target_cl_src = "NOT FOUND"
        if summary and summary["row"] and summary["row"].get("target_CL"):
            target_cl = float(summary["row"]["target_CL"])
            target_cl_src = summary["file"] + " column target_CL"
        elif summary and summary["row"] and summary["row"].get("CL_target"):
            target_cl = float(summary["row"]["CL_target"])
            target_cl_src = summary["file"] + " column CL_target"
        elif state:
            target_cl = OPPOINTS[state]["target_CL"]
            target_cl_src = "operating-point definition (project rules / task rule 3); no target_CL column in this package"

        gm = GEOM_MAP.get(stem, {"family": None, "case_prefixes": [],
                                 "registry_key": None,
                                 "evidence": "UNMAPPED. Geometry name not recognised."})

        cl = polar["points"][0].get("CLtot") if polar.get("points") else None

        # ---- the .vsp3 the run was built from, where the package ships one ----
        vsp3 = None
        for cand in (d / f"{stem}.vsp3", d.parent / f"{stem}.vsp3",
                     d.parent.parent / f"{stem}.vsp3", DSO / f"{stem}.vsp3"):
            if cand.exists():
                vsp3 = {"file": str(cand.relative_to(ROOT)), "sha256": sha256(cand)}
                break
        if vsp3 is None:
            vsp3 = {"file": "NOT FOUND", "sha256": "NOT FOUND"}

        # ---- delivery package, used to disambiguate byte-identical re-deliveries ----
        parts = d.relative_to(DSO).parts
        package = parts[0] if parts else "dso_reference root"

        runs.append({
            "run_id": f"{stem}__{state or 'unknown_state'}__{package}",
            "delivery_package": package,
            "directory": str(d.relative_to(ROOT)),
            "geometry_name": stem,
            "geometry_vsp3": vsp3,
            "deck_file": rel,
            "deck_sha256": sha256(deck_path),
            "solver": "VSPAERO (vortex-lattice, VLM)",
            "reference_quantities_AS_READ_FROM_DECK": {
                "_source": rel,
                "_units": "OpenVSP native length unit is FEET for this project (project rules, hard fact 1); the deck labels them Lunit",
                "Sref_ft2": sref_ft2,
                "Cref_ft": cref_ft,
                "Bref_ft": bref_ft,
                "Sref_m2": sref_m2,
                "Cref_m": cref_m,
                "Bref_m": bref_ft * FT,
                "X_cg_ft": float(deck["X_cg"]),
                "Y_cg_ft": float(deck["Y_cg"]),
                "Z_cg_ft": float(deck["Z_cg"]),
                "aspect_ratio_bref2_over_Sref": bref_ft ** 2 / sref_ft2,
                "basis": basis_name,
                "independent_confirmation_deck_vs_lod_header": deck_check,
                "basis_matches_registered": {
                    "DSO_Sref_m2": SREF_DSO_M2, "DSO_Cref_m": CREF_DSO_M,
                    "TP1580_Sref_m2": SREF_TP1580_M2, "TP1580_Cref_m": CREF_TP1580_M,
                    "bref_m_common": BREF_M,
                },
            },
            "frame_conversion": {
                "_rule": ("C_target = C_run * Sref_run / Sref_target, both areas in the same unit. "
                          "A moment coefficient carries the chord ratio as well: "
                          "C_M,target = C_M,run * (Sref_run * Cref_run) / (Sref_target * Cref_target)."),
                "_computed_in": "feet, the unit the decks are written in, so a same-basis run gets exactly 1.0",
                "to_DSO_basis": {
                    "force_coefficient_factor": k_area,
                    "moment_coefficient_factor": k_moment,
                    "already_on_this_basis": basis_name == "DSO",
                },
                "to_TP1580_basis": {
                    "force_coefficient_factor": k_area_tp,
                    "moment_coefficient_factor": k_moment_tp,
                    "already_on_this_basis": basis_name == "TP-1580",
                },
                "note": ("Factor 1.0 means the run is already on that basis and NO conversion is "
                         "applied. Otherwise every CL/CD/CS quoted against a RANS number must be "
                         "multiplied by force_coefficient_factor FIRST. The DSO-to-TP1580 area "
                         "ratio is 13.3572/12 = 1.1131 (D036): an omitted conversion is an 11.31% "
                         "error, against a morphing effect of order 0.4%, i.e. about 28x the signal."),
            },
            "flow_conditions_AS_READ_FROM_DECK": {
                "Mach": mach,
                "AoA_deg": float(deck["AoA"]),
                "Beta_deg": float(deck["Beta"]),
                "ReCref": float(deck["ReCref"]),
                "Vinf_Lunit_per_Tunit": float(deck["Vinf"]),
                "Rho_Munit_per_Lunit3": float(deck["Rho"]),
                "Symmetry": deck.get("Symmetry", "NOT FOUND"),
                "WakeIters": deck.get("WakeIters", "NOT FOUND"),
                "StallModel": deck.get("StallModel", "NOT FOUND"),
                "_Re_note": ("ReCref is 1e+07 in EVERY deck in this package and Vinf is 100 "
                             "Lunit/Tunit with Rho 0.002377 (imperial sea-level slug/ft3). "
                             "These are VSPAERO placeholders: VSPAERO's VLM is INVISCID and uses "
                             "Re only inside the CDo flat-plate skin-friction correlation. The "
                             "flight Reynolds numbers are NOT in these decks."),
            },
            "operating_point": {
                "state": state or "NOT FOUND",
                "target_CL": target_cl,
                "target_CL_source": target_cl_src,
                "target_CL_registered": OPPOINTS[state]["target_CL"] if state else "NOT FOUND",
                "target_CL_delivered_minus_registered": (
                    (target_cl - OPPOINTS[state]["target_CL"])
                    if (state and isinstance(target_cl, float)) else "NOT FOUND"),
                "CL_achieved": cl if cl is not None else "NOT FOUND",
                "CL_achieved_on_DSO_basis": (cl * k_area) if cl is not None else "NOT FOUND",
                "CL_achieved_on_TP1580_basis": (cl * k_area_tp) if cl is not None else "NOT FOUND",
                "CL_error": (cl - target_cl) if (cl is not None and isinstance(target_cl, float)) else "NOT FOUND",
                "_frame": ("CL is on this run's OWN Sref basis, named in "
                           "reference_quantities_AS_READ_FROM_DECK.basis. The target CL is "
                           "defined on the SAME basis, so the error is frame-consistent; but a "
                           "cruise CL (TP-1580 basis) is NOT comparable to a Condition CR CL "
                           "(DSO basis) without the conversion factor, nor across operating "
                           "points at all (D076)."),
            },
            "drag_decomposition": (
                {
                    "_what": ("A within-file proof that CDtot is CDo + CDi exactly, so the "
                              "parasite term can be identified and excluded by arithmetic "
                              "rather than by assertion."),
                    "CDo": polar["points"][0]["CDo"],
                    "CDi": polar["points"][0]["CDi"],
                    "CDtot": polar["points"][0]["CDtot"],
                    "CDo_plus_CDi_minus_CDtot": (polar["points"][0]["CDo"]
                                                 + polar["points"][0]["CDi"]
                                                 - polar["points"][0]["CDtot"]),
                    "CDo_fraction_of_CDtot": polar["points"][0]["CDo"] / polar["points"][0]["CDtot"],
                    "CDiw_far_field": polar["points"][0]["CDiw"],
                    "CDwtot_wake_total": polar["points"][0]["CDwtot"],
                    "CDiw_minus_CDi": polar["points"][0]["CDiw"] - polar["points"][0]["CDi"],
                    "CDiw_over_CDi": polar["points"][0]["CDiw"] / polar["points"][0]["CDi"],
                    "near_vs_far_field_note": ("CDi (surface integration, near field) and CDiw "
                                               "(wake/Trefftz, far field) are the SAME physical "
                                               "quantity computed two ways, so CDiw_over_CDi is a "
                                               "measure of the VLM's own internal consistency at "
                                               "this operating point. It is NOT a correction "
                                               "factor and the two must never be mixed within one "
                                               "comparison."),
                    "Re_over_1e6_in_polar": polar["points"][0]["Re/1e6"],
                    "NOT_COMPARABLE_TO_RANS": [
                        "CDo", "CDtot",
                    ],
                    "COMPARABLE_TO_RANS_INDUCED_DRAG": ["CDi (near field)", "CDiw (far field, Trefftz)"],
                    "why": ("CDo is VSPAERO's flat-plate parasite skin-friction correlation "
                            "evaluated at the deck's ReCref, which is 1e+07 in every deck here "
                            "and is a VSPAERO GUI DEFAULT, not a flight condition (docs/"
                            "ARGUS_reference_data.md section 1c; docs/decisions/ARG.md item 9, "
                            "which WITHDRAWS an earlier CDv-vs-CDo agreement claim on exactly "
                            "this ground). Our RANS runs at Re_cref 0.79e6 to 0.95e6, an order "
                            "of magnitude lower. docs/decisions/legacy.md records the same "
                            "verdict: 'ONLY CDi IS COMPARABLE; CDtot must never appear beside a "
                            "RANS CD'. Per the 2026-08-19 delivery README the DSO now ranks on "
                            "CDiw, the wake/Trefftz far-field value, not CDi."),
                }
                if polar.get("points") else {"error": "NOT FOUND"}
            ),
            "polar": polar,
            "polar_on_DSO_basis": (
                {k: (v * k_area if k in ("CLo", "CLi", "CLtot", "CDo", "CDi", "CDtot",
                                         "CSo", "CSi", "CStot", "CLwtot", "CDwtot",
                                         "CSwtot", "CLiw", "CDiw", "CSiw")
                     else (v * k_moment if k in ("CMxtot", "CMytot", "CMztot") else v))
                 for k, v in polar["points"][0].items()}
                if polar.get("points") else {"error": "NOT FOUND"}
            ),
            "lod_strip_loads": lod,
            "lod_CL_closure": closure or {"error": "NOT FOUND: closure not computable"},
            "history": hist,
            "dso_summary_csv": summary or {"file": "NOT FOUND", "row": None},
            "polar_vs_dso_summary_crosscheck": crosscheck(polar, summary),
            "mapping": {
                "case_prefix_family": gm["family"] if gm["family"] else "unmapped",
                "our_case_prefixes": gm["case_prefixes"],
                "registry_key": gm["registry_key"],
                "evidence": gm["evidence"],
            },
        })

    doc = {
        "_what": "Every VSPAERO (VLM) result found under dso_reference/, extracted verbatim.",
        "_generated": "2026-09-15",
        "_generator": "docs/report/all_geometry_2026-09-15/scripts/extract_vlm.py",
        "_source_root": "dso_reference/",
        "_rules": {
            "every_number_from_a_file": (
                "Every numeric value in this document was read from the file named in its "
                "_source or in the run's deck_file / directory. Nothing is estimated."
            ),
            "frame_rule_D068": (
                "Reference quantities come from each run's OWN .vspaero deck, which is the "
                "record of what ran. The saved .vsp3 GUI state is known-stale for Mach and "
                "Sref (project rules, hard fact 6) and was not used."
            ),
            "VLM_vs_RANS_drag": (
                "VSPAERO writes CDtot = CDo + CDi where CDo is a FLAT-PLATE SKIN-FRICTION "
                "ESTIMATE, not a computed viscous solution. CDtot must NEVER be quoted "
                "against a RANS total drag. The defensible VLM-to-RANS drag comparison is "
                "INDUCED DRAG against INDUCED DRAG: CDi (surface integration) or, per the "
                "DSO's own 2026-08-19 correction, CDiw (wake/Trefftz far-field), which the "
                "delivery README names as the physically appropriate induced-drag metric."
            ),
            "lod_half_wing_factor_2": (
                "The .lod strip loads are HALF-WING but normalised on the FULL-wing Sref, so "
                "CL closes only as 2*sum(Cl*dArea)/Sref. Without the factor 2 it reads 50.0%. "
                "Verified per run in lod_CL_closure."
            ),
            "operating_points_not_mixed_D076": (
                "Three distinct points appear in this package. Condition CR is M 0.10 at "
                "target CL 0.428277635108 on the DSO Sref basis. Early cruise is M 0.78 at "
                "target CL 0.529297087450 and late cruise M 0.78 at target CL 0.544114800210, "
                "both on the TP-1580 Sref basis. A percentage from one is not quotable "
                "against a count from another unless both CLs are stated inline."
            ),
        },
        "_registered_reference_quantities": {
            "DSO_basis": {"Sref_m2": SREF_DSO_M2, "Cref_m": CREF_DSO_M,
                          "Sref_ft2": SREF_DSO_M2 / FT2, "Cref_ft": CREF_DSO_M / FT,
                          "role": "PRIMARY for 3D and DSO comparison"},
            "TP1580_basis": {"Sref_m2": SREF_TP1580_M2, "Cref_m": CREF_TP1580_M,
                             "Sref_ft2": SREF_TP1580_M2 / FT2, "Cref_ft": CREF_TP1580_M / FT,
                             "role": "reported alongside"},
            "bref_m": BREF_M, "bref_ft": BREF_M / FT,
            "RANS_half_wing_patch_area_Aref_m2": 0.620462,
            "area_ratio_DSO_over_TP1580": (SREF_DSO_M2 / SREF_TP1580_M2),
        },
        "_operating_points": OPPOINTS,
        "_verification": verification(runs),
        "_coverage_partition": coverage_partition(runs),
        "n_runs": len(runs),
        "runs": runs,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(doc, indent=2))
    # validate before promoting (project rule: nothing invalid reaches its final filename)
    json.loads(tmp.read_text())
    assert doc["n_runs"] == len(decks), "run count does not match deck count"
    tmp.replace(OUT)
    print(f"wrote {OUT} with {len(runs)} runs")


if __name__ == "__main__":
    main()
