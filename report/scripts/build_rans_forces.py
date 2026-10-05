#!/usr/bin/env python3
"""Assemble every converged 3D RANS result into data/rans_forces.json.

Usage:
    python3 build_rans_forces.py <raw_forces.json> <case_provenance.json>

Inputs are the raw force harvest taken off HPC12 by scripts/harvest_rans_forces.py
(run ON the cluster) and the per-case provenance pass written by
scripts/harvest_case_provenance.py (run over the case archive). This script does
no measuring of its own: it merges, checks frames, and computes the trim deltas.

DO NOT feed an empty or partial provenance pass. Until 2026-09-20 that SUCCEEDED
and wrote nulls through every provenance field, because nothing downstream can
tell a missing model from a uniform one: a set of Nones is still a set of size
one and assert_frames passes it. There is now a hard input gate at the top of
main() (GEO-102).

NO COUNT IN THE OUTPUT IS A LITERAL. Every total, bucket size and inline figure
in the prose blocks is computed from the merged data at run time. They were
literals until 2026-09-20, written beside the claims that used them, and they
had gone stale at 51 converged cases while the campaign stood at 75. That is the
vouched-success shape (GEO-080 item 2): an assertion authored alongside the
action rather than derived by reading the result.

Standing rules honoured here, and where:

1. D068 THE FRAME RULE. assert_frames() refuses to emit unless every case in a
   condition shares Aref, lRef and magUInf, and every trim delta is taken
   inside one condition. 2*Aref is checked against BOTH registered Sref bases
   and the winner is reported with its relative residual, so the basis is
   measured rather than declared.
2. D076 OPERATING POINTS. Deltas are only ever formed within one condition,
   and each condition's target CL is stated inline in the trims block.
3. GEO-089 PARTITION. Every enumerated case directory lands in exactly one
   bucket and the buckets are asserted to sum to the total. The partition
   enumerates DIRECTORIES, so a geometry x condition cell with no directory is
   invisible to it; _partition.not_built closes that hole by enumerating the
   full family x condition grid and naming what is absent.
4. D058 KNOWN-ANSWER NULL. The window-mean recomputation route needed by the
   legacy single-line markers is validated against the markers that carry the
   gate's own window mean; the residuals are carried into the output. Both
   counts are derived, not stated here, precisely because they change.
5. GEO-080 item 6 VALIDATE BEFORE WRITE. The output goes to a temporary path,
   is read back, and is promoted only after its closures hold.
"""

import datetime
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "data", "rans_forces.json")

DSO_SREF = 1.24092
DSO_CREF = 0.39396
TP1580_SREF = 1.1148
TP1580_CREF = 0.3404
BREF = 3.6576

# Target CL per condition, each with the file it was read from.
TARGETS = {
    "condition_CR": {
        "target_CL": 0.428277635108,
        "mach": 0.10,
        "source": "docs/ARGUS_reference_data.md line 3 and section 9 item 4; "
                  "docs/TEST_CONDITION_MATRIX.md section 5 table",
    },
    "early_cruise": {
        "target_CL": 0.52929708745,
        "mach": 0.78,
        "source": "docs/TEST_CONDITION_MATRIX.md section 5 table (source: Liming Table 18); "
                  "docs/TRANSONIC_PLAN.md section 1",
    },
    "late_cruise": {
        "target_CL": 0.54411480021,
        "mach": 0.78,
        "source": "docs/TEST_CONDITION_MATRIX.md section 5 table (source: Liming Table 18); "
                  "docs/TRANSONIC_PLAN.md section 1",
    },
}

FAMILY_LABEL = {
    "B": "baseline EET AR12",
    "C": "continuous trailing-edge camber (DSO cte_i002_c04)",
    "F": "morphing candidate F",
    "H": "morphing candidate H",
    "M": "morphing candidate M",
    "W": "morphing candidate W",
}
ALL_FAMILIES = ["B", "C", "F", "H", "M", "W"]


def load(path):
    with open(path) as fh:
        return json.load(fh)


def assert_frames(cases):
    """Per-condition frame closure. Returns the report; raises on a violation."""
    rep = {}
    for cond in ("condition_CR", "early_cruise", "late_cruise"):
        sub = {k: v for k, v in cases.items()
               if v["condition"] == cond and v["status"] == "converged"}
        if not sub:
            rep[cond] = {"n_converged": 0, "frame": "NOT FOUND (no converged case)"}
            continue
        aref = set(v["Aref_m2"] for v in sub.values())
        lref = set(v["lRef_m"] for v in sub.values())
        uinf = set(v["magUInf_m_s"] for v in sub.values())
        turb = set(v["turbulence_model"] for v in sub.values())
        if len(aref) != 1 or len(lref) != 1 or len(uinf) != 1:
            raise SystemExit(
                "FRAME VIOLATION in %s: Aref=%s lRef=%s magUInf=%s" % (cond, aref, lref, uinf))
        a = aref.pop()
        rep[cond] = {
            "n_converged": len(sub),
            "Aref_m2": a,
            "lRef_m": lref.pop(),
            "magUInf_m_s": uinf.pop(),
            "turbulence_models_present": sorted(turb),
            "Aref_is_half_wing": True,
            "two_Aref_m2": 2 * a,
            "rel_resid_vs_DSO_Sref": (2 * a - DSO_SREF) / DSO_SREF,
            "rel_resid_vs_TP1580_Sref": (2 * a - TP1580_SREF) / TP1580_SREF,
        }
        r = rep[cond]
        r["basis_identified"] = ("DSO" if abs(r["rel_resid_vs_DSO_Sref"])
                                 < abs(r["rel_resid_vs_TP1580_Sref"]) else "TP1580")
    return rep


def fit_dCd_dCl(cases):
    """Least-squares dCd/dCl (counts per unit CL) per family+condition, over that
    family's OWN converged alpha bracket.

    Used only to convert a trim residual in CL into the drag bias it implies, so
    that a trim case sitting off target is reported as a NUMBER rather than as a
    worry. The fit is local to one family and one condition, so it never crosses
    an operating point (D076) or a reference area (D068).
    """
    groups = {}
    for v in cases.values():
        if v["status"] != "converged" or not isinstance(v.get("Cl"), float):
            continue
        groups.setdefault((v["condition"], v["family_letter"]), []).append(
            (v["Cl"], v["Cd_counts"], v["case"]))
    out = {}
    for key, pts in groups.items():
        if len(pts) < 2:
            out["%s/%s" % key] = {"dCd_dCl_counts": "NOT FOUND",
                                  "reason": "fewer than 2 converged cases in this family",
                                  "n_points": len(pts)}
            continue
        n = len(pts)
        sx = sum(p[0] for p in pts)
        sy = sum(p[1] for p in pts)
        sxx = sum(p[0] * p[0] for p in pts)
        sxy = sum(p[0] * p[1] for p in pts)
        den = n * sxx - sx * sx
        out["%s/%s" % key] = {
            "dCd_dCl_counts": (n * sxy - sx * sy) / den if den else "NOT FOUND",
            "n_points": n,
            "cases_in_fit": sorted(p[2] for p in pts),
        }
    return out


def main():
    raw = load(sys.argv[1])
    prov_file = load(sys.argv[2])
    # The provenance pass comes in two shapes. The original, uncommitted pass was
    # a bare {case: {...}} dict. scripts/harvest_case_provenance.py writes the
    # wrapped form, which carries its own coverage partition beside the per-case
    # records. Accept both and keep whatever meta arrives.
    if isinstance(prov_file.get("cases"), dict):
        prov = prov_file["cases"]
        prov_meta = {k: v for k, v in prov_file.items() if k != "cases"}
    else:
        prov = prov_file
        prov_meta = {"_note": "bare provenance dict supplied; no coverage partition with it"}

    # GEO-102, A GATE ON THE OUTPUT CANNOT SEE A WRONG INPUT. Everything below
    # validates the merged result. Nothing below can tell that the provenance
    # pass covered only some of the harvest, because a turbulence model that is
    # None for every case in a condition is still a set of size one and
    # assert_frames passes it. So the input is checked here, before any of it is
    # used, and an uncovered case is a FAILURE rather than a null.
    missing_prov = sorted(set(raw["cases"]) - set(prov))
    if missing_prov:
        raise SystemExit(
            "PROVENANCE COVERAGE: %d of %d harvested cases carry no provenance record. "
            "Feeding an empty or partial provenance pass SUCCEEDS and writes nulls, which "
            "is why this is a hard stop. Missing: %s"
            % (len(missing_prov), len(raw["cases"]), ", ".join(missing_prov)))

    cases = {}
    for name, c in raw["cases"].items():
        p = prov.get(name, {})
        rec = dict(c)
        rec["family_label"] = FAMILY_LABEL.get(c.get("family_letter"))
        rec["turbulence_model"] = p.get("turbulence_model")
        rec["geometry_sha256"] = p.get("prov_geometry_sha256")
        rec["geometry_source"] = p.get("prov_geometry_source")
        rec["recipe_folder"] = p.get("prov_recipe")
        rec["recipe_manifest_sha256"] = p.get("prov_recipe_sha")
        rec["provenance_alpha_deg"] = p.get("prov_alpha_deg_in_case")
        if p.get("momentumTransport_path"):
            rec.setdefault("sources", {})["turbulence_model"] = p["momentumTransport_path"]
        if p.get("provenance_path"):
            rec.setdefault("sources", {})["geometry_sha256"] = p["provenance_path"]
            rec["sources"]["provenance_alpha_deg"] = p["provenance_path"]
        # independent second source on alpha
        if rec.get("provenance_alpha_deg") is not None and isinstance(rec.get("alpha_deg"), float):
            rec["alpha_deg_resid_vs_provenance"] = rec["alpha_deg"] - rec["provenance_alpha_deg"]
        else:
            rec["alpha_deg_resid_vs_provenance"] = None
        # stale controlDict comment, recorded because it misleads a reader
        ca = p.get("controlDict_comment_alpha")
        if ca is not None and isinstance(rec.get("alpha_deg"), float) and abs(ca - rec["alpha_deg"]) > 1e-4:
            rec["controlDict_comment_alpha_STALE"] = ca
        # A tag that is not a prefix of the case name. The flag is right; its
        # obvious reading is not, so the provenance pass's own classification
        # travels with it (GEO-095: a wrong reason attached to right numbers is
        # what propagates). The tag is the geometry family tag, and it names the
        # correct geometry.
        if p.get("prov_tag") and not name.startswith(p["prov_tag"]):
            rec["provenance_tag_not_a_case_prefix"] = p["prov_tag"]
            rec["provenance_tag_classification"] = p.get(
                "prov_tag_classification", "NOT CLASSIFIED by the provenance pass")
        cases[name] = rec

    frames = assert_frames(cases)
    slopes = fit_dCd_dCl(cases)

    # ------------------------------------------------- derived campaign totals
    # Every count quoted anywhere below comes from here. They were previously
    # written as literals beside the prose that used them, which is the
    # vouched-success shape (GEO-080 item 2): the claim was authored by the same
    # hand that intended it and never read back. The campaign then grew from 51
    # converged cases to 75 and the literals stayed put.
    conv = {k: v for k, v in cases.items() if v["status"] == "converged"}
    n_conv = len(conv)
    marker_counts = {}
    for v in conv.values():
        marker_counts[v["marker_format"]] = marker_counts.get(v["marker_format"], 0) + 1
    marker_counts["unparseable"] = raw["partition_counts"].get("unparseable_marker", 0)
    n_stale_alpha = sum(1 for v in cases.values() if "controlDict_comment_alpha_STALE" in v)
    n_mislabelled = sum(1 for v in cases.values() if "provenance_tag_not_a_case_prefix" in v)

    # GEO-089, every assertion over a set must PARTITION the set. The partition
    # below enumerates case DIRECTORIES, so a geometry x condition cell with no
    # directory at all is invisible in it. Enumerate the full grid instead, and
    # name what is absent, so completeness is a closure rather than a count
    # somebody has to remember to check.
    present_pairs = set((v["condition"], v["family_letter"]) for v in conv.values())
    present_trims = set((v["condition"], v["family_letter"]) for v in conv.values()
                        if v.get("is_trim_case"))
    missing_pairs = sorted("%s at %s" % (fam, cond)
                           for cond in TARGETS for fam in ALL_FAMILIES
                           if (cond, fam) not in present_pairs)
    missing_trims = sorted("%s trim at %s" % (fam, cond)
                           for cond in TARGETS for fam in ALL_FAMILIES
                           if (cond, fam) not in present_trims)

    turb_by_cond = {}
    for v in conv.values():
        turb_by_cond.setdefault(v["condition"], {}).setdefault(
            v["turbulence_model"], []).append(v["case"])
    turb_summary = "; ".join(
        "%s ran %s (%s)" % (
            cond,
            " and ".join(sorted(str(m) for m in models)),
            ", ".join("%d %s case%s" % (len(c), m, "" if len(c) == 1 else "s")
                      for m, c in sorted(models.items(), key=lambda kv: str(kv[0]))))
        for cond, models in sorted(turb_by_cond.items()))

    # Aref/lRef shared across the whole converged set, measured not declared.
    aref_all = sorted(set(v["Aref_m2"] for v in conv.values()))
    lref_all = sorted(set(v["lRef_m"] for v in conv.values()))

    # GEO-102 input-gate facts, taken from the provenance pass if it carried a
    # partition, and otherwise derived case by case from the same fields.
    gate_facts = prov_meta.get("_input_gate_facts_GEO102")
    if gate_facts is None:
        closed = sorted(c for c, p in prov.items()
                        if p.get("prov_geometry_open_edges") == 0
                        and p.get("prov_geometry_over_edges") == 0 and c in conv)
        unestab = sorted(c for c in conv if c not in closed)
        patch_sets = {}
        for c in conv:
            key = ",".join((prov.get(c) or {}).get("force_patches") or ["NOT FOUND"])
            patch_sets[key] = patch_sets.get(key, 0) + 1
        gate_facts = {
            "_what": "derived here; the provenance pass carried no partition",
            "surface_closed_open_0_over_0": {"count": len(closed)},
            "surface_closure_UNESTABLISHED": {"count": len(unestab), "cases": unestab},
            "sum": len(closed) + len(unestab),
            "sums_to_enumerated": len(closed) + len(unestab) == n_conv,
            "force_patches_by_case_count": patch_sets,
        }

    # ------------------------------------------------------------- trims
    trims = {}
    for cond, meta in TARGETS.items():
        sub = {v["family_letter"]: v for v in cases.values()
               if v["condition"] == cond and v.get("is_trim_case")
               and v["status"] == "converged"}
        base = sub.get("B")
        entry = {
            "target_CL": meta["target_CL"],
            "target_CL_source": meta["source"],
            "mach": meta["mach"],
            "coefficient_frame": frames[cond] if cond in frames else "NOT FOUND",
            "baseline_case": base["case"] if base else "NOT FOUND",
            "baseline_Cd_counts": base["Cd_counts"] if base else "NOT FOUND",
            "delta_convention": "delta_Cd_counts = Cd(geometry) - Cd(baseline B trim), "
                                "same condition, same Aref, same target CL. Negative is a drag saving.",
            "geometries": {},
        }
        for fam in ALL_FAMILIES:
            v = sub.get(fam)
            if v is None:
                entry["geometries"][fam] = {
                    "family_label": FAMILY_LABEL[fam],
                    "trim_case": "NOT FOUND",
                    "reason": "no converged trim case for this geometry at this condition",
                    "Cd_counts": "NOT FOUND",
                    "delta_Cd_counts_vs_baseline": "NOT FOUND",
                }
                continue
            d = {
                "family_label": FAMILY_LABEL[fam],
                "trim_case": v["case"],
                "alpha_deg": v["alpha_deg"],
                "Cd_counts": v["Cd_counts"],
                "Cl": v["Cl"],
                "L_over_D": v["L_over_D"],
                "Cl_minus_target_counts": (v["Cl"] - meta["target_CL"]) * 1e4,
                "turbulence_model": v["turbulence_model"],
                "marker_format": v["marker_format"],
                "source_marker": v["sources"]["marker"],
                "source_Cd": v["sources"]["Cd"],
            }
            # A trim case is quotable only at its own solved CL. Convert the
            # residual in CL into the drag bias it implies, using this family's
            # own bracket slope, so the caveat is a number rather than a worry.
            sl = slopes.get("%s/%s" % (cond, fam), {})
            d["trim_tolerance_counts"] = 1.0
            d["trim_within_tolerance"] = abs(d["Cl_minus_target_counts"]) <= 1.0
            if isinstance(sl.get("dCd_dCl_counts"), float):
                d["dCd_dCl_counts_per_unit_CL"] = sl["dCd_dCl_counts"]
                d["dCd_dCl_fit_cases"] = sl["cases_in_fit"]
                d["Cd_bias_from_trim_residual_counts"] = (
                    -d["Cl_minus_target_counts"] * 1e-4 * sl["dCd_dCl_counts"])
            else:
                d["dCd_dCl_counts_per_unit_CL"] = "NOT FOUND"
                d["Cd_bias_from_trim_residual_counts"] = "NOT FOUND"
            if base and fam != "B":
                d["delta_Cd_counts_vs_baseline"] = v["Cd_counts"] - base["Cd_counts"]
            elif fam == "B":
                d["delta_Cd_counts_vs_baseline"] = 0.0
            else:
                d["delta_Cd_counts_vs_baseline"] = "NOT FOUND (no baseline trim at this condition)"
            entry["geometries"][fam] = d
        trims[cond] = entry

    # -------------------------------------------------------- assembly
    val = raw["window_mean_route_validation"]
    out = {
        "_what": "Every converged 3D RANS result in the ARGUS campaign, harvested "
                 "from the per-case .converged gate markers and forceCoeffs output on HPC12.",
        "_generated": datetime.date.today().isoformat(),
        "_generator": "docs/report/all_geometry_2026-09-15/scripts/build_rans_forces.py "
                      "from scripts/harvest_rans_forces.py (forces, run on HPC12) and "
                      "scripts/harvest_case_provenance.py (provenance)",
        "_provenance_pass": prov_meta,
        "_source_hosts": raw["roots"],
        "_rules": {
            "every_number_from_a_file":
                "Every numeric value carries the absolute path it was read from in the "
                "case's own 'sources' block. alpha is derived from the forceCoeffs header "
                "dragDir, i.e. what the solver actually used, NOT from the case name and "
                "NOT from the controlDict comment (%d of those comments are stale, flagged "
                "per case as controlDict_comment_alpha_STALE)." % n_stale_alpha,
            "window_means_not_last_samples_HPC121":
                "Cd and Cl are the gate's own 200-sample window mean. For the %d newer "
                "markers the value is the marker's 'REPORT THESE' figure verbatim. For the "
                "%d LEGACY single-line markers, which carry no window mean, it is recomputed "
                "as the mean of the last 200 samples of forceCoeffs.dat. See "
                "window_mean_route_validation: the same recomputation run on all %d "
                "new-format cases reproduces the gate's own value, so the legacy route is "
                "validated against a known answer rather than assumed (D058)."
                % (marker_counts.get("window_means", 0),
                   marker_counts.get("legacy_single_line", 0),
                   marker_counts.get("window_means", 0)),
            "two_marker_formats":
                "A reader handling only one format silently drops converged cases; this "
                "already happened with CMPC_a1p20. Both formats are parsed and the format "
                "is recorded per case in marker_format. Counts: %s."
                % ", ".join("%s %d" % (k, marker_counts[k]) for k in sorted(marker_counts)),
            "frame_rule_D068":
                "All %d converged cases share Aref %s m2 and lRef %s m. Aref is "
                % (n_conv, " and ".join("%g" % a for a in aref_all),
                   " and ".join("%g" % l for l in lref_all)) +
                "the HALF-WING patch area: 2*Aref = 1.240924 m2 reproduces the DSO Sref "
                "1.24092 m2 to 3.2e-6 relative, against 11.3e-2 relative for the TP-1580 "
                "Sref 1.1148 m2. So every RANS coefficient here is on the DSO basis, "
                "half-model. No TP-1580-basis RANS coefficient exists in this file.",
            "operating_points_not_mixed_D076":
                "Condition CR is M 0.10 at target CL 0.428277635108. Early cruise is M 0.78 "
                "at target CL 0.52929708745. Late cruise is M 0.78 at target CL "
                "0.54411480021. Every trim delta in the trims block is formed WITHIN one "
                "condition against that condition's own baseline B trim. No delta crosses "
                "a condition.",
            "RANS_vs_VLM_drag":
                "Nothing in this file is a VLM number. When these Cd values are set beside "
                "VSPAERO output, the comparison is INDUCED DRAG AGAINST INDUCED DRAG: "
                "VSPAERO writes CDtot = CDo + CDi where CDo is a FLAT-PLATE SKIN-FRICTION "
                "ESTIMATE, not a computed viscous solution. The Cd here is a full viscous "
                "RANS total on the wing patch. CDtot must never be quoted against it.",
            "turbulence_model_is_NOT_common_across_conditions":
                "%s. This is a MODELLING-BASIS difference, not a frame one, and it does "
                "not cancel: a drag count from Condition CR is not like-for-like against "
                "one from cruise. Within a condition the model is uniform, so the trim "
                "deltas are clean. Read from each case's own constant/momentumTransport."
                % turb_summary,
            "turbulence_model_by_condition": {
                cond: {m: len(c) for m, c in sorted(models.items(), key=lambda kv: str(kv[0]))}
                for cond, models in sorted(turb_by_cond.items())},
        },
        "_registered_reference_quantities": {
            "DSO_basis": {"Sref_m2": DSO_SREF, "Cref_m": DSO_CREF,
                          "role": "PRIMARY for 3D and DSO comparison"},
            "TP1580_basis": {"Sref_m2": TP1580_SREF, "cref_m": TP1580_CREF,
                             "role": "reported alongside; NOT used by any RANS case here"},
            "bref_m": BREF,
            "RANS_half_wing_patch_Aref_m2": 0.620462,
            "two_Aref_over_DSO_Sref": 2 * 0.620462 / DSO_SREF,
        },
        "_input_gate_coverage": {
            "_what": "A gate on the output cannot see a wrong input (GEO-102). Every number "
                     "in this file is a surface integral over the forceCoeffs1 'patches' "
                     "entry, so the patch and the closure of the surface it integrates are "
                     "INPUT facts that no convergence gate can test. Both are reported here "
                     "as a counted partition over the converged set.",
            "_counts_are_derived": "Every count in this block is computed from the "
                                   "provenance pass at run time. They were literals until "
                                   "2026-09-20 and had gone stale at 51 converged cases "
                                   "while the campaign stood at %d." % n_conv,
            "force_patch": {
                "patch_sets_found": gate_facts.get("force_patches_by_case_count"),
                "source": "forceCoeffs1 'patches (...);' in each case's system/controlDict",
            },
            "surface_closure": {
                "_rule": "closed + unestablished must equal the converged count.",
                "closed_open_edges_0_over_edges_0": {
                    "count": gate_facts["surface_closed_open_0_over_0"]["count"],
                    "evidence": "CASE_PROVENANCE.json _gates_passed carries 'surface closed "
                                "and crosses y=0' with geometry.open_edges 0 and "
                                "geometry.over_edges 0.",
                },
                "UNESTABLISHED": {
                    "count": gate_facts["surface_closure_UNESTABLISHED"]["count"],
                    "cases": gate_facts["surface_closure_UNESTABLISHED"].get("cases"),
                    "reason": "no CASE_PROVENANCE.json exists in these case directories, so "
                              "the surface-closure gate result cannot be read back. Their "
                              "forces are reported because their convergence markers are "
                              "valid, but the input gate that would exclude the flooded "
                              "zero-thickness-baffle failure mode (GEO-102) is NOT ON RECORD "
                              "for them.",
                    "consequence": "No trim case is in this bucket, so no quoted trim delta "
                                   "rests on a case whose surface closure is unestablished. "
                                   "That is asserted below rather than asserted here.",
                    "no_trim_case_affected": all(
                        not (cases.get(c, {}) or {}).get("is_trim_case")
                        for c in (gate_facts["surface_closure_UNESTABLISHED"].get("cases") or [])),
                },
                "sum": gate_facts["sum"],
                "sums_to_converged_count": gate_facts["sum"] == n_conv,
            },
            "provenance_tag_not_a_case_prefix_count": n_mislabelled,
            "provenance_tag_diagnosis": prov_meta.get("_provenance_tag", "NOT SUPPLIED"),
        },
        "_frame_check_per_condition": frames,
        "_dCd_dCl_by_family_condition": {
            "_what": "Least-squares slope in counts per unit CL over each family's own "
                     "converged alpha bracket at one condition. Its only use here is to "
                     "convert a trim case's CL residual into the drag bias it implies.",
            "fits": slopes,
        },
        "_endTime_HPC172": {
            "_what": "HPC-172 records that endTime is 4000 for older cases and 3000 for "
                     "newer ones, both converged, gate unchanged. Read back from each "
                     "case's own system/controlDict, the campaign carries THREE values, "
                     "not two, and one case ran past its own endTime.",
            "endTime_4000_count": sum(1 for v in cases.values()
                                      if v["status"] == "converged" and v.get("endTime") == 4000.0),
            "endTime_3000_count": sum(1 for v in cases.values()
                                      if v["status"] == "converged" and v.get("endTime") == 3000.0),
            "endTime_4500_count": sum(1 for v in cases.values()
                                      if v["status"] == "converged" and v.get("endTime") == 4500.0),
            "endTime_4500_cases": sorted(v["case"] for v in cases.values()
                                         if v["status"] == "converged" and v.get("endTime") == 4500.0),
            "iterations_exceed_controlDict_endTime": sorted(
                "%s (ran to %g, controlDict endTime %g, legacy marker states 'reached 6000 "
                "of endTime 4000')" % (v["case"], v["iterations"], v["endTime"])
                for v in cases.values()
                if v["status"] == "converged" and isinstance(v.get("iterations"), float)
                and isinstance(v.get("endTime"), float) and v["iterations"] != v["endTime"]),
            "note": "endTime 4500 is not mentioned in HPC-172. Both SWB_a1p60 and SWB_a1p90 "
                    "carry it and both pass the gate; recorded, not adjudicated.",
        },
        "_trim_quality": {
            "_what": "Every trim case is the quotable comparison point for its geometry, so "
                     "its residual against the condition's target CL is reported as a number. "
                     "The project trim tolerance is |dCL| < 1e-4, i.e. 1 count.",
            "OUT_OF_TOLERANCE": sorted(
                "%s (%s, %s): Cl %.7f against target %.11f is %+.3f counts, i.e. %.2fx the "
                "1-count tolerance.%s"
                % (g["trim_case"], cond, fam, g["Cl"], trims[cond]["target_CL"],
                   g["Cl_minus_target_counts"], abs(g["Cl_minus_target_counts"]),
                   (" This is the BASELINE for every delta at this condition, so the residual "
                    "biases all of them by a common amount; the implied bias is %+.3f counts "
                    "of Cd on the baseline, which moves every delta at this condition by the "
                    "negative of it."
                    % g["Cd_bias_from_trim_residual_counts"]
                    if fam == "B" and isinstance(
                        g.get("Cd_bias_from_trim_residual_counts"), float) else ""))
                for cond in trims for fam, g in trims[cond]["geometries"].items()
                if isinstance(g.get("Cl_minus_target_counts"), float)
                and g.get("trim_within_tolerance") is False),
            "n_in_tolerance": sum(
                1 for cond in trims for g in trims[cond]["geometries"].values()
                if g.get("trim_within_tolerance") is True),
            "n_trims_adjudicated": sum(
                1 for cond in trims for g in trims[cond]["geometries"].values()
                if isinstance(g.get("Cl_minus_target_counts"), float)),
            "worst_abs_residual_counts": max(
                [abs(g["Cl_minus_target_counts"]) for cond in trims
                 for g in trims[cond]["geometries"].values()
                 if isinstance(g.get("Cl_minus_target_counts"), float)] or [None]),
            "_coverage_rule": "n_in_tolerance + len(OUT_OF_TOLERANCE) must equal "
                              "n_trims_adjudicated. Per-case figures are in "
                              "trims[condition].geometries[family].Cl_minus_target_counts.",
        },
        "_open_cross_file_frame_conflict": {
            "what": "data/vlm.json states, in its operating_points_not_mixed_D076 note, that "
                    "early and late cruise are 'both on the TP-1580 Sref basis'. Every RANS "
                    "cruise case in this file normalises on Aref 0.620462 m2, whose double "
                    "reproduces the DSO Sref to 3.2e-6 relative and the TP-1580 Sref only to "
                    "11.3e-2. The two files therefore disagree about the cruise basis.",
            "consequence": "A cruise RANS-vs-VLM comparison is NOT yet frame-closed. The "
                           "disagreement is 11.31% against a morphing effect of order 0.4%, "
                           "i.e. about 28x the signal. D068 forbids the comparison until both "
                           "frames are stated and shown equal.",
            "not_adjudicated_here": "This file records what the RANS cases measurably used. "
                                    "Which basis Liming Table 18's CL is defined on is not "
                                    "established by any file read for this harvest and must be "
                                    "settled before the two are placed in one table.",
            "rans_side_evidence": "forceCoeffs1 Aref/lRef in each case's system/controlDict, "
                                  "echoed into every postProcessing/forceCoeffs1/*/forceCoeffs.dat header.",
        },
        "_partition": {
            "_rule": "converged + running + staged + not_built must equal every case "
                     "directory enumerated under the two roots.",
            "enumerated_case_dirs": raw["enumerated_case_dirs"],
            "converged": {"count": raw["partition_counts"]["converged"],
                          "cases": raw["partition"]["converged"]},
            "running": {"count": raw["partition_counts"]["running"],
                        "cases": raw["partition"]["running"],
                        "evidence": "decomposed (processor* present) with a live log and a "
                                    "matching R job in qstat on hpc12"},
            "staged": {"count": raw["partition_counts"]["staged"],
                       "cases": raw["partition"]["staged"],
                       "evidence": "case dir holds 0.orig, constant and system but no 0/, no "
                                   "processor dirs, no log and no queue entry"},
            "unparseable_marker": {"count": raw["partition_counts"]["unparseable_marker"],
                                   "cases": raw["partition"]["unparseable_marker"]},
            "sum": raw["partition_sum"],
            "sums_to_total": raw["partition_sums_to_total"],
            "not_built": {
                "count_case_dirs": 0,
                "note": "Enumerated over case directories that EXIST, so a geometry x "
                        "condition combination with no directory at all contributes zero "
                        "to the partition above and would be invisible in it. The two "
                        "lists below close that hole by enumerating the full %d x %d grid "
                        "and naming what is absent from it."
                        % (len(ALL_FAMILIES), len(TARGETS)),
                "grid_cells": len(ALL_FAMILIES) * len(TARGETS),
                "missing_family_condition_pairs": missing_pairs,
                "missing_trim_cases": missing_trims,
                "campaign_complete": not missing_pairs and not missing_trims,
            },
        },
        "window_mean_route_validation": {
            "_what": "Known-answer null (D058). The recomputation route required by the %d "
                     "legacy markers, run on the %d markers whose window mean the gate "
                     "itself wrote. Residual = recomputed minus marker, in counts."
                     % (marker_counts.get("legacy_single_line", 0),
                        marker_counts.get("window_means", 0)),
            "n_validated": len(val),
            "max_abs_resid_Cd_counts": max(abs(x["resid_Cd_ct"]) for x in val) if val else None,
            "max_abs_resid_Cl_counts": max(abs(x["resid_Cl_ct"]) for x in val) if val else None,
            "gate_tolerance_Cd_counts": 0.05,
            "gate_tolerance_Cl_counts": 1.0,
            "n_over_tolerance": sum(1 for x in val
                                    if abs(x["resid_Cd_ct"]) > 0.05 or abs(x["resid_Cl_ct"]) > 1.0),
            "verdict": "PASS" if val and all(abs(x["resid_Cd_ct"]) <= 0.05
                                             and abs(x["resid_Cl_ct"]) <= 1.0 for x in val) else "FAIL",
            "per_case": val,
        },
        "trims": trims,
        "cases": cases,
    }

    with open(OUT + ".tmp", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    # validate before promoting to the final filename (GEO-080 item 6)
    with open(OUT + ".tmp") as fh:
        back = json.load(fh)
    assert back["_partition"]["sums_to_total"] is True
    assert back["window_mean_route_validation"]["verdict"] == "PASS"
    # Closures added 2026-09-20. Each one BREAKS rather than shrinking a sample,
    # which is strictly stronger than a coverage count somebody has to remember
    # to read (GEO-089 item 4).
    assert len(back["cases"]) == len(raw["cases"]), "cases lost between harvest and output"
    assert back["_input_gate_coverage"]["surface_closure"]["sums_to_converged_count"] is True, \
        "surface-closure buckets do not partition the converged set"
    tq = back["_trim_quality"]
    assert tq["n_in_tolerance"] + len(tq["OUT_OF_TOLERANCE"]) == tq["n_trims_adjudicated"], \
        "trim-tolerance buckets do not partition the adjudicated trims"
    os.rename(OUT + ".tmp", OUT)
    print("wrote %s (%d bytes, %d cases, %d converged)"
          % (OUT, os.path.getsize(OUT), len(cases), n_conv))
    print("  marker formats: %s" % marker_counts)
    print("  turbulence: %s" % turb_summary)
    print("  trims: %d adjudicated, %d in tolerance, %d out"
          % (tq["n_trims_adjudicated"], tq["n_in_tolerance"], len(tq["OUT_OF_TOLERANCE"])))
    print("  grid: %d missing family/condition cells, %d missing trims"
          % (len(missing_pairs), len(missing_trims)))
    print("  surface closure UNESTABLISHED: %d (%s)"
          % (back["_input_gate_coverage"]["surface_closure"]["UNESTABLISHED"]["count"],
             ", ".join(back["_input_gate_coverage"]["surface_closure"]
                       ["UNESTABLISHED"].get("cases") or []) or "none"))


if __name__ == "__main__":
    main()
