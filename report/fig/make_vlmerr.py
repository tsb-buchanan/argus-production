#!/usr/bin/env python3
"""
RANS versus VLM disagreement, quantified.  ARGUS all-geometry report, 2026-09-15.

Answers: can VSPAERO (VLM) be trusted as an optimisation tool, and how far off is it.

EVERY NUMBER IN EVERY FIGURE IS READ FROM ONE OF:
    data/rans_forces.json    (RANS, HPC12 harvest of .converged window means)
    data/vlm.json            (VSPAERO decks, .polar, .lod)
    data/dso_vlm_claims.json (DSO report claims; used only for the noise-floor band)
Nothing is typed in from memory.  Where a quantity does not exist in any file the
script writes NOT FOUND into the audit dump and the panel says so on its face.

FRAME (D068).  Every coefficient plotted here is on the DSO basis,
Sref 1.24092 m2 / Cref 0.39396 m, bref 3.6576 m.
  RANS  : forceCoeffs Aref 0.620462 m2 is the HALF-WING patch; 2*Aref reproduces the
          DSO Sref to 3.2e-6 relative (11.3e-2 against TP-1580), so the half-model
          coefficient IS the full-wing coefficient.  Verified in rans_forces.json.
  VSPAERO: every Condition CR deck carries Sref 13.3572 ft2 = 1.2409245 m2, the DSO
          basis, so frame_conversion.to_DSO_basis.force_coefficient_factor is exactly 1.
  The script asserts both rather than trusting them.

OPERATING POINT (D076).  Condition CR only: M 0.10, target C_L 0.428277635108.
No cruise panel exists; see the GAPS block printed at the end for why.

DRAG PHYSICS (task rule 4).  VSPAERO CDtot = CDo + CDi with CDo a flat-plate
skin-friction correlation evaluated at the deck's ReCref = 1e7, which is a GUI default
and is not our Reynolds number.  CDtot is NEVER plotted against a RANS drag here.
Only CDi (near field, surface) and CDiw (far field, Trefftz) are used on the VLM side.
The RANS side has no induced-drag decomposition at all (NOT FOUND, see GAPS), so every
drag panel compares an INCREMENT at fixed C_L, not an absolute level.
"""

import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)
DATA = os.path.join(REPORT, "data")
OUT = HERE

RANS_JSON = os.path.join(DATA, "rans_forces.json")
VLM_JSON = os.path.join(DATA, "vlm.json")
DSO_CLAIMS_JSON = os.path.join(DATA, "dso_vlm_claims.json")

CONDITION = "condition_CR"
CT = 1.0e4  # 1 drag count = 1e-4 in a force coefficient

# ---------------------------------------------------------------- style
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 10.5,
    "axes.labelsize": 10.5,
    "axes.titlesize": 10.5,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 8.8,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "figure.dpi": 120,
})

# Okabe-Ito, colour-blind safe.  Every series ALSO carries a distinct marker and a
# distinct dash pattern or hatch, so the figures survive greyscale printing.
OKABE = {
    "black":  "#000000",
    "orange": "#E69F00",
    "sky":    "#56B4E9",
    "green":  "#009E73",
    "blue":   "#0072B2",
    "verm":   "#D55E00",
    "purple": "#CC79A7",
}
GEOM_STYLE = {
    "B": dict(c=OKABE["black"],  m="o", ls="-",             hatch=""),
    "C": dict(c=OKABE["blue"],   m="s", ls="--",            hatch="//"),
    "M": dict(c=OKABE["green"],  m="^", ls="-.",            hatch="\\\\"),
    "F": dict(c=OKABE["verm"],   m="D", ls=(0, (1, 1)),     hatch="xx"),
    "W": dict(c=OKABE["purple"], m="v", ls=(0, (5, 1, 1, 1)), hatch=".."),
    "H": dict(c=OKABE["orange"], m="P", ls=(0, (3, 1, 1, 1, 1, 1)), hatch="++"),
}
MORPHS = ["C", "M", "F", "W", "H"]

AUDIT = {}          # everything the figures display, dumped to JSON beside them
GAPS = []           # stated on the figures and returned to the caller


def note(msg):
    GAPS.append(msg)


# ---------------------------------------------------------------- load
with open(RANS_JSON) as fh:
    RANS = json.load(fh)
with open(VLM_JSON) as fh:
    VLM = json.load(fh)
with open(DSO_CLAIMS_JSON) as fh:
    DSO_CLAIMS = json.load(fh)

# ---------------------------------------------------------------- frame assertions
rframe = RANS["_frame_check_per_condition"][CONDITION]
assert rframe["basis_identified"] == "DSO", rframe
assert rframe["rel_resid_vs_DSO_Sref"] < 1e-5, rframe
assert rframe["Aref_is_half_wing"] is True, rframe
RANS_AREF = rframe["Aref_m2"]
RANS_2AREF = rframe["two_Aref_m2"]
RANS_TURB = rframe["turbulence_models_present"]

TARGET_CL = RANS["trims"][CONDITION]["target_CL"]
MACH = RANS["trims"][CONDITION]["mach"]

AUDIT["frame"] = {
    "rans_Aref_m2": RANS_AREF,
    "rans_two_Aref_m2": RANS_2AREF,
    "rans_rel_resid_vs_DSO_Sref": rframe["rel_resid_vs_DSO_Sref"],
    "rans_rel_resid_vs_TP1580_Sref": rframe["rel_resid_vs_TP1580_Sref"],
    "rans_turbulence_models": RANS_TURB,
    "target_CL": TARGET_CL,
    "mach": MACH,
    "basis": "DSO  Sref 1.24092 m2 / Cref 0.39396 m / bref 3.6576 m",
}

# ---------------------------------------------------------------- RANS side
rtrim = RANS["trims"][CONDITION]["geometries"]
rcases = RANS["cases"]

rans = {}
for letter in ["B"] + MORPHS:
    g = rtrim[letter]
    case = g["trim_case"]
    rec = rcases[case]
    assert abs(rec["Aref_m2"] - RANS_AREF) < 1e-12, (case, rec["Aref_m2"])
    rans[letter] = {
        "case": case,
        "Cd_counts": g["Cd_counts"],
        "Cl": g["Cl"],
        "alpha_deg": g["alpha_deg"],
        "dCd_counts": g["delta_Cd_counts_vs_baseline"],
        "Cl_minus_target_counts": g["Cl_minus_target_counts"],
        "Cd_bias_from_trim_residual_counts": g["Cd_bias_from_trim_residual_counts"],
        "trim_within_tolerance": g["trim_within_tolerance"],
        "turbulence_model": g["turbulence_model"],
        "drift_Cd_ct": rec.get("drift_Cd_ct"),
        "drift_Cl_ct": rec.get("drift_Cl_ct"),
        "window_samples": rec.get("window_samples_used"),
    }

# RANS lift slope: least squares over each family's OWN converged alpha bracket at this
# condition.  Derived here because rans_forces.json carries dCd/dCl but no dCl/dalpha.
lift_slope = {}
for letter in ["B"] + MORPHS:
    pts = sorted(
        (v["alpha_deg"], v["Cl"], k)
        for k, v in rcases.items()
        if v.get("status") == "converged"
        and v.get("condition") == CONDITION
        and v.get("family_letter") == letter
    )
    if len(pts) < 2:
        lift_slope[letter] = {"dCl_dalpha_per_deg": "NOT FOUND", "n": len(pts)}
        continue
    a = np.array([p[0] for p in pts])
    cl = np.array([p[1] for p in pts])
    A = np.vstack([a, np.ones_like(a)]).T
    (slope, icept), *_ = np.linalg.lstsq(A, cl, rcond=None)
    resid = cl - (slope * a + icept)
    lift_slope[letter] = {
        "dCl_dalpha_per_deg": float(slope),
        "alpha0_deg": float(-icept / slope),
        "rms_resid_Cl": float(np.sqrt(np.mean(resid ** 2))),
        "n": len(pts),
        "cases": [p[2] for p in pts],
        "alpha_range_deg": [float(a.min()), float(a.max())],
    }
AUDIT["rans_lift_slope_fit"] = lift_slope
AUDIT["rans_trim"] = rans

# ---------------------------------------------------------------- VLM side
vruns = {}
for r in VLM["runs"]:
    if r["operating_point"]["state"] != CONDITION:
        continue
    vruns.setdefault(r["geometry_name"], r)

# Map by the extractor's own checksum mapping (registry + CASE_PROVENANCE), never by name.
vlm = {}
vlm_unmapped = []
for name, r in vruns.items():
    prefixes = r["mapping"].get("our_case_prefixes") or []
    sw = [p for p in prefixes if p.startswith("SW")]
    if not sw:
        vlm_unmapped.append((name, r["mapping"].get("case_prefix_family")))
        continue
    letter = sw[0][2:]
    assert letter in GEOM_STYLE, (name, sw)
    fc = r["frame_conversion"]["to_DSO_basis"]["force_coefficient_factor"]
    assert fc == 1.0, (name, fc)          # already DSO basis; no conversion applied
    assert r["reference_quantities_AS_READ_FROM_DECK"]["basis"] == "DSO", name
    dd = r["drag_decomposition"]
    assert abs(dd["CDo_plus_CDi_minus_CDtot"]) < 1e-9, name
    vlm[letter] = {
        "geometry_name": name,
        "package": r["delivery_package"],
        "sha256": r["geometry_vsp3"]["sha256"],
        "alpha_deg": r["flow_conditions_AS_READ_FROM_DECK"]["AoA_deg"],
        "alpha_deg_csv": float(r["dso_summary_csv"]["row"]["alpha_trim_deg"]),
        "CL": r["operating_point"]["CL_achieved"],
        "CL_error": r["operating_point"]["CL_error"],
        "CDi": dd["CDi"],
        "CDiw": dd["CDiw_far_field"],
        "CDo": dd["CDo"],
        "CDtot": dd["CDtot"],
        "CDo_frac_of_CDtot": dd["CDo_fraction_of_CDtot"],
        "ReCref_deck": r["flow_conditions_AS_READ_FROM_DECK"]["ReCref"],
        "mapping_family": r["mapping"]["case_prefix_family"],
    }

missing = [l for l in ["B"] + MORPHS if l not in vlm]
assert not missing, "VLM run missing for %s" % missing

# The SECOND baseline.  baseline_corrected is the post-S3 aerofoil-interpolation
# lineage the DSO ranks against; it is NOT the surface our RANS meshes.
alt = vruns.get("baseline_corrected")
assert alt is not None
ALT_B = {
    "geometry_name": "baseline_corrected",
    "sha256": alt["geometry_vsp3"]["sha256"],
    "alpha_deg": alt["flow_conditions_AS_READ_FROM_DECK"]["AoA_deg"],
    "CL": alt["operating_point"]["CL_achieved"],
    "CL_error": alt["operating_point"]["CL_error"],
    "CDi": alt["drag_decomposition"]["CDi"],
    "CDiw": alt["drag_decomposition"]["CDiw_far_field"],
    "mapping_family": alt["mapping"]["case_prefix_family"],
}
AUDIT["vlm_runs_mapped"] = vlm
AUDIT["vlm_alternative_baseline"] = ALT_B
AUDIT["vlm_unmapped_condition_CR"] = vlm_unmapped

# ---------------------------------------------------------------- increments
def deltas(base_CDi, base_CDiw):
    out = {}
    for l in MORPHS:
        out[l] = {
            "dCDi_ct": (vlm[l]["CDi"] - base_CDi) * CT,
            "dCDiw_ct": (vlm[l]["CDiw"] - base_CDiw) * CT,
        }
    return out


D_MAP = deltas(vlm["B"]["CDi"], vlm["B"]["CDiw"])          # vs checksum-mapped B
D_ALT = deltas(ALT_B["CDi"], ALT_B["CDiw"])                # vs baseline_corrected

BASELINE_SHIFT_CDi_ct = (vlm["B"]["CDi"] - ALT_B["CDi"]) * CT
BASELINE_SHIFT_CDiw_ct = (vlm["B"]["CDiw"] - ALT_B["CDiw"]) * CT

AUDIT["vlm_deltas_vs_mapped_B"] = D_MAP
AUDIT["vlm_deltas_vs_baseline_corrected"] = D_ALT
AUDIT["vlm_baseline_choice_shift_counts"] = {
    "CDi": BASELINE_SHIFT_CDi_ct,
    "CDiw": BASELINE_SHIFT_CDiw_ct,
    "meaning": ("every VLM increment moves by this much when the reference is switched "
                "from baseline_wing_only_refined (the meshed surface) to "
                "baseline_corrected (the surface the DSO ranks against)"),
}

# Independent check that baseline_corrected IS the DSO's reference: recompute the
# 19 Aug published CDiw reductions from vlm.json alone.
AUDIT["reproduce_DSO_19aug_CDiw_percent_vs_baseline_corrected"] = {
    l: (vlm[l]["CDiw"] / ALT_B["CDiw"] - 1.0) * 100.0 for l in MORPHS
}
AUDIT["reproduce_superseded_CDi_percent_vs_baseline_corrected"] = {
    l: (vlm[l]["CDi"] / ALT_B["CDi"] - 1.0) * 100.0 for l in MORPHS
}

# ---------------------------------------------------------------- noise floor
def find_key(obj, key):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == key and isinstance(v, str):
                return v
            got = find_key(v, key)
            if got is not None:
                return got
    elif isinstance(obj, list):
        for v in obj:
            got = find_key(v, key)
            if got is not None:
                return got
    return None


floor_txt = find_key(DSO_CLAIMS, "our_independent_agreement")
assert floor_txt, "noise-floor statement not found in dso_vlm_claims.json"
m = re.search(r"noise floor of ([0-9.]+) to ([0-9.]+) counts", floor_txt)
assert m, floor_txt
FLOOR_LO, FLOOR_HI = float(m.group(1)), float(m.group(2))
m2 = re.search(r"acceptance band .*? is ([0-9.]+) counts", floor_txt)
ACCEPT = float(m2.group(1)) if m2 else None
AUDIT["noise_floor"] = {
    "counts_lo": FLOOR_LO, "counts_hi": FLOOR_HI,
    "preregistered_acceptance_band_counts": ACCEPT,
    "source_field": "dso_vlm_claims.json claims[].our_independent_agreement",
    "source_text": floor_txt,
}

# The morphing signal: the span the RANS actually measures over the six geometries.
rans_d = np.array([rans[l]["dCd_counts"] for l in ["B"] + MORPHS], dtype=float)
SIGNAL_LO, SIGNAL_HI = float(rans_d.min()), float(rans_d.max())
AUDIT["rans_morphing_signal_span_counts"] = [SIGNAL_LO, SIGNAL_HI]

# Gate drift, the within-dataset convergence noise indicator.
AUDIT["rans_gate_drift_Cd_counts"] = {l: rans[l]["drift_Cd_ct"] for l in ["B"] + MORPHS}

# ---------------------------------------------------------------- rank statistics
def ranks(values):
    order = sorted(values, key=lambda k: values[k])
    return {k: i + 1 for i, k in enumerate(order)}


def spearman_kendall(a, b, keys):
    ra = np.array([a[k] for k in keys], float)
    rb = np.array([b[k] for k in keys], float)
    n = len(keys)
    d = ra - rb
    rho = 1.0 - 6.0 * np.sum(d ** 2) / (n * (n ** 2 - 1))
    conc = disc = 0
    swapped = []
    for i in range(n):
        for j in range(i + 1, n):
            s = (ra[i] - ra[j]) * (rb[i] - rb[j])
            if s > 0:
                conc += 1
            elif s < 0:
                disc += 1
                swapped.append((keys[i], keys[j]))
    tau = (conc - disc) / (conc + disc) if (conc + disc) else float("nan")
    return float(rho), float(tau), conc, disc, swapped


RANS_RANK = ranks({l: rans[l]["dCd_counts"] for l in MORPHS})
VLM_RANK_NEAR = ranks({l: D_MAP[l]["dCDi_ct"] for l in MORPHS})
VLM_RANK_FAR = ranks({l: D_MAP[l]["dCDiw_ct"] for l in MORPHS})

rho_n, tau_n, cn, dn, swap_n = spearman_kendall(VLM_RANK_NEAR, RANS_RANK, MORPHS)
rho_f, tau_f, cf, df_, swap_f = spearman_kendall(VLM_RANK_FAR, RANS_RANK, MORPHS)
AUDIT["ranking"] = {
    "rans_rank_by_dCd_total": RANS_RANK,
    "vlm_rank_by_dCDi_near": VLM_RANK_NEAR,
    "vlm_rank_by_dCDiw_far": VLM_RANK_FAR,
    "near_field": dict(spearman=rho_n, kendall=tau_n, concordant=cn, discordant=dn,
                       swapped_pairs=swap_n),
    "far_field": dict(spearman=rho_f, kendall=tau_f, concordant=cf, discordant=df_,
                      swapped_pairs=swap_f),
    "rank_is_baseline_invariant": ("an increment ranking is invariant to the choice of "
                                   "VLM baseline: switching it shifts every increment by "
                                   "the same constant"),
}

# Unresolved RANS rank pairs: separations below the measured noise floor.
unresolved = []
for i in range(len(MORPHS)):
    for j in range(i + 1, len(MORPHS)):
        a, b = MORPHS[i], MORPHS[j]
        sep = abs(rans[a]["dCd_counts"] - rans[b]["dCd_counts"])
        if sep < FLOOR_HI:
            unresolved.append({"pair": [a, b], "separation_counts": sep})
AUDIT["rans_rank_pairs_below_noise_floor"] = unresolved

# ---------------------------------------------------------------- errors
err = {}
for l in MORPHS:
    rd = rans[l]["dCd_counts"]
    en = D_MAP[l]["dCDi_ct"] - rd
    ef = D_MAP[l]["dCDiw_ct"] - rd
    err[l] = {
        "rans_dCd_total_ct": rd,
        "vlm_dCDi_near_ct": D_MAP[l]["dCDi_ct"],
        "vlm_dCDiw_far_ct": D_MAP[l]["dCDiw_ct"],
        "abs_err_near_ct": en,
        "abs_err_far_ct": ef,
        "rel_err_near_pct": 100.0 * en / rd,
        "rel_err_far_pct": 100.0 * ef / rd,
        "rel_err_ill_conditioned": abs(rd) < FLOOR_HI,
        "vlm_over_rans_near": D_MAP[l]["dCDi_ct"] / rd,
    }
AUDIT["errors_vs_mapped_B"] = err

# Trim-residual sensitivity on the RANS side (file-supplied bias, not a model).
bias_shift = {
    l: (rans[l]["Cd_bias_from_trim_residual_counts"]
        - rans["B"]["Cd_bias_from_trim_residual_counts"])
    for l in MORPHS
}
AUDIT["rans_delta_shift_if_trim_residual_removed_counts"] = bias_shift

# VLM trim-residual sensitivity.  MODEL, not a file read: CDi ~ CL^2 at fixed geometry,
# so dCDi/dCL = 2 CDi / CL.  Marked as a model inside the result (GEO-095).
vlm_trim_sens = {}
for l in ["B"] + MORPHS:
    v = vlm[l]
    vlm_trim_sens[l] = {
        "CL_error_counts": v["CL_error"] * CT,
        "implied_CDi_bias_counts_MODEL_CL_squared": 2.0 * v["CDi"] / v["CL"] * (-v["CL_error"]) * CT,
    }
AUDIT["vlm_trim_residual_sensitivity_MODEL"] = vlm_trim_sens

# ---------------------------------------------------------------- lift comparison
lift = {}
for l in ["B"] + MORPHS:
    da = vlm[l]["alpha_deg"] - rans[l]["alpha_deg"]
    s = lift_slope[l]["dCl_dalpha_per_deg"]
    lift[l] = {
        "alpha_vlm_deg": vlm[l]["alpha_deg"],
        "alpha_rans_deg": rans[l]["alpha_deg"],
        "d_alpha_deg": da,
        "implied_dCL_at_fixed_alpha": da * s if isinstance(s, float) else "NOT FOUND",
        "implied_dCL_pct_of_target": (100.0 * da * s / TARGET_CL) if isinstance(s, float) else "NOT FOUND",
    }
alt_da = ALT_B["alpha_deg"] - rans["B"]["alpha_deg"]
lift["B_alternative_baseline_corrected"] = {
    "alpha_vlm_deg": ALT_B["alpha_deg"],
    "alpha_rans_deg": rans["B"]["alpha_deg"],
    "d_alpha_deg": alt_da,
    "implied_dCL_at_fixed_alpha": alt_da * lift_slope["B"]["dCl_dalpha_per_deg"],
    "implied_dCL_pct_of_target": 100.0 * alt_da * lift_slope["B"]["dCl_dalpha_per_deg"] / TARGET_CL,
}
AUDIT["lift"] = lift

# ---------------------------------------------------------------- gaps
note("RANS induced drag is NOT FOUND. rans_forces.json carries no CDi, no CDiw and no "
     "pressure/viscous split for any 3D case, and qoi_methods.md records that no consumer "
     "of the archived argusTrefftz planes exists. Every drag panel therefore compares an "
     "INCREMENT at fixed C_L (VLM induced-drag increment against RANS TOTAL-drag "
     "increment), not CDi against CDi. A true induced-against-induced correlation cannot "
     "be drawn from the present data set.")
note("VLM lift slope dCL/dalpha is NOT FOUND. All 23 VSPAERO runs carry exactly one polar "
     "point and no stability deck (.stab), so no alpha derivative can be read or fitted. "
     "The lift panel therefore compares TRIM ALPHA at the common target C_L, which needs "
     "no derivative, and the lift-slope panel shows the RANS side alone.")
note("No cruise RANS-vs-VLM comparison is drawn, for two independent reasons. (1) Geometry: "
     "of the 14 cruise VSPAERO runs only rigid_baseline maps to any of our cases by "
     "checksum, and it maps to baseline_corrected, which is NOT the meshed surface; the "
     "other six cruise geometries are a separate set with no derived STL. (2) Frame: "
     "rans_forces.json records an open, unadjudicated conflict, the cruise RANS cases "
     "measurably normalising on the DSO Sref while the cruise decks carry the TP-1580 "
     "Sref, an 11.31 percent gap against a 0.4 percent signal. D068 forbids the table "
     "until it is settled.")
note("Late cruise has no converged trim case for any geometry, so no late-cruise point "
     "of any kind exists.")

FIGS = []


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    size = os.path.getsize(path)
    FIGS.append((path, size))
    print("wrote %s  %d bytes" % (path, size))
    return path


# Label offsets chosen so no two annotations collide; F and W sit within 0.6 ct of
# each other on both drag panels, which is the whole point of the F/W tie.
LABEL_OFF = {
    "dCDi_ct":  {"C": (7, -4), "M": (7, 3), "F": (-15, -11), "W": (4, -13), "H": (8, -2)},
    "dCDiw_ct": {"C": (-16, -11), "M": (7, 2), "F": (-15, -10), "W": (4, -13), "H": (7, -2)},
}


# ================================================================ FIGURE 1
def _corr_panel(ax, key, xlabel, title, shift_ct):
    xs_f = [D_MAP[l][key] for l in MORPHS]
    xs_o = [D_ALT[l][key] for l in MORPHS]
    ys = [rans[l]["dCd_counts"] for l in MORPHS]
    lo = min(xs_f + xs_o + ys + [0.0])
    hi = max(xs_f + xs_o + ys + [0.0])
    pad = 0.10 * (hi - lo) + 0.4
    lim = (lo - pad, hi + pad)
    ax.axhspan(-FLOOR_HI, FLOOR_HI, color="0.88", zorder=0)
    ax.axhline(0.0, color="0.75", lw=0.7, zorder=0)
    ax.axvline(0.0, color="0.75", lw=0.7, zorder=0)
    ax.plot(lim, lim, color="0.40", ls="-", lw=1.0, zorder=1)
    for l in MORPHS:
        s = GEOM_STYLE[l]
        ax.plot(D_ALT[l][key], rans[l]["dCd_counts"], marker=s["m"], ms=6.2,
                mfc="none", mec=s["c"], mew=1.1, ls="none", zorder=2)
        ax.plot(D_MAP[l][key], rans[l]["dCd_counts"], marker=s["m"], ms=7.2,
                mfc=s["c"], mec="black", mew=0.6, ls="none", zorder=3)
        ax.annotate(l, (D_MAP[l][key], rans[l]["dCd_counts"]),
                    textcoords="offset points", xytext=LABEL_OFF[key][l],
                    fontsize=9.5, zorder=4)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"RANS $\Delta C_{D}$ (total, viscous), counts")
    ax.set_title(title.split(" ")[0], loc="left", fontsize=10.0)
    ax.tick_params(direction="in", top=True, right=True)


def figure_correlation():
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.0), layout="constrained")
    ax_near, ax_far, ax_alpha, ax_slope = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    _corr_panel(ax_near, "dCDi_ct", r"VLM $\Delta C_{Di}$ (near field), counts",
                r"(a) VLM near-field $C_{Di}$", BASELINE_SHIFT_CDi_ct)
    _corr_panel(ax_far, "dCDiw_ct", r"VLM $\Delta C_{Diw}$ (far field), counts",
                r"(b) VLM far-field $C_{Diw}$", BASELINE_SHIFT_CDiw_ct)

    # ---- (c) trim alpha at the common target C_L
    ax = ax_alpha
    xs = [vlm[l]["alpha_deg"] for l in ["B"] + MORPHS] + [ALT_B["alpha_deg"]]
    ys = [rans[l]["alpha_deg"] for l in ["B"] + MORPHS]
    lo = min(xs + ys) - 0.10
    hi = max(xs + ys) + 0.16
    ax.plot([lo, hi], [lo, hi], color="0.40", lw=1.0, zorder=1)
    ax.plot(ALT_B["alpha_deg"], rans["B"]["alpha_deg"], marker="o", ms=6.2,
            mfc="none", mec="black", mew=1.1, ls="none", zorder=3)
    ax.annotate("B$^{\\ast}$", (ALT_B["alpha_deg"], rans["B"]["alpha_deg"]),
                textcoords="offset points", xytext=(-24, -4), fontsize=9.5)
    for l in ["B"] + MORPHS:
        s = GEOM_STYLE[l]
        off = {"B": (7, -4), "C": (7, -3), "M": (-15, -4), "F": (7, -5),
               "W": (-16, -4), "H": (7, -1)}[l]
        ax.plot(vlm[l]["alpha_deg"], rans[l]["alpha_deg"], marker=s["m"], ms=7.2,
                mfc=s["c"], mec="black", mew=0.6, ls="none", zorder=3)
        ax.annotate(l, (vlm[l]["alpha_deg"], rans[l]["alpha_deg"]),
                    textcoords="offset points", xytext=off, fontsize=9.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"VLM trim $\alpha$, deg")
    ax.set_ylabel(r"RANS trim $\alpha$, deg")
    worst = max(MORPHS, key=lambda l: abs(lift[l]["d_alpha_deg"]))
    ax.set_title("(c)",
                 loc="left", fontsize=8.8)
    ax.tick_params(direction="in", top=True, right=True)
    # ---- (d) lift slope, RANS side only
    ax = ax_slope
    letters = ["B"] + MORPHS
    vals = np.array([lift_slope[l]["dCl_dalpha_per_deg"] for l in letters])
    yerr = np.array([
        lift_slope[l]["rms_resid_Cl"]
        / (lift_slope[l]["alpha_range_deg"][1] - lift_slope[l]["alpha_range_deg"][0])
        for l in letters])
    x = np.arange(len(letters), dtype=float)
    for i, l in enumerate(letters):
        s = GEOM_STYLE[l]
        ax.errorbar(x[i], vals[i], yerr=yerr[i], marker=s["m"], ms=7.2, mfc=s["c"],
                    mec="black", mew=0.6, ecolor="0.35", elinewidth=0.9, capsize=3,
                    ls="none", zorder=3)
    span = vals.max() - vals.min()
    ax.set_xticks(x)
    ax.set_xticklabels(letters)
    ax.set_xlim(-0.6, len(letters) - 0.4)
    ax.set_ylim(vals.min() - yerr.max() - 0.35 * span,
                vals.max() + yerr.max() + 1.55 * span)
    ax.set_xlabel("geometry")
    ax.set_ylabel(r"RANS $\mathrm{d}C_L/\mathrm{d}\alpha$, per deg")
    ax.set_title("(d)", loc="left", fontsize=10.0)
    ax.tick_params(direction="in", top=True, right=True)
    ax.grid(axis="y", ls=":", lw=0.6, color="0.82")
    ax.set_axisbelow(True)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    # ONE LINE IN THE FIGURE; the rest belongs in the caption (project decision, 2026-09-22).
    # This was a seven-line boxed paragraph explaining why no VLM line is drawn. The
    # reader needs to know there is no VLM counterpart; the reason it is absent is prose.
    pass  # "no VLM counterpart" is stated in prose (22 Sep)

    handles = [
        Line2D([], [], color="0.40", lw=1.0, label="1:1"),
        Line2D([], [], marker="o", ls="none", mfc="0.35", mec="black", ms=6.4,
               label="VLM referenced to mapped B (the meshed surface)"),
        Line2D([], [], marker="o", ls="none", mfc="none", mec="0.35", mew=1.1, ms=6.0,
               label="VLM referenced to baseline_corrected (the DSO's own reference)"),
        Patch(facecolor="0.88", label="RANS noise floor $\\pm$%.3f ct" % FLOOR_HI),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False,
               fontsize=8.4, handlelength=1.8, columnspacing=1.4)
    fig.suptitle("", fontsize=9.4)
    return save(fig, "fig_vlmerr_correlation.pdf")


# ================================================================ FIGURE 2
def figure_ranking():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.3), sharey=True,
                             layout="constrained")
    for ax, vrank, vname, stats, swaps, title in (
        (axes[0], VLM_RANK_NEAR, "VLM $\\Delta C_{Di}$\n(near field)",
         (rho_n, tau_n, len(swap_n)), swap_n, "(a) near-field $C_{Di}$"),
        (axes[1], VLM_RANK_FAR, "VLM $\\Delta C_{Diw}$\n(far field, Trefftz)",
         (rho_f, tau_f, len(swap_f)), swap_f, "(b) far-field $C_{Diw}$"),
    ):
        for l in MORPHS:
            s = GEOM_STYLE[l]
            ax.plot([0, 1], [vrank[l], RANS_RANK[l]], color=s["c"], ls=s["ls"],
                    lw=1.8, marker=s["m"], ms=7.5, mfc=s["c"], mec="black", mew=0.6,
                    clip_on=False, zorder=3)
            ax.annotate(l, (0, vrank[l]), textcoords="offset points",
                        xytext=(-17, -4), fontsize=10, annotation_clip=False)
            ax.annotate(l, (1, RANS_RANK[l]), textcoords="offset points",
                        xytext=(10, -4), fontsize=10, annotation_clip=False)
        ax.set_xlim(-0.12, 1.06)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([vname, "RANS $\\Delta C_D$\n(total, viscous)"], fontsize=9.2)
        ax.set_ylim(5.5, 0.5)
        ax.set_yticks([1, 2, 3, 4, 5])
        ax.set_title(title, loc="left", fontsize=9.4)  # tau and swaps are in the caption
        ax.tick_params(axis="y", direction="in", right=True, pad=22)
        ax.tick_params(axis="x", length=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        pass  # swapped pairs are in the caption and the table (22 Sep)
    axes[0].set_ylabel("rank  (1 = least drag)")
    tie = ", ".join("%s/%s %.3f ct" % (p["pair"][0], p["pair"][1], p["separation_counts"])
                    for p in unresolved)
    fig.suptitle("", fontsize=9.2)
    return save(fig, "fig_vlmerr_ranking.pdf")


# ================================================================ FIGURE 3
def figure_counts():
    fig = plt.figure(figsize=(7.2, 6.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0])
    ax_bar = fig.add_subplot(gs[0, :])
    ax_abs = fig.add_subplot(gs[1, 0])
    ax_rel = fig.add_subplot(gs[1, 1])

    x = np.arange(len(MORPHS), dtype=float)
    w = 0.26

    # ---- (a) the three increments
    ax_bar.axhspan(-FLOOR_HI, FLOOR_HI, color="0.87", zorder=0)
    ax_bar.axhspan(SIGNAL_LO, SIGNAL_HI, facecolor="none", edgecolor="0.45",
                   lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax_bar.bar(x - w, [rans[l]["dCd_counts"] for l in MORPHS], width=w,
               color=OKABE["black"], edgecolor="black", linewidth=0.6, zorder=3,
               label="RANS $\\Delta C_D$ (total, viscous)")
    ax_bar.bar(x, [D_MAP[l]["dCDi_ct"] for l in MORPHS], width=w,
               color=OKABE["blue"], hatch="//", edgecolor="black", linewidth=0.6,
               zorder=3, label="VLM $\\Delta C_{Di}$ (near field)")
    ax_bar.bar(x + w, [D_MAP[l]["dCDiw_ct"] for l in MORPHS], width=w,
               color=OKABE["orange"], hatch="xx", edgecolor="black", linewidth=0.6,
               zorder=3, label="VLM $\\Delta C_{Diw}$ (far field)")
    ax_bar.bar(x, [D_ALT[l]["dCDi_ct"] for l in MORPHS], width=w, facecolor="none",
               edgecolor=OKABE["blue"], linewidth=1.2, linestyle=":", zorder=4,
               label="VLM $\\Delta C_{Di}$ referenced to baseline_corrected")
    ax_bar.axhline(0.0, color="black", lw=0.8, zorder=5)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(MORPHS)
    ax_bar.set_xlabel("geometry")
    ax_bar.set_ylabel("drag increment vs B, counts")
    ax_bar.set_title("(a)", loc="left", fontsize=9.6)
    ax_bar.tick_params(direction="in", top=True, right=True)
    h, _ = ax_bar.get_legend_handles_labels()
    h += [Patch(facecolor="0.87", label="RANS noise floor $\\pm$%.3f ct" % FLOOR_HI),
          Patch(facecolor="none", edgecolor="0.45", ls=(0, (4, 3)),
                label="RANS signal span %.2f to %.2f ct" % (SIGNAL_LO, SIGNAL_HI))]
    ylo = min(min(D_ALT[l]["dCDi_ct"] for l in MORPHS), SIGNAL_LO) - 1.0
    yhi = max(SIGNAL_HI, max(rans[l]["dCd_counts"] for l in MORPHS))
    ax_bar.set_ylim(ylo, yhi + 0.58 * (yhi - ylo))
    ax_bar.legend(handles=h, loc="upper left", frameon=False, fontsize=8.0,
                  ncol=2, handlelength=1.7, columnspacing=1.1, labelspacing=0.32)
    ax_bar.text(0.985, 0.035,
                "",
                transform=ax_bar.transAxes, ha="right", va="bottom", fontsize=7.8,
                color="0.20", bbox=dict(fc="white", ec="0.75", lw=0.5, pad=2.0))

    # ---- (b) absolute error
    e_near = [err[l]["abs_err_near_ct"] for l in MORPHS]
    e_far = [err[l]["abs_err_far_ct"] for l in MORPHS]
    ax_abs.axhspan(-FLOOR_HI, FLOOR_HI, color="0.87", zorder=0,
                   label="noise floor $\\pm$%.3f ct" % FLOOR_HI)
    ax_abs.bar(x - w / 2, e_near, width=w, color=OKABE["blue"], hatch="//",
               edgecolor="black", linewidth=0.6, label="near field", zorder=3)
    ax_abs.bar(x + w / 2, e_far, width=w, color=OKABE["orange"], hatch="xx",
               edgecolor="black", linewidth=0.6, label="far field", zorder=3)
    ax_abs.axhline(0.0, color="black", lw=0.8, zorder=5)
    ax_abs.axhline(SIGNAL_HI - SIGNAL_LO, color="0.35", ls=(0, (4, 3)), lw=1.0,
                   zorder=2, label="full signal span %.2f ct" % (SIGNAL_HI - SIGNAL_LO))
    ax_abs.set_xticks(x)
    ax_abs.set_xticklabels(MORPHS)
    ax_abs.set_xlabel("geometry")
    ax_abs.set_ylabel("VLM $-$ RANS, counts")
    ax_abs.set_title("(b)", loc="left", fontsize=10.0)
    ax_abs.tick_params(direction="in", top=True, right=True)
    lo_b = min(e_far + e_near)
    hi_b = SIGNAL_HI - SIGNAL_LO
    ax_abs.set_ylim(lo_b - 0.08 * (hi_b - lo_b), hi_b + 0.42 * (hi_b - lo_b))
    ax_abs.legend(loc="upper left", frameon=False, fontsize=7.8, ncol=1,
                  handlelength=1.6, columnspacing=1.0, labelspacing=0.3)

    # ---- (c) relative error, as the ratio of predicted to measured increment
    r_near = [err[l]["vlm_dCDi_near_ct"] / err[l]["rans_dCd_total_ct"] for l in MORPHS]
    r_far = [err[l]["vlm_dCDiw_far_ct"] / err[l]["rans_dCd_total_ct"] for l in MORPHS]
    ill = [err[l]["rel_err_ill_conditioned"] for l in MORPHS]
    cn_ = ["0.72" if i else OKABE["blue"] for i in ill]
    cf_ = ["0.72" if i else OKABE["orange"] for i in ill]
    ax_rel.axhline(1.0, color="0.25", ls="-", lw=1.0, zorder=2)
    ax_rel.axhline(0.0, color="black", lw=0.8, zorder=5)
    ax_rel.bar(x - w / 2, r_near, width=w, color=cn_, hatch="//", edgecolor="black",
               linewidth=0.6, zorder=3)
    ax_rel.bar(x + w / 2, r_far, width=w, color=cf_, hatch="xx", edgecolor="black",
               linewidth=0.6, zorder=3)
    allr = r_near + r_far
    span = max(allr) - min(allr)
    ax_rel.set_ylim(min(allr) - 0.24 * span, max(allr) + 0.52 * span)
    for i in range(len(MORPHS)):
        for xx, v in ((x[i] - w / 2, r_near[i]), (x[i] + w / 2, r_far[i])):
            up = v >= 0
            ax_rel.annotate("%+.2f" % v, (xx, v), textcoords="offset points",
                            xytext=(0, 4 if up else -4), ha="center",
                            va="bottom" if up else "top",
                            fontsize=6.9, rotation=90)
    ax_rel.set_xticks(x)
    ax_rel.set_xticklabels(MORPHS)
    ax_rel.set_xlabel("geometry")
    ax_rel.set_ylabel(r"VLM $\Delta$ / RANS $\Delta$")
    ax_rel.set_title("(c)", loc="left", fontsize=10.0)
    ax_rel.tick_params(direction="in", top=True, right=True)
    ax_rel.text(0.975, 0.975,
                "",
                transform=ax_rel.transAxes, va="top", ha="right", fontsize=7.4,
                bbox=dict(fc="white", ec="0.75", lw=0.6, pad=2.5))

    fig.suptitle("", fontsize=9.2)
    return save(fig, "fig_vlmerr_counts.pdf")


# ================================================================ run
p1 = figure_correlation()
p2 = figure_ranking()
p3 = figure_counts()

AUDIT["_gaps"] = GAPS
AUDIT["_figures"] = [{"path": p, "bytes": s} for p, s in FIGS]
AUDIT["_what"] = ("Every value displayed by fig_vlmerr_*.pdf, dumped so a reader can "
                  "check the figure against the source JSONs without rerunning it.")
AUDIT["_sources"] = {"rans": RANS_JSON, "vlm": VLM_JSON, "dso_claims": DSO_CLAIMS_JSON}

audit_path = os.path.join(OUT, "fig_vlmerr_values.json")
with open(audit_path, "w") as fh:
    json.dump(AUDIT, fh, indent=1, sort_keys=True, default=str)
print("wrote %s  %d bytes" % (audit_path, os.path.getsize(audit_path)))

# verify on disk, not on exit status
bad = [p for p, s in FIGS if (not os.path.isfile(p)) or s == 0]
if bad:
    print("FAILED to write: %s" % bad, file=sys.stderr)
    sys.exit(1)

print("\n--- GAPS ---")
for g in GAPS:
    print("*", g)
print("\n--- headline numbers ---")
print("RANS dCd (ct):", {l: round(rans[l]["dCd_counts"], 4) for l in MORPHS})
print("VLM dCDi (ct) vs mapped B:", {l: round(D_MAP[l]["dCDi_ct"], 4) for l in MORPHS})
print("VLM dCDiw (ct) vs mapped B:", {l: round(D_MAP[l]["dCDiw_ct"], 4) for l in MORPHS})
print("VLM dCDi (ct) vs baseline_corrected:", {l: round(D_ALT[l]["dCDi_ct"], 4) for l in MORPHS})
print("baseline shift CDi / CDiw (ct): %.4f / %.4f" % (BASELINE_SHIFT_CDi_ct, BASELINE_SHIFT_CDiw_ct))
print("near field  rho %.3f tau %.3f swaps %s" % (rho_n, tau_n, swap_n))
print("far  field  rho %.3f tau %.3f swaps %s" % (rho_f, tau_f, swap_f))
print("abs err near (ct):", {l: round(err[l]["abs_err_near_ct"], 4) for l in MORPHS})
print("abs err far  (ct):", {l: round(err[l]["abs_err_far_ct"], 4) for l in MORPHS})
print("rel err near (%):", {l: round(err[l]["rel_err_near_pct"], 1) for l in MORPHS})
print("d_alpha (deg):", {l: round(lift[l]["d_alpha_deg"], 4) for l in ["B"] + MORPHS})
print("d_alpha B* (deg): %.4f" % lift["B_alternative_baseline_corrected"]["d_alpha_deg"])
print("RANS dCl/dalpha (/deg):", {l: round(lift_slope[l]["dCl_dalpha_per_deg"], 6) for l in ["B"] + MORPHS})
print("DSO 19Aug CDiw %% reproduced:",
      {l: round(v, 4) for l, v in AUDIT["reproduce_DSO_19aug_CDiw_percent_vs_baseline_corrected"].items()})
print("superseded CDi %% reproduced:",
      {l: round(v, 4) for l, v in AUDIT["reproduce_superseded_CDi_percent_vs_baseline_corrected"].items()})
print("unresolved RANS pairs:", unresolved)
print("noise floor: %.3f to %.3f ct, acceptance band %s ct" % (FLOOR_LO, FLOOR_HI, ACCEPT))
