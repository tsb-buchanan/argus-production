#!/usr/bin/env python3
"""
ARGUS all-geometry report, 2026-09-15: SPANWISE LOADING, RANS against VSPAERO VLM.

Self-contained. Reads only:
    ../data/vlm.json                              .lod strip loads, 69 strips per run
    ../data/rans_forces.json                      window-mean C_L, alpha, frame per case
    ../data/spanwise/<CASE>_spanwise_dense.csv    dense RANS sectional loads  [preferred]
    ../data/spanwise/<CASE>_spanwise_dense.json   its provenance and CL closure
    ../data/sections/<CASE>_eta<NNN>.csv          5-station Cp cuts           [fallback]
Writes vector PDFs into this directory.

RULES ENFORCED IN CODE, NOT IN PROSE
------------------------------------
 1. EVERY NUMBER COMES FROM A FILE. No aerodynamic value is typed in this script. Reference
    areas, reference chords, target C_L, alpha, strip loads and sectional loads are all read.
    The only literals are plotting constants and definitions (a drag count is 1e-4).

 2. FRAME RULE (D068). The plotted quantity is the SECTIONAL LOAD c_l*c/Cref. c_l is on the
    LOCAL chord on both sides; Cref is the DSO reference chord, read from each VLM deck and
    from the RANS spanwise provenance, and the two are ASSERTED equal before anything is
    differenced. Every VLM run used here is checked to carry Sref 13.3572 ft2 (DSO basis)
    from its own deck. The cruise VLM runs carry Sref 12 ft2 (TP-1580) AND are a different
    geometry set, so they are excluded by an explicit partition, never silently dropped.

 3. THE HALF-WING FACTOR OF TWO (frame-rule instance 1). The .lod strips are HALF-WING but
    normalised on the FULL-wing Sref, so CL closes only as 2*sum(Cl*dArea)/Sref. This script
    recomputes that closure from the strip table for every run it plots, BOTH with and
    without the factor, and refuses to plot a run whose with-factor closure exceeds a
    pre-registered tolerance. The RANS side carries the same factor in its own closure,
    2*int(c_l c dy)/Sref, computed by the extractor and re-read here.

 4. OPERATING POINTS NOT MIXED (D076). Condition CR (M 0.10, target C_L 0.428277635108) and
    early cruise (M 0.78, target C_L 0.529297087450) get separate files. Each figure states
    its own Mach and target C_L inline. Late cruise has no trim case and appears nowhere.

 5. WINDOW MEANS, NEVER LAST SAMPLES (HPC-121). Every RANS C_L quoted for a closure is the
    value rans_forces.json took from that case's .converged gate marker.

 6. PARTITIONS CLOSE (GEO-089). Plotted + excluded-with-reason == enumerated, asserted, for
    the VLM runs and for the RANS cases at each condition.

WHAT IS PLOTTED, AND WHY THE SECOND ROW EXISTS
----------------------------------------------
Row 1 is the absolute load c_l*c/Cref against eta = y/(b/2). The six geometries differ by a
few per cent there and the curves lie on top of one another: the absolute plot shows the
planform, not the morphing. Row 2 is the DIFFERENCE from the baseline of the same method, on
the same eta grid, which is where the morphing effect actually lives.

The baseline is interpolated onto each candidate's own eta grid before differencing, because
the VLM panel centroids MOVE when the surface is morphed (up to 2.1e-03 ft here). The offset
and its bound on the difference are printed by this script and belong in the caption.
"""

import csv
import glob
import json
import math
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

FIG_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(FIG_DIR, "..", "data"))
VLM_JSON = os.path.join(DATA_DIR, "vlm.json")
RANS_JSON = os.path.join(DATA_DIR, "rans_forces.json")
SPANWISE_DIR = os.path.join(DATA_DIR, "spanwise")
SECTIONS_DIR = os.path.join(DATA_DIR, "sections")

# Pre-registered tolerance on the .lod half-wing closure. The extraction phase measured
# max |relative error| 4.578e-06 over all 23 runs and attributed it to the .lod's 5-decimal
# print of Cl and dArea. 1e-04 is two decades above that: it passes print round-off and
# fails a dropped or doubled factor, which is the thing the gate is for.
LOD_CLOSURE_TOL = 1.0e-4

# ----------------------------------------------------------------------------------
# style: serif, colour-blind-safe (Okabe-Ito), distinct dashes for greyscale survival
# ----------------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 9.0,
    "axes.labelsize": 9.5,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 7.6,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 1.1,
    "legend.frameon": True,
    "legend.framealpha": 1.0,
    "legend.edgecolor": "0.4",
    "legend.fancybox": False,
    "legend.borderpad": 0.4,
    "pdf.fonttype": 42,
    "savefig.dpi": 600,
})

ORDER = ["B", "C", "F", "H", "M", "W"]

STYLE = {
    "B": dict(color="#000000", ls="-",                           marker="o"),
    "C": dict(color="#0072B2", ls=(0, (5, 1.6)),                 marker="s"),
    "F": dict(color="#D55E00", ls=(0, (1.2, 1.2)),               marker="^"),
    "H": dict(color="#009E73", ls=(0, (4, 1.4, 1, 1.4)),         marker="D"),
    "M": dict(color="#CC79A7", ls=(0, (3, 1.2, 1, 1.2, 1, 1.2)), marker="v"),
    "W": dict(color="#E69F00", ls=(0, (7, 1.6, 1, 1.6)),         marker="P"),
}

GREY = "0.35"


# ==================================================================================
# loaders
# ==================================================================================
def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def _vlm_strip_series(r):
    """One VLM run reduced to the plotted quantities, with its own closure recomputed."""
    rq = r["reference_quantities_AS_READ_FROM_DECK"]
    assert rq["basis"] == "DSO", "%s is not on the DSO basis" % r["run_id"]
    st = r["lod_strip_loads"]["strip_table"]
    col = {c: k for k, c in enumerate(st["_columns"])}
    A = np.asarray(st["rows"], dtype=float)

    Yavg = A[:, col["Yavg_ft"]]
    chord = A[:, col["Chord_ft"]]
    dArea = A[:, col["dArea_ft2"]]
    cl = A[:, col["Cl"]]

    # HALF-WING FACTOR OF TWO, recomputed here from the strip table rather than read from
    # the extractor's own field. Deriving the assertion from the state is the point
    # (GEO-080): a closure quoted by the party that wrote it vouches for itself.
    Sref_ft2 = rq["Sref_ft2"]
    CL_x2 = float(2.0 * np.sum(cl * dArea) / Sref_ft2)
    CL_x1 = 0.5 * CL_x2
    CLtot = float(r["polar"]["points"][0]["CLtot"])
    rel_x2 = (CL_x2 - CLtot) / CLtot
    assert abs(rel_x2) <= LOD_CLOSURE_TOL, (
        "%s: .lod closure WITH the factor 2 is %+.3e, outside the pre-registered %.0e. "
        "Not plotted." % (r["run_id"], rel_x2, LOD_CLOSURE_TOL))

    bref_ft = rq["Bref_ft"]
    Cref_ft = rq["Cref_ft"]
    return dict(
        geometry=r["geometry_name"], run_id=r["run_id"],
        delivery_package=r["delivery_package"],
        eta=Yavg / (0.5 * bref_ft), load=cl * chord / Cref_ft, cl=cl, chord_ft=chord,
        Yavg_ft=Yavg, dArea_ft2=dArea, n_strips=int(r["lod_strip_loads"]["n_strips"]),
        Sref_ft2=Sref_ft2, Cref_ft=Cref_ft, bref_ft=bref_ft, basis=rq["basis"],
        CL_with_factor_2=CL_x2, CL_without_factor_2=CL_x1, CLtot_polar=CLtot,
        rel_err_with_factor_2=rel_x2, ratio_without_factor_2=CL_x1 / CLtot,
        CL_target=r["operating_point"]["target_CL"],
        CL_error=r["operating_point"]["CL_error"],
        sha256=r["geometry_vsp3"].get("sha256"),
        mapping_note=r["mapping"].get("case_prefix_family", ""))


def vlm_condition_cr(vlm):
    """The Condition-CR VLM runs that map onto one of our six RANS families.

    Returns (series, reference_baseline, partition). The partition enumerates EVERY run in
    vlm.json and puts each into exactly one named bucket, because a verdict over an unstated
    subset is not a verdict (GEO-087).

    TWO BASELINE VLM RUNS EXIST AND THEY ARE NOT THE SAME WING.

      `baseline_wing_only_refined`  IS the surface every SWB/CMPB/LCB RANS case meshes
                                    (registry baseline::geometry_3d_laddercap, status
                                    SANCTIONED_PRODUCTION_BASELINE, verified bit-for-bit
                                    from this .vsp3 and from no other).
      `baseline_corrected`          is the corrected-aerofoil-interpolation lineage. Its
                                    laddercap STL carries NOT_USED_overruled_by_ARG-099 and
                                    "DO NOT MESH ... No case has ever used it".

    ALL FIVE CANDIDATES ARE ON THE CORRECTED LINEAGE. So a candidate differenced against
    `baseline_wing_only_refined` carries MORPHING PLUS LINEAGE, and the two do not separate.
    D071's question decides it: does the defect differ between the two things being
    differenced? Here it does, so it enters the number and must be audited rather than
    recorded. This function returns BOTH baselines; the figure differences the VLM
    candidates against the LINEAGE-CONSISTENT one and draws the lineage offset separately,
    so the reader sees the term the RANS panel cannot remove.
    """
    series, buckets = {}, {"plotted": [], "reference_baseline_other_lineage": [],
                           "not_condition_CR": [], "unmapped": [],
                           "duplicate_redelivery_byte_identical": []}
    bcorr = None
    for r in vlm["runs"]:
        gid = r["geometry_name"] + "/" + r["operating_point"]["state"]
        if r["operating_point"]["state"] != "condition_CR":
            buckets["not_condition_CR"].append(gid)
            continue

        if r["geometry_name"] == "baseline_corrected":
            s = _vlm_strip_series(r)
            if bcorr is None:
                bcorr = s
                buckets["reference_baseline_other_lineage"].append(gid)
            else:
                # NOT ASSUMED IDENTICAL: the .vsp3 checksum and every strip row are checked
                # before the second delivery is discarded as a re-delivery.
                assert s["sha256"] == bcorr["sha256"], "%s: same name, different .vsp3" % gid
                assert np.array_equal(s["cl"], bcorr["cl"]) and \
                       np.array_equal(s["chord_ft"], bcorr["chord_ft"]) and \
                       np.array_equal(s["Yavg_ft"], bcorr["Yavg_ft"]), \
                       "%s: same .vsp3, different strips" % gid
                buckets["duplicate_redelivery_byte_identical"].append(gid)
            continue

        prefixes = r["mapping"].get("our_case_prefixes") or []
        letters = sorted({p[2] for p in prefixes if p.startswith("SW")})
        if not letters:
            buckets["unmapped"].append(gid)
            continue
        L = letters[0]
        assert L not in series, "two VLM runs claim family %s" % L

        rq = r["reference_quantities_AS_READ_FROM_DECK"]
        series[L] = _vlm_strip_series(r)
        buckets["plotted"].append(gid)

    n = len(vlm["runs"])
    got = sum(len(v) for v in buckets.values())
    assert got == n, "VLM partition does not close: %d bucketed != %d runs" % (got, n)

    # D068 again, on the differencing itself: nothing is subtracted until the two frames
    # are shown equal. The reference baseline is included in the check, because it is one
    # of the two sides of every VLM difference drawn.
    allv = list(series.values()) + ([bcorr] if bcorr else [])
    crefs = {round(v["Cref_ft"], 6) for v in allv}
    brefs = {round(v["bref_ft"], 6) for v in allv}
    srefs = {round(v["Sref_ft2"], 6) for v in allv}
    assert len(crefs) == 1 and len(brefs) == 1 and len(srefs) == 1, (
        "VLM runs do not share a frame: Cref %s, bref %s, Sref %s" % (crefs, brefs, srefs))
    return series, bcorr, buckets


def read_dense_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for k in rows[0]:
        out[k] = np.array([float(r[k]) for r in rows])
    return out


def read_section_csv(path):
    """One committed 5-station Cp cut: header dict plus (x/c, z/c, cp)."""
    hdr, xs, zs, cps = {}, [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                parts = line[1:].strip().split(",", 1)
                if len(parts) == 2:
                    hdr[parts[0].strip()] = parts[1].strip().strip('"')
                continue
            if line.startswith("x_over_c"):
                continue
            f = line.split(",")
            xs.append(float(f[0])); zs.append(float(f[1])); cps.append(float(f[2]))
    return hdr, np.array(xs), np.array(zs), np.array(cps)


def rans_spanwise(case, semi_m):
    """Sectional load for one case. DENSE route preferred, 5-station route as fallback.

    Which route ran is RETURNED, not assumed, and the caller reports it. The fallback is a
    genuinely coarser measurement, not a degraded version of the same one: five stations
    spanning eta 0.20 to 0.95 cannot see the root, so its C_L closure is far from unity and
    that is a property of the sampling, not a defect in the solution.

    ETA IS REBUILT HERE ON ONE SEMISPAN, `semi_m`, AND THAT IS A D068 FIX, NOT A TIDY-UP.
    The dense extractor placed its bands on the SAMPLED SURFACE's own maximum y, which is a
    face CENTROID and sits about 3.4 mm inboard of the reference tip; the committed
    5-station cuts used the registered bref/2. Left alone the two RANS routes would carry
    eta axes 0.19% apart from each other and from the VLM's Yavg/(bref/2), and a 0.19%
    stretch is a real shift against a load gradient of order 0.5 per unit eta. y in metres
    is the frame-free quantity both routes record, so eta is rebuilt from it.
    """
    dense = os.path.join(SPANWISE_DIR, "%s_spanwise_dense.csv" % case)
    prov = os.path.join(SPANWISE_DIR, "%s_spanwise_dense.json" % case)
    if os.path.exists(dense) and os.path.exists(prov):
        d = read_dense_csv(dense)
        p = load_json(prov)
        eta = d["y_m"] / semi_m
        return dict(route="dense_surface_integration", eta=eta,
                    semispan_from_surface_m=p.get("semispan_from_surface_m"),
                    semispan_used_m=semi_m,
                    load=d["load_cl_c_over_cref"], cl=d["cl"], cn=d["cn"], ca=d["ca"],
                    chord_m=d["chord_m"], n_stations=len(eta),
                    Cref_m=p["Cref_m_DSO"], Sref_m2=p["Sref_m2_DSO"],
                    alpha_deg=p["alpha_deg"], q=p["q_from_file"],
                    CL_from_strips=p.get("CL_from_strips"),
                    CL_forceCoeffs=p.get("CL_forceCoeffs"),
                    CL_closure_rel=p.get("CL_closure_rel"),
                    eta_min=float(eta.min()), eta_max=float(eta.max()),
                    vtk=p.get("vtk"), vtk_sha256=p.get("vtk_sha256"),
                    time=p.get("time"), skipped=p.get("skipped", []))

    files = sorted(glob.glob(os.path.join(SECTIONS_DIR, "%s_eta*.csv" % case)))
    if not files:
        return None

    # The quadrature is the one already committed and closure-tested; importing it keeps a
    # single implementation (D056) and makes the two routes differ ONLY in their input path.
    sys.path.insert(0, os.path.normpath(os.path.join(FIG_DIR, "..", "..", "..", "..", "scripts")))
    from spanwise_dense import integrate_station  # noqa: E402

    eta, load, cl_l, ch_l, cn_l, ca_l, alphas = [], [], [], [], [], [], []
    Cref_m = None
    for f in files:
        hdr, xc, zc, cp = read_section_csv(f)
        dd = [float(v) for v in hdr["dragDir"].split()]
        alpha = math.degrees(math.atan2(dd[2], dd[0]))
        c = float(hdr["local_chord_m"])
        r = integrate_station(xc, cp, zc, alpha)
        if r is None:
            continue
        cn, ca, cl, _nu, _nl = r
        eta.append(float(hdr["eta"])); cl_l.append(cl); ch_l.append(c)
        cn_l.append(cn); ca_l.append(ca); alphas.append(alpha)
    if not eta:
        return None
    return dict(route="five_station_section_cuts", eta=np.array(eta),
                semispan_from_surface_m=None, semispan_used_m=semi_m,
                cl=np.array(cl_l), cn=np.array(cn_l), ca=np.array(ca_l),
                chord_m=np.array(ch_l), n_stations=len(eta),
                Cref_m=None, Sref_m2=None, alpha_deg=alphas[0], q=None,
                CL_from_strips=None, CL_forceCoeffs=None, CL_closure_rel=None,
                eta_min=min(eta), eta_max=max(eta), vtk=None, vtk_sha256=None,
                time=None, skipped=[], load=None)


def rans_series(rans, condition, cref_m, semi_m):
    """RANS sectional loads for every converged TRIM case at one condition.

    cref_m is passed in from the VLM decks where a VLM comparison exists, so the two sides'
    normalising chord is ONE number used twice rather than two numbers assumed equal. Where
    there is no VLM side it is the registered DSO Cref read from the RANS provenance.
    """
    series, buckets = {}, {"plotted": [], "not_this_condition": [], "not_trim": [],
                           "not_converged": [], "no_spanwise_data": []}
    for name, c in sorted(rans["cases"].items()):
        if c.get("condition") != condition:
            buckets["not_this_condition"].append(name)
            continue
        if c.get("status") != "converged":
            buckets["not_converged"].append(name)
            continue
        if not c.get("is_trim_case"):
            buckets["not_trim"].append(name)
            continue
        sp = rans_spanwise(name, semi_m)
        if sp is None:
            buckets["no_spanwise_data"].append(name)
            continue
        L = c["family_letter"]
        if sp["Cref_m"] is not None:
            assert abs(sp["Cref_m"] - cref_m) / cref_m < 1e-4, (
                "%s normalises on Cref %.6f m, the comparison frame is %.6f m"
                % (name, sp["Cref_m"], cref_m))
            load = sp["load"]
        else:
            load = sp["cl"] * sp["chord_m"] / cref_m
        series[L] = dict(case=name, label=c.get("family_label"), eta=sp["eta"], load=load,
                         cl=sp["cl"], chord_m=sp["chord_m"], route=sp["route"],
                         n_stations=sp["n_stations"], alpha_deg=c["alpha_deg"],
                         Cl_forces=c["Cl"], Cd_ct=c["Cd"] * 1.0e4,
                         turb=c["turbulence_model"], Aref_m2=c["Aref_m2"],
                         lRef_m=c["lRef_m"], U=c["magUInf_m_s"],
                         CL_from_strips=sp["CL_from_strips"],
                         CL_closure_rel=sp["CL_closure_rel"],
                         alpha_from_surface=sp["alpha_deg"], vtk=sp["vtk"],
                         vtk_sha256=sp["vtk_sha256"], time=sp["time"],
                         skipped=sp["skipped"], Cref_m=cref_m,
                         semispan_from_surface_m=sp["semispan_from_surface_m"],
                         semispan_used_m=sp["semispan_used_m"])
        buckets["plotted"].append(name)

    n = len(rans["cases"])
    got = sum(len(v) for v in buckets.values())
    assert got == n, "RANS partition does not close: %d bucketed != %d cases" % (got, n)

    # alpha is recorded twice by two different readers (the force file at harvest time and
    # the force file at extraction time); if they disagree the two are not the same solution.
    for L, s in series.items():
        if s["alpha_from_surface"] is not None:
            d = abs(s["alpha_from_surface"] - s["alpha_deg"])
            assert d < 1e-6, ("%s: alpha %.9f from the spanwise extraction against %.9f "
                              "from rans_forces.json" % (s["case"], s["alpha_from_surface"],
                                                         s["alpha_deg"]))
    return series, buckets


# ==================================================================================
# differencing
# ==================================================================================
def trim_note(rans, condition, letters):
    """One line naming each plotted geometry's trim residual, read from rans_forces.json.

    THE FIGURE IS AT A FIXED C_L, so the residual against that target is the quantity that
    says whether "at the same C_L" is literally true. It is stated rather than absorbed
    (D051: establish what a number is load-bearing for before leaving it out). The
    tolerance is the project's own, also read from the file.
    """
    g = rans.get("trims", {}).get(condition, {}).get("geometries", {})
    parts, worst, worst_L = [], 0.0, None
    for L in ORDER:
        if L not in letters or L not in g:
            continue
        v = g[L].get("Cl_minus_target_counts")
        if not isinstance(v, (int, float)):
            continue
        parts.append("%s %+.2f" % (L, v))
        if abs(v) > abs(worst):
            worst, worst_L = v, L
    if not parts:
        return ""
    tail = ""
    if abs(worst) > 1.0:
        tail = ("  %s is %.2fx the 1-count trim tolerance and is the baseline for every "
                "difference in this figure." % (worst_L, abs(worst)))
    return ("Trim residual $C_L$ - target, counts (1 ct = 1e-4): "
            + ", ".join(parts) + "." + tail)


def difference_from_baseline(series, baseline=None, skip=("B",)):
    """Candidate minus baseline on the CANDIDATE's own eta grid.

    `baseline` defaults to series["B"]; pass a different one where the lineage-consistent
    reference is not the series' own B (the VLM case, see vlm_condition_cr).

    The baseline is INTERPOLATED, because morphing moves the VLM panel centroids and can
    move a RANS band's chord extent. The maximum eta offset and the slope-implied bound on
    the interpolation error are returned so the caption can state them rather than the
    figure implying an exactness it does not have.
    """
    b = baseline if baseline is not None else series.get("B")
    if b is None:
        return {}, None
    order = np.argsort(b["eta"])
    be, bl = b["eta"][order], b["load"][order]
    out, max_off, max_bound = {}, 0.0, 0.0
    slope = float(np.abs(np.gradient(bl, be)).max())
    for L, s in series.items():
        if L in skip:
            continue
        e = s["eta"]
        out[L] = dict(eta=e, dload=s["load"] - np.interp(e, be, bl))
        if len(e) == len(be):
            off = float(np.abs(np.sort(e) - be).max())
            max_off = max(max_off, off)
            max_bound = max(max_bound, slope * off)
    return out, dict(max_eta_offset=max_off, max_interp_error_bound=max_bound,
                     baseline=b["geometry"] if "geometry" in b else b.get("case"))


def truncation_null(vseries, eta_lo, eta_hi):
    """How much of the RANS closure gap is the RANS eta window, measured not asserted.

    The RANS strip integral spans only the stations that were cut, so its C_L closure is
    short by whatever the unsampled root and tip carry. That deficit is NOT separately
    knowable from the RANS data alone. It IS knowable from the VLM strips, whose full-span
    integral reproduces the .polar C_L to 5e-06: integrating the SAME .lod over the SAME
    eta window and comparing against the full-span value measures the truncation directly,
    on a case whose answer is independently known (D058: a null needs an answer you know).

    This does not correct the RANS number. It bounds how much of the gap is geometry-free
    truncation and how much is left for the pressure-only quadrature to explain.
    """
    out = {}
    for L, v in vseries.items():
        m = (v["eta"] >= eta_lo) & (v["eta"] <= eta_hi)
        full = float(2.0 * np.sum(v["cl"] * v["dArea_ft2"]) / v["Sref_ft2"])
        trunc = float(2.0 * np.sum(v["cl"][m] * v["dArea_ft2"][m]) / v["Sref_ft2"])
        out[L] = dict(full=full, truncated=trunc, rel=(trunc - full) / full,
                      n_strips_kept=int(m.sum()), n_strips_total=int(len(m)))
    return out


# ==================================================================================
# figures
# ==================================================================================
def panel(ax, series, keys, ykey, xkey="eta", markevery=None, lw=1.2):
    for L in ORDER:
        if L not in series or L not in keys:
            continue
        s = series[L]
        st = STYLE[L]
        ax.plot(s[xkey], s[ykey], color=st["color"], ls=st["ls"], lw=lw,
                marker=(st["marker"] if markevery else None),
                markevery=markevery, ms=2.8, mew=0.6, mfc="none", label=L)


def fig_condition_cr(rs, vs, bcorr, rd, vd, lineage, rmeta, vmeta, out_pdf):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), sharex=True)
    (a1, a2), (a3, a4) = axes

    for ax in (a1, a2, a3, a4):
        ax.set_xlim(0.0, 1.0)
        ax.grid(True, lw=0.35, color="0.88")
        ax.set_axisbelow(True)

    # --- row 1: absolute load -----------------------------------------------------
    panel(a1, rs, set(rs), "load")
    panel(a2, vs, set(vs), "load", markevery=6)
    # The second VLM baseline deck is not drawn: one baseline, B (22 Sep). Its series is
    # still loaded and asserted above, so the frame checks it carries still run.
    a1.set_title("(a) RANS", loc="left")
    a2.set_title("(b) VLM", loc="left")
    a1.set_ylabel(r"$c_\ell\,c\,/\,C_{\mathrm{ref}}$")

    y1 = np.concatenate([rs[L]["load"] for L in rs] + [vs[L]["load"] for L in vs]
                        + ([bcorr["load"]] if bcorr is not None else []))
    lo, hi = float(y1.min()), float(y1.max())
    pad = 0.06 * (hi - lo)
    for ax in (a1, a2):
        ax.set_ylim(lo - pad, hi + pad)

    # --- row 2: difference from baseline ------------------------------------------
    for L in ORDER:
        if L in rd:
            st = STYLE[L]
            a3.plot(rd[L]["eta"], rd[L]["dload"], color=st["color"], ls=st["ls"], lw=1.25)
        if L in vd:
            st = STYLE[L]
            a4.plot(vd[L]["eta"], vd[L]["dload"], color=st["color"], ls=st["ls"], lw=1.25,
                    marker=st["marker"], markevery=6, ms=2.8, mew=0.6, mfc="none")
    # lineage offset curve not drawn: one baseline, B (22 Sep)
    a3.set_title("(c) RANS $\\Delta$ from B", loc="left")
    a4.set_title("(d) VLM $\\Delta$ from B", loc="left")
    a3.set_ylabel(r"$\Delta\,(c_\ell\,c\,/\,C_{\mathrm{ref}})$")
    for ax in (a3, a4):
        ax.axhline(0.0, color="0.5", lw=0.7)
        ax.set_xlabel(r"$\eta = y/(b/2)$")

    # Row-2 limits come from the CANDIDATE differences only. The grey lineage curve is
    # allowed to run off the bottom of (d) rather than compressing the signal both panels
    # exist to show; its extreme value is annotated on the axis so nothing is hidden.
    y2 = np.concatenate([rd[L]["dload"] for L in rd] + [vd[L]["dload"] for L in vd])
    lo2, hi2 = float(y2.min()), float(y2.max())
    pad2 = 0.16 * (hi2 - lo2)
    for ax in (a3, a4):
        ax.set_ylim(lo2 - pad2, hi2 + pad2)
    if lineage is not None:
        k = int(np.argmax(np.abs(lineage["dload"])))
        pass  # second-baseline annotation removed: one baseline, B (22 Sep)

    handles = [Line2D([], [], color=STYLE[L]["color"], ls=STYLE[L]["ls"], lw=1.3,
                      marker=STYLE[L]["marker"], ms=3.0, mew=0.6, mfc="none",
                      label="%s  %s" % (L, rmeta["labels"].get(L, "")))
               for L in ORDER if L in rs or L in vs]
    fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.015),
               columnspacing=1.4, handlelength=2.8)
    if vmeta.get("trim_note"):
        pass  # trim note lives in the prose (22 Sep)

    fig.tight_layout(rect=(0, 0.0, 1, 1.0))
    # The header is placed AFTER tight_layout, as figure text just above the axes, because
    # suptitle reserves a band whose height tight_layout then cannot reclaim.
    fig.text(0.5, 1.005,
             "",
             ha="center", va="bottom", fontsize=9.0)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def fig_cruise(rs, rd, meta, out_pdf):
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9))
    a1, a2 = axes
    for ax in (a1, a2):
        ax.set_xlim(0.0, 1.0)
        ax.grid(True, lw=0.35, color="0.85")
        ax.set_axisbelow(True)
        ax.set_xlabel(r"$\eta = y/(b/2)$")

    panel(a1, rs, set(rs), "load")
    a1.set_title("(a)", loc="left")
    a1.set_ylabel(r"$c_\ell\,c\,/\,C_{\mathrm{ref}}$")

    for L in ORDER:
        if L in rd:
            st = STYLE[L]
            a2.plot(rd[L]["eta"], rd[L]["dload"], color=st["color"], ls=st["ls"], lw=1.2)
    a2.axhline(0.0, color=GREY, lw=0.7)
    a2.set_title(r"(b) $\Delta$ from B", loc="left")
    a2.set_ylabel(r"$\Delta\,(c_\ell\,c\,/\,C_{\mathrm{ref}})$")

    if meta.get("missing"):
        print("ABSENT FROM THE CRUISE SPANWISE FIGURE, NAME IT IN THE CAPTION: %s"
              % ", ".join(meta["missing"]))
    if meta.get("trim_note"):
        pass  # trim note lives in the prose (22 Sep)

    handles = [Line2D([], [], color=STYLE[L]["color"], ls=STYLE[L]["ls"], lw=1.3,
                      label="%s  %s" % (L, meta["labels"].get(L, "")))
               for L in ORDER if L in rs]
    fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, -0.03),
               columnspacing=1.1, handlelength=2.8)
    fig.tight_layout(rect=(0, 0.0, 1, 1.0))
    fig.text(0.5, 1.01,
             "",
             ha="center", va="bottom", fontsize=9.0)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


# ==================================================================================
def main():
    vlm = load_json(VLM_JSON)
    rans = load_json(RANS_JSON)

    report = {}

    # ---------------- Condition CR ------------------------------------------------
    vs, bcorr, vbuckets = vlm_condition_cr(vlm)
    print("VLM runs bucketed:", {k: len(v) for k, v in vbuckets.items()})
    for k in ("unmapped", "reference_baseline_other_lineage",
              "duplicate_redelivery_byte_identical"):
        for g in vbuckets[k]:
            print("    %-38s %s" % (k, g))

    Cref_ft = next(iter(vs.values()))["Cref_ft"]
    Sref_ft2 = next(iter(vs.values()))["Sref_ft2"]
    # The RANS side normalises in metres. ONE conversion, stated, from the deck's own foot
    # value: 1 ft = 0.3048 m exactly (project rules, hard fact 1, a definition not a measurement).
    Cref_m = Cref_ft * 0.3048
    # ONE semispan for both sides of the comparison, from the VLM deck's own Bref.
    semi_m = 0.5 * next(iter(vs.values()))["bref_ft"] * 0.3048

    rs_cr, rbuckets = rans_series(rans, "condition_CR", Cref_m, semi_m)
    print("RANS condition_CR bucketed:", {k: len(v) for k, v in rbuckets.items()})
    for n in rbuckets["no_spanwise_data"]:
        print("    no_spanwise_data   %s" % n)

    print()
    print("=== .lod HALF-WING FACTOR OF 2 CLOSURE, recomputed here from the strip table ===")
    worst_x2, worst_x1 = 0.0, []
    for L, v in [(L, vs[L]) for L in ORDER if L in vs] + \
                ([("B*", bcorr)] if bcorr is not None else []):
        print("  %-2s %-28s  CL x2 %.8f  polar CLtot %.8f  rel %+.3e | without x2 ratio %.7f"
              % (L, v["geometry"], v["CL_with_factor_2"], v["CLtot_polar"],
                 v["rel_err_with_factor_2"], v["ratio_without_factor_2"]))
        worst_x2 = max(worst_x2, abs(v["rel_err_with_factor_2"]))
        worst_x1.append(v["ratio_without_factor_2"])
    print("  MAX |rel err| WITH factor 2 = %.3e ; without it the ratio spans %.7f to %.7f"
          % (worst_x2, min(worst_x1), max(worst_x1)))

    print()
    print("=== RANS spanwise route and its own CL closure (2*int c_l c dy / Sref) ===")
    for L in ORDER:
        if L not in rs_cr:
            continue
        s = rs_cr[L]
        cl_s = s["CL_from_strips"]
        rel = s["CL_closure_rel"]
        print("  %s %-10s route %-28s stations %3d  alpha %.6f  CL_strips %s  vs forces "
              "%.7f  gap %s"
              % (L, s["case"], s["route"], s["n_stations"], s["alpha_deg"],
                 ("%.6f" % cl_s) if cl_s is not None else "NOT COMPUTED",
                 s["Cl_forces"], ("%+.2f%%" % (100 * rel)) if rel is not None else "n/a"))

    rd_cr, rioff = difference_from_baseline(rs_cr)
    # VLM candidates are differenced against the LINEAGE-CONSISTENT baseline; B itself is
    # differenced too, and that curve IS the lineage offset drawn in grey.
    vd_all, vioff = difference_from_baseline(vs, baseline=bcorr, skip=())
    lineage = vd_all.pop("B", None)
    vd_cr = vd_all

    # THE TRUNCATION NULL: how much of the RANS closure gap is the RANS eta window alone,
    # measured on the VLM strips whose full-span answer is independently known.
    eta_lo = min(rs_cr[L]["eta"].min() for L in rs_cr)
    eta_hi = max(rs_cr[L]["eta"].max() for L in rs_cr)
    trunc = truncation_null(vs, eta_lo, eta_hi)
    print()
    print("=== TRUNCATION NULL: the same .lod integrated over the RANS window "
          "eta %.4f to %.4f ===" % (eta_lo, eta_hi))
    for L in ORDER:
        if L in trunc:
            t = trunc[L]
            print("  %s  full-span CL %.6f, window-only %.6f -> %+.3f%%  (%d of %d strips)"
                  % (L, t["full"], t["truncated"], 100 * t["rel"],
                     t["n_strips_kept"], t["n_strips_total"]))
    tr = [trunc[L]["rel"] for L in trunc]
    print("  VLM truncation costs %+.3f%% to %+.3f%%. The RANS gaps below are to be read "
          "against that, not against zero." % (100 * min(tr), 100 * max(tr)))

    # The legend names each family by the DSO GEOMETRY it is, taken from the VLM mapping
    # (which is by .vsp3 sha256, not by name). rans_forces.json's family_label reads
    # "morphing candidate F" for four of the six, which names nothing.
    labels = {L: vs[L]["geometry"] for L in vs}
    for L in rs_cr:
        labels.setdefault(L, rs_cr[L]["label"])
    rmeta = dict(labels=labels)
    turbs = sorted({rs_cr[L]["turb"] for L in rs_cr})
    assert len(turbs) == 1, "Condition CR is not on one turbulence model: %s" % turbs
    vmeta = dict(target_CL=next(iter(vs.values()))["CL_target"],
                 Sref_ft2=Sref_ft2, Cref_ft=Cref_ft,
                 bref_ft=next(iter(vs.values()))["bref_ft"], turb=turbs[0],
                 trim_note=trim_note(rans, "condition_CR", set(rs_cr)))

    out1 = os.path.join(FIG_DIR, "spanwise_loading_condition_CR.pdf")
    fig_condition_cr(rs_cr, vs, bcorr, rd_cr, vd_cr, lineage, rmeta, vmeta, out1)

    print()
    print("=== PEAK DIFFERENCE FROM BASELINE, Condition CR ===")
    for L in ORDER:
        if L == "B":
            continue
        line = "  %s" % L
        for tag, dd in (("RANS", rd_cr), ("VLM ", vd_cr)):
            if L in dd:
                d = dd[L]["dload"]
                k = int(np.argmax(np.abs(d)))
                line += "   %s peak %+.4f at eta %.3f" % (tag, d[k], dd[L]["eta"][k])
            else:
                line += "   %s NOT FOUND" % tag
        print(line)
    if lineage is not None:
        k = int(np.argmax(np.abs(lineage["dload"])))
        print("  LINEAGE OFFSET  B(as meshed) - B(corrected), VLM: peak %+.4f at eta %.3f, "
              "rms %.4f" % (lineage["dload"][k], lineage["eta"][k],
                            float(np.sqrt((lineage["dload"] ** 2).mean()))))

    # HOW CLOSELY THE TWO METHODS AGREE ON THE MORPHING EFFECT ITSELF. The VLM difference
    # is interpolated onto the RANS eta grid over the OVERLAP only; the two are not compared
    # outside the range both sampled, because an extrapolated VLM strip is not a measurement.
    agree = {}
    print("  RANS minus VLM, on the DIFFERENCE curves (overlap only):")
    for L in ORDER:
        if L not in rd_cr or L not in vd_cr:
            continue
        e = rd_cr[L]["eta"]
        ve, vl = vd_cr[L]["eta"], vd_cr[L]["dload"]
        o = np.argsort(ve)
        m = (e >= ve[o].min()) & (e <= ve[o].max())
        g = rd_cr[L]["dload"][m] - np.interp(e[m], ve[o], vl[o])
        agree[L] = dict(rms=float(np.sqrt((g ** 2).mean())), max_abs=float(np.abs(g).max()),
                        eta_lo=float(e[m].min()), eta_hi=float(e[m].max()),
                        n_points=int(m.sum()))
        print("    %s  rms %.4f, max %.4f over eta %.3f to %.3f (%d stations)"
              % (L, agree[L]["rms"], agree[L]["max_abs"], agree[L]["eta_lo"],
                 agree[L]["eta_hi"], agree[L]["n_points"]))
    print("  baseline interpolation: max eta offset RANS %.2e, VLM %.2e; "
          "implied bound on the difference RANS %.2e, VLM %.2e"
          % (rioff["max_eta_offset"] if rioff else float("nan"),
             vioff["max_eta_offset"] if vioff else float("nan"),
             rioff["max_interp_error_bound"] if rioff else float("nan"),
             vioff["max_interp_error_bound"] if vioff else float("nan")))

    report["condition_CR"] = dict(
        vlm_buckets={k: v for k, v in vbuckets.items()},
        rans_routes={L: rs_cr[L]["route"] for L in rs_cr},
        rans_cases={L: rs_cr[L]["case"] for L in rs_cr},
        rans_stations={L: rs_cr[L]["n_stations"] for L in rs_cr},
        rans_turbulence=turbs,
        vlm_geometry={L: vs[L]["geometry"] for L in vs},
        vlm_reference_baseline=(bcorr["geometry"] if bcorr else None),
        vlm_n_strips={L: vs[L]["n_strips"] for L in vs},
        lod_closure={L: dict(with_factor_2=v["CL_with_factor_2"],
                             polar_CLtot=v["CLtot_polar"],
                             rel_err=v["rel_err_with_factor_2"],
                             ratio_without_factor_2=v["ratio_without_factor_2"])
                     for L, v in list(vs.items()) + ([("B_corrected", bcorr)] if bcorr else [])},
        rans_closure={L: dict(CL_from_strips=rs_cr[L]["CL_from_strips"],
                              CL_forceCoeffs=rs_cr[L]["Cl_forces"],
                              rel=rs_cr[L]["CL_closure_rel"]) for L in rs_cr},
        rans_eta_window=[eta_lo, eta_hi],
        semispan_reconciliation={
            L: dict(from_surface_m=rs_cr[L]["semispan_from_surface_m"],
                    used_m=rs_cr[L]["semispan_used_m"],
                    difference_mm=(1000.0 * (rs_cr[L]["semispan_used_m"]
                                             - rs_cr[L]["semispan_from_surface_m"]))
                    if rs_cr[L]["semispan_from_surface_m"] else None)
            for L in rs_cr},
        lod_closure_max_abs_rel_with_factor_2=worst_x2,
        lod_ratio_without_factor_2_range=[min(worst_x1), max(worst_x1)],
        vlm_truncation_null=trunc,
        peak_difference={L: dict(
            rans=(float(rd_cr[L]["dload"][int(np.argmax(np.abs(rd_cr[L]["dload"])))]),
                  float(rd_cr[L]["eta"][int(np.argmax(np.abs(rd_cr[L]["dload"])))]))
            if L in rd_cr else None,
            vlm=(float(vd_cr[L]["dload"][int(np.argmax(np.abs(vd_cr[L]["dload"])))]),
                 float(vd_cr[L]["eta"][int(np.argmax(np.abs(vd_cr[L]["dload"])))]))
            if L in vd_cr else None) for L in ORDER if L != "B"},
        lineage_offset=(dict(
            peak=float(lineage["dload"][int(np.argmax(np.abs(lineage["dload"])))]),
            eta_at_peak=float(lineage["eta"][int(np.argmax(np.abs(lineage["dload"])))]),
            rms=float(np.sqrt((lineage["dload"] ** 2).mean())),
            what="VLM B(baseline_wing_only_refined, the surface every RANS B case meshes) "
                 "MINUS B(baseline_corrected, the lineage all five candidates are on). "
                 "The RANS difference panel carries this term and cannot remove it, "
                 "because no RANS case has ever been meshed on the corrected baseline "
                 "(registry: NOT_USED_overruled_by_ARG-099, 'DO NOT MESH').")
            if lineage is not None else None),
        rans_minus_vlm_on_the_difference=agree,
        interp=dict(rans=rioff, vlm=vioff),
        figure=out1)

    # ---------------- cruise, both conditions ------------------------------------
    # Late cruise was excluded by design while it had no trims; both cruise conditions
    # now carry welded trims with dense loads (ARG-196), so both are drawn the same way.
    outs = [out1]
    for cond in ("early_cruise", "late_cruise"):
        rs_cc, cbuckets = rans_series(rans, cond, Cref_m, semi_m)
        print()
        print("RANS %s bucketed:" % cond, {k: len(v) for k, v in cbuckets.items()})
        for n in cbuckets["no_spanwise_data"]:
            print("    no_spanwise_data   %s" % n)
        if rs_cc and "B" in rs_cc:
            print("=== RANS spanwise, %s ===" % cond)
            for L in ORDER:
                if L not in rs_cc:
                    continue
                s = rs_cc[L]
                print("  %s %-10s route %-28s stations %3d  alpha %.6f  CL_strips %s  vs forces "
                      "%.7f  gap %s"
                      % (L, s["case"], s["route"], s["n_stations"], s["alpha_deg"],
                         ("%.6f" % s["CL_from_strips"]) if s["CL_from_strips"] is not None
                         else "NOT COMPUTED", s["Cl_forces"],
                         ("%+.2f%%" % (100 * s["CL_closure_rel"]))
                         if s["CL_closure_rel"] is not None else "n/a"))
            rd_cc, cioff = difference_from_baseline(rs_cc)
            # target C_L for this condition, READ from a file, not typed
            target_cl = vlm["_operating_points"][cond]["target_CL"]
            cturbs = sorted({rs_cc[L]["turb"] for L in rs_cc})
            assert len(cturbs) == 1, "%s is not on one turbulence model: %s" % (cond, cturbs)
            # The family absent from this condition, named rather than omitted (GEO-089). It is
            # DERIVED from the letters present at Condition CR, not typed.
            missing = [("%s (%s)" % (L, labels.get(L, "")))
                       for L in ORDER if L in rs_cr and L not in rs_cc]
            meta = dict(labels=labels, missing=missing,
                        target_CL=target_cl, Cref_m=Cref_m, turb=cturbs[0],
                        trim_note=trim_note(rans, cond, set(rs_cc)))
            out2 = os.path.join(FIG_DIR, "spanwise_loading_%s.pdf" % cond)
            fig_cruise(rs_cc, rd_cc, meta, out2)
            outs.append(out2)
            print("=== PEAK DIFFERENCE FROM BASELINE, %s ===" % cond)
            for L in ORDER:
                if L in rd_cc:
                    d = rd_cc[L]["dload"]
                    k = int(np.argmax(np.abs(d)))
                    print("  %s  RANS peak %+.4f at eta %.3f" % (L, d[k], rd_cc[L]["eta"][k]))
            report[cond] = dict(
                rans_routes={L: rs_cc[L]["route"] for L in rs_cc},
                rans_cases={L: rs_cc[L]["case"] for L in rs_cc},
                rans_stations={L: rs_cc[L]["n_stations"] for L in rs_cc},
                rans_turbulence=sorted({rs_cc[L]["turb"] for L in rs_cc}),
                rans_closure={L: dict(CL_from_strips=rs_cc[L]["CL_from_strips"],
                                      CL_forceCoeffs=rs_cc[L]["Cl_forces"],
                                      rel=rs_cc[L]["CL_closure_rel"]) for L in rs_cc},
                peak_difference={L: (float(rd_cc[L]["dload"][int(np.argmax(np.abs(rd_cc[L]["dload"])))]),
                                     float(rd_cc[L]["eta"][int(np.argmax(np.abs(rd_cc[L]["dload"])))]))
                                 for L in rd_cc},
                missing_families=missing,
                missing_case_reason={n: rans["cases"][n].get("status")
                                     for n in cbuckets["not_converged"]
                                     if rans["cases"][n].get("is_trim_case")},
                vlm_counterpart="NONE. Every VSPAERO cruise run in vlm.json is on the TP-1580 "
                                "Sref basis AND is a different geometry set (all six unmapped "
                                "by .vsp3 sha256). No frame-equal VLM comparison exists.",
                figure=out2)
        else:
            print("%s: no baseline spanwise data, figure NOT produced" % cond)
            report[cond] = dict(figure=None,
                                          reason="no baseline spanwise data available")

    # ---------------- verify the files exist, on disk, non-zero -------------------
    print()
    print("=== FILES WRITTEN (checked on disk, not inferred from the exit code) ===")
    for p in outs:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            print("  OK   %s  %d bytes" % (p, os.path.getsize(p)))
        else:
            print("  FAIL %s  MISSING OR EMPTY" % p)

    jf = os.path.join(FIG_DIR, "spanwise_caption_values.json")
    with open(jf, "w") as fh:
        fh.write(json.dumps(report, indent=2, default=float) + "\n")
    print("  OK   %s  %d bytes" % (jf, os.path.getsize(jf)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
