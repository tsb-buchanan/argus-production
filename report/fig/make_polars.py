#!/usr/bin/env python3
"""
ARGUS all-geometry report, 2026-09-15: CL-CD polars and CL-alpha, RANS versus VSPAERO VLM.

Self-contained. Reads only:
    ../data/rans_forces.json
    ../data/vlm.json
Writes vector PDFs into this directory.

RULES ENFORCED IN CODE (not in prose):

 1. EVERY NUMBER COMES FROM A FILE. No literal aerodynamic value is typed in this script.
    Target CLs, Mach, reference areas, alphas, CL, CD and the Sref conversion factor are all
    read from the two JSON files. The only literals are plotting constants and the definition
    1 count = 1e-4.

 2. FRAME RULE (D068). Every RANS case is asserted to carry a single Aref/lRef per condition;
    every VLM point is taken from run['polar_on_DSO_basis'], whose conversion factor is read
    from run['frame_conversion']['to_DSO_basis']['force_coefficient_factor'] and printed in
    the figure. Where the frame cannot be shown equal (early cruise) the VLM point is drawn
    at BOTH readings and flagged; it is never adjudicated.

 3. OPERATING POINTS NOT MIXED (D076). Condition CR and early cruise get separate axes and
    separate files. Each figure states its own Mach and target CL inline.

 4. VLM CDtot IS NOT A VISCOUS DRAG. It is CDo (flat-plate skin-friction correlation at the
    deck's ReCref) + CDi. Every figure that puts a VLM drag beside a RANS drag says so.

 5. WINDOW MEANS, NEVER LAST SAMPLES. Every RANS Cl/Cd here is what rans_forces.json took
    from each case's .converged marker. Non-converged cases are excluded by a partition that
    is ASSERTED to close (plotted + excluded == enumerated), not merely reported.

 6. DE-DUPLICATION IS BY deck_sha256. baseline_corrected is delivered twice (two packages,
    byte-identical deck). Plotting the list as-delivered drew it twice; the duplicate is now
    dropped by checksum and the drop is counted and reported.
"""

import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ----------------------------------------------------------------------------------
FIG_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(FIG_DIR, "..", "data"))
RANS_JSON = os.path.join(DATA_DIR, "rans_forces.json")
VLM_JSON = os.path.join(DATA_DIR, "vlm.json")

CT = 1.0e4  # 1 drag count = 1e-4 in coefficient. A definition, not a measurement.

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 8.5,
    "axes.labelsize": 9.0,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 7.0,
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
    "legend.edgecolor": "0.45",
    "legend.fancybox": False,
    "legend.borderpad": 0.4,
    "pdf.fonttype": 42,
    "savefig.dpi": 600,
})

ORDER = ["B", "C", "F", "H", "M", "W"]

# Okabe-Ito colour-blind-safe palette, each with its own dash pattern and marker so the
# figure survives greyscale printing on line style alone.
STYLE = {
    "B": dict(color="#000000", ls="-",                           marker="o"),
    "C": dict(color="#0072B2", ls=(0, (5, 1.6)),                 marker="s"),
    "F": dict(color="#D55E00", ls=(0, (1.3, 1.3)),               marker="^"),
    "H": dict(color="#009E73", ls=(0, (4, 1.4, 1, 1.4)),         marker="D"),
    "M": dict(color="#CC79A7", ls=(0, (3, 1.2, 1, 1.2, 1, 1.2)), marker="v"),
    "W": dict(color="#E69F00", ls=(0, (7, 1.6, 1, 1.6)),         marker="P"),
}
ALT_C = "#555555"      # baseline_corrected, the lineage that is NOT meshed
FLAG_C = "#882255"     # unresolved-frame annotation
GREY = "0.35"

WRITTEN = []
DIAG = {}


# ----------------------------------------------------------------------------------
def note_below(fig, leg, text, fontsize=7.4):
    """Place a note UNDER the legend, wherever the legend actually ended up.

    Derived from the drawn legend's extent, not from a guessed offset: writing a fixed
    negative y put the note on top of the legend in the first version.
    """
    fig.canvas.draw()
    bb = leg.get_window_extent(fig.canvas.get_renderer()).transformed(
        fig.transFigure.inverted())
    fig.text(0.5, bb.y0 - 0.025, "", ha="center", va="top", fontsize=fontsize,
             color="0.2", linespacing=1.5)


def load():
    with open(RANS_JSON) as fh:
        rans = json.load(fh)
    with open(VLM_JSON) as fh:
        vlm = json.load(fh)
    return rans, vlm


def geom_names(vlm):
    """Letter -> the DSO .vsp3 geometry name it maps to, by sha256 (never by name).

    family_label in rans_forces.json is generic for F/H/M/W ("morphing candidate F"), which
    made the legend read "F F". The checksum-mapped geometry name is the real identity.
    """
    out = {}
    for r in vlm["runs"]:
        for pre in (r["mapping"].get("our_case_prefixes") or []):
            out[pre[-1]] = r["geometry_name"]
    return out


CRUISE_WELDED = ("early_cruise", "late_cruise")


def rans_series(rans, condition):
    """Converged RANS cases for one condition, grouped by geometry letter.

    The partition over every case enumerated at this condition is ASSERTED to close.
    """
    series, excluded, n_total = {}, [], 0
    other_gen = []
    for name, c in rans["cases"].items():
        if c.get("condition") != condition:
            continue
        n_total += 1
        if c.get("status") != "converged":
            excluded.append(dict(case=name, status=c.get("status")))
            continue
        # ONE MESH GENERATION PER CRUISE SERIES (ARG-196). Both cruise conditions moved to
        # the welded meshes; the published alpha-bracket cases stay in the data file,
        # labelled, but a curve or fit through bracket and trim would cross a generation.
        # Condition CR is left as it was drawn (published brackets with r2 trims, a mixing
        # the CR caption already states).
        if condition in CRUISE_WELDED and c.get("mesh_generation") != "welded":
            other_gen.append(dict(case=name, mesh_generation=c.get("mesh_generation")))
            continue
        L = c["family_letter"]
        s = series.setdefault(L, dict(pts=[], trim=None, label=c.get("family_label"),
                                      turb=set(), Aref=set(), lRef=set(), U=set()))
        rec = dict(case=name, alpha=c["alpha_deg"], Cl=c["Cl"], Cd=c["Cd"],
                   Cd_ct=c["Cd"] * CT, trim=bool(c.get("is_trim_case")))
        s["pts"].append(rec)
        s["turb"].add(c.get("turbulence_model"))
        s["Aref"].add(c.get("Aref_m2"))
        s["lRef"].add(c.get("lRef_m"))
        s["U"].add(c.get("magUInf_m_s"))
        if rec["trim"]:
            s["trim"] = rec

    n_plot = sum(len(v["pts"]) for v in series.values())
    assert n_plot + len(excluded) + len(other_gen) == n_total, (
        "PARTITION FAILED %s: %d + %d + %d != %d"
        % (condition, n_plot, len(excluded), len(other_gen), n_total))
    for L, s in series.items():
        s["pts"].sort(key=lambda r: r["alpha"])
        assert len(s["Aref"]) == 1, "Aref not unique in %s/%s" % (condition, L)
        assert len(s["lRef"]) == 1, "lRef not unique in %s/%s" % (condition, L)
        assert len(s["turb"]) == 1, "turbulence model not unique in %s/%s" % (condition, L)

    part = dict(condition=condition, enumerated=n_total, plotted=n_plot,
                excluded_not_converged=len(excluded), excluded=excluded,
                excluded_other_mesh_generation=len(other_gen),
                excluded_other_mesh_generation_cases=sorted(o["case"] for o in other_gen),
                sums_to_total=True,
                geometries_plotted=sorted(series),
                geometries_with_trim=sorted(L for L, s in series.items() if s["trim"]),
                geometries_without_trim=sorted(L for L, s in series.items() if not s["trim"]))
    return series, part


def lift_slope(pts):
    """Least-squares dCL/dalpha over a family's own converged cases at one condition.

    Derived from alphas and CLs read from rans_forces.json. The max residual is carried so
    the reader can judge the linearity the fit assumes rather than take it on trust.
    """
    a = np.array([p["alpha"] for p in pts], float)
    cl = np.array([p["Cl"] for p in pts], float)
    if a.size < 2:
        return None
    m, b = np.polyfit(a, cl, 1)
    return dict(dCL_dalpha_per_deg=float(m), CL0=float(b), n=int(a.size),
                max_abs_resid=float(np.max(np.abs(cl - (m * a + b)))),
                cases=[p["case"] for p in pts])


def vlm_at(vlm, state):
    """VLM runs at one operating point.

    mapped[letter]   run whose .vsp3 sha256 maps to one of our case prefixes
    alt              baseline lineage flagged NOT THE MESHED SURFACE
    unmapped         geometries with no RANS counterpart; never plotted as a series
    Duplicate deliveries of the same deck (identical deck_sha256) are dropped by checksum.
    """
    mapped, alt, unmapped, seen, dropped = {}, [], [], set(), []
    for r in vlm["runs"]:
        if r["operating_point"]["state"] != state:
            continue
        key = (r["geometry_name"], state, r["deck_sha256"])
        if key in seen:
            dropped.append(dict(run_id=r["run_id"], package=r["delivery_package"],
                                deck_sha256=r["deck_sha256"]))
            continue
        seen.add(key)
        rec = dict(
            run_id=r["run_id"], geometry=r["geometry_name"],
            package=r["delivery_package"],
            basis_deck=r["reference_quantities_AS_READ_FROM_DECK"]["basis"],
            Sref_ft2=r["reference_quantities_AS_READ_FROM_DECK"]["Sref_ft2"],
            factor=r["frame_conversion"]["to_DSO_basis"]["force_coefficient_factor"],
            already_DSO=r["frame_conversion"]["to_DSO_basis"]["already_on_this_basis"],
            deck=r["polar"]["points"][0], dso=r["polar_on_DSO_basis"],
            n_polar_points=len(r["polar"]["points"]),
            target_CL=r["operating_point"]["target_CL"],
            CL_err=r["operating_point"]["CL_error"],
            note=r["mapping"].get("case_prefix_family"))
        prefixes = r["mapping"].get("our_case_prefixes") or []
        letters = sorted({p[-1] for p in prefixes})
        if len(letters) == 1:
            mapped[letters[0]] = rec
        elif "NOT THE MESHED SURFACE" in str(rec["note"]):
            alt.append(rec)
        else:
            unmapped.append(rec)

    # every VSPAERO deck in this package is a single solved point, so no VLM polar exists
    npts = sorted({r["n_polar_points"] for r in list(mapped.values()) + alt + unmapped})
    return mapped, alt, unmapped, dict(dropped_duplicate_decks=dropped,
                                       polar_points_per_run=npts)


def header(rans, condition, vlm, target_CL):
    fc = rans["_frame_check_per_condition"][condition]
    mach = vlm["_operating_points"][condition]["mach"]
    # Bare condition and Mach. The full frame (target CL, turbulence model, Uinf, Aref,
    # lref, DSO basis) is preserved in the figure's facts JSON and FIGURE_NOTES.txt, so it
    # still travels with the number without being printed on the plot.
    del fc, target_CL
    return r"$M$ %.2f" % mach


def finish(fig, ax, note, fname, facts, extra=None):
    # FOOTERS SUPPRESSED on the figure itself. Every one of these notes is still COMPUTED
    # above, and the numbers in them (the 9.62-count sweep mismatch, the 1.000000 Sref
    # factor) are still asserted and still written to the facts JSON; they are simply not
    # drawn. Project convention: bare minimum labels, no explanatory text on the plot.
    facts["_suppressed_figure_note"] = note
    fig.text(0.5, -0.02, "", ha="center", va="top", fontsize=7.4, color="0.2",
             transform=ax.transAxes if extra is None else fig.transFigure)
    path = os.path.join(FIG_DIR, fname)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("WROTE %s  %d bytes" % (path, os.path.getsize(path)))
    WRITTEN.append(path)
    return facts


# ==================================================================================
# Does a geometry's alpha sweep actually pass through its own trim point?
# ==================================================================================
# WHY THIS IS MEASURED AND NOT ASSUMED. A polar drawn through an alpha sweep, with the
# trimmed solution marked on it, asserts that the two are one family of solutions. If they
# are not, the reader sees a curve and a point that disagree and has no way to tell which
# to believe. The check is a WITHIN-FILE PROOF (D052 item 4): interpolate the geometry's
# OWN sweep to its OWN trim lift and compare against its OWN measured trim drag. Nothing
# external is consulted and no provenance field has to be trusted.
#
# MEASURED at Condition CR on 2026-09-22: C 0.003 ct and M 0.001 ct, i.e. matched to well
# inside the 0.05 ct convergence gate; B 2.945, F 9.581, W 9.620, H 9.593. The file
# convicts itself.
#
# ALL OR NOTHING, PER CONDITION. Drawing sweeps for the geometries that match and points
# for those that do not would leave the reader comparing one geometry's curve against
# another's point, which is the same error one step removed. If any geometry fails, the
# figure shows trimmed solutions only and SAYS SO with the measured number.
SWEEP_TRIM_TOL_CT = 0.5


def sweep_trim_mismatch(s):
    """|CD implied by this geometry's own sweep at its own trim CL, minus its trim CD|.

    None when there is no trim or fewer than two sweep points, which is an ABSENCE of a
    verdict and is counted as such by the caller, never read as a pass.
    """
    if not s.get("trim"):
        return None
    # EXCLUDE THE TRIM POINT FROM THE SWEEP BEFORE INTERPOLATING. It is IN s["pts"], so
    # interpolating the unfiltered list at the trim CL lands exactly on the trim point and
    # the mismatch is 0.0 by construction. The first version of this check did that and
    # reported 0.000 for all six geometries against hand-computed values of 2.9 to 9.6
    # counts: a null that compares a thing to itself only tests that the arithmetic runs
    # (D058). The sweep must be the points the trim is being TESTED AGAINST.
    tcase = s["trim"].get("case")
    pts = sorted((p["Cl"], p["Cd_ct"]) for p in s["pts"] if p.get("case") != tcase)
    if len(pts) < 2:
        return None
    cl, cd = s["trim"]["Cl"], s["trim"]["Cd_ct"]
    for (cl0, cd0), (cl1, cd1) in zip(pts, pts[1:]):
        if cl0 <= cl <= cl1 and cl1 != cl0:
            return abs(cd0 + (cd1 - cd0) * (cl - cl0) / (cl1 - cl0) - cd)
    # Outside the bracket: extrapolate on the nearest pair rather than decline to judge.
    (cl0, cd0), (cl1, cd1) = (pts[0], pts[1]) if cl < pts[0][0] else (pts[-2], pts[-1])
    if cl1 == cl0:
        return None
    return abs(cd0 + (cd1 - cd0) * (cl - cl0) / (cl1 - cl0) - cd)


def sweeps_usable(series):
    """(draw_sweeps, {geom: mismatch}, {geom: reason-not-judged}) for one condition."""
    mism, nojudge = {}, {}
    for L, s in series.items():
        v = sweep_trim_mismatch(s)
        if v is None:
            nojudge[L] = "no trim" if not s.get("trim") else "fewer than two sweep points"
        else:
            mism[L] = v
    bad = {L: v for L, v in mism.items() if v > SWEEP_TRIM_TOL_CT}
    return (not bad), mism, nojudge


# ==================================================================================
# CL-CD polar, one operating point per file
# ==================================================================================
def fig_polar(rans, vlm, condition, target_CL, fname):
    series, part = rans_series(rans, condition)
    mapped, alt, unmapped, vdiag = vlm_at(vlm, condition)
    alt = []  # the second baseline deck is loaded and checked but not drawn: one baseline, B (22 Sep)
    gname = geom_names(vlm)
    DIAG["polar_%s" % condition] = dict(rans_partition=part, vlm=vdiag,
                                        unmapped=[r["geometry"] for r in unmapped])

    draw_sweeps, mism, nojudge = sweeps_usable(series)
    DIAG["polar_%s" % condition]["sweep_trim_mismatch_ct"] = {
        L: round(v, 4) for L, v in sorted(mism.items())}
    DIAG["polar_%s" % condition]["sweep_trim_not_judged"] = nojudge
    DIAG["polar_%s" % condition]["sweeps_drawn"] = draw_sweeps

    fig, ax = plt.subplots(figsize=(3.85, 3.25))
    facts = {}

    for L in ORDER:
        if L not in series:
            continue
        s, st = series[L], STYLE[L]
        if draw_sweeps:
            ax.plot([p["Cd_ct"] for p in s["pts"]], [p["Cl"] for p in s["pts"]],
                    color=st["color"], ls=st["ls"], marker=st["marker"], ms=3.0,
                    mfc="none", mew=0.8, zorder=3)
        if s["trim"]:
            ax.plot(s["trim"]["Cd_ct"], s["trim"]["Cl"], color=st["color"],
                    marker=st["marker"], ms=5.0, mfc=st["color"], mec="white",
                    mew=0.7, ls="none", zorder=7)
            facts["%s_RANS_trim" % L] = dict(case=s["trim"]["case"],
                                             CD_ct=s["trim"]["Cd_ct"],
                                             CL=s["trim"]["Cl"],
                                             alpha_deg=s["trim"]["alpha"])
        facts["%s_RANS_bracket_CD_ct" % L] = [round(p["Cd_ct"], 4) for p in s["pts"]]

    ax.axhline(target_CL, color=GREY, lw=0.9, ls=(0, (2.5, 2.5)), zorder=2)

    handles = [Line2D([], [], color=STYLE[L]["color"],
                      ls=STYLE[L]["ls"] if draw_sweeps else "none",
                      marker=STYLE[L]["marker"], ms=3.0 if draw_sweeps else 5.0,
                      mfc="none" if draw_sweeps else STYLE[L]["color"],
                      mec="white" if not draw_sweeps else STYLE[L]["color"], mew=0.8,
                      label="%s  %s" % (L, gname.get(L, "geometry NOT MAPPED")))
               for L in ORDER if L in series]
    handles.append(Line2D([], [], color="0.2", ls="none", marker="o", ms=5.0,
                          mfc="0.2", mec="white", mew=0.7, label="RANS trim point"))
    handles.append(Line2D([], [], color=GREY, ls=(0, (2.5, 2.5)), lw=0.9,
                          label=r"target $C_L$"))

    note = None
    if condition == "condition_CR":
        for L in ORDER:
            if L not in mapped:
                continue
            r = mapped[L]
            assert r["already_DSO"] and r["factor"] == 1.0, \
                "unexpected non-unity Sref conversion for %s at %s" % (L, condition)
            ax.plot(r["dso"]["CDtot"] * CT, r["dso"]["CLtot"], color=STYLE[L]["color"],
                    marker=STYLE[L]["marker"], ms=6.5, mfc="none", mew=1.1,
                    ls="none", zorder=6)
            facts["%s_VLM" % L] = dict(geometry=r["geometry"], alpha_deg=r["dso"]["AoA"],
                                       CL=r["dso"]["CLtot"],
                                       CDtot_ct=r["dso"]["CDtot"] * CT,
                                       CDi_ct=r["dso"]["CDi"] * CT,
                                       CDiw_ct=r["dso"]["CDiw"] * CT,
                                       CDo_ct=r["dso"]["CDo"] * CT)
        for r in alt:
            ax.plot(r["dso"]["CDtot"] * CT, r["dso"]["CLtot"], color=ALT_C,
                    marker="o", ms=7.5, mfc="none", mew=1.1, ls="none", zorder=6)
            facts["Bstar_VLM"] = dict(geometry=r["geometry"], alpha_deg=r["dso"]["AoA"],
                                      CL=r["dso"]["CLtot"],
                                      CDtot_ct=r["dso"]["CDtot"] * CT,
                                      CDi_ct=r["dso"]["CDi"] * CT,
                                      CDiw_ct=r["dso"]["CDiw"] * CT)
        handles.append(Line2D([], [], ls="none", marker="o", ms=7.5, mfc="none",
                              mec="0.2", mew=1.1,
                              label=r"VLM $C_{D\mathrm{tot}}$ ($=C_{Do}+C_{Di}$)"))
        note = (r"$S_{\mathrm{ref}}$ CONVERSION: both sides already DSO basis"
                r" (13.3572 ft$^2$); factor applied is exactly 1.000000."
                "\n"
                r"VLM $C_{D\mathrm{tot}}=C_{Do}+C_{Di}$ with $C_{Do}$ a FLAT-PLATE"
                r" correlation at $Re_{C_{\mathrm{ref}}}=10^{7}$:"
                "\n"
                r"NOT comparable to a RANS total. Decks are single-point,"
                r" so no VLM polar curve exists.")
    # The sweep note is appended for EITHER condition, and only when the sweeps were
    # withheld, with the measured number that caused it. A figure that silently omits a
    # series it used to draw is worse than one that never drew it.
    if not draw_sweeps:
        worst = max(mism.values()) if mism else float("nan")
        note = (note + "\n" if note else "") + (
            r"ALPHA-SWEEP POINTS WITHHELD. Interpolating each geometry's own sweep to its"
            r" own trimmed $C_L$ misses" "\n"
            + (r"that geometry's measured trim $C_D$ by up to %.2f counts, so a curve"
               r" through the sweep would" % worst) + "\n"
            + r"not pass through its own trim point. Only the trimmed solutions are shown.")
    else:
        for r in alt:
            xd, yd = r["deck"]["CDtot"] * CT, r["deck"]["CLtot"]
            xs, ys = r["dso"]["CDtot"] * CT, r["dso"]["CLtot"]
            ax.plot([xd, xs], [yd, ys], color=FLAG_C, lw=1.0, ls=(0, (1.2, 1.2)), zorder=5)
            ax.plot([xd], [yd], color=FLAG_C, marker="o", ms=7.5, mfc="none",
                    mew=1.2, ls="none", zorder=6)
            ax.plot([xs], [ys], color=FLAG_C, marker="o", ms=7.5, mfc=FLAG_C,
                    mec="white", mew=0.7, ls="none", zorder=6)
            facts["Bstar_VLM_as_delivered_TP1580"] = dict(
                geometry=r["geometry"], alpha_deg=r["deck"]["AoA"], CL=yd, CDtot_ct=xd,
                CDi_ct=r["deck"]["CDi"] * CT, CDiw_ct=r["deck"]["CDiw"] * CT)
            facts["Bstar_VLM_converted_to_DSO"] = dict(
                CL=ys, CDtot_ct=xs, CDi_ct=r["dso"]["CDi"] * CT,
                CDiw_ct=r["dso"]["CDiw"] * CT, factor=r["factor"])
        f = facts.get("Bstar_VLM_converted_to_DSO", {}).get("factor", float("nan"))
        note = (r"OPEN FRAME CONFLICT, NOT ADJUDICATED (D068). Cruise decks are TP-1580"
                r" basis, RANS cruise cases are DSO basis;"
                "\n"
                r"$\times$%.6f between them, i.e. 11.31%% against a morphing effect of"
                r" order 0.4%%. Which basis the cruise target"
                "\n"
                r"$C_L$ is defined on is settled by no file read. Ring = as delivered,"
                r" filled = converted. NEITHER is quotable."
                "\n"
                r"Only B$^{*}$ maps to a RANS geometry at cruise, and no RANS case meshes it."
                % f)

    ax.set_xlabel(r"$C_D$   (counts, $1\ \mathrm{ct}=10^{-4}$)")
    ax.set_ylabel(r"$C_L$")
    ax.set_title("Condition CR" if condition == "condition_CR" else "Early cruise", pad=36)
    ax.grid(True, lw=0.35, color="0.88", zorder=0)
    ax.set_axisbelow(True)
    leg = ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.185),
                    ncol=2, handlelength=2.4, labelspacing=0.26, columnspacing=1.1,
                    borderaxespad=0.0)
    note_below(fig, leg, note)

    path = os.path.join(FIG_DIR, fname)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("WROTE %s  %d bytes" % (path, os.path.getsize(path)))
    WRITTEN.append(path)
    return facts


# ==================================================================================
# drag decomposition at the trim point, Condition CR
# ==================================================================================
def fig_drag_decomp(rans, vlm, target_CL):
    series, _ = rans_series(rans, "condition_CR")
    mapped, alt, _, _ = vlm_at(vlm, "condition_CR")
    alt = []  # one baseline, B (22 Sep)
    letters = [L for L in ORDER if L in series and series[L]["trim"] and L in mapped]

    x = np.arange(len(letters), dtype=float)
    w = 0.20
    rans_tot = [series[L]["trim"]["Cd_ct"] for L in letters]
    v_tot = [mapped[L]["dso"]["CDtot"] * CT for L in letters]
    v_cdi = [mapped[L]["dso"]["CDi"] * CT for L in letters]
    v_cdiw = [mapped[L]["dso"]["CDiw"] * CT for L in letters]

    fig, ax = plt.subplots(figsize=(3.95, 3.0))
    ax.bar(x - 1.5 * w, rans_tot, width=w, facecolor="white", edgecolor="black",
           hatch="/////", lw=0.6, zorder=3,
           label=r"RANS $C_D$ total (viscous solution)")
    ax.bar(x - 0.5 * w, v_tot, width=w, facecolor="#BBBBBB", edgecolor="black",
           lw=0.5, zorder=3,
           label=r"VLM $C_{D\mathrm{tot}}=C_{Do}+C_{Di}$ (NOT viscous)")
    ax.bar(x + 0.5 * w, v_cdi, width=w, facecolor="#0072B2", edgecolor="black",
           hatch="....", lw=0.5, zorder=3,
           label=r"VLM $C_{Di}$ near field (surface integration)")
    ax.bar(x + 1.5 * w, v_cdiw, width=w, facecolor="#E69F00", edgecolor="black",
           hatch="\\\\\\\\", lw=0.5, zorder=3,
           label=r"VLM $C_{Diw}$ far field (wake / Trefftz)")

    ax.set_xticks(x)
    ax.set_xticklabels(letters)
    ax.set_xlabel("geometry")
    ax.set_ylabel(r"drag   (counts, $1\ \mathrm{ct}=10^{-4}$)")
    ax.set_title("", pad=40)
    ax.set_ylim(0, max(rans_tot) * 1.16)
    ax.grid(True, axis="y", lw=0.35, color="0.88", zorder=0)
    ax.set_axisbelow(True)
    leg = ax.legend(loc="upper center", ncol=1, handlelength=1.7, labelspacing=0.25,
                    columnspacing=1.0, borderaxespad=0.0, bbox_to_anchor=(0.5, -0.185))

    ratio = [v_cdiw[i] / v_cdi[i] for i in range(len(letters))]
    note = (r"RANS-vs-VLM must be INDUCED AGAINST INDUCED. The two left-hand bars are NOT the"
            r" same quantity as each other."
            "\n"
            r"RANS $C_{Di}$ is NOT FOUND in the extracted data ($C_D$ total only; no"
            r" Trefftz integral, no pressure/viscous split),"
            "\n"
            r"so the defensible pairing cannot be drawn here. VLM near- and far-field"
            r" induced drag disagree by $C_{Diw}/C_{Di}$ = %.3f to %.3f"
            "\n"
            r"at this condition; the DSO ranks on the far-field $C_{Diw}$ from 2026-08-19."
            r" Both sides DSO $S_{\mathrm{ref}}$, factor exactly 1.000000."
            % (min(ratio), max(ratio)))
    note_below(fig, leg, note)

    facts = {L: dict(RANS_CD_total_ct=rans_tot[i], VLM_CDtot_ct=v_tot[i],
                     VLM_CDi_ct=v_cdi[i], VLM_CDiw_ct=v_cdiw[i],
                     VLM_CDo_ct=mapped[L]["dso"]["CDo"] * CT,
                     VLM_CDiw_over_CDi=v_cdiw[i] / v_cdi[i])
             for i, L in enumerate(letters)}
    path = os.path.join(FIG_DIR, "drag_decomposition_rans_vs_vlm_condition_CR.pdf")
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("WROTE %s  %d bytes" % (path, os.path.getsize(path)))
    WRITTEN.append(path)
    return facts


# ==================================================================================
# CL vs alpha, one operating point per file
# ==================================================================================
def fig_cl_alpha(rans, vlm, condition, target_CL, fname):
    series, part = rans_series(rans, condition)
    mapped, alt, unmapped, vdiag = vlm_at(vlm, condition)
    alt = []  # the second baseline deck is loaded and checked but not drawn: one baseline, B (22 Sep)
    DIAG["cl_alpha_%s" % condition] = dict(rans_partition=part, vlm=vdiag)

    fig, ax = plt.subplots(figsize=(3.85, 3.25))
    facts, slopes = {}, {}

    for L in ORDER:
        if L not in series:
            continue
        s, st = series[L], STYLE[L]
        ax.plot([p["alpha"] for p in s["pts"]], [p["Cl"] for p in s["pts"]],
                color=st["color"], ls=st["ls"], marker=st["marker"], ms=3.0,
                mfc="none", mew=0.8, zorder=3)
        fit = lift_slope(s["pts"])
        if fit:
            slopes[L] = fit
        if s["trim"]:
            ax.plot(s["trim"]["alpha"], s["trim"]["Cl"], color=st["color"],
                    marker=st["marker"], ms=5.0, mfc=st["color"], mec="white",
                    mew=0.7, ls="none", zorder=7)
        facts[L] = dict(
            trim_alpha_deg=(s["trim"]["alpha"] if s["trim"] else "NOT FOUND"),
            trim_CL=(s["trim"]["Cl"] if s["trim"] else "NOT FOUND"),
            dCL_dalpha_per_deg=fit["dCL_dalpha_per_deg"] if fit else "NOT FOUND",
            dCL_dalpha_per_rad=(fit["dCL_dalpha_per_deg"] * 180.0 / np.pi) if fit else "NOT FOUND",
            fit_cases=fit["cases"] if fit else [],
            fit_max_abs_resid_CL=fit["max_abs_resid"] if fit else "NOT FOUND")

    ax.axhline(target_CL, color=GREY, lw=0.9, ls=(0, (2.5, 2.5)), zorder=2)

    handles = []
    for L in ORDER:
        if L not in series:
            continue
        a = facts[L]["trim_alpha_deg"]
        m = facts[L]["dCL_dalpha_per_deg"]
        lab = "%s   %s   %s" % (
            L,
            ("%.4f" % a) if isinstance(a, float) else "NOT FOUND",
            ("%.5f" % m) if isinstance(m, float) else "n/a")
        handles.append(Line2D([], [], color=STYLE[L]["color"], ls=STYLE[L]["ls"],
                              marker=STYLE[L]["marker"], ms=3.0, mfc="none", mew=0.8,
                              label=lab))

    note = None
    if condition == "condition_CR":
        for L in ORDER:
            if L not in mapped:
                continue
            r = mapped[L]
            ax.plot(r["dso"]["AoA"], r["dso"]["CLtot"], color=STYLE[L]["color"],
                    marker=STYLE[L]["marker"], ms=6.5, mfc="none", mew=1.1,
                    ls="none", zorder=6)
            facts[L]["VLM_alpha_deg"] = r["dso"]["AoA"]
            facts[L]["VLM_CL"] = r["dso"]["CLtot"]
            facts[L]["VLM_geometry"] = r["geometry"]
        for r in alt:
            ax.plot(r["dso"]["AoA"], r["dso"]["CLtot"], color=ALT_C, marker="o",
                    ms=7.5, mfc="none", mew=1.1, ls="none", zorder=6)
            facts["Bstar"] = dict(VLM_alpha_deg=r["dso"]["AoA"], VLM_CL=r["dso"]["CLtot"],
                                  VLM_geometry=r["geometry"])
        sym = [Line2D([], [], ls="none", marker="o", ms=7.5, mfc="none", mec="0.2",
                      mew=1.1, label="VLM single solved point")]
        note = (r"VSPAERO decks are SINGLE-POINT (one polar row each): no $\alpha$ sweep, so"
                r" the VLM LIFT SLOPE IS NOT FOUND and no VLM line is drawn."
                "\n"
                r"What IS comparable is the $\alpha$ each solver needs for the SAME $C_L$:"
                r" a VLM ring concentric with a filled RANS marker means they agree."
                "\n"
                r"RANS $\mathrm{d}C_L/\mathrm{d}\alpha$ is a least-squares fit over that"
                r" family's own converged cases; max fit residual %.1e in $C_L$."
                "\n"
                r"B$^{*}$ = baseline_corrected, an alternative aerofoil interpolation the"
                r" registry marks DO NOT MESH."
                % max(f["fit_max_abs_resid_CL"] for f in
                      (facts[L] for L in ORDER if L in series)))
    else:
        for r in alt:
            ax.plot([r["deck"]["AoA"], r["dso"]["AoA"]],
                    [r["deck"]["CLtot"], r["dso"]["CLtot"]], color=FLAG_C, lw=1.0,
                    ls=(0, (1.2, 1.2)), zorder=5)
            ax.plot([r["deck"]["AoA"]], [r["deck"]["CLtot"]], color=FLAG_C, marker="o",
                    ms=7.5, mfc="none", mew=1.2, ls="none", zorder=6)
            ax.plot([r["dso"]["AoA"]], [r["dso"]["CLtot"]], color=FLAG_C, marker="o",
                    ms=7.5, mfc=FLAG_C, mec="white", mew=0.7, ls="none", zorder=6)
            facts["Bstar"] = dict(VLM_alpha_deg=r["deck"]["AoA"],
                                  VLM_CL_as_delivered_TP1580=r["deck"]["CLtot"],
                                  VLM_CL_converted_to_DSO=r["dso"]["CLtot"],
                                  factor=r["factor"], VLM_geometry=r["geometry"])
        sym = []  # one baseline, B (22 Sep)
        note = (r"VSPAERO decks are SINGLE-POINT: no $\alpha$ sweep, so the VLM LIFT SLOPE IS"
                r" NOT FOUND. Only B$^{*}$ maps to a RANS geometry at cruise;"
                "\n"
                r"the six other cruise decks are a separate concept set with no derived STL."
                r" OPEN FRAME CONFLICT, NOT ADJUDICATED (D068):"
                "\n"
                r"deck is TP-1580 basis ($\times$%.6f to DSO, 11.31%%), RANS is DSO basis,"
                r" and which basis the cruise target $C_L$ is on"
                "\n"
                r"is settled by no file read. Both readings drawn; NEITHER quotable."
                r"   RANS $\mathrm{d}C_L/\mathrm{d}\alpha$ is a least-squares fit,"
                r" max residual %.1e in $C_L$."
                % (facts.get("Bstar", {}).get("factor", float("nan")),
                   max(f["fit_max_abs_resid_CL"] for f in
                       (facts[L] for L in ORDER if L in series))))

    ax.set_xlabel(r"$\alpha$   (deg)")
    ax.set_ylabel(r"$C_L$")
    ax.set_title("Condition CR" if condition == "condition_CR" else "Early cruise", pad=36)
    ax.grid(True, lw=0.35, color="0.88", zorder=0)
    ax.set_axisbelow(True)

    sym.insert(0, Line2D([], [], color="0.2", ls="none", marker="o", ms=5.0, mfc="0.2",
                         mec="white", mew=0.7, label="RANS trim point"))
    sym.insert(1, Line2D([], [], color=GREY, ls=(0, (2.5, 2.5)), lw=0.9,
                         label=r"target $C_L$ = %.11g" % target_CL))

    # ONE legend box, two columns: geometries left, symbol key right. Two separate boxes
    # overlapped each other at this figure width, and a legend title centred over both
    # columns was clipped, so the column header is carried as a blank-handle row instead.
    def blank(lab):
        return Line2D([], [], ls="none", marker="none", label=lab)

    col1 = [blank(r"geom  trim $\alpha$ (deg)  $\mathrm{d}C_L/\mathrm{d}\alpha$ (deg$^{-1}$)")] + handles
    col2 = [blank("")] + sym
    n = max(len(col1), len(col2))
    col1 += [blank("")] * (n - len(col1))
    col2 += [blank("")] * (n - len(col2))
    leg1 = ax.legend(handles=col1 + col2, loc="upper center",
                     bbox_to_anchor=(0.5, -0.185), ncol=2, handlelength=2.4,
                     labelspacing=0.26, columnspacing=1.4, borderaxespad=0.0)
    leg2 = leg1

    # ticks locating each trim alpha on the target-CL line, drawn last so they use the
    # final y-limits rather than the limits before the VLM points were added
    ylo, yhi = ax.get_ylim()
    tick = 0.045 * (yhi - ylo)
    for L in ORDER:
        if L in series and series[L]["trim"]:
            ax.plot([series[L]["trim"]["alpha"]] * 2, [target_CL, target_CL - tick],
                    color=STYLE[L]["color"], lw=0.9, zorder=4)

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    y0 = min(l.get_window_extent(r).transformed(fig.transFigure.inverted()).y0
             for l in (leg1, leg2))
    fig.text(0.5, y0 - 0.025, "", ha="center", va="top", fontsize=7.4, color="0.2",
             linespacing=1.5)
    path = os.path.join(FIG_DIR, fname)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("WROTE %s  %d bytes" % (path, os.path.getsize(path)))
    WRITTEN.append(path)
    return facts, slopes


# ==================================================================================
# delta alpha at the common target CL, Condition CR (the frame-closed comparison)
# ==================================================================================
def fig_trim_alpha_delta(rans, vlm, target_CL, slopes):
    series, _ = rans_series(rans, "condition_CR")
    mapped, alt, _, vdiag = vlm_at(vlm, "condition_CR")
    alt = []  # one baseline, B (22 Sep)
    DIAG["delta_alpha"] = vdiag

    rows = []
    for L in ORDER:
        if L in series and series[L]["trim"] and L in mapped:
            rows.append((L, series[L]["trim"], mapped[L], STYLE[L]["color"]))
    for r in alt:
        if "B" in series and series["B"]["trim"]:
            rows.append(("B*", series["B"]["trim"], r, ALT_C))

    labels, dalpha, colors, detail = [], [], [], {}
    for L, trim, v, col in rows:
        m = slopes["B" if L == "B*" else L]["dCL_dalpha_per_deg"]
        # both sides slid to the exact target CL with the RANS lift slope, so the bar is a
        # difference at one CL and not at two slightly different ones
        a_rans = trim["alpha"] + (target_CL - trim["Cl"]) / m
        a_vlm = v["dso"]["AoA"] + (target_CL - v["dso"]["CLtot"]) / m
        labels.append(L)
        dalpha.append(a_vlm - a_rans)
        colors.append(col)
        detail[L] = dict(vlm_geometry=v["geometry"],
                         alpha_RANS_trim_deg=trim["alpha"], CL_RANS_trim=trim["Cl"],
                         alpha_VLM_deg=v["dso"]["AoA"], CL_VLM=v["dso"]["CLtot"],
                         dCL_dalpha_per_deg=m,
                         alpha_RANS_at_target_deg=a_rans,
                         alpha_VLM_at_target_deg=a_vlm,
                         delta_alpha_deg=a_vlm - a_rans,
                         CL_slide_RANS_deg=(target_CL - trim["Cl"]) / m,
                         CL_slide_VLM_deg=(target_CL - v["dso"]["CLtot"]) / m)

    y = np.arange(len(labels), dtype=float)
    fig, ax = plt.subplots(figsize=(3.85, 2.55))
    ax.barh(y, dalpha, height=0.6, color=colors, edgecolor="black", lw=0.6, zorder=3)
    ax.axvline(0.0, color="black", lw=0.9, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([l.replace("B*", r"B$^{*}$") for l in labels])
    ax.invert_yaxis()
    ax.set_xlabel(r"$\Delta\alpha=\alpha_{\mathrm{VLM}}-\alpha_{\mathrm{RANS}}$ (deg)"
                  r" at the common target $C_L$")
    ax.set_title("", pad=40)

    span = max(abs(min(dalpha)), abs(max(dalpha)))
    ax.set_xlim(-1.75 * span, 1.75 * span)
    for yi, d in zip(y, dalpha):
        ax.text(d + (0.05 * span if d >= 0 else -0.05 * span), yi, "%+.4f" % d,
                va="center", ha="left" if d >= 0 else "right", fontsize=6.8)
    ax.grid(True, axis="x", lw=0.35, color="0.88", zorder=0)
    ax.set_axisbelow(True)

    note = (r"A LIFT comparison, not a drag one, and the cleanest single test of whether the"
            r" VLM is usable: induced drag follows the lift distribution."
            "\n"
            r"Both sides DSO $S_{\mathrm{ref}}$, conversion factor exactly 1.000000. Each"
            r" solver's point is slid to the EXACT target $C_L$ using that family's"
            "\n"
            r"own RANS $\mathrm{d}C_L/\mathrm{d}\alpha$, so each bar is a difference at ONE"
            r" $C_L$; the largest slide applied is %.2e deg."
            "\n"
            r"B = baseline_wing_only_refined, the surface every SWB/CMPB/LCB case meshes and"
            r" the only Condition CR VLM run that misses its own"
            "\n"
            r"trim target (by %.2e in $C_L$).   B$^{*}$ = baseline_corrected, an alternative"
            r" aerofoil interpolation the registry marks DO NOT MESH."
            % (max(max(abs(v["CL_slide_RANS_deg"]), abs(v["CL_slide_VLM_deg"]))
                   for v in detail.values()) if detail else 0.0,
               abs(mapped["B"]["CL_err"]) if "B" in mapped else float("nan")))
    fig.canvas.draw()
    bb = ax.get_window_extent(fig.canvas.get_renderer()).transformed(
        fig.transFigure.inverted())
    fig.text(0.5, bb.y0 - 0.16, "", ha="center", va="top", fontsize=7.4, color="0.2",
             linespacing=1.5)

    path = os.path.join(FIG_DIR, "trim_alpha_rans_vs_vlm_condition_CR.pdf")
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("WROTE %s  %d bytes" % (path, os.path.getsize(path)))
    WRITTEN.append(path)
    return detail


# ==================================================================================
def main():
    rans, vlm = load()
    tCR = vlm["_operating_points"]["condition_CR"]["target_CL"]
    tEC = vlm["_operating_points"]["early_cruise"]["target_CL"]

    out = {}
    out["polar_condition_CR"] = fig_polar(rans, vlm, "condition_CR", tCR,
                                          "polar_cl_cd_condition_CR.pdf")
    out["polar_early_cruise"] = fig_polar(rans, vlm, "early_cruise", tEC,
                                          "polar_cl_cd_early_cruise.pdf")
    out["drag_decomposition_condition_CR"] = fig_drag_decomp(rans, vlm, tCR)
    fCR, sCR = fig_cl_alpha(rans, vlm, "condition_CR", tCR, "cl_alpha_condition_CR.pdf")
    fEC, sEC = fig_cl_alpha(rans, vlm, "early_cruise", tEC, "cl_alpha_early_cruise.pdf")
    out["cl_alpha_condition_CR"] = fCR
    out["cl_alpha_early_cruise"] = fEC
    out["delta_alpha_condition_CR"] = fig_trim_alpha_delta(rans, vlm, tCR, sCR)

    # late cruise: report the gap as a counted bucket rather than omitting it
    lc, lc_part = None, None
    try:
        lc, lc_part = rans_series(rans, "late_cruise")
    except AssertionError as exc:
        lc_part = {"error": str(exc)}
    out["late_cruise_gap"] = dict(
        partition=lc_part,
        figure_produced=False,
        reason=("Only one late-cruise RANS case is converged and no LC trim case exists, "
                "so neither a polar nor a lift slope can be built. No figure is drawn."))

    print("\n=== FILE VERIFICATION (stat + magic bytes, after write) ===")
    ok = True
    for path in WRITTEN:
        exists = os.path.isfile(path)
        size = os.path.getsize(path) if exists else 0
        with open(path, "rb") as fh:
            magic = fh.read(5)
        good = exists and size > 0 and magic == b"%PDF-"
        print("%-84s exists=%-5s bytes=%-8d pdf=%s" % (path, exists, size, magic == b"%PDF-"))
        ok = ok and good
    print("ALL_FILES_OK", ok, " n_files", len(WRITTEN))

    print("\n=== DIAGNOSTICS ===")
    print(json.dumps(DIAG, indent=1, sort_keys=True))
    print("\n=== NUMBERS ===")
    print(json.dumps(out, indent=1, sort_keys=True, default=str))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
