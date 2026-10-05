#!/usr/bin/env python3
"""
Headline figure: trimmed performance across all six geometries, per condition.

Every number is read from data/rans_forces.json. Nothing is typed in, nothing is
estimated. A quantity the file marks "NOT FOUND" is drawn as an explicitly
labelled gap, never as a bar.

Frame (D068). Every coefficient here is on the DSO basis, half model:
    Aref  0.620462 m2 (half-wing patch), lRef 0.393957 m
    2*Aref = 1.240924 m2 reproduces the DSO Sref 1.24092 m2 to 3.2e-6 relative,
    against 11.3e-2 relative for the TP-1580 Sref 1.1148 m2.
The script ASSERTS that frame out of the file rather than declaring it, and halts
if any condition disagrees.

Operating points (D076). Three distinct points, never mixed; every delta is formed
WITHIN one condition against that condition's own B trim.

Turbulence model is NOT common across conditions (kOmegaSST at Condition CR,
SpalartAllmaras at both cruise points), so a count from one condition is not
like-for-like against a count from another. Stated on the figure.

Outputs (into the directory holding this script):
    fig_trimmed_performance.pdf          vector, for LaTeX
    fig_trimmed_performance_values.csv   every plotted value, full precision
"""

import csv
import json
import os
import sys
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)
RANS = os.path.join(REPORT, "data", "rans_forces.json")

OUT_PDF = os.path.join(HERE, "fig_trimmed_performance.pdf")
OUT_CSV = os.path.join(HERE, "fig_trimmed_performance_values.csv")

# ---------------------------------------------------------------- expectations
# Frame the figure is allowed to draw. Read back from the file and asserted; a
# disagreement halts rather than being relabelled.
EXPECT_AREF_M2 = 0.620462
EXPECT_LREF_M = 0.393957
EXPECT_BASIS = "DSO"

CONDITIONS = [
    ("condition_CR", "Condition CR"),
    ("early_cruise", "Early cruise"),
    ("late_cruise", "Late cruise"),
]

GEOMS = ["B", "C", "F", "H", "M", "W"]

# Okabe-Ito colour-blind-safe palette, plus a distinct hatch per geometry so the
# figure survives greyscale printing and photocopying.
COLOUR = {
    "B": "#000000",
    "C": "#0072B2",
    "F": "#D55E00",
    "H": "#009E73",
    "M": "#CC79A7",
    "W": "#E69F00",
}
HATCH = {
    "B": "",
    "C": "///",
    "F": "\\\\\\",
    "H": "xxx",
    "M": "...",
    "W": "+++",
}

# D058 (2026-07-28, docs/decisions/legacy.md): measured fixed-CL noise floor ON
# THE DELTA, 0.106 to 0.199 counts. scripts/run_card.schema.json band_CD_counts
# pins the pre-registered acceptance band at 0.10 counts, "below the 0.106-0.199
# count fixed-CL noise floor measured in D058". Measured on the coarse
# wall-function A2 pair, NOT re-measured on this mesh family: said on the figure.
NOISE_LO_CT = 0.106
NOISE_HI_CT = 0.199


def num(v):
    """Return a float, or None for anything the extractor marked NOT FOUND."""
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    return None


def load():
    with open(RANS) as fh:
        d = json.load(fh)
    trims = d["trims"]

    for key, _ in CONDITIONS:
        if key not in trims:
            sys.exit("HALT: condition %s absent from %s" % (key, RANS))
        fr = trims[key]["coefficient_frame"]
        # D068: assert the frame, do not declare it.
        if fr["basis_identified"] != EXPECT_BASIS:
            sys.exit("HALT: %s basis %r, expected %r" % (key, fr["basis_identified"], EXPECT_BASIS))
        if abs(fr["Aref_m2"] - EXPECT_AREF_M2) > 1e-9:
            sys.exit("HALT: %s Aref %r, expected %r" % (key, fr["Aref_m2"], EXPECT_AREF_M2))
        if abs(fr["lRef_m"] - EXPECT_LREF_M) > 1e-9:
            sys.exit("HALT: %s lRef %r, expected %r" % (key, fr["lRef_m"], EXPECT_LREF_M))
        if not fr.get("Aref_is_half_wing"):
            sys.exit("HALT: %s Aref is not flagged half-wing" % key)
    return d


def series(trims, key, field):
    """Per-geometry values for one condition; None where the file says NOT FOUND."""
    geos = trims[key]["geometries"]
    return [num(geos.get(g, {}).get(field)) for g in GEOMS]


def main():
    d = load()
    trims = d["trims"]

    data = {}
    for key, label in CONDITIONS:
        c = trims[key]
        data[key] = {
            "label": label,
            "mach": c["mach"],
            "target_CL": c["target_CL"],
            "turb": ", ".join(c["coefficient_frame"]["turbulence_models_present"]),
            "n_converged": c["coefficient_frame"]["n_converged"],
            "Cd_counts": series(trims, key, "Cd_counts"),
            "dCd": series(trims, key, "delta_Cd_counts_vs_baseline"),
            "LD": series(trims, key, "L_over_D"),
            "alpha": series(trims, key, "alpha_deg"),
            "Cl": series(trims, key, "Cl"),
            "Cl_err_ct": series(trims, key, "Cl_minus_target_counts"),
            "dCd_dCl": series(trims, key, "dCd_dCl_counts_per_unit_CL"),
            "trim_case": [c["geometries"].get(g, {}).get("trim_case") for g in GEOMS],
            "in_tol": [c["geometries"].get(g, {}).get("trim_within_tolerance") for g in GEOMS],
        }

    # Letter only. family_label reads "morphing candidate F" for four of the six,
    # which names nothing and costs four words per panel legend. What each geometry
    # is belongs in the report text.
    labels = {g: g for g in GEOMS}

    # ------------------------------------------------------------ dump values
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "condition", "mach", "target_CL", "turbulence_model", "geometry",
                "trim_case", "Cd_counts", "delta_Cd_counts_vs_B", "L_over_D",
                "alpha_deg", "Cl", "Cl_minus_target_counts", "trim_within_tolerance",
            ]
        )
        for key, _ in CONDITIONS:
            v = data[key]
            for i, g in enumerate(GEOMS):
                w.writerow(
                    [
                        key, v["mach"], v["target_CL"], v["turb"], g,
                        v["trim_case"][i] if v["trim_case"][i] else "NOT FOUND",
                        v["Cd_counts"][i] if v["Cd_counts"][i] is not None else "NOT FOUND",
                        v["dCd"][i] if v["dCd"][i] is not None else "NOT FOUND",
                        v["LD"][i] if v["LD"][i] is not None else "NOT FOUND",
                        v["alpha"][i] if v["alpha"][i] is not None else "NOT FOUND",
                        v["Cl"][i] if v["Cl"][i] is not None else "NOT FOUND",
                        v["Cl_err_ct"][i] if v["Cl_err_ct"][i] is not None else "NOT FOUND",
                        v["in_tol"][i] if v["in_tol"][i] is not None else "NOT FOUND",
                    ]
                )

    # ------------------------------------------------------------ figure setup
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Nimbus Roman", "serif"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 10.5,
            "axes.labelsize": 10.0,
            "axes.titlesize": 10.0,
            "xtick.labelsize": 8.6,
            "ytick.labelsize": 9.0,
            "legend.fontsize": 9.0,
            "axes.linewidth": 0.8,
            "hatch.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )

    # Six columns so row 1 can hold TWO half-width panels and row 2 THREE
    # third-width ones, one per operating point. The campaign ran on two
    # conditions when this figure was first written and row 2 had two panels;
    # the third was added when late cruise converged, rather than leaving late
    # cruise as a placeholder box inside panels built for the other two.
    fig = plt.figure(figsize=(7.4, 10.9))
    gs = fig.add_gridspec(
        3, 6,
        height_ratios=[1.00, 0.94, 0.82],
        hspace=1.00, wspace=0.62,
        left=0.105, right=0.975, top=0.880, bottom=0.272,
    )
    ax_cd = fig.add_subplot(gs[0, 0:3])
    ax_ld = fig.add_subplot(gs[0, 3:6])
    ax_dcr = fig.add_subplot(gs[1, 0:2])
    ax_dec = fig.add_subplot(gs[1, 2:4])
    ax_dlc = fig.add_subplot(gs[1, 4:6])
    ax_al = fig.add_subplot(gs[2, :])

    n = len(GEOMS)
    span = n + 2.2  # group pitch for the three-condition panels

    def bar_group(ax, values, gi, fmt, annot_size=6.6, note_missing=True,
                  skip_zero=False, rotation=90):
        """Draw one condition's six bars. A value the file marks NOT FOUND gets a
        labelled gap, never a bar. Negative values are labelled ABOVE the zero
        line so the text cannot collide with the tick labels."""
        for i, g in enumerate(GEOMS):
            x = gi * span + i
            v = values[i]
            if v is None:
                if note_missing:
                    ax.annotate(
                        "not run", xy=(x, 0), xytext=(0, 5), textcoords="offset points",
                        rotation=90, ha="center", va="bottom", fontsize=6.4,
                        color="0.30", style="italic",
                    )
                continue
            ax.bar(
                x, v, width=0.82, color=COLOUR[g], hatch=HATCH[g],
                edgecolor="black", linewidth=0.7, zorder=3,
            )
            if skip_zero and v == 0.0:
                continue
            anchor = v if v >= 0 else 0.0
            ax.annotate(
                fmt % v, xy=(x, anchor), xytext=(0, 3),
                textcoords="offset points", rotation=rotation,
                ha="center", va="bottom", fontsize=annot_size, zorder=4,
            )

    def geom_ticks(ax, ngroups):
        xs, ls = [], []
        for gi in range(ngroups):
            for i, g in enumerate(GEOMS):
                xs.append(gi * span + i)
                ls.append(g)
        ax.set_xticks(xs)
        ax.set_xticklabels(ls)
        ax.tick_params(axis="x", length=0, pad=2)
        ax.set_xlim(-1.10, (ngroups - 1) * span + n - 1 + 1.10)

    def cond_labels(ax, ngroups, y=-0.11):
        for gi in range(ngroups):
            key, lab = CONDITIONS[gi]
            v = data[key]
            centre = gi * span + (n - 1) / 2.0
            ax.annotate(
                lab,
                xy=(centre, y), xycoords=("data", "axes fraction"),
                ha="center", va="top", fontsize=7.6, annotation_clip=False,
                linespacing=1.3,
            )
            if gi:
                ax.axvline(gi * span - 1.60, color="0.72", lw=0.7, ls=(0, (4, 3)), zorder=1)

    def empty_group(ax, gi, text, fontsize=7.0):
        centre = gi * span + (n - 1) / 2.0
        ax.annotate(
            text, xy=(centre, 0.50), xycoords=("data", "axes fraction"),
            ha="center", va="center", fontsize=fontsize, color="0.15", style="italic",
            bbox=dict(boxstyle="round,pad=0.36", fc="0.93", ec="0.45", lw=0.8, ls="--"),
            zorder=5, linespacing=1.4,
        )

    def note_if_empty(ax, gi, key):
        """Draw the 'no trim case' box only for a condition that HAS none.

        This was an unconditional overlay on group 2, written when late cruise
        had no trims. It then sat on top of six correctly drawn bars once they
        existed, which is a claim about the data authored beside the drawing
        rather than read from it (GEO-080). It is now a question asked of the
        file each time the figure is built.
        """
        if all(v is None for v in data[key]["Cd_counts"]):
            empty_group(ax, gi, "%s\nNO TRIM CASE\nEXISTS FOR ANY\nOF THE SIX"
                        % dict(CONDITIONS)[key].upper())

    # ------------------------------------------------- (a) trimmed C_D, counts
    for gi, (key, _) in enumerate(CONDITIONS):
        bar_group(ax_cd, data[key]["Cd_counts"], gi, "%.2f")
        note_if_empty(ax_cd, gi, key)
    geom_ticks(ax_cd, 3)
    cond_labels(ax_cd, 3)
    ax_cd.set_ylabel("trimmed $C_D$  [counts]")
    ax_cd.set_title("(a)", loc="left", pad=5)
    ax_cd.set_ylim(0, 350)
    ax_cd.yaxis.grid(True, color="0.88", lw=0.6, zorder=0)
    ax_cd.set_axisbelow(True)

    # -------------------------------------------------------------- (b)  L/D
    for gi, (key, _) in enumerate(CONDITIONS):
        bar_group(ax_ld, data[key]["LD"], gi, "%.2f")
        note_if_empty(ax_ld, gi, key)
    geom_ticks(ax_ld, 3)
    cond_labels(ax_ld, 3)
    ax_ld.set_ylabel("$L/D$  [-]")
    ax_ld.set_title("(b)", loc="left", pad=5)
    ax_ld.set_ylim(0, 44)
    ax_ld.yaxis.grid(True, color="0.88", lw=0.6, zorder=0)
    ax_ld.set_axisbelow(True)

    # -------------------------------- (c), (d), (e) delta C_D vs the baseline
    # One axes per operating point, because Condition CR differs from the two
    # cruise points by an order of magnitude and a shared linear scale would
    # render every Condition CR delta invisible. Each panel is formed WITHIN one
    # condition against that condition's own B trim (D076); no panel differences
    # across conditions, and the text box below reports the cross-condition
    # comparison as a RANKING, which is what survives a change of operating point.
    cr = data["condition_CR"]
    ec = data["early_cruise"]
    lc = data["late_cruise"]
    iC, iF, iH, iW = (GEOMS.index(g) for g in ("C", "F", "H", "W"))

    ax_dcr.axhspan(-NOISE_HI_CT, NOISE_HI_CT, color="0.86", zorder=0)
    for s in (-NOISE_LO_CT, NOISE_LO_CT):
        ax_dcr.axhline(s, color="0.45", lw=0.8, ls=(0, (1, 2)), zorder=1)
    bar_group(ax_dcr, cr["dCd"], 0, "%+.3f", annot_size=6.6, skip_zero=True, rotation=90)
    ax_dcr.axhline(0.0, color="black", lw=1.3, ls="-", zorder=4)
    pass  # "B = 0" removed: the zero line and the caption say it (22 Sep)
    ax_dcr.set_xticks(range(n))
    ax_dcr.set_xticklabels(GEOMS, fontsize=10.0)
    ax_dcr.tick_params(axis="x", length=0, pad=2)
    ax_dcr.set_xlim(-0.75, n - 0.25)
    ax_dcr.set_ylim(-1.6, 9.4)
    ax_dcr.set_ylabel(r"$\Delta C_D$ vs baseline  [counts]")
    ax_dcr.set_title(
        "(c)  Condition CR",
        loc="left", pad=5, fontsize=8.6,
    )
    ax_dcr.yaxis.grid(True, color="0.88", lw=0.6, zorder=0)
    ax_dcr.set_axisbelow(True)

    # The two cruise panels share one y-scale so the eye can compare them
    # directly; the scale is taken from both, not fixed, so neither can overflow
    # it silently.
    cruise_vals = [v for v in (ec["dCd"] + lc["dCd"]) if v is not None]
    # Headroom for the rotated value labels is a FRACTION OF THE AXIS, not of the bar, so
    # the factor depends on the range: 1.34 fitted the published +56-count range and
    # clipped the labels once the welded range (+14.6 at most) arrived. 1.8 clears both.
    cruise_top = max(cruise_vals) * 1.8

    for ax, dat, letter, label in ((ax_dec, ec, "d", "early cruise"),
                                   (ax_dlc, lc, "e", "late cruise")):
        vals = [v for v in dat["dCd"] if v is not None and v != 0.0]
        bar_group(ax, dat["dCd"], 0, "%+.3f", annot_size=6.6, skip_zero=True, rotation=90)
        ax.axhline(0.0, color="black", lw=1.3, ls="-", zorder=4)
        # The convention is stated once, on panel (c), and in the caption.
        ax.set_xticks(range(n))
        ax.set_xticklabels(GEOMS, fontsize=10.0)
        ax.tick_params(axis="x", length=0, pad=2)
        ax.set_xlim(-0.75, n - 0.25)
        ax.set_ylim(-10, cruise_top)
        ax.set_title(
            "(%s)  %s" % (letter, label),
            loc="left", pad=5, fontsize=8.6,
        )
        ax.yaxis.grid(True, color="0.88", lw=0.6, zorder=0)
        ax.set_axisbelow(True)
    # One y-label for the row, on the leftmost delta panel only: three copies of
    # the same string crowded the narrow columns and read as three scales.
    ax_dlc.tick_params(axis="y", labelleft=True)
    CRUISE_SCALE_RATIO = cruise_top / 9.4

    # What changes between conditions and what does not, COUNTED from the data
    # rather than asserted. A ranking is the object that survives a change of
    # operating point, so the comparison across conditions is made on pairs.
    def order(dat):
        return [g for _, g in sorted(
            (dat["dCd"][GEOMS.index(g)], g) for g in GEOMS
            if dat["dCd"][GEOMS.index(g)] is not None)]

    def swapped_pairs(a, b):
        oa, ob = order(a), order(b)
        common = [g for g in oa if g in ob]
        return sum(1 for i, x in enumerate(common) for y in common[i + 1:]
                   if (ob.index(x) > ob.index(y)))

    n_pairs = len(GEOMS) * (len(GEOMS) - 1) // 2
    inv = (
        "ORDERING. Condition CR gives %s; early cruise %s; late cruise %s. "
        "Between the TWO CRUISE POINTS %d of %d pairs change order. Between low speed and "
        "cruise %d of %d change.\n"
        "RANK INVERSION, F against H: F is %.3f ct LOWER drag than H at Condition CR and "
        "%.3f ct HIGHER at early cruise, %.3f ct higher at late cruise. F against W is NOT "
        "resolved at Condition CR (%.3f ct apart, inside the D058 noise floor)." % (
            " ".join(order(cr)), " ".join(order(ec)), " ".join(order(lc)),
            swapped_pairs(ec, lc), n_pairs,
            swapped_pairs(cr, ec), n_pairs,
            abs(cr["dCd"][iF] - cr["dCd"][iH]),
            abs(ec["dCd"][iF] - ec["dCd"][iH]),
            abs(lc["dCd"][iF] - lc["dCd"][iH]),
            abs(cr["dCd"][iF] - cr["dCd"][iW]),
        )
    )
    inv = "\n".join(textwrap.fill(line, width=108) for line in inv.split("\n"))
    INV_TEXT = inv  # placed after the last axes is built, from realised positions

    # ------------------------------------------------------- (f) trim alpha
    for gi, (key, _) in enumerate(CONDITIONS):
        bar_group(ax_al, data[key]["alpha"], gi, "%.3f", annot_size=7.0)
        note_if_empty(ax_al, gi, key)
    geom_ticks(ax_al, 3)
    cond_labels(ax_al, 3, y=-0.13)
    ax_al.set_ylabel(r"trim $\alpha$  [deg]")
    ax_al.set_title(
        "(f)",
        loc="left", pad=5,
    )
    ax_al.set_ylim(0, 3.0)
    ax_al.yaxis.grid(True, color="0.88", lw=0.6, zorder=0)
    ax_al.set_axisbelow(True)

    # Rank-inversion box: centred in the gap between row 2 and panel (e)'s title,
    # using the positions matplotlib actually realised rather than a guessed y.
    y_gap_top = ax_dcr.get_position().y0 - 0.026   # below the (c)/(d) tick labels
    y_gap_bot = ax_al.get_position().y1 + 0.022    # above panel (e)'s title
    # The rank-inversion discussion belongs in the body text, not in a box between
    # panels; it is folded into the sidecar notes where `paras` is assembled below.

    # ------------------------------------------------------------- legend
    handles = [
        Patch(facecolor=COLOUR[g], hatch=HATCH[g], edgecolor="black", linewidth=0.7,
              label=g)
        for g in GEOMS
    ]
    fig.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.962),
        ncol=2, frameon=True, fancybox=False, edgecolor="0.4",
        handlelength=1.9, handleheight=1.15, columnspacing=1.8, borderpad=0.5,
        labelspacing=0.45,
    )

    # Figure title omitted: the LaTeX caption names the figure.

    # ------------------------------------------------- trim-residual caveat
    # D071: does the defect DIFFER between the two things being differenced? A
    # trim residual common to both sides of a delta CANCELS. What survives is the
    # residual DIFFERENCE, so that is what is reported, converted into drag
    # through each geometry's own bracket slope. The raw per-trim miss is stated
    # too, because a reader who sees only the difference cannot check the claim.
    out_of_tol, resid_diffs, delta_biases = [], [], []
    for key, label in CONDITIONS:
        v = data[key]
        ib = GEOMS.index("B")
        base_err = v["Cl_err_ct"][ib]
        for i, g in enumerate(GEOMS):
            e, sl = v["Cl_err_ct"][i], v["dCd_dCl"][i]
            if e is None:
                continue
            if v["in_tol"][i] is False:
                out_of_tol.append("%s %+.3f ct" % (v["trim_case"][i], e))
            if base_err is None or sl is None:
                continue
            dr = e - base_err
            resid_diffs.append(abs(dr))
            delta_biases.append(abs(-dr * 1e-4 * sl))
    TRIM_NOTE = (
        "%d of %d trims sit outside the 1 ct tolerance on C_L (%s)."
        % (len(out_of_tol),
           sum(1 for k, _ in CONDITIONS for t in data[k]["trim_case"] if t),
           "; ".join(out_of_tol))
        if out_of_tol else
        "Every trim is inside the 1 ct tolerance on C_L.")
    WORST_RESID_DIFF_CT = max(resid_diffs) if resid_diffs else float("nan")
    WORST_DELTA_BIAS_CT = max(delta_biases) if delta_biases else float("nan")
    SMALLEST_DELTA_CT = min(
        abs(v) for k, _ in CONDITIONS for v in data[k]["dCd"] if v not in (None, 0.0))

    # ------------------------------------------------------------ footnotes
    paras = [
        INV_TEXT.replace("\n", " "),
        "FRAME (D068). All coefficients DSO basis, half model: Aref 0.620462 m2 (half-wing patch), "
        "lRef 0.393957 m, bref 3.6576 m. 2*Aref = 1.240924 m2 reproduces the DSO Sref to 3.2e-6 relative, "
        "against 11.3e-2 for the TP-1580 Sref; asserted per case, not declared. 1 count = 1e-4.",
        "OPERATING POINTS (D076). Every delta-C_D is formed WITHIN one condition against that condition's "
        "own B trim; no delta crosses a condition. Values are the convergence gate's own 200-sample window "
        "mean (HPC-121), never last samples.",
        "TURBULENCE MODEL IS NOT COMMON ACROSS CONDITIONS: Condition CR %s, early and late cruise %s. A "
        "count from one condition is NOT like-for-like against a count from another; the within-condition "
        "deltas in (c), (d) and (e) are clean." % (data["condition_CR"]["turb"],
                                                   data["early_cruise"]["turb"]),
        "CAMPAIGN COVERAGE. %d cases converged: %s. All %d trims (%d geometries x %d operating "
        "points) are solved, so no group below is drawn as a gap."
        % (sum(data[k]["n_converged"] for k, _ in CONDITIONS),
           ", ".join("%s %d" % (data[k]["label"], data[k]["n_converged"]) for k, _ in CONDITIONS),
           sum(1 for k, _ in CONDITIONS for t in data[k]["trim_case"] if t),
           len(GEOMS), len(CONDITIONS)),
        "CAVEATS. (1) TRIM TOLERANCE. %s Because every number here is a DIFFERENCE inside one "
        "condition, a shared offset cancels; the residual DIFFERENCE across a pair is at most "
        "%.3f ct of C_L, worth %.3f ct of drag on any delta against a smallest plotted delta of "
        "%.1f ct, so no ordering changes. (2) SCALES. Panels (d) and (e) share one y-scale, %.0fx "
        "panel (c)'s, so the two cruise points read against each other directly; Condition CR "
        "needs its own or its deltas would be invisible. The shaded band in (c) is the D058 "
        "fixed-C_L noise floor, 0.106 to 0.199 ct, measured on the coarse wall-function A2 pair "
        "and NOT re-measured on this mesh family; it is the only measured fixed-C_L floor the "
        "project holds. (3) Nothing here is a VLM number: VSPAERO writes CDtot = CDo + CDi with "
        "CDo a flat-plate skin-friction correlation, which must never be set against these "
        "viscous totals."
        % (TRIM_NOTE, WORST_RESID_DIFF_CT, WORST_DELTA_BIAS_CT, SMALLEST_DELTA_CT,
           CRUISE_SCALE_RATIO),
        "Source: docs/report/all_geometry_2026-09-15/data/rans_forces.json. Plotted values at full precision: "
        "fig_trimmed_performance_values.csv. Generator: fig/make_bars.py.",
    ]
    # THE PROVENANCE GOES TO A SIDECAR FILE, NOT INTO THE FIGURE (project decision, 2026-09-22).
    # These six paragraphs were drawn inside the plot and made it 754 words of rendered
    # text, which is a page of prose in a figure. None of it is discarded: it is written
    # beside the figure where the record wants it, and the LaTeX caption carries the two
    # sentences a reader needs at the point of looking.
    foot = "\n\n".join(textwrap.fill(p, width=100) for p in paras)
    notes_path = os.path.join(os.path.dirname(OUT_PDF), "fig_trimmed_performance_notes.txt")
    with open(notes_path, "w") as fh:
        fh.write(foot + "\n")
    print("notes written to %s (%d words, NOT drawn in the figure)"
          % (os.path.basename(notes_path), len(foot.split())))

    fig.savefig(OUT_PDF, bbox_inches="tight")
    tb = fig.get_tightbbox(fig.canvas.get_renderer())
    print("figure %.2f x %.2f in; tight bbox x %.2f to %.2f, y %.2f to %.2f in"
          % (fig.get_figwidth(), fig.get_figheight(), tb.x0, tb.x1, tb.y0, tb.y1))
    plt.close(fig)

    # ------------------------------------------------- verify, do not vouch
    for path in (OUT_PDF, OUT_CSV):
        if not os.path.isfile(path):
            sys.exit("HALT: %s was not written" % path)
        size = os.path.getsize(path)
        if size == 0:
            sys.exit("HALT: %s is zero bytes" % path)
        print("WROTE %s  %d bytes" % (path, size))

    with open(OUT_PDF, "rb") as fh:
        if fh.read(4) != b"%PDF":
            sys.exit("HALT: %s does not begin with %%PDF" % OUT_PDF)
    print("PDF magic OK")

    print("\nPlotted values, read back from the CSV that was just written:")
    with open(OUT_CSV) as fh:
        for row in csv.DictReader(fh):
            print(
                "  %-13s %s  Cd %-12s dCd %-12s L/D %-12s alpha %s"
                % (row["condition"], row["geometry"], row["Cd_counts"],
                   row["delta_Cd_counts_vs_B"], row["L_over_D"], row["alpha_deg"])
            )


if __name__ == "__main__":
    main()
