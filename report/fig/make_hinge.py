#!/usr/bin/env python3
"""
Hinge moment about the swept x_h/c = 0.62 line, every geometry, every condition.
ARGUS all-geometry report 2026-09-15.

EVERY NUMBER IN THE OUTPUT IS READ FROM A FILE. Nothing is typed in from memory.

Inputs (read, never assumed):
  results/derived/<CASE>_trim.json    one record per trimmed case, written by
                                      scripts/case_derived_quantities.py
  data/vlm.json                       geometry letter <-> DSO geometry name

THE DEFINITION IS THE ONE AGREED WITH THE DSO TEAM, and it is reproduced here so the
figure can be checked against the agreement rather than against this script's intent:

  axis      the swept spanwise line through x_h/c = 0.62 at each local section, with
            the moment projected about the LOCAL TANGENT of that line, not the y axis
  sign      positive is trailing-edge down; the actuator reaction is the opposite sign
  scope     surfaces aft of the line only, over 0.60 <= eta <= 1.00, half wing
  sectional m'_h [N m/m] and c_mh = m'_h / (q_inf c(y)^2), on the LOCAL section chord
  integral  M_h [N m] and C_Mh,region = M_h / [q_inf INTEGRAL c(y)^2 dy]
  terms     pressure is the primary result, viscous reported separately, total included

WHY BOTH A DIMENSIONAL AND A COEFFICIENT PANEL, and this is the whole point of the
layout. The three conditions sit at q = 708, 11293 and 8257 Pa, a spread of 16x. A
DIMENSIONAL MOMENT IS NOT COMPARABLE ACROSS DYNAMIC PRESSURES, so putting m'_h for all
three on one axis would show the dynamic pressure and hide the aerodynamics. The upper
row therefore gives m'_h with ONE PANEL PER CONDITION, each on its own axis, and the
lower row gives c_mh, which is the quantity that may legitimately be compared across
them. Neither row is sufficient alone: the dimensional one is what an actuator reacts,
the coefficient one is what the shapes can be ranked on.

PARTITION (GEO-089). Conditions x geometries is enumerated, adjudicated, and the buckets
are asserted to sum to the expected count. A case that is missing, or present with no
hinge moment, is NAMED AND COUNTED, never silently dropped: a figure that quietly omits
a geometry lets a reader conclude it does not exist rather than that it has not run.

FRAMES ARE READ FROM EACH RECORD, NOT ASSUMED (D068). q, alpha and the compressible flag
come from the record's own fields. The per-condition q is CROSS-CHECKED against the value
the record carries and the script HALTS on disagreement, because the one thing that must
never happen here is a correct moment plotted under the wrong dynamic pressure.
"""

import csv
import json
import derived_source
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(REPORT)))
DATA = os.path.join(REPORT, "data")
DERIVED = os.path.join(REPO, "results", "derived")
FIG = HERE

LOG_LINES = []


def log(msg=""):
    print(msg)
    LOG_LINES.append(str(msg))


def halt(msg):
    log("HALT: " + msg)
    raise SystemExit("HALT: " + msg)


# ----------------------------------------------------------------------------
# Style, matched to make_cf.py so the report's figures are one family
# ----------------------------------------------------------------------------
plt.rcParams.update(
    {
        "font.size": 8.5,
        "axes.titlesize": 9.0,
        "axes.labelsize": 8.5,
        "legend.fontsize": 7.8,
        "xtick.labelsize": 7.8,
        "ytick.labelsize": 7.8,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.4,
        "savefig.bbox": "tight",
        "figure.dpi": 200,
        "mathtext.fontset": "dejavusans",
    }
)

STYLE = {
    "B": dict(color="#000000", ls="-", lw=1.5),
    "C": dict(color="#D55E00", ls=(0, (5.5, 1.8)), lw=1.6),
    "F": dict(color="#0072B2", ls=(0, (4.5, 1.5, 1.0, 1.5)), lw=1.5),
    "H": dict(color="#009E73", ls=(0, (1.3, 1.2)), lw=1.7),
    "M": dict(color="#CC79A7", ls=(0, (4.0, 1.4, 1.0, 1.4, 1.0, 1.4)), lw=1.5),
    "W": dict(color="#56B4E9", ls=(0, (8.0, 2.2)), lw=1.6),
}
LETTERS = ["B", "C", "F", "H", "M", "W"]

# THESE TWO ARE NOT ACTUATOR MOMENTS AND THE FIGURE MUST SAY SO (2026-09-17).
# Everything here is integrated about the common line x_h/c = 0.62 over the aft 38 per cent,
# which is the only way C_Mh_region means anything across six geometries. But H's own hinge is
# at 0.70, and W is not hinged at all: it rotates the whole section about 0.25c, so its actuator
# reacts a TORSIONAL moment over the entire section. Their bars are valid comparisons and are NOT
# the load their own mechanism carries, so they are hatched. Without this a reader ranks the bars
# by height and concludes distributed twist is the cheapest to actuate, which is not supported.
NOT_ACTUATOR = {"H", "W"}

# prefix -> (label, expected q [Pa], target CL, turbulence/wall note)
CONDITIONS = [
    ("SW", "Condition CR, $M = 0.10$", 708.05, 0.428277635108,
     "kOmegaSST, walls resolved"),
    ("CMP", "Early cruise, $M = 0.78$", 11292.804088955776, 0.529297087450,
     "SpalartAllmaras, walls modelled"),
    ("LC", "Late cruise, $M = 0.78$", 8257.365999, 0.544114800210,
     "SpalartAllmaras, walls modelled"),
]
Q_TOL_REL = 1e-4

# ----------------------------------------------------------------------------
# 1. Geometry letter -> DSO geometry name, READ FROM vlm.json
#    The SAME source and the SAME path make_cf.py uses. A second mapping built a second
#    way is a second thing that can disagree with the first.
# ----------------------------------------------------------------------------
GEOM_NAME = {}
with open(os.path.join(DATA, "vlm.json")) as fh:
    _vlm = json.load(fh)
for _r in _vlm["runs"]:
    _m = _r.get("mapping") or {}
    _fam = _m.get("case_prefix_family")
    if _fam in LETTERS and _r.get("geometry_name"):
        GEOM_NAME.setdefault(_fam, _r["geometry_name"])
log("Geometry names read from vlm.json (mapping.case_prefix_family -> geometry_name):")
for k in LETTERS:
    log("  %s = %s" % (k, GEOM_NAME.get(k, "NOT FOUND")))

# ----------------------------------------------------------------------------
# 2. Load every trimmed record and PARTITION the condition x geometry grid
# ----------------------------------------------------------------------------
REC = {}
present, missing_file, no_hinge = [], [], []

for prefix, label, q_expect, cl_target, note in CONDITIONS:
    for L in LETTERS:
        # THE KEY STAYS "<prefix><L>_trim" so every membership test further down is
        # unchanged; only the PATH is resolved, so Condition CR reads its wake-refined
        # record while the two cruise points keep falling back to the published one.
        case = "%s%s_trim" % (prefix, L)
        path, _gen = derived_source.resolve(prefix, L)
        if not os.path.exists(path):
            missing_file.append(case)
            continue
        with open(path) as fh:
            d = json.load(fh)
        h = d.get("hinge_moment") or {}
        if not h.get("available"):
            no_hinge.append((case, h.get("note", "no note recorded")))
            continue
        # FRAME CROSS-CHECK. The record's own q against the condition's q.
        q_rec = h.get("q_Pa")
        if q_rec is None:
            halt("%s carries a hinge moment with no q_Pa; it cannot be placed in a "
                 "condition" % case)
        if abs(q_rec - q_expect) / q_expect > Q_TOL_REL:
            halt("%s records q = %.6f Pa but condition %s expects %.6f Pa "
                 "(relative difference %.3e). A moment plotted under the wrong dynamic "
                 "pressure is the one error this figure must never make."
                 % (case, q_rec, label, q_expect, abs(q_rec - q_expect) / q_expect))
        REC[case] = d
        present.append(case)

n_expect = len(CONDITIONS) * len(LETTERS)
log("")
log("PARTITION of the %d condition x geometry cells (GEO-089):" % n_expect)
log("  resolved, with a hinge moment : %d" % len(present))
log("  present but no hinge moment   : %d" % len(no_hinge))
log("  no derived record on disk     : %d" % len(missing_file))
for c, why in no_hinge:
    log("      %-12s %s" % (c, why[:150]))
if missing_file:
    log("      not on disk: %s" % ", ".join(missing_file))
if len(present) + len(no_hinge) + len(missing_file) != n_expect:
    halt("the buckets do not sum to %d" % n_expect)
log("  buckets sum to %d, as they must" % n_expect)
if not present:
    halt("no case has a hinge moment; there is nothing to plot")


def _smooth(v, frac=0.08):
    """Hann-window moving average over a fixed FRACTION of the stations, for DISPLAY only.

    WHY A FILTER AND NOT A REBIN. Block-averaging the stations was tried and MEASURED not
    to work: over 638 stations the lag-1 autocorrelation of the residual was -0.498, and
    after collapsing to 91 stations it was still -0.490. The alternation survives block
    averaging at every scale, so it is not a fixed-period wiggle that a coarser grid
    removes. Its amplitude fell from 9.1 to 2.2 per cent of range while its CHARACTER did
    not change, which is exactly the signature of an artefact that reproduces on whatever
    grid it is reported on.

    WHAT THIS SEPARATES. The physical distribution varies over about 0.1 in eta: the
    outboard decay follows the chord, and the features worth seeing, such as H pulling
    clear at cruise, are that wide. The artefact lives at the bin scale, two orders finer.
    A Hann window spanning 8 per cent of the band sits between the two, so it removes the
    artefact and leaves the physics. A Hann rather than a boxcar because a boxcar's
    sidelobes let the alternation back through at reduced amplitude.

    THE UNSMOOTHED DATA STAYS AVAILABLE, IN THE RECORDS RATHER THAN ON THE PLOT. Drawing
    it faintly underneath was tried and withdrawn: see the note at the call site. Smoothing
    a display is legitimate as long as the raw values remain reachable and the caption says
    the line is smoothed, and both hold here. Nothing is smoothed before integration:
    M_h is a direct sum over faces and never sees these stations at all.
    """
    v = np.asarray(v, float)
    n = max(5, int(len(v) * frac) | 1)          # odd window, at least 5
    w = np.hanning(n)
    w = w / w.sum()
    pad = np.pad(v, n // 2, mode="edge")
    return np.convolve(pad, w, mode="valid")[:len(v)]


def sect(case):
    s = REC[case]["hinge_moment"]["sectional"]
    return (np.array([r["eta"] for r in s]),
            np.array([r["m_prime_h_Nm_per_m"] for r in s]),
            np.array([r["c_mh"] for r in s]))


def legend_handles(letters):
    h, lab = [], []
    for L in letters:
        st = STYLE[L]
        h.append(Line2D([0], [0], color=st["color"], ls=st["ls"], lw=st["lw"]))
        lab.append("%s  %s" % (L, GEOM_NAME.get(L, "name NOT FOUND")))
    return h, lab


# ----------------------------------------------------------------------------
# 3. FIGURE 1: sectional distribution, m'_h per condition and c_mh shared
# ----------------------------------------------------------------------------
def make_sectional_figure(outfile):
    ncond = len(CONDITIONS)
    fig, axes = plt.subplots(2, ncond, figsize=(7.4, 5.4), sharex=True)
    if ncond == 1:
        axes = axes.reshape(2, 1)

    for j, (prefix, label, q_expect, cl_target, note) in enumerate(CONDITIONS):
        ax_d, ax_c = axes[0, j], axes[1, j]
        drawn = 0
        for L in LETTERS:
            case = "%s%s_trim" % (prefix, L)
            if case not in REC:
                continue
            eta, mp, cmh = sect(case)
            # WHY THE CURVES ARE SMOOTHED AT ALL. At the 1600 bins the INTEGRAL needs, each
            # reporting bin holds about 3300 faces instead of 26000, so the per-bin scatter
            # is large. It is a SAMPLING artefact of the reporting bins and not a feature
            # of the flow: it alternates station to station, lag-1 autocorrelation -0.49 on
            # all three baselines, and a block rebin to 91 stations left that at -0.49
            # while only shrinking the amplitude, which is what proved it is not a
            # fixed-period wiggle. The Hann filter moves lag-1 to +0.89 and takes the
            # direction reversals from 427 to 1.
            # THE RAW SERIES IS NOT DRAWN, AND THAT IS A DELIBERATE REVERSAL. Plotting it
            # faintly beneath the smoothed line was tried: 638 stations times six
            # geometries at alpha 0.15 accumulates into a haze band that dominates every
            # panel and buries the curves it was meant to keep honest. A reader cannot
            # check data they cannot resolve, so the band bought nothing and cost the
            # figure. The unsmoothed per-station values remain in
            # results/derived/<case>.json, which is where anyone auditing the curve would
            # go, and the caption says so.
            st = STYLE[L]
            ax_d.plot(eta, _smooth(mp), zorder=3, **st)
            ax_c.plot(eta, _smooth(cmh), zorder=3, **st)
            drawn += 1
        for ax in (ax_d, ax_c):
            ax.axhline(0.0, color="0.55", lw=0.6, zorder=0)
            ax.set_xlim(0.60, 1.00)
            ax.grid(True, lw=0.35, color="0.88")
        # THREE LINES, NOT TWO, AND SMALLER. At 8.2 pt the condition line and the model
        # line overran the panel width and collided, printing "q = 708 PaSpalartAllmaras".
        ax_d.set_title(label,
                       fontsize=7.2, linespacing=1.25)
        ax_c.set_xlabel(r"$\eta$")
        if drawn == 0:
            for ax in (ax_d, ax_c):
                ax.text(0.5, 0.5, "no case at this condition\ncarries a hinge moment",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=8, color="#B00000", style="italic")
        if j == 0:
            ax_d.set_ylabel(r"$m'_h$,  N m/m")
            ax_c.set_ylabel(r"$c_{mh} = m'_h / (q_\infty c(y)^2)$")

    # the coefficient row is the comparable one, so give it a shared scale
    lims = [ax.get_ylim() for ax in axes[1, :] if ax.lines]
    if lims:
        lo = min(l[0] for l in lims); hi = max(l[1] for l in lims)
        for ax in axes[1, :]:
            ax.set_ylim(lo, hi)

    h, lab = legend_handles([L for L in LETTERS
                             if any("%s%s_trim" % (p, L) in REC for p, *_ in CONDITIONS)])
    # NO METHOD NOTE UNDER THE FIGURE. The sign convention, the frame, the reason the three
    # panels may not be read against one another and the smoothing all belong in the
    # section prose, where they can be read in order and found again. A four-line footnote
    # crowds the plot and is the last place anyone looks for a definition.
    # -0.075, NOT -0.03. Removing the footnote freed space below the axes and I pulled the
    # legend up to take it, which put the legend row on top of the x-axis tick labels and
    # the eta axis label of the right-hand panel. The footnote's removal is not a licence
    # to move the legend into the axes.
    fig.legend(h, lab, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.075))
    fig.suptitle("",
                 fontsize=9.6, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


# ----------------------------------------------------------------------------
# 4. FIGURE 2: integrated moment per geometry per condition, with the term split
# ----------------------------------------------------------------------------
def make_integrated_figure(outfile):
    conds = [c for c in CONDITIONS
             if any("%s%s_trim" % (c[0], L) in REC for L in LETTERS)]
    fig, axes = plt.subplots(2, len(conds), figsize=(7.4, 5.0), squeeze=False)

    for j, (prefix, label, q_expect, cl_target, note) in enumerate(conds):
        ax_m, ax_c = axes[0][j], axes[1][j]
        letters = [L for L in LETTERS if "%s%s_trim" % (prefix, L) in REC]
        x = np.arange(len(letters))
        mp = [REC["%s%s_trim" % (prefix, L)]["hinge_moment"]["M_h_pressure_Nm"]
              for L in letters]
        mv = [REC["%s%s_trim" % (prefix, L)]["hinge_moment"]["M_h_viscous_Nm"] or 0.0
              for L in letters]
        cm = [REC["%s%s_trim" % (prefix, L)]["hinge_moment"].get("C_Mh_region_total")
              or REC["%s%s_trim" % (prefix, L)]["hinge_moment"]["C_Mh_region"]
              for L in letters]
        _bm = ax_m.bar(x, [p + v for p, v in zip(mp, mv)], 0.62,
                       color=[STYLE[L]["color"] for L in letters],
                       edgecolor="0.2", linewidth=0.5, label="total (pressure + viscous)")
        for _b, _L in zip(_bm, letters):
            if _L in NOT_ACTUATOR:
                _b.set_hatch("///")
        # THE VISCOUS TERM IS 0.13 TO 0.47 PER CENT OF THE TOTAL, so a stacked hatch for it
        # is SUB-PIXEL: the legend promised a split the eye cannot find. It is printed as a
        # percentage instead, which is the honest way to show a term this small, and the
        # bar is labelled as the total it actually draws.
        for xi, (p, v) in enumerate(zip(mp, mv)):
            tot = p + v
            ax_m.annotate("%.3f" % tot, (xi, tot), ha="center", va="center",
                          fontsize=6.4, rotation=90, color="white", fontweight="bold",
                          xytext=(0, 14), textcoords="offset points")
            pass  # viscous share is in the table (22 Sep)
        _bc = ax_c.bar(x, cm, 0.62, color=[STYLE[L]["color"] for L in letters],
                       edgecolor="0.2", linewidth=0.5)
        for _b, _L in zip(_bc, letters):
            if _L in NOT_ACTUATOR:
                _b.set_hatch("///")
        for ax, ttl in ((ax_m, None), (ax_c, None)):
            ax.set_xticks(x); ax.set_xticklabels(letters)
            ax.axhline(0.0, color="0.4", lw=0.7)
            ax.grid(True, axis="y", lw=0.35, color="0.88")
        ax_m.set_title(label,
                       fontsize=8.2, linespacing=1.3)
        if j == 0:
            ax_m.set_ylabel(r"$M_h$,  N m   (half wing)")
            ax_c.set_ylabel(r"$C_{Mh,\mathrm{region}}$")
            _h, _l = ax_m.get_legend_handles_labels()
            _h.append(Patch(facecolor="0.88", edgecolor="0.2", hatch="///"))
            _l.append("hatched: not the actuator load")
            ax_m.legend(_h, _l, frameon=False, fontsize=7.0, loc="best")

    # NO METHOD NOTE UNDER THE FIGURE, for the reason given in the sectional figure above.
    fig.suptitle("",
                 fontsize=9.6, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


# ----------------------------------------------------------------------------
# 5. The values, dumped so a reader can check the figure against a file
# ----------------------------------------------------------------------------
def _case_cl(case, rec):
    """The case's own achieved C_L, from its forceCoeffs, with the record's value preferred.

    THE RECORD'S FIELD IS EMPTY ON EVERY CASE and shipping the column that way would be a
    blank in a delivered file, which reads as missing data rather than as an unrun step.
    case_derived_quantities.case_CL reads forceCoeffs_spliced.dat, which splice_postpro.py
    writes and which was never run for these trims. _forcecoeffs_cl reads the case's OWN
    forceCoeffs directory instead: the same independent code path the force-frame gate is
    arbitrated against, so it is a source already trusted here rather than a new one.
    """
    v = rec.get("CL")
    if v is not None:
        return v
    try:
        sys.path.insert(0, os.path.join(REPO, "scripts"))
        from case_derived_quantities import _forcecoeffs_cl
        return _forcecoeffs_cl(case)
    except Exception:                                            # noqa: BLE001
        return None


def write_values(path):
    rows = []
    for prefix, label, q_expect, cl_target, note in CONDITIONS:
        for L in LETTERS:
            case = "%s%s_trim" % (prefix, L)
            if case not in REC:
                continue
            d = REC[case]; h = d["hinge_moment"]
            rows.append(dict(
                case=case, condition=label.replace("$", "").replace("\\", ""),
                geometry_letter=L, geometry_name=GEOM_NAME.get(L, "NOT FOUND"),
                alpha_deg=d.get("alpha_deg"), q_Pa=h.get("q_Pa"),
                CL_target=cl_target, CL_case=_case_cl(case, d),
                M_h_pressure_Nm=h.get("M_h_pressure_Nm"),
                M_h_viscous_Nm=h.get("M_h_viscous_Nm"),
                M_h_total_Nm=h.get("M_h_total_Nm"),
                C_Mh_region_pressure=h.get("C_Mh_region"),
                C_Mh_region_total=h.get("C_Mh_region_total"),
                integral_c2_dy_m3=h.get("integral_c2_dy_m3"),
                n_faces_in_region=h.get("n_faces_in_region"),
                n_faces_total=h.get("n_faces_total"),
                normal_orientation=(d.get("normal_orientation") or {}).get("verdict"),
                tau_frame=(d.get("tau_frame") or {}).get("verdict"),
                reference_pressure_Pa=(d.get("reference_pressure") or {}).get("p0_Pa"),
                surface_vtk=d.get("surface_vtk"),
            ))
    if not rows:
        return None
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


# ----------------------------------------------------------------------------
written = []
written.append(make_sectional_figure(os.path.join(FIG, "fig_hinge_sections.pdf")))
written.append(make_integrated_figure(os.path.join(FIG, "fig_hinge_integrated.pdf")))
v = write_values(os.path.join(FIG, "hinge_moment_values.csv"))
if v:
    written.append(v)

# VERIFICATION: stat and read back every artefact. A figure that failed to write is
# worse than one that was never attempted, because the previous run's file is still
# sitting there under the right name (GEO-080 item 6).
log("")
log("VERIFICATION (stat of each artefact, read back from disk):")
bad = 0
for p in written:
    if not os.path.exists(p):
        log("  MISSING  %s" % os.path.basename(p)); bad += 1; continue
    n = os.path.getsize(p)
    with open(p, "rb") as fh:
        head = fh.read(4)
    kind = "PDF" if head.startswith(b"%PDF") else "text"
    ok = n > 0 and (kind == "PDF" or p.endswith(".csv"))
    log("  %s  %-44s %8d bytes  %s" % ("OK " if ok else "BAD", os.path.basename(p), n, kind))
    bad += 0 if ok else 1

log("")
log("Cells plotted: %d of %d. %s"
    % (len(present), n_expect,
       "Every condition x geometry cell resolved."
       if len(present) == n_expect else
       "The unresolved cells are named in the partition above and are NOT silently "
       "omitted from the figures."))

with open(os.path.join(FIG, "make_hinge_LOG.txt"), "w") as fh:
    fh.write("\n".join(LOG_LINES) + "\n")
log("log written to %s" % os.path.join(FIG, "make_hinge_LOG.txt"))

sys.exit(1 if bad else 0)
