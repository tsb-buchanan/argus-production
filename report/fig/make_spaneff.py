#!/usr/bin/env python3
"""
Trefftz-plane induced drag and span efficiency at Condition CR.
ARGUS all-geometry report 2026-09-15.

EVERY NUMBER IN THE OUTPUT IS READ FROM A FILE. Nothing is typed in from memory.

Input (read, never assumed):
  docs/report/.../data/span_efficiency.json   written by scripts/build_span_efficiency.py,
                                              which refuses to emit a file whose geometries
                                              do not share one mesh generation and one
                                              trimmed operating point.

-------------------------------------------------------------------------------
WHY A DECAY CURVE AND NOT A BAR CHART OF ONE NUMBER
-------------------------------------------------------------------------------
Crossflow kinetic energy is conserved in an inviscid wake, so the Trefftz integral has NO
PLATEAU to quote: it decays downstream at a rate set by how fast the grid dissipates the
trailing vortex system. A single station would be a number without the thing that
qualifies it. Panel (a) therefore plots C_Di against station and panel (b) the span
efficiency the same integral implies, and the SPREAD across stations is the uncertainty on
C_Di rather than scatter to be averaged away.

-------------------------------------------------------------------------------
THE LINE AT e = 1 IS NOT A TOLERANCE
-------------------------------------------------------------------------------
e = CL^2 / (pi * AR * C_Di) cannot exceed 1 for a planar wake. It is a self-certifying
gate: a curve above it is reporting that its own wake is unresolved, not reporting a
result. The data file refuses to carry such a case at all, so the line is drawn here as
the reader's check on that refusal rather than as a limit anyone is near.

FRAME, drawn on the figure rather than left to the caption (D065/D068):
  C_Di on the HALF wing, Aref 0.620462 m2, matching each case's own forceCoeffs.
  e on the DSO-basis full-wing aspect ratio, bref 3.6576 m, Sref 1.24092 m2.
  Condition CR, U_inf 34.0 m/s, every geometry trimmed to CL = 0.428277635108.
"""

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)
DATA = os.path.join(REPORT, "data")

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def halt(m):
    log("HALT: " + m)
    raise SystemExit("HALT: " + m)


plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 9.0, "axes.labelsize": 8.5,
    "legend.fontsize": 7.8, "xtick.labelsize": 7.8, "ytick.labelsize": 7.8,
    "axes.linewidth": 0.7, "savefig.bbox": "tight", "figure.dpi": 200,
    "mathtext.fontset": "dejavusans",
})

STYLE = {
    "B": dict(color="#000000"), "C": dict(color="#D55E00"),
    "F": dict(color="#0072B2"), "H": dict(color="#009E73"),
    "M": dict(color="#CC79A7"), "W": dict(color="#56B4E9"),
}
LETTERS = ["B", "C", "F", "H", "M", "W"]
# Letter only; the descriptions belong in the report text.
LABEL = {g: g for g in ["B", "C", "F", "H", "M", "W"]}


def main():
    path = os.path.join(DATA, "span_efficiency.json")
    if not os.path.exists(path):
        halt("no %s. Run scripts/build_span_efficiency.py first; it refuses to write one "
             "until every geometry has a trimmed evaluation on one mesh generation." % path)
    d = json.load(open(path))
    g = d["geometries"]

    have = [L for L in LETTERS if L in g]
    absent = [L for L in LETTERS if L not in g]
    log("geometries in the data file: %d of %d (%s)" % (len(have), len(LETTERS),
                                                        ", ".join(have)))
    # A DECLARED SUBSET IS ACCEPTABLE; AN UNDECLARED ONE IS NOT. build_span_efficiency.py
    # writes a `subset` block ONLY when a subset was requested with a recorded reason, so
    # the presence of that block is itself the evidence that the narrowing was deliberate.
    # Without it, missing geometries mean the data file is incomplete by accident and the
    # figure must not be drawn (GEO-087).
    sub = d.get("subset")
    if absent and not sub:
        halt("the data file is missing %s and carries no declared subset. This figure is "
             "captioned as the whole field, so it will not be drawn over an undeclared "
             "subset (GEO-087). Re-run build_span_efficiency.py with --allow-subset and "
             "--reason if the narrowing is intended." % ", ".join(absent))
    if absent and sub:
        declared = set(sub.get("excluded", []))
        if set(absent) != declared:
            halt("the data file declares %s excluded but %s are actually missing. The "
                 "declaration and the content disagree."
                 % (sorted(declared), sorted(absent)))
        log("DECLARED SUBSET: %d of %d geometries, excluding %s"
            % (len(have), len(LETTERS), ", ".join(sorted(declared))))
        log("  reason: %s" % sub.get("reason", "(none recorded)"))

    fr = d["frame"]
    log("frame: %s mesh, Aref %.6f m2, AR %.4f, %s, CL target %.12f"
        % (fr["mesh_generation"], fr["Aref_m2"], fr["AR_full"], fr["condition"],
           fr["CL_target"]))

    # The data file already asserted one generation and one station set. Re-assert here
    # rather than trust it: this script is the last thing between the numbers and the
    # report, and a guard that reads a file inherits that file's correctness (GEO-083).
    gens = sorted({g[L]["mesh_generation"] for L in have})
    if len(gens) != 1:
        halt("the data file carries %d mesh generations, %s" % (len(gens), gens))
    stations = sorted({tuple(round(s["x_m"], 3) for s in g[L]["stations"]) for L in have})
    if len(stations) != 1:
        halt("the data file carries %d station sets, %s" % (len(stations), stations))
    xs = list(stations[0])
    log("stations x = %s m" % ", ".join("%.3f" % x for x in xs))

    # EVERY GEOMETRY IS DRAWN, WITH ITS VALUE, IN ONE STYLE. Project decision, 2026-09-23. The
    # figure is captioned as the whole field, and a tick with no curve on it reads as a
    # missing result rather than as a sourcing note. The guards above are NOT weakened:
    # they still run over the current-evaluation set. This only adds the entries the data
    # file declares under published_evaluation, and asserts their preconditions here.
    #
    # THE MIX IS BOUNDED, NOT ASSUMED AWAY (D071's caveat: cancellation is exact only for a
    # linear differential). Measured on the two geometries evaluated BOTH ways, the mesh
    # generation shifts e by 0.0024 (M: 0.9698 published, 0.9722 wake-refined) and 0.0019
    # (B: 0.9654, 0.9673). Nothing in the field ordering moves. The geometries where the
    # generation DID matter, W H and F at -0.30 in e, are all on the current evaluation.
    # The per-row generation stays in fig_spaneff_values.csv, which is where the frame
    # belongs (D068 requires it to travel with the number in the RECORD).
    pub = d.get("published_evaluation", {}).get("geometries", {})
    gall = dict(g)
    extra = []
    for L in absent:
        if L not in pub:
            halt("%s has no evaluation in either the current or the published block, so "
                 "this figure cannot carry the field it is captioned as." % L)
        # THE EXEMPTION ASSERTS ITS OWN PRECONDITION (GEO-092). An untrimmed leg carries a
        # uniform freestream-w residual across the whole plane, because span_efficiency
        # removes w = U*sin(alpha) with the alpha it was handed. Falling through to the
        # draw loop because the trimmed flag merely went unchecked is the failure shape.
        if not pub[L].get("trimmed"):
            halt("%s's published entry is not a trimmed leg (trimmed=%r). An untrimmed "
                 "evaluation must not be drawn beside trimmed ones."
                 % (L, pub[L].get("trimmed")))
        if not pub[L].get("stations"):
            halt("%s's published entry carries no stations" % L)
        gall[L] = pub[L]
        extra.append(L)
    drawn = [L for L in LETTERS if L in gall]
    if extra:
        log("also drawn, from the declared published evaluation: %s" % ", ".join(extra))
    # PARTITION CLOSES (GEO-089): current + published must be every letter, asserted.
    if len(drawn) != len(LETTERS):
        halt("drawing %d of %d geometries (%s); the caption claims the whole field"
             % (len(drawn), len(LETTERS), ", ".join(drawn)))

    fig, axes = plt.subplots(1, 2, figsize=(6.69, 3.1))
    ax_cdi, ax_e = axes

    for L in drawn:
        e = gall[L]
        st = sorted(e["stations"], key=lambda s: s["x_m"])
        x = [s["x_m"] for s in st]
        ax_cdi.plot(x, [s["CDi_counts"] for s in st], "-o", ms=3.0, lw=1.2,
                    label=LABEL[L], **STYLE[L])
        ax_e.plot(x, [s["span_efficiency"] for s in st], "-o", ms=3.0, lw=1.2,
                  label=LABEL[L], **STYLE[L])
        log("  %s  CDi %7.3f -> %7.3f ct over %.3f m (decay %5.2f%%),  e_last %.4f,  "
            "CL %+6.3f ct of target"
            % (L, e["CDi_counts_first"], e["CDi_counts_last"], x[-1] - x[0],
               e["decay_first_to_last_percent"], e["e_last"],
               e["CL_minus_target_counts"]))

    ax_cdi.set_xlabel(r"Trefftz plane station $x$ [m]")
    ax_cdi.set_ylabel(r"$C_{Di}$ [counts]")
    ax_cdi.set_title("(a)")

    ax_e.axhline(1.0, color="#7a1616", lw=0.8, ls="--", zorder=1, label="$e = 1$")
    ax_e.set_xlabel(r"Trefftz plane station $x$ [m]")
    ax_e.set_ylabel(r"span efficiency $e$")
    ax_e.set_title("(b)")
    ax_e.legend(loc="lower right", frameon=False)
    # Keep e = 1 visible even though nothing approaches it: the reader checks the gate by
    # seeing the distance to it, and an axis autoscaled to the data alone hides that.
    lo = min(s["span_efficiency"] for L in drawn for s in gall[L]["stations"])
    ax_e.set_ylim(min(lo - 0.02, 0.86), 1.02)

    for ax in axes:
        ax.grid(True, lw=0.4, color="#e8e7e3", zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    pass  # the e = 1 line is a legend entry now (22 Sep)

    ax_cdi.legend(loc="upper right", frameon=False, handlelength=1.6,
                  borderaxespad=0.2, labelspacing=0.3)

    fig.text(0.5, -0.045,
             "",
             ha="center", fontsize=6.8, color="#52514e")

    fig.tight_layout()
    out = os.path.join(HERE, "fig_spaneff.pdf")
    fig.savefig(out)
    log("wrote %s" % out)

    # The values behind the figure, so the prose quotes a file rather than the plot.
    csv = os.path.join(HERE, "fig_spaneff_values.csv")
    with open(csv, "w") as fh:
        fh.write("id,case,mesh_generation,CDi_counts_first,CDi_counts_last,"
                 "decay_percent,e_last,CL_minus_target_counts\n")
        for L in drawn:
            e = gall[L]
            fh.write("%s,%s,%s,%.4f,%.4f,%.3f,%.5f,%.4f\n"
                     % (L, e["case"], e["mesh_generation"], e["CDi_counts_first"],
                        e["CDi_counts_last"], e["decay_first_to_last_percent"],
                        e["e_last"], e["CL_minus_target_counts"]))
    log("wrote %s" % csv)

    with open(os.path.join(HERE, "make_spaneff_LOG.txt"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
