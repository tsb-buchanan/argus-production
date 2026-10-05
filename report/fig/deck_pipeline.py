#!/usr/bin/env python3
"""Slide: the automated pipeline, one case from geometry to a quotable number. Boxes and
arrows in TU Delft colours; gates are the grey boxes with the heavy border. Same content as
fig/pipeline_diagram.tex, drawn at slide size."""
import deck_common as dc
from matplotlib.patches import FancyBboxPatch

dc.style()
fig, ax = dc.plt.subplots()
ax.set_xlim(-0.32, 12.1); ax.set_ylim(0, 4.55); ax.axis("off"); ax.grid(False)
FILL = {"box": "#DCF3FA", "gate": "#E4E4E4", "rec": "#D6EBDD"}
EDGE = {"box": dc.TUD_CYAN, "gate": "#444444", "rec": "#2E7D5B"}
LW = {"box": 1.4, "gate": 2.4, "rec": 1.4}
BW, BH = 1.74, 1.02


def box(x, y, text, kind="box"):
    ax.add_patch(FancyBboxPatch((x - BW / 2, y - BH / 2), BW, BH, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=FILL[kind], ec=EDGE[kind], lw=LW[kind]))
    ax.text(x, y, text, ha="center", va="center", fontsize=13.5, linespacing=1.12)


def arrow(x0, y0, x1, y1):
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#333333",
                                                       shrinkA=0, shrinkB=0, mutation_scale=14))


rows = [
    [("OpenVSP geometry\n.vsp3, STL export", "box"),
     ("Scale ft to m,\nstrip non-wing\nsolids", "box"),
     ("surfaceCheck:\nclosed, oriented", "gate"),
     ("snappyHexMesh:\nlevel 11,\nprism layers", "box"),
     ("checkMesh:\nquality limits", "gate"),
     ("decomposePar\n192 subdomains", "box")],
    [("Bracket: three\nangles solved", "box"),
     ("Convergence gate:\ndrift and span", "gate"),
     ("Lift-slope fit\nto target $C_L$", "box"),
     ("Trim case at\nthe fitted angle", "box"),
     ("Gate again,\nlift on target", "gate"),
     ("Run card, JSON:\nchecksum, mesh,\nforces", "rec")],
    [("Harvest to\nrans_forces.json", "rec"),
     ("Tables and figures\nregenerated", "rec"),
     ("Report gates:\ntables match data", "gate")],
]
ys = [3.85, 2.28, 0.71]
xs = [1.05 + i * 2.0 for i in range(6)]
# row 1 left to right, row 2 right to left, row 3 left to right (a snake, like the TikZ)
for i, (t, k) in enumerate(rows[0]):
    box(xs[i], ys[0], t, k)
    if i:
        arrow(xs[i - 1] + BW / 2, ys[0], xs[i] - BW / 2, ys[0])
for i, (t, k) in enumerate(rows[1]):
    x = xs[5 - i]
    box(x, ys[1], t, k)
    if i:
        arrow(xs[5 - i + 1] - BW / 2, ys[1], x + BW / 2, ys[1])
for i, (t, k) in enumerate(rows[2]):
    box(xs[i], ys[2], t, k)
    if i:
        arrow(xs[i - 1] + BW / 2, ys[2], xs[i] - BW / 2, ys[2])
arrow(xs[5], ys[0] - BH / 2, xs[5], ys[1] + BH / 2)
arrow(xs[0], ys[1] - BH / 2, xs[0], ys[2] + BH / 2)
# side labels, the same three the report diagram carries
for y, t in zip(ys, ["geometry", "solve and trim", "record"]):
    ax.text(-0.2, y, t, ha="center", va="center", fontsize=12.5, style="italic", color="#555555",
            rotation=90)
dc.save(fig, "deck_pipeline", dc.WIDE_W, dc.WIDE_H)
