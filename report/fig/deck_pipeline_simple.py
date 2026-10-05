#!/usr/bin/env python3
"""The pipeline end to end in eight boxes, one row, TU Delft colours."""
import deck_common as dc
from matplotlib.patches import FancyBboxPatch

dc.style()
steps = ["OpenVSP\ngeometry", "Geometry\nprocessing", "OpenFOAM\ncase setup", "HPC\nsubmission",
         "Mesh", "Solve\nand trim", "Post-\nprocessing", "Automated\nreport"]
W, H = dc.WIDE_W, 1.6
fig, ax = dc.plt.subplots()
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off"); ax.grid(False)
bw, bh = 1.28, 1.0
gap = (W - 8 * bw) / 9.0
for i, t in enumerate(steps):
    x = gap + i * (bw + gap)
    ax.add_patch(FancyBboxPatch((x, (H - bh) / 2), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.1",
                                fc="#DCF3FA", ec=dc.TUD_CYAN, lw=1.8))
    ax.text(x + bw / 2, H / 2, t, ha="center", va="center", fontsize=15)
    if i:
        ax.annotate("", (x, H / 2), (x - gap, H / 2),
                    arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#333333", shrinkA=0, shrinkB=0, mutation_scale=14))
dc.save(fig, "deck_pipeline_simple", W, H)
