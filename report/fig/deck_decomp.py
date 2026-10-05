#!/usr/bin/env python3
"""Drag decomposition at Condition CR: the RANS total beside what the VLM writes, C_Dtot =
C_Do + C_Di, its near-field C_Di and its far-field C_Diw, per geometry, in counts. RANS from
rans_forces.json trims; VLM from vlm.json polar_on_DSO_basis through the report's mapping."""
import deck_common as dc
from matplotlib.patches import Patch

dc.style()
t = dc.trims("condition_CR")
SP = dc.jload("fig/spanwise_caption_values.json")["condition_CR"]
V = dc.vlm()


def polar(L):
    r = [x for x in V["runs"] if x["geometry_name"] == SP["vlm_geometry"][L]
         and x["operating_point"]["state"] == "condition_CR"][0]
    return r["polar_on_DSO_basis"]


fig, ax = dc.plt.subplots()
w = 0.2
spec = [("RANS $C_D$ total", lambda L: t[L]["Cd_counts"], 1.0, ""),
        (r"VLM $C_{Dtot} = C_{Do} + C_{Di}$", lambda L: 1e4 * polar(L)["CDtot"], 0.55, ""),
        (r"VLM $C_{Di}$ near (surface)", lambda L: 1e4 * polar(L)["CDi"], 0.4, "//"),
        (r"VLM $C_{Diw}$ far (wake)", lambda L: 1e4 * polar(L)["CDiw"], 0.25, "..")]
for i, L in enumerate(dc.ORDER):
    for j, (lab, f, a, h) in enumerate(spec):
        b = ax.bar(i + (j - 1.5) * w, f(L), width=w * 0.95, color=dc.COLOURS[L], alpha=a,
                   hatch=h, edgecolor=dc.COLOURS[L])
        dc.bar_labels(ax, b, fmt="%.0f", size=11)
ax.set_xticks(range(len(dc.ORDER))); ax.set_xticklabels(dc.ORDER)
ax.set_xlabel("geometry, Condition CR")
ax.set_ylabel("drag  (counts)")
ax.set_ylim(0, 265)
ax.legend([Patch(fc="0.3"), Patch(fc="0.65"), Patch(fc="0.8", ec="0.4", hatch="//"), Patch(fc="0.92", ec="0.4", hatch="..")],
          [s[0] for s in spec], loc="upper left", ncol=2, fontsize=12.5, columnspacing=1.0)
dc.save(fig, "deck_decomp", dc.RESULT_W, dc.RESULT_H)
