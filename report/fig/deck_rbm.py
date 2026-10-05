#!/usr/bin/env python3
"""Slide: root bending. VLM and RANS percentage increase per candidate, 6.8 per cent screen."""
import deck_common as dc

dc.style()
rows = {r["letter"]: r for r in dc.cload("fig/rootbending_values.csv")}
fig, ax = dc.plt.subplots()
w = 0.36
for i, L in enumerate(dc.CAND):
    v = float(rows[L]["vlm_increase_pct_recomputed"]); r = float(rows[L]["rans_increase_pct"])
    b1 = ax.bar(i - w / 2, v, width=w, color=dc.COLOURS[L], alpha=0.45, hatch="//",
                edgecolor=dc.COLOURS[L])
    b2 = ax.bar(i + w / 2, r, width=w, color=dc.COLOURS[L])
    dc.bar_labels(ax, b1, fmt="%.2f", size=13); dc.bar_labels(ax, b2, fmt="%.2f", size=13)
ax.axhline(6.8, color="0.25", lw=1.4, ls="--")
ax.set_xticks(range(len(dc.CAND))); ax.set_xticklabels(dc.CAND)
ax.set_ylim(0, 9.0)
ax.set_xlabel("candidate, Condition CR")
ax.set_ylabel("root bending increase against B  (%)")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
ax.legend([Patch(fc="0.75", ec="0.45", hatch="//"), Patch(fc="0.45"),
           Line2D([], [], color="0.25", lw=1.4, ls="--")],
          ["VLM", "RANS", "6.8 % screen"], loc="upper left", ncol=3)
dc.save(fig, "deck_rbm", dc.RESULT_W, dc.RESULT_H)
