#!/usr/bin/env python3
"""Slide: actuator loads on their own lines. H about its flap hinge (x/c 0.70, eta 0.71 to
0.97) and W about the quarter chord over its twist span, each beside the baseline on the same
line and beside the campaign's 0.62 comparison line. Source results/derived/actuator_moments_hw.json."""
import json
import os
import deck_common as dc
from matplotlib.patches import Patch

dc.style()
J = os.path.normpath(os.path.join(dc.R, "..", "..", "..", "results", "derived", "actuator_moments_hw.json"))
d = json.load(open(J))
rows = [("H", "H_flap_0p70", "H, flap hinge x/c 0.70"), ("W", "W_torsion_0p25", "W, torsion axis x/c 0.25")]
fig, ax = dc.plt.subplots()
w = 0.36
for i, (g, line, lab) in enumerate(rows):
    v = d["against_baseline"][line][g]
    b1 = ax.bar(i - w / 2, v["baseline_M_total_Nm"], width=w, color="0.55")
    b2 = ax.bar(i + w / 2, v["M_total_Nm"], width=w, color=dc.COLOURS[g])
    dc.bar_labels(ax, b1, fmt="%.3f", size=13); dc.bar_labels(ax, b2, fmt="%.3f", size=13)
ax.axhline(0, color="black", lw=1.0)
ax.set_xticks(range(len(rows))); ax.set_xticklabels([r[2] for r in rows])
ax.set_ylabel("M about the line  (N m), TE down +")
ax.set_xlabel("geometry and line, Condition CR")
lo = min(min(d["against_baseline"][l][g]["M_total_Nm"], d["against_baseline"][l][g]["baseline_M_total_Nm"]) for g, l, _ in rows)
hi = max(max(d["against_baseline"][l][g]["M_total_Nm"], d["against_baseline"][l][g]["baseline_M_total_Nm"]) for g, l, _ in rows)
ax.set_ylim(min(lo, 0) * 1.25 - 0.1, max(hi, 0) * 1.25 + 0.1)
ax.legend([Patch(fc="0.55"), Patch(fc="0.25")], ["B on the same line", "candidate"], loc="lower left")
dc.save(fig, "deck_actuator", dc.RESULT_W, dc.RESULT_H)
