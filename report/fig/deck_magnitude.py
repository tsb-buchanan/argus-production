#!/usr/bin/env python3
"""Slide: magnitude. RANS total-drag increment beside both VLM induced-drag increments per
candidate, one panel, from fig_vlmerr_values.json 'errors_vs_mapped_B'."""
import deck_common as dc

dc.style()
e = dc.vlmerr()["errors_vs_mapped_B"]
fig, ax = dc.plt.subplots()
w = 0.26
spec = [("rans_dCd_total_ct", 1.0, "", "RANS total"),
        ("vlm_dCDi_near_ct", 0.45, "//", "VLM near field (surface)"),
        ("vlm_dCDiw_far_ct", 0.25, "..", "VLM far field (wake)")]
for i, L in enumerate(dc.CAND):
    for j, (k, a, h, lab) in enumerate(spec):
        b = ax.bar(i + (j - 1) * w, e[L][k], width=w * 0.95, color=dc.COLOURS[L], alpha=a,
                   hatch=h, edgecolor=dc.COLOURS[L])
        dc.bar_labels(ax, b, fmt="%+.2f", size=12)
ax.axhline(0, color="black", lw=1.0)
ax.set_xticks(range(len(dc.CAND))); ax.set_xticklabels(dc.CAND)
ax.set_xlabel("candidate, Condition CR")
ax.set_ylabel(r"drag increment against B  (counts)")
from matplotlib.patches import Patch
ax.legend([Patch(fc="0.35"), Patch(fc="0.75", ec="0.4", hatch="//"), Patch(fc="0.9", ec="0.4", hatch="..")],
          [s[3] for s in spec], loc="upper left", ncol=1)
ax.set_ylim(-2.7, 5.8)
dc.save(fig, "deck_magnitude", dc.RESULT_W, dc.RESULT_H)
