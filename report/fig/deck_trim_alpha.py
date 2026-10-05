#!/usr/bin/env python3
"""Slide: trim angle. RANS against VLM trim alpha at the common Condition CR target C_L."""
import numpy as np
import deck_common as dc

dc.style()
lift = dc.vlmerr()["lift"]
fig, ax = dc.plt.subplots()
xs, ys = [], []
for L in dc.ORDER:
    e = lift[L]
    xs.append(e["alpha_vlm_deg"]); ys.append(e["alpha_rans_deg"])
    ax.plot(e["alpha_vlm_deg"], e["alpha_rans_deg"], "o", ms=13, color=dc.COLOURS[L],
            mec="white", mew=1.2, label=L)
    ax.annotate("%s  %+.3f°" % (L, e["d_alpha_deg"]), (e["alpha_vlm_deg"], e["alpha_rans_deg"]),
                xytext=(12, -20 if L == "M" else (12 if L == "C" else -4)), textcoords="offset points", fontsize=14)
lo = min(xs + ys) - 0.12; hi = max(xs + ys) + 0.12
ax.plot([lo, hi], [lo, hi], color="0.4", lw=1.2, ls="--", label="1:1")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_aspect("equal")
ax.set_xlabel(r"VLM trim $\alpha$  (deg)")
ax.set_ylabel(r"RANS trim $\alpha$  (deg)")
ax.legend(loc="upper left", ncol=2)
dc.save(fig, "deck_trim_alpha", dc.RESULT_W, dc.RESULT_H)
