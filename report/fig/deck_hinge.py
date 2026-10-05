#!/usr/bin/env python3
"""Slide: hinge moment. One panel of the region coefficient at the three operating points,
grouped by condition, one bar per geometry, labelled with its change against that
condition's own B. Source fig/hinge_moment_values.csv (as tab:hinge)."""
import numpy as np
import deck_common as dc

dc.style()
rows = dc.cload("fig/hinge_moment_values.csv")
conds = ["Condition CR", "Early cruise", "Late cruise"]
fig, ax = dc.plt.subplots()
w = 0.13
for ci, cname in enumerate(conds):
    sub = {r["geometry_letter"]: r for r in rows if r["condition"].startswith(cname)}
    base = float(sub["B"]["M_h_total_Nm"])   # the report's percentages are on the moment
    for gi, L in enumerate(dc.ORDER):
        v = -float(sub[L]["C_Mh_region_total"])
        x = ci + (gi - 2.5) * w
        b = ax.bar(x, v, width=w * 0.92, color=dc.COLOURS[L], label=L if ci == 0 else None)
        if L != "B":
            pct = 100.0 * (float(sub[L]["M_h_total_Nm"]) - base) / abs(base)
            ax.annotate("%+.1f%%" % abs(pct), (x, v), xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=12, rotation=90)
ax.set_xticks(range(3)); ax.set_xticklabels([c if c != "Condition CR" else "Condition CR" for c in conds])
ax.set_ylim(0, 0.062)
ax.set_ylabel(r"$-C_{M_h}$ about $x_h/c = 0.62$,  $0.60 \leq \eta \leq 1$")
ax.set_xlabel("operating point")
ax.legend(ncol=6, loc="upper left", handlelength=1.0, columnspacing=1.0)
dc.save(fig, "deck_hinge", dc.RESULT_W, dc.RESULT_H)
