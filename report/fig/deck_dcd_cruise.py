#!/usr/bin/env python3
"""Slide: trimmed drag at early and late cruise. Two panels, one shared scale."""
import deck_common as dc

dc.style()
fig, axes = dc.plt.subplots(1, 2, sharey=True)
allv = []
for ax, cond, lab in zip(axes, ["early_cruise", "late_cruise"], ["early cruise", "late cruise"]):
    t = dc.trims(cond)
    vals = [t[L]["delta_Cd_counts_vs_baseline"] for L in dc.CAND]
    allv += vals
    bars = ax.bar(dc.CAND, vals, color=[dc.COLOURS[L] for L in dc.CAND], width=0.62)
    ax.axhline(0, color="black", lw=1.0)
    dc.bar_labels(ax, bars, fmt="%+.1f")
    ax.set_xlabel("candidate, " + lab)
axes[0].set_ylabel(r"$\Delta C_D$ against B  (counts)")
# The floor follows the data: a fixed 0 hid W's saving (-0.77 ct) once the welded trims landed,
# which is the headline of this slide. A negative bar gets room for its label below it.
lo, hi = min(allv), max(allv)
pad = 0.08 * (hi - min(lo, 0.0))
axes[0].set_ylim(min(0.0, lo) - (2.5 * pad if lo < 0 else 0.0), hi * 1.15)
dc.save(fig, "deck_dcd_cruise", dc.RESULT_W, dc.RESULT_H)
