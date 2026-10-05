#!/usr/bin/env python3
"""Slide: trimmed drag, Condition CR. Delta C_D per candidate against B, one panel."""
import deck_common as dc

dc.style()
t = dc.trims("condition_CR")
fig, ax = dc.plt.subplots()
vals = [t[L]["delta_Cd_counts_vs_baseline"] for L in dc.CAND]
bars = ax.bar(dc.CAND, vals, color=[dc.COLOURS[L] for L in dc.CAND], width=0.62)
ax.axhline(0, color="black", lw=1.0)
dc.bar_labels(ax, bars)
ax.set_ylabel(r"$\Delta C_D$ against B  (counts)")
ax.set_xlabel("candidate, Condition CR")
lo, hi = min(vals), max(vals)
ax.set_ylim(lo - 0.25 * (hi - lo) - 0.1, hi + 0.25 * (hi - lo) + 0.1)
dc.save(fig, "deck_dcd_cr", dc.RESULT_W, dc.RESULT_H)
