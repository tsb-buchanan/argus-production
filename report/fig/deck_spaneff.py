#!/usr/bin/env python3
"""Slide: induced drag and span efficiency. RANS e at the last Trefftz plane, e = 1 marked.
Every geometry is drawn, with its value, in one style.

WHY ONE STYLE. An entry served from the published evaluation used to be drawn hatched and
half-opaque. Project decision, 2026-09-23: the slide carries the result, not a commentary on
which evaluation produced it. The generation difference is RECORDED in
data/span_efficiency.json, which is where the frame belongs (D068 requires the frame to
travel with the number in the record, not that it be drawn on every bar).

AND IT IS BOUNDED, NOT ASSUMED AWAY. Measured on the two geometries evaluated BOTH ways,
the generation shifts e by 0.0024 (M: 0.9698 published, 0.9722 wake-refined) and 0.0019
(B: 0.9654, 0.9673). The axis below spans 0.12 in e, so the shift is about 2 per cent of
the axis and changes no ordering. The geometries where generation DID matter, W H and F at
-0.30 in e, are all served from the current evaluation."""
import deck_common as dc

dc.style()
s = dc.spaneff()
fig, ax = dc.plt.subplots()
xs = list(range(len(dc.ORDER)))
pub = s.get("published_evaluation", {}).get("geometries", {})
missing = [L for L in dc.ORDER
           if L not in s["subset"]["emitted"] and L not in pub]
if missing:
    raise SystemExit("HALT: no span efficiency for %s in either block; the slide would "
                     "carry a bar-less tick and claim the whole field" % ", ".join(missing))
for i, L in enumerate(dc.ORDER):
    e = (s["geometries"][L] if L in s["subset"]["emitted"] else pub[L])["e_last"]
    b = ax.bar(i, e, color=dc.COLOURS[L], width=0.62)
    dc.bar_labels(ax, b, fmt="%.4f")
ax.axhline(1.0, color="0.3", lw=1.3, ls="--", label=r"$e = 1$")
ax.set_xticks(xs); ax.set_xticklabels(dc.ORDER)
ax.set_ylim(0.90, 1.02)
ax.set_xlabel("geometry, Condition CR, last Trefftz plane")
ax.set_ylabel(r"RANS span efficiency $e$")
ax.legend(loc="upper left")
dc.save(fig, "deck_spaneff", dc.RESULT_W, dc.RESULT_H)
