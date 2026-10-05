#!/usr/bin/env python3
"""Surface delta maps: candidate minus baseline for c_p and signed c_f on the planform, upper
and lower, one PNG per quantity per condition. Same bins, same difference and same symmetric
99.5th-percentile scale as scripts/render_delta_maps_of12.py; no titles, letters only, one
colour bar. Writes fig/deck/deck_delta_maps.json with each figure's scale for the builder."""
import json
import os
import numpy as np
import deck_common as dc
from deck_maps_common import planform_mesh, load_map

dc.style()
COND = {"CR": ("SW", "Condition CR"), "EC": ("CMP", "early cruise")}
QTY = {"cp": ("cp", r"$\Delta c_p$", 1.0), "cfs": ("cf_s", r"$\Delta c_{f,s} \times 10^{3}$", 1e3)}
record = {}
for ctag, (prefix, clabel) in COND.items():
    B, bprov = load_map(prefix + "B_trim")
    X, Y = planform_mesh(B)
    cands = {g: load_map("%s%s_trim" % (prefix, g)) for g in dc.CAND}
    for g, (m, p) in cands.items():
        for k in ("condition", "mach", "U_inf", "n_eta", "n_xc"):
            if p.get(k) != bprov.get(k):
                raise SystemExit("HALT: %s %s %r differs from baseline %r" % (g, k, p.get(k), bprov.get(k)))
    for qtag, (field, qlabel, scale) in QTY.items():
        deltas = {(g, s): scale * (cands[g][0]["%s_%s" % (field, s)] - B["%s_%s" % (field, s)])
                  for g in dc.CAND for s in ("upper", "lower")}
        allv = np.concatenate([d[np.isfinite(d)].ravel() for d in deltas.values()])
        vmax = float(np.percentile(np.abs(allv), 99.5))
        fig, axes = dc.plt.subplots(2, 5, figsize=(dc.WIDE_W, 4.5))
        fig.set_layout_engine("none")
        fig.subplots_adjust(left=0.035, right=0.995, top=0.995, bottom=0.2, wspace=0.04, hspace=0.0)
        for row, s in enumerate(("upper", "lower")):
            for col, g in enumerate(dc.CAND):
                ax = axes[row, col]
                pm = ax.pcolormesh(Y, X, np.ma.masked_invalid(deltas[(g, s)]), cmap="coolwarm",
                                   vmin=-vmax, vmax=vmax, shading="nearest")
                ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                if row == 0:
                    ax.text(0.03, 0.95, g, transform=ax.transAxes, fontsize=20, fontweight="bold",
                            color=dc.COLOURS[g], va="top")
                if col == 0:
                    ax.set_ylabel(s, fontsize=17)
        cax = fig.add_axes([0.3, 0.125, 0.4, 0.035])
        cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
        cb.set_label("%s, candidate minus B, %s" % (qlabel, clabel), fontsize=16)
        cb.ax.tick_params(labelsize=14)
        name = "deck_delta_%s_%s" % (qtag, ctag)
        dc.save(fig, name, dc.WIDE_W, 4.5)
        record[name] = dict(scale=vmax, quantity=field, condition=clabel,
                            iterations={g: cands[g][1]["iteration"] for g in dc.CAND},
                            baseline_iteration=bprov["iteration"])
with open(os.path.join(dc.OUT, "deck_delta_maps.json"), "w") as fh:
    json.dump(record, fh, indent=1)
