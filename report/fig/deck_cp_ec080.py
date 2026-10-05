#!/usr/bin/env python3
"""Slide: why the cruise penalty appears. c_p at eta 0.80, early cruise, all six, one panel.

Upper surface solid, lower dashed, 200-bin medians as fig_cp_sections_earlycruise draws.
The sonic c_p* and the separated runs come from fig/cp_sections_measured.json (make_cp.py).
"""
import numpy as np
import deck_common as dc

dc.style()
meas = dc.jload("fig/cp_sections_measured.json")
fig, ax = dc.plt.subplots()
cp_star = None
for L in dc.ORDER:
    case = "CMP%s_trim" % L
    path, _gen = dc.section_file(case, 80)
    x, cp, _cfs, up = dc.read_section_csv(path)
    xu, cu = dc.binned_median(x[up], cp[up])
    xl, cl = dc.binned_median(x[~up], cp[~up])
    ax.plot(xu, cu, color=dc.COLOURS[L], lw=2.2, label=L)
    ax.plot(xl, cl, color=dc.COLOURS[L], lw=1.6, ls="--")
    sk = "EC|%s|eta080" % L
    cp_star = meas["shock"][sk]["cp_star"]
    for run in meas["separation"][sk]["runs"]:
        ax.axvspan(run["x_start"], run["x_end"], color=dc.COLOURS[L], alpha=0.16, lw=0)
ax.axhline(cp_star, color="0.35", lw=1.2, ls=":", label=r"$c_p^{*}$")
ax.set_xlim(0, 1)
ax.set_xlabel(r"$x/c$ at $\eta = 0.80$,  early cruise")
ax.plot([], [], color="0.3", lw=2.2, label="upper")
ax.plot([], [], color="0.3", lw=1.6, ls="--", label="lower")
ax.set_ylim(1.05, -1.42)
ax.set_ylabel(r"$c_p$  (inverted)")
ax.legend(ncol=5, loc="lower right", handlelength=1.6, columnspacing=0.9)
dc.save(fig, "deck_cp_ec080", dc.RESULT_W, dc.RESULT_H)
