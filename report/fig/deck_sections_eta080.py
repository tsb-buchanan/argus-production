#!/usr/bin/env python3
"""Slide: surface pressure and skin friction at eta 0.80, one PNG per condition, c_p over
c_f,s, all six geometries (upper solid, lower dashed), separated runs shaded. Section files
resolved per case through fig/section_source.py, as the report figures are."""
import numpy as np
import deck_common as dc

dc.style()
meas = dc.jload("fig/cp_sections_measured.json")
for tag, prefix in [("CR", "SW"), ("EC", "CMP")]:
    fig, (a1, a2) = dc.plt.subplots(2, 1, sharex=True, height_ratios=[1.35, 1.0])
    for L in dc.ORDER:
        path, _ = dc.section_file("%s%s_trim" % (prefix, L), 80)
        x, cp, cfs, up = dc.read_section_csv(path)
        xu, cu = dc.binned_median(x[up], cp[up]); xl, cl = dc.binned_median(x[~up], cp[~up])
        a1.plot(xu, cu, color=dc.COLOURS[L], lw=2.0, label=L)
        a1.plot(xl, cl, color=dc.COLOURS[L], lw=1.4, ls="--")
        xu, fu = dc.binned_median(x[up], cfs[up]); xl, fl = dc.binned_median(x[~up], cfs[~up])
        a2.plot(xu, fu, color=dc.COLOURS[L], lw=2.0)
        a2.plot(xl, fl, color=dc.COLOURS[L], lw=1.4, ls="--")
        for run in meas["separation"]["%s|%s|eta080" % (tag, L)]["runs"]:
            for a in (a1, a2):
                a.axvspan(run["x_start"], run["x_end"], color=dc.COLOURS[L], alpha=0.16, lw=0)
    a1.set_ylim(1.05, -1.42)
    a2.axhline(0, color="black", lw=0.9)
    a1.set_xlim(0, 1)
    a1.set_ylabel(r"$c_p$  (inverted)")
    a2.set_ylabel(r"$c_{f,s}$")
    a2.set_xlabel(r"$x/c$ at $\eta = 0.80$,  %s" % ("Condition CR" if tag == "CR" else "early cruise"))
    a1.legend(ncol=6, loc="lower right", handlelength=1.4, columnspacing=0.9)
    dc.save(fig, "deck_sections_%s_eta080" % tag, dc.HALF_W, dc.HALF_H)
