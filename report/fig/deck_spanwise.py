#!/usr/bin/env python3
"""Slide: spanwise load. Candidate minus baseline at Condition CR, RANS (solid) and VLM
(dashed), one panel. Same series construction as make_spanwise.py: RANS from
data/spanwise/<case>_spanwise_dense.csv (load = c_l c / C_ref), VLM from the .lod strip table
in data/vlm.json (load = Cl * Chord / Cref, eta = Yavg / (b/2)), each candidate minus its own
method's reference baseline interpolated onto the candidate's eta grid. The VLM reference is
the deck the report's difference panel uses (spanwise_caption_values.json
'vlm_reference_baseline'), read from that file rather than chosen here.
"""
import os
import numpy as np
import deck_common as dc

dc.style()
cap = dc.jload("fig/spanwise_caption_values.json")["condition_CR"]
V = dc.vlm()


def vlm_series(name):
    r = [x for x in V["runs"] if x["geometry_name"] == name
         and x["operating_point"]["state"] == "condition_CR"][0]
    rq = r["reference_quantities_AS_READ_FROM_DECK"]
    st = r["lod_strip_loads"]["strip_table"]
    col = {c: k for k, c in enumerate(st["_columns"])}
    A = np.asarray(st["rows"], dtype=float)
    eta = A[:, col["Yavg_ft"]] / (0.5 * rq["Bref_ft"])
    load = A[:, col["Cl"]] * A[:, col["Chord_ft"]] / rq["Cref_ft"]
    o = np.argsort(eta)
    return eta[o], load[o]


def rans_series(case):
    rows = dc.cload("data/spanwise/%s_spanwise_dense.csv" % case)
    eta = np.array([float(r["eta"]) for r in rows])
    load = np.array([float(r["load_cl_c_over_cref"]) for r in rows])
    o = np.argsort(eta)
    return eta[o], load[o]


rb = rans_series(cap["rans_cases"]["B"])
vb = vlm_series(cap["vlm_reference_baseline"])
fig, ax = dc.plt.subplots()
for L in dc.CAND:
    e, l = rans_series(cap["rans_cases"][L])
    ax.plot(e, l - np.interp(e, *rb), color=dc.COLOURS[L], lw=2.2, label=L)
    e, l = vlm_series(cap["vlm_geometry"][L])
    ax.plot(e, l - np.interp(e, *vb), color=dc.COLOURS[L], lw=1.8, ls="--")
ax.axhline(0, color="black", lw=0.8)
ax.plot([], [], color="0.3", lw=2.2, label="RANS")
ax.plot([], [], color="0.3", lw=1.8, ls="--", label="VLM")
ax.set_xlim(0, 1)
ax.set_xlabel(r"$\eta = y/(b/2)$,  Condition CR")
ax.set_ylabel(r"$\Delta(c_\ell\, c / C_\mathrm{ref})$ against B")
ax.legend(ncol=4, loc="lower right", handlelength=1.8)
dc.save(fig, "deck_spanwise", dc.RESULT_W, dc.RESULT_H)
