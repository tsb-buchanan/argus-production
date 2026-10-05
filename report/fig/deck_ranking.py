#!/usr/bin/env python3
"""Slide: ordering. Rank lines, VLM on the left and RANS on the right, one PNG per VLM
metric, from fig_vlmerr_values.json 'ranking' (the data behind fig_vlmerr_ranking.pdf)."""
import deck_common as dc

dc.style()
rk = dc.vlmerr()["ranking"]
for key, name, lab in [("vlm_rank_by_dCDi_near", "near", "VLM near field\n(surface)"),
                       ("vlm_rank_by_dCDiw_far", "far", "VLM far field\n(wake)")]:
    fig, ax = dc.plt.subplots()
    for L in dc.CAND:
        v, r = rk[key][L], rk["rans_rank_by_dCd_total"][L]
        ax.plot([0, 1], [v, r], "-o", color=dc.COLOURS[L], lw=2.6, ms=11, mec="white", mew=1.2)
        ax.annotate(L, (0, v), xytext=(-14, 0), textcoords="offset points", ha="right",
                    va="center", fontsize=16, color=dc.COLOURS[L], fontweight="bold")
        ax.annotate(L, (1, r), xytext=(14, 0), textcoords="offset points", ha="left",
                    va="center", fontsize=16, color=dc.COLOURS[L], fontweight="bold")
    ax.set_xlim(-0.45, 1.45); ax.set_ylim(5.6, 0.4)
    ax.set_xticks([0, 1]); ax.set_xticklabels([lab, "RANS"])
    ax.set_yticks([1, 2, 3, 4, 5])
    ax.set_ylabel("rank by drag increment (1 = least)")
    ax.grid(False)
    dc.save(fig, "deck_ranking_%s" % name, dc.NARROW_W, dc.RESULT_H)
