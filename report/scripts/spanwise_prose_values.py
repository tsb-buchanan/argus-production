#!/usr/bin/env python3
"""Print the numbers the spanwise prose quotes, read from fig/spanwise_caption_values.json.

WHY THIS EXISTS. Two passages quote values straight off the spanwise figure:

  sections/s6_low_order.tex    the r.m.s. of RANS minus VLM on the DIFFERENCE curves, the
                              range of RANS peak amplitudes, the r.m.s. as a percentage of
                              the peak, and the largest disagreement in the station of peak
                              change (section 6.3 of the restructured report).
  sections/s4_condition_cr.tex the eta at which each geometry's RANS difference peaks, in the
                              caption of fig:geom:spanwise (section 4.3).

Every one of those moves when the spanwise extraction is regenerated, which happens each
time a geometry is re-trimmed. Recomputing them by hand is how a half-updated paragraph
gets written: most of it right, one figure stale, and nothing to flag it. Run this after
every `python3 fig/make_spanwise.py` and paste what it prints.

THE PERCENTAGE IS RECOMPUTED FROM THE TWO QUANTITIES, never carried. r.m.s. and peak are
both in the file; a percentage stored beside them is a third number that can go stale
independently of the two it comes from.
"""
import json
import sys
from pathlib import Path

# THE DATA ROOT IS FOUND, NOT ASSUMED. This was parent.parent, which is the report
# root only when data/, fig/ and scripts/ are siblings. The delivery repo lifts data/
# to its own root, so the assumption held in one layout and failed silently in the
# other. data_root() walks up for the directory whose data/ carries rans_forces.json
# and RAISES if there is none, rather than returning a plausible wrong tree.
# In the working repo and in an assembled build directory it returns exactly what
# parent.parent returned, so this changes no behaviour there.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_paths import data_root  # noqa: E402
HERE = data_root()


def main():
    src = HERE / "fig/spanwise_caption_values.json"
    if not src.exists():
        sys.exit("HALT: no %s. Run fig/make_spanwise.py first." % src)
    d = json.loads(src.read_text())["condition_CR"]

    rm = d["rans_minus_vlm_on_the_difference"]
    pk = d["peak_difference"]
    letters = sorted(rm)

    print("Condition CR, RANS against VLM on the DIFFERENCE curves")
    print("%-3s %10s %10s %9s %10s %10s %9s"
          % ("id", "rms", "max_abs", "n", "RANS peak", "rms/peak", "d(eta_pk)"))
    rows = []
    for L in letters:
        rms = rm[L]["rms"]
        peak, eta_r = pk[L]["rans"]
        eta_v = pk[L]["vlm"][1]
        frac = 100.0 * rms / peak
        deta = abs(eta_r - eta_v)
        rows.append((L, rms, rm[L]["max_abs"], rm[L]["n_points"], peak, frac, deta))
        print("%-3s %10.6f %10.6f %9d %10.6f %9.2f%% %9.4f"
              % (L, rms, rm[L]["max_abs"], rm[L]["n_points"], peak, frac, deta))

    rmin, rmax = min(r[1] for r in rows), max(r[1] for r in rows)
    pmin, pmax = min(r[4] for r in rows), max(r[4] for r in rows)
    fmin, fmax = min(r[5] for r in rows), max(r[5] for r in rows)
    dmax = max(r[6] for r in rows)
    nset = sorted({r[3] for r in rows})

    print("\n--- for sections/s6_low_order.tex (6.3 Load redistribution) ---")
    print("  r.m.s. range          %.4f to %.4f" % (rmin, rmax))
    print("  RANS peak range       %.3f to %.3f" % (pmin, pmax))
    print("  r.m.s. as %% of peak   %.1f to %.1f per cent" % (fmin, fmax))
    print("  peak station agrees to within %.3f" % dmax)
    print("  stations per curve    %s" % (nset[0] if len(nset) == 1 else nset))

    print("\n--- for sections/s4_condition_cr.tex (fig:geom:spanwise caption) ---")
    by_eta = {}
    for L in letters:
        by_eta.setdefault(round(pk[L]["rans"][1], 4), []).append(L)
    for eta in sorted(by_eta, reverse=True):
        print("  eta %.4f: %s" % (eta, ", ".join(sorted(by_eta[eta]))))

    # The geometries whose loads are not yet on the current mesh generation, named rather
    # than left for the reader to work out, because the sentence above spans all of them.
    prov = HERE / "data/spanwise"
    times = {}
    for L in letters:
        case = d["rans_cases"][L]
        j = prov / ("%s_spanwise_dense.json" % case)
        if j.exists():
            times[L] = json.loads(j.read_text()).get("time")
    print("\n--- extraction time per geometry (a mixed set is a frame mix, D068) ---")
    for L in sorted(times):
        print("  %s %-9s t=%s" % (L, d["rans_cases"][L], times[L]))
    distinct = sorted({str(t) for t in times.values()})
    if len(distinct) > 1:
        print("  NOTE: %d distinct extraction times %s. Geometries at the older time are "
              "still on the published mesh." % (len(distinct), distinct))
    return 0


if __name__ == "__main__":
    sys.exit(main())
