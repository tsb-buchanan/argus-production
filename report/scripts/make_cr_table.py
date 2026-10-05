#!/usr/bin/env python3
"""Emit the Condition CR performance table rows as LaTeX, read from the force data.

WHY GENERATE RATHER THAN TYPE. The rows in sections/06_performance.tex were hand-written,
and a hand-typed number has already gone wrong once in this work: an SWH x2/x3 pair was
mistyped into a report table and only caught because the value was later recomputed from
the file. Record-keeping rule 3 says plots regenerate from a script; a table carrying the
same numbers deserves the same treatment.

RANK IS COMPUTED, NOT CARRIED. The published table ranked C, M, F, W, H as 1 to 5. On the
current numbers that order changes, and a stale rank column beside fresh deltas is exactly
the kind of half-updated table a reader trusts because most of it is right.

C AND M CARRY THEIR OWN BASELINE. Each geometry's delta is read from the data file, where
it is stored with the baseline it was formed against. The alternative, recomputing every
delta against a single baseline, would give C +2.83 instead of -0.149, a sign reversal
caused entirely by the baseline having moved (D071).
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
DATA = HERE / "data/rans_forces.json"
ORDER_ID = ["B", "C", "M", "F", "W", "H"]
CASE = {"B": "SWB", "C": "SWC", "M": "SWM", "F": "SWF", "W": "SWW", "H": "SWH"}
TARGET_CL = 0.428277635108


def main():
    d = json.loads(DATA.read_text())
    cr = d["trims"]["condition_CR"]["geometries"]

    rows = []
    for g in ORDER_ID:
        e = cr[g]
        rows.append(dict(
            g=g, case=CASE[g] + r"\_trim", alpha=e["alpha_deg"], cd=e["Cd_counts"],
            dcd=e["delta_Cd_counts_vs_baseline"],
            dcl_ct=(e["Cl"] - TARGET_CL) * 1e4, ld=e["L_over_D"],
        ))

    ranked = sorted([r for r in rows if r["g"] != "B"], key=lambda r: r["dcd"])
    rank = {r["g"]: i + 1 for i, r in enumerate(ranked)}

    out = []
    for r in rows:
        rk = "--" if r["g"] == "B" else str(rank[r["g"]])
        out.append("%s & %s & %.4f & %.3f & $%+.3f$ & $%+.3f$ & %.3f & %s \\\\"
                   % (r["g"], r["case"], r["alpha"], r["cd"], r["dcd"],
                      r["dcl_ct"], r["ld"], rk))
    print("\n".join(out))
    print()
    print("%% rank order by dCd, cheapest first: %s"
          % ", ".join(r["g"] for r in ranked), file=sys.stderr)
    trimmed = [r for r in rows if abs(r["dcl_ct"]) <= 1.0]
    print("%% %d of %d rows within 1.0 count of target CL"
          % (len(trimmed), len(rows)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
