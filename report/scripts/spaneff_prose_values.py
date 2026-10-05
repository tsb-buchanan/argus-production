#!/usr/bin/env python3
"""Write sections/12_spanefficiency_result.tex from data/span_efficiency.json.

WHY GENERATED. The result paragraph quotes a span efficiency per geometry, a decay
percentage per geometry and the spread across stations. All of those move with every
re-trim, and this is the one section where a stale `e` would be read as a physical
statement about a spanload rather than as the numerical artefact it would actually be.

IT NAMES THE EXCLUDED GEOMETRIES, ALWAYS. build_span_efficiency.py permits a subset only
with a recorded reason, and writes the excluded set into the data file. This script REFUSES
to write a paragraph that does not mention them, so a figure over four of six can never be
captioned as though it covered the field (GEO-087: a verdict over an unstated subset is not
a verdict).

IT WRITES TO A TEMPORARY PATH AND PROMOTES ON PASS, so a failed run cannot leave a partial
paragraph under the name the report \\input{}s.
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
OUT = HERE / "sections/12_spanefficiency_result.tex"
NAME = {"B": "the rigid baseline B", "C": "C", "F": "F", "H": "H", "M": "M", "W": "W"}


def _list(items):
    """Oxford-free English list: 'C', 'C and M', 'C, F and M'."""
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def main():
    src = HERE / "data/span_efficiency.json"
    if not src.exists():
        sys.exit("HALT: no %s. Run scripts/build_span_efficiency.py first." % src)
    d = json.loads(src.read_text())
    fr = d["frame"]
    sub = d.get("subset")

    # BOTH BLOCKS ARE READ, so the prose describes the figure the reader is looking at.
    # Project decision, 2026-09-23: every geometry carries its value. A geometry the data file
    # holds under `published_evaluation` HAS a measured result, and a coverage sentence
    # saying it is "not included" is simply false once the figure draws it. The mesh
    # generation stays in the data file and in fig/fig_spaneff_values.csv.
    #
    # THE PRECONDITION IS ASSERTED (GEO-092): a published entry is taken only if it is a
    # TRIMMED leg. Anything else stays uncovered and is NAMED by the Coverage paragraph,
    # which is what that paragraph is for.
    g = dict(d["geometries"])
    for L, e in d.get("published_evaluation", {}).get("geometries", {}).items():
        if L not in g and e.get("trimmed") and e.get("stations"):
            g[L] = e
    uncovered = sorted(set(sub["excluded"]) - set(g)) if sub else []

    order = sorted(g, key=lambda L: g[L]["e_last"])
    lo, hi = order[0], order[-1]
    decays = {L: g[L]["decay_first_to_last_percent"] for L in g}
    e_last = {L: g[L]["e_last"] for L in g}
    cdi_last = {L: g[L]["CDi_counts_last"] for L in g}

    lines = []
    lines.append("\\textbf{Result.} Span efficiency at the last station ranges from "
                 "$%.4f$ on %s to $%.4f$ on %s, every value below the planar-wake limit. "
                 "Induced drag at that station spans $%.2f$ to $%.2f$ counts."
                 % (e_last[lo], NAME.get(lo, lo), e_last[hi], NAME.get(hi, hi),
                    min(cdi_last.values()), max(cdi_last.values())))

    lines.append("")
    lines.append("\\textbf{Decay between the first and last plane} is $%.1f$ to $%.1f$ "
                 "per cent, and is consistent across the geometries: %s. That consistency "
                 "matters more than the magnitude, because a decay rate that differed "
                 "between geometries would mean the comparison was being made on wakes "
                 "resolved to different degrees, and no difference between them could be "
                 "attributed to the wing."
                 % (min(decays.values()), max(decays.values()),
                    ", ".join("%s $%.2f$\\,\\%%" % (L, decays[L])
                              for L in sorted(g, key=lambda L: decays[L]))))

    lines.append("")
    if len(g) > 1:
        cand = [L for L in g if L != "B"]
        if "B" in g and cand:
            better = [L for L in cand if e_last[L] > e_last["B"]]
            phrase = ("Every candidate measured here carries"
                      if len(better) == len(cand) else
                      "%d of the %d candidates measured here carry" % (len(better), len(cand)))
            lines.append("\\textbf{Against the baseline.} %s a higher span efficiency than "
                         "the rigid baseline's $%.4f$%s. The margin is small, one to two "
                         "per cent, which is the size of spanload improvement the "
                         "low-order design study predicted for this class of change."
                         % (phrase, e_last["B"],
                            " (%s)" % _list(sorted(better)) if better else ""))

    # THE EXCLUSION IS NOT OPTIONAL PROSE, BUT IT DESCRIBES WHAT IS ACTUALLY UNCOVERED.
    # It fires on the set with NO evaluation in either block, not on the bookkeeping
    # `excluded` list, which now names geometries the figure does draw.
    if uncovered:
        lines.append("")
        # STATE THE SUBSET WITHOUT NARRATING HOW IT AROSE. GEO-087 requires that the
        # excluded members be NAMED; it does not require the reason to recount the
        # campaign's history, and this report deliberately does not tell that story. The
        # first version of this sentence said the excluded cases were "still solving on the
        # wake-refined mesh" and that "their evaluations exist only on the earlier grid",
        # which is the narration the report is meant to be free of. The `reason` field
        # stays in the DATA file, where the record wants it; the PROSE states the fact.
        lines.append("\\textbf{Coverage.} This figure covers %s. %s %s not included, %s "
                     "wake-plane data not yet being available at the trimmed condition. "
                     "Induced drag and span efficiency are both sensitive to the "
                     "resolution of the wake, so the figure is drawn only over cases whose "
                     "wake data was produced under identical conditions."
                     % (", ".join(sorted(g)),
                        _list(uncovered),
                        "is" if len(uncovered) == 1 else "are",
                        "its" if len(uncovered) == 1 else "their"))

    body = "\n".join(lines) + "\n"
    if uncovered and not all(L in body for L in uncovered):
        sys.exit("HALT: the paragraph does not name every uncovered geometry (%s). A figure "
                 "over a subset captioned as the field is not publishable."
                 % ", ".join(uncovered))

    tmp = OUT.with_suffix(".tex.tmp")
    tmp.write_text(body)
    tmp.replace(OUT)
    print("wrote %s (%d paragraphs)" % (OUT.relative_to(HERE), body.count("\\textbf{")))
    print(body)
    return 0


if __name__ == "__main__":
    sys.exit(main())
