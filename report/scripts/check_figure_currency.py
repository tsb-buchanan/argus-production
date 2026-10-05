#!/usr/bin/env python3
"""Report how current every figure the document actually includes is.

WHY IT EXISTS. latexmk -halt-on-error proves every figure RESOLVES. It says nothing about
whether the file it resolved to is the one the current data would produce. Four days of
stale delta maps compiled clean every single time, because a stale figure is a perfectly
valid figure of the wrong data. That is the input-gate lesson (GEO-102) applied to the
report: the compile is an output gate and cannot see a stale input.

IT RESOLVES PATHS THE WAY LaTeX DOES, and that is the whole reason this file exists rather
than a one-line find. My first version assumed every \\includegraphics argument was
repo-relative and reported FIVE MISSING figures that were present and current, because
argus_rans_validation_report.tex declares \\graphicspath{{fig/}} and those five omit the prefix. A gate that
reports a defect in a healthy artefact is the always-fires species of broken gate, and it
burns exactly the trust the gate exists to provide. So the search path is PARSED from the
document rather than assumed, and a figure counts as missing only after every candidate
path has been tried.

BUCKETS ARE A PARTITION AND IT IS ASSERTED. current + stale + missing must equal the
number of distinct figures included, or the report is over a set that was silently
narrowed (GEO-087).

STALE IS NOT AUTOMATICALLY A DEFECT, so this PRINTS rather than fails by default. Several
figures are provenance diagrams that no solver result feeds; they are correctly older than
the data. Use --fail-stale in a pipeline that wants the stricter reading, and --since to
set the cutoff.
"""
import argparse
import collections
import datetime
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent

# Figures that legitimately predate the run data, with the REASON pinned. A bare list of
# names to skip would let a genuinely stale figure hide behind an exemption, so each entry
# has to say why it does not depend on solver output (the same rule as pinning an
# exclusion's expected values rather than just naming the field).
STATIC_OK = {
    "concept_displacement_3d.png":
        "geometry concept diagram, drawn from the OML; no solver field feeds it",
    "fig_lineage_staircase.pdf":
        "provenance staircase of geometry lineage; carries no run result",
}


def search_paths(main):
    """The directories LaTeX will look in, parsed from the document itself."""
    paths = [HERE]
    m = re.search(r"\\graphicspath\{(.+?)\}\s*$", main.read_text(), re.M)
    if m:
        for d in re.findall(r"\{([^}]*)\}", m.group(1)):
            paths.append(HERE / d)
    return paths


def resolve(arg, paths):
    for base in paths:
        for cand in (base / arg, base / (arg + ".pdf"), base / (arg + ".png")):
            if cand.exists():
                return cand
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--since", default="2026-09-22",
                    help="figures older than this date are reported as stale")
    ap.add_argument("--fail-stale", action="store_true",
                    help="exit non-zero if any non-exempt figure is stale")
    a = ap.parse_args()
    cut = datetime.datetime.strptime(a.since, "%Y-%m-%d")

    mainf = HERE / "argus_rans_validation_report.tex"
    paths = search_paths(mainf)
    tex = sorted(HERE.glob("sections/*.tex")) + [mainf]

    inc = collections.OrderedDict()
    for p in tex:
        for m in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", p.read_text()):
            inc.setdefault(m.group(1), set()).add(p.name)

    cur, stale, static, missing = [], [], [], []
    for arg, users in inc.items():
        q = resolve(arg, paths)
        if q is None:
            missing.append((arg, users, "")); continue
        mt = datetime.datetime.fromtimestamp(q.stat().st_mtime)
        row = (q.name, users, mt.strftime("%Y-%m-%d %H:%M"))
        if mt >= cut:
            cur.append(row)
        elif q.name in STATIC_OK:
            static.append(row)
        else:
            stale.append(row)

    print("FIGURE CURRENCY AUDIT  (cutoff %s)" % a.since)
    print("  search path: %s" % ", ".join(str(p.relative_to(HERE)) or "." for p in paths))
    print("  %d distinct figures included across %d tex files\n" % (len(inc), len(tex)))

    for name, rows in (("MISSING", missing), ("STALE", stale),
                       ("STATIC, exempt with reason", static), ("CURRENT", cur)):
        print("== %s : %d ==" % (name, len(rows)))
        for f, users, mt in sorted(rows):
            note = "   [%s]" % STATIC_OK[f] if f in STATIC_OK else ""
            print("   %-44s %-17s <- %s%s" % (f, mt or "NOT FOUND",
                                              ",".join(sorted(users)), note))
        print()

    # ---- IS THE BUILT PDF OLDER THAN ITS OWN INPUTS? --------------------------------
    # THE GATE ABOVE CANNOT SEE THIS. It compares each figure against a DATE CUTOFF, so it
    # reported "0 stale" at a moment when argus_rans_validation_report.pdf was older than two of the files it is
    # built from. That is the output-gate blind spot one level up: every figure was current
    # and the artefact assembled from them was not.
    #
    # AND AN MTIME DIFFERENCE IS NOT PROOF OF A CONTENT CHANGE. matplotlib stamps a
    # CreationDate into every PDF it writes, so re-running a generator always moves the
    # bytes even when the drawing is identical. Measured instance: fig_spaneff.pdf
    # regenerated to the same 24990 bytes with identical rendered text and only the
    # CreationDate differing. So this REPORTS the ordering and says what it does not prove,
    # rather than demanding a rebuild the content may not need.
    pdf = HERE / "argus_rans_validation_report.pdf"
    print("== BUILT ARTEFACT ==")
    if not pdf.exists():
        print("   argus_rans_validation_report.pdf ABSENT\n")
    else:
        pmt = datetime.datetime.fromtimestamp(pdf.stat().st_mtime)
        inputs = []
        for arg in inc:
            q = resolve(arg, paths)
            if q:
                inputs.append((q.name, datetime.datetime.fromtimestamp(q.stat().st_mtime)))
        for f in sorted(HERE.glob("sections/*.tex")) + [mainf]:
            inputs.append((f.name, datetime.datetime.fromtimestamp(f.stat().st_mtime)))
        newer = sorted([(n, m) for n, m in inputs if m > pmt], key=lambda r: r[1],
                       reverse=True)
        print("   argus_rans_validation_report.pdf  %s" % pmt.strftime("%Y-%m-%d %H:%M:%S"))
        if not newer:
            print("   NEWER THAN ALL %d inputs.\n" % len(inputs))
        else:
            print("   OLDER THAN %d of %d inputs, most recent first:" % (len(newer), len(inputs)))
            for n, m in newer[:6]:
                print("      %-44s %s" % (n, m.strftime("%Y-%m-%d %H:%M:%S")))
            print("   THIS IS AN ORDERING, NOT A DIFF. A generator re-run moves a figure's")
            print("   bytes via its embedded CreationDate even when the drawing is")
            print("   unchanged. Compare rendered content before concluding a rebuild is")
            print("   needed: pdftotext the old and new figure and diff those.\n")

    total = len(cur) + len(stale) + len(static) + len(missing)
    ok = total == len(inc)
    print("PARTITION: %d included = %d current + %d stale + %d static + %d missing  %s"
          % (len(inc), len(cur), len(stale), len(static), len(missing),
             "OK" if ok else "*** DOES NOT RECONCILE ***"))
    if not ok:
        return 3
    if missing:
        return 2
    if stale and a.fail_stale:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
