#!/usr/bin/env python3
"""Every table row a generator emits must appear verbatim in the .tex that carries it.

WHY THIS EXISTS, and it is a defect of mine from 2026-09-24. After switching M's Condition
CR result from the published mesh to the wake-refined trim I declared "no stale M numbers
left anywhere in the tex" on the strength of a grep for FIVE patterns I happened to
remember: 191.246, +0.240, 1.108, 22.394, 0.9698. The trim ANGLE was not among them, so
tab:qv:lift kept M at the published 1.252018 with three further derived values formed from
it, and the report rebuilt clean. A separate review found it by reading the generator output.

THAT IS A VERDICT OVER A SET I CHOSE (GEO-087). A grep for remembered values can only ever
confirm what the author already suspected; it cannot find the row nobody thought of. The
check that works is the other direction: take what the GENERATORS produce, and require the
document to contain it. A row the generator emits and the tex lacks is then a failure by
construction rather than something a reader has to notice.

HOW IT MATCHES. Generators print their bodies after a line naming the table, e.g.
"%% table body for tab:qv:lift". Each following row that starts with a geometry letter and
contains "&" is a row. Whitespace is collapsed before comparison, because LaTeX alignment
padding is not content, and a padding-only difference reported as a mismatch is the
always-fires species of broken gate.

IT REPORTS BEFORE IT EXITS (GEO-089 item 3). Every missing row is printed, with the tex row
it most likely replaces, before any non-zero return. A coverage check that hides its own
diagnosis turns a useful finding into an unactionable alarm.
"""
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
PY = sys.executable

# generator -> the .tex files that may carry its rows
GENERATORS = {
    "scripts/make_rank_table.py": ["sections/t5_results.tex", "sections/t6_vlm_verdict.tex"],
    "scripts/make_cr_table.py": ["sections/t5_results.tex"],
}

ROW = re.compile(r"^\s*([A-Z])(\$\^\{?\*\}?\$)?\s*&")
MARKER = re.compile(r"table body for (tab:[A-Za-z0-9:_-]+)")


def norm(s):
    """Collapse whitespace; LaTeX column padding is not content."""
    return " ".join(s.replace("\\\\", "").split())


def rows_of(text):
    return [norm(l) for l in text.splitlines() if ROW.match(l) and "&" in l]


def emitted_by_table(text):
    """{table label: [rows]} from a generator's stdout, split on its own markers.

    KEYED BY TABLE, NOT LUMPED. A generator emits bodies for several tables, and some of
    those tables live in files this check does not carry, or were restructured away. The
    first version of this compared EVERY emitted row against BOTH target files and called
    sixteen of them stale when one was. That is the always-fires gate: it reports a defect
    in a healthy document and burns the trust the check exists to provide.
    """
    out, cur = {}, None
    for line in text.splitlines():
        m = MARKER.search(line)
        if m:
            cur = m.group(1); out.setdefault(cur, [])
            continue
        if cur and ROW.match(line) and "&" in line:
            out[cur].append(norm(line))
        elif cur and line.strip() and not line.startswith(" ") and "&" not in line:
            cur = None                      # body ended
    return {k: v for k, v in out.items() if v}


def main():
    texts = {p: (HERE / p).read_text() for p in sorted({q for v in GENERATORS.values() for q in v})}
    rows_in = {p: set(rows_of(t)) for p, t in texts.items()}

    checked, missing, not_in_doc, not_carried = 0, [], [], []
    for gen, targets in sorted(GENERATORS.items()):
        r = subprocess.run([PY, str(HERE / gen)], capture_output=True, text=True, cwd=str(HERE))
        if r.returncode != 0:
            print("HALT: %s exited %d" % (gen, r.returncode)); return 3
        for label, rows in sorted(emitted_by_table(r.stdout).items()):
            # WHICH FILE CARRIES THIS TABLE, decided by its \label, not by guessing.
            host = [p for p in targets if ("\\label{%s}" % label) in texts[p]]
            if not host:
                not_in_doc.append((label, len(rows))); continue
            for row in rows:
                checked += 1
                if any(row in rows_in[p] for p in host):
                    continue
                # STALE vs NOT CARRIED are different faults and must not share a bucket.
                # STALE: the document HAS this row, same leading letter and same column
                # count, with different content. That is the M defect: a wrong number
                # sitting in a live table.
                # NOT CARRIED: the generator emits a row the document simply does not
                # have, e.g. the B* baseline-corrected row that tab:qv:lift deliberately
                # ends before. Unused generator output is not a stale document.
                letter, ncol = row.split("&")[0].strip(), len(row.split("&"))
                twin = [c for p in host for c in rows_in[p]
                        if c.split("&")[0].strip() == letter and len(c.split("&")) == ncol]
                (missing if twin else not_carried).append((label, host, row, twin))

    print("GENERATED-TABLE RECONCILIATION")
    print("  rows checked against the table that declares their own label: %d" % checked)
    print("  tables emitted but NOT present in these .tex files: %d  %s"
          % (len(not_in_doc), ", ".join("%s(%d rows)" % t for t in not_in_doc) or "none"))
    print("  rows emitted but NOT CARRIED by the document (not a fault): %d  %s"
          % (len(not_carried), ", ".join("%s/%s" % (l, r.split("&")[0].strip())
                                         for l, _, r, _ in not_carried) or "none"))
    print("  STALE rows: %d\n" % len(missing))

    for label, host, row, twin in missing:
        print("  STALE in %s, table %s" % (", ".join(host), label))
        print("    generator : %s" % row)
        for c in twin:
            print("    document  : %s" % c)
        print()

    if checked == 0:
        print("*** NOTHING WAS COMPARED. That is not the same as agreeing.")
        return 3
    print("VERDICT: %s" % ("OK, every generated row is in the document"
                           if not missing else "*** %d STALE ***" % len(missing)))
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
