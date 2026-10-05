#!/usr/bin/env python3
"""Read the split decision log as one document, in date order (GEO-089).

THE SPLIT. docs/decisions.md became docs/decisions/{legacy,GEO,ARG,HPC}.md so that
no workstream ever has to commit a file containing another workstream's work. That was the
last remaining reason for a shared-index hazard, and it cost four near-misses in
one week. NO ENTRY WAS RENUMBERED AND NO CITATION WAS BROKEN.

THIS IS THE READ PATH, SO IT IS GATED LIKE ANYTHING ELSE HERE. A migration that
silently drops an entry is indistinguishable from one that does not, unless
something asserts otherwise. So:

    EVERY ASSERTION OVER A SET MUST PARTITION THE SET (standing rule). Compared
    plus deliberately-excluded-with-reason must equal all members; anything in
    neither is a FAILURE, not a skip.

applied to the migration itself. MANIFEST_AT_SPLIT.txt holds the 138 heading lines
frozen at the moment of the split, and --check asserts every one of them appears
EXACTLY ONCE in the concatenated output. Not "at least once": a duplicated heading
is as much a migration defect as a dropped one, and duplication is what actually
happened the last time entries moved between files.

The manifest is a KNOWN ANSWER captured before the move, not derived from the
result afterwards, which is the whole reason it can detect a bad move (D056).
"""
import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEC = REPO / "docs" / "decisions"
MANIFEST = DEC / "MANIFEST_AT_SPLIT.txt"
# legacy first only as a default; the real ordering is by each entry's own date.
FILES = ["legacy.md", "GEO.md", "ARG.md", "HPC.md"]
HEADING = re.compile(r"^## (\d{4}-\d{2}-\d{2}) ", re.M)


def entries(path):
    """Split one file into (date, heading_line, body) blocks. Text before the
    first heading is the file's own preamble and is not an entry."""
    text = path.read_text()
    marks = [(m.start(), m.group(1)) for m in HEADING.finditer(text)]
    out = []
    for i, (pos, date) in enumerate(marks):
        end = marks[i + 1][0] if i + 1 < len(marks) else len(text)
        block = text[pos:end]
        out.append((date, block.splitlines()[0], block, path.name))
    return out


def collect():
    all_entries = []
    missing = []
    for name in FILES:
        p = DEC / name
        if not p.exists():
            missing.append(name)
            continue
        all_entries.extend(entries(p))
    # Stable date order; ties keep file order, which is deterministic.
    all_entries.sort(key=lambda e: e[0])
    return all_entries, missing


def check(all_entries, missing, verbose=True):
    ok = True
    headings = [e[1] for e in all_entries]
    counts = {}
    for h in headings:
        counts[h] = counts.get(h, 0) + 1

    if not MANIFEST.exists():
        print("MANIFEST ABSENT: %s. The migration cannot be checked, so it is not "
              "asserted." % MANIFEST.relative_to(REPO))
        return False

    frozen = [l.rstrip("\n") for l in MANIFEST.read_text().splitlines() if l.strip()]
    lost = [h for h in frozen if counts.get(h, 0) == 0]
    dupes = [h for h in frozen if counts.get(h, 0) > 1]
    # New entries written since the split are EXPECTED and are counted, not
    # silently tolerated: the partition is frozen + added = all.
    added = [h for h in headings if h not in set(frozen)]
    added_dupes = sorted({h for h in added if counts[h] > 1})

    if verbose:
        print("DECISION-LOG MIGRATION CHECK (GEO-089)")
        print("  files read              %s" % ", ".join(
            "%s (%d)" % (n, sum(1 for e in all_entries if e[3] == n)) for n in FILES
            if (DEC / n).exists()))
        print("  frozen at split         %d" % len(frozen))
        print("  present now             %d" % (len(frozen) - len(lost)))
        print("  added since the split   %d" % len(added))
        print("  total headings          %d = %d frozen + %d added"
              % (len(headings), len(frozen) - len(lost), len(added)))
    if missing:
        ok = False
        print("  MISSING FILES           %s" % ", ".join(missing))
    if lost:
        ok = False
        print("  *** %d FROZEN HEADING(S) LOST IN THE MIGRATION ***" % len(lost))
        for h in lost[:10]:
            print("      %s" % h[:100])
    if dupes or added_dupes:
        ok = False
        print("  *** %d DUPLICATED HEADING(S) ***" % (len(dupes) + len(added_dupes)))
        for h in (dupes + added_dupes)[:10]:
            print("      %s" % h[:100])
    if len(headings) != (len(frozen) - len(lost)) + len(added):
        ok = False
        print("  *** PARTITION FAILURE: headings do not equal frozen + added ***")
    if verbose:
        print()
        print("MIGRATION CHECK: %s" % ("PASS, every frozen heading present exactly once"
                                       if ok else "FAIL"))
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="gate only, no document")
    ap.add_argument("--headings", action="store_true", help="list headings in order")
    args = ap.parse_args()

    all_entries, missing = collect()
    ok = check(all_entries, missing, verbose=True)
    if not ok:
        return 1
    if args.check:
        return 0
    print()
    if args.headings:
        for date, head, _body, src in all_entries:
            print("%-14s %s" % ("[%s]" % src.replace(".md", ""), head))
        return 0
    print("=" * 78)
    print("ARGUS decision log, all four files concatenated in date order")
    print("=" * 78)
    for _date, _head, body, src in all_entries:
        print(body.rstrip("\n"))
        print("<!-- source: docs/decisions/%s -->\n" % src)
    return 0


if __name__ == "__main__":
    sys.exit(main())
