#!/usr/bin/env python3
"""Refuse to let the report narrate the September mesh work.

WHY THIS EXISTS, AND WHY IT IS NOT THE DECK'S COPY. The deck has had this check since it
was found to be narrating the mesh work on four slides and five speaker notes while its own
docstring claimed it did not. The REPORT never got the same check, and it had the same
defect: the span-efficiency coverage paragraph, which an earlier edit generated, read "both are
still solving on the wake-refined mesh. Their evaluations exist only on the earlier grid".
A gate that exists for one artefact and not its sibling is not a gate on the requirement,
it is a gate on the artefact someone happened to think about.

WHAT IT PERMITS. Naming an excluded case is REQUIRED (GEO-087). Explaining that exclusion by
recounting which grid produced what is NOT, and is what this refuses. The distinction is
between "C and M are not included, their wake-plane data not yet being available" and "C and
M are still on the earlier grid". Both state the subset; only one tells the story.

FALSE POSITIVES ARE EXPECTED AND ARE HANDLED BY CONTEXT, NOT BY WEAKENING THE LIST. The word
"defect" appears legitimately about a defective STL surface and about a defect in a
delivered record, and "superseded" about a bin count. Those are ALLOWED_CONTEXT below, each
pinned to the file and phrase it occurs in, so a NEW occurrence of the same word still
fires. A blanket removal of the term would be the fifth species of broken gate: an exemption
that lets through exactly the input it was not written for.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent

# Terms that would tell the reader a mesh was changed mid-campaign.
FORBIDDEN = [
    "mesh defect", "mesh artefact", "mesh artifact", "starved mesh", "starved meshes",
    "re-mesh", "remesh", "re-meshed", "remeshed", "mesh generation", "mesh generations",
    "earlier grid", "previous grid", "old mesh", "previous mesh", "corrected mesh",
    "wake-refined", "wake refined", "the re-run", "before the re-trim",
]

# Occurrences that are about something else entirely, pinned to file AND surrounding
# phrase so that a NEW use of the same word is still caught.
# File names are those of the reordered report (22 September 2026, afternoon pass): Appendix A
# took the section-cut provenance, t5 (results by quantity) took hinge, root bending and span
# efficiency. The
# twelve pre-restructure files were in sec_new/superseded/, which the glob below never
# reached, so their pins were removed rather than left pointing at nothing. That directory
# has since been renamed aside pending deletion (working repo) and removed (delivery), so
# the path named here no longer exists ANYWHERE; it is kept only to say why those pins are
# absent. Raised in review: a rename pass rewrote sec_new to sections inside this
# comment and turned a historical statement into a claim about a live path.
ALLOWED_CONTEXT = [
    ("sections/tA_sections.tex", "A defective surface would be ragged"),
    ("sections/tA_sections.tex", "The third was a defect and"),
    ("sections/tA_sections.tex", "the superseded set carrying the geometric split"),
    ("sections/t5_results.tex", "defect in the delivered record"),
    ("sections/t5_results.tex", "is superseded"),
    # The span-efficiency caption names the SCRIPT's behaviour, not the campaign's history.
    ("sections/t5_results.tex", "refuses to mix mesh generations"),
]


def allowed(rel, line):
    return any(f == rel and phrase in line for f, phrase in ALLOWED_CONTEXT)


def main():
    files = sorted(HERE.glob("sections/*.tex")) + [HERE / "argus_rans_validation_report.tex"]
    hits, scanned = [], 0
    for p in files:
        if not p.exists():
            continue
        scanned += 1
        rel = str(p.relative_to(HERE))
        # MATCH ACROSS LINE BREAKS. This scanned line by line, so a forbidden term split
        # by LaTeX's wrapping was invisible: "two mesh\ngenerations" passed cleanly while
        # "two mesh generations" on one line was caught. Found 2026-09-24 when my own text
        # evaded it by accident. The gate saw LINES; the document is not lines. Each line is
        # now tested joined to the one after it, with whitespace collapsed, so a term
        # straddling the break is caught and still reported on its first line.
        lines = p.read_text().split("\n")
        for n, line in enumerate(lines, 1):
            nxt = lines[n] if n < len(lines) else ""
            low = " ".join((line + " " + nxt).split()).lower()
            for term in FORBIDDEN:
                if term in low:
                    if allowed(rel, line):
                        continue
                    hits.append((rel, n, term, line.strip()[:120]))

    for rel, n, term, text in hits:
        print("MESH NARRATION %s:%d  %r" % (rel, n, term))
        print("    %s" % text)
    print("scanned %d files, %d forbidden terms, %d hit(s), %d pinned exception(s)"
          % (scanned, len(FORBIDDEN), len(hits), len(ALLOWED_CONTEXT)))
    if hits:
        print("\nThis report reports the final numbers and does not narrate how they were "
              "arrived at. Name an excluded case without explaining which grid produced "
              "what.")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
