#!/usr/bin/env python3
"""prune_legacy_fo.py: remove function-object blocks that argusPostPro supersedes.

    python3 scripts/prune_legacy_fo.py --self-test
    python3 scripts/prune_legacy_fo.py <controlDict> [<controlDict> ...]
    python3 scripts/prune_legacy_fo.py --dry-run <controlDict>

WHY THIS EXISTS. The L11 recipes define their own function objects, and
install_postpro.py then inserts `#include "argusPostPro"` WITHOUT removing them. It is
additive by design, which is safe in isolation and catastrophic in aggregate: every
output is sampled and written TWICE. The wing surface is 510 MB a write, so campaign-wide
that was 121 GB of byte-for-byte duplicate carrying zero information, on a filesystem
whose operator was emailing us about being full.

WHAT IS REMOVED, AND WHAT IS NOT
    wingSurface       -> argusWingSurface
    wallShearStress   -> argusWallShearStress
    yPlus             -> argusYPlus
    forces1           -> argusForces

`forceCoeffs1` IS NEVER REMOVED AND CANNOT BE. argusPostPro HAS NO EQUIVALENT, and six
scripts read postProcessing/forceCoeffs1: splice_postpro.py, plot_convergence.py,
at_target_cl.py, trim_from_bracket.py, stationarity.py and collect_postpro.sh. Removing
it would break the convergence gate, the trim and every quotable number in the project.
It is on a REFUSE list rather than merely absent from the remove list, because "not
listed" is a fact about today's list and "refused" is a fact about the code.

BRACE-MATCHED, NOT REGEX-TO-THE-NEXT-BRACE. A function-object body contains nested braces
(`surfaces ( wing { ... } )`), so a lazy match to the first `}` truncates the block and
leaves a syntactically broken dictionary that OpenFOAM reports thirty lines later as a
mismatched brace somewhere else entirely.
"""
import argparse
import re
import sys
from pathlib import Path

# argusPostPro supersedes these. Order is irrelevant; presence is not assumed.
SUPERSEDED = ("wingSurface", "wallShearStress", "yPlus", "forces1")

# Never removable, whatever a caller asks for.
REFUSE = ("forceCoeffs1",)


class PruneError(RuntimeError):
    pass


def remove_block(text, name):
    """Remove a top-level `name { ... }` block from a functions dictionary.

    Returns (new_text, removed_bool). Brace-matched from the block's opening brace, so
    nested dictionaries and inline `{ }` inside lists are handled.
    """
    if name in REFUSE:
        raise PruneError(
            "%s is on the REFUSE list: argusPostPro has no equivalent and six scripts "
            "read postProcessing/%s. Removing it breaks the convergence gate, the trim "
            "and every quotable number." % (name, name))

    m = re.search(r"^([ \t]*)%s[ \t]*\r?\n[ \t]*\{" % re.escape(name), text, re.M)
    if not m:
        # also accept `name {` on one line
        m = re.search(r"^([ \t]*)%s[ \t]*\{" % re.escape(name), text, re.M)
        if not m:
            return text, False

    open_idx = text.index("{", m.start())
    depth, i = 0, open_idx
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    if depth != 0:
        raise PruneError("unbalanced braces while removing %s" % name)

    end = i + 1
    while end < len(text) and text[end] in " \t":
        end += 1
    if end < len(text) and text[end] == "\n":
        end += 1
    return text[:m.start()] + text[end:], True


def prune(path, dry_run=False):
    p = Path(path)
    txt = p.read_text()
    removed = []
    for name in SUPERSEDED:
        txt, did = remove_block(txt, name)
        if did:
            removed.append(name)

    # DERIVED CHECK, not an assumption: forceCoeffs1 must still be there afterwards.
    if "forceCoeffs1" not in txt and "forceCoeffs1" in p.read_text():
        raise PruneError("%s: forceCoeffs1 disappeared during pruning" % path)

    if removed and not dry_run:
        p.write_text(txt)
    return removed, txt


def self_test():
    ok = True

    def check(name, cond, note=""):
        nonlocal ok
        print("  [%s] %s%s" % ("ok " if cond else "FAIL", name,
                               ("  " + note) if note else ""))
        ok = ok and cond

    # A block whose body contains NESTED braces, which is the case that breaks a lazy regex.
    src = """functions
{
    forces1
    {
        type forces;
    }
    forceCoeffs1
    {
        type forceCoeffs;
    }
    wingSurface
    {
        type surfaces;
        surfaces ( wing { type patch ; patches ( wing ) ; } );
    }
    keepMe
    {
        type something;
    }
}
"""
    out, did = remove_block(src, "wingSurface")
    check("nested-brace block removed whole", did and "surfaces (" not in out)
    check("the block AFTER it survives", "keepMe" in out and "type something;" in out)
    check("braces still balanced", out.count("{") == out.count("}"),
          "%d open, %d close" % (out.count("{"), out.count("}")))

    out2, did2 = remove_block(src, "forces1")
    check("forces1 removed", did2 and "type forces;" not in out2)
    check("forceCoeffs1 NOT removed by the forces1 pass",
          "forceCoeffs1" in out2 and "type forceCoeffs;" in out2)

    out3, did3 = remove_block(src, "notPresent")
    check("absent name is a no-op", (not did3) and out3 == src)

    try:
        remove_block(src, "forceCoeffs1")
        check("forceCoeffs1 is REFUSED", False, "it was allowed")
    except PruneError:
        check("forceCoeffs1 is REFUSED", True)

    # removing every superseded name leaves exactly forceCoeffs1 and keepMe
    t = src
    for n in SUPERSEDED:
        t, _ = remove_block(t, n)
    check("only forceCoeffs1 and keepMe remain",
          "forceCoeffs1" in t and "keepMe" in t and "wingSurface" not in t
          and "forces1" not in t)
    check("still balanced after all removals", t.count("{") == t.count("}"))

    print("  self-test: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("controlDicts", nargs="*")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not a.controlDicts:
        ap.error("give at least one controlDict, or --self-test")
    rc = 0
    for c in a.controlDicts:
        try:
            removed, _ = prune(c, a.dry_run)
            print("  %s: %s%s" % (c, ", ".join(removed) if removed else "nothing to remove",
                                  "  (dry run)" if a.dry_run and removed else ""))
        except (PruneError, OSError) as e:
            print("  %s: HALT %s" % (c, e))
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
