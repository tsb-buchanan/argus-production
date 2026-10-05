#!/usr/bin/env python3
"""Resolve which derived record supplies each (condition, geometry), and say which it was.

WHY THIS EXISTS. The hinge-moment and root-bending figures build their record name as
"<prefix><letter>_trim.json", which at Condition CR is SWB_trim.json and friends: the
PUBLISHED records, taken from the published meshes at the published trim alphas. Condition
CR has since been re-run on the wake-refined meshes, and those records are SWB_r2_trim.json
and friends. So the two figures were drawing published-mesh loads beside wake-refined drag
and span efficiency, with nothing saying so.

THE REMEDY IS A RESOLVER, NOT AN OVERWRITE. Writing the r2 numbers into SWB_trim.json would
destroy the published record, which is the thing every superseded figure and table is still
checked against. The published records stay exactly as they are; this only chooses.

THE RULE. For a prefix whose condition was re-run, prefer the re-run record when it exists,
else fall back to <prefix><L>_trim.json. Condition CR (SW) was re-run on the wake-refined
meshes, records <prefix><L>_r2_trim.json. Early (CMP) and late (LC) cruise were re-run on
the welded meshes (2026-10-05, ARG-196), records <prefix><L>_welded.json.

IT RECORDS WHAT IT CHOSE. A figure that silently mixes generations is the failure this
module exists to prevent, so resolve() returns the generation tag with the path and
write_provenance() puts the whole mapping on disk. The frame travels with the number.

MEASURED DIFFERENCE, so a caller can state the residual rather than assume it away. Between
the published and wake-refined records at Condition CR:
    hinge moment M_h total   +0.0036 to +0.0180 N m   (0.4 to 2.0 per cent)
    root bending RBM total   -0.139 to -0.388 N m     (0.1 to 0.27 per cent)
and the split is itself informative: F, H and W, whose meshes were built from unwelded
surfaces and lost shell refinement, move about FOUR TIMES as far as B, C and M, which did
not. That is the seam defect showing up in a surface integral.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)

# THE DERIVED DIRECTORY IS SEARCHED FOR, NOT COUNTED TO. The first version computed it as
# three directory levels above the report directory, which is true in the working repo and
# WRONG everywhere else: in the delivery layout it resolves to /home/results/derived, i.e.
# outside the repository entirely. Here that path does not exist so the guard fires, but on
# a machine where /home/results/derived happened to exist the figures would have read a
# stranger's data and said nothing. Counting levels is a guess about the tree above you.
#
# THE SENTINEL IS A FILE, NOT THE DIRECTORY, for the same reason as report_paths.py: an
# empty or half-copied results/derived would satisfy a directory test and the caller would
# then fail later with a confusing missing-key error instead of a clear missing-root one.
#
# RESOLUTION IS LAZY. Raising at import time would break anything that merely imports this
# module to read REMESHED_PREFIXES, including in a repo where the derived records are
# deliberately not shipped.
SENTINEL = os.path.join("results", "derived", "SWB_trim.json")


def derived_dir():
    """Nearest ancestor carrying results/derived/<baseline record>. Raises if there is none."""
    d = HERE
    while True:
        if os.path.isfile(os.path.join(d, SENTINEL)):
            return os.path.join(d, "results", "derived")
        parent = os.path.dirname(d)
        if parent == d:
            raise SystemExit(
                "HALT: cannot locate results/derived. Looked for %r in %s and every parent. "
                "The delivery repository does not ship results/derived, so these figures "
                "regenerate only in a tree that carries it." % (SENTINEL, HERE))
        d = parent


# Conditions whose records were re-run, and the (record suffix, generation tag) of the re-run.
REMESHED = {"SW": ("_r2_trim", "r2"), "CMP": ("_welded", "welded"), "LC": ("_welded", "welded")}
REMESHED_PREFIXES = set(REMESHED)


def _path(name):
    return os.path.join(derived_dir(), name + ".json")


def resolve(prefix, letter):
    """-> (path, generation) for one (condition, geometry). generation is r2, welded or published.

    Falls back rather than failing: a caller that wants the record to exist should check
    the returned path, because "which generation" and "does it exist" are different
    questions and conflating them is how a missing case becomes an invisible one.
    """
    if prefix in REMESHED:
        suffix, gen = REMESHED[prefix]
        p = _path("%s%s%s" % (prefix, letter, suffix))
        if os.path.exists(p):
            return p, gen
    return _path("%s%s_trim" % (prefix, letter)), "published"


def write_provenance(pairs, path=None):
    """Record which generation served each (prefix, letter). pairs: [(prefix, letter), ...]"""
    served = {}
    for prefix, letter in pairs:
        p, gen = resolve(prefix, letter)
        served.setdefault(gen, []).append("%s%s" % (prefix, letter))
    doc = {
        "purpose": "which derived record served each condition and geometry",
        "rule": "the re-run record (SW: <prefix><L>_r2_trim; CMP, LC: <prefix><L>_welded) where it exists, else <prefix><L>_trim",
        "remeshed_prefixes": sorted(REMESHED_PREFIXES),
        "served": {k: sorted(v) for k, v in served.items()},
        "counts": {k: len(v) for k, v in served.items()},
        "mixed": len(served) > 1,
        "note": (
            "Condition CR is on the wake-refined (r2) meshes and both cruise conditions on "
            "the welded meshes. Within one condition every geometry shares a generation, so "
            "a delta against B does not cross one; a table setting Condition CR beside "
            "cruise crosses a generation as well as an operating point."
        ),
    }
    out = path or os.path.join(HERE, "derived_provenance.json")
    with open(out, "w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return doc


if __name__ == "__main__":
    LETTERS = ["B", "C", "F", "H", "M", "W"]
    d = write_provenance([(p, L) for p in ("SW", "CMP", "LC") for L in LETTERS])
    for gen, names in sorted(d["served"].items()):
        print("  %-10s %s" % (gen, ", ".join(names)))
    print("  counts: %s   mixed=%s" % (d["counts"], d["mixed"]))
