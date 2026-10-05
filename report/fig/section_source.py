#!/usr/bin/env python3
"""Resolve, per case, WHICH section dataset supplies it, and keep its q/dragDir row with it.

WHY THIS EXISTS AS ONE MODULE. Two figure scripts (make_cp.py, make_cf.py) need the same
answer to the same question. A rule that lives in two files goes stale in one of them, and
this project has the scar: the same increments table existed in three places and a sweep
found two still carrying superseded numbers after the third was corrected. So the rule
lives here and both import it.

THE RULE. A case is served from data/sections_welded if it exists there, else from
data/sections_r2, else from data/sections_nflag. Project decision, 2026-09-22: show all six
geometries, using the newer evaluation wherever it exists rather than withholding four
because two are not ready. sections_welded (2026-10-05, ARG-196) holds the twelve cruise
trims re-run on the welded meshes (130M-cell cap); it is checked FIRST because only cruise cases
are in it, so it cannot shadow a Condition CR case.

THE TRAP THIS MODULE EXISTS TO PREVENT, and it is measured, not hypothetical. The two
datasets do NOT agree on alpha_deg, and therefore do not agree on dragDir:

    case   nflag alpha    r2 alpha     delta
    SWB    1.742369       1.747609     +0.005240 deg
    SWW    1.490454       1.520545     +0.030091 deg
    SWH    1.522701       1.552310     +0.029609 deg
    SWF    1.446357       1.475456     +0.029099 deg

cf_s is mult * (tau . dragDir) / q, so dragDir enters the skin friction directly. Taking
the section CSV from one dataset and the Q_PINF row from the other applies the wrong wind
axes: exactly the D019 defect, where a 0.0083 deg axis lag biased cd by cl times the lag.
The deltas above are up to 3.6x that. SO THE PAIR TRAVELS TOGETHER OR NOT AT ALL, and
resolve() returns both from one generation rather than letting a caller assemble them.

WHAT IT DOES NOT DO. It does not decide whether a MIXED-GENERATION FIGURE is acceptable;
that is the caller's business and the project has ruled on it. It records what was used, so the
mixing is a stated fact rather than an invisible one, and writes that record to
fig/section_provenance.json for anything downstream that needs to quote it.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)
DATA = os.path.join(REPORT, "data")

WELDED = os.path.join(DATA, "sections_welded")
R2 = os.path.join(DATA, "sections_r2")
PUBLISHED = os.path.join(DATA, "sections_nflag")
# Every dataset, in resolution order. Callers that load all Q_PINF tables iterate THIS rather
# than naming directories, so a new generation cannot be resolved here and missed there.
DATASETS = (WELDED, R2, PUBLISHED)
GENERATION = {WELDED: "welded", R2: "r2", PUBLISHED: "published"}

# Measured bound on what the generation difference does to a surface coefficient, from the
# B and W cases evaluated BOTH ways. Quoted so a caller can state the residual on a
# mixed-generation difference rather than assert it away (D071's caveat: cancellation is
# exact only for a linear differential; for a viscous one, bound it).
GENERATION_DCP_P95 = {"SWB": 0.0126, "SWW": 0.0130}


def _has(dirpath, case):
    """True if this dataset actually carries this case. Tests a station file, not the
    directory: an empty or half-written directory must not count as coverage."""
    return os.path.exists(os.path.join(dirpath, "%s_eta060.csv" % case))


def resolve(case):
    """-> (section_dir, qpi_table_path, generation_tag) for one case, all one generation.

    `case` is the full case name as it appears in filenames, e.g. 'SWB_trim'.
    """
    for d in DATASETS:
        if _has(d, case):
            return d, os.path.join(d, "Q_PINF_TABLE.csv"), GENERATION[d]
    raise KeyError("no section data for %s in any dataset" % case)


def partition(cases):
    """Bucket every case and ASSERT the buckets sum to the input (GEO-089).

    A claim over a set must partition the set: served-from-r2 plus served-from-published
    plus unavailable must equal every case asked about. Anything in none of them is a
    failure, not a skip.
    """
    out = {"welded": [], "r2": [], "published": [], "unavailable": []}
    for c in cases:
        try:
            out[resolve(c)[2]].append(c)
        except KeyError:
            out["unavailable"].append(c)
    total = sum(len(v) for v in out.values())
    if total != len(cases):
        raise AssertionError("partition lost cases: %d bucketed, %d asked" % (total, len(cases)))
    return out


def write_provenance(cases, path=None):
    """Record which generation served each case, so the mixing is a stated fact.

    THIS IS THE RECORD, NOT THE REPORT. The report does not narrate how the evaluations
    were produced; this file is where a later reader, or a script, establishes it.
    """
    p = partition(cases)
    doc = {
        "purpose": "which section dataset served each case in the Cp/Cf figures",
        "rule": "sections_welded, else sections_r2, else sections_nflag, first that carries the case",
        "served_from_welded": sorted(p["welded"]),
        "served_from_r2": sorted(p["r2"]),
        "served_from_published": sorted(p["published"]),
        "unavailable": sorted(p["unavailable"]),
        "counts": {k: len(v) for k, v in p.items()},
        "mixed": sum(1 for k in ("welded", "r2", "published") if p[k]) > 1,
        "residual_bound_note": (
            "Where a difference is taken between cases served from different generations, "
            "the generation itself differs between the two sides and does not cancel "
            "(D071). Measured bound on that residual, from the cases evaluated both ways: "
            "p95 |dCp| " + ", ".join("%s %.4f" % (k, v)
                                     for k, v in sorted(GENERATION_DCP_P95.items())) + "."
        ),
        "axes_note": (
            "alpha_deg and dragDir differ between the datasets by up to 0.0301 deg, so each "
            "case's Q_PINF row is taken from the SAME dataset as its section CSV. Mixing "
            "them would reproduce the D019 wind-axis bias."
        ),
    }
    out = path or os.path.join(HERE, "section_provenance.json")
    with open(out, "w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return doc


if __name__ == "__main__":
    CR = ["SWB_trim", "SWC_trim", "SWF_trim", "SWH_trim", "SWM_trim", "SWW_trim"]
    EC = ["CMPB_trim", "CMPC_trim", "CMPF_trim", "CMPH_trim", "CMPM_trim", "CMPW_trim"]
    LC = ["LCB_trim", "LCC_trim", "LCF_trim", "LCH_trim", "LCM_trim", "LCW_trim"]
    d = write_provenance(CR + EC + LC)
    print("served from welded    : %s" % ", ".join(d["served_from_welded"]))
    print("served from r2        : %s" % ", ".join(d["served_from_r2"]))
    print("served from published : %s" % ", ".join(d["served_from_published"]))
    print("unavailable           : %s" % (", ".join(d["unavailable"]) or "none"))
    print("counts                : %s   mixed=%s" % (d["counts"], d["mixed"]))
