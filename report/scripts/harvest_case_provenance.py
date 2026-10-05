#!/usr/bin/env python3
"""Derive the per-case provenance pass that build_rans_forces.py takes as argv[2].

WHY THIS EXISTS. build_rans_forces.py has always taken TWO inputs: the raw force
harvest (scripts/harvest_rans_forces.py, run on the cluster) and a provenance
dict. The provenance pass was never committed, so when the campaign grew from 51
to 75 converged cases the structured file could not be rebuilt. This script is
that missing half, written so the pass is repeatable by someone who was not here.

EVERY FIELD IS READ BACK FROM A FILE IN THE CASE DIRECTORY. Nothing is authored
beside the action that was supposed to produce it (GEO-080 item 2). The path each
value came from travels with it, so a reader can re-derive any single number
without trusting this script.

SOURCE OF TRUTH IS THE ARCHIVE, NOT SCRATCH. The default root is the U: case
archive. That is deliberate: if the report's data file can be rebuilt from the
archive alone, the archive is demonstrably sufficient to carry the campaign
forward, which is a property worth establishing rather than assuming. The root
actually read is recorded in the output.

STANDING RULES HONOURED HERE, AND WHERE:

1. GEO-089 PARTITION. Every enumerated case directory lands in exactly ONE named
   bucket for each field group, and the buckets are ASSERTED to sum to the
   enumerated count. A case that yields nothing is a FAILURE bucket, never a skip.
2. GEO-087 A GATE MUST ASSERT IT SAW EVERYTHING IT CLAIMS TO HAVE CHECKED. The
   enumerated count and the adjudicated count are asserted equal, not merely
   printed side by side.
3. GEO-080 item 6 VALIDATE BEFORE WRITE. The output goes to a temporary path, is
   read back and checked, and is promoted to its real filename only on pass.
4. D052 REPORT THE UNTESTED SET AS A COUNTED, NAMED BUCKET. Cases with no
   CASE_PROVENANCE.json are named individually with the consequence stated, not
   dropped.
5. GEO-102 A GATE ON THE OUTPUT CANNOT SEE A WRONG INPUT. The surface-closure and
   force-patch facts are INPUT facts no convergence gate can test, so they are
   carried here as a counted partition for build_rans_forces.py to report rather
   than being hardcoded there.

Usage:
    python3 harvest_case_provenance.py [CASE_ROOT] [-o OUT.json]
"""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOT = "/mnt/u/argus_hpc12_backup/cases"
DEFAULT_OUT = os.path.join(os.path.dirname(HERE), "data", "case_provenance_current.json")

# RAS { model <name>; } in constant/momentumTransport. OpenFOAM-org 12 spells the
# key "model"; org 7 spelled it "RASModel". Both are accepted and which one was
# found is recorded, because a silent miss here would null the turbulence model
# for a whole condition and assert_frames would not catch it.
RE_RAS_MODEL = re.compile(r"^\s*(model|RASModel)\s+([A-Za-z][A-Za-z0-9]*)\s*;", re.M)
RE_SIMTYPE = re.compile(r"^\s*simulationType\s+([A-Za-z][A-Za-z0-9]*)\s*;", re.M)
# The controlDict carries a human comment naming the alpha the direction vectors
# were built from. It is frequently STALE: it is not what the solver used. It is
# harvested precisely so build_rans_forces.py can flag the disagreement.
RE_CD_ALPHA = re.compile(r"alpha\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*deg")
# forceCoeffs1 { ... patches (wing); ... }: the patch set the forces integrate.
RE_PATCHES = re.compile(r"^\s*patches\s*\(([^)]*)\)\s*;", re.M)
# Case names are <condition prefix><family letter>_<suffix>: SWB_trim,
# CMPC_a1p20, LCM_trim. SW is Condition CR, CMP early cruise, LC late cruise.
RE_CASE_NAME = re.compile(r"^(SW|CMP|LC)([BCFHMW])_(.+)$")
CONDITION_PREFIX = {"SW": "condition_CR", "CMP": "early_cruise", "LC": "late_cruise"}


def classify_tag(case, tag):
    """Say WHAT a CASE_PROVENANCE tag is, not merely that it looks odd.

    build_rans_forces.py flags a tag that is not a prefix of the case name. That
    test is factually right and its obvious reading is wrong: the tag is the
    GEOMETRY FAMILY tag, minted under the Condition-CR ('SW') naming and carried
    unchanged into the cruise recipes that were cloned from it. So SWM on
    LCM_trim names the right geometry, not the wrong case.

    A wrong reason attached to right numbers is what propagates (GEO-095), so
    the classification travels with the flag rather than being left to a reader.
    """
    m = RE_CASE_NAME.match(case)
    if tag is None:
        return "NO TAG FIELD"
    if case.startswith(tag):
        return "prefix of the case name"
    if m and tag == "SW" + m.group(2):
        return ("geometry family tag (SW + family letter) carried from the "
                "Condition-CR naming into a cloned cruise recipe; names the "
                "correct geometry family for this case")
    return "UNEXPLAINED: not a case prefix and not SW+family letter"


def read(path):
    with open(path, errors="replace") as fh:
        return fh.read()


def harvest_one(case_dir, case):
    """Return (record, notes). Every value carries the path it was read from."""
    rec = {}
    notes = {}

    # ---------------------------------------------------------- CASE_PROVENANCE
    ppath = os.path.join(case_dir, "CASE_PROVENANCE.json")
    if not os.path.isfile(ppath):
        notes["provenance"] = "ABSENT"
    else:
        try:
            p = json.loads(read(ppath))
        except ValueError as exc:
            notes["provenance"] = "UNPARSEABLE: %s" % exc
        else:
            notes["provenance"] = "OK"
            rec["provenance_path"] = ppath
            rec["prov_tag"] = p.get("tag")
            rec["prov_tag_classification"] = classify_tag(case, p.get("tag"))
            rec["prov_alpha_deg_in_case"] = p.get("alpha_deg_in_case")
            geo = p.get("geometry") or {}
            rec["prov_geometry_sha256"] = geo.get("sha256")
            rec["prov_geometry_source"] = geo.get("source")
            # GEO-102 input facts. Carried, not adjudicated here.
            rec["prov_geometry_triangles"] = geo.get("triangles")
            rec["prov_geometry_open_edges"] = geo.get("open_edges")
            rec["prov_geometry_over_edges"] = geo.get("over_edges")
            rec["prov_gates_passed"] = p.get("_gates_passed")
            rcp = p.get("recipe") or {}
            rec["prov_recipe"] = rcp.get("folder")
            rec["prov_recipe_sha"] = rcp.get("manifest_sha256")
            rec["prov_recipe_file_count"] = rcp.get("files")
            rec["prov_speed_m_s"] = p.get("speed_m_s")
            rec["prov_axis_lag_deg"] = p.get("axis_lag_deg")
            rec["prov_case_field"] = p.get("case")

    # -------------------------------------------------------- momentumTransport
    mpath = os.path.join(case_dir, "constant", "momentumTransport")
    if not os.path.isfile(mpath):
        notes["turbulence"] = "ABSENT"
    else:
        txt = read(mpath)
        m = RE_RAS_MODEL.search(txt)
        s = RE_SIMTYPE.search(txt)
        if m is None:
            notes["turbulence"] = "NO MODEL KEY MATCHED"
        else:
            notes["turbulence"] = "OK"
            rec["turbulence_model"] = m.group(2)
            rec["turbulence_model_key"] = m.group(1)
            rec["momentumTransport_path"] = mpath
        if s is not None:
            rec["simulationType"] = s.group(1)

    # -------------------------------------------------------------- controlDict
    cpath = os.path.join(case_dir, "system", "controlDict")
    if not os.path.isfile(cpath):
        notes["controlDict"] = "ABSENT"
    else:
        txt = read(cpath)
        notes["controlDict"] = "OK"
        rec["controlDict_path"] = cpath
        a = RE_CD_ALPHA.search(txt)
        if a is not None:
            rec["controlDict_comment_alpha"] = float(a.group(1))
        pt = RE_PATCHES.search(txt)
        if pt is not None:
            rec["force_patches"] = [w for w in pt.group(1).split() if w]

    # ------------------------------------------------------------ recipe name
    rpath = os.path.join(case_dir, "RECIPE_FILES.json")
    if os.path.isfile(rpath):
        try:
            r = json.loads(read(rpath))
        except ValueError:
            notes["recipe_files"] = "UNPARSEABLE"
        else:
            notes["recipe_files"] = "OK"
            rec["recipe_name"] = r.get("name")
            rec["recipe_files_listed"] = len(r.get("files") or [])
            rec["recipe_files_path"] = rpath
            # Does the recipe DECLARE a CONDITION.json that the archive lacks?
            # This is open gap (2) in the handover, derived rather than asserted.
            rec["recipe_declares_CONDITION_json"] = "CONDITION.json" in (r.get("files") or [])
    else:
        notes["recipe_files"] = "ABSENT"

    rec["archive_has_CONDITION_json"] = os.path.isfile(
        os.path.join(case_dir, "CONDITION.json"))
    return rec, notes


def main():
    argv = [a for a in sys.argv[1:]]
    out_path = DEFAULT_OUT
    if "-o" in argv:
        i = argv.index("-o")
        out_path = argv[i + 1]
        del argv[i:i + 2]
    root = argv[0] if argv else DEFAULT_ROOT

    if not os.path.isdir(root):
        raise SystemExit("case root not found: %s" % root)

    cases = sorted(d for d in os.listdir(root)
                   if os.path.isdir(os.path.join(root, d)) and not d.startswith("_"))
    enumerated = len(cases)

    prov = {}
    # Buckets. Every case lands in exactly one bucket per group (GEO-089).
    buckets = {
        "provenance": {"OK": [], "ABSENT": [], "UNPARSEABLE": []},
        "turbulence": {"OK": [], "ABSENT": [], "NO_MODEL_KEY": []},
        "controlDict": {"OK": [], "ABSENT": []},
        "recipe_files": {"OK": [], "ABSENT": [], "UNPARSEABLE": []},
    }
    for case in cases:
        rec, notes = harvest_one(os.path.join(root, case), case)
        prov[case] = rec
        for group, note in notes.items():
            if note == "OK":
                buckets[group]["OK"].append(case)
            elif note == "ABSENT":
                buckets[group]["ABSENT"].append(case)
            elif note.startswith("UNPARSEABLE"):
                buckets[group]["UNPARSEABLE"].append(case)
            else:
                buckets[group]["NO_MODEL_KEY"].append(case)
        # A group whose file is absent never produced a note for the other
        # groups, so assert every group saw this case exactly once.
        for group in buckets:
            seen = sum(1 for b in buckets[group].values() if case in b)
            if seen != 1:
                raise SystemExit(
                    "PARTITION VIOLATION: %s lands in %d buckets of group %s"
                    % (case, seen, group))

    coverage = {}
    for group, b in buckets.items():
        counts = {k: len(v) for k, v in b.items()}
        total = sum(counts.values())
        coverage[group] = {
            "counts": counts,
            "sum": total,
            "enumerated": enumerated,
            "sums_to_enumerated": total == enumerated,
            "named": {k: sorted(v) for k, v in b.items() if v and k != "OK"},
        }
        if total != enumerated:
            raise SystemExit(
                "COVERAGE VIOLATION in %s: adjudicated %d against enumerated %d"
                % (group, total, enumerated))

    # ---- the tag question, CLASSIFIED rather than merely flagged
    tag_classes = {}
    for c, r in prov.items():
        if "prov_tag_classification" not in r:
            continue
        tag_classes.setdefault(r["prov_tag_classification"], []).append(c)
    mislabelled = sorted(
        "%s (tag %s)" % (c, r["prov_tag"]) for c, r in prov.items()
        if r.get("prov_tag") and not c.startswith(r["prov_tag"]))
    unexplained_tags = sorted(
        c for c, r in prov.items()
        if str(r.get("prov_tag_classification", "")).startswith("UNEXPLAINED"))

    # ---- GEO-102 input-gate coverage, a counted partition over the enumerated set
    closed, unestablished = [], []
    for c, r in prov.items():
        if r.get("prov_geometry_open_edges") == 0 and r.get("prov_geometry_over_edges") == 0:
            closed.append(c)
        else:
            unestablished.append(c)
    patch_sets = {}
    for c, r in prov.items():
        key = ",".join(r.get("force_patches") or ["NOT FOUND"])
        patch_sets.setdefault(key, []).append(c)

    out = {
        "_what": "Per-case provenance and solver metadata for every archived ARGUS "
                 "OF-12 case, read back from the case's own files. This is argv[2] "
                 "of build_rans_forces.py.",
        "_generator": "docs/report/all_geometry_2026-09-15/scripts/harvest_case_provenance.py",
        "_case_root": root,
        "_root_note": "Paths in this file point at the root actually read. Where that "
                      "is the U: archive rather than HPC12 scratch, the path recorded "
                      "is the archive path, because that is the file the value was "
                      "read from. The raw force harvest carries scratch paths for the "
                      "same cases; the two path families in the merged output are not "
                      "an inconsistency, they are each value naming its own source.",
        "_enumerated_case_dirs": enumerated,
        "_coverage": coverage,
        "_provenance_tag": {
            "_what": "build_rans_forces.py flags a CASE_PROVENANCE 'tag' that is not a "
                     "prefix of the case name. That test is right and its obvious reading "
                     "is wrong, so the DIAGNOSIS is carried here beside the flag: the tag "
                     "is the GEOMETRY FAMILY tag, minted under the Condition-CR 'SW' "
                     "naming and carried unchanged into the cruise recipes cloned from "
                     "it. SWM on LCM_trim names the right geometry, not the wrong case. "
                     "It affects no number: alpha, geometry sha256 and recipe in the same "
                     "file are what every quoted value uses, and they are per-case.",
            "classification_counts": {k: len(v) for k, v in sorted(tag_classes.items())},
            "n_not_a_case_prefix": len(mislabelled),
            "not_a_case_prefix": mislabelled,
            "UNEXPLAINED": {
                "_rule": "A tag that is neither a case prefix nor SW+family letter is a "
                         "FAILURE, not a skip. Zero is the expected count.",
                "count": len(unexplained_tags),
                "cases": unexplained_tags,
            },
        },
        "_input_gate_facts_GEO102": {
            "_what": "Facts about what the solver was HANDED, which no convergence gate "
                     "can test. Counted here so build_rans_forces.py reports a derived "
                     "partition instead of a hardcoded one.",
            "surface_closed_open_0_over_0": {"count": len(closed)},
            "surface_closure_UNESTABLISHED": {
                "count": len(unestablished),
                "cases": sorted(unestablished),
                "reason": "no CASE_PROVENANCE.json in the archived case, so the "
                          "surface-closure gate result cannot be read back",
            },
            "sum": len(closed) + len(unestablished),
            "sums_to_enumerated": len(closed) + len(unestablished) == enumerated,
            "force_patches_by_case_count": {k: len(v) for k, v in sorted(patch_sets.items())},
            "force_patches_sum": sum(len(v) for v in patch_sets.values()),
        },
        "cases": prov,
    }

    tmp = out_path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    # validate before promoting (GEO-080 item 6)
    back = json.loads(read(tmp))
    assert len(back["cases"]) == enumerated, "case count changed on round trip"
    for group, cov in back["_coverage"].items():
        assert cov["sums_to_enumerated"] is True, "coverage broken in %s" % group
    assert back["_input_gate_facts_GEO102"]["sums_to_enumerated"] is True
    os.rename(tmp, out_path)

    print("wrote %s (%d bytes, %d cases from %s)"
          % (out_path, os.path.getsize(out_path), enumerated, root))
    for group, cov in sorted(coverage.items()):
        print("  %-14s %s" % (group, cov["counts"]))
        for k, v in sorted(cov["named"].items()):
            print("      %s: %s" % (k, ", ".join(v)))
    print("  tag classification: %s"
          % {k: len(v) for k, v in sorted(tag_classes.items())})
    if unexplained_tags:
        print("      UNEXPLAINED TAGS (expected 0): %s" % ", ".join(unexplained_tags))
    print("  surface closure: %d closed, %d UNESTABLISHED (%s)"
          % (len(closed), len(unestablished), ", ".join(sorted(unestablished)) or "none"))


if __name__ == "__main__":
    main()
