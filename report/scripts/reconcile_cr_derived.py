#!/usr/bin/env python3
"""Recompute every Condition CR field that is a FUNCTION of other fields in the same record.

WHY THIS EXISTS. update_cr_from_r2.py patched Cd_counts and Cl into the force files and did
not touch the fields DERIVED from them. The result was a record that contradicted itself:

    H  Cd 188.065  Cl 0.42573121   Cl_minus_target_counts -0.132   trim_within_tolerance true

Cl is 25.46 counts off target and the field beside it said 0.132, and the boolean said the
trim was inside a 1.0-count gate it had never been checked against. Five of six rows
disagreed with themselves. Nothing surfaced it, because every consumer read the derived
field rather than recomputing it (GEO-080: a claim about a value must come from the value,
not from the hand that wrote it).

WHAT IS DERIVED AND WHAT IS NOT.
  DERIVED, recomputed here:
    Cl_minus_target_counts      (Cl - target) * 1e4
    trim_within_tolerance       |Cl_minus_target_counts| <= trim_tolerance_counts
    L_over_D                    Cl / (Cd_counts * 1e-4)
  NOT DERIVED, never touched:
    Cd_counts, Cl, alpha_deg, delta_Cd_counts_vs_baseline, baseline_Cd_counts_used, and
    every provenance field. Those are measurements and pointers to measurements.

IT REPORTS BEFORE IT WRITES, and it writes only on --apply. A dry run that silently fixed
things would be the same class of defect one level up.
"""
import argparse
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
TARGET_CL = 0.428277635108
FILES = ["data/rans_forces.json", "data/rans_forces_current.json"]
ORDER = ["B", "C", "M", "W", "H", "F"]


def derive(e):
    """Return the derived fields this record's own measurements imply."""
    out = {}
    if "Cl" in e:
        out["Cl_minus_target_counts"] = (e["Cl"] - TARGET_CL) * 1e4
    if "Cl" in e and "Cd_counts" in e and e["Cd_counts"]:
        out["L_over_D"] = e["Cl"] / (e["Cd_counts"] * 1e-4)
    # THE TRIM-RESIDUAL BIAS IS DERIVED TOO, and four of six rows carried a stale one.
    # It is the correction that would be ADDED to Cd to move it to the target lift, so it
    # is MINUS the slope times the residual. B stored -0.015729, which is its slope times
    # the +0.516 ct residual that was itself stale; its actual residual is -0.094 ct and
    # the bias is -0.003. Same shape as the other derived fields: patched once, never
    # recomputed, and read by fig/make_vlmerr.py without question.
    slope = e.get("dCd_dCl_counts_per_unit_CL")
    if slope is not None and "Cl" in e:
        out["Cd_bias_from_trim_residual_counts"] = -slope * (e["Cl"] - TARGET_CL)
    tol = e.get("trim_tolerance_counts")
    if tol is not None and "Cl_minus_target_counts" in out:
        out["trim_within_tolerance"] = abs(out["Cl_minus_target_counts"]) <= tol
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--apply", action="store_true",
                    help="write the corrections; without it this only reports")
    a = ap.parse_args()

    total_changed = 0
    for rel in FILES:
        path = HERE / rel
        if not path.exists():
            print("ABSENT %s" % rel)
            continue
        doc = json.loads(path.read_text())
        geoms = doc.get("trims", {}).get("condition_CR", {}).get("geometries", {})
        if not geoms:
            print("%s: no condition_CR trims block" % rel)
            continue
        print("\n=== %s ===" % rel)
        changed = 0
        for L in [x for x in ORDER if x in geoms] + \
                 [x for x in sorted(geoms) if x not in ORDER]:
            e = geoms[L]
            want = derive(e)
            for k, v in want.items():
                have = e.get(k)
                same = (have == v if isinstance(v, bool)
                        else have is not None and abs(float(have) - float(v)) < 1e-9)
                if same:
                    continue
                print("  %s.%-24s %-12s -> %s" % (L, k, str(have), _fmt(v)))
                changed += 1
                if a.apply:
                    e[k] = v
        print("  fields disagreeing with their own record: %d" % changed)
        total_changed += changed
        if a.apply and changed:
            path.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
            print("  WROTE %s" % rel)

    if not a.apply:
        print("\nDRY RUN. %d field(s) would change. Re-run with --apply." % total_changed)
        return 1 if total_changed else 0

    # READ BACK. The file on disk is the thing that matters, not what was intended.
    print("\n=== read back from disk ===")
    bad = 0
    for rel in FILES:
        path = HERE / rel
        if not path.exists():
            continue
        # SAME GUARD AS THE MAIN LOOP. This read-back did not have it and crashed on the
        # file that legitimately carries no trims block, turning a successful write into a
        # traceback. A verification pass that cannot survive the inputs the pass it
        # verifies accepts is not a verification pass.
        geoms = json.loads(path.read_text()).get("trims", {}) \
                    .get("condition_CR", {}).get("geometries", {})
        for L in sorted(geoms):
            for k, v in derive(geoms[L]).items():
                have = geoms[L].get(k)
                ok = (have == v if isinstance(v, bool)
                      else have is not None and abs(float(have) - float(v)) < 1e-9)
                if not ok:
                    print("  STILL WRONG %s %s.%s" % (rel, L, k))
                    bad += 1
    print("  residual disagreements: %d" % bad)
    return 1 if bad else 0


def _fmt(v):
    return str(v) if isinstance(v, bool) else "%.6f" % v


if __name__ == "__main__":
    sys.exit(main())
