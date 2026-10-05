#!/usr/bin/env python3
"""check_mesh_caps.py: FAIL a mesh whose achieved value equals one of its own caps.

THE PATTERN THIS GATE EXISTS FOR: A CAP THAT STOPS AND RETURNS A PLAUSIBLE RESULT IS
INDISTINGUISHABLE FROM CONVERGENCE AT THE EXIT CODE. snappyHexMesh hits a limit, writes a
normal-looking summary, and exits 0. Nothing objects. Two instances so far and BOTH WERE
FOUND BY LOOKING RATHER THAN BY BEING TOLD:
  1. maxLocalCells 4000000, which would have silently capped the level-10 grid-family mesh
     and broken the r = 2 refinement ratio the whole GCI rests on, with no error anywhere.
  2. nLayerIter 50, which truncated layer addition on BOTH production meshes while coverage
     was still falling at -0.032 pp/iter, making 74.488% and 73.577% truncation points
     rather than converged values.
The sweep that found both was a one-off, and A ONE-OFF SWEEP PROTECTS ONLY THE MESHES THAT
EXISTED ON THE DAY IT RAN. An instruction in a traps file is prose, and prose passes every
automated check. This is the gate version.

HOW IT DECIDES. For each cap it reads THE VALUE FROM THE DICTIONARY ACTUALLY USED (not a
remembered default) and THE ACHIEVED VALUE FROM THE CORRESPONDING LOG, and fails on
EQUALITY. Equality, not "close to": an iteration counter that stops exactly at its limit did
not choose to stop there. Where a count is not reported in the log the cap is reported as
NOT-TESTABLE and counted as such, never silently as a pass (D052: a sweep that reports only
what it could test overstates its own coverage).

THE OVERRIDE IS DELIBERATE AND REQUIRES A REASON. The grid family legitimately runs at
nLayerIter 50 on all three levels, because converging layers would mean rebuilding the
frozen production point distribution and because converging only two of three levels would
leave the family varying in cell size AND layer-iteration state at once, which makes the
observed order of convergence absorb both and the GCI uninterpretable. That is a decision,
so it must be SAID OUT LOUD in the run card rather than silently tolerated:
    --allow nLayerIter --reason "grid family: consistent cap, see D084"
An override without a reason string is refused.

Usage:
  python3 scripts/check_mesh_caps.py <case_dir> [<case_dir> ...]
        [--allow CAP [--allow CAP ...]] [--reason "why"] [--write-run-card]
"""
import argparse
import gzip
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# THREE CLASSES, READ OUT OF THE OPENFOAM-org 7 SOURCE, NOT INFERRED FROM THE NAMES.
# Every classification below cites the file and line that decides it. An earlier version of
# this file put nRelaxIter in FIXED_COUNT on the strength of its name and was WRONG.
#
# CLASS 1, TRUE CAP: bounds a loop that CAN BREAK EARLY on its own success criterion, so
# reaching the bound means the process was CUT SHORT and the result is not converged.
#   nLayerIter   snappyLayerDriver.C  - layer addition loop, breaks when no more layers move
#   nRelaxIter   medialAxisMeshMover.C:2050-2093
#                    const label nSnap = readLabel(meshQualityDict.lookup("nRelaxIter"));
#                    for (label iter = 0; iter < 2*nSnap; iter++) {
#                        if (iter == nSnap) { ...setErrorReduction(0.0); }
#                        if (meshMover_.scaleMesh(...)) { meshOk = true; break; }
#                    }
#                    return meshOk;
#                THE BOUND IS 2*nRelaxIter, NOT nRelaxIter, and at iter == nRelaxIter the
#                mover CHANGES BEHAVIOUR (error reduction to zero) rather than stopping.
#                Exhausting it returns meshOk FALSE, i.e. the mesh was never successfully
#                moved. That is a truncation and a real finding.
#                snapParameters.C:36 reads the SAME KEY in snapControls into nSnap_, so one
#                dictionary name feeds two different consumers.
#   maxLocalCells / maxGlobalCells  - stop refinement that would otherwise continue
#
# CLASS 2, FIXED COUNT: always executes exactly N times BY DESIGN, so achieved == limit is
# the normal and ONLY possible outcome. Flagging these would fire on EVERY MESH EVER BUILT
# and train the reader to dismiss the gate, which is the same end state as no gate.
#   nSmoothPatch      snappySnapDriver.C:668   smoothIter < snapParams.nSmoothPatch()
#   nSolveIter        snappySnapDriver.C:1719  iter < snapParams.nSmoothDispl()   (no break)
#   nFeatureSnapIter  snappySnapDriver.C:2358  nFeatIter = snapParams.nFeatureSnap()
#   nCellsBetweenLevels / minRefinementCells  - not loop counters at all
#
# CLASS 3, THRESHOLD: not a bound on anything. Crossing it CHANGES BEHAVIOUR mid-run.
#   nRelaxedIter  snappyLayerDriver.C:3157-3162
#                     const dictionary& meshQualityDict =
#                         (iteration < layerParams.nRelaxedIter() ? motionDict
#                                                                : motionDict.subDict("relaxed"));
#                 Above it, snappy switches to the RELAXED meshQuality constraints. Default
#                 is labelMax (layerParameters.C:143), i.e. never relax, and it is ABSENT
#                 from our dictionaries, so we never relax. Reaching it is neither pass nor
#                 fail; it is a change of regime and must be REPORTED, never gated on.
THRESHOLD = {"nRelaxedIter"}

# a TRUE CAP whose effective bound is a multiple of the dictionary value
CAP_MULTIPLIER = {"nRelaxIter": 2}

FIXED_COUNT = {"nFeatureSnapIter", "nSolveIter", "nSmoothPatch",
               "nCellsBetweenLevels", "minRefinementCells"}

# cap name -> (regex for the dict value, list of regexes for achieved values in the log)
# Each achieved regex must capture ONE integer.
CAPS = {
    "nLayerIter":        (r"nLayerIter\s+(\d+)",
                          [r"^Layer addition iteration"]),          # counted, not captured
    "maxLocalCells":     (r"maxLocalCells\s+(\d+)",     [r"Cells:\s*(\d+)"]),
    "maxGlobalCells":    (r"maxGlobalCells\s+(\d+)",    [r"Cells:\s*(\d+)"]),
    "nRelaxIter":        (r"nRelaxIter\s+(\d+)",        []),
    "nSmoothPatch":      (r"nSmoothPatch\s+(\d+)",      []),
    "nSolveIter":        (r"nSolveIter\s+(\d+)",        []),
    "nFeatureSnapIter":  (r"nFeatureSnapIter\s+(\d+)",  [r"Snapping to features in (\d+) iterations"]),
    "minRefinementCells": (r"minRefinementCells\s+(\d+)", []),
    "nCellsBetweenLevels": (r"nCellsBetweenLevels\s+(\d+)", []),
}
# caps whose achieved value is a COUNT OF LOG LINES rather than a captured number
COUNTED = {"nLayerIter": r"^Layer addition iteration"}


def read_text(p: Path):
    if p.suffix == ".gz":
        return gzip.open(p, "rt", errors="ignore").read()
    return p.read_text(errors="ignore")


def achieved(cap, log, dict_txt):
    """(value, how) or (None, reason-not-testable)."""
    if cap in COUNTED:
        n = len(re.findall(COUNTED[cap], log, re.M))
        return (n, f"counted {n} '{COUNTED[cap].strip('^')}' lines") if n else \
               (None, "no layer-addition iterations in log")
    pats = CAPS[cap][1]
    if not pats:
        return None, "no achieved value is reported in the log for this cap"
    best = None
    for pat in pats:
        m = [int(x) for x in re.findall(pat, log)]
        if m:
            best = max(m) if best is None else max(best, max(m))
    return (best, "from log") if best is not None else (None, "pattern not found in log")


def check_case(case: Path, allow=(), reason=None):
    d = case / "system/snappyHexMeshDict"
    logs = [case / "log.snappyHexMesh"]
    log_p = next((p for p in logs if p.exists()), None)
    if not d.exists():
        return dict(case=str(case), status="SKIP", why="no snappyHexMeshDict", rows=[])
    if log_p is None:
        return dict(case=str(case), status="NOT-TESTABLE",
                    why="no log.snappyHexMesh on disk; caps cannot be checked and MUST NOT "
                        "be assumed un-hit", rows=[])
    dt, log = read_text(d), read_text(log_p)
    rows, failed, nottest = [], [], 0
    for cap, (dpat, _) in CAPS.items():
        m = re.search(dpat, dt)
        if not m:
            continue
        limit = int(m.group(1)) * CAP_MULTIPLIER.get(cap, 1)
        if cap in CAP_MULTIPLIER:
            how_mult = f"effective bound {limit} = {CAP_MULTIPLIER[cap]}x dict value"
        got, how = achieved(cap, log, dt)
        if cap in THRESHOLD:
            rows.append((cap, limit, got, "threshold", "regime change, never a gate"))
            continue
        if got is None:
            rows.append((cap, limit, None, "NOT-TESTABLE", how)); nottest += 1
            continue
        if got >= limit:
            if cap in FIXED_COUNT:
                rows.append((cap, limit, got, "info (fixed)", how))
                continue
            state = "ALLOWED" if cap in allow else "AT CAP"
            rows.append((cap, limit, got, state, how))
            if cap not in allow:
                failed.append(cap)
        else:
            rows.append((cap, limit, got, "ok", how))
    status = "FAIL" if failed else "PASS"
    return dict(case=str(case), status=status, rows=rows, failed=failed,
                not_testable=nottest, allowed=list(allow), reason=reason)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--allow", action="append", default=[],
                    help="cap name permitted to sit at its limit; requires --reason")
    ap.add_argument("--reason", default=None)
    ap.add_argument("--write-run-card", action="store_true")
    a = ap.parse_args()
    if a.allow and not a.reason:
        print("check_mesh_caps: --allow requires --reason. A cap that is deliberately hit "
              "must SAY SO OUT LOUD in the run card; silent tolerance is the defect this "
              "gate exists to stop.", file=sys.stderr)
        return 2

    rc = 0
    for c in a.cases:
        r = check_case(Path(c), allow=tuple(a.allow), reason=a.reason)
        print(f"\n{r['case']}: {r['status']}" + (f"  ({r.get('why','')})" if r.get("why") else ""))
        for cap, lim, got, state, how in r["rows"]:
            g = "-" if got is None else str(got)
            flag = {"AT CAP": "  <-- FAIL, achieved == cap",
                    "info (fixed)": "  <-- fixed count, by design, NOT a failure",
                    "ALLOWED": "  <-- at cap, ALLOWED with reason",
                    "threshold": "  <-- THRESHOLD, regime change not failure",
                    "NOT-TESTABLE": "  <-- NOT TESTABLE, not a pass"}.get(state, "")
            print(f"    {cap:22s} cap {lim:<10} achieved {g:<10} {state:13s}{flag}")
        if r["rows"]:
            n_ok = sum(1 for x in r["rows"] if x[3] == "ok")
            print(f"    coverage: {n_ok} tested-ok, {len(r.get('failed',[]))} at cap, "
                  f"{r.get('not_testable',0)} NOT-TESTABLE of {len(r['rows'])} caps present")
        if r["status"] in ("FAIL", "NOT-TESTABLE"):
            rc = 1
        if a.write_run_card and r["rows"]:
            p = Path(c) / "run_card.json"
            card = json.loads(p.read_text()) if p.exists() else {}
            card.setdefault("mesh", {})["caps"] = {
                "checked": [x[0] for x in r["rows"]],
                "at_cap": r.get("failed", []) + [x[0] for x in r["rows"] if x[3] == "ALLOWED"],
                "not_testable": [x[0] for x in r["rows"] if x[3] == "NOT-TESTABLE"],
                "override_allowed": r["allowed"], "override_reason": r["reason"]}
            p.write_text(json.dumps(card, indent=2) + "\n")
            print(f"    run card updated: {p}")
    print(f"\nGATE: {'PASS' if rc == 0 else 'FAIL / NOT-TESTABLE'}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
