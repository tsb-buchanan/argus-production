#!/usr/bin/env python3
"""audit_3d_campaign.py: every 3D case, its state, why it stopped, and what survives.

REQUIREMENT, 2026-08-12: a matrix of every planned 3D case, which ones failed and why, what
needs rerunning, and what results exist.

EVERY CASE LANDS IN EXACTLY ONE NAMED BUCKET AND THE BUCKETS SUM TO THE TOTAL (GEO-089). A
campaign audit that reports only the cases it could classify overstates its own coverage; the
unclassifiable ones are a bucket, not an omission.

THE STOP REASON IS READ FROM THE RUN, NEVER INFERRED FROM THE ITERATION COUNT. "Stopped at 1009
of 4000" is consistent with a time limit, a crash, and a cancel, and those need completely
different responses:
    TIME LIMIT   the recipe is fine, the allocation was wrong    -> resubmit bigger/longer
    FPE / FATAL  the solve diverged                              -> fix numerics or mesh
    COMPLETED    reached endTime                                 -> check convergence separately
Reaching endTime is NOT convergence and being killed is NOT divergence (PROJECT_STATE 0.5.4).

MESH LOCATION MATTERS FOR THE RERUN COST, so it is reported. After `snappy -parallel -overwrite`
the fine mesh exists ONLY in processor*/constant/polyMesh; constant/polyMesh still holds the
background blockMesh. A case whose fine mesh is decomposed N ways cannot simply be resubmitted
at M ranks -- that needs reconstructParMesh or a full remesh, and the cost difference is hours.
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]

REMOTE_SCAN = r'''
for root in /scratch/$USER/argus/phaseA /scratch/$USER/argus/layers \
            /scratch/$USER/argus/prod /scratch/$USER/argus/cases \
            /scratch/$USER/argus/campaign_v4/cases /scratch/$USER/argus/gridconv \
            /scratch/$USER/argus/legB_wallresolved /scratch/$USER/argus/trefftz; do
  [ -d "$root" ] || continue
  for d in "$root"/*/; do
    [ -d "$d" ] || continue
    [ -d "$d/system" ] || continue
    n=$(basename "$d")
    np=$(ls -d "$d"/processor* 2>/dev/null | wc -l)
    # NEWEST log generation THAT ACTUALLY CONTAINS A SOLVER LOG.
    # Taking the newest directory unconditionally reported L11_base/L11_morph as MESH-ONLY:
    # a resubmission on 2026-08-11 created an EMPTY log/ and died, so the newest generation
    # holds nothing while the generation before it holds 1009 healthy iterations ending on a
    # SLURM time limit. "No log in the newest directory" and "never solved" are different
    # facts and must not collapse onto the same label.
    g=""; for c in $(ls -dt "$d"/log "$d"/log.prev.* 2>/dev/null); do
      [ -f "$c/simpleFoam.log" ] && { g="$c"; break; }
    done
    [ -n "$g" ] || g=$(ls -dt "$d"/log "$d"/log.prev.* 2>/dev/null | head -1)
    gen=$(basename "${g:-none}")
    sf="$g/simpleFoam.log"
    # INPUT GATE: DOES THIS LOG BELONG TO THIS CASE? (GEO-102)
    #
    # A solver log states its own case on its "Case :" line. Three phaseA cases
    # (A1_er133, A1b_n18er125, A2_relax) each held a BYTE-IDENTICAL simpleFoam.log,
    # sha256 b0a02b59, whose Case line reads prod/WR_L9_n23 -- a fourth case entirely.
    # One crash, copied into three unrelated directories. The audit read all three and
    # reported "CRASHED (fpe) 3", and PROJECT_STATE recorded "THREE DIVERGENCES ... one
    # cause, not three". It was one EVENT, misattributed, and all three of those cases
    # have never solved at all: processor0 holds only time 0.
    #
    # A gate on the output cannot see a wrong input. So check the input: if the log names
    # a different case, it is FOREIGN and tells us nothing about this one.
    logcase=$(grep -m1 '^Case  *:' "$sf" 2>/dev/null | sed 's/.*:  *//' | sed 's|/*$||')
    foreign=0
    if [ -n "$logcase" ] && [ "$(basename "$logcase")" != "$(basename "$d")" ]; then
      foreign=1
    fi
    it=$(grep -oP '^Time = \K[0-9]+' "$sf" 2>/dev/null | tail -1)
    endt=$(grep -oP '^endTime\s+\K[0-9]+' "$d"/system/controlDict 2>/dev/null | head -1)
    # STOP REASON, read from the tail of the solver log and the slurm output
    reason=""
    if [ -f "$sf" ]; then
      t=$(tail -40 "$sf" 2>/dev/null)
      # ORDER MATTERS: the most specific failure wins. And "End" is matched ANCHORED,
      # because the bare word appears in plenty of benign lines -- an unanchored match
      # labelled three heap-aborted CTE cases "completed" at ~2330 of 5000 iterations,
      # which is exactly the failure ARG-104 recorded. A crash relabelled as a
      # completion is the worst error this table can make.
      case "$t" in
        *"TIME LIMIT"*)               reason=timelimit ;;
        *"double free"*|*"corruption"*|*"glibc"*) reason=heap-abort ;;
        *"FOAM FATAL"*)               reason=foamfatal ;;
        *"Floating point exception"*) reason=fpe ;;
        *"received signal"*)          reason=signal ;;
        *"CANCELLED"*)                reason=cancelled ;;
        *"Aborted"*|*"srun: error"*)  reason=aborted ;;
      esac
      # A FOREIGN LOG OVERRIDES EVERY VERDICT ABOVE. Whatever it says happened, it did
      # not happen here. Report it as its own bucket rather than silently dropping the
      # case, so the partition still sums (GEO-089).
      if [ "$foreign" = "1" ]; then
        reason="foreign-log:$(basename "$logcase")"
        it=""
      fi
      # converged-and-stopped is a POSITIVE outcome and is detected separately
      if [ -z "$reason" ]; then
        grep -q '^End' "$sf" 2>/dev/null && reason=ended
        grep -qi 'SIMPLE solution converged' "$sf" 2>/dev/null && reason=residualControl
      fi
      # the SLURM stream carries kills the solver log never sees
      so=$(ls -t "$d"/log_slurm/*.out "$d"/log_slurm/*.err "$d"/*.out "$d"/*.err 2>/dev/null | head -2)
      if [ -n "$so" ]; then
        st=$(tail -25 $so 2>/dev/null)
        case "$st" in
          *"TIME LIMIT"*) reason=timelimit ;;
          *"DUE TO NODE FAILURE"*) reason=nodefail ;;
        esac
      fi
    fi
    cells=$(grep -oP '^\s+cells:\s+\K[0-9]+' "$g"/checkMesh.log 2>/dev/null | head -1)
    fine=$(grep -aoP 'nCells:\s*\K[0-9]+' "$d"/processor0/constant/polyMesh/owner 2>/dev/null | head -1)
    bg=$(grep -aoP 'nCells:\s*\K[0-9]+' "$d"/constant/polyMesh/owner 2>/dev/null | head -1)
    times=$(ls -d "$d"/processor0/[0-9]*/ 2>/dev/null | xargs -n1 basename 2>/dev/null \
            | grep -E '^[0-9]+$' | sort -g | tr '\n' ',' )
    lt=$(ls -d "$d"/processor0/[0-9]*/ 2>/dev/null | xargs -n1 basename 2>/dev/null \
         | grep -E '^[0-9]+$' | sort -g | tail -1)
    flds=""
    [ -n "$lt" ] && flds=$(ls "$d"/processor0/"$lt" 2>/dev/null | tr '\n' ' ')
    fc=$(ls "$d"/postProcessing/force*/*/force*.dat 2>/dev/null | tail -1)
    frows=0; [ -n "$fc" ] && frows=$(grep -vc '^#' "$fc" 2>/dev/null)
    sz=$(du -sh "$d" 2>/dev/null | cut -f1)
    echo "CASE|$root|$n|$np|${it:-}|${endt:-}|${reason:-}|${cells:-}|${fine:-}|${bg:-}|${times}|${frows}|${sz}|${lt:-}|${gen}|${flds}"
  done
done
'''


def _i(v):
    """int or None. A field that will not parse is UNREAD, never silently zero."""
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--remote", default="delftblue")
    ap.add_argument("--out", default="results/campaign_3d_audit.json")
    ap.add_argument("--from-cache", action="store_true")
    a = ap.parse_args()

    cache = REPO / "results/campaign_3d_scan.txt"
    if a.from_cache and cache.exists():
        raw = cache.read_text()
    else:
        r = subprocess.run(["ssh", a.remote, "bash -s"], input=REMOTE_SCAN, text=True,
                           capture_output=True, timeout=900)
        raw = r.stdout
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(raw)
        if r.returncode != 0:
            print("ssh exited %d; stderr tail:\n%s" % (r.returncode, r.stderr[-400:]))

    rows = []
    for l in raw.splitlines():
        if not l.startswith("CASE|"):
            continue
        parts = l.split("|", 15)
        (_, root, name, np_, it, endt, reason, cells, fine, bg, times, frows, sz) = parts[:13]
        latest_t = parts[13] if len(parts) > 13 else ""
        gen = parts[14] if len(parts) > 14 else ""
        flds = (parts[15] if len(parts) > 15 else "").split()
        rows.append(dict(root=root.rsplit("/", 1)[-1], name=name, ranks=int(np_ or 0),
                         iter=int(it) if it else None,
                         endTime=int(endt) if endt else None,
                         reason=reason or None,
                         cells=_i(cells), fine_per_rank=_i(fine), bg_cells=_i(bg),
                         saved_times=[t for t in times.split(",") if t],
                         force_rows=int(frows or 0), size=sz,
                         latest_time=latest_t or None, fields_at_latest=flds,
                         log_generation=gen or None))

    # ---- CLASSIFY. Every case gets exactly one bucket.
    def bucket(r):
        if r["ranks"] == 0 and not r["fine_per_rank"]:
            return "NOT-BUILT"
        # FOREIGN LOG FIRST, and it is NOT folded into MESH-ONLY. Both mean "never solved
        # here", but one of them also means a log from another case is sitting in this
        # directory misleading anything that reads it, which is a defect to clear rather
        # than a state to record.
        if str(r["reason"]).startswith("foreign-log"):
            return "MESH-ONLY, FOREIGN LOG PRESENT (%s)" % str(r["reason"]).split(":", 1)[1]
        if r["iter"] is None:
            return "MESH-ONLY (no solve log)"
        if r["reason"] in ("fpe", "foamfatal", "signal", "heap-abort", "aborted"):
            return "CRASHED (%s)" % r["reason"]
        if r["reason"] == "nodefail":
            return "NODE FAILURE"
        if r["reason"] == "residualControl":
            return "CONVERGED on residualControl"
        if r["reason"] == "timelimit":
            return "KILLED ON TIME LIMIT"
        if r["reason"] == "cancelled":
            return "CANCELLED"
        if r["endTime"] and r["iter"] >= r["endTime"]:
            return "REACHED endTime"
        if r["reason"] == "ended":
            return "ENDED EARLY, cause unread"
        return "STOPPED SHORT, reason unread"

    def has_data(r):
        """Fields AND forces actually on disk. Requirement, 2026-08-12: a case is complete only
        when its processor data and results are on disk.
        A log line saying End is a CLAIM about a run; the field data is the run's OUTPUT.
        L11 is the worked case: both members ran 1009 healthy iterations and a later
        resubmission reset 0/, so the logs read fine and NOTHING survives to post-process."""
        need = {"U", "p"}
        return bool(need <= set(r["fields_at_latest"])) and r["force_rows"] > 0

    for r in rows:
        r["bucket"] = bucket(r)
        r["has_data"] = has_data(r)
        if r["bucket"] in ("REACHED endTime", "CONVERGED on residualControl") \
                and not r["has_data"]:
            r["bucket"] = "RAN, BUT NO DATA ON DISK"
        # rerun cost: can it be resubmitted at a DIFFERENT rank count cheaply?
        r["remesh_needed_to_change_ranks"] = bool(
            r["fine_per_rank"] and r["bg_cells"] and r["fine_per_rank"] > 4 * r["bg_cells"])

    order = sorted({r["bucket"] for r in rows},
                   key=lambda b: (not b.startswith("REACHED"),
                                  not b.startswith("CONVERGED"), b))
    print("3D CAMPAIGN AUDIT  --  %d cases found on the cluster\n" % len(rows))
    print("%-9s %-26s %5s %8s %8s %-14s %10s %6s"
          % ("root", "case", "ranks", "iter", "endTime", "stopped by", "cells", "frows"))
    for b in order:
        g = [r for r in rows if r["bucket"] == b]
        if not g:
            continue
        print("\n--- %s  (%d)" % (b, len(g)))
        for r in sorted(g, key=lambda q: (q["root"], q["name"])):
            print("%-9s %-26s %5d %8s %8s %-14s %10s %6d"
                  % (r["root"], r["name"][:26], r["ranks"],
                     r["iter"] if r["iter"] is not None else "-",
                     r["endTime"] or "-", r["reason"] or "-",
                     r["cells"] or "-", r["force_rows"]))

    tally = {}
    for r in rows:
        tally[r["bucket"]] = tally.get(r["bucket"], 0) + 1
    print("\n=== BUCKETS ===")
    for b in order:
        if b in tally:
            print("  %-32s %d" % (b, tally[b]))
    s = sum(tally.values())
    print("  %-32s %d  %s" % ("TOTAL", s,
                              "" if s == len(rows) else "*** GAP: buckets do not sum ***"))

    o = REPO / a.out
    o.write_text(json.dumps(dict(cases=rows, buckets=tally), indent=2) + "\n")
    print("\nwrote %s" % o.relative_to(REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
