#!/usr/bin/env python3
"""splice_postpro.py: turn the per-restart force fragments into one history per case,
and write an index that says what is plottable and what is not.

WHY SPLICING IS NOT CONCATENATION. A case that restarts writes
postProcessing/forceCoeffs1/0/, then /600/, then /4500/. Each segment starts AT its
restart time, so the segments OVERLAP at the seam: the last row of one and the first
row of the next describe the same iteration. Concatenating gives duplicate
iterations, and any drift computed across a duplicate is nonsense. Splicing keeps the
LATER value at a repeated iteration, because that is the one the continuing run
actually used.

WHY THE INDEX MATTERS MORE THAN THE SPLICE. Every consumer of this tree asks the same
two questions -- how far did it get, and is it converged -- and the honest answer for
most cases is "not converged", which a bare .dat file does not say. The index answers
both per case, against the pinned gate, so a plot script cannot silently treat a
mid-transient case as a result.
"""
import json
import sys
from pathlib import Path

# hpc/v5_recipe.json convergence_gate, pinned.
WINDOW_SAMPLES = 200
CD_DRIFT_MAX_CT = 0.05      # raised from 0.01 on 2026-09-01, see drift() below
CL_DRIFT_MAX_CT = 1.0
TARGET_CL_CR = 0.428277635108
TARGET_CL_COMPRESSIBLE = 0.52929708745
ENDTIME = 4500


def read_segments(case_dir, fo="forceCoeffs1", stem="forceCoeffs.dat"):
    """Every (iteration -> row) pair across all restart segments, later wins."""
    pts, cols = {}, None
    for f in sorted((case_dir / "postProcessing" / fo).glob(f"*/{stem}"),
                    key=lambda p: float(p.parent.name)):
        try:
            txt = f.read_text(errors="replace").splitlines()
        except OSError:
            continue
        hdr = [l for l in txt if l.startswith("#")]
        if hdr:
            c = hdr[-1].lstrip("#").split()
            if "Cl" in c:
                cols = c
        if not cols:
            continue
        ic, idd = cols.index("Cl"), cols.index("Cd")
        icm = cols.index("Cm") if "Cm" in cols else None
        for r in txt:
            if not r or r.startswith("#"):
                continue
            p = r.split()
            if len(p) <= max(ic, idd):
                continue
            try:
                it = int(float(p[0]))
            except ValueError:
                continue
            pts[it] = (float(p[idd]), float(p[ic]),
                       float(p[icm]) if icm is not None and len(p) > icm else None)
    return dict(sorted(pts.items()))


def drift(hist):
    """DRIFT OF THE MEAN over the trailing window, in counts. Half-to-half.

    THIS IS THE THIRD PLACE THE SAME QUANTITY IS COMPUTED and it was the third different
    formula. It took the ENDPOINT DIFFERENCE, last sample minus first sample of the
    window. Both endpoints are single draws from a limit cycle that never decays, so that
    statistic carries TWICE the oscillation amplitude as noise and says almost nothing
    about whether the mean has moved. Measured 2026-09-02 on SWB_a1p60: endpoint
    difference +0.061 ct against a true mean drift of 0.0019 ct, a factor of 32.

    The consequence was not academic. The solve-time gate had finalised three cases as
    CONVERGED while this one wrote "not converged" into INDEX.json, and INDEX.json is
    what decides whether a case is QUOTABLE and whether its surfaces get pulled. Two
    gates on one quantity disagreeing means at least one is wrong and nothing says which.

    Now identical to the solve-time gate: mean of the second half minus mean of the
    first, against 0.05 counts on Cd and 1.0 on Cl.
    """
    its = sorted(hist)
    if len(its) < 4:
        return None, None, None
    win = its[-WINDOW_SAMPLES:] if len(its) >= WINDOW_SAMPLES else its
    h = len(win) // 2
    def mean(idx, k):
        return sum(hist[i][k] for i in idx) / len(idx)
    dcd = (mean(win[h:], 0) - mean(win[:h], 0)) * 1e4
    dcl = (mean(win[h:], 1) - mean(win[:h], 1)) * 1e4
    return dcd, dcl, (win[0], win[-1], len(win))


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "results/postpro_latest")
    cases = sorted(p for p in root.iterdir() if p.is_dir())
    index, n_conv = {}, 0

    print(f"  {'case':<14}{'samples':>8}{'iter':>7}{'Cd (ct)':>10}{'Cl':>10}"
          f"{'dCd':>9}{'dCl':>9}  gate")
    for c in cases:
        hist = read_segments(c)
        if not hist:
            index[c.name] = {"status": "no force history"}
            print(f"  {c.name:<14}{'-':>8}{'-':>7}   no force history")
            continue
        its = sorted(hist)
        last = its[-1]
        dcd, dcl, win = drift(hist)

        # SPLICED HISTORY AS ONE FILE, so a plot script never has to know about restarts.
        out = c / "forceCoeffs_spliced.dat"
        with out.open("w") as fh:
            fh.write("# spliced across restart segments; later value wins at a repeat\n")
            fh.write("# Time\tCd\tCl\tCm\n")
            for it in its:
                cd, cl, cm = hist[it]
                fh.write(f"{it}\t{cd:.10e}\t{cl:.10e}\t"
                         f"{'' if cm is None else format(cm, '.10e')}\n")

        enough = win is not None and win[2] >= WINDOW_SAMPLES
        passed = (enough and dcd is not None
                  and abs(dcd) <= CD_DRIFT_MAX_CT and abs(dcl) <= CL_DRIFT_MAX_CT)
        if passed:
            n_conv += 1
        gate = ("PASS" if passed
                else ("short window" if not enough else "not converged"))
        tgt = TARGET_CL_COMPRESSIBLE if c.name.startswith("CMPB") else TARGET_CL_CR
        surfaces = sorted(p.name for p in (c / "surfaces").glob("*")) if (c / "surfaces").is_dir() else []

        index[c.name] = {
            "samples": len(its), "last_iteration": last, "reached_endtime": last >= ENDTIME,
            "Cd": hist[last][0], "Cl": hist[last][1], "Cm": hist[last][2],
            "Cl_minus_target": hist[last][1] - tgt, "target_Cl": tgt,
            "drift_window": {"from": win[0], "to": win[1], "samples": win[2]} if win else None,
            "dCd_counts": dcd, "dCl_counts": dcl,
            "gate": {"cd_max_ct": CD_DRIFT_MAX_CT, "cl_max_ct": CL_DRIFT_MAX_CT,
                     "window_samples": WINDOW_SAMPLES, "verdict": gate},
            "surface_times": surfaces,
            "spliced_file": str(out.relative_to(root)),
            "quotable": bool(passed and last >= ENDTIME),
        }
        # A DRIFT OVER A WINDOW SHORTER THAN THE GATE'S IS NOT A DRIFT. Those windows
        # span the start-up transient, so the number is enormous and meaningless
        # (-6979 counts was printed for a case at iteration 305). Printing it invites
        # someone to quote it. Show a dash and let the verdict carry the information.
        ds = f"{dcd:>+9.3f}{dcl:>+9.2f}" if enough else f"{'-':>9}{'-':>9}"
        print(f"  {c.name:<14}{len(its):>8}{last:>7}{hist[last][0]*1e4:>10.2f}"
              f"{hist[last][1]:>10.5f}{ds}  {gate}")

    (root / "INDEX.json").write_text(json.dumps(index, indent=2))
    print(f"\n  index written: {root/'INDEX.json'}")
    print(f"  {n_conv} of {len(cases)} case(s) pass the pinned drift gate")
    # A CASE IS ONLY 'quotable' IF IT BOTH REACHED endTime AND PASSED THE GATE.
    q = [k for k, v in index.items() if v.get("quotable")]
    print(f"  quotable (reached endTime AND converged): {', '.join(q) if q else 'none'}")


if __name__ == "__main__":
    main()
