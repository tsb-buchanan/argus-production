#!/usr/bin/env python3
"""watch_run.py: live monitor for a running OpenFOAM case -- residuals, forces, y+.

WHY (project decision, 2026-08-08). Watching a transonic case converge means answering three questions at
once: are the residuals falling, has the force settled, and is the wall spacing where it was
designed to be. Those live in three different files, and grepping them one at a time while a
run is in progress is how I have repeatedly read a force off an unconverged case today.

WHAT IT READS, and it prefers the function-object output over the log:
  postProcessing/residuals/*/residuals*.dat   per-equation initial residuals
  postProcessing/forces_coeffs/*/*.dat        Cd, Cl, Cm
  postProcessing/yplus/*/yPlus.dat            min/max/mean y+ on each patch
  log.*                                       fallback, and the only source of `bounding`

RESTARTS ARE CONCATENATED, NOT OVERWRITTEN. OpenFOAM writes residuals.dat, then residuals_1.dat
on the next start, and so on. Reading only the newest file silently drops everything before the
last restart; reading them in glob order gets them in the WRONG order once there are ten
(residuals_10 sorts before residuals_2). Files are ordered by their numeric suffix and
concatenated, with duplicate times from an overlapping restart dropped.

THE CONVERGENCE TEST IS D023, not the residual. A residual measures how well the current
linear system was solved; it says nothing about whether the FORCE has stopped moving. This
project has twice quoted a drag that was still drifting under a converged-looking residual, so
D023 -- |mean(2nd half) - mean(1st half)| over the last N samples, in drag counts -- is printed
beside it and is the number that decides quotability.

BOUNDING IS REPORTED AS A DISTRIBUTION. A total count cannot tell a startup transient from
ongoing clipping, and those are opposite verdicts. The final-fifth count is what matters, and
the relative magnitude matters too: a k undershoot of 1e-4 against a mean of 7 is noise, while
the same count on the pressure field is a diverging solution held up by a limiter.
"""
import argparse
import glob
import os
import pathlib
import re
import statistics
import sys
import time

import numpy as np


def _ordered(pattern):
    """Function-object files in RESTART order.

    OpenFOAM writes a function object's output to postProcessing/<name>/<START TIME>/<file>.
    A warm restart at t=15000 therefore produces .../0/forceCoeffs.dat AND
    .../15000/forceCoeffs.dat -- two IDENTICALLY NAMED files in DIFFERENT time directories.
    Sorting on a `_N` suffix in the FILENAME gives both the same key, so glob order decided
    which was read last, and the trim read the COLD run's final cl after a warm continue had
    already moved it. Two different alphas returned identical cl to five decimals, the secant
    slope came out zero, and the driver died on a division by zero.

    Order on the TIME DIRECTORY first, then the filename suffix. Both are numeric and both
    matter: _1/_2 handle a restart into the same directory, the directory handles a restart
    into a new one.
    """
    def key(p):
        q = pathlib.Path(p)
        try:
            tdir = float(q.parent.name)
        except ValueError:
            tdir = -1.0
        m = re.search(r"_(\d+)\.dat$", q.name)
        return (tdir, int(m.group(1)) if m else 0)
    return sorted(glob.glob(pattern), key=key)


def read_dat(case, sub, name):
    """Concatenate a function object's files across restarts, dropping duplicate times."""
    files = _ordered(str(pathlib.Path(case) / "postProcessing" / sub / "*" / name))
    if not files:
        return None, []
    hdr, rows, seen = [], [], {}
    for f in files:
        for line in open(f, errors="replace"):
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if "Time" in s or s.count(" ") > 2:
                    hdr = s.lstrip("# ").split()
                continue
            v = s.split()
            try:
                t = float(v[0])
            except ValueError:
                continue
            # KEEP NON-NUMERIC COLUMNS. yPlus.dat carries the PATCH NAME as column 1
            # ("Time patch min max average"), so a blanket float() raises on a
            # perfectly well-formed file. Non-numeric entries pass through as strings.
            row = []
            for x in v:
                try:
                    row.append(float(x))
                except ValueError:
                    row.append(x)
            # LATER FILES WIN ON A REPEATED TIME. A rerun restarts at t = 0 and writes
            # forceCoeffs_1.dat alongside the old forceCoeffs.dat, so keeping the FIRST
            # occurrence of each time keeps the DEAD RUN and silently discards the live
            # one. This tool reported Cd -10831 from a crashed GAMG run while its healthy
            # replacement sat in the next file. Files are read in restart order, so the
            # last writer of a given time is the current one.
            if t in seen:
                rows[seen[t]] = row
            else:
                seen[t] = len(rows)
                rows.append(row)
    return hdr, rows


def d023(cd, n=2000):
    w = cd[-n:] if len(cd) >= n else cd
    h = len(w) // 2
    return abs(statistics.mean(w[h:]) - statistics.mean(w[:h])) * 1e4 if h else float("nan")


def bounding(case):
    """Distribution, not a total: where the events fall decides what they mean."""
    logs = _ordered(str(pathlib.Path(case) / "log.*"))
    log = None
    for f in sorted(glob.glob(str(pathlib.Path(case) / "log.*")), key=os.path.getmtime):
        if re.search(r"log\.(solver|rho\w*Foam|simpleFoam)$", f):
            log = f
    if not log:
        return {}, 0, 0
    t = open(log, errors="replace").read()
    tot = t.count("\nTime = ")
    it, tail, fields = 0, 0, {}
    for l in t.splitlines():
        if l.startswith("Time = "):
            it += 1
        elif "bounding" in l:
            f = l.split("bounding")[1].split(",")[0].strip()
            fields[f] = fields.get(f, 0) + 1
            if tot and it > 0.8 * tot:
                tail += 1
    return fields, tail, tot


def snapshot(case, ref_cd=None, ref_cl=None):
    out = {"case": pathlib.Path(case).name}
    hdr, rows = read_dat(case, "residuals", "residuals*.dat")
    if rows:
        rows = [[x for x in r if isinstance(x, float)] for r in rows]
        w = max(len(r) for r in rows)
        a = np.array([r for r in rows if len(r) == w])
        out["iters"] = int(a[-1, 0])
        out["res"] = {hdr[i]: a[-1, i] for i in range(1, min(len(hdr), a.shape[1]))} if hdr else {}
        out["res_hist"] = a
        out["res_names"] = hdr or []
    hdr, rows = read_dat(case, "forces_coeffs", "forceCoeffs*.dat")
    if rows:
        rows = [[x for x in r if isinstance(x, float)] for r in rows]
        w = max(len(r) for r in rows)
        a = np.array([r for r in rows if len(r) == w])
        cols = {n: i for i, n in enumerate(hdr)} if hdr else {}
        ci = cols.get("Cd", 2); li = cols.get("Cl", 3)
        cd, cl = a[:, ci].tolist(), a[:, li].tolist()
        out["cd"], out["cl"] = cd[-1], cl[-1]
        out["d023"] = d023(cd)
        out["cd_hist"], out["t_hist"] = cd, a[:, 0].tolist()
        if ref_cd:
            out["cd_err"] = 100 * (cd[-1] / ref_cd - 1)
        if ref_cl:
            out["cl_err"] = 100 * (cl[-1] / ref_cl - 1)
    hdr, rows = read_dat(case, "yplus", "yPlus*.dat")
    if rows:
        num = [x for x in rows[-1] if isinstance(x, float)]
        out["yplus"] = num[-3:] if len(num) >= 3 else None
        out["yplus_patch"] = next((x for x in rows[-1] if isinstance(x, str)), "")
    out["bounding"], out["bnd_tail"], out["n_log"] = bounding(case)
    return out


def render(snaps, plot=None):
    print("\033[2J\033[H", end="")
    print("  %-18s %7s  %-10s %-10s %-8s  %-9s  %s"
          % ("case", "iters", "Cd", "Cl", "D023", "y+ avg", "bounding(tail)"))
    print("  " + "-" * 92)
    for s in snaps:
        b = "%d(%d)" % (sum(s.get("bounding", {}).values()), s.get("bnd_tail", 0))
        cd = ("%.6f" % s["cd"]) if "cd" in s else "-"
        if "cd_err" in s:
            cd += " %+.1f%%" % s["cd_err"]
        yp = ("%.3f" % s["yplus"][2]) if s.get("yplus") else "-"
        gate = "QUOTABLE" if (s.get("d023", 9e9) < 0.2 and s.get("bnd_tail", 1) == 0) else ""
        print("  %-18s %7s  %-18s %-10s %-8s  %-9s  %-10s %s"
              % (s["case"], s.get("iters", s.get("n_log", "-")), cd,
                 ("%.5f" % s["cl"]) if "cl" in s else "-",
                 ("%.3f" % s["d023"]) if "d023" in s else "-", yp, b, gate))
        if s.get("res"):
            print("      residuals: " + "  ".join("%s %.2e" % (k, v) for k, v in s["res"].items()))
    print("\n  D023 < 0.2 counts AND zero bounding in the final fifth = quotable.")
    print("  Ctrl-C to stop.  --plot writes a PNG each refresh.")
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        for s in snaps:
            if s.get("res_hist") is not None:
                a, names = s["res_hist"], s.get("res_names") or []
                for i in range(1, a.shape[1]):
                    # LABEL BY FIELD NAME. Labelling by column index made the plot
                    # unreadable: five unlabelled curves and no way to tell which equation
                    # was stalling, which is the whole question a residual plot answers.
                    nm = names[i] if i < len(names) else "col%d" % i
                    ax[0].semilogy(a[:, 0], np.abs(a[:, i]), lw=1.2,
                                   label="%s%s" % (nm, "" if len(snaps) == 1
                                                   else "  " + s["case"]))
            if s.get("cd_hist"):
                ax[1].plot(s["t_hist"], s["cd_hist"], lw=1.3, label=s["case"])
        ax[0].set_ylabel("initial residual"); ax[0].grid(alpha=.3, which="both")
        ax[0].legend(frameon=False, fontsize=8, ncol=3)
        ax[1].set_ylabel("$C_d$"); ax[1].set_xlabel("iteration"); ax[1].grid(alpha=.3)
        ax[1].legend(frameon=False, fontsize=8)
        for a_ in ax:
            for sp in ("top", "right"):
                a_.spines[sp].set_visible(False)
        fig.tight_layout(); fig.savefig(plot, dpi=120); plt.close(fig)
        print("  plot -> %s" % plot)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--ref-cd", type=float, default=None)
    ap.add_argument("--ref-cl", type=float, default=None)
    ap.add_argument("--interval", type=float, default=20.0)
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--plot", default=None, help="write a PNG of residuals and Cd each refresh")
    a = ap.parse_args()
    try:
        while True:
            render([snapshot(c, a.ref_cd, a.ref_cl) for c in a.cases], a.plot)
            if a.once:
                return 0
            time.sleep(a.interval)
    except KeyboardInterrupt:
        print("\n  stopped.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
