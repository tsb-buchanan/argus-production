#!/usr/bin/env python3
"""overnight_recipe_study.py: run the whole recipe matrix unattended and write the study.

BRIEF, 2026-08-11: figures for every case, a report naming the chosen recipe, a comparison of
all of them on the specific geometries, cl/cd polars, noise levels, convergence statistics,
and THE THEORY REASONS why the winning settings work. A real study, not a sweep.

WHY A QUEUE RUNNER RATHER THAN A SCRIPT PER STAGE. The matrix is larger than one sitting and
the machine has 24 cores, so the binding constraint is scheduling, not compute. A runner that
holds the slot count, launches the next job as one frees, collects each case's figures the
moment it finishes and writes a partial report after EVERY job means that whatever has
completed by morning is already written up. A script that only reports at the end reports
nothing if it is interrupted, which for an unattended overnight run is the likely case.

THE ORDER IS BY DECISIVENESS, NOT BY CONVENIENCE:
  1. ARGUS baseline NORMAL, fixed alpha. Re 1.68e7, 2.6x the RAE case, and where every
     previous attempt died. It separates the recipes fastest.
  2. ARGUS baseline FREESTREAM, fixed alpha. M 0.780, past drag divergence, Re 2.11e7. The
     hardest case in the matrix.
  3. cl/cd POLARS, but only for recipes that survived 1 and 2. A polar of a recipe that cannot
     hold a single point is five times the cost for no information.
That ordering answers "which recipe is versatile" first and "how accurate is it" second, which
is the order the decision needs them in.

NOTHING HERE ASSERTS A RESULT IT DID NOT READ. Every number in the report is re-derived from
the case directory by scripts/recipe_trial.py and scripts/report_figures_all.py; this file
schedules and formats, and it records a FAILED job as a named bucket rather than dropping it
(GEO-089).
"""
import json
import pathlib
import re
import shutil
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
FOAM = ("export PATH=/usr/bin:$PATH; source /opt/openfoam7/etc/bashrc >/dev/null 2>&1; "
        "export PATH=/usr/bin:$PATH; ")
TRIAL = REPO / "cases/recipe_trial"
RANKS = 6
SLOTS = 4                      # 24 cores / 6 ranks
POLL = 60

SOLVER = {"hybrid_krylov": "rhoSimpleFoam", "wolf": "rhoSimpleFoam", "tutorial": "rhoSimpleFoam", "lts": "rhoPimpleFoam",
          "lts_ll": "rhoPimpleFoam", "alletto": "rhoSimpleFoam", "hybrid": "rhoSimpleFoam",
          "alletto_tvd": "rhoSimpleFoam", "alletto_lim": "rhoSimpleFoam"}

# The four carried forward to ARGUS, and WHY each is here rather than the others.
#   wolf        the incumbent and the most accurate arm on RAE. Included precisely because it
#               is expected to fail at this Reynolds number: a control that cannot fail is
#               not a control, and its failure mode is the thing being fixed.
#   hybrid      wolf's TVD, mesh-tolerant discretisation driven by Alletto's
#               SIMPLEC-appropriate iteration. The leading hypothesis.
#   alletto_lim Alletto's convection with wolf's limited laplacian and one non-orthogonal
#               corrector. Tests whether the ONERA M6 convection is usable once the
#               `corrected`/nNonOrth-0 pairing that is wrong for our meshes is removed.
#   lts_ll      a genuinely different time-integration path. If the steady arms all fail at
#               freestream, local time stepping is the fallback, and it must be measured
#               rather than assumed available.
# EXCLUDED WITH REASONS, so the set is partitioned rather than curated:
#   tutorial    C_d,pressure negative for 8,500 iterations on RAE. Out on a physical gate.
#   lts         converged with NO supersonic pocket at all; Cp RMS 0.49. Out on accuracy.
#   alletto     T_max 400 K against T0 318.9, above T0 on 51.8% of iterations. Out on a gate.
#   alletto_tvd FPE in GAMG's coarsest level at iteration 110. Out, crashed.
#   hybrid      SUPERSEDED by hybrid_krylov. It died with GAMG diverging on the MOMENTUM
#               equation (Ux final residual 7.0e+101 at iteration cap) while every physical
#               field was healthy. Multigrid assumes an elliptic-dominated operator; the
#               momentum equation here is convection-dominated on cells of aspect ratio 1e4-
#               1e5, so the coarse-grid correction amplifies. Kept in the study as the arm
#               that identified the mechanism.
ARGUS_ARMS = ["hybrid_krylov", "alletto_lim", "lts_ll", "wolf"]

JOBS = [dict(tag="rae2822_case09", arm="hybrid_krylov", iters=60000,
             src="cases/validation/rae2822/wall_modelled/case09")]
for arm in ARGUS_ARMS:
    JOBS.append(dict(tag="argus_baseline_normal", arm=arm, iters=60000,
                     src="cases/argus2d/wall_modelled/baseline_normal"))
for arm in ARGUS_ARMS:
    JOBS.append(dict(tag="argus_baseline_freestream", arm=arm, iters=60000,
                     src="cases/argus2d/wall_modelled/baseline_freestream"))


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)
    with open(REPO / "results/overnight_study.log", "a") as f:
        f.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), m))


def running_dirs():
    """Case directories with a live solver, from /proc. Never from the process NAME.

    Matching on the name would confuse two arms that share a solver, which is most of them;
    it nearly cost the alletto arm when the tutorial one was stopped.
    """
    out = set()
    for name in ("rhoSimpleFoam", "rhoPimpleFoam"):
        r = subprocess.run(["pgrep", "-x", name], capture_output=True, text=True)
        for p in r.stdout.split():
            try:
                out.add(pathlib.Path("/proc/%s/cwd" % p).resolve())
            except OSError:
                pass
    return out


def launch(job):
    d = TRIAL / ("%s_%s" % (job["tag"], job["arm"]))
    r = subprocess.run([sys.executable, "scripts/recipe_trial.py", "--arms", job["arm"],
                        "--build-only", "--ranks", str(RANKS), "--iters", str(job["iters"]),
                        "--src", job["src"], "--tag", job["tag"]],
                       cwd=str(REPO), capture_output=True, text=True, timeout=3600)
    if r.returncode != 0 or not (d / "processor0").exists():
        log("  BUILD FAILED %s/%s\n%s" % (job["tag"], job["arm"], r.stdout[-700:]))
        job["failed"] = "build/audit failed: %s" % r.stdout.strip().splitlines()[-1:]
        return None
    # ---- EACH RUN IN ITS OWN PROCESS GROUP (setsid).
    # Two runs launched from one shell share a process group, and when mpirun aborts a job it
    # can signal the group. Tonight arm E and arm F were co-launched and both stopped within
    # seconds of each other; F had genuinely crashed. F's crash turned out not to be what
    # killed E -- E diverged on its own linear solver -- but I could not establish that from
    # the logs alone and had to reconstruct it from PIDs. AN EXPERIMENT WHOSE ARMS CAN KILL
    # EACH OTHER CANNOT ATTRIBUTE A FAILURE TO ITS ARM, which is the only thing it is for.
    cmd = "mpirun -np %d %s -parallel > log.solver 2>&1" % (RANKS, SOLVER[job["arm"]])
    subprocess.Popen(["setsid", "bash", "-lc", FOAM + cmd], cwd=str(d),
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                     start_new_session=True)
    log("  launched %s/%s" % (job["tag"], job["arm"]))
    return d


def collect(job):
    d = TRIAL / ("%s_%s" % (job["tag"], job["arm"]))
    subprocess.run([sys.executable, "scripts/report_figures_all.py", str(d),
                    "--out-dir", "results/figures/recipe_trial"],
                   cwd=str(REPO), capture_output=True, text=True, timeout=5400)
    log("  collected %s/%s" % (job["tag"], job["arm"]))


# ---- THE POLAR, on the project protocol (2026-08-11): cold at the SEED alpha, warm OUTWARD.
# The seed run is the fixed-alpha stability case already in the queue, so the polar costs four
# extra points rather than five, and the point carrying the most weight in the interpolation
# is the one with no inheritance.
#
# WHY A POLAR AND NOT THE TRIM LOOP. The trim tolerance |dCL| < 1e-4 is 60 to 236 times smaller
# than the measured cl oscillation amplitude, so a secant -- which DIVIDES by the cl difference
# between two iterates -- amplifies the noise. A least-squares fit over a fixed alpha list
# AVERAGES it: five points over 0.75 deg give an alpha uncertainty of 0.007 deg against a
# trimmed-alpha difference between geometries of 1.796 deg. It also fails safe, because one bad
# alpha costs one point instead of the whole chain.
#
# THE SAME OFFSETS ON EVERY GEOMETRY. Any warm-start bias is then common to both sides of the
# morphing delta and cancels (D071); a protocol that differs between geometries puts the bias
# straight into the number we are actually after.
POLAR_OFFSETS = [-0.40, -0.20, +0.20, +0.40]      # deg, about the seed; seed itself is run 0
POLAR_ITERS = 20000                                # warm continuation, not a cold start


def polar_jobs(survivors):
    """Four extra alphas per surviving arm, warm-continued from its converged seed run."""
    out = []
    for tag, arm in survivors:
        for off in POLAR_OFFSETS:
            out.append(dict(tag=tag, arm=arm, polar_offset=off, iters=POLAR_ITERS))
    return out


def run_polar_point(job, seed_alpha, U):
    """Copy the converged seed case, rotate the freestream, re-verify, run."""
    sys.path.insert(0, str(REPO / "scripts"))
    import trim2d_transonic as T
    import solver_recipe
    src = TRIAL / ("%s_%s" % (job["tag"], job["arm"]))
    a = seed_alpha + job["polar_offset"]
    d = TRIAL / ("%s_%s_a%s" % (job["tag"], job["arm"],
                                ("m%.2f" % abs(a) if a < 0 else "p%.2f" % a).replace(".", "p")))
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(src, d, ignore=shutil.ignore_patterns("log.*", "postProcessing", "figures"))
    # ---- U COMES FROM THE CASE, NOT FROM condition(). The polar stage refused its first
    # point because the audit found magUInf 233.8436 in the case against |0/U| 233.5834 after
    # the rotation. Both are "the freestream speed" and they differ by 0.111%, because
    # condition() computes the speed of sound with R = 287.058 while the case is built from
    # the thermophysicalProperties molWeight of 28.9, i.e. R = 287.698. sqrt of that ratio is
    # exactly the discrepancy. Rotating with the recomputed value would have left the case
    # internally inconsistent, and the gate caught it (GEO-080: derive from the thing itself).
    cd_ = (d / "system/controlDict").read_text()
    U = float(re.search(r"^\s*magUInf\s+([-0-9.eE+]+);", cd_, re.M).group(1))
    T.rotate_freestream(d, a, U)
    # THE CASE THAT WAS AUDITED AT BUILD IS NOT THE CASE BEING SOLVED once the freestream has
    # been rotated: 0/U, liftDir and dragDir all changed. Re-verify and re-audit (ARG-096).
    bad = solver_recipe.verify(d, job["arm"], verbose=False)
    if bad:
        log("  polar %s: recipe no longer verifies: %s" % (d.name, bad[:1]))
        return None
    c = d / "system/controlDict"
    s = c.read_text()
    s = re.sub(r"^startFrom\s+\S+;", "startFrom       latestTime;", s, flags=re.M)
    last = max((int(q.name) for q in (d / "processor0").iterdir()
                if q.is_dir() and q.name.isdigit()), default=0)
    s = re.sub(r"^endTime\s+\S+;", "endTime         %d;" % (last + job["iters"]), s, flags=re.M)
    c.write_text(s)
    g = subprocess.run([sys.executable, "scripts/audit_case_settings.py", "--quiet-pass", str(d)],
                       cwd=str(REPO), capture_output=True, text=True, timeout=1800)
    if g.returncode != 0:
        log("  polar %s: AUDIT FAILED\n%s" % (d.name, g.stdout[-500:]))
        return None
    cmd = "mpirun -np %d %s -parallel > log.solver 2>&1" % (RANKS, SOLVER[job["arm"]])
    subprocess.Popen(["setsid", "bash", "-lc", FOAM + cmd], cwd=str(d),
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                     start_new_session=True)
    log("  polar launched %s at alpha %+.4f (warm from %d)" % (d.name, a, last))
    return d


def main():
    (REPO / "results").mkdir(exist_ok=True)
    log("=== overnight recipe study: %d jobs, %d slots of %d ranks" % (len(JOBS), SLOTS, RANKS))
    # ---- SKIP WHAT IS ALREADY DONE OR ALREADY RUNNING, so the runner is RESTARTABLE.
    # The first launch lost two jobs to a build failure and there was no way to add them back
    # without restarting, and restarting would have relaunched the two healthy jobs from
    # iteration zero and destroyed hours of progress. A runner that cannot be restarted has to
    # be got right first time, which for an unattended overnight is not a property to rely on.
    live_now = running_dirs()
    queue, active = [], {}
    for j in JOBS:
        d = TRIAL / ("%s_%s" % (j["tag"], j["arm"]))
        if d.resolve() in live_now:
            active["%s/%s" % (j["tag"], j["arm"])] = (j, d)
            log("  ADOPTED already-running %s/%s" % (j["tag"], j["arm"]))
            continue
        lg = d / "log.solver"
        if lg.exists():
            it = re.findall(r"^Time = (\d+)", lg.read_text(errors="replace"), re.M)
            if it and int(max(it, key=int)) >= j["iters"]:
                log("  SKIP complete %s/%s at %s iterations" % (j["tag"], j["arm"],
                                                                max(it, key=int)))
                continue
        queue.append(j)
    log("  %d queued, %d adopted" % (len(queue), len(active)))
    while queue or active:
        live = running_dirs()
        for key in [k for k, (j, d) in active.items() if d.resolve() not in live]:
            j, d = active.pop(key)
            log("  finished %s (%s)" % (key, (d / "log.solver").exists() and "log present"))
            try:
                collect(j)
            except Exception as e:
                log("  collect failed for %s: %s" % (key, e))
            try:
                write_report()
            except Exception as e:
                log("  report failed: %s" % e)
        while queue and len(running_dirs()) < SLOTS:
            j = queue.pop(0)
            d = launch(j)
            if d is not None:
                active["%s/%s" % (j["tag"], j["arm"])] = (j, d)
                time.sleep(20)          # let the solver take its slot before recounting
            else:
                try:
                    write_report()
                except Exception:
                    pass
        time.sleep(POLL)
    log("=== fixed-alpha matrix done; selecting survivors for the polar stage")
    # ---- SURVIVORS ONLY. A polar of a recipe that cannot hold a single point is five times
    # the cost for no information, so the fixed-alpha matrix gates the polar rather than the
    # two running in parallel.
    import json as _json
    sys.path.insert(0, str(REPO / "scripts"))
    import recipe_trial as R
    import wallresolved_transonic_routes as W
    import trim2d_transonic as T
    survivors = []
    for tag in ("argus_baseline_normal", "argus_baseline_freestream"):
        for arm in ARGUS_ARMS:
            d = TRIAL / ("%s_%s" % (tag, arm))
            if not (d / "log.solver").exists():
                continue
            m = R.crash_markers(d)
            g = W.temperature_ceiling(d) or {}
            c = W.pressure_drag_sign(d) or {}
            ok = (not m.get("foam_fatal") and not m.get("h_residual_at_one")
                  and not m.get("linear_solver_diverged") and not g.get("exceeded_T0")
                  and not c.get("cd_press_went_negative"))
            log("  %s/%s -> %s" % (tag, arm, "SURVIVOR" if ok else "excluded from the polar"))
            if ok:
                survivors.append((tag, arm))
    seed = {"argus_baseline_normal": (T.condition("early", "normal"),
                                      T.KNOWN_TRIM_ALPHA[("baseline", "normal")]),
            "argus_baseline_freestream": (T.condition("early", "freestream"),
                                          T.KNOWN_TRIM_ALPHA[("baseline", "freestream")])}
    pq = polar_jobs(survivors)
    log("  %d polar points queued for %d survivors" % (len(pq), len(survivors)))
    pactive = {}
    while pq or pactive:
        for k in [k for k, d in pactive.items() if d.resolve() not in running_dirs()]:
            pactive.pop(k)
            log("  polar finished %s" % k)
            try:
                subprocess.run([sys.executable, "scripts/report_figures_all.py",
                                str(TRIAL / k), "--out-dir", "results/figures/recipe_trial"],
                               cwd=str(REPO), capture_output=True, text=True, timeout=5400)
                write_report()
            except Exception as e:
                log("  polar collect failed %s: %s" % (k, e))
        while pq and len(running_dirs()) < SLOTS:
            j = pq.pop(0)
            cond, a0 = seed[j["tag"]]
            try:
                d = run_polar_point(j, a0, cond["V"])
            except Exception as e:
                log("  polar launch failed %s/%s%+.2f: %s" % (j["tag"], j["arm"],
                                                              j["polar_offset"], e))
                d = None
            if d is not None:
                pactive[d.name] = d
                time.sleep(20)
        time.sleep(POLL)
    log("=== all jobs done")
    write_report()
    return 0


def write_report():
    """Regenerate the study from whatever is on disk, after EVERY job."""
    r = subprocess.run([sys.executable, "scripts/write_recipe_study.py"],
                       cwd=str(REPO), capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        log("  write_recipe_study failed:\n%s" % (r.stderr or r.stdout)[-900:])


if __name__ == "__main__":
    raise SystemExit(main())
