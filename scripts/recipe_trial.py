#!/usr/bin/env python3
"""recipe_trial.py: PHASE 0. Decide the compressible recipe by measurement on RAE 2822 case 9.

THE QUESTION. Every ARGUS 2D case dies at the corrected Reynolds number with the energy
residual pegged at 1.0 and T clamped against both limitTemperature stops. The suspect is
`transonic yes`, which makes rho.relax() unreachable in OpenFOAM-7's rhoSimpleFoam (see
scripts/solver_recipe.py for the source reading). The suspect is NOT settled by argument, so
three recipes run on the one case in this project that has experimental data.

THREE ARMS, IDENTICAL MESH, IDENTICAL GEOMETRY, IDENTICAL COLD START from 0_org:
    A wolf      rhoSimpleFoam, SIMPLEC, transonic yes,  no field relaxation
    B tutorial  rhoSimpleFoam, SIMPLE,  transonic off,  fields{p 0.7, rho 0.01}
    C lts       rhoPimpleFoam, localEuler, no relaxation factors at all

WHAT IS HELD COMMON, and why it matters (D071: does the defect differ between the two things
being differenced?). Mesh, geometry, thermophysical model, wall treatment, turbulence model,
the momentum and energy divSchemes, the temperature limiter and THE PRESSURE CLAMPS. The
clamps in particular are corrected in all three rather than left at the historical
pMinFactor 0.5, which sits 2.3% below the AGARD measured suction peak; leaving it in arm A
alone would have credited B and C with fixing a defect orthogonal to the question.

THE DECISION RULE IS PRE-REGISTERED (the approved plan, 2026-08-10) and is repeated in
verdict() so it cannot drift to fit the numbers:
    1. any arm that trips a physical gate is OUT regardless of its forces
    2. among survivors: shock position within one tap spacing (0.025 c) of 0.550, then Cp RMS,
       then C_d
    3. robustness breaks ties, and A must win on accuracy by a clear margin to be retained,
       being already known to be marginal at Re 6.5e6 and to fail outright at 1.7e7
    4. if C matches B on accuracy, C is preferred for the 3D path

ALL THREE ARE REPORTED IN FULL, INCLUDING ANY THAT FAILS (GEO-089). The arms are asserted to
partition into converged / ran-to-endTime / crashed / not-launched, with a reason on every
member of the last two buckets. A comparison that quotes only the winner is not a comparison.
"""
import argparse
import json
import pathlib
import re
import shutil
import subprocess
import sys
import time

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
FOAM = ("export PATH=/usr/bin:$PATH; source /opt/openfoam7/etc/bashrc >/dev/null 2>&1; "
        "export PATH=/usr/bin:$PATH; ")
# ---- THE TRIAL RUNS ON MORE THAN ONE GEOMETRY, and the SOURCE is what selects which.
# RAE 2822 case 9 is the only case with EXPERIMENT, so it screens for accuracy. The ARGUS
# baseline at Re 1.68e7 is the case every recipe has historically DIED on, so it screens for
# stability. A recipe is only in the toolbox if it passes both, which is why one driver runs
# both rather than two drivers drifting apart.
SRC = REPO / "cases/validation/rae2822/wall_modelled/case09"
TRIAL = REPO / "cases/recipe_trial"
TAG = "rae2822_case09"
ARMS = ("wolf", "tutorial", "lts", "lts_ll", "alletto", "hybrid", "alletto_tvd", "alletto_lim", "hybrid_krylov",
        "hybrid_tight", "hybrid_relax", "hybrid_unbounded", "wolf_gamg", "alletto_smooth", "wolf_damped", "wolf_hdamp", "simple_tvd", "alletto_elim",
        "sa_elim", "sa_lim", "sa_wolf", "sa_alletto",
        "alletto_k03", "alletto_flin2", "alletto_vanleer",
        "sa_k03", "sa_flin2", "sa_vanleer",
        "sa_flin2_t", "sa_flin2_k05")
SOLVER = {"wolf": "rhoSimpleFoam", "tutorial": "rhoSimpleFoam",
          "lts": "rhoPimpleFoam", "lts_ll": "rhoPimpleFoam",
          "alletto": "rhoSimpleFoam", "hybrid": "rhoSimpleFoam",
          "alletto_tvd": "rhoSimpleFoam", "alletto_lim": "rhoSimpleFoam",
          "hybrid_krylov": "rhoSimpleFoam",
          "hybrid_tight": "rhoSimpleFoam", "hybrid_relax": "rhoSimpleFoam",
          "hybrid_unbounded": "rhoSimpleFoam",
          "wolf_gamg": "rhoSimpleFoam",
          "alletto_smooth": "rhoSimpleFoam",
          "wolf_damped": "rhoSimpleFoam", "wolf_hdamp": "rhoSimpleFoam",
          "simple_tvd": "rhoSimpleFoam",
          "alletto_elim": "rhoSimpleFoam", "sa_elim": "rhoSimpleFoam",
          "sa_lim": "rhoSimpleFoam", "sa_wolf": "rhoSimpleFoam", "sa_alletto": "rhoSimpleFoam", "alletto_k03": "rhoSimpleFoam",
          "alletto_flin2": "rhoSimpleFoam", "alletto_vanleer": "rhoSimpleFoam",
          "sa_k03": "rhoSimpleFoam", "sa_flin2": "rhoSimpleFoam",
          "sa_vanleer": "rhoSimpleFoam",
          "sa_flin2_t": "rhoSimpleFoam", "sa_flin2_k05": "rhoSimpleFoam"}
# AGARD AR-138 case 9, Cook/McDonald/Firmin. Quoted here only for the printed comparison; the
# Cp RMS and shock position come from the tap file itself via report_figures_all.
AGARD = dict(CD=0.0168, CN=0.803, CM=-0.099, shock=0.550, M=0.730, Re=6.5e6, alpha=2.79)


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def sh(cmd, cwd, timeout=3600):
    return subprocess.run(["bash", "-lc", FOAM + cmd], cwd=str(cwd),
                          capture_output=True, text=True, timeout=timeout)


def case_dir(tag, arm):
    """cases/recipe_trial/<turbulence model>/<recipe>/<case> (project decision, 2026-08-11).

    61 flat directories with names like `argus_baseline_freestream_hybrid_krylov_am0p33`, in
    which the TURBULENCE MODEL appeared nowhere -- it was implied only by whether the recipe
    happened to start with `sa_`. Putting the model in the PATH means it cannot be misread,
    which is the same reason cases/argus2d is <wall treatment>/<geometry>/<convention>.
    """
    model = "SpalartAllmaras" if arm.startswith("sa_") else "kOmegaSST"
    return TRIAL / model / arm / tag


def build(arm, ranks, iters, potential=False):
    """Copy the mesh, reset to a COLD 0/, apply the recipe, decompose, audit. Halt on any."""
    d = case_dir(TAG, arm)
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for sub in ("constant", "system"):
        shutil.copytree(SRC / sub, d / sub)
    # THE COLD START COMES FROM 0/, AND NOT FROM 0_org, WHICH IS A TRAP.
    #
    # I wrote this the other way round first, on the reasonable-sounding ground that 0_org is
    # "the untouched initial condition". IT IS UNTOUCHED BY US AND IT IS THE WRONG FLOW: it is
    # the wolfdynamics NACA 0012 tutorial's condition, M 0.699 / p 73,048 Pa / T 283.24 K / U
    # 236 m/s, which build_rae2822_wolf.py copies in wholesale and then never updates, because
    # it writes the operating point into 0/ only. Case 9 is M 0.730 / 100,000 Pa / 288.15 K /
    # 248.69 m/s. Reset from 0_org and the case runs a DIFFERENT AEROFOIL'S OPERATING POINT
    # while every dictionary, every force coefficient and every figure still says case 9.
    #
    # CAUGHT BY THE PRE-RUN AUDIT, on its first real use, on all three arms at once
    # (magUInf 248.6936 vs |0/U| 236.0662; rhoInf 1.206270 vs p/(RT) 0.896431). That is the
    # gate doing precisely what ARG-096 item 10 says an INPUT gate is for: nothing downstream
    # could have seen it, because a converged solution to the wrong freestream is a perfectly
    # well-formed result. The three committed RAE cases still carry the stale 0_org.
    shutil.copytree(SRC / "0", d / "0")
    for stale in ("system/blockMeshDict",):
        if (d / stale).exists():
            (d / stale).unlink()

    sys.path.insert(0, str(REPO / "scripts"))
    import solver_recipe
    # THE TURBULENCE MODEL IS SWAPPED BEFORE THE RECIPE IS WRITTEN, because the recipe's
    # divSchemes and solver groups name the turbulence field: an SA recipe verified against a
    # case still carrying k and omega would pass on entries the solver never reads.
    if arm.startswith("sa_"):
        solver_recipe.to_spalart_allmaras(d)
    solver_recipe.apply(d, arm, iters=iters)          # writes AND verifies, raises otherwise


    if potential:
        # ---- POTENTIALFOAM INITIALISATION (project decision, 2026-08-11: try
        # potentialFoam first). Every wolf-family failure on ARGUS dies in the STARTUP -- 1,773, 1,958,
        # 2,022 and 2,880 iterations -- never reaching a developed flow, and heavier
        # relaxation only moved the death later rather than preventing it. That is the
        # signature of the INITIAL CONDITION, not of the converged state, and no amount of
        # relaxation tuning addresses it.
        #
        # WHAT THE COLD START ACTUALLY DOES: it applies no-slip INSTANTANEOUSLY to a uniform
        # stream. The first iteration therefore contains an infinite shear layer at the wall
        # and a pressure field that knows nothing about the body. At Re 1.68e7 the near-wall
        # cells are 1.05e-04 c, so that impulsive shear is resolved onto a very fine grid and
        # produces exactly the energy-equation excursion we measure.
        #
        # potentialFoam solves a Laplacian for the velocity potential, so it delivers a
        # DIVERGENCE-FREE, irrotational field that already flows AROUND the aerofoil with
        # roughly the right stagnation and acceleration pattern. The RANS solve then only has
        # to grow a boundary layer on a sensible field instead of inventing the whole flow.
        #
        # ONLY U IS TAKEN. potentialFoam's pressure is kinematic and inviscid, and this is a
        # COMPRESSIBLE case where p is absolute in Pa; writing its p would be a frame error of
        # the kind this project keeps paying for. p and T stay at their uniform freestream
        # values, which are correct, and only the velocity field is improved.
        fv = d / "system/fvSolution"
        s = fv.read_text()
        if "potentialFlow" not in s:
            s += "\npotentialFlow\n{\n    nNonOrthogonalCorrectors 5;\n}\n"
            fv.write_text(s)
        # ---- Phi's BCs ARE NOT COSMETIC AND `calculated` IS NOT ONE OF THEM.
        # My first version used `calculated` on every patch and potentialFoam aborted:
        # "cannot be called for a calculatedFvPatchField ... you are probably trying to solve
        # for a field with a default boundary condition." Phi is SOLVED, so it needs real
        # conditions, and they have to mirror how the flux is specified:
        #   wall and inlet   the normal flux is IMPOSED by U, so Phi is zeroGradient
        #   outlet           the flux is free and the potential is what anchors the problem,
        #                    so Phi is fixedValue 0. It also supplies the reference level that
        #                    an all-Neumann Laplacian would otherwise leave singular.
        # Written from the mesh's OWN patch list rather than a guessed set of names.
        b = (d / "constant/polyMesh/boundary").read_text(errors="replace")
        patches = re.findall(r"\n    (\w+)\s*\n    \{[^}]*?type\s+(\w+);", b, re.S)
        bf = ""
        for name, ptype in patches:
            if ptype == "empty":
                bf += "    %s\n    {\n        type            empty;\n    }\n" % name
            elif name == "outlet":
                bf += ("    %s\n    {\n        type            fixedValue;\n"
                       "        value           uniform 0;\n    }\n" % name)
            else:
                bf += "    %s\n    {\n        type            zeroGradient;\n    }\n" % name
        (d / "0/Phi").write_text(
            "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            "    class volScalarField;\n    object Phi;\n}\n\n"
            "dimensions      [0 2 -1 0 0 0 0];\ninternalField   uniform 0;\n\n"
            "boundaryField\n{\n%s}\n" % bf)
        v = d / "system/fvSolution"
        sv = v.read_text()
        if '"Phi"' not in sv and "\n    Phi\n" not in sv:
            sv = sv.replace("solvers\n{\n", "solvers\n{\n    Phi\n    {\n"
                            "        solver          GAMG;\n        smoother        "
                            "GaussSeidel;\n        tolerance       1e-08;\n"
                            "        relTol          0.01;\n    }\n", 1)
            v.write_text(sv)
        r0 = sh("potentialFoam -writePhi > log.potentialFoam 2>&1", d, 3600)
        # DERIVED FROM THE FIELD, NOT FROM THE EXIT STATUS. A potentialFoam that ran and
        # changed nothing leaves a uniform 0/U and a cold start wearing a warm label.
        if "nonuniform" not in (d / "0/U").read_text(errors="replace"):
            raise SystemExit("%s: potentialFoam left 0/U uniform, so it did nothing" % arm)
        # ---- AND ITS FLUX MUST GO. potentialFoam also writes 0/phi, a VOLUMETRIC flux with
        # dimensions [0 3 -1]; rhoSimpleFoam is compressible and needs a MASS flux, [1 0 -1].
        # Left in place, the solver refused the case outright with
        #     [U[0 1 -2]] + [(rho*MRFZoneList:acceleration)[1 -2 -2]]
        # which is the momentum equation reporting that its two halves no longer have the same
        # dimensions. Caught by the audit's one-step solver check, which exists precisely for
        # the failures the enumerated checks did not anticipate.
        # Deleting it costs nothing: rhoSimpleFoam recomputes phi from rho and U on startup,
        # so the improved VELOCITY field is inherited and the inconsistent flux is not. That
        # is what "only U is taken" has to mean in practice rather than in intent.
        for stale in ("0/phi", "0/Phi"):
            if (d / stale).exists():
                (d / stale).unlink()
        log("  %-9s potentialFoam initialised 0/U (now nonuniform)" % arm)

    # ---- BEFORE decomposePar, AND THAT ORDERING IS THE WHOLE POINT.
    # Run after decomposition, potentialFoam writes 0/U while the parallel solve reads
    # processor*/0/U, so the improved field is silently discarded. The run then came back
    # BIT-IDENTICAL to the cold one -- dead at the same iteration 1,958 with the same 290
    # energy-residual events -- which is the same signature as the ARG-156 "damped" runs that
    # matched the undamped ones exactly. An initialisation that is ignored is indistinguishable
    # from one that did not help, and only the identical iteration COUNT gave it away.
    dp = d / "system/decomposeParDict"
    dp.write_text(re.sub(r"numberOfSubdomains\s+\d+;", "numberOfSubdomains %d;" % ranks,
                         dp.read_text()))
    r = sh("rm -rf processor*; decomposePar -force > log.decomposePar 2>&1", d, 3600)
    n = len(list(d.glob("processor*")))
    if n != ranks:
        raise SystemExit("%s: %d subdomains vs %d requested\n%s" % (arm, n, ranks, r.stdout))

    au = [sys.executable, "scripts/audit_case_settings.py", str(d)]
    if arm.startswith("sa_"):
        au += ["--allow-rasmodel", "SpalartAllmaras", "--reason",
               "Authorised by project decision 2026-08-11: kOmegaSST clips k/omega 45,001-62,527 times per "
               "60,000 iterations on every arm; omega ~ 6nu/(beta y^2) is stiff at the wall "
               "and SA has no omega equation. Pending ARG-161; nothing from an SA run is "
               "quotable until that entry exists."]
    a = subprocess.run(au, cwd=str(REPO), capture_output=True, text=True, timeout=1800)
    (d / "log.audit").write_text(a.stdout + a.stderr)
    if a.returncode != 0:
        # REPORT FIRST, EXIT AFTER (GEO-089 item 3). An audit failure that hides which check
        # failed is an unactionable alarm.
        print(a.stdout)
        raise SystemExit("%s: audit failed, not launching" % arm)
    log("  %-9s built, recipe verified, %d ranks, audit clean" % (arm, ranks))
    return d


def continue_run(arm, iters):
    """Restart an arm from its latest time and carry it to `iters`.

    WHY (review of the LTS residual plot, 2026-08-10): the forces suggested the case needs to run
    longer. Correct, and the plot says something sharper than "longer". The residuals are FLAT at 2e-4 from iteration 2,500 onward while C_d wanders
    +/-8 counts with a period of roughly 5,000 iterations. THAT IS A LIMIT CYCLE, NOT A DECAYING
    TRANSIENT, so a residual threshold will never fire and 15,000 iterations covers only about
    two periods. An oscillating force converges to its MEAN (ARG-157), and a mean over two
    periods is not a mean -- it is a sample of the phase the run happened to stop in. This is
    the same aliasing that made D023 unquotable for these cases: measured 1.5 to 16.5 counts on
    one run depending only on where the window ended.

    THE RECIPE IS RE-VERIFIED BEFORE THE RESTART, not assumed to have survived. Continuation
    edits controlDict, and a run that silently resumed under different numerics would be
    indistinguishable in the log from one that did not.
    """
    d = case_dir(TAG, arm)
    sys.path.insert(0, str(REPO / "scripts"))
    import solver_recipe
    bad = solver_recipe.verify(d, arm, verbose=False)
    if bad:
        raise SystemExit("%s: recipe no longer verifies, refusing to continue:\n  %s"
                         % (arm, "\n  ".join(bad)))
    c = d / "system/controlDict"
    s = c.read_text()
    s = re.sub(r"^startFrom\s+\S+;", "startFrom       latestTime;", s, flags=re.M)
    s = re.sub(r"^endTime\s+\S+;", "endTime         %d;" % iters, s, flags=re.M)
    c.write_text(s)
    got = re.search(r"^startFrom\s+(\S+);", c.read_text(), re.M).group(1)
    end = re.search(r"^endTime\s+(\S+);", c.read_text(), re.M).group(1)
    if got != "latestTime" or int(end) != iters:
        raise SystemExit("%s: continuation edit did not land (%s / %s)" % (arm, got, end))
    have = sorted(int(p.name) for p in (d / "processor0").iterdir()
                  if p.is_dir() and p.name.isdigit())
    log("  %-9s continuing from t=%d to %d" % (arm, have[-1] if have else 0, iters))
    return d


def crash_markers(case):
    """Read the run's own outcome out of its log. Nothing is asserted by the launcher.

    THE h RESIDUAL HITTING 1.0 IS THE CRASH TEST (project decision, 2026-08-10): an h residual
    of 1 means the run has crashed; hitting limits means it exploded. A residual
    of exactly 1 means the initial residual equalled the normalisation, i.e. the update was as
    large as the field. Counted, not thresholded.
    """
    case = pathlib.Path(case)
    # READ THE COMPRESSED LOG TOO. scripts/reclaim_storage.py gzips solver logs to free space
    # in the WSL image, and a reader that only knows log.solver would report "never launched"
    # for a run that completed perfectly -- a storage decision silently becoming a WRONG
    # RESULT rather than a missing one. Concatenated because a continued run leaves
    # log.solver.0_15000 beside log.solver and the crash markers span both.
    # ORDER-INDEPENDENT, AFTER TWO ORDERING ASSUMPTIONS FAILED IN A ROW.
    #   attempt 1, sort by NAME: "log.solver" sorts before "log.solver.0_15000", so the
    #     ARCHIVED first segment landed last and its "Time = 15000" overwrote the
    #     continuation's 60,000 -- a naming convention I invented deciding a physical number.
    #   attempt 2, sort by MTIME: correct until scripts/reclaim_storage.py gzipped the logs
    #     and stamped both copies with the compression time. The SAME BUG CAME BACK from an
    #     unrelated operation, which is the tell that the fix was a different assumption
    #     rather than the removal of one.
    # So no ordering is assumed at all: every segment is parsed, iteration count is the MAX
    # over all of them, wall clock is the SUM, and the crash markers are a union. All three
    # are order-free by construction, and the marker counts were already order-free.
    parts = []
    for lg in sorted(case.glob("log.solver*")):
        if lg.suffix == ".gz":
            import gzip
            parts.append(gzip.open(lg, "rt", errors="replace").read())
        else:
            parts.append(lg.read_text(errors="replace"))
    if not parts:
        return dict(launched=False, why="no log.solver")
    t = "\n".join(parts)
    out = dict(launched=True)
    out["converged"] = "SIMPLE solution converged" in t
    # ANCHORED TO THE ACTUAL FAILURE MARKERS, NOT TO THE WORDS.
    # The first version matched "floating point exception" and therefore matched
    #     sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
    # which OpenFOAM prints in the BANNER OF EVERY RUN. So every arm, including the one that
    # completed 60,000 iterations cleanly, was reported as having raised a FOAM FATAL and the
    # partition came back with zero survivors. A GATE THAT ALWAYS FIRES is as useless as one
    # that cannot, and this is the fifth substring-instead-of-entry failure tonight.
    out["foam_fatal"] = bool(re.search(
        r"^-->\s*FOAM FATAL|Foam::error::printStack|Foam::sigFpe::sigHandler"
        r"|\*\*\* Process received signal \*\*\*|exited on signal", t, re.M))
    # THE FIRST ITERATION IS EXCLUDED, AND THAT IS NOT A LOOSENING.
    # OpenFOAM normalises the initial residual by the field's own variation. On a UNIFORM
    # initial field that normalisation is degenerate and h reads exactly 1 on iteration 1 for
    # every cold start; it is an artefact of the initial condition, not of the solve. Counting
    # it marked the alletto arm as "crashed, not slow" on the strength of one sample taken
    # before the solver had done anything. The project criterion is about a residual that reaches
    # 1 DURING a run, which is what is counted here.
    h = [float(x) for x in re.findall(r"for h, Initial residual = ([-0-9.eE+]+)", t)]
    out["h_residual_at_one"] = int(sum(1 for x in h[1:] if x >= 1.0))
    out["h_residual_at_one_including_startup"] = int(sum(1 for x in h if x >= 1.0))
    out["h_residual_final"] = h[-1] if h else None
    out["clamp_events"] = len(re.findall(r"pressureControl:", t))
    out["clamp_min_events"] = len(re.findall(r"pressureControl: p min", t))
    out["clamp_max_events"] = len(re.findall(r"pressureControl: p max", t))
    out["bounding_events"] = len(re.findall(r"bounding (k|omega|epsilon|h|e),", t))
    # ---- A DIVERGING LINEAR SOLVE IS ITS OWN FAILURE MODE, distinct from the physics.
    # Arm E died with "Solving for Ux ... Final residual = 7.0e+101, No Iterations 1000" while
    # every physical field was healthy. Reported as a crash it looks like an unstable recipe;
    # it is actually GAMG diverging on a convection-dominated, aspect-ratio-1e5 momentum
    # matrix, and the fix is the linear solver, not the numerics. A gate that cannot tell
    # those apart sends you to change the wrong thing.
    solves = re.findall(r"Solving for (\w+), Initial residual = ([-0-9.eE+]+), "
                        r"Final residual = ([-0-9.eE+]+), No Iterations (\d+)", t)
    div = [(f, float(i), float(fi), int(n)) for f, i, fi, n in solves if float(fi) > float(i)]
    out["linear_solver_diverged"] = len(div)
    out["linear_solver_diverged_worst"] = (
        max(div, key=lambda q: q[2])[:1] + (max(q[2] for q in div),) if div else None)
    out["linear_solver_maxiter"] = sum(1 for f, i, fi, n in solves if int(n) >= 1000)
    # SUM the wall clock over segments (each restarts its own ExecutionTime at zero) and take
    # the MAX iteration, never the last of either.
    seg_clocks = []
    for seg in parts:
        c = re.findall(r"ClockTime = (\d+)", seg)
        seg_clocks.append(int(c[-1]) if c else 0)
    out["clock_s"] = sum(seg_clocks) or None
    out["clock_s_per_segment"] = seg_clocks
    it = [int(x) for x in re.findall(r"^Time = (\d+)", t, re.M)]
    out["iterations"] = max(it) if it else 0
    yp = re.findall(r"patch aerofoil y\+ : min = ([-0-9.eE+]+), max = ([-0-9.eE+]+), "
                    r"average = ([-0-9.eE+]+)", t)
    if yp:
        out["yplus"] = dict(zip(("min", "max", "avg"), (float(x) for x in yp[-1])))
    return out


def tail_stats(case, frac=0.25):
    """Drift of the MEAN and the amplitude, over the last `frac` of the force history.

    D023 IS NOT USED HERE AND THAT IS DELIBERATE. Over a 2,000-sample window it ALIASES on an
    oscillating force: measured 1.5 to 16.5 counts on one and the same run depending only on
    where the window ended. Two separate numbers instead -- how far the mean is still moving,
    and how wide the oscillation is -- because they answer different questions and collapsing
    them into one convergence index is what made D023 unquotable for these cases.
    """
    case = pathlib.Path(case)
    f = sorted(case.glob("postProcessing/forces_coeffs/*/forceCoeffs*.dat"))
    if not f:
        return None
    rows = []
    for l in open(f[-1], errors="replace"):
        if l.startswith("#"):
            continue
        v = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", l)]
        if len(v) >= 5:
            rows.append(v)
    if len(rows) < 40:
        return None
    w = min(len(r) for r in rows)
    a = np.array([r[:w] for r in rows])
    n = len(a)
    lo = int(n * (1 - frac))
    out = {}
    for name, col in (("cd", 2), ("cl", 3)):
        s = a[lo:, col]
        half = len(s) // 2
        out[name] = dict(
            final=float(s[-1]),
            tail_mean=float(s.mean()),
            # drift: first half of the tail against the second half, in COUNTS for cd
            drift=float(abs(s[half:].mean() - s[:half].mean())),
            amplitude=float(s.max() - s.min()))
    out["samples"] = int(n)
    out["tail_from_iter"] = float(a[lo, 0])
    return out


def measure(arm, d):
    """Everything about one arm, DERIVED from its own case directory."""
    sys.path.insert(0, str(REPO / "scripts"))
    import wallresolved_transonic_routes as W
    rec = dict(arm=arm, case=str(d.relative_to(REPO)), solver=SOLVER[arm])
    rec["run"] = crash_markers(d)
    rec["tail"] = tail_stats(d)
    try:
        rec["forces"] = W.forces(d)
    except Exception as e:
        rec["forces"] = {"error": "%s: %s" % (type(e).__name__, e)}
    rec["gate_temperature"] = W.temperature_ceiling(d)
    rec["gate_pressure_drag"] = W.pressure_drag_sign(d)
    return rec


def figures(cases, iters):
    """Cp/Cf/y+/shock/Cp-RMS against AGARD, and the figures, via the existing pipeline."""
    r = subprocess.run([sys.executable, "scripts/report_figures_all.py"]
                       + [str(c) for c in cases] + ["--out-dir", "results/figures/recipe_trial"],
                       cwd=str(REPO), capture_output=True, text=True, timeout=7200)
    print(r.stdout[-3000:])
    if r.returncode != 0:
        log("  report_figures_all returned %d\n%s" % (r.returncode, r.stderr[-1500:]))
    out = {}
    for c in cases:
        j = REPO / "results/figures/recipe_trial" / (pathlib.Path(c).name + ".json")
        if j.exists():
            out[pathlib.Path(c).name] = json.loads(j.read_text())
    return out


def verdict(recs):
    """Apply the PRE-REGISTERED rule. The rule is quoted here so it cannot drift to fit."""
    surv, out_ = [], []
    for r in recs:
        why = []
        if not r["run"].get("launched"):
            why.append("never launched: %s" % r["run"].get("why"))
        if r["run"].get("foam_fatal"):
            why.append("solver raised a FOAM FATAL")
        if r["run"].get("h_residual_at_one", 0) > 0:
            why.append("energy residual reached 1.0 on %d iterations (crashed, not slow)"
                       % r["run"]["h_residual_at_one"])
        g = r.get("gate_temperature") or {}
        if g.get("exceeded_T0"):
            why.append("T exceeded the adiabatic ceiling %.1f K at iteration %s"
                       % (g.get("T0", float("nan")), g.get("exceeded_at")))
        p = r.get("gate_pressure_drag") or {}
        if p.get("cd_press_went_negative"):
            why.append("C_d,pressure went negative at iteration %s (d'Alembert violation)"
                       % p.get("cd_press_neg_at"))
        (out_ if why else surv).append((r, why))
    return surv, out_


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--arms", nargs="*", default=list(ARMS), choices=list(ARMS))
    ap.add_argument("--ranks", type=int, default=6)
    ap.add_argument("--src", default=None,
                    help="case directory to take the mesh, geometry and operating point from")
    ap.add_argument("--tag", default=None, help="prefix for the trial case directories")
    ap.add_argument("--iters", type=int, default=15000,
                    help="15000 matches the existing case09 run, so arm A is comparable "
                         "against a number already in the record")
    ap.add_argument("--build-only", action="store_true")
    ap.add_argument("--potential", action="store_true",
                    help="initialise 0/U with potentialFoam before the solve")
    ap.add_argument("--continue-to", type=int, default=None,
                    help="restart the named arms from their latest time and carry them to this "
                         "iteration; the recipe is re-verified first")
    ap.add_argument("--collect-only", action="store_true",
                    help="skip build and solve; re-measure what is on disk")
    a = ap.parse_args()

    global SRC, TAG
    if a.src:
        SRC = REPO / a.src
    if a.tag:
        TAG = a.tag
    if not (SRC / "constant/polyMesh/points").exists():
        raise SystemExit("no mesh in %s" % SRC)
    TRIAL.mkdir(parents=True, exist_ok=True)
    built, failed = {}, {}
    if a.continue_to:
        for arm in a.arms:
            try:
                built[arm] = continue_run(arm, a.continue_to)
            except SystemExit as e:
                failed[arm] = str(e)
                log("  %-9s NOT CONTINUED: %s" % (arm, e))
    elif not a.collect_only:
        for arm in a.arms:
            try:
                built[arm] = build(arm, a.ranks, a.iters, a.potential)
            except SystemExit as e:
                failed[arm] = str(e)
                log("  %-9s NOT BUILT: %s" % (arm, e))
    else:
        for arm in a.arms:
            d = case_dir(TAG, arm)
            (built if d.exists() else failed).__setitem__(arm, d if d.exists() else "absent")
    if a.build_only:
        return 0 if not failed else 1

    if not a.collect_only:
        procs = {}
        for arm, d in built.items():
            cmd = ("mpirun -np %d %s -parallel > log.solver 2>&1"
                   % (a.ranks, SOLVER[arm]))
            procs[arm] = subprocess.Popen(["bash", "-lc", FOAM + cmd], cwd=str(d),
                                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            log("  %-9s launched: %s" % (arm, cmd))
        t0 = time.time()
        while any(p.poll() is None for p in procs.values()):
            time.sleep(120)
            st = []
            for arm, p in procs.items():
                m = crash_markers(built[arm])
                st.append("%s %s%s" % (arm, m.get("iterations", 0),
                                       "" if p.poll() is None else "(done)"))
            log("  %5.0f min  %s" % ((time.time() - t0) / 60, "   ".join(st)))
        for arm, p in procs.items():
            log("  %-9s exit %s" % (arm, p.returncode))

    recs = [measure(arm, built[arm]) for arm in a.arms if arm in built]
    figs = figures([built[arm] for arm in a.arms if arm in built], a.iters)
    for r in recs:
        r["figures"] = figs.get(pathlib.Path(r["case"]).name, {})

    surv, out_ = verdict(recs)

    print("\n" + "=" * 100)
    print("PHASE 0 -- RAE 2822 case 9, M %.3f, Re %.2e, alpha %.2f deg. AGARD: C_D %.4f, "
          "C_N %.3f, shock x/c %.3f" % (AGARD["M"], AGARD["Re"], AGARD["alpha"],
                                        AGARD["CD"], AGARD["CN"], AGARD["shock"]))
    print("=" * 100)
    hdr = ("arm", "solver", "iters", "clock", "C_d", "d%", "C_l", "l%", "shock",
           "dshock", "CpRMS", "h=1", "clamp", "bnd")
    print("%-9s %-14s %7s %7s %9s %7s %8s %7s %7s %7s %7s %5s %6s %6s" % hdr)
    for r in recs:
        f = r.get("forces") or {}
        fg = r.get("figures") or {}
        run = r["run"]
        cd, cl = f.get("cd"), f.get("cl")
        sk = fg.get("shock_cfd")
        print("%-9s %-14s %7s %7s %9s %7s %8s %7s %7s %7s %7s %5s %6s %6s" % (
            r["arm"], r["solver"], run.get("iterations", "-"),
            "%.0fm" % (run["clock_s"] / 60) if run.get("clock_s") else "-",
            "%.6f" % cd if cd is not None else "-",
            "%+.1f" % (100 * (cd / AGARD["CD"] - 1)) if cd else "-",
            "%.5f" % cl if cl is not None else "-",
            "%+.1f" % (100 * (cl / AGARD["CN"] - 1)) if cl else "-",
            "%.3f" % sk if sk else "-",
            "%+.3f" % (sk - AGARD["shock"]) if sk else "-",
            "%.4f" % (fg.get("cp_rms_upper") or fg.get("cp_rms_vs_ref"))
            if (fg.get("cp_rms_upper") or fg.get("cp_rms_vs_ref")) else "-",
            run.get("h_residual_at_one", "-"), run.get("clamp_events", "-"),
            run.get("bounding_events", "-")))

    print("\nTAIL BEHAVIOUR (last 25%% of the force history; D023 deliberately not used, it "
          "aliases on an oscillating force)")
    for r in recs:
        t = r.get("tail")
        if not t:
            print("  %-9s no force history" % r["arm"]); continue
        print("  %-9s from it %6.0f   C_d mean %.6f  drift %.2f counts  amplitude %.2f counts"
              "   |  C_l mean %.5f  drift %.5f  amplitude %.5f"
              % (r["arm"], t["tail_from_iter"], t["cd"]["tail_mean"],
                 1e4 * t["cd"]["drift"], 1e4 * t["cd"]["amplitude"],
                 t["cl"]["tail_mean"], t["cl"]["drift"], t["cl"]["amplitude"]))

    print("\nPHYSICAL GATES")
    for r in recs:
        g, p = r.get("gate_temperature") or {}, r.get("gate_pressure_drag") or {}
        print("  %-9s T_max %.1f K vs T0 %.1f K (%.3fx)%s   |   C_d,press %+.6f  C_d,visc "
              "%+.6f%s" % (
                  r["arm"], g.get("T_max_reached", float("nan")), g.get("T0", float("nan")),
                  g.get("T_max_over_T0", float("nan")),
                  "  *** EXCEEDED ***" if g.get("exceeded_T0") else "",
                  p.get("cd_press_final", float("nan")), p.get("cd_visc_final", float("nan")),
                  "  *** WENT NEGATIVE ***" if p.get("cd_press_went_negative") else ""))

    print("\nPARTITION OF THE %d ARMS (GEO-089: nothing may be in neither bucket)" % len(ARMS))
    print("  SURVIVED THE PHYSICAL GATES : %s" % ", ".join(r["arm"] for r, _ in surv) or "none")
    for r, why in out_:
        print("  OUT  %-9s %s" % (r["arm"], "; ".join(why)))
    for arm, why in failed.items():
        print("  NOT BUILT  %-9s %s" % (arm, why))
    accounted = len(surv) + len(out_) + len(failed)
    print("  accounted %d of %d requested%s"
          % (accounted, len(a.arms), "" if accounted == len(a.arms) else "   *** GAP ***"))

    print("\nPRE-REGISTERED DECISION RULE (from the approved plan; NOT re-derived from these "
          "numbers)")
    print("  1. any arm tripping a physical gate is OUT regardless of its forces")
    print("  2. among survivors: shock within one tap spacing (0.025 c) of 0.550, then Cp RMS,"
          " then C_d")
    print("  3. A must win on accuracy by a CLEAR MARGIN to be retained")
    print("  4. if C matches B on accuracy, C is preferred for the 3D path")

    p = REPO / "results/recipe_trial_case09.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(
        agard=AGARD, iters=a.iters, ranks=a.ranks, arms=recs,
        survived=[r["arm"] for r, _ in surv],
        out=[{"arm": r["arm"], "why": w} for r, w in out_],
        not_built=failed), indent=2, default=str) + "\n")
    print("\n  written %s" % p.relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
