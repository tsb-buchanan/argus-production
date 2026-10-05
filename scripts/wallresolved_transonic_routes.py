#!/usr/bin/env python3
"""wallresolved_transonic_routes.py: the two untried routes to a wall-resolved transonic solve.

WHY (project decision, 2026-08-09 night). Eight attempts at wall-resolved transonic have failed (ARG-155)
and every one shared a signature: k and omega clipped negative in the first hundred
iterations, C_d crossing zero, then divergence. Our own hpc/slurm/mesh_solve.slurm already
names the cause -- "applying no-slip at iteration 1 across a 4.5 um first cell is a shear rate
of 9.1e6 1/s, so nut*S^2 production explodes". Every attempt so far varied the MESH or the
SCHEMES. Neither touches the cold start.

These two do:

  ROUTE 1, MACH RAMP. We HAVE a working wall-resolved solution -- naca0012_te60 at M 0.15,
  C_d -0.02% against CFL3D. Walk the Mach number up in small steps, warm-continuing each from
  the last converged field. The solver never sees a cold start and never sees a transonic
  field it has not been walked into. The shock forms gradually on an already-developed
  boundary layer instead of appearing on iteration one.

  ROUTE 2, TRANSIENT MARCH. rhoPimpleFoam with a small timestep and adjustTimeStep. A steady
  solver takes an effectively infinite step damped only by relaxation; a transient one has a
  timestep, which is exactly the limiter the cold start lacks. Slower per iteration -- the
  tutorial spent 159,724 steps to reach 64% of a run -- but it cannot take the unbounded first
  step that kills the steady solver.

BOTH ARE JUDGED ON THE SAME KNOWN ANSWER. RAE 2822 case 9: C_d 0.0168, C_N 0.803, shock at
x/c 0.550, and ATTACHED. A route that produces separation on this case has failed, whatever
its residuals look like -- our wall-MODELLED setup already gets 0.0% reversed flow here, so
the bar is set and it is not a guess.

WHAT IS RECORDED EITHER WAY. Both routes write their trajectory, so a failure says WHERE it
left the physical solution rather than only that it did. C_d crossing zero is the marker; it
is unphysical and every previous attempt hit it before diverging.
"""
import argparse
import json
import math
import os
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
TARGET = dict(CD=0.0168, CN=0.803, shock=0.550, attached=True)
# AGARD f8621 corrected angles. Case 7 carries the WEAKEST shock of the three and sits
# furthest from the buffet boundary, which is why it is the discriminator for whether the
# wall-resolved drift on case 9 is a numerical failure or an incipient-unsteadiness one.
AGARD = {"case07": dict(CD=0.0107, CN=0.658, alpha=2.21, M=0.725),
         "case06": dict(CD=0.0127, CN=0.743, alpha=2.54, M=0.725),
         "case09": dict(CD=0.0168, CN=0.803, alpha=2.79, M=0.730)}


def sh(cmd, cwd=None, timeout=None):
    return subprocess.run(["bash", "-lc", FOAM + cmd], cwd=str(cwd) if cwd else None,
                          capture_output=True, text=True, timeout=timeout)


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def forces(case):
    sys.path.insert(0, str(REPO / "scripts"))
    from watch_run import read_dat, d023, bounding
    hdr, rows = read_dat(str(case), "forces_coeffs", "forceCoeffs*.dat")
    if not rows:
        return None
    rows = [[x for x in r if isinstance(x, float)] for r in rows]
    w = max(len(r) for r in rows)
    a = np.array([r for r in rows if len(r) == w])
    c = {n: i for i, n in enumerate(hdr)} if hdr else {}
    cd, cl = a[:, c.get("Cd", 2)], a[:, c.get("Cl", 3)]
    b, tail, tot = bounding(str(case))
    # THE UNPHYSICAL MARKER. Every failed attempt crossed C_d = 0 before diverging, so it is
    # reported as a state rather than left to be spotted in a plot afterwards.
    n = len(cd)
    neg = np.where(cd[n // 5:] < 0)[0]
    return dict(iters=int(a[-1, 0]), cd=float(cd[-1]), cl=float(cl[-1]),
                d023=float(d023(cd.tolist())), bnd_tail=int(tail),
                cd_went_negative=bool(len(neg)),
                cd_neg_at=float(a[n // 5 + neg[0], 0]) if len(neg) else None)


def temperature_ceiling(case, margin=0.05):
    """WHEN T_max FIRST EXCEEDS THE ADIABATIC STAGNATION TEMPERATURE. Earliest marker we have.

    T0 = T_inf (1 + (gamma-1)/2 M^2). In an adiabatic flow with no work input NOTHING can be
    hotter than that, so T_max > T0 is impossible, not merely large. Case 6 runs at M 0.725,
    T_inf 288.15, so T0 = 318.44 K.

    AND THE fvOptions limitTemperature ENTRY IS HIDING IT. The case ships `limitT` with
    min 100 / max 900, whose own comment reads "900-1000 with 800 becomes unbounded" -- i.e.
    the ceiling was tuned upward until runs stopped crashing. 900 K is 2.8x T0. The limiter
    converts a hard, immediate failure into a slow silent corruption that still writes
    plausible forces, which is the most expensive shape of defect this project has (GEO-080).
    Measured, sustained above T0 + 5%:
        case06  from it  1,763   limiter saturates  8,533   crash 12,834
        case07  from it 11,681   limiter saturates 13,459   crash 21,498
    """
    case = pathlib.Path(case)
    fs = sorted(case.glob("postProcessing/minmaxdomain/*/fieldMinMax.dat"))
    if not fs:
        return None
    # T_inf AND Mach FROM case_normalisation, WHICH HANDLES A WARM-STARTED CASE.
    # 0/T is `nonuniform List<scalar>` after mapFields, so reading `Tinlet` or an
    # `internalField uniform` value returns nothing and this check silently reported None --
    # a gate that goes quiet on exactly the runs it exists to judge.
    sys.path.insert(0, str(REPO / "scripts"))
    import case_normalisation as cn
    try:
        f_ = cn.derive(str(case))["flow"]
    except SystemExit:
        return None
    Tinf, M = f_["T_inf"], f_["Mach"]
    T0 = Tinf * (1 + 0.2 * M * M)
    t, tmax = [], []
    for l in open(fs[-1], errors="replace"):
        if l.startswith("#"):
            continue
        q = l.split("\t")
        if len(q) < 8 or q[1].strip() != "T":
            continue
        t.append(float(q[0]))
        tmax.append(float(re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", q[5])[0]))
    if not t:
        return None
    t, tmax = np.array(t), np.array(tmax)
    thr = T0 * (1 + margin)
    onset = None
    for i in np.where(tmax > thr)[0]:
        # SUSTAINED, not a single spike: a lone excursion at the stagnation point is noise.
        if (tmax[i:i + 200] > thr).mean() > 0.8:
            onset = float(t[i]); break
    return dict(T_inf=Tinf, Mach=float(M), T0=float(T0), threshold=float(thr),
                T_max_reached=float(tmax.max()), T_max_over_T0=float(tmax.max() / T0),
                exceeded_T0=onset is not None, exceeded_at=onset,
                iterations=float(t[-1]))


def pressure_drag_sign(case):
    """WHERE C_d,PRESSURE FIRST GOES NEGATIVE. A one-directional physical bound (D060).

    PRESSURE DRAG ON A LIFTING AEROFOIL WITH A SHOCK CANNOT BE NEGATIVE. Wave drag is
    strictly positive and form drag is non-negative, so C_d,press < 0 is a d'Alembert
    violation, not a large error. It is therefore a SIGN test and strictly stronger than any
    tolerance on the total.

    WHY IT MATTERS MORE THAN THE TOTAL, measured on case 7 (refined shock band):
        it  2000  press +0.006525  visc 0.004469  total 0.010994
        it  8000  press +0.003507  visc 0.004644  total 0.008152
        it 12000  press -0.001770  visc 0.004716  total 0.002946   <-- IMPOSSIBLE
    The run crashed at 21,498. The pressure drag went negative around 11,000, i.e. TEN
    THOUSAND ITERATIONS EARLIER, while the total C_d was still smooth to 3.4 counts
    peak-to-peak and every residual looked flat. The blow-up is the last symptom, not the
    fault, and the total C_d hides the crossing because the viscous part stays healthy and
    masks it.
    """
    case = pathlib.Path(case)
    cd_ = (case / "system/controlDict").read_text()
    g = lambda k: float(re.search(r"%s\s+([-0-9.eE+]+);" % k, cd_).group(1))
    try:
        A, rho, U = g("Aref"), g("rhoInf"), g("magUInf")
        dd = np.array([float(x) for x in
                       re.search(r"^\s*dragDir\s+\(([^)]*)\);", cd_, re.M).group(1).split()])
    except AttributeError:
        return None
    f = sorted(case.glob("postProcessing/forces/*/force*.dat"))
    if not f:
        return None
    rows = []
    for l in open(f[-1], errors="replace"):
        if l.startswith("#"):
            continue
        v = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", l)]
        if len(v) >= 19:
            rows.append(v[:19])
    if not rows:
        return None
    a = np.array(rows)
    q = 0.5 * rho * U * U
    cdp = (a[:, 1:4] @ dd) / (q * A)
    cdv = (a[:, 4:7] @ dd) / (q * A)
    neg = np.where(cdp < 0)[0]
    # ---- THE VERDICT IS THE SUSTAINED FRACTION AFTER STARTUP, NOT THE FIRST OCCURRENCE.
    # A cold start from a uniform freestream is not a flow yet: the first few hundred
    # iterations routinely put C_d,pressure below zero while the field is still nothing like
    # an aerofoil solution. Judged on first occurrence, this gate FIRES ON EVERY COLD START
    # and cannot tell them apart from a real d'Alembert violation. Measured on this project:
    #     hybrid_krylov, ARGUS normal   negative for 16 iterations of 60,000 (it 132-147),
    #                                   0.0% after startup, then 60,000 iterations at 0.01
    #                                   counts of drift -- a perfectly healthy run
    #     tutorial, RAE case 9          negative for 46.2% OF THE WHOLE RUN, it 4 to 22,076 --
    #                                   a genuine violation
    # First occurrence rates those identically. The fraction after the first 20% separates
    # them completely, and it is the same startup exclusion forces() already uses.
    # THE FIRST-OCCURRENCE ITERATION IS KEPT as a diagnostic, because ARG-158's value was that
    # the crossing PRECEDES the visible failure by thousands of iterations. It is information
    # about when to look, not the verdict.
    n = len(cdp)
    tail = cdp[n // 5:]
    frac = float((tail < 0).mean())
    return dict(cd_press_final=float(cdp[-1]), cd_visc_final=float(cdv[-1]),
                cd_press_went_negative=bool(frac > 0),
                cd_press_neg_fraction_after_startup=frac,
                cd_press_neg_fraction_whole_run=float((cdp < 0).mean()),
                cd_press_neg_iterations=int(len(neg)),
                cd_press_first_neg_at=float(a[neg[0], 0]) if len(neg) else None,
                cd_press_neg_at=float(a[neg[0], 0]) if len(neg) else None,
                last_physical_iteration=float(a[neg[0] - 1, 0]) if len(neg) else float(a[-1, 0]))


def sep_fraction(case):
    """Reversed-flow fraction on the upper surface. The experiment says 0 for case 9."""
    case = pathlib.Path(case)
    sh("reconstructPar -latestTime > log.recon 2>&1", case, 3600)
    sh("rhoSimpleFoam -postProcess -func wallShearStress -latestTime > log.wss 2>&1", case, 3600)
    (case / "system/surfaces").write_text(
        'type surfaces;\nlibs ("libsampling.so");\nsurfaceFormat raw;\n'
        'interpolationScheme cell;\nfields (p wallShearStress yPlus);\n'
        'surfaces ( airfoil { type patch; patches (".*(wall|aerofoil|airfoil).*"); '
        'interpolate false; } );\n')
    sh("postProcess -func surfaces -latestTime > log.sample 2>&1", case, 3600)
    g = sorted(case.glob("postProcessing/surfaces/*/wallShearStress_*.raw"))
    if not g:
        return None
    sys.path.insert(0, str(REPO / "scripts"))
    from plot_case_report import split_surfaces
    q = np.loadtxt(g[-1])
    (xu, tu), _ = split_surfaces(q[:, 0], q[:, 1], q[:, 3])
    return round(100 * float((tu > 0).mean()), 2)


def set_mach(case, M, Re, T=288.15, p_inf=1.0e5):
    """Retarget a case to a new Mach at the same Reynolds number."""
    case = pathlib.Path(case)
    tp = case / "constant/thermophysicalProperties"
    tt = tp.read_text()
    R = 8314.47 / float(re.search(r"^\s*molWeight\s+([-0-9.eE+]+);", tt, re.M).group(1))
    a = math.sqrt(1.4 * R * T)
    U = M * a
    rho = p_inf / (R * T)
    mu = rho * U * 1.0 / Re
    tp.write_text(re.sub(r"^(\s*mu\s+)[-0-9.eE+]+;", r"\g<1>%.8e;" % mu, tt, flags=re.M))
    u = case / "0/U"
    ut = u.read_text()
    m = re.search(r"^Uinlet\s+\(([^)]*)\);", ut, re.M)
    # FAIL, DO NOT FALL THROUGH. The direction comes from Uinlet and everything below reuses
    # it; without the macro this used to raise NameError twenty lines later, in a message that
    # named neither the file nor the reason (GEO-092: a bypass that does not assert its own
    # precondition is a bypass without conditions).
    if not m:
        raise SystemExit("%s/0/U carries no `Uinlet (...)` macro: cannot retarget the Mach "
                         "without knowing the freestream direction" % case)
    v = [float(x) for x in m.group(1).split()]
    n = math.sqrt(sum(q * q for q in v)) or 1.0
    ut = re.sub(r"^(Uinlet\s+)\([^)]*\);",
                r"\g<1>(%.6f %.6f %.6f);" % tuple(U * q / n for q in v), ut, flags=re.M)
    u.write_text(ut)
    for fn in ("system/controlDict", "system/functionObject0"):
        q = case / fn
        if q.exists():
            s = q.read_text()
            s = re.sub(r"^(\s*magUInf\s+)[-0-9.eE+]+;", r"\g<1>%.6f;" % U, s, flags=re.M)
            s = re.sub(r"^(\s*rhoInf\s+)[-0-9.eE+]+;", r"\g<1>%.7f;" % rho, s, flags=re.M)
            q.write_text(s)
    # every U file the solver will read, including the restart fields
    #
    # THE ZERO VECTOR IS EXEMPT, AND THAT EXEMPTION IS THE WHOLE SAFETY OF THIS LOOP.
    # A blanket rewrite of every `uniform (a b c)` also hits the wall's own entry: OF-7's
    # noSlipFvPatchVectorField::write emits `value uniform (0 0 0)` into every restart field,
    # so retargeting the Mach would have set the AEROFOIL SURFACE to the freestream velocity.
    # It happens to be harmless here only because noSlip re-imposes Zero on construction and
    # ignores the dictionary value -- i.e. we would be relying on a BC to undo our own
    # corruption. Change that patch to `fixedValue` and the same line silently produces a wall
    # blowing at 248 m/s, which converges beautifully to the wrong problem (ARG-096: a gate on
    # the output cannot see a wrong input). The freestream vector is never zero and the wall
    # vector always is, so the magnitude separates them exactly.
    new = tuple(U * q / n for q in v)
    def _retarget(mo):
        if all(abs(float(x)) < 1e-30 for x in mo.group(2, 3, 4)):
            return mo.group(0)                     # a wall, not the freestream
        return "%s(%.6f %.6f %.6f)" % ((mo.group(1),) + new)
    touched = []
    for f in list(case.glob("processor*/*/U")) + \
             [q / "U" for q in case.iterdir()
              if q.is_dir() and q.name.replace(".", "").isdigit() and q.name != "0"]:
        if not f.exists():
            continue
        s = f.read_text(errors="replace")
        s2 = re.sub(r"(uniform\s+)\(\s*([-0-9.eE+]+)\s+([-0-9.eE+]+)\s+([-0-9.eE+]+)\s*\)",
                    _retarget, s)
        if s2 != s:
            f.write_text(s2)
        touched.append(f)
    # DERIVED FROM THE FILE, NOT ASSERTED BESIDE THE EDIT (GEO-080). Read back every wall
    # block and confirm it is still at rest. The check costs nothing and does not care which
    # BC type is in use, so it survives the change that would make the exemption load-bearing.
    for f in touched:
        for blk in re.findall(r"\n\s{4}\w*(?:wall|aerofoil|airfoil)\w*\s*\n\s{4}\{(.*?)\n\s{4}\}",
                              f.read_text(errors="replace"), re.S | re.I):
            for mo in re.finditer(r"uniform\s+\(\s*([-0-9.eE+]+)\s+([-0-9.eE+]+)"
                                  r"\s+([-0-9.eE+]+)\s*\)", blk):
                if any(abs(float(x)) > 1e-12 for x in mo.group(1, 2, 3)):
                    raise SystemExit("set_mach corrupted a wall BC in %s: %s"
                                     % (f, mo.group(0)))
    return U, rho, mu


def route_mach_ramp(base, ranks, steps, iters_per, out):
    """ROUTE 1: walk M up from a converged wall-resolved subcritical solution."""
    d = REPO / "cases/tutorials/wr_ramp_case09"
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(REPO / "cases/tutorials" / base, d, ignore=shutil.ignore_patterns(
        "processor*", "postProcessing", "log.*", "[1-9]*", "qoi_normalisation.json"))
    log("ROUTE 1 MACH RAMP from %s (a CONVERGED wall-resolved solution)" % base)
    log("  steps: %s" % ", ".join("%.3f" % m for m in steps))
    dp = d / "system/decomposeParDict"
    dp.write_text(re.sub(r"numberOfSubdomains\s+\d+;", "numberOfSubdomains %d;" % ranks,
                         dp.read_text()))
    sh("rm -rf processor*; decomposePar -force > log.decomposePar 2>&1 && "
       "mpirun -np %d renumberMesh -parallel -overwrite > log.renumber 2>&1" % ranks, d, 3600)
    n_proc = len(list(d.glob("processor*")))
    if n_proc != ranks:
        log("  ABORT: decomposePar made %d subdomains, asked for %d" % (n_proc, ranks))
        return []
    log("  decomposed into %d, matching -np %d" % (n_proc, ranks))
    hist, total = [], 0
    for i, M in enumerate(steps):
        U, rho, mu = set_mach(d, M, 6.5e6)
        total += iters_per
        cd_ = d / "system/controlDict"
        s = cd_.read_text()
        s = re.sub(r"^endTime.*$", "endTime         %d;" % total, s, flags=re.M)
        s = re.sub(r"^startFrom.*$",
                   "startFrom       %s;" % ("startTime" if i == 0 else "latestTime"),
                   s, flags=re.M)
        s = re.sub(r"^writeInterval.*$", "writeInterval   %d;" % max(iters_per // 3, 500),
                   s, flags=re.M)
        cd_.write_text(s)
        # AUDIT EVERY STEP, NOT ONLY THE FIRST. set_mach rewrites mu, 0/U, magUInf, rhoInf and
        # every restart field, so each step is a fresh opportunity for the exact defect that
        # cost a doubled C_d: a controlDict number that no longer matches the case around it.
        # The audit re-derives magUInf from 0/U and rhoInf from p/(RT), which is precisely the
        # consistency set_mach has to maintain, and it runs in a tempdir so it cannot disturb
        # the run. Nine steps at ~40 s is 6 minutes against a multi-hour ramp.
        g = subprocess.run([sys.executable, "scripts/audit_case_settings.py", "--quiet-pass",
                            str(d.relative_to(REPO))], cwd=str(REPO),
                           capture_output=True, text=True, timeout=1800)
        if g.returncode != 0:
            log("  AUDIT FAILED at M %.3f, not launching\n%s" % (M, g.stdout[-800:]))
            break
        log("  M %.3f : U %.3f  rho %.5f  mu %.3e  audit clean" % (M, U, rho, mu))
        sh("mpirun -np %d rhoSimpleFoam -parallel > log.M%.3f 2>&1" % (ranks, M), d, 86400)
        sh("cat log.M* > log.solver 2>/dev/null", d, 300)
        f = forces(d)
        if f is None:
            log("  M %.3f : no forces" % M); break
        hist.append(dict(M=M, **f))
        log("  M %.3f -> Cd %+.6f  Cl %+.5f  D023 %.3f  %s"
            % (M, f["cd"], f["cl"], f["d023"],
               "Cd WENT NEGATIVE at t=%.0f" % f["cd_neg_at"] if f["cd_went_negative"] else "ok"))
        (REPO / out).write_text(json.dumps({"route": "mach_ramp", "history": hist},
                                           indent=2, default=str) + "\n")
        if abs(f["cd"]) > 1 or not np.isfinite(f["cd"]):
            log("  DIVERGED at M %.3f -- the ramp got to %.3f" % (M, steps[max(i - 1, 0)]))
            break
    if hist and abs(hist[-1]["cd"]) < 1:
        s = sep_fraction(d)
        log("  reversed flow on the upper surface: %s%%  (experiment says 0 for case 9)" % s)
        hist[-1]["reversed_pct"] = s
        (REPO / out).write_text(json.dumps({"route": "mach_ramp", "history": hist},
                                           indent=2, default=str) + "\n")
    return hist



def route_transient(base, ranks, iters, out, warm_from=None, end_time=0.30,
                    damp=False, max_co=1.0):
    """ROUTE 2: transient march to steady state with rhoPimpleFoam.

    THE COLD START IS AN UNBOUNDED FIRST STEP. A steady solver advances by an effectively
    infinite step damped only by relaxation, which is why no-slip applied across a micron
    first cell blows k and omega negative on iteration one. A TRANSIENT solver has a
    timestep, and that is precisely the limiter the steady one lacks.

    STARTED SMALL AND LET TO GROW. deltaT 1e-7 with adjustTimeStep and maxCo 5: the early
    steps are tiny enough to survive the initial transient, and the step grows by itself once
    the field settles. The wolfdynamics tutorial's failure mode was the opposite -- a FIXED
    4e-08 for the whole run, which is why it burned 159,724 steps to reach 64%.

    FOUR THINGS rhoPimpleFoam NEEDS THAT rhoSimpleFoam DOES NOT, all found the hard way
    (ARG-145 amendment 1): a `rho` entry in solvers, div(phiv,p) in divSchemes, Final-corrector
    solver entries, and PIMPLE.residualControl in SINGLE-VALUE form rather than SIMPLE's
    dictionary form.
    """
    d = REPO / "cases/tutorials/wr_trans_case09"
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(REPO / "cases/tutorials" / base, d, ignore=shutil.ignore_patterns(
        "processor*", "postProcessing", "log.*", "[1-9]*", "qoi_normalisation.json"))
    log("ROUTE 2 TRANSIENT MARCH from %s%s"
        % (base, " (warm, mapped from %s)" % warm_from if warm_from else " (cold)"))
    if warm_from and not warm_map(d, warm_from):
        return []


    c = d / "system/controlDict"
    s = c.read_text()
    s = re.sub(r"^application.*$", "application     rhoPimpleFoam;", s, flags=re.M)
    s = re.sub(r"^startFrom.*$", "startFrom       startTime;", s, flags=re.M)
    # ~%.0f convective times at c/U = 4.02 ms. A COLD start needs tens of them to develop the
    # boundary layer; a WARM one inherits a developed field and only has to relax the near-wall
    # layer the new mesh added, so it is given a shorter horizon rather than the same one.
    s = re.sub(r"^endTime.*$", "endTime         %g;" % end_time, s, flags=re.M)
    s = re.sub(r"^deltaT.*$", "deltaT          1e-07;", s, flags=re.M)
    s = re.sub(r"^writeControl.*$", "writeControl    runTime;", s, flags=re.M)
    s = re.sub(r"^writeInterval.*$", "writeInterval   %g;" % (end_time / 6), s, flags=re.M)
    if "adjustTimeStep" not in s:
        # maxCo 5 WAS TOO AGGRESSIVE AND ALSO UNREACHABLE. The first attempt showed max
        # Courant oscillating 0.36 -> 46 -> 0.36 while adjustTimeStep drove deltaT to 1e-165
        # chasing it, which is what a NaN in one cell looks like from the outside. Start
        # conservative and let it grow.
        s = s.replace("writeControl", "adjustTimeStep  yes;\nmaxCo           %g;\n"
                      "maxDeltaT       1e-5;\n\nwriteControl" % max_co, 1)
    c.write_text(s)

    f = d / "system/fvSchemes"
    s = f.read_text()
    s = re.sub(r"(ddtSchemes\s*\{[^}]*?default\s+)steadyState;", r"\1Euler;", s, flags=re.S)
    if "div(phiv,p)" not in s:
        s = re.sub(r"(divSchemes\s*\n\{\s*\n)",
                   r"\1    div(phiv,p)     Gauss upwind;   // rhoPimpleFoam only\n", s, count=1)
    f.write_text(s)

    v = d / "system/fvSolution"
    s = v.read_text()
    s = re.sub(r"\bSIMPLE\b", "PIMPLE", s)
    blk = re.search(r"PIMPLE\s*\{.*?\n\}", s, re.S)
    s = s[:blk.start()] + """PIMPLE
{
    nOuterCorrectors 2;
    nCorrectors      2;
    nNonOrthogonalCorrectors 1;
    transonic        yes;
    // OF-7 wants SINGLE VALUES here; the dictionary form belongs in
    // outerCorrectorResidualControl and is a hard parse error in this block.
    residualControl
    {
        "(U|k|omega|e|h)" 1e-6;
        p                 1e-6;
    }
}""" + s[blk.end():]
    # rho and the Final correctors, which rhoSimpleFoam never asks for
    #
    # THE PRESENCE TEST IS SCOPED TO THE solvers BLOCK AND ANCHORED TO A LINE START.
    # `if fld not in s` is a SUBSTRING test over the whole file, and this file contains
    # `urf_rho 0.5;` in relaxationFactors plus two commented `rho` lines. So "rho" was found,
    # the rho solver entry was never written, and rhoPimpleFoam refused the case with
    # "keyword rho is undefined in dictionary fvSolution.solvers". A substring is not an entry
    # and the whole file is not the block -- the same shape as the presence-vs-validity failure
    # that put the one-step solver check into the audit in the first place.
    sb = re.search(r"solvers\s*\n\{\s*\n", s)
    blk_txt = s[sb.end():] if sb else s
    add = ""
    for fld, sol, pre in (("rho", "PCG", "DIC"), ("pFinal", "PBiCGStab", "DILU"),
                          ("UFinal", "PBiCGStab", "DILU"), ("hFinal", "PBiCGStab", "DILU"),
                          ("eFinal", "PBiCGStab", "DILU"), ("kFinal", "PBiCGStab", "DILU"),
                          ("omegaFinal", "PBiCGStab", "DILU"), ("rhoFinal", "PCG", "DIC")):
        if not re.search(r"^\s*\"?[\w()|]*\b%s\b[\w()|]*\"?\s*$" % re.escape(fld),
                         blk_txt, re.M):
            add += ("        %s\n        {\n            solver          %s;\n"
                    "            preconditioner  %s;\n            tolerance       1e-08;\n"
                    "            relTol          0;\n        }\n\n" % (fld, sol, pre))
    if add:
        s = re.sub(r"(solvers\s*\n\{\s*\n)", r"\1" + add, s, count=1)
    v.write_text(s)
    # DERIVED FROM THE FILE (GEO-080): every field the solver will ask for is now present as
    # an ENTRY, not merely as a character sequence somewhere in the dictionary.
    sb2 = re.search(r"solvers\s*\n\{\s*\n", s)
    for fld in ("rho", "pFinal", "UFinal"):
        if not re.search(r"^\s*\"?[\w()|]*\b%s\b[\w()|]*\"?\s*$" % re.escape(fld),
                         s[sb2.end():], re.M):
            raise SystemExit("fvSolution.solvers still has no entry for %s" % fld)

    dp = d / "system/decomposeParDict"
    dp.write_text(re.sub(r"numberOfSubdomains\s+\d+;", "numberOfSubdomains %d;" % ranks,
                         dp.read_text()))
    sh("rm -rf processor*; decomposePar -force > log.decomposePar 2>&1", d, 3600)
    n_proc = len(list(d.glob("processor*")))
    if n_proc != ranks:
        log("  ABORT: %d subdomains vs -np %d" % (n_proc, ranks)); return []
    g = subprocess.run([sys.executable, "scripts/audit_case_settings.py", "--quiet-pass",
                        "cases/tutorials/wr_trans_case09"], cwd=str(REPO),
                       capture_output=True, text=True, timeout=1800)
    if g.returncode != 0:
        log("  AUDIT FAILED, not launching\n%s" % g.stdout[-500:]); return []
    log("  audit clean; marching with adjustTimeStep, maxCo 5")
    sh("mpirun -np %d rhoPimpleFoam -parallel > log.solver 2>&1" % ranks, d, 86400)
    f2 = forces(d)
    if f2:
        f2["reversed_pct"] = sep_fraction(d)
        log("  Cd %+.6f  Cl %+.5f  reversed %s%%  (target 0.0168 / 0.803 / 0%%)"
            % (f2["cd"], f2["cl"], f2["reversed_pct"]))
        (REPO / out.replace(".json", "_transient.json")).write_text(
            json.dumps({"route": "transient", "result": f2}, indent=2, default=str) + "\n")
    return [f2] if f2 else []


def warm_map(d, warm_from):
    """mapFields a converged solution into case `d`'s 0/, and PROVE it landed.

    Shared by the steady and the transient route, because the cold start is one hypothesis and
    the time-integration scheme is another, and testing them together tells you nothing about
    either.

    THREE THINGS ARE ASSERTED, none of them from mapFields' exit status:
      1. the source HAS a reconstructed time > 0. mapFields falls back to 0/ otherwise, which
         makes a cold start wearing a warm label and nothing downstream can tell.
      2. 0/U comes out NONUNIFORM. A map that silently produced a uniform field returns 0.
      3. the WALL TREATMENT is still the low-Re set. The source is wall-MODELLED, and
         nutUSpaldingWallFunction arriving on a y+ 1 mesh would give a beautifully convergent
         run of the experiment we are not doing, with every force and figure looking right.
    """
    src = REPO / "cases/tutorials" / warm_from
    if not [q for q in src.iterdir() if q.is_dir() and q.name.replace(".", "").isdigit()
            and q.name != "0"]:
        sh("reconstructPar -latestTime > log.recon 2>&1", src, 7200)
    times = sorted((float(q.name) for q in src.iterdir()
                    if q.is_dir() and q.name.replace(".", "").isdigit() and float(q.name) > 0),
                   reverse=True)
    if not times:
        log("  ABORT: %s has no reconstructed time to map from" % warm_from); return False
    log("  mapping from %s t=%g (its converged solution)" % (warm_from, times[0]))
    sh("mapFields ../%s -sourceTime %g -consistent > log.mapFields 2>&1"
       % (warm_from, times[0]), d, 7200)
    if "nonuniform" not in (d / "0/U").read_text(errors="replace"):
        log("  ABORT: 0/U still uniform after mapFields -- the map did not take\n%s"
            % (d / "log.mapFields").read_text(errors="replace")[-600:])
        return False
    # MATCH THE PATCH OR THE GROUP: the pre-map files address the wall through its GROUP
    # (`wall`) and mapFields writes the resolved PATCH name (`aerofoil`). Keying on the group
    # alone reported MISSING for all three on a map that had preserved every BC correctly.
    want = dict(nut="nutLowReWallFunction", k="fixedValue", omega="omegaWallFunction")
    got = {}
    for fld in want:
        mm = re.search(r"\n\s{4}(?:wall|aerofoil)\s*\n\s{4}\{.*?type\s+(\w+);",
                       (d / ("0/" + fld)).read_text(errors="replace"), re.S)
        got[fld] = mm.group(1) if mm else "MISSING"
    if got != want:
        log("  ABORT: mapFields changed the wall treatment: got %s, want %s" % (got, want))
        return False
    log("  warm start real; wall BCs survived the map: %s" % got)
    return True


def apply_damping(d):
    """The two turbulence-transport changes, shared by the steady and transient routes.

    Extracted because the first version of this lived inline in one route and was inserted
    into the WRONG one by an anchor that matched twice (ARG-156 item 11a): cold+damped came
    back bit-identical to cold+undamped. A shared function cannot land in the wrong place.
    """
    # ---- DAMP THE TURBULENCE TRANSPORT, AND NOTHING ELSE.
    # Shot B (warm, undamped) reached Cd 0.016583 / Cl 0.823852 at iteration 1501 against
    # targets 0.0168 / 0.803, then drifted DOWN over 7,500 iterations -- Cl 0.824, 0.798,
    # 0.749, 0.716, 0.615 -- and diverged, with 10,983 k/omega bounding events. A slow
    # monotone loss of lift with the turbulence field being clipped is the boundary layer
    # progressively losing its eddy viscosity, not a startup transient and not the domain.
    #
    # So the two things that touch that mechanism move, and the momentum, energy and
    # pressure discretisation is left exactly as validated:
    #   div(phi,k), div(phi,omega)  linearUpwind -> upwind   (their own commented option)
    #   k, omega relaxation         0.5 -> 0.3
    # UPWIND ON THE TURBULENCE ONLY IS NOT A LOSS OF ACCURACY WHERE IT COUNTS. k and omega
    # are not reported quantities; C_d, C_l and the shock position are, and those ride on
    # div(phi,U), div(phid,p) and the energy scheme, all untouched.
    f = d / "system/fvSchemes"
    t = f.read_text()
    t = t.replace("turbulence      bounded Gauss linearUpwind limitedT;",
                  "turbulence      bounded Gauss upwind;   // damped: see ARG-156")
    f.write_text(t)
    # ANCHORED TO A LINE START, WHICH IS THE WHOLE POINT. This dict carries the commented
    # alternative `//turbulence      bounded Gauss upwind;` two lines above the live entry,
    # so a substring test finds the wanted text IN A COMMENT and reports success on a file
    # it did not change. The negative control caught exactly that: with the anchor removed
    # the edit no-opped and the check still passed. SECOND INSTANCE TONIGHT of a substring
    # standing in for an entry, after `urf_rho` swallowing the `rho` solver check -- an
    # OpenFOAM dictionary is full of commented-out alternatives, so `in text` is never the
    # right test in one.
    if not re.search(r"^\s*turbulence\s+bounded Gauss upwind;", f.read_text(), re.M):
        raise SystemExit("damped turbulence scheme did not land in fvSchemes")
    v = d / "system/fvSolution"
    t = v.read_text()
    t = re.sub(r'^(\s*"\(k\|omega\)"\s+)\$urf1;', r"\g<1>0.3;", t, flags=re.M)
    v.write_text(t)
    if not re.search(r'^\s*"\(k\|omega\)"\s+0\.3;', v.read_text(), re.M):
        raise SystemExit("damped k/omega relaxation did not land in fvSolution")
    log("  damped: turbulence divSchemes -> upwind, k/omega relaxation -> 0.3")


def route_wolf(base, ranks, iters, out, warm_from=None, damp=False):
    """ROUTES 3 and 4: the wall-resolved mesh on the VALIDATED WOLFDYNAMICS DOMAIN.

    WHY THIS SUPERSEDES THE RAMP (found 2026-08-10, ARG-156). Every one of the eight failures
    in ARG-155, and the ramp's own first step, ran on the te60 topology, whose farfield is
    -1.7 to 4.0 chords streamwise and +/-2 chords normal. The converged wall-MODELLED case
    runs +/-15. A subcritical case tolerates a 2-chord box because disturbances decay; a shock
    reflecting off it with `transonic yes` does not. THE MACH NUMBER WAS NEVER THE VARIABLE.

    So the honest experiment is the SINGLE-VARIABLE one that was never run: take the case that
    converged -- same domain, same geometry, same operating point, same schemes, same SIMPLEC,
    same 0.5 uniform relaxation -- and change ONLY the wall-normal distribution and the nut/k
    wall BCs. zCells 120 -> 200, zGrading 800 -> 1.02e5, per-cell growth 1.0578 -> 1.0597.

    TWO STARTS, because the cold start is a separate hypothesis from the domain:
      route "wolf"     COLD. If the domain was the whole problem, this alone converges.
      route "wolfmap"  WARM, mapFields from the converged wall-modelled solution. The
                       interpolated field is already right everywhere except the near-wall
                       layer the new mesh exists to resolve, so there is no cold start at all.
    Running both separates "the domain was wrong" from "the cold start was wrong". Running
    only the warm one would leave us unable to say which.
    """
    tag = ("wolfmap" if warm_from else "wolf") + ("damp" if damp else "")
    # THE CASE ID COMES FROM THE BASE, NOT A LITERAL. This was hardcoded "case09", so building
    # case07 and pointing the route at it would have written into the case09 directory and
    # reported case07 results under a case09 name -- the two would have been indistinguishable
    # afterwards, which is the expensive kind of wrong.
    mid = re.search(r"(case\d+)", base)
    if not mid:
        raise SystemExit("cannot read a case id out of base %r" % base)
    d = REPO / "cases/tutorials" / ("wr_%s_%s" % (tag, mid.group(1)))
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(REPO / "cases/tutorials" / base, d, ignore=shutil.ignore_patterns(
        "processor*", "postProcessing", "log.*", "[1-9]*", "qoi_normalisation.json"))
    log("ROUTE %s: wall-resolved on the wolfdynamics domain, from %s" % (tag.upper(), base))

    cd_ = d / "system/controlDict"
    s = cd_.read_text()
    s = re.sub(r"^endTime.*$", "endTime         %d;" % iters, s, flags=re.M)
    s = re.sub(r"^startFrom.*$", "startFrom       startTime;", s, flags=re.M)
    cd_.write_text(s)

    if warm_from and not warm_map(d, warm_from):
        return []

    if damp:
        apply_damping(d)

    dp = d / "system/decomposeParDict"
    dp.write_text(re.sub(r"numberOfSubdomains\s+\d+;", "numberOfSubdomains %d;" % ranks,
                         dp.read_text()))
    sh("rm -rf processor*; decomposePar -force > log.decomposePar 2>&1 && "
       "mpirun -np %d renumberMesh -parallel -overwrite > log.renumber 2>&1" % ranks, d, 7200)
    n_proc = len(list(d.glob("processor*")))
    if n_proc != ranks:
        log("  ABORT: decomposePar made %d subdomains, asked for %d" % (n_proc, ranks))
        return []
    g = subprocess.run([sys.executable, "scripts/audit_case_settings.py", "--quiet-pass",
                        str(d.relative_to(REPO))], cwd=str(REPO),
                       capture_output=True, text=True, timeout=1800)
    if g.returncode != 0:
        log("  AUDIT FAILED, not launching\n%s" % g.stdout[-800:]); return []
    log("  decomposed into %d, audit clean; solving %d iterations" % (n_proc, iters))
    sh("mpirun -np %d rhoSimpleFoam -parallel > log.solver 2>&1" % ranks, d, 172800)
    f = forces(d)
    if not f:
        log("  no forces -- see %s/log.solver" % d.name); return []
    f["reversed_pct"] = sep_fraction(d)
    f["route"] = tag
    tc = temperature_ceiling(d)
    if tc:
        f.update({"T0": tc["T0"], "T_max_reached": tc["T_max_reached"],
                  "T_exceeded_T0_at": tc["exceeded_at"]})
        if tc["exceeded_T0"]:
            log("  *** T_max EXCEEDED THE ADIABATIC CEILING T0 = %.1f K at it %.0f "
                "(peak %.0f K = %.1fx T0) -- UNPHYSICAL ***"
                % (tc["T0"], tc["exceeded_at"], tc["T_max_reached"], tc["T_max_over_T0"]))
        else:
            log("  T_max stayed under T0 = %.1f K (peak %.1f K)"
                % (tc["T0"], tc["T_max_reached"]))
    ps = pressure_drag_sign(d)
    if ps:
        f.update(ps)
        if ps["cd_press_went_negative"]:
            log("  *** C_d,PRESSURE WENT NEGATIVE at it %.0f -- UNPHYSICAL from there on. "
                "Last physical iteration %.0f ***"
                % (ps["cd_press_neg_at"], ps["last_physical_iteration"]))
        else:
            log("  C_d,pressure stayed positive throughout (final %.6f, viscous %.6f)"
                % (ps["cd_press_final"], ps["cd_visc_final"]))
    tgt = AGARD.get(mid.group(1), TARGET)
    log("  %s  Cd %+.6f (target %.4f)  Cl %+.5f (C_N target %.3f)  D023 %.3f  reversed %s%%"
        % (mid.group(1), f["cd"], tgt["CD"], f["cl"], tgt["CN"], f["d023"], f["reversed_pct"]))
    (REPO / out).write_text(json.dumps({"route": tag, "base": base, "warm_from": warm_from,
                                        "case": mid.group(1), "target": tgt, "result": f},
                                       indent=2, default=str) + "\n")
    return [f]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--route", choices=("ramp", "transient", "both", "wolf", "wolfmap",
                                        "transientwarm", "wolfmapdamp", "wolfdamp"),
                    default="both")
    ap.add_argument("--max-co", type=float, default=1.0,
                    help="PIMPLE target Courant; 5 was unreachable on this mesh")
    ap.add_argument("--end-time", type=float, default=0.30,
                    help="transient horizon in seconds; c/U = 4.02 ms per convective time")
    ap.add_argument("--base", default="rae2822_case09",
                    help="RAE geometry on the te60 topology -- the only\n                    wall-resolved topology that has ever converged here. The\n                    ramp starts at M 0.15, where it behaves like te60.")
    ap.add_argument("--warm-from", default="wolfrae_case09",
                    help="converged case to mapFields from, for --route wolfmap")
    ap.add_argument("--iters", type=int, default=15000)
    ap.add_argument("--ranks", type=int, default=10)
    ap.add_argument("--iters-per-step", type=int, default=4000)
    ap.add_argument("--out", default="results/wallresolved_routes.json")
    a = ap.parse_args()
    steps = [0.15, 0.30, 0.45, 0.55, 0.62, 0.67, 0.70, 0.72, 0.73]
    if a.route in ("ramp", "both"):
        route_mach_ramp(a.base, a.ranks, steps, a.iters_per_step, a.out)
    if a.route in ("transient", "both"):
        route_transient(a.base, a.ranks, a.iters_per_step, a.out,
                        end_time=a.end_time)
    if a.route == "wolf":
        route_wolf(a.base, a.ranks, a.iters, a.out)
    if a.route == "transientwarm":
        route_transient(a.base, a.ranks, a.iters_per_step, a.out,
                        warm_from=a.warm_from, end_time=a.end_time, damp=True,
                        max_co=a.max_co)
    if a.route == "wolfmap":
        route_wolf(a.base, a.ranks, a.iters, a.out, warm_from=a.warm_from)
    if a.route == "wolfmapdamp":
        route_wolf(a.base, a.ranks, a.iters, a.out, warm_from=a.warm_from, damp=True)
    if a.route == "wolfdamp":
        route_wolf(a.base, a.ranks, a.iters, a.out, damp=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
