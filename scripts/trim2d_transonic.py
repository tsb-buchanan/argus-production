#!/usr/bin/env python3
"""trim2d_transonic.py: trim an ARGUS 2D section to a target cl at a compressible condition.

WHY (project decision, 2026-08-09). scripts/trim2d.py already does exactly this -- polar, slope fit,
iterate alpha to |cl - target| <= 1e-4, write trim_history.json -- but it drives the
INCOMPRESSIBLE build_2d_case.py O-grid. The transonic recipe validated in ARG-151 is a
different case-construction path, so this is the compressible sibling. The iteration logic,
the tolerance and the history format are deliberately the same, so the run cards match.

THE CONDITION (Liming, 2026-08-09; recorded ARG-152):
    early cruise   M 0.78   CL 0.52929708745   74,500 kg at 10 km   rho 0.41271
    late  cruise   M 0.78   CL 0.54411480021   56,000 kg at 12 km   rho 0.31083
Liming asked for early cruise first: the root-bending constraint is active there and it is
where the trailing-edge/twist ordering differs most.

TWO MACH CONVENTIONS ARE RUN AS A PAIR (project decision). The wing has 27 deg of quarter-chord
sweep, so a section sees M_n = 0.78 cos 27 = 0.695 -- which lands BETWEEN RAE case 01 (0.676)
and case 07 (0.725), i.e. inside the range the compressible chain is validated over. The
freestream value 0.78 is OUTSIDE it. Running both measures what the convention is worth
instead of assuming it, and makes the extrapolation explicit rather than silent.

REYNOLDS IS DERIVED, NOT ASSUMED. Reference-data 9.6e carried only "Re ~2.5e7, order of
magnitude only". From report Table 13 plus the standard atmosphere and the aircraft-scale
mean chord S_ref/b = 122.6/38.36 = 3.196 m:
    early  V 233.6 m/s  Re 2.11e7 freestream, 1.68e7 normal to the sweep
    late   V 230.2 m/s  Re 1.61e7 freestream, 1.28e7 normal
This is Liming's STAGE 2 (flight Re for the concept comparison). His stage 1, TM X-71996 at
tunnel Re with the grit treatment stated, is a 3D baseline task and is not this script.

THE SIGN CHECK THAT MATTERS (D060). The same sections at M 0.1 trim to alpha 0.917 / 0.769 /
0.428 deg at cl 0.4777. Here the cl is HIGHER (0.529) and the Mach is HIGHER, and
compressibility RAISES the lift slope, so the trimmed alpha must come out LOWER. A higher
alpha is a sign failure, not scatter, and the run is rejected rather than reported.
"""
import argparse
import json
import math
import pathlib
import re
import os
import shutil
import subprocess
import sys
import time

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
FOAM = ("export PATH=/usr/bin:$PATH; source /opt/openfoam7/etc/bashrc >/dev/null 2>&1; "
        "export PATH=/usr/bin:$PATH; ")
WOLF = REPO / "cases/tutorials/wolf_naca0012_M07"

SECTIONS = {
    "baseline":   "geometry/derived/eet_section/eet_cruise_faired.dat",
    "mbr_eta070": "geometry/derived/eet_section/morphed/mbr_c01_l011_eta0p70.dat",
    "mbr_eta080": "geometry/derived/eet_section/morphed/mbr_c01_l011_eta0p80.dat",
    "mbr_eta090": "geometry/derived/eet_section/morphed/mbr_c01_l011_eta0p90.dat",
}
# ---- THE RESPONSE CURVE (ARG-161). Ten camber amplitudes spanning REFLEX to DROOP, run to
# build dcl(A/c) and dcd(A/c) as FUNCTIONS so any candidate's spanwise schedule maps onto them
# without new solves. The `rc` names carry the SIGNED amplitude, not a span station: rcm2117 is
# A/c -0.021174 and rcp3497 is +0.034968. `m`/`p` for minus/plus, because a name that hides the
# sign of a reflex section is a name that will be misread in a table.
#
# WHY REFLEX IS HERE AT ALL: te_b03_c01, the current selected early-cruise candidate, runs A/c
# from -0.021174 to +0.017733. HALF ITS SCHEDULE IS TRAILING-EDGE UP, and every 2D section this
# project generated before ARG-161 was positive droop, so that half had never been tested.
#
# THREE EXACT ANTISYMMETRIC PAIRS: +/-0.005318, +/-0.010633, +/-0.017733. Thin-aerofoil theory
# is exactly antisymmetric in A, so a matched pair converts "is the response linear" from a
# curve-fit judgement into a difference of two measurements with a known answer of zero.
_RC = {"rcm2117": 0.00, "rcm1773": 0.10, "rcm1063": 0.20, "rcm0532": 0.30,
       "rcp0022": 0.40, "rcp0532": 0.50, "rcp1063": 0.60, "rcp1773": 0.70,
       "rcp2598": 0.80, "rcp3497": 0.90}
for _n, _e in _RC.items():
    SECTIONS[_n] = ("geometry/derived/eet_section/morphed/resp_curve_eta%s.dat"
                    % ("%.2f" % _e).replace(".", "p"))
# rcp0022 (A/c +0.000218) IS THE NULL: geometrically the baseline to 0.02% of chord, predicted
# dcl 0.0029 at M 0.695. It is the only cell whose answer is known independently of everything
# else being tested, so it runs FIRST and the rest are not readable until it passes.
# report Table 13; rho and CL are Liming's, V and Re derived here from the standard atmosphere
STATES = {
    "early": dict(M=0.78, CL=0.52929708745, alt=10000, rho=0.41271, T=223.15),
    "late":  dict(M=0.78, CL=0.54411480021, alt=12000, rho=0.31083, T=216.65),
}
SWEEP_DEG = 27.0            # quarter-chord sweep, TP-1580 / reference-data section 2
CHORD_AC = 3.196            # aircraft-scale mean chord, S_ref/b = 122.6/38.36
LOWSPEED_ALPHA = {"baseline": 0.917, "mbr_eta070": -0.410,
                  "mbr_eta080": -1.000, "mbr_eta090": -0.818}



def free_gb(path=REPO):
    st = os.statvfs(str(path))
    return st.f_bavail * st.f_frsize / 1024 ** 3


def strip_case(case, keep_logs=True):
    """Drop everything heavy, keep the evidence.

    A 4-iteration trim leaves 3 dead cases behind, and 8 trim points would leave 24 at
    ~200 MB each. Only the CONVERGED iteration is a result; the others are the path taken,
    and their logs and dictionaries record that in ~20 MB rather than 200.
    """
    # WHAT IS BULK AND WHAT IS RESULT, measured rather than assumed:
    #     processor0/            216 MB   <- bulk, and reproducible from the reconstructed field
    #     latest time directory    8 MB   <- the U and p fields; needed for the flow plots
    #     postProcessing/         23 MB   <- forces, residuals, yPlus. THE RESULT.
    # My first version deleted postProcessing to save space. It is the CHEAPEST thing in the
    # case and the only part that cannot be regenerated once the fields are gone, which is
    # why there was no Cp or Cf to plot. Only processor* and SUPERSEDED time directories go.
    case = pathlib.Path(case)
    for q in case.glob("processor*"):
        shutil.rmtree(q, ignore_errors=True)
    times = sorted([q for q in case.iterdir()
                    if q.is_dir() and q.name not in ("0", "0_org", "constant", "system")
                    and q.name.replace(".", "").isdigit()], key=lambda q: float(q.name))
    for q in times[:-1]:                       # keep the LATEST; it carries U and p
        shutil.rmtree(q, ignore_errors=True)
    if not keep_logs:
        for q in case.glob("log.*"):
            q.unlink(missing_ok=True)


def sh(cmd, cwd=None, timeout=None):
    return subprocess.run(["bash", "-lc", FOAM + cmd], cwd=str(cwd) if cwd else None,
                          capture_output=True, text=True, timeout=timeout)


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def condition(state, convention):
    """(M, V, rho, mu, Re) for a state, on the chosen Mach convention.

    NORMAL CONVENTION: both the Mach AND the chord are taken normal to the quarter-chord, so
    Re_n = rho (V cos L) (c cos L) / mu. Taking the Mach normal but leaving the chord
    streamwise would be a frame error -- two halves of one transformation, applied to one.
    """
    s = STATES[state]
    T, rho = s["T"], s["rho"]
    # 287.058 here is the THIRD appearance of a gas constant in this file and it is the only
    # one that is allowed to differ from the case's own (ARG-161). It sets the TARGET Re, which
    # build() then holds EXACTLY by sizing mu as rho*U*c/Re with the case-consistent U, so the
    # achieved Reynolds number is cond["Re"] whatever R is used here. What it does mean is that
    # the target itself is "Re at M 0.695 with R = 287.058" rather than with the case's 287.698,
    # a 0.11% difference in the LABEL and none in the solve. V below is likewise nominal; the
    # velocity actually imposed comes from case_speed(). Do not quote V from here.
    a = math.sqrt(1.4 * 287.058 * T)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    cosL = math.cos(math.radians(SWEEP_DEG))
    if convention == "normal":
        M, V, c = s["M"] * cosL, s["M"] * a * cosL, CHORD_AC * cosL
    else:
        M, V, c = s["M"], s["M"] * a, CHORD_AC
    return dict(M=M, V=V, rho=rho, T=T, mu=mu, a=a, chord_ac=c,
                Re=rho * V * c / mu, CL=s["CL"], convention=convention, state=state)


def build(name, section, cond, alpha, ranks, span=0.2, outdir=None, recipe=None):
    """Case construction on the ARG-151 recipe, with an arbitrary section and condition."""
    sys.path.insert(0, str(REPO / "scripts"))
    import build_rae2822_transonic as B

    # PATH SAYS WHAT THE CASE IS (project decision, 2026-08-10). cases/tutorials reached 67
    # directories with five mesh generations sharing three names, which is how the one good
    # wall-resolved result was destroyed. The tree is
    # cases/argus2d/<wall treatment>/<geometry>/<convention>, so what is running is legible
    # from the path alone.
    dst = REPO / (outdir or "cases/tutorials") / name
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(WOLF, dst, ignore=shutil.ignore_patterns(
        "processor*", "postProcessing", "log.*", "0", "[1-9]*", "polyMesh",
        "*Zone.Identifier", "sol_logs", "gnuplot"))
    shutil.copytree(WOLF / "0_org", dst / "0", dirs_exist_ok=True)

    # ---- geometry. Same pre-flip as the RAE build: transformPoints inverts the section,
    # which is invisible on a symmetric aerofoil and puts a cambered one upside down.
    raw = []
    for l in (REPO / section).read_text().splitlines():
        s = l.strip()
        if s and not s.startswith("#"):
            v = s.replace(",", " ").split()
            if len(v) >= 2:
                try:
                    raw.append((float(v[0]), float(v[1])))
                except ValueError:
                    pass
    loop = np.array(raw)
    if np.hypot(*(loop[0] - loop[-1])) < 1e-9:
        loop = loop[:-1]
    keep = [0] + [i for i in range(1, len(loop)) if np.hypot(*(loop[i] - loop[i - 1])) > 1e-9]
    loop = (loop[keep] * np.array([1.0, -1.0]))[::-1]

    g = dst / "constant/geometry"
    g.mkdir(parents=True, exist_ok=True)
    for old in g.glob("*.obj"):
        old.unlink()
    npts, rs = B.write_obj(loop, g / "SECTION.obj", span)

    # ---- thermophysical state. R from the case's own molWeight, never a constant here.
    tp = dst / "constant/thermophysicalProperties"
    tt = tp.read_text()
    R = 8314.47 / float(re.search(r"^\s*molWeight\s+([-0-9.eE+]+);", tt, re.M).group(1))
    T = cond["T"]
    # our sections are chord 1 m; hold M and Re by choosing p (hence rho) and mu
    aS = math.sqrt(1.4 * R * T)
    U = cond["M"] * aS
    p_inf = cond["rho"] * R * T
    mu = cond["rho"] * U * 1.0 / cond["Re"]
    tt = re.sub(r"^(\s*mu\s+)[-0-9.eE+]+;", r"\g<1>%.8e;" % mu, tt, flags=re.M)
    # ---- transport MUST be `const`, or the mu above is dead (ARG-160).
    # The wolfdynamics case this is copied from ships `transport sutherland`, which reads ONLY
    # As and Ts. Every ARGUS case built here therefore ran at Sutherland's value for its own
    # T (1.457e-05 at 223.15 K) and Re 5.9e6 / 6.6e6, against the 1.68e7 / 2.11e7 this line
    # computes and the reference data documents: A FACTOR 2.9 TO 3.2 TOO LOW, i.e. ~26% on
    # skin friction and ~12% on total C_d. THE TARGET WAS ALWAYS RIGHT AND NEVER TOOK EFFECT.
    tt = re.sub(r"transport\s+\w+;", "transport       const;", tt)
    tp.write_text(tt)
    got = re.search(r"transport\s+(\w+);", tp.read_text()).group(1)
    if got != "const":
        raise SystemExit("transport model did not land: %r" % got)
    # ---- and the temperature limiter must sit near the adiabatic ceiling, not at 900 K.
    T0 = cond["T"] * (1 + 0.2 * cond["M"] ** 2) if "T" in cond else None
    fo = dst / "system/fvOptions"
    if T0 and fo.exists():
        t2 = re.sub(r"(limitT\s*\{[^}]*?max\s+)[-0-9.eE+]+;", r"\g<1>%.6g;" % (1.3 * T0),
                    fo.read_text(), flags=re.S)
        fo.write_text(t2)

    # ---- blockMeshDict: swap the obj, shoulder station, explicit polyLine edges
    # ---- THE ENERGY RESIDUAL WAS NEVER BEING WRITTEN.
    # system/residuals asks for `fields (p U e k omega)`, but thermophysicalProperties selects
    # `energy sensibleEnthalpy`, so the solved variable is h and `e` does not exist. The
    # function object skips it SILENTLY: every residual plot this project has made for a
    # compressible case is missing the one equation that actually fails. The temperature
    # runaway had to be found through fieldMinMax instead, iterations later than necessary.
    rq = dst / "system/residuals"
    if rq.exists():
        t0 = rq.read_text()
        t0 = re.sub(r"fields\s*\([^)]*\);", "fields (p U h k omega);", t0)
        rq.write_text(t0)
        if "fields (p U h k omega);" not in rq.read_text():
            raise SystemExit("residuals field list did not land")

    bm = dst / "system/blockMeshDict_1"
    # ---- WALL-NORMAL SPACING FOR THE CORRECTED REYNOLDS NUMBER (ARG-160 amendment 1).
    # The inherited wolf distribution (zCells 120, zGrading 800, first cell 1.022e-03 c) was
    # sized for Re 6.5e6. At the ARGUS flight Re of 2.12e7 it predicts y+ 369, ABOVE the 300
    # ceiling for a wall function: the first cell then sits outside the log layer and the law
    # of the wall does not apply there. y+ SCALES AS Re^0.9, so correcting the Reynolds number
    # moved the wall treatment out of validity WITHOUT ANY MESH FILE CHANGING -- which is why
    # this is set from the condition rather than left as a constant.
    # 160 cells at zGrading 2618 gives first cell 2.77e-04 c and y+ ~100, the same regime the
    # RAE validation cases achieve (mean 91-102), so the wall model does the same job in the
    # deliverable as in the case that validates it.
    _t = bm.read_text()
    # ---- CHORDWISE, WHERE THE SHOCK ACTUALLY SITS.
    # Measured on the converged ARGUS baseline: the shock lands at x/c 0.224 and the mesh is
    # COARSEST there (ds/c 0.0071), because the forward block is graded toward the nose and the
    # shoulder is at x/c 0.30. The shock spans 3.0 cells against a published 0.002-0.004 c
    # target, and a shock hunting between cells in a coarse band is a numerical limit cycle --
    # which is exactly what the 1,154-iteration oscillation in C_d is.
    # xUCells 120 with leadGrading 0.5 was the best-conditioned option in the RAE sweep
    # (min cell volume 5.17e-11 against the 4.23e-11 baseline, i.e. BETTER), and it avoids the
    # tiny nose cells that caused a kOmegaSST::F2 divide-by-zero at xUCells 180 / lead 0.2.
    _t = re.sub(r"^(\s*xUCells\s+)\d+;", r"\g<1>%d;" % XUCELLS, _t, flags=re.M)
    _t = re.sub(r"^(\s*xMCells\s+)\d+;", r"\g<1>%d;" % XMCELLS, _t, flags=re.M)
    _t = re.sub(r"^(\s*leadGrading\s+)[-0-9.eE+]+;", r"\g<1>%g;" % LEADGRADING, _t, flags=re.M)
    _t = re.sub(r"^(\s*zCells\s+)\d+;", r"\g<1>%d;" % ZCELLS, _t, flags=re.M)
    _t = re.sub(r"^(\s*zGrading\s+)[-0-9.eE+]+;", r"\g<1>%g;" % ZGRADING, _t,
                flags=re.M)
    bm.write_text(_t)
    # ---- nNonOrthogonalCorrectors. One corrector is thin for a mesh at 38-45 deg max
    # non-orthogonality, and the correction is what keeps the pressure equation consistent on
    # a skewed cell. Cheap: about 20% more work per outer iteration.
    fvs = dst / "system/fvSolution"
    _u = re.sub(r"(nNonOrthogonalCorrectors\s+)\d+;", r"\g<1>%d;" % NNONORTH, fvs.read_text())
    _u = re.sub(r"^(urf1\s+)[-0-9.eE+]+;", r"\g<1>%g;" % URF, _u, flags=re.M)
    fvs.write_text(_u)
    if not re.search(r"^urf1\s+%g;" % URF, fvs.read_text(), re.M):
        raise SystemExit("relaxation urf1 did not land")
    if not re.search(r"nNonOrthogonalCorrectors\s+%d;" % NNONORTH, fvs.read_text()):
        raise SystemExit("nNonOrthogonalCorrectors did not land")
    if not (re.search(r"^\s*zCells\s+%d;" % ZCELLS, bm.read_text(), re.M)
            and re.search(r"^\s*zGrading\s+%g;" % ZGRADING, bm.read_text(), re.M)):
        raise SystemExit("wall-normal spacing did not land in blockMeshDict_1")
    t = bm.read_text().replace('"NACA0012.obj"', '"SECTION.obj"')
    x_sh = 0.3
    ups, los = rs[rs[:, 1] >= 0], rs[rs[:, 1] < 0]
    zu = float(np.interp(x_sh, np.sort(ups[:, 0]), ups[np.argsort(ups[:, 0]), 1]))
    zl = float(np.interp(x_sh, np.sort(los[:, 0]), los[np.argsort(los[:, 0]), 1]))
    # THE LEADING AND TRAILING EDGE SEEDS MUST COME FROM THE SECTION TOO.
    # I set only the shoulder pair and left xLead/zLead/xTrail/zTrail at the NACA's
    # (0,0) and (1,0). That is harmless for a section whose TE sits on the chord line and
    # wrong for one that does not: the morphed eta 0.70 TE is at (1.0015, +0.01479), so the
    # TE VERTEX was being projected from a point 0.015 c away from the surface it is meant
    # to land on -- 41 negative-volume cells, all in one radial column at x/c 0.997 on the
    # drooped side. The baseline survived only because its TE is at +0.00125, ten times
    # closer. ANOTHER SEED THAT IS RIGHT FOR THE REFERENCE AEROFOIL AND WRONG FOR OURS.
    te_i = int(np.argmax(rs[:, 0]))
    le_i = int(np.argmin(rs[:, 0]))
    for k, v in (("xUpper", x_sh), ("zUpper", zu), ("xLower", x_sh), ("zLower", zl),
                 ("xLead", float(rs[le_i, 0])), ("zLead", float(rs[le_i, 1])),
                 ("xTrail", float(rs[te_i, 0])), ("zTrail", float(rs[te_i, 1]))):
        t = re.sub(r"^(\s*%s\s+)[-0-9.eE+]+;" % k, r"\g<1>%.6f;" % v, t, flags=re.M)
    t = re.sub(r"^width1\s+[-0-9.eE+]+;", "width1 %.6f;" % (-span / 2), t, flags=re.M)
    t = re.sub(r"^width2\s+[-0-9.eE+]+;", "width2  %.6f;" % (span / 2), t, flags=re.M)
    ed = re.search(r"^edges\s*\n\((.*?)\n\);", t, re.S | re.M)
    kept = "\n".join(l for l in ed.group(1).splitlines() if "(aerofoil)" not in l)
    t = t[:ed.start()] + ("edges\n(\n%s\n%s);"
                          % (kept, B.polyline_edges(rs, x_sh, span))) + t[ed.end():]
    bm.write_text(t)

    # ---- 0/ fields. Incidence enters through the freestream direction, in the x-y plane
    # that transformPoints leaves behind.
    ar = math.radians(alpha)
    u = dst / "0/U"
    u.write_text(re.sub(r"^(Uinlet\s+)\([^)]*\);",
                        r"\g<1>(%.6f %.6f 0);" % (U * math.cos(ar), U * math.sin(ar)),
                        u.read_text(), flags=re.M))
    for f, key, val in (("0/p", "pOut", p_inf), ("0/T", "Tinlet", T)):
        q = dst / f
        q.write_text(re.sub(r"^(%s\s+)[-0-9.eE+]+;" % key, r"\g<1>%.6e;" % val,
                            q.read_text(), flags=re.M))

    for fn in ("system/controlDict", "system/functionObject0"):
        q = dst / fn
        if not q.exists():
            continue
        s = q.read_text()
        # BOTH reference lengths come from the MEASURED section. The faired EET baseline
        # spans 1.0014 after derotation, not exactly 1, and leaving lRef at the nominal 1.0
        # puts the moment coefficient on a different chord from the forces.
        chord_m = float(rs[:, 0].max() - rs[:, 0].min())
        s = re.sub(r"^(\s*Aref\s+)[-0-9.eE+]+;", r"\g<1>%.8f;" % (chord_m * span), s, flags=re.M)
        s = re.sub(r"^(\s*lRef\s+)[-0-9.eE+]+;", r"\g<1>%.8f;" % chord_m, s, flags=re.M)
        s = re.sub(r"^(\s*magUInf\s+)[-0-9.eE+]+;", r"\g<1>%.6f;" % U, s, flags=re.M)
        s = re.sub(r"^(\s*rhoInf\s+)[-0-9.eE+]+;", r"\g<1>%.7f;" % cond["rho"], s, flags=re.M)
        s = re.sub(r"^(\s*liftDir\s+)\([^)]*\);",
                   r"\g<1>(%.8f %.8f 0);" % (-math.sin(ar), math.cos(ar)), s, flags=re.M)
        s = re.sub(r"^(\s*dragDir\s+)\([^)]*\);",
                   r"\g<1>(%.8f %.8f 0);" % (math.cos(ar), math.sin(ar)), s, flags=re.M)
        q.write_text(s)
    d = dst / "system/decomposeParDict"
    d.write_text(re.sub(r"numberOfSubdomains\s+\d+;", "numberOfSubdomains %d;" % ranks,
                        d.read_text()))
    cd = dst / "system/controlDict"
    s = cd.read_text()
    s = re.sub(r"^purgeWrite.*$", "purgeWrite      2;", s, flags=re.M)
    s = re.sub(r"^writeInterval.*$", "writeInterval   5000;", s, flags=re.M)
    cd.write_text(s)

    r = sh("blockMesh -dict system/blockMeshDict_1 > log.blockMesh 2>&1 && "
           "transformPoints -yawPitchRoll '(0 0 90)' > log.transform 2>&1 && "
           "checkMesh > log.checkMesh 2>&1", dst, 3600)

    # REFERENCE LENGTHS COME FROM THE MESH THAT WAS BUILT, not the section that was asked for.
    # Setting them from the resampled ordinates gave lRef 1.0015 against a meshed chord of
    # 0.9972 -- a 0.43% mismatch that failed the audit on every morphed section. The resample
    # and the projection produce slightly different chords, and comparing across them is the
    # frame error D068 exists for. Derive AFTER, from the artefact (GEO-080).
    if (dst / "constant/polyMesh/points").exists():
        sys.path.insert(0, str(REPO / "scripts"))
        import case_normalisation as _cn
        try:
            _n = _cn.derive(str(dst))
            c_mesh, s_mesh = _n["geometry"]["chord"], _n["geometry"]["span"]
            for fn in ("system/controlDict", "system/functionObject0"):
                q = dst / fn
                if not q.exists():
                    continue
                s2 = q.read_text()
                s2 = re.sub(r"^(\s*Aref\s+)[-0-9.eE+]+;", r"\g<1>%.8f;" % (c_mesh * s_mesh),
                            s2, flags=re.M)
                s2 = re.sub(r"^(\s*lRef\s+)[-0-9.eE+]+;", r"\g<1>%.8f;" % c_mesh,
                            s2, flags=re.M)
                q.write_text(s2)
            (dst / _cn.NAME).unlink(missing_ok=True)   # stale the moment we edited the dict
        except SystemExit:
            pass
    if not (dst / "constant/polyMesh/points").exists():
        return None, "blockMesh failed"
    cm = (dst / "log.checkMesh").read_text(errors="replace")
    neg = re.search(r"negative volume cells:?\s*(\d+)", cm)
    if neg and int(neg.group(1)) > 0:
        return None, "%s negative-volume cells" % neg.group(1)

    # ---- APPLY THE SETTLED RECIPE, LAST, AND VERIFY IT FROM THE CASE (ARG-161).
    # Everything above patches a TEMPLATE, which carries whatever numerics the template was
    # built with. That was fine while this script owned its own recipe; it is not fine now that
    # a recipe has been SELECTED on evidence (`flin2`: energy convection Gauss filteredLinear2
    # 0.2 0). A case that silently ran the template's numerics while the report said `flin2`
    # would be the run-card literal failure in a new place -- a claim authored beside the
    # action rather than derived from the thing it describes.
    #
    # solver_recipe.apply() WRITES THEN VERIFIES with foamDictionary -expand and never returns
    # having only written, so this cannot half-apply. The SA conversion is separate because it
    # touches constant/turbulenceProperties and the 0/ fields, not the schemes.
    # ORDER MATTERS AND I GOT IT BACKWARDS FIRST TIME. to_spalart_allmaras() rewrites the
    # turbulence properties and the 0/ fields, INCLUDING the solver and relaxation field names
    # (k|omega -> nuTilda). apply() then writes the recipe's own fvSchemes/fvSolution whole and
    # verifies them. Applying first and converting second leaves the converter's edits on top
    # of the verified dictionaries, so verify() fails against the spec it just wrote -- which
    # is what happened: both null-gate cases returned "recipe sa_flin2 did not verify".
    # scripts/recipe_trial.py:140-141 has the correct order and I did not read it first.
    if recipe:
        import solver_recipe as _sr
        if recipe not in _sr.RECIPES:
            return None, "recipe %s is not registered in solver_recipe.RECIPES" % recipe
        if recipe.startswith("sa_"):
            _sr.to_spalart_allmaras(dst, verbose=False)
        _sr.apply(dst, recipe, verbose=False)      # writes AND verifies; raises otherwise
        (dst / "RECIPE").write_text(recipe + "\n")
        # ---- THE RECIPE'S OWN RELAXATION STANDS. DO NOT OVERRIDE IT WITH URF.
        #
        # solver_recipe.apply() writes system/fvSolution WHOLE, so the recipe's
        # `U 0.9; nuTilda 0.9;` replaces the URF = 0.3 set 130 lines above. THAT IS
        # CORRECT BEHAVIOUR AND MUST NOT BE "FIXED". I did fix it, and it was wrong.
        #
        # URF = 0.3 was measured by ARG-160 against the SHIPPED WOLFDYNAMICS TEMPLATE's
        # numerics. `flin2` is a different scheme set, selected on evidence in
        # LOSS_AND_RECIPE_VERDICT.md, and it was validated AS A PACKAGE WITH ITS OWN
        # RELAXATION. Forcing 0.3 onto it is D077's shape exactly: a value established in
        # one configuration reused in another where the thing that made it true no longer
        # holds.
        #
        # MEASURED, and the margin is not subtle. polar() cases run with the recipe's own
        # 0.9 settle to D023 0.0012 to 0.158 counts with bnd_tail 0. The same baseline run
        # through trim() with URF forced to 0.3 gave D023 95.0 counts at alpha -0.1758 and
        # 46.4 at -0.7049, with cl swinging 6,000 to 8,000 counts peak to peak. THREE
        # ORDERS OF MAGNITUDE WORSE, in the direction I had claimed was an improvement.
        #
        # If a recipe ever needs different relaxation, register a recipe variant in
        # solver_recipe.RECIPES so the change is named, versioned and testable. Do not
        # patch it in here, where it silently applies to every recipe at once.
    # ---- COLLAPSE DUPLICATED `type` ENTRIES IN 0/ (ARG-161). RUNS LAST, DELIBERATELY:
    # to_spalart_allmaras() rewrites the 0/ fields, so an earlier pass was undone by it and
    # the audit still caught the duplicate at build time. Anything that touches 0/ must
    # come BEFORE this, and 0_org/0.orig are normalised too because run() restores from them.
    # The shipped wolfdynamics template declares BOTH nutkWallFunction and
    # nutUSpaldingWallFunction in 0/nut's `wall` block. OpenFOAM takes the LAST and says
    # nothing, so the runs were always Spalding and no result is wrong. But the file reads as
    # the other one, and a regex looking for `type\s+(\w+);` finds the DEAD entry first --
    # including this project's own verifier. KEEP THE LAST, which is what the solver used, so
    # behaviour is bit-identical and the file stops lying about it. The template itself is a
    # DELIVERED artefact and is not edited; this normalises the COPY.
    for _d in ("0", "0_org", "0.orig"):
        if not (dst / _d).is_dir():
            continue
        for _f in sorted((dst / _d).iterdir()):
            if not _f.is_file():
                continue
            _t = _f.read_text(errors="replace")


            def _collapse(mo):
                body = mo.group(2)
                ty = re.findall(r"^[ \t]*type[ \t]+[\w.]+[ \t]*;[ \t]*\n", body, re.M)
                if len(ty) < 2:
                    return mo.group(0)
                body = re.sub(r"^[ \t]*type[ \t]+[\w.]+[ \t]*;[ \t]*\n", "", body, count=len(ty) - 1,
                              flags=re.M)
                return mo.group(1) + body + mo.group(3)
            _n = re.sub(r"(\w+\s*\n?\s*\{)([^{}]*)(\})", _collapse, _t)
            if _n != _t:
                _f.write_text(_n)

    return dst, "ok"



def finalise(case, tag):
    """Reconstruct, sample the surface, and write every plot the report needs.

    WHY THIS IS A STEP AND NOT AN AFTERTHOUGHT (requirement, 2026-08-09): Cp and Cf data
    to plot. The wolfdynamics case ships forces, residuals and yPlus function objects but NOT
    a surface sampler, so a converged run leaves no Cp or Cf on disk at all. Sampling has to
    happen while the fields still exist -- and it must happen BEFORE strip_case removes the
    processor directories, because the surface sample needs the reconstructed field.

    Produces, per case:
        Cp(x/c) and Cf(x/c)      surface distributions
        y+(x/c)                  the wall regime, measured
        U and p contours         the final field, to see the shock and any separation
        forces history           cl and cd against iteration, with the D023 window marked
    """
    case = pathlib.Path(case)
    sh("reconstructPar -latestTime > log.recon 2>&1", case, 3600)
    sh("rhoSimpleFoam -postProcess -func wallShearStress -latestTime > log.wss 2>&1", case, 3600)
    (case / "system/surfaces").write_text(
        'type surfaces;\nlibs ("libsampling.so");\nsurfaceFormat raw;\n'
        'interpolationScheme cell;\nfields (p wallShearStress yPlus);\n'
        'surfaces ( airfoil { type patch; patches (".*(wall|aerofoil|airfoil).*"); '
        'interpolate false; } );\n')
    sh("postProcess -func surfaces -latestTime > log.sample 2>&1", case, 3600)
    got = sorted(case.glob("postProcessing/surfaces/*/*.raw"))
    r = subprocess.run([sys.executable, "scripts/plot_case_report.py",
                        "--case", str(case.relative_to(REPO)), "--tag", tag],
                       cwd=str(REPO), capture_output=True, text=True, timeout=1800)
    for l in (r.stdout or "").splitlines():
        if l.strip():
            log("      " + l.strip())
    if r.returncode != 0:
        log("      figures FAILED: %s" % (r.stderr or "")[-300:])
    return [q.name for q in got]



def rotate_freestream(case, alpha, U):
    """Point the freestream at a new alpha WITHOUT rebuilding, and warm-continue.

    WHY (project decision, 2026-08-09). Rebuilding the case for every alpha throws away a converged
    solution that is one alpha step away and re-runs the risky cold start each time. That cold
    start is what killed iteration 2 of the normal point: 2,937 bounding events on k, the field
    growing 357 -> 765, then a floating-point exception at iteration 1,911 -- while iteration 4
    ran clean at essentially the SAME alpha (0.288 vs 0.294 deg) on a BYTE-IDENTICAL mesh. Same
    mesh, same condition, opposite outcome: a marginal startup, not a physical instability.

    scripts/trim2d.py already warm-continues for exactly this reason and I failed to carry it
    over. Starting from the previous solution removes the cold start entirely and reaches a
    given D023 in far less wall-clock, which is also what unblocks the 1e-4 trim tolerance.

    Rotates: the freestream vector in every U file (0/, latest time, and each processor), and
    the forceCoeffs lift/drag directions, which must follow the flow or C_d picks up a
    component of lift (D019).
    """
    case = pathlib.Path(case)
    ar = math.radians(alpha)
    vec = "(%.6f %.6f 0)" % (U * math.cos(ar), U * math.sin(ar))
    n = 0
    for f in list(case.glob("0/U")) + list(case.glob("processor*/*/U")) \
            + [q / "U" for q in case.iterdir()
               if q.is_dir() and q.name.replace(".", "").isdigit() and q.name != "0"]:
        if not f.exists():
            continue
        s = f.read_text(errors="replace")
        s2 = re.sub(r"^(Uinlet\s+)\([^)]*\);", r"\g<1>%s;" % vec, s, flags=re.M)
        s2 = re.sub(r"(uniform\s+)\(\s*[-0-9.eE+]+\s+[-0-9.eE+]+\s+[-0-9.eE+]+\s*\)",
                    r"\g<1>%s" % vec, s2)
        if s2 != s:
            f.write_text(s2)
            n += 1
    for fn in ("system/controlDict", "system/functionObject0"):
        q = case / fn
        if not q.exists():
            continue
        s = q.read_text()
        s = re.sub(r"^(\s*liftDir\s+)\([^)]*\);",
                   r"\g<1>(%.8f %.8f 0);" % (-math.sin(ar), math.cos(ar)), s, flags=re.M)
        s = re.sub(r"^(\s*dragDir\s+)\([^)]*\);",
                   r"\g<1>(%.8f %.8f 0);" % (math.cos(ar), math.sin(ar)), s, flags=re.M)
        s = re.sub(r"^startFrom.*$", "startFrom       latestTime;", s, flags=re.M)
        q.write_text(s)
    return n


def read_cl_cd(case):
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
    return dict(cl=float(cl[-1]), cd=float(cd[-1]), d023=float(d023(cd.tolist())),
                iters=int(a[-1, 0]), bnd_tail=int(tail))


# CHORDWISE. leadGrading STAYS AT THE TUTORIAL'S 0.2: raising it to 0.5 to even out the
# forward block took max non-orthogonality from 38.3 to 69.1 and max skewness from 0.49 to
# 1.29, and the run diverged with h exploding. THE GRADING RATIO, NOT THE CELL COUNT, IS WHAT
# WRECKS THE ORTHOGONALITY -- so add cells and leave the ratio alone.
# MEDIUM GRID LEVEL, which is the ORIGINAL chordwise resolution with the wall spacing
# corrected for the Reynolds number. Chordwise refinement to 120/240 was tried and REJECTED:
# it took max non-orthogonality from 45.5 to 69.1 and skewness from 1.14 to 1.29, and the run
# diverged with h exploding. The medium level has IDENTICAL quality to the original mesh
# (45.5 / 1.137), so nothing was given up by fixing y+.
XUCELLS, XMCELLS, LEADGRADING = 80, 80, 0.2
NNONORTH = 2      # the wolf recipe ships 1; these meshes run 38-45 max non-orthogonality
URF = 0.3         # the wolf recipe ships 0.5. At the CORRECT Re 1.68e7 that diverges: the
                  # physical viscous damping is 2.85x weaker than at the Sutherland Re these
                  # cases used to run at, so the same factor is far more aggressive. 0.3 was
                  # measured to recover a case that 0.5 could not (ARG-160 controls D vs E).
# ---- SIZED ON THE MEASURED y+ MAXIMUM, NOT THE PREDICTED MEAN (project decision, 2026-08-11).
# The old pair gave first cell 2.757e-04 c and a MEASURED y+ of min 1.69 / mean 80.4 /
# MAX 639.6 at the correct normal-convention Re 1.681e7, i.e. 2.1x outside
# nutUSpaldingWallFunction's band. At the FREESTREAM condition, Re 2.1145e7, y+ scales as
# Re^0.9 and the maximum would be 786.
#
# THE COMMENT IT REPLACES SAID "y+ ~80-100" AND WAS TRUE OF THE MEAN. The flat-plate estimate
# it came from returned 80.41 against a measured mean of 80.38, reproducing the mean to 0.04%
# while being blind to the maximum, because a flat-plate correlation has no leading edge in
# it. A band is a constraint on EVERY wall face, so the maximum is the quantity that has to
# meet it; the audit now reads it from the run rather than predicting it.
#
# SIZED FOR THE HARDER CONDITION so one mesh serves both conventions: first cell
# 1.052e-04 c puts y+max at 300 at freestream and ~245 at normal. Going to 200 cells rather
# than holding 160 also IMPROVES the stretching, per-cell growth 1.0507 -> 1.0448, for a 25%
# cell increase (76,800 -> 96,000). A finer wall spacing on the same cell count would have
# made the mesh both finer AND more stretched, which is how you trade one defect for another.
# REVERTED, and the reason is a lesson not a typo. I refined this to first cell 1.052e-04 c
# to bring a MEASURED y+ maximum of 639.6 down to 300 -- but that 639.6 was the last write
# before that case DIED with a floating point exception at iteration 800, i.e. a measurement
# taken off a diverging field. The runs that actually converge on the refined mesh measure
# y+ mean 30 / max 46, so the ORIGINAL spacing was giving mean ~79 / max ~121 all along:
# squarely in the wall-function band, and exactly what the flat-plate estimate predicted.
#
# THE REFINEMENT MADE THINGS WORSE IN THREE WAYS. y+ ~30 is the worst place to sit with a wall
# function -- too fine for the log law to be valid, too coarse to resolve the sublayer -- and
# it tripled the near-wall aspect ratio on a mesh checkMesh already WARNS about (7,578 max,
# 1,680 cells flagged, against 779 on the RAE mesh that every recipe survives).
#
# "MEASURED BEATS PREDICTED" ONLY IF THE MEASUREMENT CAME FROM A VALID SOLUTION. The audit now
# refuses to use y+ from a run that crashed, which is the gate this needed and did not have.
# 200 cells rather than the original 160 is kept: same first cell, gentler growth (1.0343 vs
# 1.0507), which is free.
ZCELLS, ZGRADING = 200, 2036.0    # first cell 2.757e-04 c -> converged y+ mean ~79, max ~121
UPWIND_ITERS = 800      # ARG-098's startup stage; never reportable
D023_TRIM_MAX = 5.0     # counts; a trim iterate must be settled to be believed
BND_TAIL_MAX = 20       # k/omega clipping in the final fifth = not converged


def ramp_reynolds(case, ranks, mu_target, stages=5, iters=600):
    """Walk mu DOWN to its target, warm-continuing, so Re rises in steps.

    WHY (project decision, 2026-08-10): reuse the routine that produced the converged cases in
    the first place. These cases DID converge -- at Re 6.6e6, which is what
    Sutherland was silently imposing. At the correct 2.12e7 they diverge in ~240 iterations of
    second-order, and an isolation test rules the mesh out: the old wall-normal distribution
    (y+ 369) and the new one (y+ 100) die at iteration 231 and 241 with 436 and 423 bounding
    events. IDENTICAL FAILURE ON BOTH MESHES AT THE SAME Re. The Reynolds number is the
    variable, so it is the one to ramp.

    THIS IS THE MACH-RAMP IDEA ON THE OTHER AXIS, and unlike that one it starts from a
    condition this exact case is KNOWN to reach: Re 6.6e6 is where all fifteen ARGUS cases
    converged. Each stage warm-continues, so the solver never meets a boundary layer it has
    not been walked into.

    mu ONLY. Nothing else moves: same mesh, same schemes, same relaxation, same freestream
    velocity and direction. Re rises purely because mu falls.
    """
    tp = case / "constant/thermophysicalProperties"
    mu0 = float(re.search(r"^\s*mu\s+([-0-9.eE+]+);", tp.read_text(), re.M).group(1))
    # geometric walk from the Re that works to the Re we want
    mus = [mu0 * (mu_target / mu0) ** ((i + 1) / stages) for i in range(stages)]
    cd = case / "system/controlDict"
    total = int(re.search(r"^endTime\s+(\d+);", cd.read_text(), re.M).group(1))
    for i, mu in enumerate(mus):
        t = re.sub(r"^(\s*mu\s+)[-0-9.eE+]+;", r"\g<1>%.8e;" % mu, tp.read_text(), flags=re.M)
        tp.write_text(t)
        got = float(re.search(r"^\s*mu\s+([-0-9.eE+]+);", tp.read_text(), re.M).group(1))
        if abs(got - mu) / mu > 1e-6:
            raise SystemExit("mu did not land: %g vs %g" % (got, mu))
        total += iters
        s2 = re.sub(r"^endTime.*$", "endTime         %d;" % total, cd.read_text(), flags=re.M)
        s2 = re.sub(r"^startFrom.*$", "startFrom       latestTime;", s2, flags=re.M)
        s2 = re.sub(r"^writeInterval.*$", "writeInterval   %d;" % iters, s2, flags=re.M)
        cd.write_text(s2)
        sh("mpirun -np %d rhoSimpleFoam -parallel > log.reramp%d 2>&1" % (ranks, i), case, 86400)
        r = read_cl_cd(case)
        log("      Re stage %d/%d: mu %.4e -> cl %s  cd %s  D023 %s"
            % (i + 1, stages, mu,
               "%.5f" % r["cl"] if r else "-", "%.6f" % r["cd"] if r else "-",
               "%.1f" % r["d023"] if r else "-"))
        if r is None or not np.isfinite(r["cl"]) or abs(r["cl"]) > 5.0:
            return False, total
    return True, total


def run(case, ranks, iters, warm=False, end_at=None):
    """Solve. On a warm continue the decomposition is REUSED, not rebuilt."""
    if not warm:
        sh("rm -rf processor*; decomposePar -force > log.decomposePar 2>&1 && "
           "mpirun -np %d renumberMesh -parallel -overwrite > log.renumber 2>&1" % ranks,
           case, 3600)
    cd = case / "system/controlDict"
    s = cd.read_text()
    s = re.sub(r"^endTime.*$", "endTime         %d;" % (end_at or iters), s, flags=re.M)
    if not warm:
        s = re.sub(r"^startFrom.*$", "startFrom       startTime;", s, flags=re.M)
    cd.write_text(s)
    # ---- NO FIRST-ORDER UPWIND STARTUP ON THIS RECIPE. IT DESTROYS THE SOLUTION.
    # One was added citing ARG-098, and it broke every ARGUS case. Measured, same
    # case, same mesh, same Sutherland Re, only the startup differing:
    #     WITH the upwind stage : T collapses to the 100 K / 900 K limiters by iteration 451
    #     WITHOUT it            : T settles at 205.9 / 244.5 K, matching the originals' 202/245
    # T0 here is 244.7 K, so the second is physically exact and the first is nonsense.
    #
    # WHY: setting `energy`, `energy1` and hence div(phi,h) and div(phi,K) to first-order
    # upwind while `div(phid,p)` STAYS second-order `Gauss limitedLinear 1` leaves the pressure
    # work inconsistent with the kinetic-energy transport in the same equation. The mismatch is
    # not a loss of accuracy, it is an unbalanced energy budget, and the temperature runs away.
    #
    # ARG-098 PRESCRIBES THE UPWIND STARTUP FOR THE INCOMPRESSIBLE 3D simpleFoam CASES, where
    # there is no `div(phid,p)` and no energy equation to unbalance. Carrying it across to the
    # compressible transonic recipe is D077's shape: a remedy valid in one configuration reused
    # where the thing that made it valid no longer holds. The wolfdynamics recipe cold-starts
    # into second order perfectly well and always did.
    sh("mpirun -np %d rhoSimpleFoam -parallel > log.solver.%s 2>&1"
       % (ranks, "warm" if warm else "cold"), case, 86400)
    sh("cat log.solver.* > log.solver 2>/dev/null", case, 300)
    return read_cl_cd(case)


# TRIMMED ALPHAS FROM THE ORIGINAL RUN (results/overnight_argus_2d.json), used as the
# STARTING POINT for the corrected-Reynolds rerun. Those runs were at the wrong Re, so these
# are not answers -- but they are far better initial guesses than the low-speed alphas, and a
# trim loop's starting point does not have to be right, only close enough not to diverge.
# Re 6.6e6 -> 2.12e7 will move the trim point by a fraction of a degree; +/-1.8 deg, which is
# what the low-speed alpha costs on the baseline at M 0.78, is what kills it.
KNOWN_TRIM_ALPHA = {
    ("mbr_eta070", "normal"): -1.4105, ("mbr_eta080", "normal"): -1.9716,
    ("mbr_eta090", "normal"): -1.4493,
    ("mbr_eta070", "freestream"): -0.8625, ("mbr_eta080", "freestream"): -0.8810,
    ("mbr_eta090", "freestream"): -0.8543,
    # BASELINE, recovered from the MACH SWEEP's own trimmed freestream direction. It was
    # absent from the trim results only because the sweep used a different case-naming
    # convention, so it fell through to the Prandtl-Glauert estimate -- and that estimate put
    # it at +0.034 deg where the case DIVERGES, while -0.176 deg runs 30,000 iterations clean.
    # A 0.21 deg difference decided it. The morphed sections never hit this because they had
    # entries here all along.
    ("baseline", "freestream"): -0.5275,   # from argus2d_baseline_early_M0p780
    ("baseline", "normal"):     -0.1758,   # from argus2d_baseline_early_M0p700 (M 0.700 vs 0.695)
}


def transonic_alpha0(alpha_lowspeed, cl_lowspeed, cond, slope_lowspeed=0.11, section=None):
    """A COMPRESSIBLE starting alpha, not the low-speed one.

    WHY (found 2026-08-10 by reading the original run's own history). The three morphed
    sections converged from their low-speed alphas ONLY because those happened to sit within
    ~0.5 deg of the transonic trim point:
        eta070  low-speed -0.410 -> trimmed -0.863
        eta080  low-speed -1.000 -> trimmed -0.881
        eta090  low-speed -0.818 -> trimmed -0.854
    THE BASELINE'S LOW-SPEED ALPHA IS +0.917, roughly 1.8 deg above where it must trim, and at
    M 0.78 that is a strong shock and separation on iteration one. It diverged every time. The
    baseline was never in the original trim results -- only the three morphed sections -- so
    this starting point had never actually been exercised.

    PRANDTL-GLAUERT ON THE LIFT SLOPE, which is the cheapest defensible estimate: the
    incompressible slope is divided by beta = sqrt(1 - M^2), so the alpha needed to hold a
    given cl shrinks by the same factor. Crude near M_crit and it does not need to be better
    than the trim loop's own secant, which converges from anywhere sane.
    """
    key = (section, cond.get("convention"))
    if key in KNOWN_TRIM_ALPHA:
        a0 = KNOWN_TRIM_ALPHA[key]
        log("    alpha0 %+.4f deg (the original run's trimmed alpha for %s/%s)"
            % (a0, section, cond.get("convention")))
        return a0
    beta = math.sqrt(max(1.0 - cond["M"] ** 2, 0.05))
    slope_c = slope_lowspeed / beta
    cl_here = cl_lowspeed / beta            # what this alpha would give compressibly
    a0 = alpha_lowspeed - (cl_here - cond["CL"]) / slope_c
    # Keep it inside the range these sections are known to behave in.
    a0 = max(-4.0, min(3.0, a0))
    log("    alpha0 %+.4f deg (Prandtl-Glauert from low-speed %+.4f at cl %.4f, beta %.4f)"
        % (a0, alpha_lowspeed, cl_lowspeed, beta))
    return a0


def case_speed(case, cond):
    """Freestream speed for `cond`, using THE CASE'S OWN gas constant.

    THERE WERE TWO GAS CONSTANTS IN THIS FILE AND THEY DISAGREED (ARG-161). build() derives
    R from the case's `molWeight` (28.9, so R = 287.698) because that is what the SOLVER will
    integrate with, and writes 0/U, mu and forceCoeffs.magUInf from it. polar() and trim()
    recomputed U from a hardcoded 287.058 and then rotated the freestream with THAT, silently
    overwriting the consistent value.

    THE COST: U 208.124 against 208.356, 0.11%. The requested Mach 0.694985 was actually
    achieved as 0.6942, and since mu was sized as rho*U*c/Re from the OTHER U, the Reynolds
    number was off by the same 0.11%. Both small; neither is a thing to leave in a delta study
    that is trying to resolve sub-count drag differences, and both were invisible because each
    routine was internally consistent.

    CAUGHT BY THE PRE-RUN AUDIT, not by inspection: `magUInf 208.3561 vs |0/U| 208.1243`. That
    check exists because a reference velocity and the velocity actually imposed are two
    different objects that nothing else reconciles. It refused to run the case, which is the
    behaviour that makes it worth having.
    """
    tp = pathlib.Path(case) / "constant/thermophysicalProperties"
    mw = float(re.search(r"^\s*molWeight\s+([-0-9.eE+]+);", tp.read_text(), re.M).group(1))
    return cond["M"] * math.sqrt(1.4 * (8314.47 / mw) * cond["T"])


def polar(tag, section, cond, alphas, ranks, iters, outdir=None, recipe=None):
    """A FIXED alpha sweep, warm-continued, instead of a secant chase to one cl.

    WHY (project decision, 2026-08-10): separate cases at several angles of attack, to give
    dcl/dalpha and cl/cd plots. Better than the trim loop here, for four reasons:

      1. THE TRIM LOOP EXTRAPOLATES, AND EXTRAPOLATION IS WHAT BLEW UP. Given one bad iterate
         it computes the next alpha from a secant through it. The baseline went +0.034 ->
         cl 3e4 -> nonsense. A FIXED alpha list cannot walk anywhere it was not sent.
      2. ONE DIVERGENT ALPHA COSTS ONE POINT, NOT THE CASE. The trim loop returns nothing if
         its sequence fails; a sweep still has the other four alphas and can interpolate.
      3. THE POLAR IS WHAT THE REPORT NEEDS ANYWAY -- dcl/dalpha, cd(cl), and the drag at ANY
         target lift, not just the one we happened to trim to.
      4. IT IS COMPARABLE TO XFOIL, which produces polars, not trim points.

    Same cost as a trim: the alphas warm-continue through one case exactly as the trim
    iterates did. The trim point is then INTERPOLATED from cl(alpha) rather than chased.
    """
    name = "argus2d_%s" % tag
    case, why = build(name, section, cond, alphas[0], ranks, outdir=outdir,
                      recipe=recipe)
    if case is None:
        return dict(reason="build: " + why, points=[])
    rel = str(case.relative_to(REPO))
    subprocess.run([sys.executable, "scripts/case_normalisation.py", "--case", rel],
                   cwd=str(REPO), capture_output=True, text=True, timeout=600)

    def gate(when):
        r = subprocess.run([sys.executable, "scripts/audit_case_settings.py", "--quiet-pass",
                            rel], cwd=str(REPO), capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            log("    AUDIT FAILED %s\n%s" % (when, r.stdout[-400:]))
            return False
        return True

    if not gate("at build"):
        return dict(reason="audit", points=[])
    U = case_speed(case, cond)   # the CASE's gas constant, not a second one (ARG-161)
    pts, total = [], 0
    for i, al in enumerate(alphas):
        if i > 0:
            rotate_freestream(case, al, U)
            if not gate("after rotating to %+.4f" % al):
                break
        total += iters
        r = run(case, ranks, iters, warm=(i > 0), end_at=total)
        if r is None:
            log("    alpha %+.4f : no forces" % al); continue
        bad = (not np.isfinite(r["cl"]) or abs(r["cl"]) > 5.0
               or r["d023"] > D023_TRIM_MAX or r["bnd_tail"] > BND_TAIL_MAX)
        log("    alpha %+.4f -> cl %.5f  cd %.6f  D023 %.2f  bnd %d %s"
            % (al, r["cl"], r["cd"], r["d023"], r["bnd_tail"], "REJECTED" if bad else "ok"))
        # RECORD THE REJECTED POINT TOO, flagged. A polar that silently drops its failures
        # looks smoother than the data is (GEO-089: every member of the set is accounted for).
        pts.append(dict(alpha_deg=al, usable=not bad, **r))
        if bad:
            log("      (kept in the record as unusable; the sweep continues)")
    good = [q for q in pts if q["usable"]]
    out = dict(points=pts, n_usable=len(good), n_total=len(pts),
               case=str(case.relative_to(REPO)))
    if len(good) >= 2:
        a = np.array([q["alpha_deg"] for q in good]); cl = np.array([q["cl"] for q in good])
        cd = np.array([q["cd"] for q in good])
        k = np.argsort(a)
        out["dcl_dalpha"] = float(np.polyfit(a[k], cl[k], 1)[0])
        tgt = cond["CL"]
        if cl.min() <= tgt <= cl.max():
            out["alpha_at_target"] = float(np.interp(tgt, cl[k], a[k]))
            out["cd_at_target"] = float(np.interp(tgt, cl[k], cd[k]))
            out["target_bracketed"] = True
            log("    TARGET cl %.6f bracketed: alpha %+.4f deg, cd %.6f"
                % (tgt, out["alpha_at_target"], out["cd_at_target"]))
        else:
            out["target_bracketed"] = False
            log("    *** target cl %.6f NOT BRACKETED by [%.5f, %.5f] -- extend the sweep ***"
                % (tgt, cl.min(), cl.max()))
    return out


def trim(tag, section, cond, alpha0, ranks, iters, tol, max_iters, slope,
         outdir=None, recipe=None):
    """One case, warm-continued through the alpha sequence.

    THE COLD START WAS THE FAILURE MODE. Rebuilding per alpha re-ran a marginal startup four
    times per point and lost one of them outright. Here the case is built ONCE and each new
    alpha rotates the freestream and continues from the previous converged field, which is
    both safer and faster -- the solution starts a few thousandths of a cl away instead of
    from uniform flow.
    """
    name = "argus2d_%s" % tag
    # RECIPE THREADED THROUGH, AND THIS WAS THE ARG-186 ROOT CAUSE. `--recipe` was
    # accepted by the CLI, honoured by polar(), and DROPPED HERE: trim() had no recipe
    # parameter and called build() without one, so EVERY TRIM EVER RUN used the shipped
    # wolfdynamics template's numerics instead of the recipe the recipe study selected
    # on evidence. The template ships pMinFactor 0.5, a static-pressure floor that
    # clamps at Cp -1.479 on a section whose suction peak needs Cp -2.0, and it does not
    # carry the flin2 energy-convection scheme at all. Ten trims returned "diverged on
    # the first iterate" against numerics nobody chose.
    case, why = build(name, section, cond, alpha0, ranks, outdir=outdir, recipe=recipe)
    if case is None:
        log("    build failed: %s" % why)
        return dict(converged=False, reason="build: " + why, history=[])
    rel = str(case.relative_to(REPO))
    subprocess.run([sys.executable, "scripts/case_normalisation.py", "--case", rel],
                   cwd=str(REPO), capture_output=True, text=True, timeout=600)

    def gate(when):
        """EVERY SETTING RE-CHECKED BEFORE EVERY SOLVE (project decision, 2026-08-10).

        The audit used to run ONCE, at build. But the trim loop rotates the freestream between
        solves -- it rewrites 0/U, liftDir and dragDir on every alpha -- so the case that gets
        audited is not the case that gets solved after iteration 0. A gate that inspects the
        input once and then lets the input change is a gate on a different object (ARG-096).

        All 17 checks, including the one-step solver run, and it HALTS rather than warns. The
        cost is ~40 s against an 8,000-iteration solve.
        """
        r = subprocess.run([sys.executable, "scripts/audit_case_settings.py", "--quiet-pass",
                            rel], cwd=str(REPO), capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            log("    AUDIT FAILED %s -- not solving\n%s" % (when, r.stdout[-500:]))
            return False
        log("    settings audit clean (%s)" % when)
        return True

    if not gate("at build"):
        return dict(converged=False, reason="audit at build", history=[])

    hist, alpha, total = [], alpha0, 0
    U = case_speed(case, cond)   # the CASE's gas constant, not a second one (ARG-161)
    for it in range(max_iters):
        if free_gb() < 25:
            log("    ABORT: %.1f GB free" % free_gb())
            break
        warm = it > 0
        if warm:
            nrot = rotate_freestream(case, alpha, U)
            log("    alpha -> %+.4f deg (warm, %d U files rotated)" % (alpha, nrot))
            if not gate("after rotating to alpha %+.4f" % alpha):
                return dict(converged=False, reason="audit after rotation", history=hist)
        total += iters
        r = run(case, ranks, iters, warm=warm, end_at=total)
        if r is None:
            return dict(converged=False, reason="no forces", history=hist)
        err = r["cl"] - cond["CL"]
        # AN UNCONVERGED ITERATE IS NOT A MEASUREMENT, AND THIS LOOP USED TO TREAT IT AS ONE.
        # The magnitude test alone passed a solve that reported D023 1230 counts and 560
        # bounding events in its final fifth: cl came back 1.75 from a field with k at 1e22.
        # The driver then extrapolated on it, reached alpha -8.5 and then -29.4 deg, and spent
        # the rest of the run solving a massively stalled aerofoil at M 0.78. BOTH NUMBERS
        # WERE PRINTED ON THE SAME LINE AS THE cl AND NEITHER WAS TESTED.
        #
        # The secant needs a cl that MEANS something. D023 and the tail bounding count are
        # exactly the evidence for that, and they are already in hand.
        bad = []
        if not np.isfinite(r["cl"]) or abs(r["cl"]) > 5.0:
            bad.append("cl %.3e" % r["cl"])
        if r["d023"] > D023_TRIM_MAX:
            bad.append("D023 %.1f > %.1f counts" % (r["d023"], D023_TRIM_MAX))
        if r["bnd_tail"] > BND_TAIL_MAX:
            bad.append("%d bounding events in the final fifth" % r["bnd_tail"])
        if bad:
            log("    REJECTED at alpha %+.4f (%s) -- halving the step"
                % (alpha, "; ".join(bad)))
            if hist:
                alpha = 0.5 * (alpha + hist[-1]["alpha_deg"])
                continue
            return dict(converged=False, reason="diverged on the first iterate", history=hist)
        hist.append(dict(alpha_deg=alpha, **r, err=err))
        log("    alpha %+.4f -> cl %.5f (target %.5f, err %+.2e)  cd %.6f  D023 %.3f  bnd_tail %d"
            % (alpha, r["cl"], cond["CL"], err, r["cd"], r["d023"], r["bnd_tail"]))
        if abs(err) <= tol:
            log("    CONVERGED; sampling surface and writing figures")
            finalise(case, tag)
            return dict(converged=True, alpha_trim_deg=alpha, target_CL=cond["CL"],
                        tolerance=tol, case=str(case.relative_to(REPO)), history=hist, **r)
        if len(hist) >= 2 and abs(hist[-1]["alpha_deg"] - hist[-2]["alpha_deg"]) > 1e-9:
            s_new = ((hist[-1]["cl"] - hist[-2]["cl"])
                     / (hist[-1]["alpha_deg"] - hist[-2]["alpha_deg"]))
            # A DEGENERATE SECANT IS A SYMPTOM, NOT A SLOPE. Two alphas returning the same cl
            # means the measurement is wrong, not that the aerofoil has no lift curve. Keep
            # the nominal slope and say so, rather than dividing by ~0 and killing the run.
            if abs(s_new) < 1e-3:
                log("    secant slope %.2e is degenerate; keeping the nominal %.3f "
                    "(two alphas returned the same cl -- check the force reader)"
                    % (s_new, slope))
            else:
                slope = s_new
        alpha = alpha - err / slope

    # Not converged. THE CASE ON DISK HOLDS THE LAST ITERATE, so finalise it and say which
    # one it is -- my previous version finalised a different directory from the best result
    # and produced a figure labelled "besteffort" showing an earlier alpha entirely.
    if hist:
        b = min(hist, key=lambda h: abs(h["err"]))
        last = hist[-1]
        log("    not converged in %d; best err %+.2e at alpha %+.4f, LAST was %+.2e at %+.4f"
            % (max_iters, b["err"], b["alpha_deg"], last["err"], last["alpha_deg"]))
        log("    finalising the case on disk, which is the LAST iterate (alpha %+.4f)"
            % last["alpha_deg"])
        finalise(case, tag)
        return dict(converged=False, reason="max iters", best_effort=True,
                    alpha_trim_deg=last["alpha_deg"], best_alpha_deg=b["alpha_deg"],
                    best_err=b["err"], target_CL=cond["CL"], tolerance=tol,
                    case=str(case.relative_to(REPO)), history=hist,
                    **{k: v for k, v in last.items() if k != "alpha_deg"})
    return dict(converged=False, reason="max iters", history=hist)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sections", nargs="*", default=list(SECTIONS))
    ap.add_argument("--state", default="early", choices=list(STATES))
    ap.add_argument("--conventions", nargs="*", default=["normal", "freestream"])
    ap.add_argument("--ranks", type=int, default=6)
    ap.add_argument("--iters", type=int, default=8000)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--max-iters", type=int, default=4)
    ap.add_argument("--slope", type=float, default=0.13, help="dcl/dalpha per deg, compressible")
    ap.add_argument("--out", default="results/overnight_argus_2d.json")
    ap.add_argument("--polar", nargs="*", type=float, default=None,
                    help="ALPHA SWEEP instead of a trim: a fixed list of alphas, warm-\n"
                         "                    continued, with the target cl INTERPOLATED. Robust where the\n"
                         "                    secant chase is not, and it yields cl(alpha) and cd(cl).")
    ap.add_argument("--recipe", default=None,
                    help="solver_recipe name applied and VERIFIED after the mesh, "
                         "e.g. sa_flin2. Omitted = the template's own numerics.")
    ap.add_argument("--outdir", default=None,
                    help="case tree root; default cases/argus2d/<wt>/<geom>/<conv>")
    ap.add_argument("--mach-sweep", nargs="*", type=float, default=None,
                    help="MACH SWEEP at fixed cl. Overrides --conventions: each M is run as "
                         "its own condition with Re scaled from the same atmosphere, so the "
                         "result is c_d(M) at constant lift -- the drag-rise curve. "
                         "TM X-71996 measures drag divergence at M 0.802 on this wing family.")
    a = ap.parse_args()

    out = []
    if a.mach_sweep:
        # A MACH SWEEP IS A SEQUENCE OF CONDITIONS, not a sequence of conventions. Re is
        # recomputed at each M from the SAME atmosphere and the same aircraft chord, because
        # in flight Re follows V and V follows M -- holding Re fixed while sweeping M would
        # be a different experiment and not the one drag rise is defined by.
        base = STATES[a.state]
        for M in a.mach_sweep:
            saved = base["M"]
            base["M"] = M
            c = condition(a.state, "freestream")
            base["M"] = saved
            c["convention"] = "M%.3f" % M
            log("condition %s M %.4f : V %.1f m/s  Re %.3e  target cl %.8f"
                % (a.state, c["M"], c["V"], c["Re"], c["CL"]))
            for sec in a.sections:
                tag = "%s_%s_M%s" % (sec, a.state, ("%.3f" % M).replace(".", "p"))
                r = trim(tag, SECTIONS[sec], c,
                         transonic_alpha0(LOWSPEED_ALPHA.get(sec, 0.0), 0.4777, c, section=sec),
                         a.ranks, a.iters, a.tol, a.max_iters, a.slope, recipe=a.recipe)
                out.append(dict(tag=tag, section=sec, **c, **r))
                (REPO / a.out).write_text(json.dumps(out, indent=2, default=str) + "\n")
        log("")
        log("DRAG RISE at fixed cl %.6f" % STATES[a.state]["CL"])
        log("   M       Re         alpha      cd         d(cd)/dM   conv")
        prev = None
        for r in out:
            if "cd" not in r:
                log("   %.3f   %-10s (%s)" % (r["M"], "", r.get("reason","?"))); continue
            slope = "" if prev is None else "%9.4f" % ((r["cd"]-prev[1])/(r["M"]-prev[0]))
            log("   %.3f   %.3e  %+7.4f  %.6f  %s  %s"
                % (r["M"], r["Re"], r["alpha_trim_deg"], r["cd"], slope,
                   "yes" if r.get("converged") else "best"))
            prev = (r["M"], r["cd"])
        log("")
        log("   M_dd by the standard criterion is where d(cd)/dM = 0.1")
        log("   TM X-71996 measures drag divergence at M 0.802 on this wing family")
        log("wrote %s" % a.out)
        return 0
    if a.polar:
        for conv in a.conventions:
            c = condition(a.state, conv)
            log("condition %s/%s : M %.4f  Re %.3e  target cl %.8f"
                % (a.state, conv, c["M"], c["Re"], c["CL"]))
            # ---- ONE ALPHA LIST FOR EVERY SECTION IN THIS CONDITION (ARG-161).
            # A DELTA STUDY COMPARES SECTIONS AT THE SAME ALPHA OR IT COMPARES NOTHING. Seeding
            # each section from its own LOWSPEED_ALPHA gives each a DIFFERENT sweep: the first
            # null-gate attempt put baseline at -0.176/+0.224 and rcp0022 at -0.883/-0.483, so
            # the "dcl between two near-identical sections" would have been 0.7 deg of incidence
            # difference and nothing to do with camber. The analytic prediction in ARG-161 is a
            # FIXED-ALPHA statement, so the frame has to be fixed alpha too.
            #
            # The list is anchored on the BASELINE's seed, which is the one with a measured
            # trimmed alpha behind it, and every section is run on it unchanged. Sections whose
            # own trim sits far from the baseline's simply produce a different cl on the same
            # alpha, which is the SIGNAL, not a problem to correct away.
            a0 = transonic_alpha0(LOWSPEED_ALPHA.get("baseline", 0.0), 0.4777, c,
                                  section="baseline")
            als = sorted(a0 + d for d in a.polar)
            log("  COMMON alpha list for every section: %s"
                % ", ".join("%+.3f" % q for q in als))
            for sec in a.sections:
                tag = "%s_%s_%s" % (sec, a.state, conv)
                log("  %s : alphas %s" % (tag, ", ".join("%+.3f" % q for q in als)))
                # ---- THE STATE AND THE RECIPE ARE IN THE PATH (ARG-161). Without the state,
                # early and late cruise write to the SAME directory for the same section and
                # convention, and the second run silently overwrites the first. They are not
                # the same case: late is a different CL and, more to the point, a 24% lower
                # Reynolds number. Without the recipe, a re-run under different numerics
                # overwrites a result it is not comparable to. Both have precedent here --
                # cases/tutorials reached 67 directories with five mesh generations sharing
                # three names, and that is how the one good wall-resolved result was lost.
                r = polar(tag, SECTIONS[sec], c, als, a.ranks, a.iters, recipe=a.recipe,
                          outdir="cases/argus2d/%s/%s/%s/%s"
                                 % (a.recipe or "template", a.state, conv, sec))
                out.append(dict(tag=tag, section=sec,
                                **{k: v for k, v in c.items()}, **r))
                (REPO / a.out).write_text(json.dumps(out, indent=2, default=str) + "\n")
        log("wrote %s" % a.out)
        return 0

    for conv in a.conventions:
        c = condition(a.state, conv)
        log("condition %s/%s : M %.4f  V %.1f m/s  Re %.3e  target CL %.8f"
            % (a.state, conv, c["M"], c["V"], c["Re"], c["CL"]))
        for sec in a.sections:
            tag = "%s_%s_%s" % (sec, a.state, conv)
            log("  %s (low-speed alpha was %+.3f deg at cl 0.4777)"
                % (tag, LOWSPEED_ALPHA.get(sec, float("nan"))))
            a0 = transonic_alpha0(LOWSPEED_ALPHA.get(sec, 0.0), 0.4777, c, section=sec)
            r = trim(tag, SECTIONS[sec], c, a0, a.ranks, a.iters, a.tol, a.max_iters,
                     a.slope, outdir=('cases/argus2d/%s/%s/%s/%s'
                                      % (a.recipe or 'wall_modelled', a.state, sec, conv)),
                     recipe=a.recipe)
            # SIGN CHECK (D060): higher cl AND higher Mach must give a LOWER trimmed alpha
            if r.get("converged"):
                ls = LOWSPEED_ALPHA.get(sec)
                if ls is not None and r["alpha_trim_deg"] > ls:
                    r["sign_check"] = ("FAIL: trimmed alpha %+.3f exceeds the low-speed %+.3f "
                                       "at a higher cl and Mach; compressibility raises the "
                                       "lift slope so it must be lower"
                                       % (r["alpha_trim_deg"], ls))
                    log("    *** %s" % r["sign_check"])
                else:
                    r["sign_check"] = "ok"
            out.append(dict(tag=tag, section=sec, **{k: v for k, v in c.items()}, **r))
            (REPO / a.out).write_text(json.dumps(out, indent=2, default=str) + "\n")

    log("")
    log("SUMMARY  state=%s  (wall-modelled, ARG-151; M 0.78 adoption pending an entry)" % a.state)
    log("  section        conv         alpha      cl        cd         D023   conv?  sign")
    for r in out:
        if r.get("converged"):
            log("  %-14s %-12s %+7.4f  %.5f  %.6f  %5.3f  yes    %s"
                % (r["section"], r["convention"], r["alpha_trim_deg"], r["cl"], r["cd"],
                   r["d023"], "ok" if r.get("sign_check") == "ok" else "FAIL"))
        else:
            log("  %-14s %-12s   --       --        --         --     no     (%s)"
                % (r["section"], r["convention"], r.get("reason", "?")))
    log("")
    log("wrote %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
