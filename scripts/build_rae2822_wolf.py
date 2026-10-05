#!/usr/bin/env python3
"""build_rae2822_wolf.py: RAE 2822 cases 7/6/9 on the WOLFDYNAMICS transonic recipe.

WHY THIS EXISTS (project decision, 2026-08-08). Our own transonic attempt diverged on all three cases.
The wolfdynamics NACA 0012 case at M 0.699 runs steady rhoSimpleFoam to a published answer, and
REPRODUCES ON OUR OpenFOAM-org 7 UNMODIFIED: C_d +0.15%, C_l -0.18% against their log, every
one of the 20 case files byte-identical apart from the rank count. So the numerics are known
good on this exact installation, and the honest move is to transplant them WHOLE rather than
keep guessing one setting at a time.

WHAT WE ARE TAKING, and it is everything except the geometry and the operating point:
    bounded Gauss on EVERY divScheme   -- subtracts div(phi)*field, which is non-zero while a
                                          compressible steady solve converges. Our unbounded
                                          schemes are the likeliest single cause of divergence.
    div(phi,U) bounded Gauss Minmod    -- TVD limiter
    laplacian/snGrad  limited 0.5      -- not `corrected`
    SIMPLEC (consistent yes)           -- with nNonOrthogonalCorrectors 1
    relaxation 0.5 UNIFORM             -- not p 0.3 / U 0.7
    hePsiThermo + sensibleEnthalpy     -- energy in h, not e
    nutUSpaldingWallFunction           -- WALL-MODELLED, y+ 30-300
    renumberMesh before the solve      -- their fvSolution relTol 0 depends on it
    domain +/-15 chords                -- against the +/-2 of our subcritical recipe

WHAT CHANGES, and BOTH change together, which is a departure from single-variable discipline:
GEOMETRY (NACA 0012 -> RAE 2822) and OPERATING POINT (M 0.699/1.55 deg -> the AGARD cases).
Both are required to reach the deliverable and there is reference Cp and shock position to
diagnose against. IF IT MISBEHAVES, THE BISECTION STEP IS RAE GEOMETRY AT WOLF CONDITIONS:
no reference data, but it isolates the geometry.

THE CAMBER FOLD COMES WITH US. Their blockMeshDict places the aerofoil block edges with
`project ... (aerofoil)`, i.e. nearest-point projection of a straight line. On the RAE 2822
that folds the upper edge through the section beyond x/c 0.886, because camber lifts the lower
skin above z = 0 near the trailing edge and the chord line runs closer to it. Verified on our
own build: 164 negative-volume cells, all at that station. A SYMMETRIC AEROFOIL CANNOT SHOW IT,
so their case gives no warning. The edges are therefore replaced with explicit polyLines here
too, exactly as in build_rae2822_transonic.py.

y+ AT THE LOWER REYNOLDS NUMBER. Their case runs Re 1.16e7 and achieves y+ 44-274. RAE Case 9
is Re 6.5e6, a factor 0.559, so the same first cell gives roughly y+ 18-178. The bottom of that
is the buffer layer, where nutkWallFunction is invalid -- but the aerofoil carries
nutUSpaldingWallFunction, a single blended law valid across it. Recorded rather than corrected,
and MEASURED from the run rather than predicted.

WALL FUNCTIONS ARE OFF-LIMITS IN THIS PROJECT WITHOUT A DECISION ENTRY (project rules). These cases
are wall-MODELLED by construction. Nothing from them is quotable until that entry is written.
"""
import argparse
import math
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
FOAM = ("export PATH=/usr/bin:$PATH; source /opt/openfoam7/etc/bashrc >/dev/null 2>&1; "
        "export PATH=/usr/bin:$PATH; ")
WOLF = REPO / "cases/tutorials/wolf_naca0012_M07"

# AGARD f8621, decoded by scripts/parse_rae2822.py. alpha is the CORRECTED angle.
CASES = {
    "case07": dict(M=0.725, alpha=2.21, Re=6.5e6, CN=0.658, CM=-0.090, CD=0.0107),
    "case06": dict(M=0.725, alpha=2.54, Re=6.5e6, CN=0.743, CM=-0.095, CD=0.0127),
    "case09": dict(M=0.730, alpha=2.79, Re=6.5e6, CN=0.803, CM=-0.099, CD=0.0168),
    "case01": dict(M=0.676, alpha=1.93, Re=5.7e6, CN=0.566, CM=-0.082, CD=0.0085),
    "case12": dict(M=0.730, alpha=2.83, Re=2.7e6, CN=0.721, CM=-0.078, CD=0.0133),
}
CHORD, T_INF, P_INF, GAMMA = 1.0, 288.15, 1.0e5, 1.4


def wall_normal_grading(n_cells, target_first, extent=15.0):
    """blockMesh zGrading (last/first) that puts the FIRST CELL at target_first.

    blockMesh's simpleGrading ratio is LAST CELL OVER FIRST, not the per-cell growth, so it
    cannot be written down from a y+ target directly. Bisect on it instead: for n cells with
    per-cell ratio r = R^(1/(n-1)), the first cell is extent*(r-1)/(r^n - 1), monotone in R.

    THE PER-CELL GROWTH IS THE QUANTITY THAT HAS TO STAY SANE, not the total. The validated
    wolf mesh runs 120 cells at zGrading 800, i.e. growth 1.0578. 200 cells at 1.02e5 gives
    growth 1.0597 and a first cell 124x smaller: the SAME stretching function extended down to
    the wall, which is why this is a single-variable change from a converged case rather than
    a new mesh.
    """
    lo, hi = 1.0 + 1e-9, 1e14
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        r = mid ** (1.0 / (n_cells - 1))
        if extent * (r - 1) / (r ** n_cells - 1) > target_first:
            lo = mid
        else:
            hi = mid
    return mid, mid ** (1.0 / (n_cells - 1))


def measure_first_cell(case):
    """First cell height AT THE WALL, read off the built mesh (GEO-080).

    Never the number we asked for. blockMesh applies grading per EDGE and the projection moves
    the surface, so the achieved spacing is not the requested one, and the requested one is
    exactly the value a report would quote if nobody looked. Returns None rather than a guess
    if the wall patch cannot be found.
    """
    pm = pathlib.Path(case) / "constant/polyMesh"
    b = (pm / "boundary").read_text(errors="replace")
    m = re.search(r"(\w+)\s*\{[^}]*?type\s+wall;[^}]*?nFaces\s+(\d+);[^}]*?startFace\s+(\d+);",
                  b, re.S)
    if not m:
        return None
    nf, sf = int(m.group(2)), int(m.group(3))
    pts = np.array(re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)",
                              (pm / "points").read_text(errors="replace")), dtype=float)
    faces = [np.array(g.split(), dtype=int) for g in
             re.findall(r"\d+\(([\d ]+)\)", (pm / "faces").read_text(errors="replace"))]
    # THE FIRST CELL HEIGHT IS THE FACE-TO-FACE SEPARATION, NOT TWICE IT.
    # My first version multiplied by 2, on the reflex that "wall to cell CENTRE, doubled, is
    # the cell height" -- true of a centre, false of the opposite FACE, which is already a
    # full cell away. It reported 1.65e-05 c and y+ 1.8 for a mesh sitting at 8.3e-06 c and
    # y+ 0.9, i.e. it would have condemned a mesh that hit its target exactly. Hence the
    # independent second route below rather than a corrected single one (D058: a null needs an
    # answer known some other way).
    #
    # ROUTE A: wall face centroid to the nearest OTHER face centroid of the same cell. With
    # aspect ratio ~7.6e4 the side faces are three orders of magnitude further off, so the
    # nearest is unambiguously the opposite face.
    # ROUTE B: wall face VERTEX to the nearest mesh point not on the wall patch, which walks
    # one step along the wall-normal grid line. Shares no arithmetic with route A.
    own = np.array((pm / "owner").read_text(errors="replace")
                   .split("(", 1)[-1].rsplit(")", 1)[0].split(), dtype=int)
    by_cell = {}
    for j, c in enumerate(own):
        by_cell.setdefault(int(c), []).append(j)
    wall_pts = set()
    for i in range(nf):
        wall_pts.update(int(k) for k in faces[sf + i])
    free = np.array(sorted(set(range(len(pts))) - wall_pts))
    hs_a, hs_b = [], []
    step = max(nf // 200, 1)
    for i in range(0, nf, step):
        wf = faces[sf + i]
        c_wall = pts[wf].mean(axis=0)
        far = [pts[faces[j]].mean(axis=0) for j in by_cell.get(int(own[sf + i]), [])
               if j != sf + i]
        if far:
            hs_a.append(min(float(np.linalg.norm(c - c_wall)) for c in far))
        d = np.linalg.norm(pts[free] - pts[wf[0]], axis=1)
        hs_b.append(float(d.min()))
    if not hs_a:
        return None
    a_med, b_med = float(np.median(hs_a)), float(np.median(hs_b))
    if abs(a_med - b_med) / max(a_med, b_med) > 0.10:
        raise SystemExit("first-cell height disagrees between two independent measurements: "
                         "face-to-face %.3e c vs point-to-point %.3e c" % (a_med, b_med))
    return float(np.min(hs_a)), a_med, float(np.max(hs_a))


def build(tag, cfg, ranks, span=0.005, zcells=None, zgrading=None, xmcells=None,
          xucells=None, leadgrading=None, transport=None, tlimit=None, outdir=None):
    sys.path.insert(0, str(REPO / "scripts"))
    import build_rae2822_transonic as B          # reuse the section, obj and polyLine work

    suffix = "_wr" if zcells else ""
    dst = REPO / (outdir or "cases/tutorials") / ("wolfrae_%s%s" % (tag, suffix))
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(WOLF, dst, ignore=shutil.ignore_patterns(
        "processor*", "postProcessing", "log.*", "0", "[1-9]*", "polyMesh",
        "*Zone.Identifier", "sol_logs", "gnuplot"))
    shutil.copytree(WOLF / "0_org", dst / "0", dirs_exist_ok=True)

    # ---- geometry: RAE 2822, resampled, in the x-z convention their dict expects
    a = B.section()
    loop = a[:-1] if np.hypot(*(a[0] - a[-1])) < 1e-9 else a
    keep = [0] + [i for i in range(1, len(loop)) if np.hypot(*(loop[i] - loop[i - 1])) > 1e-9]
    loop = loop[keep]
    g = dst / "constant/geometry"
    for old in g.glob("*.obj"):
        old.unlink()
    # PRE-FLIP THE SECTION, because transformPoints inverts it.
    # Their pipeline runs `transformPoints -yawPitchRoll '(0 0 90)'` after blockMesh, which
    # rotates the section's normal axis and puts the aerofoil UPSIDE DOWN in the final frame.
    # For their SYMMETRIC NACA 0012 that is completely invisible. On the RAE 2822 it inverts
    # the +1.26% camber so it OPPOSES the incidence: measured camber at x/c 0.757 came out
    # -0.01249 against the source's +0.01264, and C_l read 0.097 where C_N should be 0.658.
    #
    # FIFTH CONVENTION TODAY THAT IS CORRECT FOR A SYMMETRIC SECTION AND WRONG FOR A CAMBERED
    # ONE, after the loop ordering, the projection resolution, the shoulder station and the
    # nearest-point edge projection. A symmetric aerofoil is its own mirror image, so every
    # orientation error in the chain cancels and nothing downstream can detect it.
    # MIRRORING REVERSES THE LOOP'S ORIENTATION, so the point order must reverse with it.
    # polyline_edges splits the loop at the leading edge and takes loop[:le+1] as the upper
    # surface; after a bare y-flip that slice IS THE LOWER SURFACE, so the upper block edge
    # gets lower-surface points and the edges cross through the section. Measured: 376
    # negative-volume cells and non-orthogonality 179.8, with the camber sign correct. A
    # MIRROR IS NOT A ROTATION -- it changes handedness, and every downstream consumer that
    # assumes a traversal direction has to be given one.
    loop = (loop * np.array([1.0, -1.0]))[::-1]
    npts, rs = B.write_obj(loop, g / "RAE2822.obj", span)

    # ---- flow state, R from THEIR thermophysicalProperties
    tp = dst / "constant/thermophysicalProperties"
    tt = tp.read_text()
    R = 8314.47 / float(re.search(r"^\s*molWeight\s+([-0-9.eE+]+);", tt, re.M).group(1))
    aS = math.sqrt(GAMMA * R * T_INF)
    U = cfg["M"] * aS
    rho = P_INF / (R * T_INF)
    mu = rho * U * CHORD / cfg["Re"]
    ar = math.radians(cfg["alpha"])
    tt = re.sub(r"^(\s*mu\s+)[-0-9.eE+]+;", r"\g<1>%.8e;" % mu, tt, flags=re.M)
    if transport:
        # ---- THE TRANSPORT MODEL DECIDES WHETHER `mu` MEANS ANYTHING (ARG-160).
        # sutherlandTransport reads ONLY As and Ts; the mu written two lines above is dead to
        # it, so every case built here ran at Sutherland's sea-level 1.79e-05 and Re 1.67e7
        # instead of the 6.5e6 this function computes. `const` uses mu and hits the target.
        # Constant mu over T 240-318 K is a ~10% transport error, which published RAE 2822
        # computations routinely accept; a 2.56x Reynolds error is not.
        tt = re.sub(r"transport\s+\w+;", "transport       %s;" % transport, tt)
        got = re.search(r"transport\s+(\w+);", tt).group(1)
        if got != transport:
            raise SystemExit("transport model did not land: %s" % got)
    tp.write_text(tt)

    # ---- blockMeshDict_1: swap the obj, fix the shoulder station, replace projected edges
    bm = dst / "system/blockMeshDict_1"
    t = bm.read_text()
    t = t.replace('"NACA0012.obj"', '"RAE2822.obj"')
    x_sh = 0.3
    ups = rs[rs[:, 1] >= 0]; los = rs[rs[:, 1] < 0]
    zu = float(np.interp(x_sh, np.sort(ups[:, 0]), ups[np.argsort(ups[:, 0]), 1]))
    zl = float(np.interp(x_sh, np.sort(los[:, 0]), los[np.argsort(los[:, 0]), 1]))
    for k, v in (("xUpper", x_sh), ("zUpper", zu), ("xLower", x_sh), ("zLower", zl)):
        t = re.sub(r"^(\s*%s\s+)[-0-9.eE+]+;" % k, r"\g<1>%.6f;" % v, t, flags=re.M)
    # ASSERT THE SIGN, since the defect was invisible on a symmetric section
    cam = float(np.interp(0.757, np.sort(ups[:, 0]), ups[np.argsort(ups[:, 0]), 1])
                + np.interp(0.757, np.sort(los[:, 0]), los[np.argsort(los[:, 0]), 1])) / 2
    if cam > 0:
        raise SystemExit("pre-flip failed: camber at x/c 0.757 is %+.5f, expected NEGATIVE "
                         "before transformPoints inverts it" % cam)

    ed = re.search(r"^edges\s*\n\((.*?)\n\);", t, re.S | re.M)
    kept = "\n".join(l for l in ed.group(1).splitlines() if "(aerofoil)" not in l)
    t = t[:ed.start()] + ("edges\n(\n%s\n%s);"
                          % (kept, B.polyline_edges(rs, x_sh, span))) + t[ed.end():]
    # the dict's half-width must equal the obj's half-span
    t = re.sub(r"^width1\s+[-0-9.eE+]+;", "width1 %.6f;" % (-span/2), t, flags=re.M)
    t = re.sub(r"^width2\s+[-0-9.eE+]+;", "width2  %.6f;" % (span/2), t, flags=re.M)
    # ---- WALL-RESOLVED VARIANT: the wall-normal distribution, and NOTHING ELSE.
    # The domain, the chordwise distribution, every scheme, SIMPLEC, the relaxation and the
    # solver are the validated wolfdynamics set that reproduced its published C_d to +0.15%
    # and matched the AGARD shock position to 0.001 c. Only zCells and zGrading move.
    if zcells:
        t = re.sub(r"^(\s*zCells\s+)\d+;", r"\g<1>%d;" % zcells, t, flags=re.M)
        t = re.sub(r"^(\s*zGrading\s+)[-0-9.eE+]+;", r"\g<1>%.6g;" % zgrading, t, flags=re.M)
    if xmcells:
        # ---- THE SHOCK BAND (project decision, 2026-08-10). xMCells is the ONLY count that sets the
        # chordwise spacing from the shoulder x/c 0.3 to the trailing edge, and that block is
        # `simpleGrading (1 1 zGrading)`, i.e. UNIFORM in x. So ds/c there is simply 0.7/xMCells
        # and nothing else touches it: 80 cells gives 0.00884, measured, against the 0.002-0.004
        # of published wall-resolved RAE 2822 grids. The mesh clusters at the LE (0.00257) and
        # the TE and leaves MID-CHORD, where the shock actually sits, as its coarsest region.
        #
        # UNIFORM IS THE RIGHT CHOICE HERE, not a bump graded onto x/c 0.5. The shock position
        # is what we are trying to predict and it MOVES between cases (0.507 to 0.58 measured
        # across 7/6/9), so grading toward an assumed location would put the fine band where we
        # already believe the answer is -- a mesh that presupposes its own result.
        t = re.sub(r"^(\s*xMCells\s+)\d+;", r"\g<1>%d;" % xmcells, t, flags=re.M)
        if not re.search(r"^\s*xMCells\s+%d;" % xmcells, t, flags=re.M):
            raise SystemExit("xMCells did not land in blockMeshDict_1")
    if xucells:
        # ---- THE FORWARD BLOCK, LE to the shoulder at x/c 0.30 (project decision, 2026-08-10).
        # Refining xMCells alone left a COARSE BAND AND A SPACING DISCONTINUITY exactly where
        # the shock travels when it destabilises. Measured on the xMCells-240 mesh:
        #     x/c 0.0-0.1  ds/c 0.00256       x/c 0.2-0.3  ds/c 0.00674  <- coarsest
        #     x/c 0.1-0.2  ds/c 0.00481       x/c 0.3-1.0  ds/c 0.00294
        # leadGrading 0.2 clusters this block toward the NOSE, so it is coarsest AT THE
        # SHOULDER, and the two blocks meet there with a 2.3x jump.
        #
        # THAT IS A FEEDBACK LOOP, NOT MERELY A COARSE PATCH. A shock drifting forward off the
        # refined band meets cells 2.3x coarser, rings harder, and is pushed further forward;
        # an abrupt spacing jump at a block interface is also a reflection site. It matches the
        # observed signature exactly: a SLOW ramp in the residuals led by Uy, which is the
        # normal-velocity oscillation a rung shock produces.
        t = re.sub(r"^(\s*xUCells\s+)\d+;", r"\g<1>%d;" % xucells, t, flags=re.M)
        if not re.search(r"^\s*xUCells\s+%d;" % xucells, t, flags=re.M):
            raise SystemExit("xUCells did not land in blockMeshDict_1")
    if leadgrading is not None:
        # ---- STOP CLUSTERING AT THE NOSE. leadGrading 0.2 shrinks the forward block's cells
        # 5x from the shoulder to the leading edge. That was harmless at xUCells 80 (LE
        # ds/c 0.00257) and became actively harmful at 180: LE ds/c 0.00114, MIN CELL VOLUME
        # 6.9e-12 against 4.2e-11, average non-orthogonality 7.8 -> 11.2, and a
        # divide-by-zero inside kOmegaSST::F2 on a case that ran fine one refinement earlier.
        #
        # REFINING A REGION IS NOT THE SAME AS CLUSTERING IT. What the shock needs is UNIFORM
        # ds/c ~ 0.003 wherever it might sit; what leadGrading gives is ever-finer cells at
        # the one place the shock never reaches. leadGrading 1 with xUCells 100 puts the whole
        # forward block at 0.0030, matching the aft block exactly and leaving no jump anywhere.
        t = re.sub(r"^(\s*leadGrading\s+)[-0-9.eE+]+;", r"\g<1>%.6g;" % leadgrading,
                   t, flags=re.M)
        if not re.search(r"^\s*leadGrading\s+%.6g;" % leadgrading, t, flags=re.M):
            raise SystemExit("leadGrading did not land in blockMeshDict_1")
    bm.write_text(t)

    # ---- 0/ fields and force axes. Their plane is x-y AFTER transformPoints, but the
    #      fields are written in the FINAL frame, so U and liftDir/dragDir are x-y.
    u = dst / "0/U"
    ut = u.read_text()
    ut = re.sub(r"^(Uinlet\s+)\([^)]*\);",
                r"\g<1>(%.6f %.6f 0);" % (U * math.cos(ar), U * math.sin(ar)), ut, flags=re.M)
    u.write_text(ut)
    for f, key, val in (("0/p", "pOut", P_INF), ("0/T", "Tinlet", T_INF)):
        q = dst / f
        s = q.read_text()
        s = re.sub(r"^(%s\s+)[-0-9.eE+]+;" % key, r"\g<1>%.6e;" % val, s, flags=re.M)
        q.write_text(s)

    if zcells:
        # ---- THE WALL TREATMENT, and it must move as a SET or not at all.
        #   nut   nutUSpaldingWallFunction -> nutLowReWallFunction   (nut = 0 at the wall)
        #   k     kqRWallFunction          -> fixedValue 1e-10       (k -> 0 at the wall)
        #   omega omegaWallFunction        -> UNCHANGED, deliberately
        #
        # omega LOOKS like the one to change and is the one to leave. OF-7's omegaWallFunction
        # is Menter's BLENDED form: it evaluates both 6 nu/(beta1 y^2) and the log-layer value
        # and blends, so at y+ ~ 1 it IS the viscous-sublayer condition. There is no
        # "omegaLowReWallFunction" to swap to, and replacing it with a fixedValue would be a
        # worse approximation than the thing it replaced.
        #
        # k IS NOT OPTIONAL, which is the part worth stating. kqRWallFunction is zero-gradient,
        # the companion of a wall FUNCTION: it lets the wall value float to the first cell's.
        # Against nut = 0 that is inconsistent, and an unpinned k next to a 8e-6 c first cell
        # is exactly the quantity ARG-155 recorded going negative in the first hundred
        # iterations on all eight previous attempts.
        q = dst / "0/nut"
        s = q.read_text()
        # the shipped file carries TWO `type` lines in the wall block (nutkWallFunction then
        # nutUSpaldingWallFunction) and OF takes the last; both go, or the dead one wins.
        s = re.sub(r"(wall\s*\n\s*\{)(.*?)(\n\s*\})",
                   lambda m: m.group(1) + "\n        type            nutLowReWallFunction;"
                             "\n        value           uniform 0;" + m.group(3), s, flags=re.S)
        q.write_text(s)
        q = dst / "0/k"
        s = q.read_text()
        s = re.sub(r"(wall\s*\n\s*\{)(.*?)(\n\s*\})",
                   lambda m: m.group(1) + "\n        type            fixedValue;"
                             "\n        value           uniform 1e-10;" + m.group(3),
                   s, flags=re.S)
        q.write_text(s)
        # DERIVED FROM THE FILES, NOT ASSERTED BESIDE THE EDIT (GEO-080).
        got = {f: re.search(r"wall\s*\n\s*\{.*?type\s+(\w+);", (dst / ("0/" + f)).read_text(),
                            re.S).group(1) for f in ("nut", "k", "omega")}
        want = dict(nut="nutLowReWallFunction", k="fixedValue", omega="omegaWallFunction")
        if got != want:
            raise SystemExit("wall-resolved BC swap did not land: got %s, want %s" % (got, want))

    for fn in ("system/controlDict", "system/functionObject0"):
        q = dst / fn
        if not q.exists():
            continue
        s = q.read_text()
        s = re.sub(r"^(\s*Aref\s+)[-0-9.eE+]+;", r"\g<1>%.8f;" % (float(rs[:,0].max()-rs[:,0].min())*span), s, flags=re.M)
        s = re.sub(r"^(\s*magUInf\s+)[-0-9.eE+]+;", r"\g<1>%.6f;" % U, s, flags=re.M)
        s = re.sub(r"^(\s*rhoInf\s+)[-0-9.eE+]+;", r"\g<1>%.7f;" % rho, s, flags=re.M)
        s = re.sub(r"^(\s*liftDir\s+)\([^)]*\);",
                   r"\g<1>(%.8f %.8f 0);" % (-math.sin(ar), math.cos(ar)), s, flags=re.M)
        s = re.sub(r"^(\s*dragDir\s+)\([^)]*\);",
                   r"\g<1>(%.8f %.8f 0);" % (math.cos(ar), math.sin(ar)), s, flags=re.M)
        q.write_text(s)
    if tlimit:
        # A limiter must sit at a physically defensible bound, not wherever keeps runs alive.
        # The tutorial ships max 900 against an adiabatic T0 of 318 K (ARG-159).
        fo = dst / "system/fvOptions"
        if fo.exists():
            t2 = fo.read_text()
            t2 = re.sub(r"(limitT\s*\{[^}]*?max\s+)[-0-9.eE+]+;", r"\g<1>%.6g;" % tlimit,
                        t2, flags=re.S)
            fo.write_text(t2)
            if not re.search(r"limitT\s*\{[^}]*?max\s+%.6g;" % tlimit, fo.read_text(), re.S):
                raise SystemExit("limitTemperature max did not land")

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

    d = dst / "system/decomposeParDict"
    d.write_text(re.sub(r"numberOfSubdomains\s+\d+;", "numberOfSubdomains %d;" % ranks,
                        d.read_text()))

    # ---- 0_org MUST TRACK 0/, OR IT IS A DIFFERENT AEROFOIL'S OPERATING POINT.
    # The copytree above brings the wolfdynamics tutorial's 0_org across verbatim: M 0.699,
    # p 73,048 Pa, T 283.24 K, U 236 m/s, its NACA 0012 condition. Everything after it writes
    # the RAE operating point into 0/ ONLY, so 0_org sat in every case looking exactly like a
    # pristine initial condition and holding the wrong flow. Resetting a case from it, or
    # building a new case off it, runs case 9's geometry at the tutorial's freestream and
    # reports it as case 9. Found when recipe_trial.py did precisely that and the pre-run
    # audit refused all three arms (ARG-096 item 10: only an INPUT gate can see this).
    #
    # ONE TRUTH PER CASE, the same rule that removed the stale blockMeshDict below.
    shutil.rmtree(dst / "0_org", ignore_errors=True)
    shutil.copytree(dst / "0", dst / "0_org")
    for fn in ("U", "p", "T"):
        if (dst / "0" / fn).read_text() != (dst / "0_org" / fn).read_text():
            raise SystemExit("0_org/%s does not match 0/%s" % (fn, fn))

    # ---- REMOVE THE DEAD blockMeshDict. The mesh is built from blockMeshDict_1 by explicit
    # -dict, so the plain file is never read, and it still declares the TUTORIAL's span. On
    # wolfrae_case09 that stale span is 0.2 and the mesh span is also 0.2, so the audit's
    # "blockMeshDict describes the mesh" check PASSED ON A FILE THAT HAD NOTHING TO DO WITH
    # THE MESH -- a correct answer by an invalid route, which is a defect and not a pass
    # (GEO-080). At any other span it fails and looks like a meshing error. One truth per case.
    stale = dst / "system/blockMeshDict"
    if stale.exists():
        stale.unlink()

    # ---- purgeWrite. The tutorial ships purgeWrite 0 / writeInterval 100, i.e. 150 full
    # field dumps over a 15,000-iteration run. That is the setting that filled 7 GB in a night.
    cdp = dst / "system/controlDict"
    s = cdp.read_text()
    s = re.sub(r"^(purgeWrite\s+)\d+;", r"\g<1>2;", s, flags=re.M)
    s = re.sub(r"^(writeInterval\s+)\d+;", r"\g<1>500;", s, flags=re.M)
    cdp.write_text(s)

    # ---- their exact mesh pipeline, including renumberMesh
    r = subprocess.run(["bash", "-lc", FOAM +
                        "blockMesh -dict system/blockMeshDict_1 > log.blockMesh 2>&1 && "
                        "transformPoints -yawPitchRoll '(0 0 90)' > log.transform 2>&1 && "
                        "checkMesh > log.checkMesh 2>&1"],
                       cwd=str(dst), capture_output=True, text=True, timeout=3600)
    cm = (dst / "log.checkMesh")
    if not (dst / "constant/polyMesh/points").exists() or not cm.exists():
        print("  %-16s MESH FAILED" % tag)
        print("\n".join("      " + l for l in
                        (dst / "log.blockMesh").read_text(errors="replace").splitlines()[-8:]))
        return None
    s = cm.read_text(errors="replace")
    nc = re.search(r"cells:\s+(\d+)", s)
    ar_ = re.search(r"Max aspect ratio = ([0-9.]+)", s) or re.search(r"Max aspect ratio: ([0-9.]+)", s)
    no = re.search(r"non-orthogonality Max: ([0-9.]+)", s)
    neg = re.search(r"negative volume cells:?\s*(\d+)", s)
    print("  %-19s M %.3f  a %.2f  U %.2f  mu %.3e  |  cells %s  aspect %s  nonortho %s  %s"
          % (dst.name, cfg["M"], cfg["alpha"], U, mu,
             nc.group(1) if nc else "?", ar_.group(1) if ar_ else "?",
             no.group(1) if no else "?",
             "Mesh OK" if "Mesh OK" in s else ("NEG VOL %s" % neg.group(1) if neg else "see log")))
    # THE ACHIEVED SPACING, MEASURED. The requested value is what a report quotes when nobody
    # looks; blockMesh grades per edge and the projection moves the surface, so the two differ.
    # The y+ here is a flat-plate ESTIMATE for sizing only -- the run's own yPlus function
    # object is the number that gets reported (D060: sign first, magnitude second).
    # ---- Aref AND lRef FROM THE MESHED CHORD, set AFTER blockMesh (GEO-080).
    # They were written from the RESAMPLED SECTION's extent before meshing, and blockMesh's
    # projection moves the surface: the section spans 1.0000 and the mesh 0.9994, so Aref came
    # out 0.06% wrong and every coefficient with it. Small, but it is a number describing the
    # mesh that was computed without looking at the mesh, and it blocked the figure pipeline's
    # normalisation guard -- which is the guard doing its job.
    pm = dst / "constant/polyMesh"
    bnd = (pm / "boundary").read_text(errors="replace")
    mw = re.search(r"(\w+)\s*\{[^}]*?type\s+wall;[^}]*?nFaces\s+(\d+);[^}]*?startFace\s+(\d+);",
                   bnd, re.S)
    if mw:
        nfw, sfw = int(mw.group(2)), int(mw.group(3))
        pts_ = np.array(re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)",
                                   (pm / "points").read_text(errors="replace")), dtype=float)
        fcs = [np.array(q.split(), dtype=int) for q in
               re.findall(r"\d+\(([\d ]+)\)", (pm / "faces").read_text(errors="replace"))]
        idx = np.unique(np.concatenate([fcs[sfw + i] for i in range(nfw)]))
        wp = pts_[idx]
        ext_ = wp.max(axis=0) - wp.min(axis=0)
        spdir = int(np.argmin(ext_))
        chord_mesh = float(max(e for k, e in enumerate(ext_) if k != spdir))
        lo_, hi_ = pts_.min(axis=0), pts_.max(axis=0)
        span_mesh = float((hi_ - lo_).min())
        for fn in ("system/controlDict", "system/functionObject0"):
            q = dst / fn
            if not q.exists():
                continue
            t3 = q.read_text()
            t3 = re.sub(r"^(\s*Aref\s+)[-0-9.eE+]+;", r"\g<1>%.8f;" % (chord_mesh * span_mesh),
                        t3, flags=re.M)
            t3 = re.sub(r"^(\s*lRef\s+)[-0-9.eE+]+;", r"\g<1>%.8f;" % chord_mesh,
                        t3, flags=re.M)
            q.write_text(t3)
        print("      Aref/lRef set from the MESH: chord %.6f x span %.6f = %.8f"
              % (chord_mesh, span_mesh, chord_mesh * span_mesh))

    h = measure_first_cell(dst)
    if h:
        kfac = math.sqrt(0.026 / cfg["Re"] ** (1 / 7.) / 2) * cfg["Re"] / CHORD
        print("      first cell measured  min %.3e  median %.3e  max %.3e c"
              "   -> y+_est %.2f - %.2f (flat plate, for SIZING only)"
              % (h[0], h[1], h[2], h[0] / 2 * kfac, h[2] / 2 * kfac))
    return dst


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cases", nargs="*", default=list(CASES))
    ap.add_argument("--ranks", type=int, default=6)
    ap.add_argument("--span", type=float, default=0.005,
                    help="2D span. ARBITRARY except through Aref, but it sets an\n                    UNAVOIDABLE aspect floor of span/first-cell: 0.2 m gave 44,823.")
    ap.add_argument("--wall-resolved", action="store_true",
                    help="refine the WALL-NORMAL distribution only, to y+ ~ 1, and swap the\n"
                         "                    nut/k wall BCs to the low-Re set. Builds\n"
                         "                    wolfrae_<tag>_wr beside the wall-modelled case.")
    ap.add_argument("--yplus", type=float, default=1.0, help="target y+ at the first CENTRE")
    ap.add_argument("--zcells", type=int, default=200)
    ap.add_argument("--xucells", type=int, default=None,
                    help="cells from the LE to the shoulder x/c 0.3, graded toward the nose.\n"
                         "                    80 leaves ds/c 0.00674 at the shoulder; 180 gives\n"
                         "                    ~0.0030, matching the aft block and removing the jump.")
    ap.add_argument("--transport", choices=("const","sutherland"), default=None,
                    help="const USES the computed mu and hits the target Re; sutherland\n"
                         "                    IGNORES it and pins Re at ~1.7e7 (ARG-160).")
    ap.add_argument("--tlimit", type=float, default=None,
                    help="limitTemperature max, K. Must be near the adiabatic T0, not 900.")
    ap.add_argument("--outdir", default=None, help="e.g. cases/validation")
    ap.add_argument("--leadgrading", type=float, default=None,
                    help="LE-block grading (last/first). 0.2 = 5x clustering at the nose;\n"
                         "                    1 = uniform, which is what a travelling shock wants.")
    ap.add_argument("--xmcells", type=int, default=None,
                    help="chordwise cells from x/c 0.3 to the TE, UNIFORM. ds/c = 0.7/N;\n"
                         "                    80 (default) gives 0.00884, 240 gives 0.00292.")
    a = ap.parse_args()
    zc = zg = None
    if a.wall_resolved:
        # sized off the FIRST case asked for; all three share Re 6.5e6 so it does not diverge
        Re = CASES[a.cases[0]]["Re"]
        kfac = math.sqrt(0.026 / Re ** (1 / 7.) / 2) * Re / CHORD
        zc = a.zcells
        zg, growth = wall_normal_grading(zc, 2 * a.yplus / kfac)
        print("  wall-resolved: zCells %d  zGrading %.4g  (per-cell growth %.4f; the "
              "wall-modelled mesh runs 1.0578)" % (zc, zg, growth))
    out = [build(t, CASES[t], a.ranks, a.span, zc, zg, a.xmcells, a.xucells, a.leadgrading,
                 a.transport, a.tlimit, a.outdir)
           for t in a.cases]
    ok = [o for o in out if o]
    print("\n  built %d of %d" % (len(ok), len(a.cases)))
    return 0 if len(ok) == len(a.cases) else 1


if __name__ == "__main__":
    raise SystemExit(main())
