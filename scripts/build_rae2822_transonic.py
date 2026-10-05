#!/usr/bin/env python3
"""build_rae2822_transonic.py: RAE 2822 cases 7, 6 and 9 on the frozen 2D recipe.

WHY THESE THREE (project decision, 2026-08-08). One case gives a point check; a pair gives a DERIVATIVE,
and a derivative cancels everything common to both sides -- the tunnel correction, the trip,
the geometry, our own mesh bias.

    7 -> 6   M 0.725 and Re 6.5e6 BOTH FIXED, alpha 2.21 -> 2.54.
             THE ONLY CLEAN alpha PAIR IN THE SET.  dC_D/dalpha = 60.6 counts/deg.
    6 -> 9   alpha nearly fixed, M 0.725 -> 0.730: a DRAG-RISE step, which is the
             sensitivity M 0.78 will actually demand of us.
    9        the strongest shock at Re 6.5e6 and the community benchmark.

READ 6 -> 9 AS AN alpha DERIVATIVE AND YOU GET 164 ct/deg, 2.7x the truth. The excess is
dM = 0.005 sitting on the drag-rise knee. That confound is the whole reason case 7 is here.

WHAT CHANGES FROM THE FROZEN NACA 0012 RECIPE (recipes/2d_compressible_aerofoil.json), and
each of these is a deliberate departure rather than an inherited default:

 1. GEOMETRY. New aerofoil means a new .obj, not a new generator. The RAE 2822 is CAMBERED and
    12.1% thick, so the blockMeshDict's projection seeds (xUpper/zUpper, xLower/zLower) cannot
    stay at the NACA's symmetric +/-0.06 -- they are taken from the section's own extrema.
 2. INCIDENCE. alpha enters through the FREESTREAM DIRECTION, not by rotating the mesh: the
    C-grid's leading-edge clustering is built around the geometric nose, and rotating the body
    moves the stagnation point off it. liftDir and dragDir in forceCoeffs rotate WITH the
    flow, or C_D picks up a component of lift (D019: axes lagging U by 0.0083 deg biased cd by
    cl x lag).
 3. SHOCK CAPTURING. The NACA case is subcritical and its unlimited linear divergence schemes
    are fine there. A shock needs LIMITED schemes or it rings; `limitedLinear 1` on the
    convective terms is the conservative choice, and vanLeer on the density.
 4. CHORDWISE CLUSTERING MOVES. The frozen recipe clusters toward the trailing edge because
    that is where the NACA case's error lived. Here the graded quantity is SHOCK POSITION at
    x/c 0.537-0.563, so the aft block is clustered toward MID-CHORD as well. The trailing-edge
    resolution is kept, because it is what made the Cp match.
 5. alpha IS THE CORRECTED ANGLE, NOT THE GEOMETRIC ONE, and the difference is ~0.4 deg. We
    run free air, so the tunnel correction is already applied to the reference. Published
    studies differ on this and it is a large part of why they scatter; ours is pinned here.

THE alpha = 0 SYMMETRY NULL IS GONE. On the NACA case upper and lower surfaces had to
coincide, which was a free correctness check. At incidence there is no equivalent, so the
C_N closure and the shock position carry that weight instead -- and C_N IS NOT C_L. The
experiment reports force NORMAL TO THE CHORD; converting needs the axial force:
    C_L = C_N cos(a) - C_A sin(a)
At 2.5 deg the two agree to 0.02%, which is a COINCIDENCE OF THE ANGLE and not a licence to
compare them directly. The comparison script does the conversion explicitly.
"""
import argparse
import json
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
BASE = REPO / "cases/tutorials/naca0012_te60"
SECTION = REPO / "validation_data/rae2822/rae2822_section.dat"
# BLUNT-TE VARIANT, 0.263% c gap against the sharp section's 0.013%.
# cgrid2d.py records the reason: "O-topology folds the cell at a CUSPED sharp trailing edge.
# Measured on RAE 2822, whose TE included angle is 8.59 deg." A cusp puts a near-degenerate
# cell at the one place a transonic wall-resolved solution is already hardest -- the wake
# origin. Blunting it is standard practice in published RAE 2822 CFD, not a fudge.
SECTION_BLUNT = REPO / "validation_data/rae2822/rae2822_section_te025.dat"

# Decoded from AGARD f8621 by scripts/parse_rae2822.py. alpha is the CORRECTED angle.
CASES = {
    "case07": dict(M=0.725, alpha=2.21, Re=6.5e6, CN=0.658, CM=-0.090, CD=0.0107),
    "case06": dict(M=0.725, alpha=2.54, Re=6.5e6, CN=0.743, CM=-0.095, CD=0.0127),
    "case09": dict(M=0.730, alpha=2.79, Re=6.5e6, CN=0.803, CM=-0.099, CD=0.0168),
}
CHORD, T_INF, P_INF, GAMMA = 1.0, 288.15, 1.0e5, 1.4
# NO R HERE ON PURPOSE: it is read from each case's thermophysicalProperties.


def section(blunt=False):
    """Read the decoded ordinates and return an ordered closed loop, plus its extrema."""
    xy = []
    src = SECTION_BLUNT if blunt else SECTION
    for l in src.read_text().splitlines():
        s = l.strip()
        if not s or s.startswith("#"):
            continue
        p = s.replace(",", " ").split()
        if len(p) >= 2:
            try:
                xy.append((float(p[0]), float(p[1])))
            except ValueError:
                pass
    a = np.array(xy)
    if a[:, 0].max() > 1.5:                       # ordinates in per-cent chord
        a = a / 100.0
    return a



def resample(loop, n_out=1200):
    """Spline-resample the measured ordinates onto a fine cosine-clustered loop.

    WHY THIS IS NEEDED. blockMesh PROJECTS its vertices onto the triSurface, so the surface
    must be FINER than the mesh that snaps to it. The decoded section carries 128 points at
    mean spacing 0.0159 c while our surface cells run 0.0042-0.0235 c: the mesh was finer than
    its own geometry, and the projection produced 150 NEGATIVE-VOLUME CELLS, 901 wrongly
    oriented faces and non-orthogonality 94. The NACA obj that meshes cleanly carries 7,996
    points for the same reason.

    THE MEASURED ORDINATES REMAIN THE SOURCE OF TRUTH. The spline only adds resolution between
    them; it must not move them, and the null below asserts that it does not.

    Cosine clustering in arc length puts points where curvature is: the leading edge and the
    sharp trailing edge, which are exactly where a uniform resample starves the projection.
    """
    # TEST SELF-INTERSECTION, WHICH IS THE ACTUAL PROPERTY, not a proxy for it.
    #
    # Two proxies were tried and both were wrong. POLAR-ANGLE MONOTONICITY assumes the section
    # is star-shaped about its centroid; a thin cambered aerofoil with a blunt trailing edge is
    # not, and it fired on the project's own faired EET baseline at 53 of 512 steps.
    # X-MONOTONICITY is worse: it would PASS the very bug this guard exists for, because
    # sorting by sign of y and then by x still gives x descending once and ascending once.
    #
    # A correctly ordered aerofoil loop is a SIMPLE POLYGON. The scrambled one is not: points
    # taken from the wrong surface make segments cross. Testing that directly is exact, costs
    # one vectorised pass over 512 points, and has no shape assumption in it at all.
    # DEGENERATE SEGMENTS FIRST. Sorting by sign of y can place the same point in both
    # halves, giving a zero-length segment. The arc-length parameter then stops increasing
    # and CubicSpline raises "`x` must be strictly increasing" -- a true rejection, but from
    # scipy and with no mention of the section. Report it here, in the section's own terms.
    seg0 = np.hypot(*np.diff(np.vstack([loop, loop[:1]]), axis=0).T)
    if (seg0 < 1e-12).any():
        k = int(np.argmin(seg0))
        raise SystemExit("section has %d zero-length segment(s), first at index %d "
                         "(x/c %.5f): the loop revisits a point"
                         % (int((seg0 < 1e-12).sum()), k, loop[k][0]))

    P = np.vstack([loop, loop[:1]])
    A, B = P[:-1], P[1:]
    n = len(A)
    def _cross(o, a, b):
        return (a[..., 0] - o[..., 0]) * (b[..., 1] - o[..., 1]) - \
               (a[..., 1] - o[..., 1]) * (b[..., 0] - o[..., 0])
    i, j = np.triu_indices(n, k=2)
    keep = ~((i == 0) & (j == n - 1))          # the closing segment legitimately touches
    i, j = i[keep], j[keep]
    d1 = _cross(A[i], B[i], A[j]) * _cross(A[i], B[i], B[j])
    d2 = _cross(A[j], B[j], A[i]) * _cross(A[j], B[j], B[i])
    hits = int(((d1 < 0) & (d2 < 0)).sum())
    if hits:
        k = np.where((d1 < 0) & (d2 < 0))[0][0]
        raise SystemExit("section self-intersects: %d crossing segment pairs, first at "
                         "index %d x/c %.4f and index %d x/c %.4f"
                         % (hits, i[k], A[i[k]][0], j[k], A[j[k]][0]))

    # cumulative arc length as the spline parameter, periodic in the closed loop
    c = np.vstack([loop, loop[:1]])
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(c[:, 0]), np.diff(c[:, 1])))])
    s /= s[-1]
    # cluster toward s = 0 and s = 1, which is the trailing edge, and toward the LE at s ~ 0.5
    u = np.linspace(0, 1, n_out, endpoint=False)
    u = u - 0.06 * np.sin(4 * np.pi * u) / (4 * np.pi) * 4
    u = np.clip(np.sort(u % 1.0), 0, 1)
    try:
        from scipy.interpolate import CubicSpline
        fx = CubicSpline(s, c[:, 0], bc_type="periodic")
        fy = CubicSpline(s, c[:, 1], bc_type="periodic")
        out = np.column_stack([fx(u), fy(u)])
    except ImportError:
        out = np.column_stack([np.interp(u, s, c[:, 0]), np.interp(u, s, c[:, 1])])

    # NULL WITH AN INDEPENDENTLY KNOWN ANSWER (D058): evaluate the spline AT THE ORIGINAL
    # PARAMETERS and require it to return the original points. The answer is known because the
    # inputs are the answer, and it is independent of how densely we sample the output.
    #
    # MY FIRST VERSION MEASURED THE DISTANCE FROM EACH ORIGINAL POINT TO THE NEAREST OUTPUT
    # POINT and failed at 9.9e-04. That is not curve error, it is OUTPUT SPACING: at 1,200
    # points the loop is sampled every ~0.0017 c, so a point lying exactly on the curve still
    # sits up to half a spacing from its nearest neighbour. Same error as comparing a requested
    # first cell against a windowed mean -- TWO QUANTITIES WITH DIFFERENT DEFINITIONS COMPARED
    # AS THOUGH THEY WERE ONE (D065), and the second time today.
    try:
        chk = np.column_stack([fx(s[:-1]), fy(s[:-1])])
        d = np.hypot(*(chk - loop).T)
        how = "spline evaluated at the source parameters"
    except NameError:
        d = np.array([np.hypot(*(out - q).T).min() for q in loop])
        how = "nearest-output fallback (linear interpolation path)"
    if d.max() > 1e-9:
        raise SystemExit("resample moved the measured ordinates by up to %.2e c (%s)"
                         % (d.max(), how))
    # AND a separate, weaker check that the OUTPUT is dense enough for the mesh to project onto
    step = np.hypot(*np.diff(np.vstack([out, out[:1]]), axis=0).T)
    if step.max() > 0.004:
        raise SystemExit("resampled spacing %.4f c is coarser than the surface cells "
                         "(0.0042 c); blockMesh will project onto facets" % step.max())
    print("      resampled %d -> %d points; ordinates exact to %.1e c (%s); "
          "max output spacing %.5f c" % (len(loop), len(out), d.max(), how, step.max()))
    return out


def write_obj(a, path, span):
    """Extrude the section into a two-triangle-deep ribbon, the form blockMesh projects onto."""
    # THE FILE IS ALREADY AN ORDERED CLOSED LOOP. Its own header says so: "Selig order:
    # TE over the upper surface to LE, then lower back to TE." USE THE ORDER GIVEN.
    #
    # My first version split on the SIGN OF y and re-sorted each half by x. That works for a
    # symmetric section and SELF-INTERSECTS A CAMBERED ONE: the RAE 2822 carries 1.26% camber
    # at x/c 0.757, so 13 of its 65 LOWER-surface points sit ABOVE y = 0 and were sorted onto
    # the upper surface. The result still meshed, and blockMesh reported success -- the damage
    # showed only in the quality numbers, aspect 96,282 against 10,017 and non-orthogonality
    # 90.2 against 30.7.
    #
    # SAME SHAPE AS D065: an ordering was GIVEN and I reconstructed one instead. Reconstruction
    # is not the object, and it fails on exactly the inputs the reconstruction rule was not
    # written for.
    loop = a
    if np.hypot(*(loop[0] - loop[-1])) < 1e-9:    # drop the repeated closing point
        loop = loop[:-1]
    keep = [0] + [i for i in range(1, len(loop))
                  if np.hypot(*(loop[i] - loop[i - 1])) > 1e-9]
    loop = loop[keep]
    # ASSERT IT IS A SIMPLE LOOP, since the failure above was silent. A star-shaped section
    # traversed correctly sweeps a monotone polar angle about an interior point; a scrambled
    # one reverses. Count the reversals rather than trusting the read.
    c = loop.mean(axis=0)
    th = np.unwrap(np.arctan2(loop[:, 1] - c[1], loop[:, 0] - c[0]))
    d = np.diff(th)
    loop = resample(loop)
    n = len(loop)
    L = ["# RAE 2822, extruded from validation_data/rae2822/rae2822_section.dat",
         "# generated by scripts/build_rae2822_transonic.py -- do not edit",
         "g airfoil"]
    for y in (-span / 2, span / 2):
        for x, z in loop:
            L.append("v %.9f %.9f %.9f" % (x, y, z))
    for i in range(n):
        j = (i + 1) % n
        a0, b0, a1, b1 = i + 1, j + 1, i + 1 + n, j + 1 + n
        L.append("f %d %d %d" % (a0, b0, b1))
        L.append("f %d %d %d" % (a0, b1, a1))
    pathlib.Path(path).write_text("\n".join(L) + "\n")
    return n, loop



def polyline_edges(loop, x_sh, span):
    """Explicit polyLine edges along the aerofoil, replacing blockMesh's `project`.

    THE DEFECT THIS FIXES, verified rather than suspected. `project 8 5 (aerofoil)` asks
    blockMesh to place the edge from the UPPER shoulder to the trailing edge by projecting the
    STRAIGHT LINE between them onto the surface, NEAREST POINT WINS. The RAE 2822 carries 1.26%
    camber, which lifts the lower skin above z = 0 near the trailing edge, so beyond x/c 0.886
    that straight line is CLOSER TO THE LOWER SURFACE than the upper one. The edge folds
    through the aerofoil, and the radial column of cells above it inverts: 164 negative-volume
    cells, all at x/c 0.886, all on the upper surface.

    A SYMMETRIC SECTION CANNOT SHOW IT. Its chord line is equidistant from both skins by
    construction, so nearest-point projection is unambiguous everywhere and the NACA case gives
    no warning. Fourth instance today of a convention that is correct for the NACA 0012 and
    wrong for the first cambered aerofoil it meets.

    THE FIX IS TO STOP ASKING. A polyLine states where the edge goes; there is no nearest-point
    decision left to get wrong, on this or any other section.
    """
    n = len(loop)
    le = int(np.argmin(loop[:, 0]))
    # A BLUNT TRAILING EDGE HAS TWO CORNERS AND A BASE; THE C-GRID HAS ONE TE VERTEX.
    # blockMesh refuses the topology outright ("Block mesh topology incorrect") if the upper
    # and lower edges end at different points. Route each polyLine along its own surface AND
    # HALF THE BASE, so both converge on the same midpoint, which is the vertex the wake
    # starts from. A sharp TE is the degenerate case of this and is unaffected.
    te_i = int(np.argmax(loop[:, 0]))
    base = np.where(loop[:, 0] > loop[:, 0].max() - 1e-3)[0]
    if len(base) > 2:                       # genuinely blunt
        mid = 0.5 * (loop[base].max(axis=0) + loop[base].min(axis=0))
    else:
        mid = loop[te_i]
    upper = np.vstack([loop[:le + 1][::-1], mid])       # LE -> upper surface -> TE midpoint
    lower = np.vstack([loop[le:], mid])                 # LE -> lower surface -> TE midpoint

    def seg(surf, x0, x1):
        m = (surf[:, 0] > min(x0, x1) + 1e-6) & (surf[:, 0] < max(x0, x1) - 1e-6)
        s = surf[m]
        return s[np.argsort(s[:, 0])] if x1 > x0 else s[np.argsort(-s[:, 0])]

    def fmt(v, pts):
        return ("    polyLine %d %d\n    (\n" % v
                + "".join("        (%.8f %.8f %.8f)\n" % (x, y, z) for x, z in pts)
                + "    )\n")

    L = []
    for y, (v_le, v_lo, v_up, v_te) in ((-span / 2, (4, 7, 8, 5)),
                                        (+span / 2, (16, 19, 20, 17))):
        L.append(fmt((v_le, v_lo), seg(lower, 0.0, x_sh)))     # LE -> lower shoulder
        L.append(fmt((v_lo, v_te), seg(lower, x_sh, 1.0)))     # lower shoulder -> TE
        L.append(fmt((v_le, v_up), seg(upper, 0.0, x_sh)))     # LE -> upper shoulder
        L.append(fmt((v_up, v_te), seg(upper, x_sh, 1.0)))     # upper shoulder -> TE
        L = [l.replace("(%.8f " % 0, "(%.8f " % 0) for l in L]  # no-op, keeps formatting stable
    # y is not carried by seg(), so splice it in per block
    out = []
    for y, chunk in ((-span / 2, L[:4]), (+span / 2, L[4:])):
        for c in chunk:
            out.append(re.sub(r"\(([-0-9.]+) ([-0-9.]+) ([-0-9.]+)\)",
                              lambda m: "(%s %.6f %s)" % (m.group(1), y, m.group(3)), c))
    return "".join(out)


def build(tag, cfg, args):
    dst = REPO / "cases/tutorials" / ("rae2822_%s" % tag)
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(BASE, dst, ignore=shutil.ignore_patterns(
        "processor*", "postProcessing", "log.*", "qoi_normalisation.json"))
    for d in list(dst.iterdir()):
        if d.is_dir() and d.name not in ("0", "constant", "system"):
            shutil.rmtree(d)

    a = section()
    span = 0.2
    gdir = dst / "constant/geometry"
    gdir.mkdir(parents=True, exist_ok=True)
    for old in gdir.glob("*.obj"):
        old.unlink()
    npts, loop = write_obj(a, gdir / "RAE2822.obj", span)

    # --- flow state. Re fixes mu; M fixes U through the speed of sound. -----------
    # R COMES FROM THE CASE'S OWN molWeight, NOT FROM A CONSTANT IN THIS FILE. I hard-coded
    # 287.058 while the thermophysicalProperties give 287.70, a 0.22% gap that put Re at
    # 6.485e6 instead of 6.5e6 and M at 0.7242 instead of 0.725, and left forceCoeffs' rhoInf
    # disagreeing with p/(RT). The pre-run audit refused to launch on exactly that, which is
    # the whole reason it exists.
    tp_txt = (dst / "constant/thermophysicalProperties").read_text()
    R_case = 8314.47 / float(re.search(r"^\s*molWeight\s+([-0-9.eE+]+);", tp_txt, re.M).group(1))
    aS = math.sqrt(GAMMA * R_case * T_INF)
    U = cfg["M"] * aS
    rho = P_INF / (R_case * T_INF)
    mu = rho * U * CHORD / cfg["Re"]
    ar = math.radians(cfg["alpha"])
    Uvec = (U * math.cos(ar), 0.0, U * math.sin(ar))

    # --- blockMeshDict: new geometry, seeds from the section's own extrema --------
    bm = dst / "system/blockMeshDict"
    t = bm.read_text()
    t = t.replace('file   "NACA0012.obj";', 'file   "RAE2822.obj";')
    # THE SHOULDER VERTICES MUST SHARE A CHORDWISE STATION. xUpper and xLower are the C-grid
    # corners where the blocks split, and the NACA dict has both at 0.3 for that reason. I first
    # set them to the section's EXTREMA LOCATIONS, which for a CAMBERED aerofoil are different
    # x (0.4378 upper, 0.3530 lower). The upper and lower blocks then span different chordwise
    # ranges, the cells twist, and blockMesh produced 150 NEGATIVE-VOLUME cells and 901 wrongly
    # oriented faces -- while still reporting success, because blockMesh builds what it is told.
    #
    # A SYMMETRIC SECTION HIDES THIS COMPLETELY: max thickness above and below are at the same
    # x, so the wrong rule and the right rule agree. Third time today that a convention valid
    # for the NACA case failed on its first cambered input (D077).
    x_sh = 0.3
    ups = loop[loop[:, 1] >= 0]
    los = loop[loop[:, 1] < 0]
    # local surface height at the shoulder station, taken from the section itself
    zu = float(np.interp(x_sh, np.sort(ups[:, 0]), ups[np.argsort(ups[:, 0]), 1]))
    zl = float(np.interp(x_sh, np.sort(los[:, 0]), los[np.argsort(los[:, 0]), 1]))
    xu = xl = x_sh
    print("      shoulder vertices at x/c %.3f: z_upper %+.6f  z_lower %+.6f "
          "(max t/c %.4f at x/c %.3f, camber peak x/c %.3f)"
          % (x_sh, zu, zl, loop[:, 1].max() - loop[:, 1].min(),
             loop[np.argmax(loop[:, 1]), 0], loop[np.argmin(loop[:, 1]), 0]))
    for key, val in (("xUpper", xu), ("zUpper", zu), ("xLower", xl), ("zLower", zl)):
        t = re.sub(r"^(\s*%s\s+)[-0-9.eE+]+;" % key, r"\g<1>%.6f;" % val, t, flags=re.M)
    # REPLACE the four aerofoil `project` edges (and their front-plane twins) with polyLines
    ed = re.search(r"^edges\s*\n\((.*?)\n\);", t, re.S | re.M)
    if not ed:
        raise SystemExit("no edges block in blockMeshDict")
    body = ed.group(1)
    kept = "\n".join(l for l in body.splitlines()
                      if "(aerofoil)" not in l)          # keep the cylinder projections
    t = t[:ed.start()] + ("edges\n(\n%s\n%s);" % (kept, polyline_edges(loop, x_sh, span))) \
        + t[ed.end():]
    bm.write_text(t)

    # --- fvSchemes: LIMITED convection, because a shock rings on unlimited linear --
    fs = dst / "system/fvSchemes"
    s = fs.read_text()
    div = re.search(r"divSchemes\s*\n\{(.*?)\n\}", s, re.S)
    new_div = """divSchemes
{
    default         none;
    // TRANSONIC: unlimited linear rings across a shock. limitedLinear 1 is the
    // conservative TVD choice here; vanLeer on density is the usual pairing.
    div(phi,U)      Gauss limitedLinearV 1;
    div(phi,e)      Gauss limitedLinear 1;
    div(phi,K)      Gauss limitedLinear 1;
    div(phi,h)      Gauss limitedLinear 1;
    div(phi,k)      Gauss limitedLinear 1;
    div(phi,omega)  Gauss limitedLinear 1;
    // UPWIND ON THE PRESSURE CONVECTION TERM, and this is the one that matters.
    // The working rhoPimpleFoam tutorial upwinds exactly this term -- `div(phiv,p) Gauss
    // upwind` -- while keeping linearUpwind on momentum and energy. It carries the startup
    // pressure wave, the "sonic boom" that fires on the first few iterations before the
    // field settles, and second order on it is what blew ours up: bounding omega at 3.7e8 on
    // iteration 1. First order here costs accuracy nowhere that matters at convergence,
    // because at steady state the term is small.
    div(phid,p)     Gauss upwind;
    div(phi,Ekp)    Gauss limitedLinear 1;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}"""
    s = s[:div.start()] + new_div + s[div.end():]
    fs.write_text(s)

    # --- fvSolution: transonic pressure equation ----------------------------------
    fv = dst / "system/fvSolution"
    v = fv.read_text()
    if re.search(r"^\s*transonic\s+", v, re.M):
        v = re.sub(r"^(\s*transonic\s+)\w+;", r"\g<1>yes;", v, flags=re.M)
    else:
        v = re.sub(r"(SIMPLE\s*\n\{\s*\n)", r"\1    transonic       yes;\n", v, count=1)
    # TRANSONIC MAKES THE PRESSURE MATRIX ASYMMETRIC, so GAMG cannot solve it.
    # `transonic yes` adds the div(phid,p) convection term to the pressure equation. GAMG is a
    # SYMMETRIC-matrix solver; handed an asymmetric one it dies in sumProd at the coarsest
    # level with a floating-point exception on iteration 1, which is exactly how all three
    # cases failed. PBiCGStab/DILU is the asymmetric equivalent.
    #
    # This is D013's lesson arriving by a second route: there GAMG STALLED on extreme
    # anisotropy, here it CANNOT APPLY at all. Same remedy, different reason, and the reason
    # matters because the subcritical NACA case runs GAMG on the same mesh perfectly well.
    pblk = re.search(r"^\s{4}p\s*\n\s{4}\{.*?\n\s{4}\}", v, re.S | re.M)
    if pblk and "GAMG" in pblk.group(0):
        v = v[:pblk.start()] + """    p
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-09;
        relTol          0.01;
    }""" + v[pblk.end():]
    fv.write_text(v)

    # ARG-146 FROZE CLASSIC p 0.3 / U 0.7. I built these from te60, which carries the
    # TUTORIAL pairing p 0.7 / U 0.3; te60_relax is the case that has the frozen one. Copying
    # "the best case" picked up a setting the recipe had already superseded.
    v = fv.read_text()
    v = re.sub(r"^(\s*p\s+)0\.7;", r"\g<1>0.3;", v, flags=re.M)
    v = re.sub(r"^(\s*U\s+)0\.3;", r"\g<1>0.7;", v, flags=re.M)
    fv.write_text(v)

    # --- 0/ fields: incidence enters through the FREESTREAM DIRECTION -------------
    u = dst / "0/U"
    ut = u.read_text()
    ut = re.sub(r"^(Uinlet\s+)\([^)]*\);",
                r"\g<1>(%.6f %.6f %.6f);" % Uvec, ut, flags=re.M)
    u.write_text(ut)
    for f, val in (("0/p", P_INF), ("0/T", T_INF)):
        p = dst / f
        txt = p.read_text()
        m = re.search(r"^(\w+)\s+[-0-9.eE+]+;", txt, re.M)
        if m and m.group(1) not in ("version",):
            txt = re.sub(r"^%s\s+[-0-9.eE+]+;" % m.group(1),
                         "%s %.6e;" % (m.group(1), val), txt, count=1, flags=re.M)
        p.write_text(txt)

    tp = dst / "constant/thermophysicalProperties"
    tt = tp.read_text()
    tt = re.sub(r"^(\s*mu\s+)[-0-9.eE+]+;", r"\g<1>%.8e;" % mu, tt, flags=re.M)
    tp.write_text(tt)

    # --- controlDict: force axes ROTATE WITH THE FLOW -----------------------------
    cd = dst / "system/controlDict"
    ct = cd.read_text()
    ct = re.sub(r"^(\s*magUInf\s+)[-0-9.eE+]+;", r"\g<1>%.6f;" % U, ct, flags=re.M)
    ct = re.sub(r"^(\s*rhoInf\s+)[-0-9.eE+]+;", r"\g<1>%.6f;" % rho, ct, flags=re.M)
    ct = re.sub(r"^(\s*liftDir\s+)\([^)]*\);",
                r"\g<1>(%.8f 0 %.8f);" % (-math.sin(ar), math.cos(ar)), ct, flags=re.M)
    ct = re.sub(r"^(\s*dragDir\s+)\([^)]*\);",
                r"\g<1>(%.8f 0 %.8f);" % (math.cos(ar), math.sin(ar)), ct, flags=re.M)
    ct = re.sub(r"^(\s*endTime\s+)[0-9]+;", r"\g<1>%d;" % args.end_time, ct, flags=re.M)
    cd.write_text(ct)

    r = subprocess.run(["bash", "-lc", FOAM + "blockMesh > log.blockMesh 2>&1 && "
                        "checkMesh -allGeometry > log.checkMesh 2>&1"],
                       cwd=str(dst), capture_output=True, text=True, timeout=3600)
    ok = (dst / "constant/polyMesh/points").exists()
    print("  %-18s M %.3f  a %.2f deg  U %.3f  mu %.3e  Re %.2e   obj %d pts   %s"
          % (tag, cfg["M"], cfg["alpha"], U, mu, cfg["Re"], npts,
             "meshed" if ok else "*** blockMesh FAILED ***"))
    if not ok:
        print("\n".join("      " + l for l in
                        (dst / "log.blockMesh").read_text(errors="replace").splitlines()[-10:]))
        return None
    cm = (dst / "log.checkMesh").read_text(errors="replace")
    nc = re.search(r"cells:\s+(\d+)", cm)
    ar_ = re.search(r"Max aspect ratio: ([0-9.]+)", cm)
    no = re.search(r"non-orthogonality Max: ([0-9.]+)", cm)
    print("      cells %s  aspect %s  nonortho %s  checkMesh %s"
          % (nc.group(1) if nc else "?", ar_.group(1) if ar_ else "?",
             no.group(1) if no else "?", "OK" if "Mesh OK" in cm else "SEE LOG"))
    return dst


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cases", nargs="*", default=list(CASES))
    ap.add_argument("--end-time", type=int, default=80000)
    a = ap.parse_args()
    made = []
    for tag in a.cases:
        d = build(tag, CASES[tag], a)
        if d:
            made.append(d)
    print("\n  built %d of %d" % (len(made), len(a.cases)))
    return 0 if len(made) == len(a.cases) else 1


if __name__ == "__main__":
    raise SystemExit(main())
