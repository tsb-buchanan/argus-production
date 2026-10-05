#!/usr/bin/env python3
"""hinge_moment_integrate.py: the ARGUS hinge moment by RANS surface integration.

THE SANCTIONED ROUTE, AND THE ONLY ONE. scripts/hinge_moment.py refuses every other:
no cmy, no transfer formula, no exception (project decision, 2026-07-29). That module
describes the route and computes nothing; this one computes it.

THE SPEC, as agreed with the DSO team and quoted from their delivery README:

    the aerodynamic moment about the swept local x_h/c = 0.62 line, projected about
    its local tangent, using only surfaces aft of the line over 0.60 <= eta <= 1.00.
    Positive is trailing-edge down. Pressure, viscous and total reported separately as:
      - sectional moment per unit span            m'_h   [N m/m]
      - c_mh = m'_h / (q_inf c(y)^2)                     [-]
      - integrated half-wing moment               M_h    [N m]
      - C_Mh,region = M_h / [q_inf integral c(y)^2 dy]   [-]

FOUR THINGS THIS GETS RIGHT ON PURPOSE, each because the obvious version is wrong.

 1. THE HINGE LINE IS SWEPT, so its tangent is NOT the y axis. x_h(y) = x_LE(y) +
    0.62 c(y), and on a swept, tapered wing that line runs at an angle. Projecting the
    moment on the y axis instead of the LOCAL TANGENT is the error a careful person
    makes first, and on this planform it is not small.

 2. "SURFACES AFT OF THE LINE" IS A PER-STATION TEST, not a global x cut. The wing is
    swept and tapered, so a single x threshold includes forward surface at the root and
    excludes aft surface at the tip. Each face is tested against x_h AT ITS OWN y.

 3. THE MOMENT ARM IS MEASURED FROM THE LINE, IN THE PLANE NORMAL TO IT. Using the
    full 3D vector from an arbitrary point on the line double-counts the spanwise
    component, which the tangential projection then partly removes: right answer by a
    wrong route, and only by accident.

 4. AREAS AND NORMALS COME FROM THE POLYGONS, not from a nominal cell size. A faceted
    surface under-reads area (D060), so an area taken any other way is both wrong and
    wrong in a known direction.

THE NULL HAS AN INDEPENDENTLY KNOWN ANSWER (D058). A flat plate at zero incidence
carries equal and opposite pressure on its two faces, so its hinge moment is EXACTLY
ZERO by symmetry, analytically, with no reference solution required. --self-test builds
one and asserts it, and also asserts three cases where the answer is known to be
NON-zero, because a null that can only pass is doing no work.

THE VISCOUS TERM IS NOW COMPUTABLE AND ITS SIGN IS THE TRAP. OpenFOAM's
`wallShearStress` is the NEGATIVE of the traction on the body: forces.C dots devRhoReff
with +Sf and wallShearStress.C with -Sf/|Sf|, so the two objects disagree in sign on the
same tensor and neither says so. Measured on the ONERA M6 patch at a ratio of exactly
-1.0000 in all three components. wall_shear_to_body_traction() does the negation and
assert_friction_opposes_flow() refuses a wall shear that integrates to a thrust, which
is impossible whatever the geometry or the model.

Usage:
    python3 scripts/hinge_moment_integrate.py --self-test
    # q IS CONDITION WT, 0.5*1.225*40.8^2 = 1019.59 Pa. An earlier version of this line
    # showed 848.9, which belongs to no condition in this project and does not reproduce
    # the committed result; a usage example is copied far more often than it is checked.
    python3 scripts/hinge_moment_integrate.py --vtk results/vtk_export/baseline_wing_surface.vtk \\
        --q 1019.59 --out results/hinge_moment_baseline.json
    # with the viscous term, on a case run with argusPostPro installed:
    python3 scripts/hinge_moment_integrate.py --vtk <case>/postProcessing/argusWingSurface/<t>/cp_wingSurface.vtk \\
        --tau-vtk <case>/postProcessing/argusWingSurface/<t>/wallShearStress_wingSurface.vtk \\
        --q <q> --alpha <deg>
"""
import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

X_H_OVER_C = 0.62
ETA_LO, ETA_HI = 0.60, 1.00
SEMISPAN_DEFAULT = 1.8288          # bref/2, m


class HingeError(RuntimeError):
    pass


# ----------------------------------------------------------------- geometry ----
def polygon_areas_normals(points, polys):
    """Area vector and centroid of each polygon, by the shoelace/fan construction.

    Sf points along the polygon's own winding. Its MAGNITUDE is the area and its
    DIRECTION is the outward normal only if the source winding is outward; that is
    asserted by the caller against the enclosed volume, never assumed here.

    VECTORISED BY VERTEX COUNT, AND THE ARITHMETIC IS UNCHANGED. This was a pure Python
    double loop: 5.28e6 polygons times four vertices is 21e6 inner iterations, each calling
    np.cross and np.linalg.norm on a single 3-vector, where the per-call overhead is far
    larger than the arithmetic. MEASURED at 65.6 s per 400k faces, which is 865 s for one
    wing surface: the whole 17-minute cost of a case was this function.

    THE PROJECT ALREADY KNEW. case_derived_quantities.poly_area_normal is the vectorised
    twin and its docstring says so outright, that a Python loop over 5.2 M faces takes
    minutes while grouping by vertex count takes seconds. The fix landed in one copy and
    never reached this one, which is precisely the failure this module's own VTK-reader
    comment warns about: duplicated code diverges on the first fix.

    THE TWO ARE STILL NOT INTERCHANGEABLE, so this is a rewrite and not a substitution.
    Their AREA VECTORS are identical, measured at max |diff| = 0.000e+00 over 400k faces.
    Their CENTROIDS ARE NOT: poly_area_normal takes the plain vertex mean, this function
    takes the AREA-WEIGHTED centroid of the fan triangles and divides by |sum(Sf)|. The
    hinge moment applies its moment arm at that centroid, so swapping the definitions would
    move the answer silently. Every line below reproduces the loop it replaces, term for
    term, and the regression is the two hinge moments the slow path already produced.
    """
    sizes = np.fromiter((len(p) for p in polys), dtype=np.int64, count=len(polys))
    cen = np.empty((len(polys), 3))
    Sf = np.empty((len(polys), 3))
    for k in np.unique(sizes):
        kk = int(k)
        idx = np.flatnonzero(sizes == k)
        conn = np.fromiter((v for i in idx for v in polys[i]),
                           dtype=np.int64, count=len(idx) * kk).reshape(len(idx), kk)
        q = points[conn]                                  # (m, k, 3)
        c0 = q.mean(axis=1)                               # (m, 3)  the fan apex
        a = q - c0[:, None, :]                            # p[i]   - c0
        b = np.roll(a, -1, axis=1)                        # p[i+1] - c0, wrapping
        t = 0.5 * np.cross(a, b)                          # (m, k, 3) fan triangle areas
        S = t.sum(axis=1)                                 # (m, 3)  == the loop's `a`
        Sf[idx] = S
        w = np.linalg.norm(t, axis=2)                     # (m, k)  == |t| per triangle
        tri_c = (q + np.roll(q, -1, axis=1) + c0[:, None, :]) / 3.0
        cw = (w[:, :, None] * tri_c).sum(axis=1)          # (m, 3)  == the loop's `cw`
        n = np.linalg.norm(S, axis=1)                     # (m,)    == |a|, the divisor
        good = n > 0
        cen[idx] = np.where(good[:, None],
                            cw / np.where(good[:, None], n[:, None], 1.0),
                            c0)
    return cen, Sf


def planform(points, nbin=200):
    """x_LE(y), x_TE(y), c(y) from the surface itself, binned in y.

    DERIVED FROM THE SURFACE, not from the registry, because the hinge line must sit
    on the geometry that was actually meshed and integrated. A planform taken from a
    sidecar can disagree with the STL and nothing would surface it (GEO-102).
    """
    y = points[:, 1]
    lo, hi = y.min(), y.max()
    edges = np.linspace(lo, hi, nbin + 1)
    yc, xle, xte = [], [], []
    for i in range(nbin):
        m = (y >= edges[i]) & (y < edges[i + 1] if i < nbin - 1 else y <= edges[i + 1])
        if m.sum() < 3:
            continue
        yc.append(0.5 * (edges[i] + edges[i + 1]))
        xle.append(points[m, 0].min())
        xte.append(points[m, 0].max())
    if len(yc) < 10:
        raise HingeError("planform: fewer than 10 populated spanwise bins")
    return np.array(yc), np.array(xle), np.array(xte)


def hinge_line(yc, xle, xte):
    """x_h(y) and the unit tangent of the swept hinge line at each station."""
    xh = xle + X_H_OVER_C * (xte - xle)
    dxh = np.gradient(xh, yc)                       # dx_h/dy, the sweep of the line
    t = np.stack([dxh, np.ones_like(dxh), np.zeros_like(dxh)], axis=1)
    t /= np.linalg.norm(t, axis=1)[:, None]
    return xh, t


# ------------------------------------------------------------------- moment ----
def hinge_moment(points, polys, p_face, tau_face=None, semispan=SEMISPAN_DEFAULT,
                 q=None, eta_lo=ETA_LO, eta_hi=ETA_HI, nbin=200):
    """Return the four agreed forms, pressure and viscous separated.

    p_face   : static pressure on each face [Pa]   (gauge; a uniform offset cancels
               on a CLOSED region but NOT on an open aft-of-line region, so the caller
               must supply gauge pressure and this is asserted in the self-test)
    tau_face : wall shear stress vector on each face [Pa], or None to omit the
               viscous term. OMITTED IS REPORTED, NEVER SILENTLY ZERO.
    """
    cen, Sf = polygon_areas_normals(points, polys)
    area = np.linalg.norm(Sf, axis=1)
    if not np.all(area > 0):
        raise HingeError("degenerate polygon with zero area")

    yc, xle, xte = planform(points, nbin=nbin)
    xh, that = hinge_line(yc, xle, xte)
    chord = xte - xle

    # per-face station lookup
    j = np.clip(np.searchsorted(yc, cen[:, 1]) - 1, 0, len(yc) - 1)
    xh_f = xh[j]
    t_f = that[j]
    c_f = chord[j]
    eta_f = cen[:, 1] / semispan

    # ITEM 2: aft-of-line is tested per station, against x_h at the face's OWN y.
    sel = (cen[:, 0] >= xh_f) & (eta_f >= eta_lo) & (eta_f <= eta_hi)
    if sel.sum() == 0:
        raise HingeError("no faces aft of the hinge line inside the eta band")

    # ITEM 3: arm measured FROM THE LINE, in the plane normal to it.
    origin = np.stack([xh_f, cen[:, 1], np.zeros_like(xh_f)], axis=1)
    r = cen - origin
    r -= (np.einsum("ij,ij->i", r, t_f))[:, None] * t_f      # strip the spanwise part

    # pressure force on a face is -p n dA ; Sf already carries n dA
    Fp = -p_face[:, None] * Sf
    Mp = np.cross(r, Fp)
    mp = np.einsum("ij,ij->i", Mp, t_f)                      # ITEM 1: local tangent

    if tau_face is not None:
        Fv = tau_face * area[:, None]
        mv = np.einsum("ij,ij->i", np.cross(r, Fv), t_f)
    else:
        mv = np.zeros_like(mp)

    Mp_tot = float(mp[sel].sum())
    Mv_tot = float(mv[sel].sum())

    out = {
        "M_h_pressure_Nm": Mp_tot,
        "M_h_viscous_Nm": (Mv_tot if tau_face is not None else None),
        "M_h_total_Nm": (Mp_tot + Mv_tot if tau_face is not None else None),
        "viscous_included": tau_face is not None,
        "n_faces_total": int(len(polys)),
        "n_faces_in_region": int(sel.sum()),
        "region": {"x_h_over_c": X_H_OVER_C, "eta": [eta_lo, eta_hi]},
        # THE SPANWISE BIN COUNT IS A RESULT, NOT A SETTING, so it is recorded with the
        # number it produced. planform() takes min(x) and max(x) WITHIN each bin, so a
        # coarse bin overstates the local chord, which moves x_h = x_LE + 0.62c against the
        # surface and changes WHICH FACES COUNT AS AFT OF THE LINE. Measured on SWB_trim:
        # nbin 100/200/400/800/1600 gives M_h 0.9375/0.8972/0.8774/0.8677/0.8630 with
        # n_faces_in_region falling 523627 to 501559 and the mean chord falling 0.2217 to
        # 0.2135, while the mean x_h holds to 0.015 per cent. The increments halve on each
        # doubling, so this is first-order convergent and Richardson puts the limit near
        # 0.858: nbin = 200 was 4.5 per cent high in magnitude.
        "n_spanwise_bins": int(nbin),
        "sign_convention": "positive is trailing-edge down, about the local tangent "
                           "of the swept x_h/c line",
    }

    # sectional and normalised forms
    if q is not None:
        band = (yc / semispan >= eta_lo) & (yc / semispan <= eta_hi)
        c2 = np.trapz((chord[band]) ** 2, yc[band])
        out["q_Pa"] = q
        out["integral_c2_dy_m3"] = float(c2)
        out["C_Mh_region"] = Mp_tot / (q * c2) if c2 > 0 else None
        if tau_face is not None:
            out["C_Mh_region_total"] = (Mp_tot + Mv_tot) / (q * c2) if c2 > 0 else None
        # sectional m'_h and c_mh, per spanwise bin
        sect = []
        for i in np.where(band)[0]:
            m = sel & (j == i)
            if m.sum() == 0:
                continue
            dy = (yc[min(i + 1, len(yc) - 1)] - yc[max(i - 1, 0)]) / 2.0
            if dy <= 0:
                continue
            mprime = float((mp[m] + mv[m]).sum()) / dy
            sect.append({"eta": float(yc[i] / semispan), "y_m": float(yc[i]),
                         "chord_m": float(chord[i]),
                         "m_prime_h_Nm_per_m": mprime,
                         "c_mh": mprime / (q * chord[i] ** 2)})
        out["sectional"] = sect
    return out


def wall_shear_to_body_traction(tau):
    """OpenFOAM's `wallShearStress` NEGATED, which is the traction ON THE BODY.

    TWO FUNCTION OBJECTS, ONE TENSOR, OPPOSITE SIGNS, AND NEITHER SAYS SO. Both start
    from devRhoReff, and then:

        forces.C:819            fT  = Sfb & devRhoReffb          <- Sf points OUT
        wallShearStress.C:86    wss = (-Sfp/magSfp) & Reffp      <- unit normal points IN

    so fT = -|Sf| * wss, face by face. `forces` is the object that reports drag, so it is
    the one giving the force on the body; `wallShearStress` is its negative.

    MEASURED, NOT INFERRED, AND THE MEASUREMENT IS THE AUTHORITY HERE. Reading the source
    first, I concluded the two agreed and was wrong. Integrating the sampled wall-shear
    over the ONERA M6 patch gives (-186.918, -2.381, -4.344) N against the forces object's
    (+186.918, +2.381, +4.344) N: a ratio of exactly -1.0000 in all three components. A
    correct answer reached by an invalid route is a defect, and so is an incorrect one
    (GEO-080 item 5); the check is what settled it, and the source then explained it.
    """
    return -np.asarray(tau, float)


def assert_friction_opposes_flow(tau_body, area, u_dir, label=""):
    """Skin friction can only RESIST the motion. A thrust means the sign is wrong.

    A D060 SIGN TEST WITH A DIRECTION PHYSICS FIXES, not a tolerance. The integrated
    viscous traction on a body in a flow has a positive component along the freestream:
    that is true whatever the geometry, the model or the mesh, so a negative one is
    impossible rather than merely surprising. It is the one check that catches the
    convention error above no matter which caller makes it.
    """
    u = np.asarray(u_dir, float)
    u = u / np.linalg.norm(u)
    D = float(np.dot((np.asarray(tau_body) * np.asarray(area)[:, None]).sum(axis=0), u))
    if D <= 0:
        raise HingeError(
            "%sthe supplied wall shear integrates to a THRUST of %.3f N along the "
            "freestream. Skin friction cannot propel a body. The sign convention is "
            "inverted: OpenFOAM's `wallShearStress` is the NEGATIVE of the traction on "
            "the body, so pass it through wall_shear_to_body_traction() first."
            % (label and label + ": ", -D))
    return D


# -------------------------------------------------------------------- nulls ----
def _box(nx=30, ny=10, chord=1.0, span=2.0, thick=0.10):
    """A closed rectangular section: upper at +t/2, lower at -t/2, plus the four ends.

    THE FLAT PLATE CANNOT TEST THE VISCOUS TERM AND THAT IS WHY THIS EXISTS. On a
    zero-thickness plate every face centre lies at z = 0, so the moment arm has no z
    component, and a chordwise traction produces NO moment about a spanwise line whatever
    its magnitude: the viscous answer is zero for a right reason and for a wrong one
    alike. Giving the section a thickness puts the two surfaces at z = +/- t/2 and makes
    the viscous moment analytic and non-zero.
    """
    xs = np.linspace(0.0, chord, nx + 1)
    ys = np.linspace(0.0, span, ny + 1)
    pts, polys, idx = [], [], {}
    for s, zv in ((0, +thick / 2), (1, -thick / 2)):
        for iy, yv in enumerate(ys):
            for ix, xv in enumerate(xs):
                idx[(s, ix, iy)] = len(pts)
                pts.append([xv, yv, zv])
    pts = np.array(pts)
    for s, sgn in ((0, +1), (1, -1)):
        for iy in range(ny):
            for ix in range(nx):
                a = idx[(s, ix, iy)]; b = idx[(s, ix + 1, iy)]
                c = idx[(s, ix + 1, iy + 1)]; d = idx[(s, ix, iy + 1)]
                polys.append([a, b, c, d] if sgn > 0 else [a, d, c, b])
    n_lift = len(polys)          # the lifting faces; the caps follow and are unloaded
    for iy in range(ny):         # leading and trailing end caps
        for (ix, w) in ((0, -1), (nx, +1)):
            a = idx[(0, ix, iy)]; b = idx[(0, ix, iy + 1)]
            c = idx[(1, ix, iy + 1)]; d = idx[(1, ix, iy)]
            polys.append([a, b, c, d] if w > 0 else [a, d, c, b])
    return pts, polys, n_lift


def _plate(nx=40, ny=20, chord=1.0, span=2.0, alpha_deg=0.0):
    """A closed zero-thickness plate: top and bottom sheets with opposite winding."""
    xs = np.linspace(0.0, chord, nx + 1)
    ys = np.linspace(0.0, span, ny + 1)
    pts, polys = [], []
    idx = {}
    for s, sgn in ((0, +1), (1, -1)):
        for iy, yv in enumerate(ys):
            for ix, xv in enumerate(xs):
                idx[(s, ix, iy)] = len(pts)
                pts.append([xv, yv, 0.0])
    pts = np.array(pts)
    for s, sgn in ((0, +1), (1, -1)):
        for iy in range(ny):
            for ix in range(nx):
                a = idx[(s, ix, iy)]; b = idx[(s, ix + 1, iy)]
                c = idx[(s, ix + 1, iy + 1)]; d = idx[(s, ix, iy + 1)]
                polys.append([a, b, c, d] if sgn > 0 else [a, d, c, b])
    return pts, polys


def self_test():
    ok = True

    def check(name, got, want, tol, note=""):
        nonlocal ok
        good = abs(got - want) <= tol
        ok = ok and good
        print("  [%s] %-46s got %+.6e  want %+.6e  tol %.0e %s"
              % ("ok " if good else "FAIL", name, got, want, tol, note))

    # 1. FLAT PLATE, ZERO INCIDENCE, EQUAL PRESSURE BOTH SIDES -> EXACTLY ZERO.
    #    Known analytically by symmetry; no reference solution needed (D058).
    pts, polys = _plate()
    p = np.full(len(polys), 1234.5)                     # same on both sheets
    r = hinge_moment(pts, polys, p, semispan=2.0, q=100.0)
    check("flat plate, uniform p, both faces", r["M_h_pressure_Nm"], 0.0, 1e-9,
          "symmetry")

    # 2. NEGATIVE CONTROL: load ONE side only. The answer must NOT be zero, or the
    #    test above is passing for the wrong reason (a gate that cannot fail).
    n = len(polys) // 2
    p2 = np.concatenate([np.full(n, 2000.0), np.full(len(polys) - n, 1000.0)])
    r2 = hinge_moment(pts, polys, p2, semispan=2.0, q=100.0)
    nz = abs(r2["M_h_pressure_Nm"]) > 1e-6
    ok = ok and nz
    print("  [%s] %-46s got %+.6e  (must be non-zero)"
          % ("ok " if nz else "FAIL", "one-sided load gives a NON-zero moment",
             r2["M_h_pressure_Nm"]))

    # 3. SIGN. Higher pressure on the UPPER surface pushes the trailing edge DOWN,
    #    which the agreed convention calls POSITIVE.
    good = r2["M_h_pressure_Nm"] > 0
    ok = ok and good
    print("  [%s] %-46s %+.6e > 0"
          % ("ok " if good else "FAIL", "upper-surface load -> TE down -> positive",
             r2["M_h_pressure_Nm"]))

    # 4. REGION. Only faces aft of x_h/c = 0.62 inside 0.60 <= eta <= 1.00 count.
    #    Plate chord 1.0, span 2.0, semispan 2.0 -> eta band is the outer 40%, and
    #    aft of the line is 38% of chord. Expected fraction of a 2-sheet grid:
    exp_frac = (1 - X_H_OVER_C) * (ETA_HI - ETA_LO)
    got_frac = r["n_faces_in_region"] / r["n_faces_total"]
    check("region selects the right face fraction", got_frac, exp_frac, 0.03,
          "(1-0.62)x(1.00-0.60)")

    # 5. AREA. The plate is two sheets of chord x span, so the total is 2*c*s exactly.
    cen, Sf = polygon_areas_normals(pts, polys)
    check("polygon areas sum to the analytic plate area",
          float(np.linalg.norm(Sf, axis=1).sum()), 2 * 1.0 * 2.0, 1e-9)

    # ---------------- the VISCOUS term, which had never been executed ----------------
    # 6. AN ANALYTIC VISCOUS MOMENT. On a box of thickness t, uniform chordwise traction
    #    tau_u on the upper surface at z = +t/2 and tau_l on the lower at z = -t/2 gives
    #    M_y = sum over faces of (r x F)_y = sum r_z F_x = (t/2)*A*(tau_u - tau_l),
    #    with A the loaded area of ONE surface. It is independent of where the hinge line
    #    sits, because the force is purely chordwise. Nothing in this file computes it.
    C, S, T = 1.0, 2.0, 0.10
    # nx = 50 AND ny = 10 ARE CHOSEN SO THE REGION BOUNDARIES FALL ON CELL EDGES.
    # x_h/c = 0.62 is an edge at 31/50 and eta = 0.60 is an edge at 6/10, so the selected
    # area is exactly (1-0.62)*(1.00-0.60) of the plan and the analytic answer below is
    # exact. At nx = 30 the hinge line lands mid-cell, the selection is 0.3667c rather
    # than 0.38c, and the check failed by 3.5% -- the FIXTURE was wrong, not the code.
    bpts, bpolys, nlift = _box(nx=50, ny=10, chord=C, span=S, thick=T)
    nhalf = nlift // 2
    tau_u, tau_l = 3.0, 1.0
    tb = np.zeros((len(bpolys), 3))
    tb[:nhalf, 0] = tau_u                       # upper sheet
    tb[nhalf:nlift, 0] = tau_l                  # lower sheet
    p0 = np.zeros(len(bpolys))
    # only the part of the box inside the region contributes, so the analytic answer is
    # scaled by the same aft-of-line and eta selection the integrator applies
    frac = (1 - X_H_OVER_C) * (ETA_HI - ETA_LO)
    want = (T / 2) * (C * S * frac) * (tau_u - tau_l)
    rb = hinge_moment(bpts, bpolys, p0, tau_face=tb, semispan=S, q=100.0)
    check("box, uniform shear tau_u=3 tau_l=1, analytic (t/2)A(du)",
          rb["M_h_viscous_Nm"], want, abs(want) * 1e-9, "(t/2)*A*(tau_u-tau_l)")

    # 7. SYMMETRY NULL, AND IT IS THE ONE THE FLAT PLATE COULD ONLY EVER GIVE. Equal
    #    shear on both surfaces must cancel exactly. Reported beside 6 so that a passing
    #    zero cannot be mistaken for a working integrator.
    ts = np.zeros((len(bpolys), 3))
    ts[:nlift, 0] = 2.0
    rs = hinge_moment(bpts, bpolys, p0, tau_face=ts, semispan=S, q=100.0)
    check("box, EQUAL shear both surfaces -> exactly zero",
          rs["M_h_viscous_Nm"], 0.0, 1e-12, "symmetry")

    # 8. THE TWO HALVES MUST AGREE WHEN THEY CARRY THE SAME FORCE. Give the viscous
    #    branch, face by face, exactly the force the pressure branch produces
    #    (F = -p n dA, so tau = -p n). The moments must then be IDENTICAL. This is the
    #    only control that can catch a sign or projection error in the VISCOUS branch
    #    alone: every other check would pass with the two branches inconsistent.
    cen, Sf = polygon_areas_normals(bpts, bpolys)
    ar = np.linalg.norm(Sf, axis=1)
    pmix = 500.0 + 300.0 * cen[:, 0] + 100.0 * cen[:, 1]
    equiv = -(pmix[:, None] * Sf) / ar[:, None]          # traction with the same force
    ra = hinge_moment(bpts, bpolys, pmix, semispan=S, q=100.0)
    rv = hinge_moment(bpts, bpolys, p0, tau_face=equiv, semispan=S, q=100.0)
    check("pressure and viscous branches agree on identical forces",
          rv["M_h_viscous_Nm"], ra["M_h_pressure_Nm"], 1e-9, "F_p == F_v face by face")

    # 9. THE FRICTION SIGN GATE MUST REFUSE A THRUST. Physics fixes the direction, so
    #    this control constructs the impossible case and passes only if it is rejected.
    try:
        assert_friction_opposes_flow(-tb[:nlift], ar[:nlift], (1, 0, 0))
        refused = False
    except HingeError:
        refused = True
    ok = ok and refused
    print("  [%s] %-46s %s"
          % ("ok " if refused else "FAIL", "a thrust-producing wall shear is REFUSED",
             "raised" if refused else "accepted it"))

    print("\n  self-test: %s" % ("PASS, 9/9 including three that must NOT be zero "
                                 "and one that must be refused"
                                 if ok else "*** FAILURES ABOVE ***"))
    return 0 if ok else 1


# ---------------------------------------------------------------------- vtk ----
# ONE READER, IMPORTED, NOT A SECOND COPY. This module carried its own, and when the
# `FIELD attributes` block turned out to be how OpenFOAM's sampled-surface writer
# actually stores cell data, only the copy in surface_integrate.py was fixed. The other
# then reported "surface carries no cp" on a file that plainly contained it: a reader
# that silently finds nothing rather than saying it cannot read. Duplicated parsers
# diverge on the first fix; the fix is to delete the duplicate, not to patch it too.
from surface_integrate import read_vtk_polydata as _read_vtk    # noqa: E402


def read_vtk_polydata(path):
    try:
        return _read_vtk(path)
    except SystemExit as e:
        raise HingeError(str(e))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--vtk")
    ap.add_argument("--tau-vtk", dest="tau_vtk",
                    help="wallShearStress surface, same patch, for the viscous term")
    ap.add_argument("--no-viscous", dest="no_viscous", action="store_true",
                    help="deliberately omit the viscous term even if the data exists")
    ap.add_argument("--alpha", type=float, default=0.0,
                    help="incidence [deg], for the friction sign gate only")
    ap.add_argument("--q", type=float, help="freestream dynamic pressure [Pa]")
    ap.add_argument("--semispan", type=float, default=SEMISPAN_DEFAULT)
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not a.vtk or a.q is None:
        ap.error("--vtk and --q are required unless --self-test")

    pts, polys, fields = read_vtk_polydata(a.vtk)
    if "cp" not in fields:
        raise HingeError("surface carries no cp; fields present: %s" % sorted(fields))
    p = fields["cp"] * a.q                                  # gauge pressure

    # THE VISCOUS TERM, IF AND ONLY IF THE WALL-SHEAR VECTOR IS THERE. Never inferred
    # from a magnitude: reverse flow is the SIGN of tau and a magnitude cannot carry it.
    tau, tau_src, drag = None, None, None
    # IF THE SHEAR SURFACE IS SITTING NEXT TO THE PRESSURE ONE, USING IT IS NOT OPTIONAL.
    # The viscous term was `null` on every case in this project for a year because nobody
    # sampled the vector. Now that cases do, the remaining way to lose it is to run this
    # script without pointing at it, and get a null that reads exactly like the old one.
    # So an unused sibling is a HALT, not a default: --no-viscous says so out loud.
    if not a.tau_vtk and not a.no_viscous:
        sib = list(Path(a.vtk).parent.glob("wallShearStress*.vtk"))
        if sib:
            raise HingeError(
                "%s sits beside this surface and was not used. The viscous term would be "
                "reported as null while the data to compute it is on disk. Pass "
                "--tau-vtk %s, or --no-viscous if omitting it is deliberate."
                % (sib[0].name, sib[0]))
    if a.tau_vtk:
        tpts, tpolys, tf = read_vtk_polydata(a.tau_vtk)
        if "wallShearStress" not in tf:
            raise HingeError("%s carries no wallShearStress; fields: %s"
                             % (a.tau_vtk, sorted(tf)))
        if len(tpolys) != len(polys):
            raise HingeError("the shear surface has %d faces against the pressure "
                             "surface's %d: these are not the same patch"
                             % (len(tpolys), len(polys)))
        tau = wall_shear_to_body_traction(tf["wallShearStress"])
        _cen, _Sf = polygon_areas_normals(pts, polys)
        drag = assert_friction_opposes_flow(tau, np.linalg.norm(_Sf, axis=1),
                                            [math.cos(math.radians(a.alpha)), 0.0,
                                             math.sin(math.radians(a.alpha))],
                                            label=Path(a.tau_vtk).name)
        tau_src = str(a.tau_vtk)

    res = hinge_moment(pts, polys, p, tau_face=tau, semispan=a.semispan, q=a.q)
    res["_source_vtk"] = str(a.vtk)
    res["_tau_vtk"] = tau_src
    if tau is None:
        res["_viscous_note"] = (
            "VISCOUS TERM OMITTED, NOT ZERO. No wall-shear VECTOR surface was supplied; "
            "cf as a magnitude cannot serve, because the moment needs the direction. "
            "Re-run with argusPostPro installed and pass --tau-vtk. Reported as null "
            "rather than 0.0 so it cannot be mistaken for a computed value.")
    else:
        res["_viscous_note"] = (
            "Viscous term INCLUDED. OpenFOAM's wallShearStress is the NEGATIVE of the "
            "traction on the body (forces.C uses +Sf, wallShearStress.C uses -Sf/|Sf|), "
            "so it is negated on read; measured against the forces object at a ratio of "
            "exactly -1.0000 in all three components. The integrated friction is "
            "%.3f N along the freestream, positive as physics requires." % drag)
        res["viscous_drag_check_N"] = drag
    print(json.dumps({k: v for k, v in res.items() if k != "sectional"}, indent=2))
    print("  sectional stations: %d" % len(res.get("sectional", [])))
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=2) + "\n")
        print("  wrote %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
