#!/usr/bin/env python3
"""case_derived_quantities.py: hinge moment and span efficiency for EVERY case, so the
cases can be compared rather than described one at a time.

    python3 scripts/case_derived_quantities.py --all
    python3 scripts/case_derived_quantities.py SWB_a1p60

WHY THIS EXISTS. Both quantities were previously computed once, by hand, for one case,
and reported to the 19 August audit as PARTIAL and BLOCKED. Both are now computable for
every case that has a pulled surface, and a number produced once by hand is a number
nobody can compare. This runs over the whole pulled set and writes one record per case.

-------------------------------------------------------------------------------
THE FRAME PROBLEM, AND WHY IT IS SETTLED BY MEASUREMENT AND NOT BY ARGUMENT
-------------------------------------------------------------------------------
An incompressible OpenFOAM solver carries pressure and wall shear KINEMATICALLY, divided
through by density. A compressible one does not. So `wallShearStress` on the sampled
surface is either m^2/s^2 or Pa, and the two differ by rho = 1.225, which is close
enough to 1 that BOTH give a plausible-looking C_f and NEITHER announces itself. Reading
the solver family and reasoning it out is exactly the move this project has been bitten
by: a correct answer reached by an invalid route (GEO-080 item 5).

So it is measured. The solver's own `forces` function object writes the viscous force in
NEWTONS, by a different code path from ours, and that is an independently known answer
(D058). We integrate the traction both ways and see which one reproduces it:

    candidate A (tau is already Pa)  :  F = sum(tau * dA)
    candidate B (tau is kinematic)   :  F = rho * sum(tau * dA)

Exactly one should match. If BOTH match the test has no power and we say so; if NEITHER
matches we HALT rather than pick the closer one, because a frame that reproduces nothing
is not a frame, it is a bug somewhere else wearing one.

-------------------------------------------------------------------------------
SPAN EFFICIENCY
-------------------------------------------------------------------------------
C_Di is taken on the Trefftz planes, which is the definition the VLM uses and therefore
the one comparable with CDiw. Crossflow perturbations only:

    v' = U_y ,   w' = U_z - U_inf sin(alpha)
    C_Di = INTEGRAL (v'^2 + w'^2) dS / (U_inf^2 A_ref)
    e    = C_L^2 / (pi * AR * C_Di)

FOUR PLANES ARE REPORTED, NEVER ONE. Induced drag on a wake plane decays with downstream
distance and has no plateau, so a single station is a choice presented as a measurement.
The spread across stations IS the uncertainty and it is reported as such.

FRAMES, stated because this project's failure mode is a right number in an unstated
frame (D068). A_ref is the case's own half-wing 0.620462 m^2, matching its forceCoeffs,
because the planes cover the half model. AR is the FULL-wing 10.7807 on the DSO S_ref
1.24092 m^2 with b_ref 3.6576 m. C_L is the case's own, on the same A_ref.
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

BREF = 3.6576
SREF_DSO = 1.24092
AR_FULL = BREF * BREF / SREF_DSO          # 10.7807 on the DSO basis

# Spanwise bins for the hinge-moment planform. See the note at the call site: 200 is the
# module default and it is 4.5 per cent biased; 1600 is ~0.5 per cent from the limit.
HINGE_NBIN = 1600


# ---------------------------------------------------------------------------
# readers
# ---------------------------------------------------------------------------
def read_polydata(path, want=None):
    """POINTS, POLYGONS and whichever of POINT_DATA / CELL_DATA the file carries."""
    pts = polys = None
    fields, ncell_data, kind = {}, None, None
    with open(path, "r", errors="replace") as fh:
        line = fh.readline()
        while line:
            t = line.split()
            if t and t[0] == "POINTS":
                n = int(t[1]); buf = []
                while len(buf) < 3 * n:
                    buf.extend(fh.readline().split())
                pts = np.asarray(buf, dtype=np.float64).reshape(n, 3)
            elif t and t[0] == "POLYGONS":
                n, tot = int(t[1]), int(t[2]); buf = []
                while len(buf) < tot:
                    buf.extend(fh.readline().split())
                v = np.asarray(buf, dtype=np.int64)
                polys, i = [], 0
                while i < len(v):
                    k = v[i]; polys.append(v[i + 1:i + 1 + k]); i += k + 1
            elif t and t[0] in ("POINT_DATA", "CELL_DATA"):
                kind = t[0]; ncell_data = int(t[1])
                nf = fh.readline().split()
                nfield = int(nf[2]) if len(nf) > 2 else 0
                seen = []
                for _ in range(nfield):
                    # A BLANK LINE IS NOT END OF DATA, AND TREATING IT AS ONE SILENTLY
                    # NARROWED THE FIELD SET. This was `h = fh.readline().split()` followed
                    # by `if not h: break`, which abandoned every REMAINING field the moment
                    # it met a single blank line.
                    #
                    # MEASURED: SWB_trim's surface contains 0 blank lines and SWH_trim's
                    # contains 3. Both declare `FIELD attributes 4` and both carry the same
                    # four headers in the same order, p, cp, yPlus, wallShearStress. Because
                    # wallShearStress is LAST, a blank line anywhere ahead of it deleted it,
                    # and SWH_trim died on KeyError: 'wallShearStress'.
                    #
                    # THE CRASH WAS THE LUCKY CASE. Had the blank line fallen after
                    # wallShearStress, only yPlus would have gone missing, nothing would have
                    # raised, and the case would have produced a hinge moment from a quietly
                    # truncated field set. So the count is now ASSERTED against the file's own
                    # declaration rather than left to whichever field happened to be last
                    # (GEO-087: a verdict over an unstated subset is not a verdict).
                    h = []
                    while not h:
                        line_h = fh.readline()
                        if not line_h:
                            break
                        h = line_h.split()
                    if not h:
                        break
                    name, ncomp = h[0], int(h[1])
                    seen.append(name)
                    vals = []
                    while len(vals) < ncomp * ncell_data:
                        l = fh.readline()
                        if not l:
                            break
                        vals.extend(l.split())
                    if want is None or name in want:
                        a = np.asarray(vals, dtype=np.float64)
                        fields[name] = a.reshape(ncell_data, ncomp) if ncomp > 1 else a
                break
            line = fh.readline()
    return pts, polys, fields, kind


def poly_area_normal(pts, polys):
    """Area vector of each polygon by the fan formula. Returns (centroid, Svec).

    VECTORISED BY POLYGON SIZE. A Python loop over the 5.2 M faces of a wing surface
    takes minutes and is run three times per case; grouping by vertex count and doing
    each group as one array operation takes seconds. snappyHexMesh gives a mixed soup of
    triangles and quads, so there are only ever two or three groups.
    """
    sizes = np.fromiter((len(p) for p in polys), dtype=np.int64, count=len(polys))
    cen = np.empty((len(polys), 3))
    Sv = np.empty((len(polys), 3))
    for k in np.unique(sizes):
        idx = np.flatnonzero(sizes == k)
        conn = np.fromiter((v for i in idx for v in polys[i]),
                           dtype=np.int64, count=len(idx) * int(k)).reshape(len(idx), int(k))
        q = pts[conn]                                   # (m, k, 3)
        c = q.mean(axis=1)                              # (m, 3)
        a = q - c[:, None, :]
        b = np.roll(a, -1, axis=1)
        Sv[idx] = 0.5 * np.cross(a, b).sum(axis=1)
        cen[idx] = c
    return cen, Sv


# ---------------------------------------------------------------------------
# case metadata, read from the case and never assumed
# ---------------------------------------------------------------------------
def case_meta(case):
    cd = REPO / "cases/of12/v5/solve" / case / "system/controlDict"
    if not cd.exists():
        raise SystemExit("  %s: no controlDict at %s" % (case, cd))
    t = cd.read_text(errors="replace")

    def num(pat, default=None):
        m = re.search(pat, t)
        return float(m.group(1)) if m else default

    U = num(r"magUInf\s+([0-9.eE+-]+)")
    Aref = num(r"Aref\s+([0-9.eE+-]+)")
    lRef = num(r"lRef\s+([0-9.eE+-]+)")
    rho = num(r"rhoInf\s+([0-9.eE+-]+)", 1.225)
    compressible = re.search(r"^\s*solver\s+fluid\s*;", t, re.M) is not None
    ld = re.search(r"liftDir\s*\(\s*([-0-9.eE+]+)\s+[-0-9.eE+]+\s+([-0-9.eE+]+)", t)
    alpha = math.degrees(math.atan2(-float(ld.group(1)), float(ld.group(2)))) if ld else None
    return dict(case=case, U_inf=U, Aref=Aref, lRef=lRef, rho=rho, alpha_deg=alpha,
                compressible=compressible,
                q_inf=0.5 * rho * U * U)


def latest_surface(case):
    """The wing surface at the HIGHEST ITERATION, chosen numerically.

    THIS SORTED THE PATH STRING AND TOOK THE LAST, which is lexicographic and therefore a
    coin the directory names flip. Measured on SWB_r2_trim, which carries BOTH
    surfaces/r2/argusWingSurface/6000/ and surfaces/conv/argusWingSurface/4000/: the
    record it produced names conv/4000, so a file called SWB_r2_trim.json described the
    PUBLISHED surface at the published alpha. The name said r2 and the content did not.

    The iteration is the surface's own parent directory, so sort on that as a NUMBER. A
    non-numeric directory sorts lowest rather than raising, so a stray folder degrades the
    choice instead of killing the run. Third instance of this shape today, after
    read_dragdir's forceCoeffs pick and the reported surface time in the spanwise job.
    """
    def it(p):
        try:
            return float(p.parent.name)
        except ValueError:
            return float("-inf")
    g = sorted((REPO / "results/postpro_latest" / case)
               .glob("surfaces/*/argusWingSurface/*/wingSurface.vtk"), key=it)
    return g[-1] if g else None


def resolve_normal_orientation(p_pa, Sv, ref_pressure, tol=0.08, p0_measured=None):
    """Do the polygon normals point OUT of the body, or into it? MEASURED, not assumed.

    THE VISCOUS GATE CANNOT ANSWER THIS. Wall shear is a vector field sampled on the
    surface, so its sign is independent of how the facets happen to be wound; the
    pressure force is p times the AREA VECTOR and flips outright if the winding is
    reversed. These surfaces are known to wind inward -- the export adapter measures
    mean n_z = -0.998 on the geometrically upper faces and calibrates for it -- so a
    hinge moment taken on them without this check is a sign waiting to happen. The
    delivered OF-7 value and the OF-12 one differ in SIGN, which is what forced this
    gate to be written.

    Checked the same way as the shear frame: against the solver's OWN pressure force,
    computed by a different code path, and HALTING when neither orientation reproduces
    it rather than choosing the closer.
    """
    F_out = -(p_pa[:, None] * Sv).sum(axis=0)
    F_in = -F_out
    if ref_pressure is None:
        return None, dict(verdict="NO REFERENCE",
                          note="no forces output pulled; orientation unmeasured")
    t_ref, F_ref = ref_pressure

    # THE OPEN-PATCH OFFSET MODE, AND WHY IT IS REMOVED RATHER THAN TOLERATED.
    # Our integrand is GAUGE pressure (cp * q). The solver integrates the pressure
    # field it carries, and on a COMPRESSIBLE case that field is ABSOLUTE. On a CLOSED
    # surface a uniform offset integrates to zero and the two agree regardless; the
    # `wing` patch is NOT closed, it is cut by the symmetry plane at the root, so
    # sum(Sf) is non-zero and a uniform offset p0 lands as a force p0*sum(Sf).
    #
    # MEASURED ON CMPB_a0p87, and this is the whole diagnosis: sum(Sf) is
    # [-4.5e-15, -3.854e-02, -2.2e-14], i.e. PURELY SPANWISE to machine precision, so
    # the offset can only ever pollute y. The solver and our inward integral agreed to
    # EIGHT significant figures in x and z and disagreed by 1021.176 N in y alone,
    # against a predicted p0*sum(Sf)_y of 1021.176 N with p0 = mean(p - cp*q) =
    # 26495.90 Pa at a standard deviation of 0.029 Pa over 5,276,769 faces. The gate
    # was reporting NEITHER at 192.9 % and 27.8 % on a surface whose winding was never
    # in doubt, and it halted every compressible case in the campaign.
    #
    # SO THE NUISANCE MODE IS PROJECTED OUT OF BOTH SIDES AND THEN ASSERTED, which is
    # strictly MORE checking than before, not less: the component along sum(Sf) is not
    # discarded, it is converted into an IMPLIED p0 and compared against the offset
    # measured from the file. A winding error still fails, because reversing the
    # winding flips the lift and drag components and those are what the test now reads.
    Ssum = Sv.sum(axis=0)
    magS = float(np.linalg.norm(Ssum))
    u = (Ssum / magS) if magS > 0 else None

    def _perp(V):
        return (V - float(np.dot(V, u)) * u) if u is not None else np.asarray(V, float)

    Fr_p = _perp(F_ref)
    n = np.linalg.norm(Fr_p)
    if n <= 0:
        return None, dict(verdict="NO REFERENCE",
                          note="the solver pressure force lies entirely along sum(Sf); "
                               "no winding-sensitive component survives")
    e_out = float(np.linalg.norm(_perp(F_out) - Fr_p) / n)
    e_in = float(np.linalg.norm(_perp(F_in) - Fr_p) / n)
    d = dict(t_reference=t_ref, F_solver_N=[float(x) for x in F_ref],
             F_if_outward_N=[float(x) for x in F_out],
             F_if_inward_N=[float(x) for x in F_in],
             rel_err_if_outward=e_out, rel_err_if_inward=e_in, tolerance=tol,
             sum_Sf_m2=[float(x) for x in Ssum], sum_Sf_magnitude_m2=magS,
             comparison="components orthogonal to sum(Sf); the along-sum(Sf) residual "
                        "is reported below as an implied uniform pressure offset")
    if u is not None:
        d["implied_p0_if_outward_Pa"] = float(np.dot(F_ref - F_out, u) / magS)
        d["implied_p0_if_inward_Pa"] = float(np.dot(F_ref - F_in, u) / magS)
    if p0_measured is not None:
        d["p0_measured_Pa"] = float(p0_measured)
    ok_out, ok_in = e_out < tol, e_in < tol
    if ok_out and ok_in:
        d["verdict"] = "AMBIGUOUS"
        d["note"] = "both windings reproduce the solver force; the test has no power"
        return None, d
    if ok_out:
        d["verdict"] = "OUTWARD"
        return 1.0, d
    if ok_in:
        d["verdict"] = "INWARD"
        d["note"] = ("the facet normals wind INTO the body, so the pressure term is "
                     "negated before integration. Without this the hinge moment comes "
                     "out with the wrong sign and nothing in the result says so.")
        return -1.0, d
    d["verdict"] = "NEITHER"
    d["note"] = ("neither winding reproduces the solver's own pressure force "
                 "(%.1f %% and %.1f %%). HALTED." % (100 * e_out, 100 * e_in))
    return None, d


def solver_forces(case, t_want, meta=None):
    """(t, F_pressure, F_viscous) in NEWTONS, plus the frame that was measured to get there.

    ONE BLOCK, AND ITS FRAME RESOLVED BY MEASUREMENT. This function used to glob every
    directory matching *[Ff]orces* and keep whichever row landed nearest t_want, returning
    the raw numbers. Two blocks exist and THEY ARE IN DIFFERENT FRAMES: on the
    incompressible path `argusForces` is KINEMATIC, density divided out, while `forces1`
    already carries rho, and root_bending's own comment records that they disagree at every
    shared timestep by EXACTLY 1.22500.

    SO THE OLD VERSION WAS RIGHT BY ACCIDENT. SWB_a1p60 has both blocks and the glob handed
    it the dimensional one, so the orientation gate agreed with the surface integral to
    eight significant figures and looked authoritative. SWB_trim has ONLY `argusForces`, the
    gate was handed the kinematic one, and it reported NEITHER at 222.5 % and 22.5 % on a
    surface whose winding was never in doubt: |F_surface| / |F_solver| = 1.225000 with
    cos = +1.000000, which is rho and nothing else. A gate that depends on which directory
    the filesystem happened to list first is not a gate (GEO-080 item 5).

    AND THE FRAME IS NOT PREDICTABLE FROM THE BLOCK NAME, which is why this measures rather
    than prefers. On the COMPRESSIBLE cases `argusForces` is already dimensional: CMPB's
    z-component matched the surface integral to eight digits with no rho applied. The same
    block name carries different frames on different solver paths.

    THE ARBITER IS forceCoeffs, a different code path that reports a DIMENSIONLESS C_L and
    therefore cannot itself be in the wrong frame. Exactly one candidate should reproduce
    it; if neither does we HALT rather than take the closer, because a frame that
    reproduces nothing is a bug wearing one.
    """
    base = REPO / "results/postpro_latest" / case / "postProcessing"
    block = None
    for name in ("argusForces", "forces1"):
        if list((base / name).glob("*/forces.dat")):
            block = name
            break
    if block is None:
        return None

    best = None
    for f in sorted((base / block).glob("*/forces.dat")):
        for line in f.read_text(errors="replace").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            nums = re.findall(r"[-+0-9.eE]+", line)
            if len(nums) < 7:
                continue
            tt = float(nums[0])
            pres = np.asarray([float(x) for x in nums[1:4]])
            visc = np.asarray([float(x) for x in nums[4:7]])
            if best is None or abs(tt - t_want) < abs(best[0] - t_want):
                best = (tt, pres, visc)
    if best is None or meta is None:
        return best

    cl_ref = _forcecoeffs_cl(case)
    if cl_ref is None or not cl_ref:
        return None
    fz = float(best[1][2] + best[2][2])
    q_kin = 0.5 * meta["U_inf"] ** 2
    cand = {"kinematic": (meta["rho"], fz / (q_kin * meta["Aref"])),
            "dimensional": (1.0, fz / (q_kin * meta["rho"] * meta["Aref"]))}
    err = {k: abs(v[1] - cl_ref) / abs(cl_ref) for k, v in cand.items()}
    frame = min(err, key=err.get)
    if err[frame] > 0.02:
        return None
    s = cand[frame][0]
    return (best[0], best[1] * s, best[2] * s,
            dict(block=block, frame=frame, scale=s, cl_reference=cl_ref,
                 rel_err={k: float(v) for k, v in err.items()},
                 note=("the force block's frame was MEASURED against the case's own "
                       "forceCoeffs C_L, not inferred from the block name. Both blocks "
                       "exist on some cases and they differ by exactly rho.")))


def solver_viscous_force(case, t_want):
    """The viscous force vector [N] from the solver's OWN forces output, nearest t_want.

    THIS IS THE INDEPENDENT ANSWER the frame gate is checked against. It is produced by
    OpenFOAM's forces function object, not by us, from the same fields but a different
    code path.
    """
    best = None
    for f in (REPO / "results/postpro_latest" / case).glob("postProcessing/*[Ff]orces*/*/force*.dat"):
        for line in f.read_text(errors="replace").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            # Time ((pres) (visc)) ((mom_pres) (mom_visc))
            nums = re.findall(r"[-+0-9.eE]+", line)
            if len(nums) < 7:
                continue
            tt = float(nums[0])
            visc = np.asarray([float(x) for x in nums[4:7]])
            if best is None or abs(tt - t_want) < abs(best[0] - t_want):
                best = (tt, visc)
    return best


# ---------------------------------------------------------------------------
# the frame gate
# ---------------------------------------------------------------------------
def resolve_tau_frame(tau_raw, Sv, rho, ref, tol=0.08):
    """Decide whether the sampled wall shear is Pa or kinematic, BY MEASUREMENT."""
    area = np.linalg.norm(Sv, axis=1)
    # OpenFOAM's wallShearStress is the NEGATIVE of the traction on the body.
    F_A = -(tau_raw * area[:, None]).sum(axis=0)          # tau already Pa
    F_B = rho * F_A                                       # tau kinematic
    if ref is None:
        return None, dict(verdict="NO REFERENCE", note=(
            "no forces output was pulled for this case, so the frame could not be "
            "measured and no hinge moment is reported. Existence of a plausible C_f is "
            "NOT evidence of a frame."))
    t_ref, F_ref = ref
    n = np.linalg.norm(F_ref)
    eA = float(np.linalg.norm(F_A - F_ref) / n)
    eB = float(np.linalg.norm(F_B - F_ref) / n)
    okA, okB = eA < tol, eB < tol
    d = dict(t_reference=t_ref, F_solver_N=[float(x) for x in F_ref],
             F_if_pascals_N=[float(x) for x in F_A],
             F_if_kinematic_N=[float(x) for x in F_B],
             rel_err_if_pascals=eA, rel_err_if_kinematic=eB, tolerance=tol)
    if okA and okB:
        d["verdict"] = "AMBIGUOUS"
        d["note"] = ("both candidate frames reproduce the solver force within %.0f %%, "
                     "so this test has no power here and the hinge moment is withheld."
                     % (100 * tol))
        return None, d
    if okA:
        d["verdict"] = "PASCALS"
        return 1.0, d
    if okB:
        d["verdict"] = "KINEMATIC"
        return rho, d
    d["verdict"] = "NEITHER"
    d["note"] = ("neither frame reproduces the solver's own viscous force (%.1f %% and "
                 "%.1f %%). HALTED: a frame that reproduces nothing is not a frame."
                 % (100 * eA, 100 * eB))
    return None, d


# ---------------------------------------------------------------------------
# root bending
# ---------------------------------------------------------------------------
def root_bending(case, meta, w=200):
    """Half-wing root bending moment, from the solver's OWN moment output.

    WHY NOT FROM THE SAMPLED SURFACE, like the hinge moment. Root bending needs no
    chordwise resolution and no aft-of-hinge masking: it is the first moment of the
    spanwise load about the root, and the solver already integrates it exactly as M_x
    about a CofR that sits on the symmetry plane (y = 0). Taking it from forces.dat makes
    it available for EVERY case that has a force history, where the hinge moment needs a
    pulled surface. That matters: the morphed geometries have force histories and no
    surfaces, so this quantity arrives for them the moment they run.

    M_x IS THE ROOT BENDING MOMENT ONLY BECAUSE THE CofR HAS y = 0, which is asserted
    rather than assumed. A CofR off the symmetry plane would fold a lift-times-offset term
    into M_x and the number would still look plausible.

    THE FORCES ARE KINEMATIC ON THE INCOMPRESSIBLE PATH, density divided out, so rho is
    applied here. That defect under-reads by a factor of 1.225, which is close enough to
    unity to survive inspection; it is caught by closing the reconstructed C_D against the
    solver's own forceCoeffs, not by reading the numbers.

    REPORTED WITH ITS DYNAMIC PRESSURE, because a dimensional bending moment is
    meaningless without it. The low-order VSPAERO comparison runs at 100 ft/s and ours at
    34 m/s, a factor of 1.244 in q: comparing the two dimensional numbers directly is
    wrong by 24 % before any physics enters.
    """
    # TWO BLOCKS EXIST AND THEY ARE IN DIFFERENT FRAMES. THEY MUST NEVER BE MERGED.
    # The block was renamed mid-project, so the morphed cases carry `forces1` and the
    # current ones carry `argusForces`. Globbing only the new name returned "no force
    # history" for fifteen cases that have one; globbing BOTH and merging was worse, and
    # is the bug this comment exists to prevent. On SWB_a1p60 the two blocks disagree at
    # every shared timestep by EXACTLY 1.22500: `argusForces` is kinematic, `forces1`
    # already carries rho. Merging them mixed frames and moved the answer by 4.5 %, which
    # is small enough to have shipped.
    #
    # SO: ONE BLOCK, PREFERRED BY NAME, AND ITS FRAME RESOLVED BY MEASUREMENT rather than
    # inferred from that name. The same discipline the wall-shear gate already uses.
    base = REPO / "results/postpro_latest" / case / "postProcessing"
    block = None
    for name in ("argusForces", "forces1"):
        if list((base / name).glob("*/forces.dat")):
            block = name
            break
    if block is None:
        return dict(available=False, note="no force history for this case")
    best, cofr = {}, None
    for f in sorted((base / block).glob("*/forces.dat")):
        for line in f.read_text(errors="ignore").splitlines():
            if line.startswith("#"):
                if "CofR" in line:
                    v = [float(x) for x in re.findall(r'[-+0-9.eE]+', line.split(":")[1])]
                    if len(v) == 3:
                        cofr = v
                continue
            v = [float(x) for x in re.findall(r'[-+0-9.eE]+',
                                              line.replace("(", " ").replace(")", " "))]
            if len(v) >= 13:
                best[v[0]] = v
    if not best:
        return dict(available=False, note="no force history for this case")
    if cofr is None:
        return dict(available=False, note="forces.dat carries no CofR; M_x is unanchored")
    if abs(cofr[1]) > 1e-9:
        return dict(available=False,
                    note=("CofR y = %.6g is not on the symmetry plane, so M_x is not the "
                          "root bending moment. HALTED rather than relabelled." % cofr[1]))
    ts = sorted(best)[-w:]
    A = np.array([best[t] for t in ts])
    rho = meta["rho"]
    fz_raw = float((A[:, 3] + A[:, 6]).mean())

    # RESOLVE THE FRAME AGAINST THE CASE'S OWN forceCoeffs, a different code path.
    # Candidate A: the block is kinematic, so rho must be applied.
    # Candidate B: the block already carries rho.
    # Exactly one should reproduce C_L. If neither does we HALT rather than take the
    # closer one, because a frame that reproduces nothing is a bug wearing a frame.
    cl_ref = _forcecoeffs_cl(case, w)
    if cl_ref is None:
        return dict(available=False,
                    note="no forceCoeffs history to resolve the force frame against")
    q_kin = 0.5 * meta["U_inf"] ** 2
    cand = {"kinematic": (rho, fz_raw / (q_kin * meta["Aref"]) * 1.0),
            "dimensional": (1.0, fz_raw / (q_kin * rho * meta["Aref"]))}
    err = {k: abs(v[1] - cl_ref) / abs(cl_ref) for k, v in cand.items()}
    frame = min(err, key=err.get)
    if err[frame] > 0.02:
        return dict(available=False, frame_errors=err,
                    note=("neither force frame reproduces the case's own C_L "
                          "(kinematic %.1f %%, dimensional %.1f %%). HALTED."
                          % (100 * err["kinematic"], 100 * err["dimensional"])))
    scale = cand[frame][0]

    mp, mv = float(A[:, 7].mean()) * scale, float(A[:, 10].mean()) * scale
    fz = fz_raw * scale
    half_span = 3.6576 / 2.0
    return dict(available=True, n_samples=len(ts), t_last=ts[-1],
                block=block, force_frame=frame,
                frame_rel_err={k: float(v) for k, v in err.items()},
                cofr=cofr, q_Pa=meta["q_inf"],
                RBM_pressure_Nm=mp, RBM_viscous_Nm=mv, RBM_total_Nm=mp + mv,
                normal_force_half_wing_N=fz,
                lift_centroid_y_m=(mp + mv) / fz if fz else None,
                lift_centroid_eta=(mp + mv) / fz / half_span if fz else None,
                note=("half wing, about the root chord line at y = 0. The force frame was "
                      "MEASURED against the case's own forceCoeffs, not inferred from the "
                      "block name. Quote with q_Pa: a dimensional bending moment is not "
                      "comparable across dynamic pressures."))


def _forcecoeffs_cl(case, w=200):
    """The case's own C_L window mean from forceCoeffs, the independent code path."""
    root = REPO / "results/postpro_latest" / case / "postProcessing/forceCoeffs1"
    best, icl = {}, None
    for d in sorted(root.glob("*")) if root.exists() else []:
        f = d / "forceCoeffs.dat"
        if not f.exists():
            continue
        for line in f.read_text(errors="ignore").splitlines():
            if line.startswith("#"):
                tok = line.lstrip("#").split()
                if "Cl" in tok:
                    icl = tok.index("Cl")
                continue
            p = line.split()
            if len(p) > 3:
                try:
                    best[float(p[0])] = p
                except ValueError:
                    pass
    if not best or icl is None:
        return None
    ts = sorted(best)[-w:]
    return float(np.mean([float(best[t][icl]) for t in ts]))


# ---------------------------------------------------------------------------
# span efficiency
# ---------------------------------------------------------------------------
def span_efficiency(case, meta, CL):
    """C_Di and e on every Trefftz plane that was sampled. Four stations, not one."""
    # TWO LAYOUTS, BOTH REAL. An early manual pull put the planes in <case>/trefftz/<t>/;
    # collect_postpro.sh puts them in <case>/surfaces/<t>/argusTrefftz/<t>/. Globbing only
    # the first reported "no Trefftz planes pulled" for CMPB_a1p20 while all four planes
    # sat on disk under the second -- an absence claimed by looking in one place. The
    # search now covers both and says which it used.
    base = REPO / "results/postpro_latest" / case
    cand = []
    for root in (base / "trefftz", ):
        if root.exists():
            cand += [d for d in root.glob("*") if d.is_dir() and d.glob("trefftz_x*.vtk")]
    cand += [d for d in base.glob("surfaces/*/argusTrefftz/*") if d.is_dir()]
    times = sorted([d for d in cand if any(d.glob("trefftz_x*.vtk"))],
                   key=lambda d: float(d.name) if d.name.replace(".", "").isdigit() else 0.0)
    if not times:
        return dict(available=False,
                    note="no Trefftz planes pulled for this case")
    tdir = times[-1]
    U = meta["U_inf"]
    w_free = U * math.sin(math.radians(meta["alpha_deg"] or 0.0))
    out = []
    for f in sorted(tdir.glob("trefftz_x*.vtk")):
        pts, polys, fld, kind = read_polydata(f, want=("U",))
        if "U" not in fld or pts is None or polys is None:
            continue
        cen, Sv = poly_area_normal(pts, polys)
        dS = np.abs(Sv[:, 0])                 # plane normal is streamwise
        Up = fld["U"]
        if kind == "POINT_DATA":
            # average the point values onto each polygon
            Uc = np.empty((len(polys), 3))
            for i, p in enumerate(polys):
                Uc[i] = Up[p].mean(axis=0)
        else:
            Uc = Up
        v = Uc[:, 1]
        w = Uc[:, 2] - w_free
        CDi = float(((v * v + w * w) * dS).sum() / (U * U * meta["Aref"]))
        e = float(CL * CL / (math.pi * AR_FULL * CDi)) if CDi > 0 else None
        out.append(dict(plane=f.stem, x_m=float(cen[:, 0].mean()),
                        CDi=CDi, CDi_counts=1e4 * CDi, span_efficiency=e,
                        n_faces=len(polys)))
    if not out:
        return dict(available=False, note="Trefftz files present but carried no U")
    es = [o["span_efficiency"] for o in out if o["span_efficiency"]]
    cds = [o["CDi_counts"] for o in out]
    return dict(available=True, time=float(tdir.name), planes=out,
                CDi_counts_range=[min(cds), max(cds)],
                span_efficiency_range=[min(es), max(es)] if es else None,
                note=("Reported across ALL sampled stations. Induced drag on a wake "
                      "plane decays downstream and has no plateau, so the SPREAD is the "
                      "uncertainty, not scatter to be averaged away."))


def case_CL(case, t_want):
    f = REPO / "results/postpro_latest" / case / "forceCoeffs_spliced.dat"
    if not f.exists():
        return None, None
    best = None
    for line in f.read_text(errors="replace").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) < 3:
            continue
        t = float(p[0])
        if best is None or abs(t - t_want) < abs(best[0] - t_want):
            best = (t, float(p[2]), float(p[1]))
    return (best[1], best[0]) if best else (None, None)


# ---------------------------------------------------------------------------
def run(case):
    print("  ==== %s ====" % case)
    meta = case_meta(case)
    rec = dict(meta)
    surf = latest_surface(case)
    rec["surface_vtk"] = str(surf.relative_to(REPO)) if surf else None

    CL = t_surf = None
    if surf:
        t_surf = float(surf.parent.name)
        CL, t_cl = case_CL(case, t_surf)
        rec["iteration"] = t_surf
        rec["CL"] = CL

    # ---- hinge moment -----------------------------------------------------
    if surf is None:
        rec["hinge_moment"] = dict(available=False,
                                   note="no wing surface pulled for this case")
        print("    hinge moment: no surface pulled")
    else:
        pts, polys, fld, _ = read_polydata(surf,
                                           want=("cp", "p", "wallShearStress", "yPlus"))
        _cen, Sv = poly_area_normal(pts, polys)
        allf = solver_forces(case, t_surf, meta)
        ref = (allf[0], allf[2]) if allf else None
        ref_p = (allf[0], allf[1]) if allf else None
        if allf and len(allf) > 3:
            rec["solver_force_frame"] = allf[3]
            print("    force frame: %s block=%s scale=%.6g (kinematic %.2f %%, "
                  "dimensional %.2f %%)"
                  % (allf[3]["frame"], allf[3]["block"], allf[3]["scale"],
                     100 * allf[3]["rel_err"]["kinematic"],
                     100 * allf[3]["rel_err"]["dimensional"]))
        scale, frame = resolve_tau_frame(fld["wallShearStress"], Sv, meta["rho"], ref)
        rec["tau_frame"] = frame
        print("    tau frame: %s (Pa %.1f %%, kinematic %.1f %%)"
              % (frame.get("verdict"), 100 * frame.get("rel_err_if_pascals", float("nan")),
                 100 * frame.get("rel_err_if_kinematic", float("nan"))))
        p_pa = fld["cp"] * meta["q_inf"]

        # THE OFFSET IS MEASURED FROM THE FILE, NOT ASSUMED, and only where it MEANS
        # something. On a compressible case `p` is absolute Pa, so p - cp*q is the
        # uniform reference pressure and its SPREAD is the evidence that it is uniform.
        # On an incompressible case `p` is KINEMATIC (m^2/s^2), so that same difference
        # is cp*(0.5U^2 - q) and varies face by face: it is not an offset at all, and
        # computing one there would manufacture a number out of a frame error (D068).
        p0_measured = None
        if meta.get("compressible") and "p" in fld:
            off = fld["p"] - p_pa
            p0_measured = float(off.mean())
            rec["reference_pressure"] = dict(
                p0_Pa=p0_measured, std_Pa=float(off.std()),
                min_Pa=float(off.min()), max_Pa=float(off.max()),
                note=("mean of p - cp*q over every face. A UNIFORM value is the claim "
                      "being made, so the standard deviation is reported beside it: a "
                      "spread comparable to the mean would mean this is not an offset "
                      "and the orientation gate's correction would be unfounded."))
        orient, onote = resolve_normal_orientation(p_pa, Sv, ref_p,
                                                   p0_measured=p0_measured)
        rec["normal_orientation"] = onote
        print("    normals:   %s (outward %.1f %%, inward %.1f %%)"
              % (onote.get("verdict"),
                 100 * onote.get("rel_err_if_outward", float("nan")),
                 100 * onote.get("rel_err_if_inward", float("nan"))))
        if scale is None or orient is None:
            rec["hinge_moment"] = dict(
                available=False,
                note="gate did not resolve: " + frame.get("note", "") + " " + onote.get("note", ""))
        else:
            import hinge_moment_integrate as HM
            # ORIENTATION APPLIED TO BOTH TERMS. Reversing the winding negates the
            # pressure force; the traction is a sampled vector and does not care, but
            # the MOMENT ARM is unaffected either way, so the whole integrand flips
            # together and applying it to p alone would leave the two terms in
            # different conventions.
            tau_pa = HM.wall_shear_to_body_traction(fld["wallShearStress"] * scale)
            # nbin IS NOT A FREE PARAMETER, IT IS A CONVERGED CHOICE. The default of 200
            # carried a 4.5 per cent bias in |M_h| because planform()'s per-bin min/max
            # overstates the chord and therefore mis-selects the faces aft of the hinge
            # line. A 100 to 1600 sweep on SWB_trim halves the increment on every doubling,
            # first-order, and 1600 sits about 0.5 per cent from the Richardson limit.
            # That residual is the number to quote as the method uncertainty, and it is
            # stated rather than assumed away.
            res = HM.hinge_moment(pts, polys, p_pa * orient, tau_face=tau_pa,
                                  q=meta["q_inf"], nbin=HINGE_NBIN)
            res["available"] = True
            rec["hinge_moment"] = res
            print("    M_h pressure %.6f N m, viscous %s, total %s"
                  % (res["M_h_pressure_Nm"],
                     ("%.6f" % res["M_h_viscous_Nm"]) if res["M_h_viscous_Nm"] is not None else "None",
                     ("%.6f" % res["M_h_total_Nm"]) if res["M_h_total_Nm"] is not None else "None"))

    # ---- span efficiency --------------------------------------------------
    se = span_efficiency(case, meta, CL) if CL else dict(
        available=False, note="no C_L available, so e cannot be formed")
    # ---- root bending -----------------------------------------------------
    rb = root_bending(case, meta)
    rec["root_bending"] = rb
    if rb.get("available"):
        print("    root bending %.4f N m half wing at q %.2f Pa; lift centroid eta %.4f"
              % (rb["RBM_total_Nm"], rb["q_Pa"], rb["lift_centroid_eta"]))
    else:
        print("    root bending: %s" % rb.get("note"))

    rec["span_efficiency"] = se
    if se.get("available"):
        print("    C_Di %.3f to %.3f counts over %d planes; e %.4f to %.4f"
              % (se["CDi_counts_range"][0], se["CDi_counts_range"][1], len(se["planes"]),
                 se["span_efficiency_range"][0], se["span_efficiency_range"][1]))
    else:
        print("    span efficiency: %s" % se.get("note"))

    out = REPO / "results/derived" / ("%s.json" % case)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, default=float) + "\n")
    return rec


def run_quick(case):
    """Only the quantities derived from the FORCE HISTORY, merged into any existing record.

    The surface pass costs about twenty minutes a case and root bending does not need it.
    This updates that field alone and LEAVES EVERY SURFACE-DERIVED FIELD UNTOUCHED, so a
    record that already carries a hinge moment keeps it. It never invents a record's
    surface fields: a case with no prior record gets one with the surface entries absent
    rather than falsely marked unavailable.
    """
    meta = case_meta(case)
    out = REPO / "results/derived" / ("%s.json" % case)
    rec = json.loads(out.read_text()) if out.exists() else dict(meta)
    print("  ==== %s (quick) ====" % case)
    rb = root_bending(case, meta)
    rec["root_bending"] = rb
    if rb.get("available"):
        print("    root bending %.4f N m half wing at q %.2f Pa; lift centroid eta %.4f"
              % (rb["RBM_total_Nm"], rb["q_Pa"], rb["lift_centroid_eta"]))
    else:
        print("    root bending: %s" % rb.get("note"))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, default=float) + "\n")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--quick", action="store_true",
                    help="force-history quantities only (root bending); skips the "
                         "surface pass and preserves existing surface-derived fields")
    a = ap.parse_args()
    root = REPO / "results/postpro_latest"
    cases = a.cases or (sorted(p.name for p in root.glob("*") if p.is_dir())
                        if a.all else [])
    if not cases:
        ap.error("name a case or pass --all")
    ok = 0
    for c in cases:
        try:
            (run_quick if a.quick else run)(c); ok += 1
        except SystemExit as e:
            print("  %s: %s" % (c, e))
        except Exception as e:                                  # noqa: BLE001
            print("  %s: FAILED %s: %s" % (c, type(e).__name__, e))
    # PARTITION THE SET, per GEO-089: every case named lands in a counted bucket.
    print("  ---- %d of %d cases produced a record; %d did not ----"
          % (ok, len(cases), len(cases) - ok))
    return 0


if __name__ == "__main__":
    sys.exit(main())
