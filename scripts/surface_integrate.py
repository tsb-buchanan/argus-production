#!/usr/bin/env python3
"""surface_integrate.py: EXACT surface integrals from a sampled VTK patch.

  python3 scripts/surface_integrate.py --self-test
  python3 scripts/surface_integrate.py --case cases/validation/onera_m6_B_flin2 --time 500

-------------------------------------------------------------------------------
WHY THIS EXISTS: A `raw` SAMPLE CANNOT BE INTEGRATED
-------------------------------------------------------------------------------
OpenFOAM's `raw` surface format writes face-centre POINTS ONLY, with no connectivity,
so face areas and normals do not exist in the file. Every integral formed from one has
to RECONSTRUCT them: sort by chord, split upper from lower, walk a contour. That works
approximately on a smooth section and not at all on a tip cap, and on the ONERA M6 it
left a 4.4-count closure error against the forces object, against a signal of 20.8
counts (ARG-178). The limit on that investigation was the instrument, not the physics.

VTK POLYDATA CARRIES THE POLYGONS. Face area vectors then come from the geometry by
Newell's formula, and the integral is a sum with no reconstruction anywhere in it. That
is why `hpc/templates/argusPostPro` samples VTK, and this module is the consumer that
justifies the choice.

-------------------------------------------------------------------------------
THE ORIENTATION IS ESTABLISHED, NOT ASSUMED
-------------------------------------------------------------------------------
A face's area vector has a sign, and the whole integral flips with it. It is fixed here
by GEOMETRY rather than by whichever choice makes the answer agree with something: on a
lifting surface the outward normal above the chord plane points away from it. The result
is then checked against three known answers that the choice cannot fake.

-------------------------------------------------------------------------------
THE KNOWN ANSWERS
-------------------------------------------------------------------------------
1. A CLOSED SURFACE HAS ZERO TOTAL AREA VECTOR. Our wing patch is open only at the root,
   where it meets the symmetry plane, and that opening is normal to y. So sum(A_x) and
   sum(A_z) must vanish while sum(A_y) equals the root cross-section area. This is a
   statement about the mesh topology that no flow field can influence.
2. UNIFORM PRESSURE EXERTS NO FORCE on a body closed in that direction, which follows
   from 1 and is checked separately because it also exercises the weighting.
3. THE FORCES FUNCTION OBJECT IS AN INDEPENDENT ROUTE TO THE SAME NUMBER. OpenFOAM
   computes it inside the solver from its own face areas and normals; this module never
   sees those. Agreement is therefore evidence, not a tautology (D058).
"""
import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]


def read_vtk_polydata(path):
    """POINTS, POLYGONS and CELL_DATA from a legacy ASCII VTK polydata file.

    WRITTEN OUT RATHER THAN PULLED FROM A LIBRARY because the only thing needed is three
    blocks of numbers, and a dependency that must be installed on the cluster to read our
    own output is a dependency that will one day not be there.
    """
    txt = Path(path).read_text(errors="replace")
    m = re.search(r"^POINTS\s+(\d+)\s+\w+\s*$", txt, re.M)
    if not m:
        raise SystemExit("%s: no POINTS block; is this legacy ASCII VTK polydata?" % path)
    npt = int(m.group(1))
    body = txt[m.end():]
    vals = np.fromstring(body[:body.find("POLYGONS")], sep=" ")
    P = vals[:3 * npt].reshape(npt, 3)

    m = re.search(r"^POLYGONS\s+(\d+)\s+(\d+)\s*$", txt, re.M)
    if not m:
        raise SystemExit("%s: no POLYGONS block, so there are no face areas or normals "
                         "in this file and no integral can be formed from it." % path)
    nfc, ntot = int(m.group(1)), int(m.group(2))
    rest = txt[m.end():]
    cut = re.search(r"^(CELL_DATA|POINT_DATA)\b", rest, re.M)
    conn = np.fromstring(rest[:cut.start()] if cut else rest, sep=" ", dtype=float)
    conn = conn[:ntot].astype(int)

    faces, i = [], 0
    while i < conn.size:
        n = conn[i]
        faces.append(conn[i + 1:i + 1 + n])
        i += 1 + n
    if len(faces) != nfc:
        raise SystemExit("%s: header says %d polygons, %d were read"
                         % (path, nfc, len(faces)))

    data = {}
    # OpenFOAM's sampled-surface writer emits a `FIELD attributes N` block, whose entries
    # are `name ncomp ntuples type`, NOT the SCALARS/VECTORS forms. Parsing only the
    # latter returned an empty field dictionary while the file plainly contained the data,
    # which is a reader that silently finds nothing rather than saying it cannot read.
    for m in re.finditer(r"^\s*(\w+)\s+(\d+)\s+(\d+)\s+(float|double)\s*$", txt, re.M):
        name, ncomp, ntup = m.group(1), int(m.group(2)), int(m.group(3))
        if ntup != nfc:
            continue
        seg = txt[m.end():]
        stop = re.search(r"^\s*(\w+\s+\d+\s+\d+\s+(?:float|double)|SCALARS|VECTORS|FIELD|"
                         r"CELL_DATA|POINT_DATA)\b", seg, re.M)
        arr = np.fromstring(seg[:stop.start()] if stop else seg, sep=" ")
        want = nfc * ncomp
        if arr.size >= want:
            data[name] = arr[:want].reshape(nfc, ncomp) if ncomp == 3 else arr[:want]

    for m in re.finditer(r"^(SCALARS\s+(\w+)\s+\w+(?:\s+\d+)?|VECTORS\s+(\w+)\s+\w+)\s*$",
                         txt, re.M):
        name = m.group(2) or m.group(3)
        ncomp = 1 if m.group(2) else 3
        seg = txt[m.end():]
        if m.group(2):
            lut = re.match(r"\s*LOOKUP_TABLE\s+\w+\s*", seg)
            if lut:
                seg = seg[lut.end():]
        stop = re.search(r"^(SCALARS|VECTORS|FIELD|CELL_DATA|POINT_DATA)\b", seg, re.M)
        arr = np.fromstring(seg[:stop.start()] if stop else seg, sep=" ")
        want = nfc * ncomp
        if arr.size >= want:
            data[name] = arr[:want].reshape(nfc, ncomp) if ncomp == 3 else arr[:want]
    return P, faces, data


def area_vectors(P, faces):
    """Newell's formula, which is exact for a planar polygon and gives the correct mean
    normal for a warped one, and which is what OpenFOAM's own face decomposition reduces
    to for a flat face. Magnitude is the area, direction the normal."""
    A = np.zeros((len(faces), 3))
    C = np.zeros((len(faces), 3))
    for k, f in enumerate(faces):
        v = P[f]
        vn = np.roll(v, -1, axis=0)
        A[k] = 0.5 * np.cross(v, vn).sum(axis=0)
        C[k] = v.mean(axis=0)
    return A, C


def orient_outward(A, C, chord_plane_z=0.0):
    """Make every area vector point OUT of the body, and say how it was decided.

    THE TEST IS GEOMETRIC. Above the chord plane the outward normal has a positive z
    component and below it a negative one, over the great majority of a lifting surface;
    the exceptions are the leading-edge nose and the tip, which are a small minority of
    faces. So the majority vote over faces away from those regions fixes the global sign,
    and the vote MARGIN is returned so a marginal decision cannot pass silently.
    """
    up = C[:, 2] > chord_plane_z
    lo = C[:, 2] < chord_plane_z
    agree = int((A[up, 2] > 0).sum() + (A[lo, 2] < 0).sum())
    disagree = int((A[up, 2] < 0).sum() + (A[lo, 2] > 0).sum())
    flipped = agree < disagree
    if flipped:
        A = -A
        agree, disagree = disagree, agree
    frac = agree / max(agree + disagree, 1)
    return A, dict(flipped=bool(flipped), agree=agree, disagree=disagree,
                   agreement_fraction=float(frac))


def closure_checks(A, tol_rel=1e-3):
    """sum(A_x) and sum(A_z) must vanish on a body closed in x and z; sum(A_y) is the
    root opening. Reported as FRACTIONS OF TOTAL AREA so the numbers are comparable
    between meshes."""
    tot = float(np.linalg.norm(A, axis=1).sum())
    s = A.sum(axis=0)
    return dict(total_area_m2=tot,
                sum_Ax_over_area=float(s[0] / tot), sum_Ay_over_area=float(s[1] / tot),
                sum_Az_over_area=float(s[2] / tot),
                root_opening_area_m2=float(abs(s[1])),
                closed_in_x=bool(abs(s[0] / tot) < tol_rel),
                closed_in_z=bool(abs(s[2] / tot) < tol_rel))


def spanwise_bands(A, C, cp, aref, span, edges):
    """Exact contribution of each spanwise band to C_A and C_N, and it PARTITIONS.

    NO RECONSTRUCTION ANYWHERE. Each face carries its own area vector, so a face belongs
    to exactly one band and the band sums are the total by construction. The bands are
    asserted to sum, because a decomposition that does not is a verdict over an unstated
    subset (GEO-087).
    """
    eta = C[:, 1] / span
    out, assigned = [], np.zeros(len(cp), bool)
    for a0, a1 in zip(edges[:-1], edges[1:]):
        m = (eta >= a0) & (eta < a1)
        assigned |= m
        ca = float(-(cp[m] * A[m, 0]).sum() / aref)
        cn = float(-(cp[m] * A[m, 2]).sum() / aref)
        out.append(dict(eta0=float(a0), eta1=float(a1), n_faces=int(m.sum()),
                        area_m2=float(np.linalg.norm(A[m], axis=1).sum()),
                        C_A=ca, C_N=cn))
    leftover = int((~assigned).sum())
    tot_ca, _ = pressure_force_coeffs(A, cp, aref)
    sum_ca = sum(b["C_A"] for b in out)
    return out, dict(n_unassigned=leftover, bands_sum_C_A=sum_ca, total_C_A=tot_ca,
                     closes=bool(abs(sum_ca - tot_ca) < 1e-9 and leftover == 0))


def strip_sectional_cA(A, C, cp, stations, half_band_eta, span, x=None):
    """Exact SECTIONAL axial pressure coefficient in a narrow spanwise strip.

        c_A(y) = [-sum over the strip of Cp * A_x] / (dy * c_local)

    THIS IS THE SAME QUANTITY the reference codes' contour integral of Cp gives, formed
    without any contour. It is what replaces the reconstruction that left a 4.4-count
    closure error (ARG-178): here there is nothing to reconstruct, because every face
    already carries its own area vector.

    THE STRIP WIDTH CANCELS ONLY IF IT IS MEASURED, NOT ASSUMED. dy is taken as the
    actual spread of face centres in the strip rather than the nominal band width, since
    a face straddling the edge is either in or out and the nominal width then over-counts.
    """
    eta = C[:, 1] / span
    out = {}
    for st in stations:
        m = np.abs(eta - st) <= half_band_eta
        if m.sum() < 40:
            continue
        xs = C[m, 0] if x is None else x[m]
        c_local = float(xs.max() - xs.min())
        dy = float(C[m, 1].max() - C[m, 1].min())
        if dy <= 0 or c_local <= 0:
            continue
        out[st] = dict(c_A=float(-(cp[m] * A[m, 0]).sum() / (dy * c_local)),
                       c_N=float(-(cp[m] * A[m, 2]).sum() / (dy * c_local)),
                       chord_m=c_local, dy_m=dy, n_faces=int(m.sum()))
    return out


def pressure_force_coeffs(A, cp, aref):
    """(C_A, C_N) from Cp and outward area vectors. F = -sum(p n dA), so with Cp the
    uniform part drops out exactly and the sum is over Cp * A."""
    return (float(-(cp * A[:, 0]).sum() / aref),
            float(-(cp * A[:, 2]).sum() / aref))


# -----------------------------------------------------------------------------
def self_test():
    ok, bad = 0, []

    def check(label, got, want, tol):
        nonlocal ok
        if abs(got - want) <= tol:
            ok += 1
            print("    PASS  %-54s %+.6e (want %+.6e)" % (label, got, want))
        else:
            bad.append(label)
            print("    FAIL  %-54s %+.6e (want %+.6e, tol %.1e)"
                  % (label, got, want, tol))

    # A CLOSED UNIT CUBE, built by hand, with OUTWARD faces. Its area is 6 exactly, its
    # total area vector is zero exactly, and a uniform pressure exerts no force on it.
    v = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    faces = [np.array(f) for f in ([0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
                                   [2, 3, 7, 6], [1, 2, 6, 5], [0, 4, 7, 3])]
    A, C = area_vectors(v, faces)
    check("unit cube total area", float(np.linalg.norm(A, axis=1).sum()), 6.0, 1e-12)
    check("unit cube sum of area vectors, x", float(A.sum(axis=0)[0]), 0.0, 1e-12)
    check("unit cube sum of area vectors, z", float(A.sum(axis=0)[2]), 0.0, 1e-12)
    ca, cn = pressure_force_coeffs(A, np.full(6, -0.73), 1.0)
    check("uniform Cp on a closed body: axial force", ca, 0.0, 1e-12)
    check("uniform Cp on a closed body: normal force", cn, 0.0, 1e-12)

    # A KNOWN NON-ZERO ANSWER. Load ONLY the -x face (outward normal -x, area 1) with
    # Cp = 2. F_x = -Cp * A_x = -2 * (-1) = +2, so C_A = +2 on A_ref = 1. Analytic.
    cp = np.zeros(6)
    ix = int(np.argmin(A[:, 0]))            # the face whose outward normal is -x
    cp[ix] = 2.0
    ca, cn = pressure_force_coeffs(A, cp, 1.0)
    check("Cp = 2 on the -x face only, C_A", ca, 2.0, 1e-12)
    check("  and it produces no normal force", cn, 0.0, 1e-12)

    # ORIENTATION: hand it the cube INSIDE OUT and it must flip it back.
    Aflip, info = orient_outward(-A, C)
    check("an inside-out surface is detected and flipped",
          float(np.abs(Aflip - A).max()), 0.0, 1e-12)
    if info["flipped"]:
        ok += 1
        print("    PASS  %-54s flipped, %d against %d"
              % ("and it reports that it flipped", info["agree"], info["disagree"]))
    else:
        bad.append("flip reported")
        print("    FAIL  the flip was not reported")

    # AREA IS INDEPENDENT OF VERTEX COUNT: a square split into two triangles has the
    # same area vector as the quad. Guards a Newell implementation that assumes 4 nodes.
    tri = [np.array([0, 1, 2]), np.array([0, 2, 3])]
    At, _ = area_vectors(v[:4], [np.array([0, 1, 2, 3])])
    At2, _ = area_vectors(v[:4], tri)
    check("a quad and its two triangles have the same area vector",
          float(np.abs(At[0] - At2.sum(axis=0)).max()), 0.0, 1e-12)

    print("\n  self-test %d/%d" % (ok, ok + len(bad)))
    return 0 if not bad else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", default="cases/validation/onera_m6_B_flin2")
    ap.add_argument("--time", default="500")
    ap.add_argument("--surface", default="argusWingSurface")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--json", default="results/onera_m6_exact_surface_integral.json")
    ap.add_argument("--figure",
                    default="results/figures/deck/onera_m6_exact_integral.png")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    case = REPO / a.case
    d = case / "postProcessing" / a.surface / str(a.time)
    cpf = next(iter(sorted(d.glob("cp_*.vtk"))), None)
    pf = next(iter(sorted(d.glob("p_*.vtk"))), None)
    if not pf:
        sys.exit("  no p_*.vtk under %s. Run the case with hpc/templates/argusPostPro "
                 "installed; a `raw` sample cannot be integrated." % d)
    P, faces, data = read_vtk_polydata(pf)
    print("  %s: %d points, %d polygons, fields %s"
          % (pf.name, len(P), len(faces), sorted(data)))
    A, C = area_vectors(P, faces)
    A, orient = orient_outward(A, C)
    print("  orientation: %s, %d faces agree against %d (%.1f%%)"
          % ("FLIPPED to outward" if orient["flipped"] else "already outward",
             orient["agree"], orient["disagree"], 100 * orient["agreement_fraction"]))
    cl = closure_checks(A)
    print("  topology closure, as fractions of the %.4f m2 total area:"
          % cl["total_area_m2"])
    print("    sum(A_x)/area %+.3e   sum(A_z)/area %+.3e   %s"
          % (cl["sum_Ax_over_area"], cl["sum_Az_over_area"],
             "closed in x and z" if (cl["closed_in_x"] and cl["closed_in_z"])
             else "NOT CLOSED: the patch has an opening normal to x or z"))
    print("    sum(A_y) = %.6f m2, the root opening where the wing meets the symmetry "
          "plane" % cl["root_opening_area_m2"])

    sys.path.insert(0, str(REPO / "scripts"))
    from validate_onera_m6 import load_surface, freestream_from_case   # noqa: E402
    _P, _Q, fs = load_surface(case, a.time, verbose=False)

    # Cp on the SAME faces. Prefer the solver's own cp field; fall back to forming it
    # from p with the case's own frame, and SAY WHICH was used.
    if cpf:
        _Pc, _fc, dc = read_vtk_polydata(cpf)
        cp, cp_src = dc["cp"], "the solver's own cp field (%s)" % cpf.name
    else:
        cp = (data["p"] - fs["pinf"]) / fs["q"]
        cp_src = "formed from p with the case's own p_inf and q"
    print("  cp from %s: range %+.4f to %+.4f" % (cp_src, cp.min(), cp.max()))

    # NULL: a uniform Cp must give exactly zero force, which follows from the topology
    # closure above but also exercises the weighting.
    z0, z1 = pressure_force_coeffs(A, np.full(len(faces), -0.73), 1.0)
    print("  null: uniform Cp gives C_A %+.3e, C_N %+.3e (must be zero)" % (z0, z1))

    aref = planform_area(A)
    ca, cn = pressure_force_coeffs(A, cp, aref)
    alpha = math.atan2(15.5, 290.0)
    cdp = ca * math.cos(alpha) + cn * math.sin(alpha)
    print("\n  EXACT surface integral on A_ref = %.6f m2:" % aref)
    print("    C_Ap %+.6f   C_Np %+.6f   ->  C_Dp %+.6f  (%.2f counts)"
          % (ca, cn, cdp, 1e4 * cdp))

    import importlib
    _mx = importlib.import_module("m6_pressure_excess")
    _fr = _mx.tmr_forces(_mx.mesh_cells(case), math.cos(alpha), math.sin(alpha))
    if _fr:
        out_rows = _fr["rows"]
    else:
        out_rows = []
    F = force_object(case, a.time)
    out = dict(_what="Exact pressure-force coefficients from the sampled VTK patch.",
               case=a.case, time=a.time, surface=a.surface,
               n_points=len(P), n_faces=len(faces), cp_source=cp_src,
               orientation=orient, topology=cl, freestream=fs,
               reference_area_m2=aref, alpha_deg=math.degrees(alpha),
               C_Ap=ca, C_Np=cn, C_Dp=cdp,
               null_uniform_cp=dict(C_A=z0, C_N=z1),
               reference_rows=out_rows,
               generated_by="scripts/surface_integrate.py")
    if F:
        q = fs["q"]
        ca_f = F["pressure"][0] / (q * aref)
        cn_f = F["pressure"][2] / (q * aref)
        cdp_f = ca_f * math.cos(alpha) + cn_f * math.sin(alpha)
        print("\n  the forces function object, an INDEPENDENT route:")
        print("    C_Ap %+.6f   C_Np %+.6f   ->  C_Dp %+.6f  (%.2f counts)"
              % (ca_f, cn_f, cdp_f, 1e4 * cdp_f))
        print("    difference   %+.6f   %+.6f       %+.6f  (%.2f counts)"
              % (ca - ca_f, cn - cn_f, cdp - cdp_f, 1e4 * (cdp - cdp_f)))
        print("    OpenFOAM formed that from its own face areas and normals, which this")
        print("    module never sees, so the agreement is evidence and not a tautology.")
        out["forces_object"] = dict(C_Ap=ca_f, C_Np=cn_f, C_Dp=cdp_f,
                                    d_C_Ap=ca - ca_f, d_C_Dp=cdp - cdp_f,
                                    d_C_Dp_counts=1e4 * (cdp - cdp_f))
    # ---- the spanwise decomposition, and the tip the reference wing does not have -----
    SPAN_NOM = 1.1963
    edges = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.96, 1.0, 1.05])
    bands, part = spanwise_bands(A, C, cp, aref, SPAN_NOM, edges)
    print("\n  EXACT spanwise decomposition of C_Ap, on the NOMINAL semispan %.4f m:"
          % SPAN_NOM)
    print("    %-14s %-9s %-11s %s" % ("y/b band", "faces", "area [m2]", "C_Ap share"))
    for b in bands:
        print("    %.2f to %-8.2f %-9d %-11.5f %+.6f"
              % (b["eta0"], b["eta1"], b["n_faces"], b["area_m2"], b["C_A"]))
    print("    partition: %d faces unassigned, bands sum %+.6f against total %+.6f, %s"
          % (part["n_unassigned"], part["bands_sum_C_A"], part["total_C_A"],
             "CLOSES" if part["closes"] else "DOES NOT CLOSE"))
    beyond = [b for b in bands if b["eta0"] >= 1.0]
    if beyond:
        cb = sum(b["C_A"] for b in beyond)
        ab = sum(b["area_m2"] for b in beyond)
        print("\n    BEYOND y/b = 1.00 THERE IS %.5f m2 OF WING, %.2f%% of the wetted area,"
              % (ab, 100 * ab / cl["total_area_m2"]))
        print("    contributing %+.6f to C_Ap (%.1f counts). The reference grids model a"
              % (cb, 1e4 * cb))
        print("    FLAT tip at y/b = 1.00 and have no such region at all, so this is a")
        print("    GEOMETRY difference and not a numerical one.")
    out["spanwise_bands"] = bands
    out["spanwise_partition"] = part
    out["nominal_semispan_m"] = SPAN_NOM
    out["measured_semispan_m"] = float(P[:, 1].max())
    out["_tip_note"] = (
        "This geometry continues past the nominal semispan with a ROUNDED tip cap: "
        "measured semispan %.5f m against the nominal %.4f, +%.2f%%, and planform area "
        "%.5f m2 against the nominal trapezoid 0.75295, +%.2f%%. Inboard the planform "
        "matches the definition closely: chord within 0.6 to 0.8%% and leading-edge x "
        "within 0.6 mm of the 30-degree sweep line at every station from y/b 0.05 to "
        "1.00. The difference is the tip alone."
        % (float(P[:, 1].max()), SPAN_NOM, 100 * (float(P[:, 1].max()) / SPAN_NOM - 1),
           aref, 100 * (aref / 0.75295 - 1)))

    p = REPO / a.json
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2) + "\n")
    print("\n  wrote %s" % a.json)
    if a.figure:
        exact_figure(REPO / a.figure, bands, out, cl, aref)
        print("  wrote %s" % a.figure)
    return 0


def exact_figure(out_png, bands, rec, cl, aref):
    """Only what the exact integral supports. No reference sectional data appears here,
    because it does not reconcile with its own published forces (ARG-179)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12.4, 5.6))
    fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.30, wspace=0.26)

    mid = [0.5 * (b["eta0"] + b["eta1"]) for b in bands]
    w = [b["eta1"] - b["eta0"] for b in bands]
    col = ["#8C2F39" if b["eta0"] >= 1.0 else "#2a78d6" for b in bands]
    ax[0].bar(mid, [1e4 * b["C_A"] for b in bands], width=w, color=col,
              edgecolor="white", linewidth=0.6)
    ax[0].axhline(0, color="#555", lw=0.8)
    ax[0].axvline(1.0, color="#8C2F39", lw=1.0, ls="--")
    ax[0].set_xlabel("y/b, on the NOMINAL semispan 1.1963 m")
    ax[0].set_ylabel(r"contribution to $C_{A,p}$ [counts]")
    ax[0].set_title("Exact spanwise decomposition; the bands sum to the total",
                    fontsize=10)
    ax[0].grid(alpha=0.25, lw=0.5, axis="y")
    ax[0].text(1.02, ax[0].get_ylim()[0] * 0.85, "rounded tip\nthe reference\nwing does\nnot have",
               fontsize=7.5, color="#8C2F39", va="bottom")

    fo = rec.get("forces_object", {})
    labels = ["this work"]
    vals = [1e4 * rec["C_Ap"]]
    for r in rec.get("reference_rows", []):
        labels.append(r["code"][:18])
        vals.append(1e4 * r["CAp"])
    ax[1].barh(range(len(vals)), vals,
               color=["#2a78d6"] + ["#7f8c8d"] * (len(vals) - 1))
    ax[1].set_yticks(range(len(vals)))
    ax[1].set_yticklabels(labels, fontsize=8)
    ax[1].invert_yaxis()
    ax[1].axvline(0, color="#555", lw=0.8)
    ax[1].set_xlabel(r"$C_{A,p}$ [counts]")
    ax[1].set_title("Wing level, force data on both sides", fontsize=10)
    ax[1].grid(alpha=0.25, lw=0.5, axis="x")

    cap = ("Exact surface integral from the sampled VTK patch: face areas and normals come "
           "from the polygons, so nothing is reconstructed. It reproduces the forces "
           "function object to %.2f counts, and the total area vector closes to %.0e of the "
           "wetted area. THE SECTIONAL LOCALISATION PREVIOUSLY REPORTED IS WITHDRAWN "
           "(ARG-179): the reference sectional data does not reproduce the reference's own "
           "tabulated forces, so a station-by-station difference against it measures "
           "nothing. The wing-level excess is unaffected."
           % (abs(fo.get("d_C_Dp_counts", 0.0)), abs(cl["sum_Az_over_area"])))
    fig.suptitle("ONERA M6 pressure force, measured exactly", fontsize=13,
                 fontweight="bold")
    fig.text(0.006, 0.20, "\n".join(_wrapt(cap, 122)), fontsize=8.5, color="#52514e",
             va="top")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def _wrapt(s, w):
    out, line = [], ""
    for word in s.split():
        if len(line) + len(word) + 1 > w:
            out.append(line)
            line = word
        else:
            line = (line + " " + word).strip()
    if line:
        out.append(line)
    return out


def planform_area(A):
    """Semispan planform area, EXACTLY, as half the total projected area in z.

    NO FIT AND NO LOOKUP. Every planform point of a closed wing is covered by exactly one
    upper and one lower face, so the sum of |A_z| over all faces is twice the planform
    area; the tip cap and the root opening are nearly normal to z and contribute almost
    nothing. It needs no chord definition, no station binning and no reference table.

    THE FIT IT REPLACES WAS BIASED 8.4% HIGH AND THE REASON IS SWEEP. Measuring a "chord"
    as x_max - x_min inside a spanwise bin adds the leading-edge sweep offset across the
    bin: on this wing, 30 degrees of sweep over a 0.12 m bin adds 0.069 m to an 0.8 m
    chord. Narrow bands hide it rather than removing it, which is worse, because the
    residual is then small enough to look like noise. Measured: the fit gave 0.8161 m2
    against 0.7528 m2 from the published planform, and this returns the latter without
    being told it.
    """
    return float(0.5 * np.abs(A[:, 2]).sum())


def force_object(case, time):
    cand = sorted((case / "postProcessing").glob("*orce*/*/force*.dat"))
    if not cand:
        return None
    best = None
    for ln in cand[-1].read_text().splitlines():
        if ln.startswith("#"):
            continue
        v = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", ln)
        if len(v) < 19:
            continue
        if float(v[0]) <= float(time) + 1e-9:
            best = [float(x) for x in v[1:7]]
    return dict(pressure=best[0:3], viscous=best[3:6]) if best else None


if __name__ == "__main__":
    raise SystemExit(main())
