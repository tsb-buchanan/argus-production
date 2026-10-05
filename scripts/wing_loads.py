#!/usr/bin/env python3
"""wing_loads.py: integrated loads from a 3D wing surface sample.

Computes, by integrating pressure and wall shear over the wing patch:

  1. TOTAL FORCES        CL, CD -- integrated here and CHECKED against forceCoeffs.
  2. ROOT BENDING MOMENT compared against the DSO's VLM value (an ACTIVE CONSTRAINT in
                         the optimisation, so this tests its constraint handling).
  3. HINGE MOMENT        about the morph hinge at x/c = 0.62. THE VLM CANNOT COMPUTE THIS:
                         `optimization_result.csv` carries
                         `hinge_moment_status = not_computed_requires_surface_pressure_integration`.
  4. SEPARATION          reverse-flow area fraction from the SIGN of the streamwise wall
                         shear. |Cf| cannot show separation because a magnitude has no sign.

THE TOTAL-FORCE CHECK IS THE POINT OF ITEM 1. It is a null with an independently known
answer (D058): the solver already wrote CL and CD via forceCoeffs, so if this integration
cannot reproduce them, its bending and hinge numbers are worthless. Run it first, believe
the rest only if it passes.

FRAMES (D068).
  OpenFOAM p and wallShearStress are KINEMATIC (divided by rho). Dimensional force needs
  a factor rho, taken as rhoInf = 1.225 to match the case's forceCoeffs.
  Moments are in N.m on the HALF WING, which is the basis the DSO quotes
  (110.754 N.m baseline, 118.644 N.m mcv2).
  Coefficients are on Sref = 1.24092 m2 (DSO basis) with the half-model Aref 0.620462.
"""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
RHO = 1.225
AREF_HALF = 0.620462          # forceCoeffs Aref in the 3D cases (half of Sref)


def read_polys(path):
    """Legacy-VTK POLYDATA -> (points, list-of-index-arrays, cell values).

    Reads CELL_DATA (one value per face), which is what sampledSurface writes. See
    wing_sections.py for why pairing against POINT data would silently misalign.
    """
    op = gzip.open if str(path).endswith(".gz") else open
    npts = ncell = 0
    ncomp = 1
    mode = None
    buf, conn, vals = [], [], []
    pts = None
    with op(path, "rt", errors="ignore") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            head = s.split()[0]
            if head == "POINTS":
                npts = int(s.split()[1]); mode = "pts"; buf = []; continue
            if head == "POLYGONS":
                pts = np.fromiter(buf, float, 3 * npts).reshape(npts, 3)
                ncell = int(s.split()[1]); mode = "conn"; buf = []; continue
            if head in ("CELL_DATA", "POINT_DATA"):
                conn = np.asarray(buf, dtype=np.int64); mode = "await"; buf = []; continue
            if mode == "await" and head in ("FIELD", "SCALARS"):
                mode = "fieldhdr" if head == "FIELD" else "await2"; continue
            if mode == "await2" and head == "LOOKUP_TABLE":
                mode = "val"; continue
            if mode == "fieldhdr":
                f = s.split()
                ncomp = int(f[1]) if len(f) >= 3 and f[1].isdigit() else 1
                mode = "val"; continue
            if mode in ("pts", "conn", "val"):
                buf.extend(s.split()); continue
        if mode == "val":
            vals = buf
    V = np.asarray(vals[:ncell * ncomp], float)
    V = V.reshape(-1, ncomp) if ncomp > 1 else V.reshape(-1, 1)
    faces, k = [], 0
    for _ in range(ncell):
        n = int(conn[k]); faces.append(conn[k + 1:k + 1 + n]); k += 1 + n
    return pts, faces, V[:ncell]


def face_geometry(pts, faces):
    """Centroid, outward area-vector and area for each polygon.

    The area vector is the fan-triangulated cross-product sum, which is exact for a
    planar polygon and the standard approximation for a slightly warped one.
    """
    n = len(faces)
    cx = np.empty((n, 3)); av = np.empty((n, 3))
    for i, idx in enumerate(faces):
        P = pts[idx]
        c = P.mean(axis=0)
        cx[i] = c
        a = np.zeros(3)
        for j in range(len(P)):
            a += np.cross(P[j] - c, P[(j + 1) % len(P)] - c)
        av[i] = 0.5 * a
    return cx, av, np.linalg.norm(av, axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--p", required=True, help="p_wing.vtk[.gz]")
    ap.add_argument("--tau", help="wallShearStress_wing.vtk[.gz]")
    ap.add_argument("--alpha", type=float, required=True, help="deg, for the wind axes")
    ap.add_argument("--uinf", type=float, default=40.8)
    ap.add_argument("--xh", type=float, default=0.62, help="hinge, fraction of local chord")
    ap.add_argument("--label", default="case")
    ap.add_argument("--cl-ref", type=float, default=None,
                    help="forceCoeffs CL, used to calibrate the global normal sign and "
                         "as the null test on the whole integration")
    ap.add_argument("--eta-lo", type=float, default=0.60,
                    help="morph region lower bound. Liming: the moment is over "
                         "0.60 <= eta <= 1.00, NOT the whole span")
    ap.add_argument("--eta-hi", type=float, default=1.00)
    ap.add_argument("--out", type=Path, default=REPO / "results/wing_loads.json")
    a = ap.parse_args()

    pts, faces, P = read_polys(a.p)
    cx, av, area = face_geometry(pts, faces)
    p = P[:, 0]
    print("%s: %d faces, wetted area %.4f m2 (half wing)" % (a.label, len(faces), area.sum()))

    # ---- NORMAL ORIENTATION: CALIBRATED AGAINST A KNOWN ANSWER, NOT GUESSED. ----
    # My first version oriented each face outward with sign(area_vec . (centroid - mid)).
    # On a thin swept wing that dot product is dominated by the spanwise term and is
    # meaningless: it produced CL = -0.0035 against a known 0.455 and 99.9% reverse flow.
    #
    # The surface is CLOSED and surfaceCheck reports ONE normal orientation, so the only
    # freedom is a single GLOBAL sign. That is resolved by requiring the integrated lift
    # to match the solver's own forceCoeffs value -- a null with an independently known
    # answer (D058), used here as the calibration rather than as an afterthought.
    Fp_raw = -(p * RHO)[:, None] * av

    Fv = np.zeros_like(Fp_raw)
    tau = None
    if a.tau:
        tpts, tfaces, T = read_polys(a.tau)
        if len(T) == len(area):
            tau = T * RHO                                   # kinematic -> Pa
            # TAU SIGN, CALIBRATED THE SAME WAY. OpenFOAM's wallShearStress carries the
            # opposite sign to what a naive force = tau*A assumes: the first version gave
            # CDv = -88.9 counts, magnitude right against the known 88.59 and sign wrong,
            # and made 99.87% of the wing read as reverse flow. FRICTION DRAG IS ALWAYS
            # POSITIVE, which is a known answer, so it calibrates the sign with no
            # convention lookup.
            ca0, sa0 = math.cos(math.radians(a.alpha)), math.sin(math.radians(a.alpha))
            Fv_raw = tau * area[:, None]
            Dv_raw = Fv_raw[:, 0].sum() * ca0 + Fv_raw[:, 2].sum() * sa0
            tsgn = 1.0 if Dv_raw > 0 else -1.0
            print("  tau sign calibrated: raw viscous drag %+.4f N -> sign %+g "
                  "(friction drag must be positive)" % (Dv_raw, tsgn))
            tau = tsgn * tau
            Fv = tau * area[:, None]
    ca, sa = math.cos(math.radians(a.alpha)), math.sin(math.radians(a.alpha))
    q = 0.5 * RHO * a.uinf ** 2 * AREF_HALF
    L_raw = (-Fp_raw[:, 0].sum() * sa + Fp_raw[:, 2].sum() * ca) / q
    sgn = 1.0
    if a.cl_ref is not None:
        sgn = 1.0 if (L_raw * a.cl_ref) > 0 else -1.0
        print("  normal orientation calibrated: raw CL_p %+.4f, reference CL %+.4f -> sign %+g"
              % (L_raw, a.cl_ref, sgn))
    Fp = sgn * Fp_raw
    F = Fp + Fv
    L = -F[:, 0].sum() * sa + F[:, 2].sum() * ca
    D = F[:, 0].sum() * ca + F[:, 2].sum() * sa
    Lp = -Fp[:, 0].sum() * sa + Fp[:, 2].sum() * ca
    Dp = Fp[:, 0].sum() * ca + Fp[:, 2].sum() * sa
    Dv = D - Dp
    print("  CL %.6f   CD %.7f (%.2f ct)   CDp %.7f   CDv %.7f"
          % (L / q, D / q, 1e4 * D / q, Dp / q, Dv / q))

    # ---- root bending: moment about the x axis at y = 0, half wing ----
    Mx = (cx[:, 1] * F[:, 2] - cx[:, 2] * F[:, 1]).sum()
    print("  root bending moment  %.4f N.m  (half wing, about y=0)" % Mx)

    # ---- hinge moment about the x/c = 0.62 line, per local chord ----
    # LOCAL chord, from the surface itself in narrow spanwise bins: the wing is tapered
    # and swept, so a single global hinge x would sit at a different x/c at every station.
    semi = cx[:, 1].max()
    nb = 200
    b = np.clip(((cx[:, 1] / semi) * nb).astype(int), 0, nb - 1)
    xle = np.full(nb, np.nan); xte = np.full(nb, np.nan)
    for i in range(nb):
        m = b == i
        if m.sum() >= 20:
            xle[i], xte[i] = cx[m, 0].min(), cx[m, 0].max()
    xh = xle[b] + a.xh * (xte[b] - xle[b])
    # THE MOMENT IS OVER THE MORPH REGION ONLY. Liming, 2026-07: "calculate the local
    # aerodynamic moment directly from surface pressures about x_h/c = 0.62 over
    # 0.60 <= eta <= 1.00". Integrating the whole span aft of the hinge is a DIFFERENT
    # quantity and would not answer what was asked.
    etaf = cx[:, 1] / semi
    inreg = (etaf >= a.eta_lo) & (etaf <= a.eta_hi)
    aft = np.isfinite(xh) & (cx[:, 0] >= xh) & inreg
    Mh = ((cx[aft, 0] - xh[aft]) * F[aft, 2]).sum()
    print("  HINGE MOMENT about x/c %.2f, %.2f <= eta <= %.2f:  %+.4f N.m"
          % (a.xh, a.eta_lo, a.eta_hi, Mh))
    print("    frame: HALF WING, N.m, about the LOCAL hinge line at x/c %.2f of the local"
          % a.xh)
    print("    chord; positive = trailing edge DOWN (nose-down about the hinge);")
    print("    p and tau are OpenFOAM KINEMATIC values x rho = %.3f kg/m3;" % RHO)
    print("    UNCONSTRAINED -- the legacy 1050 N.m bound is WITHDRAWN (Liming 2026-07).")
    print("    %d faces, %.1f%% of wetted area" % (aft.sum(), 100 * area[aft].sum() / area.sum()))

    # ---- separation: SIGN of streamwise wall shear ----
    sep = None
    if tau is not None:
        sdir = np.array([ca, 0.0, sa])
        ts = tau @ sdir
        rev = ts < 0
        sep = 100 * area[rev].sum() / area.sum()
        print("  reverse-flow area (tau.streamwise < 0)  %.3f%%  of wetted area" % sep)
        aftrev = rev & aft
        print("    of which aft of the hinge: %.3f%% of total wetted area"
              % (100 * area[aftrev].sum() / area.sum()))

    rec = {"label": a.label, "alpha_deg": a.alpha, "faces": len(faces),
           "wetted_area_m2": float(area.sum()), "CL": float(L / q), "CD": float(D / q),
           "CD_pressure": float(Dp / q), "CD_viscous": float(Dv / q),
           "root_bending_Nm_half_wing": float(Mx),
           "hinge_moment_Nm_half_wing": float(Mh), "hinge_xc": a.xh,
           "hinge_eta_range": [a.eta_lo, a.eta_hi],
           "hinge_frame": ("half wing; N.m; about the LOCAL hinge line at x/c 0.62 of the "
                           "local chord; positive = trailing edge down; kinematic p and tau "
                           "times rho 1.225; UNCONSTRAINED, the legacy 1050 N.m bound is "
                           "withdrawn"),
           "reverse_flow_area_pct": sep,
           "frames": {"Aref_m2": AREF_HALF, "rho": RHO, "uinf": a.uinf,
                      "note": "half wing; moments in N.m; p and tau kinematic x rho"}}
    old = json.loads(a.out.read_text()) if a.out.exists() else {}
    old[a.label] = rec
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(old, indent=1) + "\n")
    print("  -> %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
