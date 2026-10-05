#!/usr/bin/env python3
"""spanwise_dense.py: dense spanwise load distribution from a 3D wing surface sample.

WHY. The RANS side of the spanwise-loading figure carried SIX stations against the VLM's
SIXTY-NINE strips, so the two panels were being read as though the sparse one were the
noisy one. The sampling was a CHOICE, not a limit: the surface holds 5.28 M faces and the
VTK is read once, so 45 stations cost the same single read as 6.

WHAT THIS ADDS OVER wing_sections.py, which it reuses for parsing:

  1. DENSE STATIONS, matched to the VLM's spanwise resolution.
  2. TRUE SECTIONAL LIFT, not the normal force. The old figure integrated
     cn = int (Cp_l - Cp_u) d(x/c) and LABELLED IT c_l. Those differ by the axial term:
         c_l = c_n cos(alpha) - c_a sin(alpha)
     At alpha 2.1028 deg the correction is small, but "small" is a measurement and not an
     assumption, so c_n, c_a and c_l are all reported and the figure can state which it
     plots (D068: the frame travels with the number).
  3. LOCAL CHORD per station, taken from the surface itself. This wing is tapered, so the
     LOAD c_l * c / cref and the sectional c_l are different curves with different shapes,
     and it is the load that drives CL and root bending.
  4. A CLOSURE, not a coverage count (GEO-089). The strips are integrated back to CL and
     compared against the case's own forceCoeffs. A missing or misplaced strip BREAKS the
     closure rather than quietly shrinking the sample.
  5. PROVENANCE. The existing per-station CSVs record neither the case nor the time they
     came from, so nothing downstream can tell which solution it is plotting. Every output
     here carries case, time, VTK size and checksum, alpha, and the extraction settings.

FRAME (D068), stated because every number below is normalised:
  Cp   = p / (0.5 Uinf^2), OpenFOAM kinematic p (already divided by rho).
  x/c  = LOCAL chord at the station, from the surface, not a global chord.
  c_l  = PRESSURE ONLY. Skin friction is excluded; it is ~0.1% of lift and is NOT
         negligible for drag, which is why no sectional drag is claimed here.
  CL closure carries the HALF-WING FACTOR OF TWO (frame-rule instance 1): the mesh is a
  half-wing on a symmetry plane while Sref is the FULL-wing reference area, so
  CL = 2 * int(c_l c dy) / Sref. Dropping the 2 reads exactly 50% of the truth.
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wing_sections import read_vtk_cell_centres, section  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
SREF_DSO = 1.24092      # m2, full wing, DSO basis (project rules, hard fact 2)
CREF_DSO = 0.39396      # m


def integrate_station(xc, cp, zc, alpha_deg, nbin=200):
    """(c_n, c_a, c_l) from a scattered point cloud of one section.

    THE UPPER/LOWER SPLIT IS PER CHORDWISE BIN, by the median z in that bin, because a
    global z split fails on a cambered section near the LE where both surfaces sit on the
    same side of the mean line. Bins with fewer than 2 points on a side are dropped from
    THAT side only, and the count that survived is returned so the caller can refuse a
    station that was too thinly sampled to integrate.
    """
    e = np.linspace(0.0, 1.0, nbin + 1)
    idx = np.clip(np.digitize(xc, e) - 1, 0, nbin - 1)
    xu, cu, zu, xl, cl_, zl = [], [], [], [], [], []
    for b in range(nbin):
        m = idx == b
        if m.sum() < 4:
            continue
        zb, cb = zc[m], cp[m]
        u = zb >= np.median(zb)
        xm = 0.5 * (e[b] + e[b + 1])
        if u.sum() >= 2:
            xu.append(xm); cu.append(cb[u].mean()); zu.append(zb[u].mean())
        if (~u).sum() >= 2:
            xl.append(xm); cl_.append(cb[~u].mean()); zl.append(zb[~u].mean())
    if len(xu) < 20 or len(xl) < 20:
        return None

    g = np.linspace(0.005, 0.995, nbin)
    cpu = np.interp(g, xu, cu)
    cpl = np.interp(g, xl, cl_)
    zzu = np.interp(g, xu, zu)
    zzl = np.interp(g, xl, zl)

    # c_n = int (Cp_l - Cp_u) d(x/c)
    cn = float(np.trapz(cpl - cpu, g))
    # c_a = int Cp_u d(z_u/c) - int Cp_l d(z_l/c)
    ca = float(np.trapz(cpu, zzu) - np.trapz(cpl, zzl))
    a = np.radians(alpha_deg)
    cl = cn * np.cos(a) - ca * np.sin(a)
    return cn, ca, cl, len(xu), len(xl)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vtk", required=True, help="p_wing.vtk for ONE case")
    ap.add_argument("--tag", required=True, help="baseline | morphed")
    ap.add_argument("--case", default="", help="case name, recorded as provenance")
    ap.add_argument("--time", default="", help="solution time, recorded as provenance")
    ap.add_argument("--alpha", type=float, required=True, help="deg, for the c_n -> c_l rotation")
    ap.add_argument("--uinf", type=float, default=40.8)
    ap.add_argument("--rho", type=float, default=1.225)
    ap.add_argument("--n", type=int, default=45, help="number of stations")
    ap.add_argument("--eta-min", type=float, default=0.06)
    ap.add_argument("--eta-max", type=float, default=0.98)
    ap.add_argument("--band", type=float, default=0.004,
                    help="half-width as a fraction of semispan; MUST match the legacy "
                         "extraction for the null below to mean anything")
    ap.add_argument("--legacy-eta", type=float, nargs="*",
                    default=[0.30, 0.38, 0.50, 0.70, 0.80, 0.90],
                    help="stations the committed CSVs already hold, recomputed as a null")
    ap.add_argument("--cl-ref", type=float, default=None,
                    help="the case's own forceCoeffs CL, for the closure")
    ap.add_argument("--out", type=Path, default=REPO / "results")
    a = ap.parse_args()

    h = hashlib.sha256()
    with open(a.vtk, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    prov = dict(case=a.case, time=a.time, vtk=os.path.basename(a.vtk),
                vtk_bytes=os.path.getsize(a.vtk), vtk_sha256=h.hexdigest()[:16],
                alpha_deg=a.alpha, uinf=a.uinf, rho=a.rho, band=a.band,
                script="scripts/spanwise_dense.py")

    print("reading %s" % a.vtk)
    P, V = read_vtk_cell_centres(a.vtk)
    semi = float(P[:, 1].max())
    prov["semispan_m"] = semi
    print("  {:,} faces, semispan {:.4f} m".format(len(P), semi))

    etas = np.linspace(a.eta_min, a.eta_max, a.n)
    rows, skipped = [], []
    for e in etas:
        r = section(P, V, float(e), semi, a.band, a.uinf)
        if r is None:
            skipped.append((float(e), "fewer than 50 faces in band"))
            continue
        xc, cp, zc = r
        # local chord: section() normalised by it, so recover it from the raw extent
        m = np.abs(P[:, 1] - e * semi) <= a.band * semi
        chord = float(P[m, 0].max() - P[m, 0].min())
        q = integrate_station(xc, cp, zc, a.alpha)
        if q is None:
            skipped.append((float(e), "too few chordwise bins with both surfaces"))
            continue
        cn, ca, cl, nu, nl = q
        rows.append(dict(eta=float(e), y_m=float(e * semi), chord_m=chord,
                         cn=cn, ca=ca, cl=cl,
                         load_cl_c_over_cref=cl * chord / CREF_DSO,
                         lift_per_span_N_per_m=0.5 * a.rho * a.uinf ** 2 * chord * cl,
                         n_faces=int(m.sum()), n_bins_upper=nu, n_bins_lower=nl))

    # EVERY STATION IS ACCOUNTED FOR (GEO-089): integrated + skipped-with-reason = asked.
    assert len(rows) + len(skipped) == len(etas), "stations do not partition"
    print("  stations: %d asked = %d integrated + %d skipped" % (len(etas), len(rows), len(skipped)))
    for e, why in skipped:
        print("    SKIPPED eta %.3f: %s" % (e, why))

    # --- closure, not a coverage count -------------------------------------------
    # CL = 2 * int c_l c dy / Sref. The 2 is the half-wing factor (frame-rule instance 1).
    et = np.array([r["eta"] for r in rows])
    clc = np.array([r["cl"] * r["chord_m"] for r in rows])
    CL_strip = float(2.0 * np.trapz(clc, et * semi) / SREF_DSO)
    prov["CL_from_strips"] = CL_strip
    prov["CL_strip_note"] = ("pressure only, and the integral spans eta %.2f to %.2f so the "
                             "root and tip regions outside that range are NOT included; the "
                             "closure gap below is therefore an upper bound on the sampling "
                             "error, not a defect on its own" % (a.eta_min, a.eta_max))
    print("  CL from strips (2x half-wing / Sref %.5f) = %.6f" % (SREF_DSO, CL_strip))
    if a.cl_ref is not None:
        prov["CL_forceCoeffs"] = a.cl_ref
        prov["CL_closure_rel"] = (CL_strip - a.cl_ref) / a.cl_ref
        print("  CL from forceCoeffs               = %.6f" % a.cl_ref)
        print("  closure gap = %+.2f%%" % (100 * prov["CL_closure_rel"]))

    # --- null: recompute the stations the committed CSVs already hold -------------
    # D057: where a new code path can reproduce an existing artefact, that is a free
    # independent null. It also ESTABLISHES the provenance those CSVs never recorded: if
    # they match, they came from this case at this time.
    nulls = []
    for e in a.legacy_eta:
        f = REPO / "results/wing_sections" / ("%s_eta%s.csv" % (a.tag, ("%.2f" % e).replace(".", "p")))
        if not f.exists():
            nulls.append(dict(eta=e, status="no committed CSV to compare"))
            continue
        r = section(P, V, e, semi, a.band, a.uinf)
        if r is None:
            nulls.append(dict(eta=e, status="band empty now"))
            continue
        new = np.sort(r[1])
        old = np.sort(np.genfromtxt(f, delimiter=",", names=True)["cp"])
        if len(new) != len(old):
            nulls.append(dict(eta=e, status="DIFFERENT SAMPLE COUNT",
                              n_new=len(new), n_old=len(old)))
            continue
        d = float(np.abs(new - old).max())
        nulls.append(dict(eta=e, status="match" if d < 1e-9 else "DIFFERS",
                          max_abs_cp_diff=d, n=len(new)))
    prov["null_vs_committed_sections"] = nulls
    print("  null against committed sections:")
    for n in nulls:
        print("    eta %.2f: %s%s" % (n["eta"], n["status"],
              ("  max|dCp| %.2e over %d pts" % (n["max_abs_cp_diff"], n["n"]))
              if "max_abs_cp_diff" in n else ""))

    a.out.mkdir(parents=True, exist_ok=True)
    cf = a.out / ("spanwise_dense_%s.csv" % a.tag)
    with open(cf, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    jf = a.out / ("spanwise_dense_%s.json" % a.tag)
    prov["n_stations"] = len(rows)
    prov["skipped"] = [{"eta": e, "why": w} for e, w in skipped]
    jf.write_text(json.dumps(prov, indent=2) + "\n")
    print("  wrote %s and %s" % (cf.name, jf.name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
