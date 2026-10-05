#!/usr/bin/env python3
"""spanwise_dense_report.py: dense spanwise sectional load from an archived wingSurface.vtk.

Driver for the 2026-09-15 all-geometry report's spanwise-loading figure. It exists because
the RANS side of that figure otherwise carries FIVE stations (the Cp section cuts) against
the VLM's SIXTY-NINE strips, and a sparse curve beside a dense one is read as the noisy one.

NOTHING NUMERICAL IS IMPLEMENTED HERE (D056 forbids a fourth implementation of a thing the
repo already has). Two existing, committed code paths are imported unmodified:

  1. `read_surface_vtk`, `derive_q`, `read_dragdir` from `section_slices_from_surface.py`,
     the extractor whose 24 known-answer tests were run before the 2026-09-15 section cuts;
  2. `integrate_station` from `spanwise_dense.py`, which is the c_n / c_a / c_l quadrature
     already closure-tested against forceCoeffs.

WHAT IS DIFFERENT FROM spanwise_dense.py, and it is a correctness difference, not a
preference. spanwise_dense.py builds Cp itself as `p / (0.5 U^2)` via `wing_sections.section`.
That is the INCOMPRESSIBLE KINEMATIC form and it is WRONG for the M 0.78 cases, whose solver
carries ABSOLUTE pressure in Pa: it would divide 26495.9 Pa by 0.5*233.93^2 and return a Cp
near unity everywhere. This driver takes the solver's OWN `cp` field from the VTK instead,
which is correct in both conditions, and cross-checks it by deriving (q, p_inf) from the
file's own p-vs-cp relation.

FRAME (D068), because every number below is normalised:
  cp      solver's own field, the freestream dynamic pressure of THAT case.
  x/c,z/c LOCAL chord of the band, taken from the surface. Not a reference chord.
  c_l     PRESSURE ONLY. Skin friction excluded (~0.1% of lift; it is NOT negligible for
          drag, which is why no sectional drag is claimed here).
  LOAD    c_l * c / Cref with Cref = 0.39396 m, the DSO basis, matching the RANS force
          coefficients (which sit on Aref 0.620462 m2 half-model, 2*Aref = DSO Sref).
  CLOSURE CL = 2 * int(c_l c dy) / Sref, Sref 1.24092 m2 full-wing DSO. The 2 is the
          half-wing factor, frame-rule instance 1.

The wing is SWEPT, so a constant-y band is not a streamwise section: x/c here is the
band's own streamwise extent, the same convention the .lod strips and the committed
section cuts use. Stated, not assumed away.
"""
import argparse
import csv
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from section_slices_from_surface import read_surface_vtk, derive_q, read_dragdir  # noqa: E402
from spanwise_dense import integrate_station                                      # noqa: E402

SREF_DSO = 1.24092      # m2, FULL wing, DSO basis (project rules, hard fact 2)
CREF_DSO = 0.39396      # m, DSO basis
SEMI_SPAN = 3.6576 / 2.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vtk", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--case-dir", default=None)
    ap.add_argument("--time", default="")
    ap.add_argument("--cl-ref", type=float, default=None,
                    help="the case's own forceCoeffs window-mean C_L, for the closure")
    ap.add_argument("--eta-min", type=float, default=0.02)
    ap.add_argument("--eta-max", type=float, default=0.99)
    ap.add_argument("--n", type=int, default=63)
    ap.add_argument("--band", type=float, default=0.004,
                    help="half-width as a fraction of semispan; matches the committed "
                         "section cuts so the two routes can be differenced")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if not os.path.isdir(a.out):
        os.makedirs(a.out)

    h = hashlib.sha256()
    with open(a.vtk, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    prov = dict(case=a.case, vtk=a.vtk, time=a.time,
                vtk_bytes=os.path.getsize(a.vtk), vtk_sha256=h.hexdigest(),
                script="scripts/spanwise_dense_report.py",
                band_half_width_frac_semispan=a.band,
                Sref_m2_DSO=SREF_DSO, Cref_m_DSO=CREF_DSO,
                frame="c_l on LOCAL chord; LOAD = c_l*c/Cref, Cref DSO 0.39396 m; "
                      "CL closure = 2*int(c_l c dy)/Sref, Sref DSO 1.24092 m2 full-wing")

    print("reading %s" % a.vtk, flush=True)
    P, _area, F, _nrm = read_surface_vtk(a.vtk)
    print("  %d faces, fields: %s" % (len(P), ", ".join(sorted(F))), flush=True)
    for need in ("p", "cp"):
        if need not in F:
            sys.exit("FATAL: %s carries no '%s' field" % (a.vtk, need))

    # WITHIN-FILE CHECK, not a declaration: the file's own p-vs-cp relation must be linear
    # and must return the q of this case. derive_q REFUSES a non-linear relation and a
    # constant cp (its own negative controls), so a pass here is evidence the cp field is
    # the one this case's forces were formed on.
    q, p_inf, rel, p99 = derive_q(F["p"], F["cp"])
    prov.update(q_from_file=q, p_inf_from_file=p_inf,
                q_fit_resid_max_over_q=rel, q_fit_resid_p99_over_q=p99)
    print("  DERIVED FROM THIS FILE: q = %.6g, p_inf = %.6g (resid max %.2e, p99 %.2e)"
          % (q, p_inf, rel, p99), flush=True)

    cp = F["cp"]
    case_dir = a.case_dir or os.path.abspath(os.path.join(os.path.dirname(a.vtk), "..", "..", ".."))
    dhat = read_dragdir(case_dir)
    alpha = float(np.degrees(np.arctan2(dhat[2], dhat[0])))
    prov.update(dragDir=[float(v) for v in dhat], alpha_deg=alpha, case_dir=case_dir)
    print("  dragDir %s -> alpha %.6f deg (from forceCoeffs.dat)"
          % (np.array2string(dhat, precision=9), alpha), flush=True)

    semi = float(P[:, 1].max())
    prov["semispan_from_surface_m"] = semi
    prov["semispan_registered_m"] = SEMI_SPAN
    print("  semispan from surface %.6f m (registered %.6f)" % (semi, SEMI_SPAN), flush=True)

    etas = np.linspace(a.eta_min, a.eta_max, a.n)
    rows, skipped = [], []
    for e in etas:
        e = float(e)
        m = np.abs(P[:, 1] - e * semi) <= a.band * semi
        if m.sum() < 50:
            skipped.append((e, "fewer than 50 faces in band (%d)" % m.sum()))
            continue
        x, z = P[m, 0], P[m, 2]
        c = float(x.max() - x.min())
        if c <= 0:
            skipped.append((e, "degenerate chord"))
            continue
        xc = (x - x.min()) / c
        zc = (z - z.mean()) / c
        r = integrate_station(xc, cp[m], zc, alpha)
        if r is None:
            skipped.append((e, "too few chordwise bins with both surfaces"))
            continue
        cn, ca, cl, nu, nl = r
        rows.append(dict(eta=e, y_m=e * semi, chord_m=c, cn=cn, ca=ca, cl=cl,
                         load_cl_c_over_cref=cl * c / CREF_DSO,
                         n_faces=int(m.sum()), n_bins_upper=nu, n_bins_lower=nl))

    # GEO-089: integrated + skipped-with-reason must equal requested. Anything in neither
    # is a failure, not a skip.
    assert len(rows) + len(skipped) == len(etas), "stations do not partition"
    print("  stations: %d asked = %d integrated + %d skipped"
          % (len(etas), len(rows), len(skipped)), flush=True)
    for e, why in skipped:
        print("    SKIPPED eta %.4f: %s" % (e, why), flush=True)
    if not rows:
        sys.exit("FATAL: no station integrated")

    et = np.array([r["eta"] for r in rows])
    clc = np.array([r["cl"] * r["chord_m"] for r in rows])
    CL_strip = float(2.0 * np.trapz(clc, et * semi) / SREF_DSO)
    prov["CL_from_strips"] = CL_strip
    prov["CL_strip_note"] = (
        "pressure only; the integral spans eta %.3f to %.3f, so the root and tip regions "
        "outside that range are NOT included. The closure gap is therefore an UPPER BOUND "
        "on the sampling error and is not by itself a defect." % (a.eta_min, a.eta_max))
    print("  CL from strips (2x half-wing / Sref %.5f) = %.6f" % (SREF_DSO, CL_strip), flush=True)
    if a.cl_ref is not None:
        prov["CL_forceCoeffs"] = a.cl_ref
        prov["CL_closure_rel"] = (CL_strip - a.cl_ref) / a.cl_ref
        print("  CL from forceCoeffs = %.6f, closure gap %+.3f%%"
              % (a.cl_ref, 100.0 * prov["CL_closure_rel"]), flush=True)

    prov["n_stations"] = len(rows)
    prov["skipped"] = [{"eta": e, "why": w} for e, w in skipped]

    cf = os.path.join(a.out, "%s_spanwise_dense.csv" % a.case)
    with open(cf, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    jf = os.path.join(a.out, "%s_spanwise_dense.json" % a.case)
    with open(jf, "w") as fh:
        fh.write(json.dumps(prov, indent=2) + "\n")
    print("  wrote %s and %s" % (os.path.basename(cf), os.path.basename(jf)), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
