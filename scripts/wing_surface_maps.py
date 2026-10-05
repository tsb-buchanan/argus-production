#!/usr/bin/env python3
"""wing_surface_maps.py: planform contour maps and a separation analysis from the wing surface.

WHY. Two things the project has never had, both derivable from surface samples ALREADY ON
DISK with no new solves:

  1. CONTOUR MAPS of Cp and Cf over the planform, baseline against morphed, and their
     DIFFERENCE. Sectional cuts show what happens at six or seventy spanwise stations; a map
     shows where the morph acts in two dimensions, including the inboard/outboard edges of
     the morphed panel which no cut passes through.
  2. A SEPARATION ANALYSIS with a SIGN TEST rather than a threshold. Reverse flow is
     wallShearStress with a NEGATIVE streamwise component. That is a one-directional
     quantity, so the sign is the test and a magnitude tolerance would be strictly weaker
     (D060). Attached flow CANNOT produce tau_s < 0, so any negative area is separation, not
     scatter.

THE STREAMWISE DIRECTION IS THE FREESTREAM, NOT GLOBAL X, and that is a frame statement
(D068). At alpha the freestream is s_hat = (cos a, 0, sin a), so tau_s = tau . s_hat. Using
global x would mis-sign faces near the leading edge by the alpha rotation. The wing is also
SWEPT, so tau_s is the freestream-aligned component and NOT the section-normal chordwise
component; those differ by the local sweep angle and the difference is recorded rather than
assumed away.

WHAT IS BINNED AND WHY. Per-face data is 5.28 M faces per field per case, which is too large
to carry off the cluster. Everything is binned onto a structured (eta, x/c) grid, SEPARATELY
for the upper and lower surface, because a planform bin contains both and averaging across
them would cancel the very Cp difference being plotted.

FRAME (D068):
  Cp  = p / (0.5 Uinf^2), OpenFOAM KINEMATIC p (already divided by rho).
  Cf  = |tau| / (0.5 Uinf^2), KINEMATIC wallShearStress, likewise.
  x/c = LOCAL chord at that spanwise station, taken from the surface; the wing is tapered.
"""
import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wing_sections import read_vtk_cell_centres  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
NU = 1.46e-5      # sea level, ARGUS_reference_data.md section 9


def sha16(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()[:16]


def planform_frame(P, n_eta):
    """eta, local chord, x/c and an upper/lower mask for every face.

    THE UPPER/LOWER SPLIT IS PER SPANWISE STRIP, against that strip's own mean camber, not
    against a global z. A swept, tapered, twisted wing puts the outboard lower surface above
    the inboard upper surface in absolute z, so a global split mislabels whole panels.
    """
    semi = P[:, 1].max()
    eta = P[:, 1] / semi
    edges = np.linspace(0.0, 1.0, n_eta + 1)
    idx = np.clip(np.digitize(eta, edges) - 1, 0, n_eta - 1)
    xc = np.full(len(P), np.nan)
    chord = np.full(len(P), np.nan)
    upper = np.zeros(len(P), dtype=bool)
    for b in range(n_eta):
        m = idx == b
        if m.sum() < 50:
            continue
        x, z = P[m, 0], P[m, 2]
        c = x.max() - x.min()
        if c <= 0:
            continue
        t = (x - x.min()) / c
        chord[m] = c
        xc[m] = t
        # camber line in 40 chordwise bins, then split about it
        cb = np.linspace(0, 1, 41)
        ci = np.clip(np.digitize(t, cb) - 1, 0, 39)
        zc = np.array([z[ci == k].mean() if (ci == k).any() else np.nan for k in range(40)])
        good = ~np.isnan(zc)
        zline = np.interp(t, (0.5 * (cb[:-1] + cb[1:]))[good], zc[good])
        u = np.zeros(m.sum(), dtype=bool)
        u[z >= zline] = True
        upper[m] = u
    return eta, chord, xc, upper, semi


def bin2d_wmean(eta, xc, val, w, mask, n_eta, n_xc):
    """AREA-WEIGHTED mean of val on an (eta, x/c) grid. Used for the separation map.

    BINNING THE MEAN OF A FIELD AND THEN THRESHOLDING IT IS NOT THE SAME AS BINNING THE
    FRACTION THAT CROSSES THE THRESHOLD, and the difference is the whole map. Reverse flow
    is 0.13% of wetted area, so almost no bin has a negative MEAN shear and a map built
    that way is blank while the separation is really there. The fraction is the quantity
    that answers "how much of this bin is separated".
    """
    ok = mask & np.isfinite(xc) & np.isfinite(val)
    ei = np.clip((eta[ok] * n_eta).astype(int), 0, n_eta - 1)
    xi = np.clip((xc[ok] * n_xc).astype(int), 0, n_xc - 1)
    flat = ei * n_xc + xi
    num = np.bincount(flat, weights=(val[ok] * w[ok]), minlength=n_eta * n_xc)
    den = np.bincount(flat, weights=w[ok], minlength=n_eta * n_xc)
    out = np.full(n_eta * n_xc, np.nan)
    nz = den > 0
    out[nz] = num[nz] / den[nz]
    return out.reshape(n_eta, n_xc)


def bin2d(eta, xc, val, mask, n_eta, n_xc):
    """Mean of val on an (eta, x/c) grid over mask. Empty bins -> nan, never zero."""
    ok = mask & np.isfinite(xc) & np.isfinite(val)
    ei = np.clip((eta[ok] * n_eta).astype(int), 0, n_eta - 1)
    xi = np.clip((xc[ok] * n_xc).astype(int), 0, n_xc - 1)
    flat = ei * n_xc + xi
    s = np.bincount(flat, weights=val[ok], minlength=n_eta * n_xc)
    c = np.bincount(flat, minlength=n_eta * n_xc)
    out = np.full(n_eta * n_xc, np.nan)
    nz = c > 0
    out[nz] = s[nz] / c[nz]
    return out.reshape(n_eta, n_xc)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", required=True, help="wingSurface/<time> directory")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--case", default="")
    ap.add_argument("--time", default="")
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--uinf", type=float, default=40.8)
    ap.add_argument("--n-eta", type=int, default=240)
    ap.add_argument("--n-xc", type=int, default=320)
    ap.add_argument("--out", type=Path, default=REPO / "results/wing_maps")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    q = 0.5 * a.uinf ** 2

    fp = Path(a.dir) / "p_wing.vtk"
    fw = Path(a.dir) / "wallShearStress_wing.vtk"
    fy = Path(a.dir) / "yPlus_wing.vtk"
    print("reading %s" % fp)
    P, p, area = read_vtk_cell_centres(fp, with_area=True)
    print("reading %s" % fw)
    Pw, tau = read_vtk_cell_centres(fw, as_vector=True)
    if len(Pw) != len(P):
        # THE TWO FILES MUST DESCRIBE THE SAME FACES. Pairing different face sets would
        # attach each shear vector to the wrong location and still plot smoothly.
        sys.exit("FATAL: p has %d faces, wallShearStress has %d; not the same surface"
                 % (len(P), len(Pw)))
    d = np.abs(P - Pw).max()
    if d > 1e-9:
        sys.exit("FATAL: face centroids differ by %.3e m between the two files" % d)

    yp = None
    if fy.exists():
        print("reading %s" % fy)
        Py, yp = read_vtk_cell_centres(fy)
        if len(Py) != len(P) or np.abs(P - Py).max() > 1e-9:
            sys.exit("FATAL: yPlus is not on the same faces as p")

    eta, chord, xc, upper, semi = planform_frame(P, a.n_eta)
    cp = p / q
    cf = np.linalg.norm(tau, axis=1) / q
    al = np.radians(a.alpha)
    shat = np.array([np.cos(al), 0.0, np.sin(al)])
    tau_s = tau @ shat                      # SIGNED streamwise shear

    # --- CALIBRATE THE SIGN CONVENTION AGAINST A KNOWN ANSWER (D058 rule 2) ----------
    # OpenFOAM's wallShearStress function object does NOT use the intuitive sign: for
    # attached flow in +x it returns a NEGATIVE x-component. Taking the intuitive sign on
    # faith reported 99.87% of the wing as reverse flow, which is a gate that ALWAYS FIRES
    # rather than a finding. So the convention is MEASURED here, not assumed, against a
    # region whose answer is known independently: the upper surface between x/c 0.20 and
    # 0.60 at eta 0.15-0.50 is inboard of the morph, at 2.1 deg on a 12% section at
    # Re ~1e6, and is certainly attached. If that region is NOT overwhelmingly one sign,
    # the calibration is invalid and the run HALTS rather than picking a sign anyway.
    ref = (valid_pre := np.isfinite(xc)) & upper & (xc >= 0.20) & (xc <= 0.60) \
        & (eta >= 0.15) & (eta <= 0.50)
    if ref.sum() < 5000:
        sys.exit("FATAL: only %d faces in the attached reference region; cannot calibrate"
                 % ref.sum())
    frac_neg = float((tau_s[ref] < 0).mean())
    if frac_neg > 0.99:
        sign, conv = -1.0, "OpenFOAM (negative for forward flow); flipped here"
    elif frac_neg < 0.01:
        sign, conv = +1.0, "intuitive (positive for forward flow); used as-is"
    else:
        sys.exit("FATAL: the attached reference region is %.1f%% negative, which is neither "
                 "convention. It is therefore not uniformly attached and cannot calibrate "
                 "a sign test." % (100 * frac_neg))
    tau_s = sign * tau_s
    print("  sign calibration: reference region %d faces, %.3f%% negative -> %s"
          % (ref.sum(), 100 * frac_neg, conv))
    cfs = tau_s / q

    # --- separation, by SIGN (D060) ------------------------------------------------
    # Attached flow cannot produce tau_s < 0, so no tolerance is needed or wanted. Faces
    # with no valid planform frame are counted in their own bucket, never dropped silently.
    valid = valid_pre
    rev = valid & (tau_s < 0.0)
    # AREA-WEIGHTED, NOT FACE-COUNTED. Face areas on this mesh vary by orders of
    # magnitude between the refined trailing edge and the coarse mid-span, and separation
    # lives precisely in the refined region, so a face-count fraction OVER-weights exactly
    # the places where separation occurs. The committed claim is phrased as an AREA
    # fraction, so it must be weighted by area to be the same quantity at all.
    A_valid = float(area[valid].sum())
    A_rev = float(area[rev].sum())
    sep = dict(sign_convention=conv, sign_applied=sign,
               calib_ref_faces=int(ref.sum()), calib_ref_frac_negative=frac_neg,
               n_faces=int(len(P)), n_valid=int(valid.sum()),
               n_unframed=int((~valid).sum()),
               n_reverse=int(rev.sum()),
               reverse_fraction_by_count=float(rev.sum() / max(valid.sum(), 1)),
               wetted_area_m2_valid=A_valid, reverse_area_m2=A_rev,
               reverse_fraction_of_area=A_rev / max(A_valid, 1e-30))
    # AFT OF THE HINGE is the number that actually answers Liming's question: does the
    # AFT CAMBER separate? Reverse flow at the root junction and the tip rim is present on
    # both wings and has nothing to do with the morph, so a whole-wing figure would let
    # those two regions answer a question they are not about.
    aft = valid & (xc >= 0.62)
    sep["aft_of_hinge_area_m2"] = float(area[aft].sum())
    sep["aft_of_hinge_reverse_area_m2"] = float(area[aft & rev].sum())
    sep["aft_of_hinge_reverse_fraction"] = float(area[aft & rev].sum() / max(area[aft].sum(), 1e-30))
    # and the same restricted to the MORPHED span, which is the sharpest form of the test
    mrph = aft & (eta >= 0.6014)
    sep["aft_of_hinge_morphed_span_reverse_fraction"] = float(
        area[mrph & rev].sum() / max(area[mrph].sum(), 1e-30))
    # BOTH NORMALISATIONS, BECAUSE THE DENOMINATOR IS THE WHOLE DISAGREEMENT (D068).
    # scripts/wing_loads.py reports this numerator over TOTAL WETTED AREA, and the report
    # quotes that as "aft of the hinge 0.008%". Over the REGION'S OWN AREA the same
    # numerator is ~0.08%, ten times larger and answering the more natural question "how
    # much of the aft camber is separated". Neither is wrong; a figure without its
    # denominator is.
    sep["aft_of_hinge_morphed_span_reverse_over_TOTAL_wetted"] = float(
        area[mrph & rev].sum() / max(A_valid, 1e-30))
    sep["aft_of_hinge_morphed_span_area_m2"] = float(area[mrph].sum())
    assert sep["n_valid"] + sep["n_unframed"] == sep["n_faces"], "faces do not partition"

    # where it is, per spanwise station: fraction of that strip's faces in reverse flow,
    # and the most UPSTREAM aft-of-midchord reversal, which is the separation onset
    nst = 40
    st_edges = np.linspace(0.0, 1.0, nst + 1)
    si = np.clip(np.digitize(eta, st_edges) - 1, 0, nst - 1)
    stations = []
    for b in range(nst):
        m = valid & (si == b)
        if m.sum() < 200:
            continue
        r = rev & m
        onset = float(np.min(xc[r & (xc > 0.5)])) if (r & (xc > 0.5)).any() else None
        stations.append(dict(eta=float(0.5 * (st_edges[b] + st_edges[b + 1])),
                             n=int(m.sum()), onset_xc=onset,
                             reverse_fraction=float(area[r].sum() / max(area[m].sum(), 1e-30)),
                             reverse_fraction_upper=float(area[r & upper].sum()
                                                          / max(area[m & upper].sum(), 1e-30))))
    sep["stations"] = stations

    fields = [("cp", cp), ("cf", cf), ("cfs", cfs)]
    if yp is not None:
        # THE FIRST-CELL HEIGHT IS THE DISCRIMINATOR, NOT y+ ITSELF. Both y+ and Cf scale
        # with the wall shear, so a streak in the flow shows up in BOTH and a y+ map alone
        # cannot say whether a streak is mesh or physics. The wall spacing implied by the
        # two sampled fields,
        #     y_wall = y+ * nu / u_tau ,   u_tau = sqrt(|tau|)   (tau is KINEMATIC here)
        # is a property of the MESH and nothing else. If it streaks, prism layers dropped
        # out; if it is smooth while Cf streaks, the streak is in the flow.
        utau = np.sqrt(np.maximum(np.linalg.norm(tau, axis=1), 1e-30))
        ywall = yp * NU / utau
        fields += [("yplus", yp), ("ywall_um", 1e6 * ywall)]
        prov_extra = dict(
            nu=NU,
            ywall_definition="y_wall = yPlus * nu / sqrt(|tau_kinematic|), in micrometres",
            ywall_request_um=4.5,
            ywall_median_um=float(np.nanmedian(1e6 * ywall[np.isfinite(xc)])),
            ywall_p90_um=float(np.nanpercentile(1e6 * ywall[np.isfinite(xc)], 90)),
            ywall_p99_um=float(np.nanpercentile(1e6 * ywall[np.isfinite(xc)], 99)),
            # ACHIEVED LAYER THICKNESS CAN ONLY BE AT OR BELOW THE REQUEST (traps item 2),
            # so a first cell ABOVE 4.5 um is layer dropout and is one-directional evidence.
            frac_above_request=float(np.mean((1e6 * ywall[np.isfinite(xc)]) > 4.5)))
    else:
        prov_extra = {}

    grids = {}
    for name, v in fields:
        for side, msk in (("upper", upper & valid), ("lower", (~upper) & valid)):
            grids["%s_%s" % (name, side)] = bin2d(eta, xc, v, msk, a.n_eta, a.n_xc)
    # area-weighted REVERSE-FLOW FRACTION per bin, which is what the separation map needs
    revf = (tau_s < 0).astype(float)
    for side, msk in (("upper", upper & valid), ("lower", (~upper) & valid)):
        grids["sepfrac_%s" % side] = bin2d_wmean(eta, xc, revf, area, msk, a.n_eta, a.n_xc)

    prov = dict(**prov_extra, case=a.case, time=a.time, tag=a.tag, alpha_deg=a.alpha,
                uinf=a.uinf,
                semispan_m=float(semi), n_eta=a.n_eta, n_xc=a.n_xc,
                p_vtk_sha256=sha16(fp), wss_vtk_sha256=sha16(fw),
                p_vtk_bytes=os.path.getsize(fp), wss_vtk_bytes=os.path.getsize(fw),
                script="scripts/wing_surface_maps.py",
                frame=("Cp = p/(0.5 U^2) kinematic; Cf = |tau|/(0.5 U^2) kinematic; "
                       "cfs is the SIGNED component along the freestream (cos a, 0, sin a); "
                       "x/c on the LOCAL chord"),
                separation=sep)

    np.savez_compressed(a.out / ("maps_%s.npz" % a.tag), **grids)
    (a.out / ("maps_%s.json" % a.tag)).write_text(json.dumps(prov, indent=2) + "\n")
    print("  %d faces, semispan %.4f m" % (len(P), semi))
    print("  aft of hinge x/c>=0.62: %.4f%% reverse (of that region's own area)"
          % (100 * sep["aft_of_hinge_reverse_fraction"]))
    print("  morphed span aft of hinge: %.4f%% of its OWN area = %.4f%% of TOTAL wetted"
          % (100 * sep["aft_of_hinge_morphed_span_reverse_fraction"],
             100 * sep["aft_of_hinge_morphed_span_reverse_over_TOTAL_wetted"]))
    print("  reverse flow: %.4f%% BY AREA (%.5f of %.4f m2), %.4f%% by face count"
          % (100 * sep["reverse_fraction_of_area"], A_rev, A_valid,
             100 * sep["reverse_fraction_by_count"]))
    ons = [s for s in stations if s["onset_xc"] is not None]
    if ons:
        print("  earliest aft separation onset: x/c %.3f at eta %.3f"
              % (min(o["onset_xc"] for o in ons),
                 min(ons, key=lambda o: o["onset_xc"])["eta"]))
    if prov_extra:
        print("  first cell implied by y+ and tau: median %.2f um, p90 %.2f, p99 %.2f "
              "(request 4.5); %.3f%% of faces ABOVE the request"
              % (prov_extra["ywall_median_um"], prov_extra["ywall_p90_um"],
                 prov_extra["ywall_p99_um"], 100 * prov_extra["frac_above_request"]))
    print("  wrote maps_%s.npz and maps_%s.json" % (a.tag, a.tag))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
