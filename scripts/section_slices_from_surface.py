#!/usr/bin/env python3
"""section_slices_from_surface.py: Cp, Cf and y+ sections cut from an archived wingSurface.vtk.

REQUIREMENT, 2026-09-11: the wing-surface VTK carries every field needed, so 2D comparison
plots along the aerofoils are slices of that VTK.

WHY THIS EXISTS BESIDE wing_sections.py RATHER THAN INSTEAD OF IT. Two things changed under
that script and both are silent:

  1. THE ARCHIVE HOLDS ONE COMBINED FILE, NOT ONE FILE PER FIELD. wing_sections.py takes
     `--base p_wing.vtk`, the OF-7-era layout. The OF-12 `argusWingSurface` function object
     writes a single wingSurface.vtk carrying FIELD attributes 4: p, cp, yPlus and
     wallShearStress. read_vtk_cell_centres stops after the FIRST field, so pointed at the
     combined file it returns `p` and CANNOT reach wallShearStress. It happens to return the
     right array, by truncation rather than by design, which is why this is worth stating:
     the old reader is not wrong on the new file, it is BLIND to three quarters of it.
  2. Cf WAS NEVER AVAILABLE AT A SECTION AT ALL, only as a planform map.

THE FRAME, AND IT IS THE WHOLE REASON THIS SCRIPT DOES NOT HARDCODE q (D068).
------------------------------------------------------------------------------
wing_sections.py and wing_surface_maps.py both compute

    Cp = p / (0.5 Uinf^2)          # OpenFOAM KINEMATIC p, incompressible

which is correct on the L11 wall-resolved cases and WRONG on the M6 compressible ones.
Measured on the archives, 2026-09-11:

    SWC_trim     p = -124.659          cp = -0.215673     -> p/cp  =  578.0 = 0.5 * 34.0^2
    CMPB_a1p35   p =  22921.6 Pa       cp = -0.316512     -> q = 11292 Pa, p_inf = 26496 Pa

The compressible solver writes TRUE forces and ABSOLUTE pressure, which ARG-188 already paid
for once: blind rho application there was wrong by a factor 2.42. A single `--uinf` cannot
express both paths, and a per-case constants table is a thing to get out of step.

SO q AND p_inf ARE DERIVED FROM THE FILE, by least squares on the case's own two fields:

    p = q * cp + p_inf

This is a WITHIN-FILE PROOF (D052/D053 item 4): the file convicts itself and needs no
external arbiter, no condition lookup and no per-case table. It returns q in whatever units
that case's `p` is in, so

    Cf = |tau| / q

is then automatically right on BOTH paths, because `wallShearStress` is in the same units as
`p` in both. The fit residual is ASSERTED, not printed: if `cp` were ever written against a
different reference than `p`, the linear relation breaks and this HALTS rather than
returning a plausible number.

THE UPPER/LOWER SPLIT IS DELIBERATELY NOT DONE HERE, and that is a portability decision
rather than a design preference. HPC12's system python is 3.6 with numpy 1.12, which has
neither `argsort(kind="stable")` nor `lstsq(rcond=None)`. The shared splitter
airfoil_surfaces.split_branches uses the former, and it is the GEOMETRY AGENT'S shared
utility: editing it to suit this cluster would be changing another workstream's file to work
around my environment, and reimplementing it is what D056 forbids in terms (three
implementations were tried and rejected; the sign-of-z one shipped briefly and inflated a
reported residual threefold on exactly the aft-cambered sections this wing carries).

So this writes the RAW BAND, sorted along the chord, and the split happens at PLOT TIME on
the laptop where the shared splitter runs unmodified. That is also wing_sections.py's stated
pattern: consume the ~600 MB VTK once on the cluster, write small CSVs, and keep the figure
script runnable with no cluster access.
"""
import argparse
import csv
import gzip
import itertools
import sys
from pathlib import Path

import numpy as np

SEMI_SPAN_DEFAULT = 3.6576 / 2.0     # bref from ARGUS_reference_data.md, half of it


def _tokens(fh):
    for line in fh:
        for t in line.split():
            yield t


def read_surface_vtk(path, with_area=False, with_normals=False):
    """Legacy ASCII VTK POLYDATA with N CELL_DATA fields -> (centroids, area, fields, normals).

    `area` is None unless asked for: it costs a second 5.28M-iteration Python loop and no
    section output uses it. Kept available because a wetted-area or force closure would.

    Parsed as a TOKEN STREAM rather than line by line, because a field header
    ("cp 1 5276970 float") and its data are not separated by any structure the line
    iterator can see: after exactly ncomp*ntuples floats, the next four tokens ARE the next
    header. Consuming exact counts makes that boundary explicit instead of guessed.

    CELL_DATA, not POINT_DATA: OpenFOAM's sampledSurface writes one value per FACE and there
    are FEWER polygons than points, so pairing values against the point array would misalign
    every value by a growing offset and still draw a plausible curve.
    """
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", errors="ignore") as fh:
        tk = _tokens(fh)
        pts = conn = None
        ncell = 0
        fields = {}
        while True:
            t = next(tk, None)
            if t is None:
                break
            if t == "POINTS":
                n = int(next(tk)); next(tk)          # count, dtype
                pts = np.fromiter(itertools.islice(tk, 3 * n), dtype=np.float64,
                                  count=3 * n).reshape(n, 3)
            elif t == "POLYGONS":
                ncell = int(next(tk)); size = int(next(tk))
                conn = np.fromiter(itertools.islice(tk, size), dtype=np.int64, count=size)
            elif t == "CELL_DATA":
                ncell = int(next(tk))
            elif t == "FIELD":
                next(tk)                              # "attributes"
                nf = int(next(tk))
                for _ in range(nf):
                    name = next(tk); ncomp = int(next(tk)); ntup = int(next(tk)); next(tk)
                    v = np.fromiter(itertools.islice(tk, ncomp * ntup), dtype=np.float64,
                                    count=ncomp * ntup)
                    fields[name] = v.reshape(ntup, ncomp) if ncomp > 1 else v
    if pts is None or conn is None:
        sys.exit("FATAL: %s carries no POINTS/POLYGONS block" % path)

    # Walk the flat connectivity ONCE to get each face's start and size. The recurrence
    # start[i+1] = start[i] + 1 + size[i] is inherently sequential; everything after it is
    # vectorised.
    starts = np.empty(ncell, dtype=np.int64)
    sizes = np.empty(ncell, dtype=np.int64)
    k = 0
    for c in range(ncell):
        starts[c] = k
        n = conn[k]
        sizes[c] = n
        k += 1 + n
    vmask = np.ones(len(conn), dtype=bool)
    vmask[starts] = False
    vidx = conn[vmask]
    gstart = np.concatenate(([0], np.cumsum(sizes)[:-1]))
    cen = np.add.reduceat(pts[vidx], gstart, axis=0) / sizes[:, None]

    # Fan triangulation about the centroid: valid for the mixed 3/4/5-sided polygons
    # snappyHexMesh produces and for non-convex ones, where a shoelace about an arbitrary
    # vertex is not.
    area = None
    if with_area:
        area = np.zeros(ncell)
        for c in range(ncell):
            q = pts[conn[starts[c] + 1: starts[c] + 1 + sizes[c]]]
            v = q - cen[c]
            area[c] = 0.5 * np.linalg.norm(np.cross(v, np.roll(v, -1, axis=0)).sum(axis=0))

    # PER-FACE NORMAL BY NEWELL'S METHOD: 0.5 * sum_i (v_i x v_{i+1}) round each face's
    # closed loop. Exact for a planar polygon and origin-independent because the loop
    # closes; robust on the non-planar and non-convex n-gons snappyHexMesh writes, where a
    # single cross product of the first three vertices is not.
    #
    # CHUNKED, and the reason is measured rather than cautious. Line 142 already
    # materialises pts[vidx], which is 652 MB on a real 5.28M-face archive surface; doing
    # the cross product on the whole array at once needs three such buffers at peak. Chunks
    # of 500k faces hold it near 50 MB, and `starts` is already sequential so the face-range
    # to vertex-range mapping is exact rather than estimated.
    #
    # NO STABLE ARGSORT AND NO lstsq HERE, deliberately: this must run on HPC12's numpy
    # 1.12, where both are unavailable (see the module docstring and derive_q). Everything
    # below is cross, arange and add.reduceat, all of which are present there.
    nrm = None
    if with_normals:
        nrm = np.empty((ncell, 3))
        chunk = 500000
        for a0 in range(0, ncell, chunk):
            a1 = min(a0 + chunk, ncell)
            v0 = int(gstart[a0])
            v1 = int(gstart[a1 - 1] + sizes[a1 - 1])
            sub = vidx[v0:v1]
            g = gstart[a0:a1] - v0
            nxt = np.arange(1, len(sub) + 1, dtype=np.int64)
            nxt[g + sizes[a0:a1] - 1] = g        # last vertex of each face wraps to its first
            nrm[a0:a1] = 0.5 * np.add.reduceat(
                np.cross(pts[sub], pts[sub[nxt]]), g, axis=0)
    return cen, area, fields, nrm


def read_dragdir(case_dir):
    """Freestream unit vector, READ FROM THE FORCE FILE THE SOLVER WROTE.

    Not from system/controlDict, and the distinction is the D080 one: the dictionary is
    what we ASKED for, the forceCoeffs.dat header is what the run actually USED. They agree
    here, and the one worth trusting is the output.

        # dragDir       : (9.99756909e-01 0.00000000e+00 2.20482143e-02)

    NOT global x. At alpha the freestream is rotated, and using global x would mis-sign
    faces near the leading edge by exactly that rotation (D068). The wing is also SWEPT, so
    this is the FREESTREAM-ALIGNED component and NOT the section-normal chordwise one; the
    two differ by the local sweep angle and that is recorded rather than assumed away.
    """
    # LATEST LEG FIRST, SORTED NUMERICALLY. A plain sorted() is LEXICOGRAPHIC, so a case
    # with legs {0, 4000, 5500} yields "0" first -- the COLD leg, at the pre-nudge alpha.
    # Measured on SWM_trim: that is 1.252018 deg against the trimmed 1.257423, and the
    # extraction then ran on wind axes 0.0054 deg from the ones the run used (D019).
    def _t(p):
        try:
            return float(p.parent.name)
        except ValueError:
            return float("-inf")       # a non-numeric dir sorts last, it does not raise
    cands = sorted(Path(case_dir).glob("postProcessing/forceCoeffs1/*/forceCoeffs.dat"),
                   key=_t, reverse=True)
    if not cands:
        sys.exit("FATAL: no forceCoeffs.dat under %s; cannot read the freestream direction"
                 % case_dir)
    for f in cands:
        with open(f, errors="ignore") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if "dragDir" in line:
                    nums = [float(t) for t in
                            line.replace("(", " ").replace(")", " ").split()
                            if _isnum(t)]
                    if len(nums) >= 3:
                        d = np.array(nums[:3], dtype=float)
                        n = float(np.linalg.norm(d))
                        if abs(n - 1.0) > 1e-6:
                            sys.exit("FATAL: dragDir in %s has norm %.9f, not a unit vector"
                                     % (f, n))
                        return d / n
    sys.exit("FATAL: no 'dragDir' header line in any forceCoeffs.dat under %s" % case_dir)


def _isnum(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def calibrate_shear_sign(tau_s, eta_f, xc_f):
    """+1 or -1 so that ATTACHED flow comes out POSITIVE, measured not assumed.

    OpenFOAM's wallShearStress function object does NOT use the intuitive sign: for
    attached flow in +x it returns a NEGATIVE x-component. wing_surface_maps.py records
    that taking the intuitive sign on faith reported 99.87% of the wing as reverse flow,
    which is a gate that ALWAYS FIRES rather than a finding.

    So the convention is measured against a region whose answer is known independently:
    INBOARD (eta <= 0.5) MID-CHORD (0.20 <= x/c <= 0.60). That is inboard of the morphed
    panel, on a 12% section at Re ~1e6 and a degree or so of incidence, and it is certainly
    attached on both surfaces. Using both surfaces is what lets this run on the cluster,
    where the upper/lower splitter is unavailable.

    HALTS if that region is not overwhelmingly one sign, rather than picking a sign anyway:
    an ambiguous calibration means the reference assumption is wrong, and a coin-flip there
    would invert every separation verdict downstream.
    """
    ref = (eta_f <= 0.5) & (xc_f >= 0.20) & (xc_f <= 0.60) & np.isfinite(tau_s)
    n = int(ref.sum())
    if n < 1000:
        sys.exit("FATAL: only %d faces in the sign-calibration region; too few to calibrate"
                 % n)
    pos = float((tau_s[ref] > 0).mean())
    frac = max(pos, 1.0 - pos)
    if frac < 0.95:
        sys.exit("FATAL: the attached reference region is only %.1f%% one sign. The sign "
                 "convention cannot be calibrated and a guess would invert every "
                 "separation verdict. Refusing." % (100.0 * frac))
    sign = 1.0 if pos > 0.5 else -1.0
    return sign, frac, n


def calibrate_surface_flag(nz_ref, cp_ref, z_ref, min_group=500, min_dcp=0.05):
    """+1 or -1 so that (mult * nz) >= 0 selects the UPPER surface, measured not assumed.

    THE FACET NORMAL SEPARATES THE TWO SURFACES; IT DOES NOT SAY WHICH IS WHICH. These
    surfaces' normals point INWARD on the OF-12 exports (measured: geometrically upper
    mid-chord faces carry mean nz -0.997), so `nz >= 0` picks the LOWER surface there. A
    hard-coded sign would be right on one export and silently inverted on the next, which is
    the same trap calibrate_shear_sign exists to avoid, and inverting this flag swaps every
    separation statistic between the two surfaces.

    SO IT IS MEASURED AGAINST A REGION WHOSE ANSWER IS KNOWN INDEPENDENTLY: mid-chord,
    x/c 0.35 to 0.65, where the upper surface carries the MORE NEGATIVE cp on any lifting
    section. That is the primary criterion.

    AND CONFIRMED BY A SECOND, INDEPENDENT ONE: at mid-chord the upper surface also sits at
    HIGHER z. The two must agree. This is deliberately not used alone, because the module
    docstring records that a sign-of-z splitter "shipped briefly and inflated a reported
    residual threefold on exactly the aft-cambered sections this wing carries" -- z ordering
    is unreliable near an aft-drooped trailing edge. At MID-CHORD it is reliable, which is
    why the cross-check is confined there and why disagreement is a refusal rather than a
    tie-break: two criteria that disagree mean the reference assumption is wrong.

    HALTS rather than guessing, in all three failure modes, because a coin-flip here would
    mislabel every reversed face and move the reversal to the wrong surface.
    """
    a = nz_ref >= 0.0
    na, nb = int(a.sum()), int((~a).sum())
    if min(na, nb) < min_group:
        sys.exit("FATAL: the normal-sign groups hold %d and %d faces in the mid-chord "
                 "reference region; too few to calibrate the upper/lower flag" % (na, nb))
    dcp = float(cp_ref[a].mean() - cp_ref[~a].mean())
    dz = float(np.median(z_ref[a]) - np.median(z_ref[~a]))
    if abs(dcp) < min_dcp:
        sys.exit("FATAL: the two normal-sign groups differ by only %+.4f in mean cp at "
                 "mid-chord (threshold %.2f). Which one is the upper surface cannot be "
                 "established, and a guess would mislabel every reversed face. Refusing."
                 % (dcp, min_dcp))
    mult = 1.0 if dcp < 0.0 else -1.0          # upper = the MORE NEGATIVE mean cp
    if (dz > 0.0) != (mult > 0.0):
        sys.exit("FATAL: the cp criterion and the geometric z criterion disagree about which "
                 "normal sign is the upper surface (mean cp gap %+.4f, median z gap %+.6f). "
                 "Two independent criteria must agree before a flag is written. Refusing."
                 % (dcp, dz))
    return mult, dcp, dz, na, nb


def derive_q(p, cp, tol=1e-4):
    """q and p_inf from the case's OWN p and cp, with the fit residual ASSERTED.

    Returns (q, p_inf, max_rel, p99_rel). Halts if the relation is not linear to `tol`,
    because a cp written against a different reference than p would otherwise hand back a q
    that silently rescales every Cf in the report.

    THE TOLERANCE IS SET AGAINST THE DATA'S PRECISION, NOT PICKED. The solver writes this
    VTK as ASCII at about six significant figures, so a relative residual of order 1e-6 is
    TEXT ROUND-OFF and nothing else. My first value here was 1e-6 and it HALTED a perfectly
    good file at 5.8e-6, which is the always-fires species of broken gate: a threshold
    below the noise floor of its own input cannot pass. 1e-4 sits ~20x above the observed
    rounding and ~1e4x below any real reference mismatch, which would show up at order 1
    rather than at 1e-6. Both the max and the 99th percentile are returned so the margin is
    visible rather than asserted.
    """
    ok = np.isfinite(p) & np.isfinite(cp)
    if ok.sum() < 1000:
        sys.exit("FATAL: too few finite (p, cp) pairs to derive q")
    if np.ptp(cp[ok]) < 1e-9:
        sys.exit("FATAL: cp is constant over the surface; cannot derive q")
    # CENTRED COVARIANCE, NOT polyfit OR lstsq, AND THE REASON IS MEASURED.
    #   lstsq(rcond=None) raises on HPC12's numpy 1.12 ("must be real number, not NoneType").
    #   polyfit RUNS there and returns a subtly WRONG slope: on the compressible surfaces,
    #   where p ~ 2.3e4 Pa and cp ~ -0.3, it reported a fit whose residual was 7.0e-3 of q
    #   while the true relation holds to 9.3e-6. Regressing uncentred data with that dynamic
    #   range loses the precision in the normal equations, and the failure PRESENTED AS BAD
    #   DATA: this gate refused two good compressible cases and named the file as the
    #   culprit. A numerical defect in the checker is indistinguishable from the defect it
    #   was written to catch, which is why the fix is an estimator with no conditioning
    #   problem rather than a looser tolerance.
    # q = cov(cp, p) / var(cp); p_inf = mean(p) - q * mean(cp). Three lines, exact for a
    # linear relation, and identical on both numpy versions.
    c0 = cp[ok] - cp[ok].mean()
    p0 = p[ok] - p[ok].mean()
    q = float((c0 * p0).sum() / (c0 * c0).sum())
    p_inf = float(p[ok].mean() - q * cp[ok].mean())
    resid = np.abs(p[ok] - (q * cp[ok] + p_inf))
    scale = max(abs(q), 1e-30)
    rel = float(resid.max() / scale)
    p99 = float(np.percentile(resid, 99.0) / scale)
    if rel > tol:
        sys.exit("FATAL: p vs cp is not linear (max residual %.3e of q, p99 %.3e). The two "
                 "fields are not on the same reference; refusing to derive a normalisation."
                 % (rel, p99))
    return float(q), float(p_inf), rel, p99


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vtk", required=True, help="archived wingSurface.vtk[.gz]")
    ap.add_argument("--case", required=True)
    ap.add_argument("--eta", type=float, nargs="+",
                    default=[0.20, 0.40, 0.60, 0.75, 0.85, 0.95])
    ap.add_argument("--band", type=float, default=0.004,
                    help="half-width of the cut as a fraction of semispan")
    ap.add_argument("--semi", type=float, default=SEMI_SPAN_DEFAULT)
    ap.add_argument("--case-dir", type=Path, default=None,
                    help="case root holding postProcessing/forceCoeffs1; "
                         "default: three levels above the VTK")
    ap.add_argument("--out", type=Path, required=True)
    # OFF BY DEFAULT so the 2026-09-15 report extraction stays byte-reproducible. When set,
    # the CSV gains a trailing `surface` column written from the CALIBRATED FACET NORMAL,
    # which is the column split_section_surfaces.py would otherwise add on the laptop from
    # the geometric branch splitter. Downstream needs no change: read_section_csv already
    # takes columns by name and already understands `surface`.
    ap.add_argument("--surface-flag", action="store_true",
                    help="write the upper/lower flag from calibrated facet normals")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    print("reading %s" % a.vtk, flush=True)
    P, _area, F, NRM = read_surface_vtk(a.vtk, with_normals=a.surface_flag)
    print("  %d faces, fields: %s" % (len(P), ", ".join(sorted(F))), flush=True)
    for need in ("p", "cp", "wallShearStress"):
        if need not in F:
            sys.exit("FATAL: %s carries no '%s' field" % (a.vtk, need))

    q, p_inf, rel, p99 = derive_q(F["p"], F["cp"])
    print("  DERIVED FROM THIS FILE: q = %.6g, p_inf = %.6g "
          "(fit residual: max %.2e of q, p99 %.2e)" % (q, p_inf, rel, p99), flush=True)

    cp = F["cp"]
    tau = F["wallShearStress"]
    cf = np.linalg.norm(tau, axis=1) / q
    yp = F.get("yPlus", np.full(len(P), np.nan))

    # <case>/postProcessing/argusWingSurface/<time>/wingSurface.vtk -> parents[3] is <case>
    case_dir = a.case_dir or Path(a.vtk).resolve().parents[3]
    dhat = read_dragdir(case_dir)
    print("  freestream dragDir %s (alpha %.6f deg), read from forceCoeffs.dat"
          % (np.array2string(dhat, precision=9), np.degrees(np.arctan2(dhat[2], dhat[0]))),
          flush=True)
    tau_s_raw = tau @ dhat

    # PASS 1: cut every band, so the sign can be calibrated on a POOLED reference region
    # before anything is written. Calibrating per station would let two stations of one
    # case disagree about which way is downstream.
    bands, summary = [], []
    for eta in a.eta:
        y0 = eta * a.semi
        m = np.abs(P[:, 1] - y0) <= a.band * a.semi
        if m.sum() < 50:
            print("  eta %.2f: only %d faces in the band, SKIPPED" % (eta, m.sum()))
            summary.append({"eta": eta, "n": int(m.sum()), "status": "too-few-faces"})
            continue
        x = P[m, 0]
        c = x.max() - x.min()
        if c <= 0:
            summary.append({"eta": eta, "n": int(m.sum()), "status": "degenerate-chord"})
            continue
        bands.append((eta, m, (x - x.min()) / c, (P[m, 2] - P[m, 2].mean()) / c, c))

    if not bands:
        sys.exit("FATAL: no station produced a band")
    eta_f = np.concatenate([np.full(b[2].shape, b[0]) for b in bands])
    xc_f = np.concatenate([b[2] for b in bands])
    ts_f = np.concatenate([tau_s_raw[b[1]] for b in bands])
    sign, frac, nref = calibrate_shear_sign(ts_f, eta_f, xc_f)
    print("  SHEAR SIGN calibrated on %d inboard mid-chord faces: %.1f%% one sign -> "
          "multiplier %+.0f (attached is positive)" % (nref, 100.0 * frac, sign), flush=True)
    cf_s = sign * tau_s_raw / q

    # THE UPPER/LOWER FLAG, calibrated on the SAME pooled band faces the shear sign used, so
    # the two calibrations cannot disagree about which faces they were established on.
    up = None
    surf_note = None
    if a.surface_flag:
        if NRM is None:
            sys.exit("FATAL: --surface-flag was set but the reader returned no normals")
        z_f = np.concatenate([b[3] for b in bands])
        nz_f = np.concatenate([NRM[b[1], 2] for b in bands])
        cpb_f = np.concatenate([cp[b[1]] for b in bands])
        ref = (xc_f >= 0.35) & (xc_f <= 0.65)
        mult, dcp, dz, na, nb = calibrate_surface_flag(nz_f[ref], cpb_f[ref], z_f[ref])
        up = (mult * NRM[:, 2]) >= 0.0
        surf_note = ("calibrated facet normal (Newell per face); upper = the normal-sign "
                     "group with the more negative mean cp at x/c 0.35-0.65, multiplier "
                     "%+.0f; mean cp gap %+.4f and median z/c gap %+.6f agree on %d mid-chord "
                     "band faces (%d/%d by normal sign)" % (mult, dcp, dz, na + nb, na, nb))
        print("  SURFACE FLAG calibrated on %d mid-chord band faces (%d/%d by normal sign): "
              "mean cp gap %+.4f, median z/c gap %+.6f -> multiplier %+.0f"
              % (na + nb, na, nb, dcp, dz, mult), flush=True)

    # PASS 2: write
    for eta, m, xc, zc, c in bands:
        cpm, cfm, cfs, ypm = cp[m], cf[m], cf_s[m], yp[m]
        o = np.argsort(xc)
        rev = float((cfs < 0).mean())

        f = a.out / ("%s_eta%03d.csv" % (a.case, round(eta * 100)))
        with open(f, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["# case", a.case])
            w.writerow(["# eta", "%.4f" % eta])
            w.writerow(["# local_chord_m", "%.6f" % c])
            w.writerow(["# q_from_file", "%.6g" % q])
            w.writerow(["# p_inf_from_file", "%.6g" % p_inf])
            w.writerow(["# dragDir", "%.9f %.9f %.9f" % tuple(dhat)])
            w.writerow(["# shear_sign_multiplier", "%+.0f" % sign])
            w.writerow(["# n_faces_in_band", int(m.sum())])
            w.writerow(["# reverse_flow_fraction", "%.6f" % rev])
            w.writerow(["# note",
                        ("raw band, chord-ordered; the `surface` column is written HERE from "
                         "calibrated facet normals, so split_section_surfaces.py must not "
                         "re-split it. " if up is not None else
                         "raw band, chord-ordered; split upper/lower at plot time. ") +
                        "cf is |tau|/q; cf_s is the SIGNED freestream-aligned "
                        "component, negative = reverse flow"])
            if up is not None:
                w.writerow(["# surface_split", surf_note])
                upm = up[m]
                w.writerow(["x_over_c", "z_over_c", "cp", "cf", "cf_s", "yplus", "surface"])
                for i in o:
                    w.writerow(["%.6f" % xc[i], "%.6f" % zc[i], "%.6f" % cpm[i],
                                "%.8f" % cfm[i], "%+.8f" % cfs[i], "%.4f" % ypm[i],
                                "upper" if upm[i] else "lower"])
            else:
                w.writerow(["x_over_c", "z_over_c", "cp", "cf", "cf_s", "yplus"])
                for i in o:
                    w.writerow(["%.6f" % xc[i], "%.6f" % zc[i], "%.6f" % cpm[i],
                                "%.8f" % cfm[i], "%+.8f" % cfs[i], "%.4f" % ypm[i]])
        print("  eta %.2f: %d faces, reverse flow %.3f%% -> %s"
              % (eta, m.sum(), 100.0 * rev, f.name), flush=True)
        summary.append({"eta": eta, "n": int(m.sum()), "chord_m": c, "status": "ok",
                        "reverse_fraction": rev})

    # EVERY REQUESTED STATION IS ACCOUNTED FOR, not just the ones that worked (GEO-089:
    # produced plus deliberately-skipped-with-reason must equal all requested).
    ok = sum(1 for s in summary if s["status"] == "ok")
    print("STATIONS: %d requested, %d written, %d skipped"
          % (len(a.eta), ok, len(a.eta) - ok))
    for s in summary:
        if s["status"] != "ok":
            print("  SKIPPED eta %.2f: %s (%d faces)" % (s["eta"], s["status"], s["n"]))
    if ok == 0:
        sys.exit("FATAL: no station produced a section")


if __name__ == "__main__":
    main()
