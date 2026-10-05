#!/usr/bin/env python3
"""plot_section_slices.py: overlay Cp and Cf sections from several cases at one span station.

Consumes the small CSVs written by section_slices_from_surface.py on the cluster, so this
runs on the laptop with no cluster access. That split is deliberate and is
wing_sections.py's pattern: the inputs are ~600 MB ASCII VTKs and the outputs are a few
hundred kB.

WHAT IT DRAWS. One figure per station: -Cp on the upper axes, Cf on the lower, x/c along the
bottom, one colour per case, upper and lower surface as solid and dashed. The aerofoil
outline is drawn under the Cp axes so the pressure features can be read against the shape
that produced them, which was requested on 2026-08-12 because the aerofoil was not visible in
the earlier image, and is the reason the field-delta contour panels were replaced.

THE UPPER/LOWER SPLIT USES THE SHARED SPLITTER, airfoil_surfaces.split_branches, imported
and never reimplemented. D056 records that three implementations were tried and rejected,
and that the obvious one (sign of z) is WRONG on these aft-cambered supercritical sections:
32 of 87 stored points at defining station 4 have z > 0 on the LOWER surface, from x/c 0.78
onward, which is exactly where aft-camber morphing does its work. The artefact and the
signal occupy the same place, so it does not look like a bug, it looks like a result.

FRAME (D068). Every curve carries the q its own case was normalised on, read from the CSV
header, and the script REFUSES to overlay cases whose q differs unless --allow-mixed-q is
given with a reason. A wall-resolved case at q = 578 (kinematic, Condition CR) and a
compressible one at q = 11292.8 Pa are both correctly normalised and are NOT comparable
curves: same symbol, different operating point, which is frame-rule instance 8 in picture
form. The refusal names both values rather than saying "mixed".
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from airfoil_surfaces import split_branches  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


def read_section_csv(path):
    """-> (meta, x/c, z/c, cp, cf, cf_s, yplus, is_upper). Header comments carry the frame.

    COLUMNS ARE TAKEN BY NAME, NOT BY POSITION, and that is the whole point of this
    function's second version. The first read every row with np.array(rows, dtype=float),
    which is only valid while every column is numeric. split_section_surfaces.py added a
    trailing STRING column, `surface`, on 2026-09-15, and this reader (2026-09-12) died on
    it with "could not convert string to float: 'upper'" for every one of the 55 files.

    The positional unpack underneath was NOT wrong: with the new schema, columns 0..5 are
    still x, z, cp, cf, cf_s, yplus. Only the float cast was. Keying on the header row fixes
    the cast and makes a future column append a no-op instead of a breakage.

    is_upper is None when the file carries no `surface` column, so callers can tell
    "no flag recorded" from "flagged lower".
    """
    meta, rows, cols = {}, [], None
    with open(path, newline="") as fh:
        for r in csv.reader(fh):
            if not r:
                continue
            if r[0].startswith("#"):
                meta[r[0].lstrip("# ").strip()] = r[1] if len(r) > 1 else ""
                continue
            if r[0] == "x_over_c":
                cols = [s.strip() for s in r]
                continue
            rows.append(r)
    if not rows:
        return None
    if cols is None:
        # Headerless legacy files: the documented positional order. The 5-column form has
        # NO cf_s and its column 4 is yplus, so the two widths cannot share one list.
        cols = (["x_over_c", "z_over_c", "cp", "cf", "yplus"] if len(rows[0]) == 5
                else ["x_over_c", "z_over_c", "cp", "cf", "cf_s", "yplus"])[:len(rows[0])]
    idx = {n: i for i, n in enumerate(cols)}

    def num(name):
        # cf_s is absent from CSVs written before the signed-shear pass; carry NaN rather
        # than silently plotting the magnitude in its place, which would draw an
        # always-attached curve and hide exactly what the column was added to show.
        if name not in idx:
            return np.full(len(rows), np.nan)
        j = idx[name]
        return np.array([float(r[j]) for r in rows])

    is_upper = None
    if "surface" in idx:
        j = idx["surface"]
        is_upper = np.array([r[j].strip().lower() == "upper" for r in rows])
    return (meta, num("x_over_c"), num("z_over_c"), num("cp"), num("cf"),
            num("cf_s"), num("yplus"), is_upper)


def branch_masks(xc, zc):
    """Boolean masks for (upper, lower) recovered from the shared splitter.

    split_branches returns COORDINATES, not indices, and Cp/Cf have to travel with the
    points. It sorts by x and returns an ordered subsequence of that sort, so the indices
    come back by a two-pointer merge against the same sort. The SPLIT itself is not
    reimplemented; only the bookkeeping the shared function does not hand back.
    """
    br = split_branches(xc, zc)
    if br is None:
        return None
    o = np.argsort(np.asarray(xc, float), kind="stable")
    xs, zs = np.asarray(xc, float)[o], np.asarray(zc, float)[o]
    out = []
    for bx, bz in br:
        m = np.zeros(len(xc), dtype=bool)
        j = 0
        for k in range(len(bx)):
            while j < len(xs) and not (xs[j] == bx[k] and zs[j] == bz[k]):
                j += 1
            if j >= len(xs):
                return None
            m[o[j]] = True
            j += 1
        out.append(m)
    return out[0], out[1]


def masks_for(xc, zc, is_upper):
    """(upper, lower) masks: the CSV's own flag where it has one, else recompute.

    PREFER THE RECORDED FLAG. split_section_surfaces.py wrote it with the same shared
    airfoil_surfaces.split_branches, and it HALTS unless its index recovery reproduces
    split_branches element for element, so the column is a proved partition rather than a
    convenience (D052 item 4: a within-file proof beats an external arbiter).

    Recomputing here is the fallback, not the default, because branch_masks would re-run
    the splitter on a 60k-row RAW BAND. split_branches was written for a clean section
    outline; a band is a finite-width slab of a SWEPT, TAPERED wing, so it carries many
    points at the same x/c with genuine spanwise z scatter. That is the shape the flag was
    added to handle, and re-deriving it at plot time would throw that work away.
    """
    if is_upper is not None:
        return is_upper, ~is_upper
    return branch_masks(xc, zc)


OUTLINE_DIR = REPO / "results/wing_outlines_of12"
OUTLINE_LETTERS = set("BCFHMW")   # the six registered geometries, per CASE_PROVENANCE.json


def _outline_frame(eta_int):
    """The ONE frame every outline at this station is re-expressed in: the baseline's.

    Any geometry would do as the datum; the baseline is chosen because it is the reference
    the whole report differences against, so the frame is named rather than incidental.
    """
    p = OUTLINE_DIR / "outline_B.json"
    if not p.exists():
        return None
    for s in json.loads(p.read_text())["sections"]:
        if round(s["eta"] * 100) == eta_int and "z_mean_m" in s:
            return s["x_min_m"], s["chord_m"], s["z_mean_m"]
    return None


def _outline_for(case, eta_int):
    """The geometry's own outline at this station, or None to fall back to the band.

    THE CASE NAME CARRIES THE GEOMETRY IN ITS PREFIX. CMPB_trim, SWB_trim and LCB_trim are
    the SAME WING at three operating points, so the geometry is the last character of the
    token before the first underscore. Checked against all twelve case names on disk, not
    assumed: CMP{B,C,F,H,M,W}_trim and SW{B,C,F,H,M,W}_trim all resolve correctly.

    IT ASSERTS ITS PRECONDITIONS INSTEAD OF FALLING THROUGH QUIETLY (GEO-092: an exemption
    creates a fall-through path, taken by exactly the inputs it was not written for). A
    prefix that yields an unregistered letter, or an artifact carrying no such station, is
    REPORTED and then declined. Only a genuinely absent artifact falls back silently, and
    that one is safe because the fallback IS the documented previous behaviour rather than
    a guess about what was wanted.
    """
    letter = case.split("_")[0][-1:].upper()
    if letter not in OUTLINE_LETTERS:
        print("  %s: prefix resolves to %r, not a registered geometry; band outline used"
              % (case, letter))
        return None
    p = OUTLINE_DIR / ("outline_%s.json" % letter)
    if not p.exists():
        return None
    rec = json.loads(p.read_text())
    for s in rec["sections"]:
        if round(s["eta"] * 100) != eta_int:
            continue
        # SIX GEOMETRIES ON ONE AXIS IS A COMPARISON, SO THEY NEED ONE FRAME (D068). Each
        # outline arrives normalised on ITS OWN chord and about ITS OWN mean z, which is the
        # right convention for a single section and the wrong one for an overlay: aft camber
        # moves the mean, so a purely aft change leaks forward as a rigid vertical offset.
        #
        # MEASURED, and this is why it is fixed rather than noted. The six wings are
        # identical forward of the hinge, so their outlines must coincide there. Un-datumed
        # they spread 0.000e+00 at eta 0.20 and 0.40, where the morph is zero, but 0.43,
        # 1.01 and 0.44 per cent of chord at eta 0.60, 0.80 and 0.95, where it is not. At
        # eta 0.80 that is 1.01 %c of pure artifact against 2.36 %c of real aft difference,
        # 43 per cent of the signal, worst between M (largest camber command) and W (none).
        #
        # The chords themselves agree to 0.039 %c, so the x re-datum is a formality; it is
        # applied anyway because a frame conversion that silently omits one axis is the
        # thing this rule exists to prevent.
        if "z_mean_m" not in s:
            print("  %s: outline_%s.json predates the frame fields; regenerate with "
                  "scripts/wing_outlines_of12.py --all. Band outline used." % (case, letter))
            return None
        fr = _outline_frame(eta_int)
        if fr is None:
            print("  %s: no baseline frame at eta%03d to re-datum onto; band outline used"
                  % (case, eta_int))
            return None
        xb0, cb, zb0 = fr
        c0, x0, z0 = s["chord_m"], s["x_min_m"], s["z_mean_m"]

        def _re(xs, zs):
            return (np.asarray(xs) * c0 + x0 - xb0) / cb, (np.asarray(zs) * c0 + z0 - zb0) / cb

        return [_re(s["upper_xc"], s["upper_zc"]), _re(s["lower_xc"], s["lower_zc"])]
    print("  %s: outline_%s.json carries no eta%03d station (has %s); band outline used"
          % (case, letter, eta_int, ", ".join("%.2f" % s["eta"] for s in rec["sections"])))
    return None


def bin_branch(xc, vals, nbin):
    """Mean of each val on a uniform x/c grid over one branch, plus the spread.

    THE RAW BAND IS NOT A SECTION AND PLOTTING IT RAW SAYS SO LOUDLY. The cut is a slab of
    finite width (default +/-0.004 semispan, about 7 mm), and this wing is SWEPT and
    TAPERED, so two faces at the same x/c sit at different y and genuinely carry different
    Cp. Drawn unbinned that spanwise variation reads as mesh noise and buries a real
    aft-camber difference of order 0.1 in Cp under a scatter of similar size.

    Returns (centres, {name: mean}, {name: (lo, hi)}) with empty bins dropped. The spread is
    the 10th-90th percentile and is RETURNED rather than discarded, because it is the honest
    width of the slab and the reader should be able to see it if they want it.
    """
    xc = np.asarray(xc, float)
    edges = np.linspace(0.0, 1.0, nbin + 1)
    idx = np.clip(np.digitize(xc, edges) - 1, 0, nbin - 1)
    keep = np.array([np.any(idx == b) for b in range(nbin)])
    centres = 0.5 * (edges[:-1] + edges[1:])[keep]
    means, spreads = {}, {}
    for name, v in vals.items():
        v = np.asarray(v, float)
        m = np.full(nbin, np.nan)
        lo = np.full(nbin, np.nan)
        hi = np.full(nbin, np.nan)
        for b in range(nbin):
            s = v[idx == b]
            if s.size:
                m[b] = s.mean()
                lo[b], hi[b] = np.percentile(s, 10.0), np.percentile(s, 90.0)
        means[name] = m[keep]
        spreads[name] = (lo[keep], hi[keep])
    return centres, means, spreads


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True, help="directory of section CSVs")
    ap.add_argument("--bins", type=int, default=240,
                    help="chordwise bins per branch; 0 plots the raw band")
    ap.add_argument("--show-spread", action="store_true",
                    help="shade the 10-90 percentile across the slab width")
    ap.add_argument("--eta", type=int, required=True, help="station, as eta*100 (e.g. 75)")
    ap.add_argument("--cases", nargs="+", required=True)
    ap.add_argument("--allow-mixed-q", default="",
                    help="reason for overlaying cases on different q; refused without one")
    ap.add_argument("--out", type=Path, default=REPO / "results/section_figs")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    data, missing = {}, []
    for c in a.cases:
        f = a.dir / ("%s_eta%03d.csv" % (c, a.eta))
        if not f.exists():
            missing.append(c)
            continue
        d = read_section_csv(f)
        if d is None:
            missing.append(c)
            continue
        data[c] = d
    if not data:
        sys.exit("FATAL: none of the requested cases has a station at eta%03d" % a.eta)

    qs = {c: float(d[0]["q_from_file"]) for c, d in data.items()}
    if len(set("%.6g" % v for v in qs.values())) > 1 and not a.allow_mixed_q:
        lines = ["FATAL: these cases are normalised on DIFFERENT q and are not comparable:"]
        for c, v in sorted(qs.items(), key=lambda kv: kv[1]):
            lines.append("    %-14s q = %-12.6g p_inf = %s" % (c, v, data[c][0]["p_inf_from_file"]))
        lines.append("  A kinematic q (0.5*U^2) and a dimensional one (0.5*rho*U^2, Pa) are")
        lines.append("  both correct and describe different operating points (D076/D068).")
        lines.append("  Pass --allow-mixed-q '<reason>' if the overlay is genuinely intended.")
        sys.exit("\n".join(lines))

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1]})
    colours = plt.cm.tab10(np.linspace(0, 1, 10))
    shown = []
    for i, (c, (meta, xc, zc, cp, cf, cfs, yp, isup)) in enumerate(sorted(data.items())):
        bm = masks_for(xc, zc, isup)
        if bm is None:
            print("  %s: branch split failed, SKIPPED" % c)
            missing.append(c)
            continue
        col = colours[i % 10]
        for m, style, lab in ((bm[0], "-", "upper"), (bm[1], "--", "lower")):
            if a.bins > 0:
                use = cfs[m] if np.isfinite(cfs[m]).any() else cf[m]
                xb, mu, sp = bin_branch(xc[m], {"cp": cp[m], "cf": use}, a.bins)
                ax0.plot(xb, -mu["cp"], style, color=col, lw=1.3,
                         label=("%s %s" % (c, lab)))
                ax1.plot(xb, mu["cf"], style, color=col, lw=1.3)
                if a.show_spread:
                    ax0.fill_between(xb, -sp["cp"][0], -sp["cp"][1], color=col, alpha=0.15,
                                     lw=0)
                    ax1.fill_between(xb, sp["cf"][0], sp["cf"][1], color=col, alpha=0.15,
                                     lw=0)
            else:
                o = np.argsort(xc[m])
                use = cfs[m] if np.isfinite(cfs[m]).any() else cf[m]
                ax0.plot(xc[m][o], -cp[m][o], style, color=col, lw=0.8,
                         label=("%s %s" % (c, lab)))
                ax1.plot(xc[m][o], use[o], style, color=col, lw=0.8)
        shown.append(c)

    # ONE OUTLINE PER GEOMETRY, NOT ONE OUTLINE FOR ALL OF THEM (review, 2026-09-15: one
    # aerofoil drawn under curves from different geometries is misleading). The first version drew data[shown[0]] alone, in neutral grey and unlabelled,
    # so six curves sat over ONE geometry's section with nothing saying whose. Which case
    # that was fell out of alphabetical order, not out of a decision.
    #
    # THAT IS THE WORST POSSIBLE OMISSION ON THIS FIGURE, because the aft section shape is
    # precisely what the geometries differ in: at eta 0.80 the trailing edge moves up to
    # 7.29 mm on a 219 mm local chord, 3.3 per cent of chord, which is plainly visible at
    # this scale. The single outline was hiding the very difference the reader is asked to
    # judge, while implying all six shared one shape.
    #
    # ONE SHARED SCALE, fixed from the first case and reused for the rest. Scaling each
    # outline to its own z range would make genuinely different shapes plot identically,
    # which is the same defect that per-panel autoscaling is in the displacement figure.
    if shown:
        ylo, yhi = ax0.get_ylim()
        order = sorted(data)
        outl = {c: _outline_for(c, a.eta) for c in sorted(shown)}
        sc = None
        for c in sorted(shown):
            meta, xc, zc, cp, cf, cfs, yp, isup = data[c]
            if sc is None:
                # SCALE FROM WHAT IS ACTUALLY DRAWN, not from the sample. When the outlines
                # are in use the band's z range is no longer the plotted one, and scaling
                # the geometry by the sample's extent would let the band back into the
                # figure through the one number nobody would think to check.
                ref = np.concatenate([z for _, z in outl[c]]) if outl.get(c) else zc
                sc = 0.35 * (yhi - ylo) / max(np.ptp(ref), 1e-9)
            col = colours[order.index(c) % 10]
            # THE OUTLINE COMES FROM THE GEOMETRY, NOT FROM THE FLOW SAMPLE. Review,
            # 2026-09-16: an outline taken from the sample made the leading and trailing edges
            # wiggle on a smooth geometry, a false impression.
            #
            # He is right and the old drawing earned it. It binned the CFD BAND's own z/c:
            # the band is +/-0.004 semispan, 7.32 mm, at EVERY station, so on a swept and
            # tapered wing it pools faces whose true z genuinely differs, and the pooled mean
            # wiggles wherever the surface curves hardest. Measured identically on both
            # sources (60 bins, median, second difference of z/c), the band is 4x rougher at
            # the root leading edge, 19x at mid-span and 59x at eta 0.95 -- and it DEGRADES
            # SEVENFOLD OUTBOARD while a planar cut stays flat, which is the signature of a
            # fixed-width slab sampling a shrinking chord.
            #
            # scripts/wing_outlines_of12.py cuts the PLACED STL, the same surface
            # snappyHexMesh was given, on an exact plane with no width and therefore nothing
            # to average. Falls back to the old band drawing only when no outline exists, so
            # a missing artifact degrades to the previous behaviour rather than to a blank.
            ol = outl.get(c)
            if ol is not None:
                for xb, zz in ol:
                    ax0.plot(xb, ylo + 0.05 * (yhi - ylo) + np.asarray(zz) * sc,
                             color=col, lw=0.9, alpha=0.8, zorder=0)
                continue
            bm = masks_for(xc, zc, isup)
            if bm is None:
                continue
            for m in bm:
                if a.bins > 0:
                    xb, mu, _ = bin_branch(xc[m], {"z": zc[m]}, a.bins)
                    zz = mu["z"]
                else:
                    o = np.argsort(xc[m]); xb, zz = xc[m][o], zc[m][o]
                ax0.plot(xb, ylo + 0.05 * (yhi - ylo) + zz * sc,
                         color=col, lw=0.8, alpha=0.75, zorder=0)

    q0 = list(qs.values())[0]
    ax0.set_ylabel(r"$-C_p$")
    ax0.set_title("eta = %.2f    q = %.6g %s"
                  % (a.eta / 100.0, q0,
                     "(kinematic, m^2/s^2)" if q0 < 1e3 else "(Pa)"))
    ax0.legend(fontsize=7, ncol=2)
    ax0.grid(alpha=0.3)
    signed = any(np.isfinite(d[5]).any() for d in data.values())
    ax1.set_ylabel(r"$C_{f,s}$ (freestream-aligned)" if signed else r"$C_f$  $|\tau|/q$")
    ax1.set_xlabel("x/c  (local chord at this station)")
    ax1.grid(alpha=0.3)
    ax1.axhline(0.0, color="k", lw=0.8)
    fig.tight_layout()
    f = a.out / ("sections_eta%03d.png" % a.eta)
    fig.savefig(f, dpi=150)
    print("wrote %s" % f)
    # COVERAGE, STATED (GEO-087): what was asked for, what was drawn, what was not.
    if signed:
        for c in shown:
            rf = data[c][0].get("reverse_flow_fraction")
            if rf is not None:
                print("  %-14s reverse flow in this band: %.3f%%" % (c, 100.0 * float(rf)))
    else:
        print("  NOTE: no cf_s column; the Cf axis is the MAGNITUDE and cannot show "
              "reverse flow. Re-extract to get the signed component.")
    print("CASES: %d requested, %d drawn, %d missing%s"
          % (len(a.cases), len(shown), len(missing),
             (": " + ", ".join(sorted(set(missing)))) if missing else ""))


if __name__ == "__main__":
    main()
