#!/usr/bin/env python3
"""render_delta_maps_of12.py: baseline plus five candidate DELTAS, upper and lower, one figure.

    python3 scripts/render_delta_maps_of12.py --quantity cf --condition CR
    python3 scripts/render_delta_maps_of12.py --all            # 2 quantities x 2 conditions
    python3 scripts/render_delta_maps_of12.py --null           # baseline against itself

WHAT THIS DRAWS, AND WHY IN THIS SHAPE. Specification, 2026-09-15: the baseline, then each
other geometry as delta Cp and delta Cf on the upper and lower surfaces, 12 panels per
variable and one figure per flight condition. So: 6 columns (B, then C F H M W) x 2 rows (UPPER, LOWER) = 12 panels, one figure
per quantity per condition, four figures in total.

COLUMN ONE IS ABSOLUTE, THE OTHER FIVE ARE DIFFERENCES, and they carry SEPARATE colour bars.
A candidate minus baseline is two orders of magnitude smaller than the field itself, so one
shared scale across all six would render the five deltas uniformly grey. The five deltas DO
share one scale with each other, because comparing them is the entire point; per-panel
autoscaling would make five different magnitudes look identical, which is the defect
render_concept_3d.py records for the displacement figure.

THE AXES ARE (eta, x/c) AND THAT IS A RECTANGLE, NOT A PLANFORM OUTLINE. The maps are binned
on the common (eta, x/c) grid, which is the only support on which a difference between two
separately meshed wings is defined at all (export_wing_vtk.py states the same and resolves it
the same way). Painting the rectangle back onto the true swept planform is a further step and
is not done here, because the bin grid is what the numbers actually live on and drawing them
on an outline they were not computed on would be a claim about resolution this does not have.
Span runs left to right, root at left. x increases UPWARD, so the leading edge is at the bottom, the
same reading order as the displacement figure.

EMPTY BINS ARE LEFT EMPTY. A bin no face landed in is nan and is drawn as blank, never as
zero. Zero is a measurement here: it means the candidate matched the baseline there.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
MAPS = REPO / "results/wing_maps_of12"

COND = {"CR": ("SWB", ["SWC", "SWF", "SWH", "SWM", "SWW"], "Condition CR, $M$ 0.10"),
        "EC": ("CMPB", ["CMPC", "CMPF", "CMPH", "CMPM", "CMPW"], "early cruise, $M$ 0.78"),
        "LC": ("LCB", ["LCC", "LCF", "LCH", "LCM", "LCW"], "late cruise, $M$ 0.78")}
LETTER = {"SWC": "C", "SWF": "F", "SWH": "H", "SWM": "M", "SWW": "W",
          "CMPC": "C", "CMPF": "F", "CMPH": "H", "CMPM": "M", "CMPW": "W",
          "LCC": "C", "LCF": "F", "LCH": "H", "LCM": "M", "LCW": "W"}
# THE FOURTH FIELD IS `signed`: DOES THIS QUANTITY'S ZERO MEAN ANYTHING? For C_p and C_f it
# does not, and a sequential scale across the data's own range is the honest choice. For
# C_{f,s} it is the whole point. Negative means the wall shear opposes the freestream, so the
# ZERO CROSSING IS THE SEPARATION LINE, and percentile limits would put the colour midpoint at
# the median of a wing that is almost entirely attached: separated and attached flow would come
# out nearly the same colour and the one feature the field exists to show would be invisible.
# A signed quantity therefore gets a diverging map on limits mirrored about zero, so the colour
# boundary and the physical boundary are the same line.
QTY = {"cp":   (r"$C_p$",     r"$\Delta C_p$",     "viridis", False),
       "cf":   (r"$C_f$",     r"$\Delta C_f$",     "viridis", False),
       "cf_s": (r"$C_{f,s}$", r"$\Delta C_{f,s}$", "RdBu_r",  True)}

# MINIMUM REVERSED FRACTION BEFORE A SEPARATION CONTOUR IS DRAWN. Measured, not chosen by
# feel: at Condition CR, where nothing separates, every geometry still carries 4 to 6 negative
# bins of 76,797 at the root junction and the blunt trailing edge, and at early cruise the
# attached cases carry 0 to 3. The real signal is 1318 bins for C and 1716 for M. Noise and
# signal are more than two orders of magnitude apart, so any cut between them serves; this one
# sits 13x above the worst noise and 17x below the smallest signal.
SEP_MIN = 0.001


def performance(tag):
    """-> (Cd counts, dCd vs that condition's baseline, L/D, alpha) for one trim case.

    THE FIGURE SHOULD CARRY THE COUNTS, NOT LEAVE THEM TO BE LOOKED UP. Specification, 2026-09-16:
    the figures must state each wing's performance, its C_L and C_D. A field
    map that does not say what the field cost is decorative: the reader sees where the
    pressure moved but not whether it was worth anything. dCd is the number the report ranks
    on, so it goes in the panel title.

    Read from data/rans_forces_current.json, the SAME harvest the performance tables use, so
    a figure and a table can never disagree. C_L is deliberately NOT per panel: every case is
    trimmed to the common target, so it is a property of the comparison and belongs in the
    caption once, not repeated six times as though it varied.
    """
    import functools
    return _perf_table()[tag]


@__import__("functools").lru_cache(maxsize=1)
def _perf_table():
    p = REPO / "docs/report/all_geometry_2026-09-15/data/rans_forces_current.json"
    cases = json.loads(p.read_text())["cases"]
    out = {}
    for cond, (base_tag, cands, _lab) in COND.items():
        b = cases["%s_trim" % base_tag]
        bcd = float(b["Cd_counts"])
        for t in [base_tag] + cands:
            v = cases["%s_trim" % t]
            out[t] = (float(v["Cd_counts"]), float(v["Cd_counts"]) - bcd,
                      float(v["L_over_D"]), float(v["alpha_deg"]))
    return out


def _cr_negative_bins():
    """(min, max) negative upper-surface C_f,s bins over the six Condition CR maps. COMPUTED: this
    range was typed into the note as "4 to 6" and went stale when the CR maps were rebuilt."""
    n = []
    for t in [COND["CR"][0]] + COND["CR"][1]:
        m, _ = load(t)
        u = m["cf_s_upper"]
        n.append(int(np.sum(u[np.isfinite(u)] < 0)))
    return min(n), max(n)


def load(tag):
    npz = MAPS / ("maps_%s_trim.npz" % tag)
    js = MAPS / ("maps_%s_trim.json" % tag)
    if not npz.exists():
        raise SystemExit("%s: no map at %s; run scripts/wing_maps_of12.py --all" % (tag, npz))
    return dict(np.load(npz)), json.loads(js.read_text())


def _qdef(quantity, prov):
    """The definition of the QUANTITY ACTUALLY PLOTTED, or a stated absence. Never silence.

    The caption used to print prov.get("cf_definition", "") whatever was on the page, so the
    C_{f,s} figure described itself with the MAGNITUDE's formula, |wallShearStress| / q, on a
    figure whose whole point is the SIGN. Worse, the "" default FAILS OPEN: a map binned before
    a definition existed yields a caption that simply omits it, and an omission reads as
    "nothing to say" rather than "not recorded". A figure that misstates how its own field was
    computed is the same defect as one that misstates its axes, which this script also had.

    cp needs no stamp: it is the solver's own cp, carried through unmodified.
    """
    if quantity == "cp":
        return "$C_p$ as written by the solver."
    v = prov.get("%s_definition" % quantity)
    if v:
        return "%s = %s." % (quantity, v)
    return ("%s: DEFINITION NOT RECORDED by the export that produced this map. Re-run "
            "scripts/of12_surface_to_export.py to stamp it." % quantity)


def robust_sym(vals, pct=99.5):
    """Symmetric limit from a percentile of |delta|, ignoring empty bins."""
    v = np.concatenate([np.abs(x[np.isfinite(x)]).ravel() for x in vals if np.isfinite(x).any()])
    if v.size == 0:
        return 1.0
    m = float(np.percentile(v, pct))
    return m if m > 0 else float(v.max() or 1.0)


def build(quantity, condition, out, dpi=170):
    base_tag, cand_tags, cond_label = COND[condition]
    label, dlabel, cmap, signed = QTY[quantity]
    B, bprov = load(base_tag)
    cands = [(t,) + load(t) for t in cand_tags]

    # FRAME CHECK BEFORE ARITHMETIC (D068): every panel must be the same operating point and
    # the same grid, or the difference is between unlike things.
    for t, _g, p in cands:
        for k in ("condition", "mach", "U_inf", "n_eta", "n_xc", "semispan_registered_m"):
            if p.get(k) != bprov.get(k):
                raise SystemExit("%s: %s is %r against the baseline's %r; refusing to difference"
                                 % (t, k, p.get(k), bprov.get(k)))

    # COVERAGE CHECK BEFORE ARITHMETIC, beside the frame check and for the same reason. A map
    # binned before a quantity existed simply lacks the key, and the delta loop below then dies
    # on a bare KeyError. It DOES stop, which is lucky, but it stops because a dict lookup blew
    # up rather than because anything checked: A CORRECT ANSWER REACHED BY AN INVALID ROUTE IS
    # A DEFECT (GEO-080 item 5). The cost is concrete under --all, where cp and cf succeed for
    # both conditions and cf_s would crash partway, leaving four of six figures written and a
    # traceback where a diagnosis belongs. Naming the cases turns it into an instruction.
    need = ["%s_%s" % (quantity, s) for s in ("upper", "lower")]
    lacking = [t for t, g, _p in [(base_tag, B, bprov)] + cands
               if any(k not in g for k in need)]
    if lacking:
        raise SystemExit(
            "%s/%s: %d of %d cases carry no %s. Re-bin them with "
            "`python3 scripts/wing_maps_of12.py --case <CASE>`; refusing to draw a panel set "
            "that silently omits cases. Missing: %s"
            % (condition, quantity, len(lacking), 1 + len(cands), quantity, ", ".join(lacking)))

    deltas = {}
    for t, g, _p in cands:
        for surf in ("upper", "lower"):
            k = "%s_%s" % (quantity, surf)
            deltas[(t, surf)] = g[k] - B[k]          # nan propagates where either is empty
    vmax = robust_sym(list(deltas.values()))

    bvals = [B["%s_%s" % (quantity, s)] for s in ("upper", "lower")]
    bfin = np.concatenate([x[np.isfinite(x)].ravel() for x in bvals])
    blo, bhi = np.percentile(bfin, [0.2, 99.8])
    if signed:
        bhi = float(max(abs(float(blo)), abs(float(bhi))))
        blo = -bhi
        # STATE THE SEPARATED FRACTION RATHER THAN LEAVING IT TO THE EYE. The map shows WHERE;
        # this says HOW MUCH, on the same faces, so the figure and any sentence written about
        # it cannot drift apart.
        for s in ("upper", "lower"):
            v = B["%s_%s" % (quantity, s)]
            g = np.isfinite(v)
            print("    baseline %s: %.3f%% of %d filled bins reversed"
                  % (s, 100.0 * float((v[g] < 0).mean()), int(g.sum())))

    # RESERVE THE BOTTOM STRIP AND THE COLOURBAR COLUMN EXPLICITLY, rather than letting a
    # layout engine place them. Under constrained_layout the caption was drawn at y=0.004 with
    # NO space reserved for it, so four lines of it ran straight across the lower row of panels
    # and over the eta axis labels; and a colorbar attached to axes[:, 0] stole its space from
    # between columns one and two, breaking the panel spacing. Both are the same failure: the
    # engine lays out the axes and then text and bars are added on top of a full canvas.
    # Fixed margins are deterministic and version-independent, which matters more here than
    # elegance.
    # THE GRID IS SIZED FROM THE CANDIDATE LIST, NOT HARDCODED AT SIX. The first version
    # wrote plt.subplots(2, 6) and then indexed cand_tags[col - 1] for col up to 5, so it
    # assumed exactly five candidates and raised IndexError on any other number. Both
    # completed conditions happen to have five, so it would never have failed in production
    # today -- and would have failed the moment it was pointed at LATE CRUISE, which has 4 of
    # 24 cases converged and will spend a long time with fewer than five candidates. Found by
    # running the frame check's POSITIVE control with a one-candidate list, which is the half
    # of a negative-control pair that is easy to skip.
    ncol = 1 + len(cand_tags)
    fig, axes = plt.subplots(2, ncol, figsize=(3.17 * ncol, 7.2), squeeze=False)
    fig.subplots_adjust(left=0.055, right=0.895, top=0.895, bottom=0.185,
                        wspace=0.09, hspace=0.12)
    # FINITE COORDINATES EVERYWHERE, INCLUDING WHERE THE VALUE IS MASKED. pcolormesh refuses
    # non-finite X or Y outright, and x_m/y_m carry NaN wherever a bin caught no face: off the
    # wing near the root fairing, and in the last x/c column at stations where the blunt
    # trailing edge falls inside the bin. Masking the FIELD is not enough, because the mesh
    # coordinates are a separate argument and are read even under a mask.
    #
    # THE FILL IS EXACT, NOT INTERPOLATED. Within one eta row the section is planar, so x_m is
    # LINEAR in x/c: x_m = x_LE(eta) + (x/c) * chord(eta). Fitting that line to the row's own
    # filled bins reproduces them and extends the row without inventing a shape. y_m is
    # constant along a row by construction, the row being a station, so its median fills it.
    # Measured centres are kept wherever they exist; the fit only ever fills gaps.
    def _planform_mesh(M):
        X, Y = np.array(M["x_m"], float), np.array(M["y_m"], float)
        n_eta, n_xc = X.shape
        xcc = (np.arange(n_xc) + 0.5) / n_xc
        sl = np.full(n_eta, np.nan)
        ic = np.full(n_eta, np.nan)
        yv = np.full(n_eta, np.nan)
        for i in range(n_eta):
            g = np.isfinite(X[i])
            if g.sum() >= 2:
                sl[i], ic[i] = np.polyfit(xcc[g], X[i][g], 1)
            gy = np.isfinite(Y[i])
            if gy.any():
                yv[i] = np.median(Y[i][gy])
        idx = np.arange(n_eta)
        for v in (sl, ic, yv):
            f = np.isfinite(v)
            if not f.any():
                raise SystemExit("planform mesh: no row has finite coordinates; refusing")
            v[~f] = np.interp(idx[~f], idx[f], v[f])
        # THE FIT IS USED EVERYWHERE, NOT ONLY IN THE GAPS, and that is the more correct
        # choice rather than a convenience. pcolormesh wants CELL CENTRES; x_m is the
        # AREA-WEIGHTED CENTROID of the faces that fell in the bin, which is a different
        # quantity and is pulled toward wherever area concentrates, worst at the leading edge
        # where curvature packs faces into a short run of x. Mixing the two made the
        # coordinate array non-monotonic, which matplotlib warned about and which would have
        # let it derive cell edges wrongly. The bin's true centre is x_LE(eta) + (x/c)*chord,
        # exactly the line fitted above, so the fit IS the geometry and the centroids were
        # never the right input.
        Xf = ic[:, None] + sl[:, None] * xcc[None, :]
        Yf = np.repeat(yv[:, None], n_xc, axis=1)
        for nm, A in (("x", Xf), ("y", Yf)):
            d = np.diff(A, axis=1 if nm == "x" else 0)
            if d.size and not (np.all(d >= 0) or np.all(d <= 0)):
                raise SystemExit("planform mesh: %s not monotonic after fit; refusing" % nm)
        if not (np.isfinite(Xf).all() and np.isfinite(Yf).all()):
            raise SystemExit("planform mesh: fill left non-finite coordinates; refusing")
        return Xf, Yf

    XmF, YmF = _planform_mesh(B)
    gmap = {t: g for t, g, _p in cands}          # candidate tag -> its own ABSOLUTE fields
    imb = imd = None
    for row, surf in enumerate(("upper", "lower")):
        for col in range(ncol):
            ax = axes[row, col]
            # DRAWN ON THE TRUE PLANFORM, NOT ON A UNIT SQUARE. Specification, 2026-09-16: on
            # the wing surface, not on a rectangle.
            # An (eta, x/c) image is a rectangle and silently implies a rectangular wing to
            # anyone who does not read the axis labels. The maps now carry x_m and y_m, the
            # physical centroid of every bin, so pcolormesh draws the real swept, tapered
            # outline. Axes match the concept-displacement figure exactly, because that is
            # the figure this was asked to look like: SPAN horizontal with the root at left,
            # CHORD vertical, equal aspect, and nothing stretched to fill a box.
            Xm, Ym = XmF, YmF
            if col == 0:
                m = B["%s_%s" % (quantity, surf)]
                imb = ax.pcolormesh(Ym, Xm, np.ma.masked_invalid(m), cmap=cmap,
                                    vmin=blo, vmax=bhi, shading="nearest")
                cd, dcd, ld, _al = performance(base_tag)
                title = "B  %s\n$C_D$ %.1f ct, $L/D$ %.2f" % (label, cd, ld)
            else:
                t = cand_tags[col - 1]
                m = deltas[(t, surf)]
                imd = ax.pcolormesh(Ym, Xm, np.ma.masked_invalid(m), cmap="coolwarm",
                                    vmin=-vmax, vmax=vmax, shading="nearest")
                _cd, dcd, ld, _al = performance(t)
                title = "%s  %s\n$\\Delta C_D$ %+.1f ct, $L/D$ %.2f" % (LETTER[t], dlabel, dcd, ld)
            # THE SEPARATION LINE, DRAWN WHERE IT ACTUALLY IS. A DELTA PANEL CANNOT SHOW
            # SEPARATION. Candidate minus baseline says the shear FELL, not that it REVERSED,
            # and telling those two apart is the single thing this campaign most needs from
            # these figures. F carries a visibly weakened boundary layer behind the same shock
            # as C and M, and reverses 0.004 per cent of its upper bins, which is the rigid
            # baseline's own value: a reader going by the blue delta patch alone would call F
            # separated and be wrong. This contour is that case's OWN C_{f,s} at zero, so
            # inside it is reversed flow and outside it is attached, whatever the delta
            # underneath happens to look like.
            if signed:
                own = (B if col == 0 else gmap[cand_tags[col - 1]])["%s_%s" % (quantity, surf)]
                fin = np.isfinite(own)
                # GATED ON SEP_MIN, BECAUSE "any negative bin" FIRES ON THE NOISE FLOOR. It
                # fired on all twelve panels at Condition CR, where nothing separates. Today
                # those contours are sub-pixel and invisible, which is LUCK RATHER THAN SAFETY:
                # printed larger they would read as separation, and a mark a reader trusts must
                # not be drawn around four cells of trailing-edge noise.
                frac = float((own[fin] < 0).mean()) if fin.any() else 0.0
                if frac >= SEP_MIN:
                    ax.contour(Ym, Xm, np.ma.masked_invalid(own), levels=[0.0],
                               colors="k", linewidths=0.7)
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            if row == 0:
                ax.set_title(title, fontsize=9)
            if col == 0:
                ax.set_ylabel(surf.upper(), fontsize=10)

    # Two bars, stacked on the right in their own axes: the absolute field and the shared
    # delta scale are different quantities and must never share a bar.
    cb_b = fig.add_axes([0.912, 0.560, 0.011, 0.320])
    cb_d = fig.add_axes([0.912, 0.200, 0.011, 0.320])
    fig.colorbar(imb, cax=cb_b).set_label(label + "  (baseline)", fontsize=9)
    fig.colorbar(imd, cax=cb_d).set_label(dlabel + "  (candidates)", fontsize=9)
    for cb in (cb_b, cb_d):
        cb.tick_params(labelsize=7)
    fig.suptitle("%s on the wing, baseline and each candidate minus baseline: %s"
                 % (label, cond_label), fontsize=13)

    nfilled = int(np.isfinite(B["%s_upper" % quantity]).sum())
    # THE THRESHOLD TRAVELS WITH THE FIGURE. A mark that is sometimes drawn and sometimes not
    # is unreadable unless the rule is stated, and an absent contour is itself a result here.
    # BREAKS ARE AUTHORED, NOT HOPED FOR. Matplotlib does not wrap figure text and no layout
    # engine polices it, so a run-on line leaves the page mid-word. This script's own concept
    # figure carries that scar already; the first version of this note reproduced it exactly,
    # truncating at "and its ABSENCE means the case does no".
    sepnote = ("" if not signed else
               "\nBlack contour is that case's own %s = 0, the separation line: inside it the "
               "flow is reversed, and its ABSENCE means the case does not separate.\n"
               "Drawn only where at least %.1f per cent of a surface's filled bins are "
               "negative: every case carries a few negative bins at the root junction and the "
               "blunt trailing edge\n(%d to %d of %d at Condition CR, where nothing separates), "
               "and a contour around those would read as separation."
               % ((label, 100 * SEP_MIN) + _cr_negative_bins() + (nfilled,)))
    fig.text(0.004, 0.004,
             "Binned on the common $(\\eta, x/c)$ grid, %d x %d, area-weighted; the two wings are "
             "meshed separately and share no cells, so this grid is the only support on which a "
             "difference exists.\n"
             "Column B is the ABSOLUTE field on its own scale; the five candidate columns are "
             "candidate minus baseline on ONE shared scale, $\\pm$%.4g at the 99.5th percentile, "
             "so their magnitudes are comparable by eye.\n"
             "Span left to right, root at left; $x$ increases upward, leading edge at the "
             "bottom. Blank bins are bins no face landed in, never zero. %s, each wing at its own "
             "incidence at the common target $C_L$, so this is at MATCHED LIFT.\n"
             "Upper %d filled bins of %d.\n%s\n"
             "Generated by scripts/render_delta_maps_of12.py from results/wing_maps_of12/."
             % (bprov["n_eta"], bprov["n_xc"], vmax, cond_label, nfilled,
                bprov["n_eta"] * bprov["n_xc"], _qdef(quantity, bprov) + sepnote),
             fontsize=7, va="bottom")

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print("  wrote %s  (delta scale +/-%.4g, %d filled bins)" % (out.name, vmax, nfilled))
    return vmax


def null_test():
    """Baseline differenced against itself must be EXACTLY zero in every filled bin."""
    bad = 0
    for cond in ("CR", "EC", "LC"):
        base_tag = COND[cond][0]
        try:
            B, _ = load(base_tag)
        except SystemExit as e:
            print("  SKIP %s: %s" % (cond, e)); continue
        for k, m in sorted(B.items()):
            d = m - m
            fin = np.isfinite(m)
            if not np.array_equal(np.isfinite(d), fin):
                print("  FAIL %s %s: filled set changed" % (cond, k)); bad += 1
            elif fin.any() and np.nanmax(np.abs(d[fin])) != 0.0:
                print("  FAIL %s %s: not exactly zero" % (cond, k)); bad += 1
            else:
                print("  PASS %s %-9s exactly 0 over %d filled bins" % (cond, k, int(fin.sum())))
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quantity", choices=sorted(QTY))
    ap.add_argument("--condition", choices=sorted(COND))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--null", action="store_true")
    ap.add_argument("--outdir", type=Path, default=REPO / "results/figures/delta_maps")
    a = ap.parse_args()

    if a.null:
        return null_test()
    jobs = ([(q, c) for q in sorted(QTY) for c in sorted(COND)] if a.all
            else [(a.quantity, a.condition)] if (a.quantity and a.condition) else [])
    if not jobs:
        return ap.error("give --quantity and --condition, or --all, or --null")
    ok, bad = 0, []
    for q, c in jobs:
        try:
            build(q, c, a.outdir / ("delta_map_%s_%s.png" % (c, q)))
            ok += 1
        except SystemExit as e:
            print("  %s %s REFUSED: %s" % (c, q, e)); bad.append("%s/%s" % (c, q))
    print("  PARTITION: %d requested = %d drawn + %d refused%s"
          % (len(jobs), ok, len(bad), (": " + ", ".join(bad)) if bad else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
