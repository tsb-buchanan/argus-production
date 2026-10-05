#!/usr/bin/env python3
"""
Cp section figures for the ARGUS all-geometry report (2026-09-15).

Inputs, all read from disk, nothing carried from memory:
  data/sections/<CASE>_trim_eta<NNN>.csv   surface band cuts (cp, cf_s, yplus, surface flag)
  data/sections/Q_PINF_TABLE.csv           per-case q, p_inf, alpha, dragDir
  data/rans_forces.json                    trim C_L, C_D, turbulence model, frame

Outputs, into this directory:
  fig_cp_sections_CR.pdf              Condition CR  (M 0.10), 6 geometries, 5 stations
  fig_cp_sections_earlycruise.pdf     early cruise  (M 0.78), 6 geometries, 5 stations
  fig_cp_sections_latecruise.pdf      late cruise   (M 0.78), 6 geometries, 5 stations
  fig_cp_shock_separation_CMPC.pdf    eta 0.80 early-cruise shock / separation detail
  cp_sections_measured.json           every number the figures display

Frame note (D068). cp here is the solver's own surface cp, normalised on each case's own
freestream dynamic pressure q read from that file's header. It carries NO reference area and
is NOT on the DSO or the TP-1580 basis. x/c is on each band's own local chord, also from the
file header. Force coefficients quoted in the annotations ARE on the DSO basis, half-model
(Aref 0.620462 m2), as recorded in rans_forces.json, and are labelled as such.

Operating points (D076). SW* is Condition CR, M 0.10, target C_L 0.428277635108.
CMP* is early cruise, M 0.78, target C_L 0.529297087450. They are DIFFERENT points and
nothing is differenced across them.

Validate-before-write: every figure is built, its numbers asserted, and only then written.
The "written" claim at the end is derived by stat-ing the files, not by the exit status.
"""

import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ----------------------------------------------------------------------------- paths
HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)
# sections_nflag: the upper/lower flag now comes from the calibrated facet normal written by
# the extractor, not from the geometric branch splitter, which mislabels faces where the two
# surfaces converge. data/sections/ is retained as the superseded 2026-09-15 set.
SECT = os.path.join(REPORT, "data", "sections_nflag")
RANS_JSON = os.path.join(REPORT, "data", "rans_forces.json")
QTAB = os.path.join(SECT, "Q_PINF_TABLE.csv")

# PER-CASE DATASET RESOLUTION. A case is served from data/sections_r2 where that exists and
# from sections_nflag otherwise, so the figures show the newer evaluation for the cases that
# have one. SECT and QTAB above remain the fallback and the module docstring's stated
# default; section_source decides per case. The q/dragDir row ALWAYS comes from the same
# dataset as the section CSV: the two datasets disagree on alpha by up to 0.0301 deg and
# dragDir enters cf_s directly, so crossing them is the D019 wind-axis defect. See
# fig/section_source.py for the measured deltas.
import section_source  # noqa: E402  (after the path constants it reads)

# ----------------------------------------------------------------------------- constants
GAMMA = 1.4  # stated because the sonic-Cp and Mach cross-checks depend on it
STATIONS = [20, 40, 60, 80, 95]
GEOMS = ["B", "C", "F", "H", "M", "W"]

# Letter only. Project convention: bare minimum labels on the figure; what each
# geometry IS belongs in the report text, not repeated in every legend.
GEOM_LABEL = {g: g for g in ["B", "C", "F", "H", "M", "W"]}

# Okabe-Ito, colour-blind safe; every geometry also has a distinct dash pattern so the
# figure survives greyscale printing.
STYLE = {
    "B": dict(color="#000000", ls="solid", lw=1.15),
    "C": dict(color="#D55E00", ls=(0, (4.0, 1.4)), lw=1.05),
    "F": dict(color="#0072B2", ls=(0, (5.0, 1.2, 1.0, 1.2)), lw=1.05),
    "H": dict(color="#009E73", ls=(0, (1.3, 1.3)), lw=1.15),
    "M": dict(color="#CC79A7", ls=(0, (6.0, 1.2, 1.0, 1.2, 1.0, 1.2)), lw=1.05),
    "W": dict(color="#E69F00", ls=(0, (2.6, 1.2)), lw=1.05),
}

# Binning of the surface band onto a common chordwise grid.
NBIN = 200
MINPTS = 8  # a bin with fewer points than this is left empty, never interpolated over
EDGES = np.linspace(0.0, 1.0, NBIN + 1)
CTR = 0.5 * (EDGES[:-1] + EDGES[1:])

# Shock flag. Stated here so the criterion travels with the number.
SHOCK_GRAD_THRESHOLD = 15.0   # d(cp)/d(x/c), per unit chord
SHOCK_WINDOW = (0.15, 0.95)   # excludes the LE suction peak and the TE recovery
SHOCK_HALFSPAN = 0.04         # chord fraction over which the cp rise is quoted

# Leading-edge cut on the DIFFERENCE panels only. Forward of this the section is steeply
# curved, the pooled span band and the upper/lower branch split both become ill-conditioned,
# and single-bin medians swing by more than the difference being measured. The size of that
# swing is measured, not asserted: see `le_roughness` in the JSON.
DIFF_XMIN = 0.03

# ----------------------------------------------------------------------------- style
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.0,
    "axes.labelsize": 8.0,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "legend.fontsize": 7.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "round",
    "pdf.fonttype": 42,          # embed TrueType, keep text as text (vector)
    "ps.fonttype": 42,
    "savefig.dpi": 600,
    "figure.dpi": 150,
})


# ----------------------------------------------------------------------------- io
def read_section(case, eta):
    """Read one band CSV. Returns (header dict, column dict). No column is invented."""
    try:
        sect_dir = section_source.resolve("%s_trim" % case)[0]
    except KeyError:
        return None, None
    path = os.path.join(sect_dir, "%s_trim_eta%03d.csv" % (case, eta))
    if not os.path.exists(path):
        return None, None
    hdr, cols, rows = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                key, _, val = line[2:].rstrip("\n").partition(",")
                hdr[key] = val
                continue
            if cols is None:
                cols = line.rstrip("\n").split(",")
                continue
            rows.append(line.rstrip("\n").split(","))
    arr = np.array(rows, dtype=object)
    out = {}
    for i, c in enumerate(cols):
        out[c] = arr[:, i].astype(str) if c == "surface" else arr[:, i].astype(float)
    hdr["_path"] = path
    return hdr, out


def _load_q_table(path):
    rows = {}
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split(",")
        for line in fh:
            v = line.rstrip("\n").split(",")
            d = dict(zip(head, v))
            rows[d["case"]] = d
    return rows


def read_q_table():
    """Per-case q, p_inf, alpha and dragDir, EACH ROW FROM THE SAME DATASET AS ITS SECTIONS.

    Not a merge with one table winning. The row is fetched through the same resolver the
    section CSVs go through, so a case can never end up with its coefficients from one
    evaluation and its wind axes from the other.
    """
    tables = {d: _load_q_table(os.path.join(d, "Q_PINF_TABLE.csv"))
              for d in section_source.DATASETS}
    rows = {}
    for case in sorted(set().union(*(set(t) for t in tables.values()))):
        try:
            qdir = os.path.dirname(section_source.resolve(case)[1])
        except KeyError:
            continue          # a q row with no section data is never plotted
        if case in tables.get(qdir, {}):
            rows[case] = dict(tables[qdir][case], _generation=section_source.resolve(case)[2])
    return rows


def bin_branch(x, y, minpts=MINPTS):
    """Median, q25 and q75 of y in each x/c bin. Empty bins stay NaN."""
    idx = np.digitize(x, EDGES) - 1
    med = np.full(NBIN, np.nan)
    q25 = np.full(NBIN, np.nan)
    q75 = np.full(NBIN, np.nan)
    cnt = np.zeros(NBIN, dtype=int)
    order = np.argsort(idx, kind="stable")
    idx_s, y_s = idx[order], y[order]
    bounds = np.searchsorted(idx_s, np.arange(NBIN + 1))
    for i in range(NBIN):
        a, b = bounds[i], bounds[i + 1]
        cnt[i] = b - a
        if cnt[i] >= minpts:
            v = y_s[a:b]
            med[i] = np.median(v)
            q25[i] = np.percentile(v, 25)
            q75[i] = np.percentile(v, 75)
    return med, q25, q75, cnt


def frac_negative(x, y, minpts=MINPTS):
    idx = np.digitize(x, EDGES) - 1
    frac = np.full(NBIN, np.nan)
    order = np.argsort(idx, kind="stable")
    idx_s, y_s = idx[order], y[order]
    bounds = np.searchsorted(idx_s, np.arange(NBIN + 1))
    for i in range(NBIN):
        a, b = bounds[i], bounds[i + 1]
        if (b - a) >= minpts:
            frac[i] = float(np.mean(y_s[a:b] < 0.0))
    return frac


def cp_star(mach):
    """Sonic pressure coefficient at freestream Mach, gamma = GAMMA."""
    t = ((1.0 + 0.5 * (GAMMA - 1.0) * mach ** 2) / (1.0 + 0.5 * (GAMMA - 1.0))) ** (
        GAMMA / (GAMMA - 1.0)
    )
    return 2.0 / (GAMMA * mach ** 2) * (t - 1.0)


def contiguous_runs(mask):
    """[(i0, i1_inclusive), ...] for each run of True in mask."""
    runs, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        elif not m and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask) - 1))
    return runs


# ----------------------------------------------------------------------------- load
def build_panel_data():
    """Load every available (case, station). Returns data dict plus a closed partition."""
    rans = json.load(open(RANS_JSON))
    qtab = read_q_table()

    conditions = {
        "CR": dict(
            prefix="SW",
            trimkey="condition_CR",
            title="Condition CR, M 0.10",
            long="Condition CR (M 0.10), target $C_L$ = 0.428277635108",
        ),
        "EC": dict(
            prefix="CMP",
            trimkey="early_cruise",
            title="Early cruise, M 0.78",
            long="Early cruise (M 0.78), target $C_L$ = 0.529297087450",
        ),
        "LC": dict(
            prefix="LC",
            trimkey="late_cruise",
            title="Late cruise, M 0.78",
            long="Late cruise (M 0.78), target $C_L$ = 0.544114800210",
        ),
    }

    data = {}
    present, absent = [], []
    for ck, cinfo in conditions.items():
        trims = rans["trims"][cinfo["trimkey"]]
        frame = trims["coefficient_frame"]
        cond = dict(
            key=ck,
            title=cinfo["title"],
            long=cinfo["long"],
            target_CL=trims["target_CL"],
            mach_nominal=trims["mach"],
            Aref_m2=frame.get("Aref_m2"),
            lRef_m=frame.get("lRef_m"),
            basis=frame.get("basis_identified"),
            turb=frame.get("turbulence_models_present"),
            geoms={},
        )
        for g in GEOMS:
            case = "%s%s" % (cinfo["prefix"], g)
            gtrim = trims["geometries"].get(g, {})
            stations = {}
            for eta in STATIONS:
                hdr, col = read_section(case, eta)
                if hdr is None:
                    absent.append((case, eta))
                    continue
                present.append((case, eta))
                q = float(hdr["q_from_file"])
                pinf = float(hdr["p_inf_from_file"])
                # Mach derived within-file from q and p_inf; zero p_inf (kinematic,
                # Condition CR) makes it undefined, which is correct, not an error.
                mach = float(np.sqrt(2.0 * q / (GAMMA * pinf))) if pinf > 1.0 else None
                up = col["surface"] == "upper"
                lo = ~up
                d = dict(
                    eta=float(hdr["eta"]),
                    chord_m=float(hdr["local_chord_m"]),
                    q=q,
                    p_inf=pinf,
                    mach_from_file=mach,
                    n_faces=int(hdr["n_faces_in_band"]),
                    rev_frac_file=float(hdr["reverse_flow_fraction"]),
                    path=hdr["_path"],
                )
                for nm, m in (("up", up), ("lo", lo)):
                    med, q25, q75, cnt = bin_branch(col["x_over_c"][m], col["cp"][m])
                    d["cp_" + nm] = med
                    d["cp_%s_q25" % nm] = q25
                    d["cp_%s_q75" % nm] = q75
                    d["n_" + nm] = cnt
                cfs_med, _, _, _ = bin_branch(col["x_over_c"][up], col["cf_s"][up])
                d["cfs_up"] = cfs_med
                # UPPER-SURFACE SUCTION PEAK BY BOTH ESTIMATORS, EACH NAMED.
                # The appendix prose quotes the RAW FACE EXTREMUM over the pooled span
                # band; the figure beside it draws the 200-BIN MEDIAN. They differ by
                # 0.02 to 0.05 in cp, always in the same direction, because a median over
                # a bin cannot reach the single most extreme face in it. NEITHER IS
                # WRONG AND THE ESTIMATOR WAS SIMPLY UNSTATED, so a reader could not find
                # the quoted peak anywhere in the plot. That is the frame rule (D065)
                # applied to an ESTIMATOR rather than to a normalisation: the number is
                # not reportable without the thing that produced it. Both are written so
                # a caption can quote the one consistent with what it sits beside.
                _rx, _rc = col["x_over_c"][up], col["cp"][up]
                _bm = d["cp_up"]
                _bok = ~np.isnan(_bm)
                _ri = int(np.argmin(_rc)) if _rc.size else None
                _bi = int(np.argmin(_bm[_bok])) if _bok.any() else None
                d["suction_peak"] = dict(
                    raw_face_cp_min=(float(_rc[_ri]) if _ri is not None else None),
                    raw_face_x_over_c=(float(_rx[_ri]) if _ri is not None else None),
                    binned_median_cp_min=(float(_bm[_bok][_bi]) if _bi is not None else None),
                    binned_median_x_over_c=(float(CTR[_bok][_bi]) if _bi is not None else None),
                    n_bins=int(NBIN),
                    estimator_note=("raw_face_* is the single most negative cp among the "
                                    "faces in the pooled band; binned_median_* is the "
                                    "minimum of the %d-bin medians the figure draws. "
                                    "Quote the one that matches the artefact." % NBIN),
                )
                # Wall treatment read off the solution, not off a recipe name. y+ ~ 0.3 is
                # an integrated near-wall solution; y+ ~ 40 is a wall-function reconstruction.
                d["yplus_median_upper"] = float(np.median(col["yplus"][up]))
                d["yplus_p99_upper"] = float(np.percentile(col["yplus"][up], 99))
                d["revfrac_up"] = frac_negative(col["x_over_c"][up], col["cf_s"][up])
                rev_u = up & (col["cf_s"] < 0.0)
                d["rev_up_n"] = int(rev_u.sum())
                d["rev_up_xmin"] = float(col["x_over_c"][rev_u].min()) if rev_u.any() else None
                d["rev_up_xmax"] = float(col["x_over_c"][rev_u].max()) if rev_u.any() else None
                stations[eta] = d
            cond["geoms"][g] = dict(
                case=case,
                stations=stations,
                Cl=gtrim.get("Cl"),
                alpha_deg=gtrim.get("alpha_deg"),
                Cd_counts=gtrim.get("Cd_counts"),
                dCd_counts=gtrim.get("delta_Cd_counts_vs_baseline"),
                trim_ok=gtrim.get("trim_within_tolerance"),
                turb=gtrim.get("turbulence_model"),
                alpha_from_section=(
                    float(qtab["%s_trim" % case]["alpha_deg"]) if "%s_trim" % case in qtab else None
                ),
                surface_time=(
                    qtab["%s_trim" % case]["surface_time"] if "%s_trim" % case in qtab else None
                ),
            )
        data[ck] = cond

    # GEO-089: the claim over the set must partition the set.
    expected = len(conditions) * len(GEOMS) * len(STATIONS)
    assert len(present) + len(absent) == expected, (
        "partition failure: %d present + %d absent != %d expected"
        % (len(present), len(absent), expected)
    )
    return data, present, absent


# ----------------------------------------------------------------------------- analysis
def detect_shock(st):
    """Max d(cp)/d(x/c) on the binned upper branch inside SHOCK_WINDOW.

    Returns a dict always; `flagged` is True only if the gradient clears the threshold AND
    the flow just upstream is supersonic by the cp* test. A gate must say what it saw, so
    the unflagged maximum is returned too, never dropped.
    """
    ok = ~np.isnan(st["cp_up"])
    xs, ys = CTR[ok], st["cp_up"][ok]
    if xs.size < 12:
        return dict(flagged=False, reason="too few populated bins", n_bins=int(xs.size))
    grad = np.gradient(ys, xs)
    win = (xs >= SHOCK_WINDOW[0]) & (xs <= SHOCK_WINDOW[1])
    if not win.any():
        return dict(flagged=False, reason="empty search window", n_bins=int(xs.size))
    i = int(np.argmax(grad[win]))
    x_sh = float(xs[win][i])
    g_sh = float(grad[win][i])
    cp_us = float(np.interp(x_sh - SHOCK_HALFSPAN, xs, ys))
    cp_ds = float(np.interp(x_sh + SHOCK_HALFSPAN, xs, ys))
    cps = cp_star(st["mach_from_file"]) if st["mach_from_file"] else None
    supersonic_upstream = (cps is not None) and (cp_us < cps)
    n_sup = int(np.sum(ys < cps)) if cps is not None else 0
    return dict(
        flagged=bool(g_sh >= SHOCK_GRAD_THRESHOLD and supersonic_upstream),
        x_over_c=x_sh,
        max_dcp_dx=g_sh,
        cp_upstream=cp_us,
        cp_downstream=cp_ds,
        delta_cp=cp_ds - cp_us,
        cp_star=cps,
        supersonic_upstream=bool(supersonic_upstream),
        n_supersonic_bins=n_sup,
        n_bins=int(xs.size),
    )


def detect_separation(st):
    """Upper-surface reverse flow from the binned median of the signed cf.

    cf_s is the freestream-aligned wall-shear component; the extractor calibrated its sign
    and negative means reverse flow. The bin medians are used, not single faces.
    """
    med = st["cfs_up"]
    ok = ~np.isnan(med)
    neg = ok & (med < 0.0)
    runs = contiguous_runs(neg)
    runs = [r for r in runs if (r[1] - r[0] + 1) >= 3]  # >= 3 bins = 0.015 c, not a single bin
    out = dict(
        any=bool(runs),
        band_reverse_fraction_all_faces=st["rev_frac_file"],
        n_reversed_upper_faces=st["rev_up_n"],
        raw_upper_x_min=st["rev_up_xmin"],
        raw_upper_x_max=st["rev_up_xmax"],
        runs=[],
    )
    for a, b in runs:
        seg = slice(a, b + 1)
        out["runs"].append(dict(
            x_start=float(EDGES[a]),
            x_end=float(EDGES[b + 1]),
            min_cf_s=float(np.nanmin(med[seg])),
            reaches_trailing_edge=bool(EDGES[b + 1] > 0.985),
        ))
    return out


def diff_vs_baseline(cond, g, eta, branch):
    """cp(candidate) - cp(baseline B) on the common x/c bin grid, one surface branch.

    Both sides are binned identically, on the same band half-width at the same eta, with
    local chords differing by at most 9.2e-4 relative. The band-pooling defect is therefore
    common to both sides of the difference and largely cancels (D071). It cancels only
    approximately, since the two geometries need not carry the same spanwise cp gradient
    inside the band, so the residual is bounded rather than asserted away.
    """
    sb = cond["geoms"]["B"]["stations"].get(eta)
    sg = cond["geoms"][g]["stations"].get(eta)
    if sb is None or sg is None:
        return None
    return sg["cp_" + branch] - sb["cp_" + branch]


def le_roughness(d):
    """Median |second difference| of a delta-cp curve, inside and outside the LE cut.

    A smooth physical signal has a small second difference; bin-level noise does not. This
    is the measurement that justifies DIFF_XMIN, so the cut is a number and not a taste.
    """
    out = {}
    for name, m in (("forward_of_cut", CTR < DIFF_XMIN),
                    ("aft_of_cut", (CTR >= DIFF_XMIN) & (CTR <= 0.90))):
        v = np.where(m, d, np.nan)
        s = v[1:-1] * 2.0 - v[:-2] - v[2:]
        s = s[~np.isnan(s)]
        out[name] = float(np.median(np.abs(s))) if s.size else None
    return out


# ----------------------------------------------------------------------------- figures
def panel_title(cond, eta):
    st = None
    for g in GEOMS:
        st = cond["geoms"][g]["stations"].get(eta)
        if st is not None:
            break
    c = st["chord_m"]
    return r"$\eta$ = %.2f    $c$ = %.3f m" % (eta / 100.0, c)


def finish_cp_axis(ax, ylim, show_y):
    ax.set_ylim(ylim[1], ylim[0])  # inverted: -Cp upward, the convention
    ax.set_xlim(-0.02, 1.02)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.grid(True, lw=0.3, color="0.85", zorder=0)
    ax.set_axisbelow(True)
    if not show_y:
        ax.tick_params(labelleft=False)


def make_condition_figure(cond, outpath, geoms_present, measured):
    """3 rows x 5 stations: cp, then candidate-minus-baseline on each surface branch."""
    is_cruise = cond["key"] == "EC"
    ncol = len(STATIONS)
    fig, axes = plt.subplots(
        # BUILT AT THE WIDTH IT IS PRINTED AT, for the reason recorded in make_cf.py: the
        # report's textwidth is 170 mm = 6.69 in and this is included at \linewidth, so a
        # 7.3 in figure was downscaled to 0.92 and every font with it. Placed two to a row
        # at 0.49\linewidth, as it was until today, the scale was 0.66 and nothing in it
        # could be read.
        3, ncol, figsize=(6.69, 6.64),
        gridspec_kw=dict(hspace=0.30, wspace=0.28,
                         left=0.085, right=0.995, top=0.895, bottom=0.135),
    )

    # Shared cp limits across the stations of this condition, from the data.
    allcp = []
    for g in geoms_present:
        for eta in STATIONS:
            st = cond["geoms"][g]["stations"].get(eta)
            if st is None:
                continue
            allcp.append(st["cp_up"])
            allcp.append(st["cp_lo"])
    lo = float(np.nanmin(np.concatenate(allcp)))
    hi = float(np.nanmax(np.concatenate(allcp)))
    pad = 0.06 * (hi - lo)
    cplim = (lo - pad, hi + pad)

    cands = [g for g in geoms_present if g != "B"]

    for j, eta in enumerate(STATIONS):
        # ---------------- row 0: cp
        ax = axes[0, j]
        st_b = cond["geoms"]["B"]["stations"].get(eta)
        if st_b is not None:
            for br in ("up", "lo"):
                m = ~np.isnan(st_b["cp_%s_q25" % br])
                ax.fill_between(CTR[m], st_b["cp_%s_q25" % br][m], st_b["cp_%s_q75" % br][m],
                                color="0.60", alpha=0.35, lw=0, zorder=1)
        if is_cruise and st_b is not None and st_b["mach_from_file"]:
            ax.axhline(cp_star(st_b["mach_from_file"]), color="0.35", lw=0.7,
                       ls=(0, (1.0, 1.6)), zorder=2)
        for g in geoms_present:
            st = cond["geoms"][g]["stations"].get(eta)
            if st is None:
                continue
            for br in ("up", "lo"):
                ax.plot(CTR, st["cp_" + br], zorder=4, **STYLE[g])
        # separation extent first, so the shading sits under the curves
        for g in geoms_present:
            sp = measured["separation"].get((cond["key"], g, eta))
            if sp and sp["any"]:
                for r in sp["runs"]:
                    ax.axvspan(r["x_start"], r["x_end"], color=STYLE[g]["color"],
                               alpha=0.13, lw=0, zorder=1)
        # shock marks, stacked so the labels of co-located shocks do not collide
        span = cplim[1] - cplim[0]
        flagged = [(g, measured["shock"][(cond["key"], g, eta)]) for g in geoms_present
                   if measured["shock"].get((cond["key"], g, eta), {}).get("flagged")]
        flagged.sort(key=lambda t: t[1]["x_over_c"])
        for k, (g, sh) in enumerate(flagged):
            y_tip = cplim[0] + (0.10 + 0.085 * k) * span
            ax.annotate(
                "", xy=(sh["x_over_c"], y_tip),
                xytext=(sh["x_over_c"], cplim[0] + 0.012 * span),
                arrowprops=dict(arrowstyle="-|>", lw=0.75, color=STYLE[g]["color"],
                                mutation_scale=5.5,
                                linestyle="solid", shrinkA=0, shrinkB=0),
                zorder=6, annotation_clip=False)
            ax.text(sh["x_over_c"] + 0.015, y_tip, "%s  %.2f" % (g, sh["x_over_c"]),
                    ha="left", va="center", fontsize=6.0,
                    color=STYLE[g]["color"], zorder=7,
                    bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6))
        finish_cp_axis(ax, cplim, show_y=(j == 0))
        ax.set_title(panel_title(cond, eta), pad=3.5)
        if j == 0:
            ax.set_ylabel(r"$c_p$   (inverted)")

        # ---------------- rows 1, 2: differences, per branch
        for r, (br, brname) in enumerate((("up", "upper"), ("lo", "lower")), start=1):
            axd = axes[r, j]
            axd.axhline(0.0, color="0.45", lw=0.6, zorder=2)
            keep = CTR >= DIFF_XMIN
            peak = 0.0
            for g in cands:
                d = diff_vs_baseline(cond, g, eta, br)
                if d is None:
                    continue
                axd.plot(CTR[keep], d[keep], zorder=4, **STYLE[g])
                w = keep & (CTR < 0.98)
                if np.any(~np.isnan(d[w])):
                    peak = max(peak, float(np.nanmax(np.abs(d[w]))))
            axd.set_xlim(DIFF_XMIN - 0.04, 1.02)
            axd.set_xticks([0.0, 0.5, 1.0])
            lim = max(1.15 * peak, 0.02)
            axd.set_ylim(-lim, lim)
            axd.grid(True, lw=0.3, color="0.85", zorder=0)
            axd.set_axisbelow(True)
            axd.tick_params(axis="y", labelsize=6.2)
            if j == 0:
                axd.set_ylabel(r"$\Delta c_p$, %s surface" % brname)
            if r == 2:
                axd.set_xlabel(r"$x/c$")

    # ---------------- header text
    turb = cond["turb"][0] if cond["turb"] else "NOT FOUND"
    # Bare title. The frame this used to state is kept in the figure's JSON (_frame)
    # and in fig/FIGURE_NOTES.txt, so the frame still travels with the number.
    head = cond["long"].split(",")[0]
    fig.suptitle("", fontsize=7.6, y=0.987, ha="center", va="top", linespacing=1.45)  # caption carries it

    # ---------------- legend
    handles = [Line2D([], [], label=GEOM_LABEL[g],
                      **dict(STYLE[g], lw=max(STYLE[g]["lw"], 1.45)))
               for g in geoms_present]
    handles.append(Patch(facecolor="0.60", alpha=0.35, label="B band IQR"))
    if is_cruise:
        handles.append(Line2D([], [], color="0.35", lw=0.8, ls=(0, (1.0, 1.6)),
                              label=r"$c_p^{*}$ at M %.4f" % cond_mach(cond)))
    sep_geoms = sorted({g for (ck, g, _), v in measured["separation"].items()
                        if ck == cond["key"] and v["any"]})
    for g in sep_geoms:
        handles.append(Patch(facecolor=STYLE[g]["color"], alpha=0.13,
                             label="%s upper-surface reverse flow" % g))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 0.004), handlelength=3.1, columnspacing=1.4,
               borderaxespad=0.0)
    return fig


def cond_mach(cond):
    for g in GEOMS:
        for eta in STATIONS:
            st = cond["geoms"][g]["stations"].get(eta)
            if st is not None and st["mach_from_file"]:
                return st["mach_from_file"]
    return float("nan")


def make_shock_detail_figure(cond, outpath, geoms_present, measured, eta=80):
    """Early-cruise eta 0.80: cp, signed cf and reverse-flow fraction on the upper surface."""
    fig, axes = plt.subplots(
        1, 3, figsize=(7.3, 2.75),
        gridspec_kw=dict(wspace=0.30, left=0.075, right=0.995, top=0.815, bottom=0.315),
    )
    xlim = (0.45, 1.01)
    sel = (CTR >= xlim[0]) & (CTR <= xlim[1])

    # (a) cp, upper surface only
    ax = axes[0]
    cps = None
    vals = []
    for g in geoms_present:
        st = cond["geoms"][g]["stations"].get(eta)
        if st is None:
            continue
        ax.plot(CTR, st["cp_up"], zorder=4, **STYLE[g])
        vals.append(st["cp_up"][sel])
        if st["mach_from_file"]:
            cps = cp_star(st["mach_from_file"])
    if cps is not None:
        ax.axhline(cps, color="0.35", lw=0.7, ls=(0, (1.0, 1.6)), zorder=2)
    v = np.concatenate(vals)
    lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
    pad = 0.08 * (hi - lo)
    ax.set_ylim(hi + pad, lo - pad)
    ax.set_xlim(*xlim)
    ax.set_ylabel(r"$c_p$, upper surface  (inverted)")
    ax.set_title("(a)", pad=3.5)
    flagged = sorted(
        ((g, measured["shock"][(cond["key"], g, eta)]) for g in geoms_present
         if measured["shock"].get((cond["key"], g, eta), {}).get("flagged")),
        key=lambda t: t[1]["x_over_c"])
    ylo, yhi = ax.get_ylim()
    for k, (g, sh) in enumerate(flagged):
        ax.axvline(sh["x_over_c"], color=STYLE[g]["color"], lw=0.9,
                   ls=(0, (2.0, 1.6)), zorder=3)
        ax.text(sh["x_over_c"], ylo + (0.955 - 0.075 * k) * (yhi - ylo),
                " %s %.2f" % (g, sh["x_over_c"]), ha="left", va="center",
                fontsize=6.2, color=STYLE[g]["color"], zorder=7,
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6))

    # (b) signed cf along the freestream direction
    ax = axes[1]
    ax.axhline(0.0, color="0.45", lw=0.6, zorder=2)
    for g in geoms_present:
        st = cond["geoms"][g]["stations"].get(eta)
        if st is None:
            continue
        ax.plot(CTR, st["cfs_up"] * 1e3, zorder=4, **STYLE[g])
    sp = measured["separation"].get((cond["key"], "C", eta))
    if sp and sp["any"]:
        for r in sp["runs"]:
            ax.axvspan(r["x_start"], r["x_end"], color=STYLE["C"]["color"],
                       alpha=0.13, lw=0, zorder=1)
    ax.set_xlim(*xlim)
    ax.set_ylabel(r"$10^{3}\,c_{f,s}$, upper surface")
    ax.set_title("(b)", pad=3.5)

    # (c) fraction of band faces in reverse flow
    ax = axes[2]
    for g in geoms_present:
        st = cond["geoms"][g]["stations"].get(eta)
        if st is None:
            continue
        ax.plot(CTR, 100.0 * st["revfrac_up"], zorder=4, **STYLE[g])
    ax.set_xlim(*xlim)
    ax.set_ylim(-4, 104)
    ax.set_ylabel("reversed upper-surface faces, %")
    ax.set_title("(c)", pad=3.5)

    for ax in axes:
        ax.set_xlabel(r"$x/c$")
        ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        ax.grid(True, lw=0.3, color="0.85", zorder=0)
        ax.set_axisbelow(True)

    turb = cond["turb"][0] if cond["turb"] else "NOT FOUND"
    yps = [cond["geoms"][g]["stations"][eta]["yplus_median_upper"]
           for g in geoms_present if eta in cond["geoms"][g]["stations"]]
    fig.suptitle("", fontsize=7.6, y=0.988, va="top", linespacing=1.45)  # caption carries it

    handles = [Line2D([], [], label=GEOM_LABEL[g],
                      **dict(STYLE[g], lw=max(STYLE[g]["lw"], 1.45)))
               for g in geoms_present]
    if cps is not None:
        handles.append(Line2D([], [], color="0.35", lw=0.9, ls=(0, (1.0, 1.6)),
                              label=r"$c_p^{*}$ at M %.4f" % cond_mach(cond)))
    handles.append(Patch(facecolor=STYLE["C"]["color"], alpha=0.13, label="C reverse-flow extent"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 0.004), handlelength=3.1, columnspacing=1.4)
    return fig


# ----------------------------------------------------------------------------- main
def main():
    data, present, absent = build_panel_data()

    measured = dict(shock={}, separation={}, maxabs_dcp={})
    for ck in ("CR", "EC", "LC"):
        cond = data[ck]
        for g in GEOMS:
            for eta in STATIONS:
                st = cond["geoms"][g]["stations"].get(eta)
                if st is None:
                    continue
                if st["mach_from_file"]:
                    measured["shock"][(ck, g, eta)] = detect_shock(st)
                measured["separation"][(ck, g, eta)] = detect_separation(st)
                if g != "B":
                    entry = {}
                    for br, brname in (("up", "upper"), ("lo", "lower")):
                        d = diff_vs_baseline(cond, g, eta, br)
                        if d is None:
                            continue
                        w = (CTR >= DIFF_XMIN) & (CTR < 0.98)
                        if np.any(~np.isnan(d[w])):
                            k = int(np.nanargmax(np.abs(np.where(w, d, np.nan))))
                            entry[brname] = dict(
                                max_abs_dcp=float(abs(d[k])),
                                at_x_over_c=float(CTR[k]),
                                signed=float(d[k]),
                                x_window=[DIFF_XMIN, 0.98],
                                le_roughness=le_roughness(d),
                            )
                    measured["maxabs_dcp"][(ck, g, eta)] = entry

    # --- within-file cross-check of the derived normalisation, asserted not assumed
    checks = {}
    for ck, expect_q in (("CR", 578.0), ("EC", 11292.8)):
        qs = sorted({data[ck]["geoms"][g]["stations"][e]["q"]
                     for g in GEOMS for e in STATIONS
                     if e in data[ck]["geoms"][g]["stations"]})
        assert len(qs) == 1 and abs(qs[0] - expect_q) < 1e-6, (ck, qs)
        checks["q_%s" % ck] = qs[0]
    m_ec = cond_mach(data["EC"])
    checks["mach_EC_from_file"] = m_ec
    checks["cp_star_EC"] = cp_star(m_ec)

    # --- wall treatment derived from the solution's own y+ field, not from a recipe name
    for ck in ("CR", "EC"):
        yp = [data[ck]["geoms"][g]["stations"][e]["yplus_median_upper"]
              for g in GEOMS for e in STATIONS if e in data[ck]["geoms"][g]["stations"]]
        checks["yplus_median_upper_range_%s" % ck] = [float(min(yp)), float(max(yp))]
    assert checks["yplus_median_upper_range_CR"][1] < 1.0, "Condition CR is not wall-resolved"
    assert checks["yplus_median_upper_range_EC"][0] > 5.0, "early cruise is not wall-modelled"
    checks["cf_cross_condition_comparable"] = False

    # --- the LE cut on the difference panels must be earned by a measurement
    rough_fwd, rough_aft = [], []
    for e in measured["maxabs_dcp"].values():
        for br in e.values():
            r = br["le_roughness"]
            if r["forward_of_cut"] is not None:
                rough_fwd.append(r["forward_of_cut"])
            if r["aft_of_cut"] is not None:
                rough_aft.append(r["aft_of_cut"])
    checks["diff_panel_x_min"] = DIFF_XMIN
    checks["le_roughness_median_forward_of_cut"] = float(np.median(rough_fwd))
    checks["le_roughness_median_aft_of_cut"] = float(np.median(rough_aft))
    checks["le_roughness_ratio"] = (
        checks["le_roughness_median_forward_of_cut"]
        / checks["le_roughness_median_aft_of_cut"])
    assert checks["le_roughness_ratio"] > 3.0, (
        "the LE cut is not justified by the data: roughness ratio %.2f"
        % checks["le_roughness_ratio"])

    # --- the shock threshold must sit in a gap, or it is an arbitrary number
    grads_flag = sorted(v["max_dcp_dx"] for v in measured["shock"].values() if v.get("flagged"))
    grads_no = sorted(v["max_dcp_dx"] for v in measured["shock"].values()
                      if not v.get("flagged") and "max_dcp_dx" in v)
    checks["shock_threshold"] = SHOCK_GRAD_THRESHOLD
    checks["flagged_gradients"] = grads_flag
    checks["largest_unflagged_gradient"] = grads_no[-1] if grads_no else None
    if grads_flag and grads_no:
        assert grads_no[-1] < SHOCK_GRAD_THRESHOLD < grads_flag[0], (
            "shock threshold does not sit in a gap: unflagged max %.2f, flagged min %.2f"
            % (grads_no[-1], grads_flag[0]))

    # --- build, then write
    written = []
    jobs = []
    for ck, fname in (("CR", "fig_cp_sections_CR.pdf"),
                      ("EC", "fig_cp_sections_earlycruise.pdf"),
                      ("LC", "fig_cp_sections_latecruise.pdf")):
        cond = data[ck]
        gp = [g for g in GEOMS if cond["geoms"][g]["stations"]]
        if not gp:
            print("SKIP %s: no station data for any geometry" % ck)
            continue
        out = os.path.join(HERE, fname)
        fig = make_condition_figure(cond, out, gp, measured)
        jobs.append((fig, out, gp, ck))

    cond = data["EC"]
    gp_ec = [g for g in GEOMS if cond["geoms"][g]["stations"]]
    if gp_ec and 80 in cond["geoms"]["B"]["stations"]:
        out = os.path.join(HERE, "fig_cp_shock_separation_CMPC.pdf")
        fig = make_shock_detail_figure(cond, out, gp_ec, measured)
        jobs.append((fig, out, gp_ec, "EC80"))

    for fig, out, gp, tag in jobs:
        fig.savefig(out, format="pdf", bbox_inches="tight")
        plt.close(fig)
        written.append((out, gp, tag))

    # --- machine-readable record of every number the figures display
    def keystr(k):
        return "%s|%s|eta%03d" % k

    summary = dict(
        _generated="make_cp.py",
        _frame=dict(
            note=("cp is on each case's own freestream q, read from that file's header. It "
                  "carries no reference area and is on neither the DSO nor the TP-1580 basis. "
                  "x/c is on the band's own local chord. Force coefficients quoted alongside "
                  "are DSO basis, half-model, Aref 0.620462 m2, from rans_forces.json."),
            gamma=GAMMA,
        ),
        _method=dict(
            binning="%d uniform x/c bins, median per bin, bins with < %d faces left empty"
                    % (NBIN, MINPTS),
            shock_criterion=("max d(cp)/d(x/c) on the binned upper branch inside x/c "
                             "[%.2f, %.2f], flagged when >= %.1f per unit chord AND cp at "
                             "x_shock - %.3f c is below cp*(M_inf)"
                             % (SHOCK_WINDOW[0], SHOCK_WINDOW[1], SHOCK_GRAD_THRESHOLD,
                                SHOCK_HALFSPAN)),
            separation_criterion=("binned median of cf_s on the upper branch below zero over "
                                  "at least 3 consecutive bins (0.015 c)"),
            difference_panels=("cp(candidate) - cp(baseline B) per branch on the shared bin "
                               "grid, plotted and scaled over x/c >= %.2f only" % DIFF_XMIN),
        ),
        _checks=checks,
        _partition=dict(
            expected_case_station=len(present) + len(absent),
            present=len(present),
            absent=len(absent),
            absent_list=sorted({c for c, _ in absent}),
            sums_to_total=True,
        ),
        conditions={},
        shock={keystr(k): v for k, v in measured["shock"].items()},
        separation={keystr(k): v for k, v in measured["separation"].items()},
        max_abs_delta_cp={keystr(k): v for k, v in measured["maxabs_dcp"].items()},
        # WRITTEN OUT, not merely computed. The first version of this built the
        # suction_peak block into each station dict and stopped there, and the station
        # dicts are never serialised, so the numbers existed in memory and in nothing
        # a reader could open. Computing a fact and reporting it as delivered is the
        # vouched-success shape (GEO-080); caught by querying the written file instead
        # of trusting the edit.
        suction_peak={
            keystr((ck, g, eta)): st["suction_peak"]
            for ck in ("CR", "EC", "LC") if ck in data
            for g, gd in data[ck].get("geoms", {}).items()
            for eta, st in gd.get("stations", {}).items()
            if "suction_peak" in st
        },
        figures=[dict(path=p, geometries=g, tag=t) for p, g, t in written],
    )
    for ck in ("CR", "EC", "LC"):
        c = data[ck]
        summary["conditions"][ck] = dict(
            title=c["title"], target_CL=c["target_CL"], mach_nominal=c["mach_nominal"],
            Aref_m2=c["Aref_m2"], lRef_m=c["lRef_m"], basis=c["basis"], turbulence=c["turb"],
            geometries={g: dict(case=c["geoms"][g]["case"],
                                Cl=c["geoms"][g]["Cl"],
                                alpha_deg=c["geoms"][g]["alpha_deg"],
                                alpha_deg_from_section_dragDir=c["geoms"][g]["alpha_from_section"],
                                Cd_counts=c["geoms"][g]["Cd_counts"],
                                delta_Cd_counts_vs_B=c["geoms"][g]["dCd_counts"],
                                trim_within_tolerance=c["geoms"][g]["trim_ok"],
                                surface_time=c["geoms"][g]["surface_time"],
                                stations_present=sorted(c["geoms"][g]["stations"].keys()),
                                yplus_median_upper={
                                    "eta%03d" % e: c["geoms"][g]["stations"][e]["yplus_median_upper"]
                                    for e in sorted(c["geoms"][g]["stations"])},
                                local_chord_m={
                                    "eta%03d" % e: c["geoms"][g]["stations"][e]["chord_m"]
                                    for e in sorted(c["geoms"][g]["stations"])})
                        for g in GEOMS},
        )
    jpath = os.path.join(HERE, "cp_sections_measured.json")
    tmp = jpath + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True, default=str)
    json.load(open(tmp))  # validate before promoting to the final name
    os.replace(tmp, jpath)
    written.append((jpath, [], "json"))

    # --- derive the written claim by stat-ing the files, never from the exit status
    print("\nFILES (derived by stat, not asserted):")
    bad = 0
    for p, _, _ in written:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            print("  OK   %10d bytes  %s" % (os.path.getsize(p), p))
        else:
            print("  FAIL missing or empty      %s" % p)
            bad += 1

    print("\nCHECKS: q_CR=%s q_EC=%s M_EC=%.6f cp*=%.4f" % (
        checks["q_CR"], checks["q_EC"], checks["mach_EC_from_file"], checks["cp_star_EC"]))
    print("SHOCK threshold %.1f; flagged %s; largest unflagged %.2f" % (
        SHOCK_GRAD_THRESHOLD, ["%.1f" % g for g in grads_flag],
        checks["largest_unflagged_gradient"]))
    print("\nFLAGGED SHOCKS")
    for k, v in sorted(measured["shock"].items()):
        if v.get("flagged"):
            print("  %-3s %s eta%03d  x/c=%.3f  dcp/dx=%6.1f  dcp(+-%.2fc)=%+.3f  cp_us=%+.3f"
                  % (k[0], k[1], k[2], v["x_over_c"], v["max_dcp_dx"], SHOCK_HALFSPAN,
                     v["delta_cp"], v["cp_upstream"]))
    print("\nUPPER-SURFACE SEPARATION (binned median cf_s < 0 over >= 3 bins)")
    for k, v in sorted(measured["separation"].items()):
        if v["any"]:
            for r in v["runs"]:
                print("  %-3s %s eta%03d  x/c %.3f to %.3f  min cf_s %+.2e  to TE: %s  "
                      "(band reverse fraction all faces %.4f)"
                      % (k[0], k[1], k[2], r["x_start"], r["x_end"], r["min_cf_s"],
                         r["reaches_trailing_edge"], v["band_reverse_fraction_all_faces"]))
    print("\nMAX |delta cp| vs B")
    for k in sorted(measured["maxabs_dcp"]):
        e = measured["maxabs_dcp"][k]
        s = "  ".join("%s %.3f@%.2f" % (b, e[b]["max_abs_dcp"], e[b]["at_x_over_c"]) for b in e)
        print("  %-3s %s eta%03d  %s" % (k[0], k[1], k[2], s))
    print("\nABSENT (case, station) pairs: %d, cases %s"
          % (len(absent), sorted({c for c, _ in absent})))

    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
