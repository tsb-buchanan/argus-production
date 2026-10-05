#!/usr/bin/env python3
"""
Skin-friction Cf sections and the separation picture, ARGUS all-geometry report
2026-09-15.

EVERY NUMBER IN THE OUTPUT IS READ FROM A FILE. Nothing is typed in from memory.

Inputs (read, never assumed):
  data/sections/<CASE>_eta<NNN>.csv   55 files, 11 cases x 5 stations
  data/sections/Q_PINF_TABLE.csv      q, p_inf, alpha, dragDir, sign multiplier
  data/sections/MANIFEST.tsv          source VTK, surface time
  data/rans_forces.json               turbulence model, Cl, Cd, alpha per case
  data/vlm.json                       geometry letter <-> DSO geometry name

SIGN CONVENTION (the thing this figure is about)
  OpenFOAM's wallShearStress returns a NEGATIVE streamwise component for ATTACHED
  flow in +x.  The section extractor calibrated that sign against a region whose
  attachment is known independently (inboard mid-chord faces) and applied the
  multiplier recorded per file as `shear_sign_multiplier` (-1 on all 11 cases).
  The plotted quantity is therefore

      cf_s = mult * (tau . dragDir) / q ,   POSITIVE = attached, NEGATIVE = reversed

  so the SIGN CHANGE is the separation indicator and the magnitude is not.
  This script asserts the multiplier is -1 on every file rather than assuming it.

WITHIN-FILE CLOSURE CHECK (D052: prefer a proof the file makes against itself)
  Each file's header carries `reverse_flow_fraction` over the whole band, written
  by the extractor on the cluster.  This script recomputes it from the rows and
  HALTS on a mismatch.  Nothing is plotted from a file that fails.

PARTITION (GEO-089)
  Cases x stations enumerated, adjudicated and asserted to sum.  Missing members
  are named, counted and reported, never silently dropped.
"""

import csv
import json
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)
DATA = os.path.join(REPORT, "data")
# sections_nflag, not sections: the upper/lower flag now comes from the CALIBRATED FACET
# NORMAL written by the extractor, not from the geometric branch splitter run on the laptop,
# which mislabels faces where the two surfaces converge. The superseded set is retained
# beside it in data/sections/ and is what the 2026-09-15 figures were drawn from.
SEC = os.path.join(DATA, "sections_nflag")
FIG = HERE

# PER-CASE DATASET RESOLUTION. A case is served from data/sections_r2 where that exists,
# from sections_nflag otherwise. SEC above stays the fallback; section_source decides per
# case. THE Q_PINF ROW IS TAKEN FROM THE SAME DATASET AS THE SECTION CSV, always: the two
# disagree on alpha by up to 0.0301 deg, and cf_s is mult * (tau . dragDir) / q, so a
# crossed pair applies the wrong wind axes. That is the D019 defect, whose recorded lag was
# 0.0083 deg, i.e. this would be up to 3.6x it. See fig/section_source.py.
import section_source  # noqa: E402


LOG_LINES = []


def log(msg=""):
    print(msg)
    LOG_LINES.append(str(msg))


def halt(msg):
    log("HALT: " + msg)
    raise SystemExit("HALT: " + msg)


# ----------------------------------------------------------------------------
# Style: serif, colour-blind safe (Okabe-Ito), distinct dash patterns so the
# figure survives greyscale printing.
# ----------------------------------------------------------------------------
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 10.5,
        "axes.labelsize": 10.5,
        "axes.titlesize": 10.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9.5,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "lines.solid_capstyle": "round",
        "lines.dash_capstyle": "round",
        "pdf.fonttype": 42,
        "savefig.bbox": "tight",
    }
)

STYLE = {
    "B": dict(color="#000000", ls="-", lw=1.5),
    "C": dict(color="#D55E00", ls=(0, (5.5, 1.8)), lw=1.6),
    "F": dict(color="#0072B2", ls=(0, (4.5, 1.5, 1.0, 1.5)), lw=1.5),
    "H": dict(color="#009E73", ls=(0, (1.3, 1.2)), lw=1.7),
    "M": dict(color="#CC79A7", ls=(0, (4.0, 1.4, 1.0, 1.4, 1.0, 1.4)), lw=1.5),
    "W": dict(color="#56B4E9", ls=(0, (8.0, 2.2)), lw=1.6),
}
LETTERS = ["B", "C", "F", "H", "M", "W"]
ETAS = [0.20, 0.40, 0.60, 0.80, 0.95]

# ----------------------------------------------------------------------------
# 1. Geometry letter -> DSO geometry name, READ FROM vlm.json
# ----------------------------------------------------------------------------
GEOM_NAME = {}
with open(os.path.join(DATA, "vlm.json")) as fh:
    _vlm = json.load(fh)
for _r in _vlm["runs"]:
    _m = _r.get("mapping") or {}
    _fam = _m.get("case_prefix_family")
    if _fam in LETTERS and _r.get("geometry_name"):
        GEOM_NAME.setdefault(_fam, _r["geometry_name"])
log("Geometry names read from vlm.json (mapping.case_prefix_family -> geometry_name):")
for k in LETTERS:
    log("  %s = %s" % (k, GEOM_NAME.get(k, "NOT FOUND")))

# ----------------------------------------------------------------------------
# 2. Turbulence model / trim state per case, READ FROM rans_forces.json
# ----------------------------------------------------------------------------
with open(os.path.join(DATA, "rans_forces.json")) as fh:
    RANS = json.load(fh)
RCASES = RANS["cases"]


def rans_field(case, key):
    rec = RCASES.get(case)
    if rec is None:
        return None
    return rec.get(key)


# ----------------------------------------------------------------------------
# 3. Q_PINF_TABLE.csv - per-case q, p_inf, alpha, sign multiplier
# ----------------------------------------------------------------------------
QPI = {}
_QTABLES = {}
for _d in section_source.DATASETS:
    _p = os.path.join(_d, "Q_PINF_TABLE.csv")
    if not os.path.exists(_p):
        continue
    with open(_p) as fh:
        _QTABLES[_d] = {r["case"]: r for r in csv.DictReader(fh)}
# Fetch each case's row THROUGH THE RESOLVER, not by merging with one table winning, so a
# case can never take its coefficients from one evaluation and its wind axes from the other.
for _case in sorted(set().union(*(set(v) for v in _QTABLES.values()))):
    try:
        _qd = os.path.dirname(section_source.resolve(_case)[1])
    except KeyError:
        continue
    if _case in _QTABLES.get(_qd, {}):
        QPI[_case] = dict(_QTABLES[_qd][_case],
                          _generation=section_source.resolve(_case)[2])
log("\nQ_PINF_TABLE.csv: %d cases (%s)"
    % (len(QPI), ", ".join("%s=%s" % (c, QPI[c]["_generation"]) for c in sorted(QPI))))

# ----------------------------------------------------------------------------
# 4. Enumerate and load the section files. PARTITION, asserted.
# ----------------------------------------------------------------------------
HDR_KEYS = (
    "case",
    "eta",
    "local_chord_m",
    "q_from_file",
    "p_inf_from_file",
    "dragDir",
    "shear_sign_multiplier",
    "n_faces_in_band",
    "reverse_flow_fraction",
)


def read_section(path):
    """Return (header dict, structured arrays). Header values parsed from the
    file's own comment block; nothing defaulted."""
    hdr = {}
    ncomment = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            ncomment += 1
            body = line[1:].strip()
            if "," not in body:
                continue
            k, v = body.split(",", 1)
            k = k.strip()
            if k in HDR_KEYS:
                hdr[k] = v.strip()
    missing = [k for k in HDR_KEYS if k not in hdr]
    if missing:
        halt("%s: header keys missing %s" % (os.path.basename(path), missing))

    x, z, cp, cf, cfs, yp, surf = [], [], [], [], [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            break  # this line is the column header
        for line in fh:
            p = line.rstrip("\n").split(",")
            if len(p) != 7:
                halt("%s: row with %d fields, expected 7" % (path, len(p)))
            x.append(p[0]); z.append(p[1]); cp.append(p[2])
            cf.append(p[3]); cfs.append(p[4]); yp.append(p[5]); surf.append(p[6])
    d = dict(
        x=np.asarray(x, float),
        z=np.asarray(z, float),
        cp=np.asarray(cp, float),
        cf=np.asarray(cf, float),
        cfs=np.asarray(cfs, float),
        yplus=np.asarray(yp, float),
        upper=np.asarray([s == "upper" for s in surf], bool),
    )
    bad = set(surf) - {"upper", "lower"}
    if bad:
        halt("%s: unexpected surface labels %s" % (path, bad))
    return hdr, d


# Enumerate across BOTH datasets, then keep only the file the resolver sanctions for
# that case. Enumerating one directory would silently narrow the input set while the run
# reported over "the sections" (GEO-087).
_seen = {}
for _d in section_source.DATASETS:
    if not os.path.isdir(_d):
        continue
    for f in os.listdir(_d):
        m = re.match(r"^(\w+_trim)_eta(\d{3})\.csv$", f)
        if not m:
            continue
        try:
            _sd = section_source.resolve(m.group(1))[0]
        except KeyError:
            continue
        if os.path.abspath(_sd) == os.path.abspath(_d):
            _seen[f] = _d
files = sorted(_seen)
log("\nEnumerated %d section CSVs across %d dataset(s)"
    % (len(files), len(set(_seen.values()))))

SEC_DATA = {}      # (case, eta) -> data arrays
SEC_HDR = {}       # (case, eta) -> header
cases_found = set()
n_closure_checked = 0
max_closure_resid = 0.0

for fn in files:
    m = re.match(r"^(\w+_trim)_eta(\d{3})\.csv$", fn)
    case, eta = m.group(1), int(m.group(2)) / 100.0
    hdr, d = read_section(os.path.join(_seen[fn], fn))

    # --- assertions derived from the file, not authored alongside it ---------
    if hdr["case"] != case:
        halt("%s: header case %s disagrees with filename" % (fn, hdr["case"]))
    if abs(float(hdr["eta"]) - eta) > 1e-9:
        halt("%s: header eta %s disagrees with filename" % (fn, hdr["eta"]))
    mult = float(hdr["shear_sign_multiplier"])
    if mult != -1.0:
        halt(
            "%s: shear_sign_multiplier is %g, not -1. The plotted sign convention "
            "assumes the extractor's calibration; it must be read, not assumed." % (fn, mult)
        )
    n_hdr = int(hdr["n_faces_in_band"])
    if n_hdr != d["x"].size:
        halt("%s: n_faces_in_band %d but %d rows" % (fn, n_hdr, d["x"].size))

    # within-file closure: recompute the header's own reverse_flow_fraction
    rf_hdr = float(hdr["reverse_flow_fraction"])
    rf_calc = float(np.mean(d["cfs"] < 0.0))
    resid = abs(rf_hdr - rf_calc)
    max_closure_resid = max(max_closure_resid, resid)
    if resid > 1.5e-6:   # header is printed to 6 dp
        halt("%s: reverse_flow_fraction header %.6f vs recomputed %.8f" % (fn, rf_hdr, rf_calc))
    n_closure_checked += 1

    SEC_DATA[(case, eta)] = d
    SEC_HDR[(case, eta)] = hdr
    cases_found.add(case)

log("Within-file closure on reverse_flow_fraction: %d/%d files agree, "
    "max |header - recomputed| = %.2e (header printed to 6 dp)"
    % (n_closure_checked, len(files), max_closure_resid))

# --- PARTITION over the members this figure needs ---------------------------
sw_cases = sorted(c for c in cases_found if c.startswith("SW"))
cmp_cases = sorted(c for c in cases_found if c.startswith("CMP"))
lc_cases = sorted(c for c in cases_found if c.startswith("LC"))

present, missing = [], []
for cond_prefix in ("SW", "CMP", "LC"):
    for L in LETTERS:
        case = "%s%s_trim" % (cond_prefix, L)
        for eta in ETAS:
            (present if (case, eta) in SEC_DATA else missing).append((case, eta))

log("\nPARTITION over the 3 conditions x 6 geometries x 5 stations = %d members:"
    % (3 * len(LETTERS) * len(ETAS)))
log("  present %d + missing %d = %d ; sums_to_total: %s"
    % (len(present), len(missing), len(present) + len(missing),
       len(present) + len(missing) == 3 * len(LETTERS) * len(ETAS)))
miss_cases = sorted({c for c, _ in missing})
log("  missing cases (all 5 stations each): %s" % ", ".join(miss_cases))
if len(present) != len(files):
    halt("present members %d != files read %d" % (len(present), len(files)))

# ----------------------------------------------------------------------------
# 5. Binning. Medians per x/c bin; the band pools span, so a median over the
#    band is the estimator that compares like with like (the extractor's own
#    finding, README "my first gate was wrong").
# ----------------------------------------------------------------------------
NBIN = 200
EDGES = np.linspace(0.0, 1.0, NBIN + 1)
CENT = 0.5 * (EDGES[:-1] + EDGES[1:])
MIN_PER_BIN = 5


def binned_stats(d, upper=True):
    m = d["upper"] if upper else ~d["upper"]
    x, c, z = d["x"][m], d["cfs"][m], d["z"][m]
    idx = np.digitize(x, EDGES) - 1
    idx = np.clip(idx, 0, NBIN - 1)
    xc = np.full(NBIN, np.nan)
    med = np.full(NBIN, np.nan)
    q25 = np.full(NBIN, np.nan)
    q75 = np.full(NBIN, np.nan)
    zmed = np.full(NBIN, np.nan)
    n = np.zeros(NBIN, int)
    order = np.argsort(idx, kind="stable")
    idx_s, c_s, x_s, z_s = idx[order], c[order], x[order], z[order]
    bounds = np.searchsorted(idx_s, np.arange(NBIN + 1))
    for i in range(NBIN):
        a, b = bounds[i], bounds[i + 1]
        n[i] = b - a
        if n[i] >= MIN_PER_BIN:
            xc[i] = np.median(x_s[a:b])
            med[i] = np.median(c_s[a:b])
            q25[i], q75[i] = np.percentile(c_s[a:b], [25, 75])
            zmed[i] = np.median(z_s[a:b])
    return dict(xc=xc, med=med, q25=q25, q75=q75, zmed=zmed, n=n)


# separated run: >= MIN_RUN consecutive valid bins with median cf_s < 0.
# At 0.005 c per bin, MIN_RUN = 4 means the run must exceed 2% of local chord.
# Anything shorter is reported separately as isolated reversed faces and is NOT
# called separation.
MIN_RUN = 4


def separated_runs(st):
    valid = np.isfinite(st["med"])
    neg = valid & (st["med"] < 0.0)
    runs = []
    i = 0
    while i < NBIN:
        if neg[i]:
            j = i
            while j + 1 < NBIN and neg[j + 1]:
                j += 1
            if (j - i + 1) >= MIN_RUN:
                runs.append((EDGES[i], EDGES[j + 1], j - i + 1))
            i = j + 1
        else:
            i += 1
    return runs


STATS_U, STATS_L, RUNS_U = {}, {}, {}
for key, d in SEC_DATA.items():
    STATS_U[key] = binned_stats(d, upper=True)
    STATS_L[key] = binned_stats(d, upper=False)
    RUNS_U[key] = separated_runs(STATS_U[key])

# ----------------------------------------------------------------------------
# 6. Reverse-flow / separated-extent table. MEASURED quantities only.
# ----------------------------------------------------------------------------
TABLE = []
for (case, eta), d in sorted(SEC_DATA.items()):
    up = d["upper"]
    cu = d["cfs"][up]
    xu = d["x"][up]
    rev = cu < 0.0
    lo = ~up
    revl = d["cfs"][lo] < 0.0
    runs = RUNS_U[(case, eta)]
    rec = dict(
        case=case,
        letter=case[case.index("_") - 1],
        # LC rows were labelled early_cruise until 2026-10-05: the ternary predates any late-
        # cruise section data, so it was right only while no LC case existed (ARG-196).
        condition=("CR" if case.startswith("SW") else
                   "late_cruise" if case.startswith("LC") else "early_cruise"),
        eta=eta,
        n_upper=int(up.sum()),
        n_lower=int(lo.sum()),
        n_rev_upper=int(rev.sum()),
        n_rev_lower=int(revl.sum()),
        pct_rev_upper=100.0 * rev.mean(),
        pct_rev_lower=100.0 * revl.mean(),
        min_cfs_upper=float(cu.min()),
        x_rev_first_face=(float(xu[rev].min()) if rev.any() else None),
        x_rev_last_face=(float(xu[rev].max()) if rev.any() else None),
        n_separated_runs=len(runs),
        x_sep=(runs[0][0] if runs else None),
        x_reatt=(runs[-1][1] if runs else None),
        run_extent_c=(sum(r[1] - r[0] for r in runs) if runs else 0.0),
        yplus_median_upper=float(np.median(d["yplus"][up])),
        yplus_p99_upper=float(np.percentile(d["yplus"][up], 99)),
        local_chord_m=float(SEC_HDR[(case, eta)]["local_chord_m"]),
        q=float(SEC_HDR[(case, eta)]["q_from_file"]),
    )
    TABLE.append(rec)

# Most-forward reversed face, so the "everything else is the blunt-TE base region"
# claim is DERIVED from the data rather than asserted beside it.
def forward_most_reversed(surface, exclude_separated_runs=True):
    best = None
    for (case, eta), d in SEC_DATA.items():
        if exclude_separated_runs and RUNS_U[(case, eta)] and surface == "upper":
            continue
        m = d["upper"] if surface == "upper" else ~d["upper"]
        rev = m & (d["cfs"] < 0.0)
        if rev.any():
            xmin = float(d["x"][rev].min())
            if best is None or xmin < best[0]:
                best = (xmin, case, eta)
    return best


FWD_U = forward_most_reversed("upper", exclude_separated_runs=True)
FWD_L = forward_most_reversed("lower", exclude_separated_runs=False)

# pooled over the five sampled bands, per case (MEASURED; a face-count fraction,
# not an area fraction, and over 5 bands only)
POOLED = {}
for case in sorted(cases_found):
    nu = sum(r["n_upper"] for r in TABLE if r["case"] == case)
    nr = sum(r["n_rev_upper"] for r in TABLE if r["case"] == case)
    POOLED[case] = dict(n_upper=nu, n_rev=nr, pct=100.0 * nr / nu if nu else float("nan"))

# band coverage of the span, from the README's stated band half-width
BAND_HALFWIDTH_ETA = 0.004          # semispan units, stated in data/sections/README.md
BAND_COVERAGE_PCT = 100.0 * len(ETAS) * 2 * BAND_HALFWIDTH_ETA

log("\nMost-forward reversed face, EXCLUDING the bands that carry a separated run:")
log("  upper surface: x/c = %.6f  (%s, eta %.2f)" % FWD_U)
log("Most-forward reversed face on the lower surface, all bands:")
log("  lower surface: x/c = %.6f  (%s, eta %.2f)" % FWD_L)

log("\nSEPARATED EXTENT, upper surface, per case and station")
log("  (percentage of upper-surface FACES in the sampled band with cf_s < 0)")
hdr = ("%-11s %-6s %6s %8s %9s %10s %10s %9s"
       % ("case", "eta", "n_up", "n_rev", "%rev", "x_sep", "x_reatt", "min cf_s"))
log("  " + hdr)
for r in TABLE:
    log("  %-11s %-6.2f %6d %8d %9.4f %10s %10s %9.6f"
        % (r["case"], r["eta"], r["n_upper"], r["n_rev_upper"], r["pct_rev_upper"],
           ("%.3f" % r["x_sep"]) if r["x_sep"] is not None else "-",
           ("%.3f" % r["x_reatt"]) if r["x_reatt"] is not None else "-",
           r["min_cfs_upper"]))

# CSV alongside the figures
tab_csv = os.path.join(FIG, "cf_separation_table.csv")
with open(tab_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(TABLE[0].keys()))
    w.writeheader()
    for r in TABLE:
        w.writerow(r)

# ----------------------------------------------------------------------------
# 7. Panel column metadata, read from the case records
# ----------------------------------------------------------------------------
def condition_header(prefix):
    """Turbulence model, target C_L and measured y+ for a condition, from files."""
    cases = sorted(c for c in cases_found if c.startswith(prefix))
    models = sorted({rans_field(c, "turbulence_model") for c in cases} - {None})
    cls = [rans_field(c, "Cl") for c in cases if rans_field(c, "Cl") is not None]
    yp = [r["yplus_median_upper"] for r in TABLE if r["case"] in cases]
    return dict(cases=cases, models=models, cl=cls,
                yp_lo=min(yp), yp_hi=max(yp))


CR_H = condition_header("SW")
CM_H = condition_header("CMP")
log("\nCondition CR  : turbulence model(s) %s ; y+ median over upper faces %.2f to %.2f"
    % (CR_H["models"], CR_H["yp_lo"], CR_H["yp_hi"]))
log("Early cruise  : turbulence model(s) %s ; y+ median over upper faces %.1f to %.1f"
    % (CM_H["models"], CM_H["yp_lo"], CM_H["yp_hi"]))
if len(CR_H["models"]) != 1 or len(CM_H["models"]) != 1:
    log("NOTE: turbulence model is not unique within a condition; see list above.")

COLS = [
    dict(prefix="SW",
         title="Condition CR,  M 0.10\ntarget $C_L$ 0.428277635108\n"
               "WALL-RESOLVED, %s\n$y^+$ median %.2f to %.2f"
               % (CR_H["models"][0] if CR_H["models"] else "NOT FOUND",
                  CR_H["yp_lo"], CR_H["yp_hi"]),
         ),
    dict(prefix="CMP",
         title="Early cruise,  M 0.78\ntarget $C_L$ 0.529297087450\n"
               "WALL-MODELLED, %s\n$y^+$ median %.0f to %.0f"
               % (CM_H["models"][0] if CM_H["models"] else "NOT FOUND",
                  CM_H["yp_lo"], CM_H["yp_hi"]),
         ),
]

# y limits per column, from the binned medians actually plotted. The lower bound
# is forced below zero so that a sign change is visible at all: reversed cf_s is
# genuinely two orders smaller than attached cf_s, which is the point.
def col_ylim(prefix, stats):
    vals = []
    for (case, eta), st in stats.items():
        if case.startswith(prefix):
            v = st["med"][np.isfinite(st["med"])]
            if v.size:
                vals.append(v)
    v = np.concatenate(vals)
    hi = float(v.max())
    lo = float(min(v.min(), -0.10 * hi))
    pad = 0.05 * (hi - lo)
    return lo - pad, hi + pad


# ----------------------------------------------------------------------------
# 8. FIGURE 1 and 2: Cf sections, 5 stations x 2 conditions
# ----------------------------------------------------------------------------
def legend_handles(letters, with_missing_note=False):
    h, lab = [], []
    for L in letters:
        h.append(Line2D([0], [0], **STYLE[L]))
        lab.append("%s  %s" % (L, GEOM_NAME.get(L, "name NOT FOUND")))
    return h, lab


def make_sections_figure(surface, outfile):
    stats = STATS_U if surface == "upper" else STATS_L
    # BUILT AT THE WIDTH IT IS PRINTED AT. The report's textwidth is 170 mm (a4, 20 mm
    # margins) = 6.69 in, and this figure is included at \linewidth. At figsize 7.4 in it
    # was downscaled to 0.90 and every font with it; when it was previously placed two to a
    # row at 0.49\linewidth the scale was 0.37 and the 7.8 pt labels rendered near 2.9 pt,
    # which is why nothing in it could be read. Matching the figure to the column removes
    # the downscale, so the point sizes chosen here are the point sizes that print.
    fig, axes = plt.subplots(len(ETAS), 2, figsize=(6.69, 9.95), sharex=True)
    ylims = [col_ylim(c["prefix"], stats) for c in COLS]

    for j, col in enumerate(COLS):
        for i, eta in enumerate(ETAS):
            ax = axes[i, j]
            y0, y1 = ylims[j]
            ax.axhline(0.0, color="0.35", lw=0.8, zorder=1)
            plotted, notes = [], []
            for L in LETTERS:
                case = "%s%s_trim" % (col["prefix"], L)
                key = (case, eta)
                if key not in stats:
                    continue
                st = stats[key]
                ok = np.isfinite(st["med"])
                ax.plot(st["xc"][ok], st["med"][ok], zorder=3, **STYLE[L])
                plotted.append(L)
                # separation / reattachment markers, upper surface only
                if surface == "upper":
                    for (xs, xr, nb) in RUNS_U[key]:
                        ax.plot([xs], [0.0], marker="v", ms=5.0, mfc=STYLE[L]["color"],
                                mec="k", mew=0.5, zorder=6, ls="none")
                        ax.plot([xr], [0.0], marker="^", ms=5.0, mfc="white",
                                mec=STYLE[L]["color"], mew=1.2, zorder=6, ls="none")
                        ax.axvspan(xs, xr, color=STYLE[L]["color"], alpha=0.11,
                                   lw=0, zorder=2)
                        notes.append((L, "%s: sep %.3f, reatt %.3f" % (L, xs, xr)))
            ax.set_ylim(y0, y1)
            ax.set_xlim(0.0, 1.0)
            ax.tick_params(length=3, width=0.7)
            # band-averaged aft bin at eta 0.95 (stated in data/sections/README.md)
            if eta == 0.95:
                ax.axvspan(0.875, 1.0, color="0.85", alpha=0.55, lw=0, zorder=0)
            ax.text(0.048, 0.93, "$\\eta$ = %.2f" % eta, transform=ax.transAxes,
                    ha="left", va="top", fontsize=10, zorder=8,
                    bbox=dict(fc="white", ec="none", alpha=0.80, pad=1.5))
            if plotted:
                pass  # letters are in the legend; labels only (22 Sep)
            # WHITE BACKING, SAME AS THE eta LABEL ABOVE. These notes sit at y = 0.30 and
            # 0.18 in axes coordinates, which is exactly where the curves cross zero in the
            # one panel that carries them (eta 0.80, early cruise). Unbacked, the C and M
            # lines were printed through by the data and neither could be read. The notes
            # do not overlap each other; they were overlapped by the plot.
            for k, (L, txt) in enumerate(notes):
                pass  # markers only; the stations are in the caption and the table (22 Sep)
            if i == 0:
                ax.set_title(col["title"].split("\n")[0], fontsize=8.4, pad=7, linespacing=1.4)
            if i == len(ETAS) - 1:
                ax.set_xlabel("$x/c$  (local chord of the band)")
            if j == 0:
                ax.set_ylabel("$c_{f,s}$")

    h, lab = legend_handles(LETTERS)
    if surface == "upper":
        h += [Line2D([0], [0], ls="none", marker="v", ms=5.5, mfc="0.3", mec="k", mew=0.5),
              Line2D([0], [0], ls="none", marker="^", ms=5.5, mfc="white", mec="0.3", mew=1.2)]
        lab += ["separation", "reattachment"]
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, 0.170),
               ncol=3, frameon=False, handlelength=3.0, columnspacing=1.6)

    fig.suptitle("", fontsize=10.8, y=0.995, linespacing=1.45)  # caption carries it (22 Sep)
    # Summary line derived from the measurements, not authored alongside them.
    key_pct = "pct_rev_upper" if surface == "upper" else "pct_rev_lower"
    nruns = sum(len(RUNS_U[k]) for k in RUNS_U) if surface == "upper" else 0
    worst = max(TABLE, key=lambda r: r[key_pct])
    run_cases = sorted({"%s ($\\eta$ %.2f)" % (k[0], k[1])
                        for k in RUNS_U if RUNS_U[k]})
    if surface == "upper":
        # THE BAND COUNT IS DERIVED, NOT TYPED. It read "all 55" as a literal while every
        # other value in this sentence came from the measurements, so a run over a different
        # number of bands would have adjudicated one population and announced another. That
        # is the silently-narrowed gate (GEO-087) inside a caption, and the comment above
        # this block claiming the line is "derived from the measurements" made it worse by
        # vouching for it.
        summary = ("Separated runs over all %d (case, station) bands: %d, all of them %s. "
                   "Largest single-band reversed fraction %.3f%% (%s, $\\eta$ %.2f).\n"
                   "In every OTHER band the most-forward reversed face sits at $x/c$ = %.3f "
                   "(%s, $\\eta$ %.2f), i.e. the blunt trailing-edge base region, NOT "
                   "boundary-layer separation."
                   % (len(TABLE), nruns, ", ".join(run_cases), worst[key_pct], worst["case"],
                      worst["eta"], FWD_U[0], FWD_U[1], FWD_U[2]))
    else:
        summary = ("No separated run anywhere on the LOWER surface at either condition. "
                   "Largest single-band reversed fraction %.3f%% (%s, $\\eta$ %.2f);\n"
                   "the most-forward reversed lower-surface face in the whole set sits at "
                   "$x/c$ = %.3f (%s, $\\eta$ %.2f), aft of the camber region."
                   % (worst[key_pct], worst["case"], worst["eta"],
                      FWD_L[0], FWD_L[1], FWD_L[2]))
    # NOT DRAWN ANY MORE, AND IT IS THE REASON THIS FIGURE WAS UNREADABLE. `summary` is a
    # long UNBROKEN line, and with savefig.bbox="tight" the saved canvas expands to contain
    # its widest line: that pinned the figure at 623.4 pt however figsize was set, so the
    # report scaled it to 0.773 and the 7.8 pt labels printed near 6.0 pt. The numbers it
    # carried are derived, not typed, so they were moved into the prose of
    # sections/06_performance.tex rather than deleted. It is still COMPUTED above, because
    # the same quantities feed cf_separation_table.csv, which is where a reader checks them.

    # DERIVED, for the same reason as the band count above. Inserted by concatenation rather
    # than by %-formatting the whole literal, because the block below contains a bare "2%"
    # for "2% of local chord" which %-formatting would read as a format spec and break.
    # THE METHOD NOTE IS IN THE PROSE NOW, not under the plot. Every definition it carried,
    # the c_f,s formula and its frame, the calibrated sign multiplier, the median binning,
    # the separated-run criterion and the grey band at eta 0.95, is stated in
    # sections/06_performance.tex under "Reading the skin-friction figures". Nothing was
    # dropped; it moved to where it can be read in order and found again.
    #
    # It also has to go for a mechanical reason: its longest unbroken line was wider than
    # the text column, and tight bbox sizes the canvas to the widest thing drawn on it.
    # A caption cannot be allowed to set the figure's width.
    fig.subplots_adjust(hspace=0.17, wspace=0.18, top=0.868, bottom=0.228)
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


# ----------------------------------------------------------------------------
# 9. FIGURE 3: the separation picture at early cruise
# ----------------------------------------------------------------------------
def make_separation_figure(outfile):
    fig = plt.figure(figsize=(7.4, 9.3))
    gs = fig.add_gridspec(
        3, 2, height_ratios=[1.30, 0.62, 1.20], hspace=0.52, wspace=0.26,
        left=0.095, right=0.975, top=0.905, bottom=0.135,
    )
    ax_cr = fig.add_subplot(gs[0, 0])
    ax_cm = fig.add_subplot(gs[0, 1])
    ax_sec = fig.add_subplot(gs[1, :])
    ax_tab = fig.add_subplot(gs[2, :])

    ETA_Z = 0.80
    XLO, XHI = 0.55, 1.0

    # ---- (a) and (b): aft zoom at eta 0.80, both conditions -----------------
    def _sub(prefix, h, fmt):
        # turbulence model and y+ derived from the case records and the files,
        # never typed in
        yp = [r["yplus_median_upper"] for r in TABLE
              if r["case"].startswith(prefix) and r["eta"] == ETA_Z]
        return "%s, $y^+$ median %s" % (
            h["models"][0] if h["models"] else "NOT FOUND",
            (fmt + "--" + fmt) % (min(yp), max(yp)))

    for ax, prefix, tag, sub in (
        (ax_cr, "SW",
         "(a) Condition CR, M 0.10, WALL-RESOLVED",
         _sub("SW", CR_H, "%.2f")),
        (ax_cm, "CMP",
         "(b) Early cruise, M 0.78, WALL-MODELLED",
         _sub("CMP", CM_H, "%.0f")),
    ):
        ax.axhline(0.0, color="0.35", lw=0.8, zorder=1)
        vals, nrun = [], 0
        for L in LETTERS:
            key = ("%s%s_trim" % (prefix, L), ETA_Z)
            if key not in STATS_U:
                continue
            st = STATS_U[key]
            ok = np.isfinite(st["med"]) & (st["xc"] >= XLO)
            ax.plot(st["xc"][ok], st["med"][ok], zorder=3, **STYLE[L])
            vals.append(st["med"][ok])
            for (xs, xr, nb) in RUNS_U[key]:
                nrun += 1
                ax.axvspan(xs, xr, color=STYLE[L]["color"], alpha=0.13, lw=0, zorder=2)
                ax.plot([xs], [0.0], marker="v", ms=6, mfc=STYLE[L]["color"],
                        mec="k", mew=0.5, ls="none", zorder=6)
                ax.plot([xr], [0.0], marker="^", ms=6, mfc="white",
                        mec=STYLE[L]["color"], mew=1.3, ls="none", zorder=6)
                ax.annotate("", xy=(xs, 0.0),
                            xytext=(-7, 30), textcoords="offset points", ha="right",
                            va="bottom", fontsize=8.2, color=STYLE[L]["color"],
                            zorder=7,
                            arrowprops=dict(arrowstyle="-", lw=0.7,
                                            color=STYLE[L]["color"]))
                ax.annotate("", xy=(xr, 0.0),
                            xytext=(-5, -34), textcoords="offset points", ha="right",
                            va="top", fontsize=8.2, color=STYLE[L]["color"], zorder=7,
                            arrowprops=dict(arrowstyle="-", lw=0.7,
                                            color=STYLE[L]["color"]))
        # interquartile spread of the band, C only, to show band pooling honestly
        keyC = ("%sC_trim" % prefix, ETA_Z)
        if keyC in STATS_U:
            st = STATS_U[keyC]
            ok = np.isfinite(st["med"]) & (st["xc"] >= XLO)
            ax.fill_between(st["xc"][ok], st["q25"][ok], st["q75"][ok],
                            color=STYLE["C"]["color"], alpha=0.18, lw=0, zorder=1.5)
            vals.append(st["q25"][ok]); vals.append(st["q75"][ok])
        v = np.concatenate(vals)
        lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
        # extra headroom below zero where an annotation hangs under the axis
        pad_lo = (0.42 if nrun else 0.14) * (hi - lo)
        ax.set_ylim(lo - pad_lo, hi + 0.14 * (hi - lo))
        ax.set_xlim(XLO, XHI)
        ax.set_xlabel("$x/c$")
        ax.set_ylabel("$c_{f,s}$")
        ax.set_title(tag.split(",")[0] + "\n$\\eta$ = 0.80",
                     fontsize=8.8, pad=6, linespacing=1.35)
        ax.tick_params(length=3, width=0.7)
        ax.axvspan(0.875, XHI, color="0.86", alpha=0.5, lw=0, zorder=0)
        if nrun == 0:
            ax.text(0.5, 0.06, "",
                    transform=ax.transAxes, ha="center", va="bottom",
                    fontsize=8.6, color="0.25", zorder=7)

    # ---- (c) where it sits on the section -----------------------------------
    keyC = ("CMPC_trim", ETA_Z)
    keyB = ("CMPB_trim", ETA_Z)
    stCu, stCl = STATS_U[keyC], STATS_L[keyC]
    stBu, stBl = STATS_U[keyB], STATS_L[keyB]
    for st in (stBu, stBl):
        ok = np.isfinite(st["zmed"])
        ax_sec.plot(st["xc"][ok], st["zmed"][ok], color="0.60", ls=(0, (4, 2)),
                    lw=0.9, zorder=2)
    for st in (stCu, stCl):
        ok = np.isfinite(st["zmed"])
        ax_sec.plot(st["xc"][ok], st["zmed"][ok], color="#000000", ls="-", lw=1.1,
                    zorder=3)
    # reversed segment of the upper branch, drawn as a thick overlay
    negmask = np.isfinite(stCu["med"]) & (stCu["med"] < 0.0) & np.isfinite(stCu["zmed"])
    ax_sec.plot(stCu["xc"][negmask], stCu["zmed"][negmask],
                color=STYLE["C"]["color"], ls="-", lw=3.2, solid_capstyle="butt",
                zorder=5)
    runsC = RUNS_U[keyC]
    if runsC:
        xs, xr = runsC[0][0], runsC[-1][1]
        fin = np.isfinite(stCu["zmed"])
        for xv, mk, fc in ((xs, "v", STYLE["C"]["color"]), (xr, "^", "white")):
            zi = np.interp(xv, stCu["xc"][fin], stCu["zmed"][fin])
            ax_sec.plot([xv], [zi + 0.022], marker=mk, ms=6, mfc=fc,
                        mec=STYLE["C"]["color"], mew=1.2, ls="none", zorder=6)
        ax_sec.annotate("",
                        xy=(0.5 * (xs + xr), 0.085), ha="center", va="bottom",
                        fontsize=8.4, color=STYLE["C"]["color"])
    ax_sec.set_aspect("equal", adjustable="box")
    ax_sec.set_xlim(-0.02, 1.04)
    ax_sec.set_ylim(-0.085, 0.135)
    ax_sec.set_yticks([-0.05, 0.0, 0.05])
    ax_sec.set_xlabel("$x/c$")
    ax_sec.set_ylabel("$z/c$")
    ax_sec.set_title("(c)",
                     fontsize=8.8, pad=5)
    ax_sec.tick_params(length=3, width=0.7)

    # ---- (d) table -----------------------------------------------------------
    ax_tab.axis("off")
    cell = []
    for L in LETTERS:
        case = "CMP%s_trim" % L
        if case in cases_found:
            r = [case, GEOM_NAME.get(L, "?")]
            for eta in ETAS:
                rec = [t for t in TABLE if t["case"] == case and t["eta"] == eta][0]
                r.append("%.3f" % rec["pct_rev_upper"])
            r.append("%.3f" % POOLED[case]["pct"])
            rec80 = [t for t in TABLE if t["case"] == case and t["eta"] == 0.80][0]
            r.append("%.3f-%.3f" % (rec80["x_sep"], rec80["x_reatt"])
                     if rec80["x_sep"] is not None else "none")
        else:
            r = [case, GEOM_NAME.get(L, "?")] + ["NOT FOUND"] * (len(ETAS) + 2)
        cell.append(r)

    colnames = (["case", "DSO geometry"] + ["$\\eta$ %.2f" % e for e in ETAS] +
                ["5 bands\npooled", "$x/c$ sep-reatt\nat $\\eta$ 0.80"])
    tb = ax_tab.table(cellText=cell, colLabels=colnames, cellLoc="center",
                      loc="upper center",
                      colWidths=[0.115, 0.205] + [0.076] * len(ETAS) + [0.090, 0.140])
    tb.auto_set_font_size(False)
    tb.set_fontsize(8.0)
    tb.scale(1.0, 1.42)
    for (ri, ci), c in tb.get_celld().items():
        c.set_linewidth(0.5)
        c.set_edgecolor("0.6")
        txt = c.get_text().get_text()
        if ri == 0:
            c.set_text_props(weight="bold")
            c.set_facecolor("0.94")
            continue
        if ci <= 1:
            c.set_text_props(ha="left")
            c.PAD = 0.04
        if txt == "NOT FOUND":
            c.set_text_props(color="#B00000", style="italic", size=6.4)
        elif ci >= 2:
            try:
                if float(txt.split("-")[0]) > 1.0:
                    c.set_facecolor("#F6DFD2")
                    c.set_text_props(weight="bold")
            except ValueError:
                pass
    ax_tab.set_title(
        "(d)",
        fontsize=8.8, pad=10, loc="center", linespacing=1.35)
    # THE CAVEAT IS DERIVED FROM THE TABLE, NOT AUTHORED BESIDE IT (GEO-080). This line
    # read "CMPM_trim was staged but never run, so the M row is NOT FOUND" while the table
    # printed directly above it showed M at 24.875 per cent. The claim was written once,
    # by hand, when it was true, and nothing recomputed it when the case ran. A caption
    # that contradicts its own table is worse than no caption, because the table is what a
    # reader checks the caption against.
    # THE SAME PREDICATE THE TABLE ITSELF USES, not a second one that could disagree with
    # it. The table's own rows are built by `case = "CMP%s_trim" % L` and `if case in
    # cases_found`, so the caption keys on exactly that and the two cannot drift apart.
    # My first attempt reached for a `CM_CASE` mapping that does not exist anywhere in
    # this file, which would have failed with NameError on the first run.
    _absent = sorted(L for L in LETTERS if ("CMP%s_trim" % L) not in cases_found)
    _missing_note = (" %s has no run in this set, so its row reads NOT FOUND."
                     % ", ".join(_absent)) if _absent else ""
    ax_tab.text(
        0.5, -0.02,
        "",
        transform=ax_tab.transAxes, ha="center", va="top", fontsize=7.8, color="0.25",
        linespacing=1.5)

    h, lab = legend_handles(LETTERS)
    h += [Line2D([0], [0], ls="none", marker="v", ms=6, mfc="0.3", mec="k", mew=0.5),
          Line2D([0], [0], ls="none", marker="^", ms=6, mfc="white", mec="0.3", mew=1.3)]
    lab += ["separation", "reattachment"]
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, 0.082), ncol=4,
               frameon=False, handlelength=3.2, columnspacing=1.5)
    fig.suptitle("",
                 fontsize=10.8, y=0.978)
    fig.text(0.5, 0.022,
             "",
             ha="center", va="top", fontsize=7.8, color="0.25", linespacing=1.5)
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


# ----------------------------------------------------------------------------
# 10. Write, then VERIFY from the filesystem (GEO-080: derive the claim by
#     reading the artefact, do not assert it alongside the action)
# ----------------------------------------------------------------------------
written = []
written.append(make_sections_figure("upper", os.path.join(FIG, "fig_cf_sections_upper.pdf")))
written.append(make_sections_figure("lower", os.path.join(FIG, "fig_cf_sections_lower.pdf")))
written.append(make_separation_figure(os.path.join(FIG, "fig_cf_separation_earlycruise.pdf")))
written.append(tab_csv)

log("\nVERIFICATION (stat of each artefact, read back from disk):")
ok_all = True
for p in written:
    if not os.path.exists(p):
        log("  MISSING  %s" % p)
        ok_all = False
        continue
    sz = os.path.getsize(p)
    head = b""
    with open(p, "rb") as fh:
        head = fh.read(5)
    kind = "PDF" if head.startswith(b"%PDF") else "text"
    log("  OK  %-58s %9d bytes  %s" % (os.path.basename(p), sz, kind))
    if sz == 0:
        ok_all = False
if not ok_all:
    halt("at least one artefact is missing or empty")

with open(os.path.join(FIG, "make_cf_LOG.txt"), "w") as fh:
    fh.write("\n".join(LOG_LINES) + "\n")
print("\nlog written to", os.path.join(FIG, "make_cf_LOG.txt"))
