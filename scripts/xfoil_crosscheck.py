#!/usr/bin/env python3
"""xfoil_crosscheck.py: task-4 CODE-TO-CODE cross-check, XFOIL vs 2D RANS.

Runs XFOIL 6.99 viscous polars on the faired EET baseline section at the
Condition CR Reynolds number (chord-based, same as the RANS 2D cases) and
overlays cl(alpha) and cd(cl) against the fine-grid OpenFOAM points.

STRICTLY code-to-code: two different physics models (e^n transitional IBL vs
fully-turbulent RANS kOmegaSST), no experimental data, no accuracy claim in
either direction. Two XFOIL curves are produced:
  1. forced transition at x/c 0.05 both surfaces (closest analogue to the
     fully-turbulent RANS),
  2. free transition Ncrit 9 (XFOIL's natural state, for information).

Usage:
  python3 scripts/xfoil_crosscheck.py --xfoil /path/to/xfoil \
      [--out results/eet2d_cr/xfoil]
"""

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
GEOM = REPO / "geometry/derived/eet_section/eet_cruise_faired.dat"

U, CHORD, NU = 34.0, 0.3404, 1.46e-5   # Condition CR (docs/ARGUS_reference_data.md s9)
RE = U * CHORD / NU
NPANEL = 320
ALFA_SEQ = (-2.0, 6.0, 0.5)
CL_BAND = (0.40, 0.53)   # operating cl range of the ARGUS 2D program

C_SURFACE = "#fcfcfb"; C_TEXT = "#0b0b0b"; C_TEXT2 = "#52514e"; C_GRID = "#e8e7e3"
C_RANS = "#1baf7a"; C_XT = "#2a78d6"; C_XF = "#a06be0"


def xfoil_polar(xfoil, workdir, polar_name, forced_xtr, aseq=ALFA_SEQ):
    """One PACC polar over `aseq`; returns parsed rows."""
    polar = workdir / polar_name
    polar.unlink(missing_ok=True)   # XFOIL appends to an existing polar file
    lines = ["PLOP", "G", "", "LOAD eet.dat", "PANE", "PPAR", f"N {NPANEL}", "", "",
             "OPER", f"VISC {RE:.0f}", "MACH 0", "ITER 300"]
    if forced_xtr is not None:
        lines += ["VPAR", f"XTR {forced_xtr} {forced_xtr}", ""]
    lines += ["PACC", polar_name, "", f"ASEQ {aseq[0]} {aseq[1]} {aseq[2]}",
              "PACC", "", "QUIT"]
    proc = subprocess.run([str(xfoil)], input="\n".join(lines) + "\n", text=True,
                          capture_output=True, cwd=workdir, timeout=600,
                          env={"PATH": "/usr/bin", "HOME": str(workdir)})
    if not polar.exists():
        raise RuntimeError(f"XFOIL produced no polar {polar_name}:\n{proc.stdout[-2000:]}")
    rows = []
    for ln in polar.read_text().splitlines():
        parts = ln.split()
        if len(parts) >= 7:
            try:
                vals = [float(v) for v in parts[:7]]
            except ValueError:
                continue
            rows.append(dict(zip(("alpha_deg", "cl", "cd", "cdp", "cm",
                                  "xtr_top", "xtr_bot"), vals)))
    return rows, polar.read_text()


def panel_sensitivity(xfoil, workdir, counts=(200, 260, 320, 360)):
    """XFOIL-side discretization check: the code-to-code analogue of the RANS
    grid study. Without it, 'the codes agree to X counts' has an unquantified
    term on the XFOIL side."""
    rows = []
    for n in counts:
        global NPANEL
        keep, NPANEL = NPANEL, n
        try:
            pol, _ = xfoil_polar(xfoil, workdir, f"pol_n{n}.txt", 0.05,
                                 aseq=(0.0, 4.0, 2.0))
        finally:
            NPANEL = keep
        for p in pol:
            rows.append({"npanel": n, **p})
    spread = {}
    for a in sorted({r["alpha_deg"] for r in rows}):
        at = [r for r in rows if r["alpha_deg"] == a]
        spread[f"alpha_{a:g}"] = {
            "cl_spread": round(max(r["cl"] for r in at) - min(r["cl"] for r in at), 5),
            "cd_spread_counts": round((max(r["cd"] for r in at)
                                       - min(r["cd"] for r in at)) * 1e4, 3)}
    return rows, spread


def rans_points():
    """Fine-grid baseline points: alpha sweep cases + grid-family fine + trims.
    Trim cases report true-wind-axis values (D019 amendment 1) via their cards."""
    pts = []
    for card_path in sorted(REPO.glob("cases/baseline/CR2d_*fine*/run_card.json")):
        card = json.loads(card_path.read_text())
        if card["case"]["condition"] != "CR":
            continue
        pts.append({"alpha_deg": card["conditions"]["alpha_deg"],
                    "cl": card["results"]["CL"], "cd": card["results"]["CD"],
                    "cdp": card["results"]["CD_pressure"], "cdv": card["results"]["CD_viscous"],
                    "cm": card["results"]["Cm"], "case": card["case"]["path"]})
    return sorted(pts, key=lambda p: p["alpha_deg"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--xfoil", required=True, help="path to the xfoil binary")
    ap.add_argument("--out", type=Path, default=Path("results/eet2d_cr/xfoil"))
    args = ap.parse_args()
    out = REPO / args.out
    out.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="xfoil_") as td:
        wd = Path(td)
        shutil.copy(GEOM, wd / "eet.dat")
        forced, forced_raw = xfoil_polar(args.xfoil, wd, "pol_forced.txt", 0.05)
        free, free_raw = xfoil_polar(args.xfoil, wd, "pol_free.txt", None)
        pan_rows, pan_spread = panel_sensitivity(args.xfoil, wd)
    (out / "xfoil_polar_forced_xtr0p05.txt").write_text(forced_raw)
    (out / "xfoil_polar_free_ncrit9.txt").write_text(free_raw)
    with (out / "xfoil_panel_sensitivity.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pan_rows[0].keys()))
        w.writeheader(); w.writerows(pan_rows)

    rans = rans_points()
    for name, rows in (("xfoil_forced_xtr0p05.csv", forced), ("xfoil_free_ncrit9.csv", free),
                       ("rans_fine_points.csv", rans)):
        with (out / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)

    # --- comparisons. Fixed-alpha is the raw diff; fixed-cl is the meaningful
    # drag comparison (the whole 2D program is a fixed-cl program).
    fa = np.array([r["alpha_deg"] for r in forced])
    fcl = np.array([r["cl"] for r in forced])
    fcd = np.array([r["cd"] for r in forced])
    at_alpha, at_cl = [], []
    for p in rans:
        if not (fa.min() <= p["alpha_deg"] <= fa.max()):
            continue
        xcl = float(np.interp(p["alpha_deg"], fa, fcl))
        xcd = float(np.interp(p["alpha_deg"], fa, fcd))
        at_alpha.append({"alpha_deg": p["alpha_deg"], "rans_cl": p["cl"], "xfoil_cl": xcl,
                         "dcl": round(xcl - p["cl"], 5), "rans_cd": p["cd"], "xfoil_cd": xcd,
                         "dcd_counts": round((xcd - p["cd"]) * 1e4, 2)})
        if fcl.min() <= p["cl"] <= fcl.max():
            xcd_cl = float(np.interp(p["cl"], fcl, fcd))   # cl monotonic over this range
            at_cl.append({"cl": p["cl"], "rans_alpha_deg": p["alpha_deg"],
                          "xfoil_alpha_deg": round(float(np.interp(p["cl"], fcl, fa)), 4),
                          "rans_cd": p["cd"], "xfoil_cd": xcd_cl,
                          "dcd_counts": round((xcd_cl - p["cd"]) * 1e4, 2)})

    def slope(alphas, cls, lo=-2.0, hi=4.0):
        a, c = np.asarray(alphas), np.asarray(cls)
        m = (a >= lo) & (a <= hi)
        if m.sum() < 2:
            return None, None
        k, b = np.polyfit(a[m], c[m], 1)
        return round(float(k), 5), round(float(-b / k), 4)

    r_slope, r_a0 = slope([p["alpha_deg"] for p in rans], [p["cl"] for p in rans])
    f_slope, f_a0 = slope(fa, fcl)
    n_slope, n_a0 = slope([r["alpha_deg"] for r in free], [r["cl"] for r in free])

    summary = {
        "label": "CODE-TO-CODE ONLY: XFOIL 6.99 (e^n transitional integral BL, inviscid "
                 "panel + viscous coupling) vs OpenFOAM-org 7 simpleFoam kOmegaSST "
                 "fully-turbulent low-Re RANS. Two different physics models; NO "
                 "experimental reference is involved and neither result is evidence of "
                 "accuracy for the other.",
        "geometry": "geometry/derived/eet_section/eet_cruise_faired.dat (D004 faired baseline)",
        "conditions": {"condition": "CR", "Re_chord": RE, "U_m_s": U, "chord_m": CHORD,
                       "mach": "0 in XFOIL; incompressible RANS (nominal M 0.10)"},
        "xfoil_setup": {"version": "6.99 (MIT source, single precision)",
                        "panels": NPANEL, "iter": 300,
                        "primary_curve": "forced transition x/c 0.05 both surfaces "
                                         "(closest analogue to fully-turbulent RANS)",
                        "secondary_curve": "free transition Ncrit 9 (informational only)",
                        "panel_sensitivity": pan_spread,
                        "panel_sensitivity_note": "spread over N = 200/260/320/360 panels; "
                                                  "the XFOIL side of the comparison is "
                                                  "panel-converged well inside 0.1 counts"},
        "n_rans_points": len(rans),
        "lift_slope_per_deg_fit_-2_to_4": {"rans": r_slope, "xfoil_forced": f_slope,
                                           "xfoil_free": n_slope},
        "alpha_at_cl0_deg": {"rans": r_a0, "xfoil_forced": f_a0, "xfoil_free": n_a0},
        "at_matched_alpha_vs_forced": at_alpha,
        "at_matched_cl_vs_forced": at_cl,
    }
    if at_cl:
        d = [x["dcd_counts"] for x in at_cl]
        band = [x for x in at_cl if CL_BAND[0] <= x["cl"] <= CL_BAND[1]]
        summary["headline_at_matched_cl"] = {
            "dcd_counts_min": min(d), "dcd_counts_max": max(d),
            "dcd_counts_mean": round(sum(d) / len(d), 2),
            "n_points": len(at_cl),
            "reading": "XFOIL(forced) minus RANS, counts, at matched cl"}
        if band:
            b = [x["dcd_counts"] for x in band]
            summary["headline_at_matched_cl"]["operating_band"] = {
                "cl_range": list(CL_BAND), "n_points": len(band),
                "dcd_counts": b, "max_abs_dcd_counts": round(max(abs(v) for v in b), 2),
                "note": "the cl range the ARGUS 2D program actually operates in "
                        "(DSO sectional targets 0.4221-0.4777 plus the grid-study point)"}
    (out / "crosscheck_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6), facecolor=C_SURFACE)
    for ax in axes:
        ax.set_facecolor(C_SURFACE); ax.grid(True, color=C_GRID, lw=0.7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(C_TEXT2)
        ax.tick_params(colors=C_TEXT2, labelsize=9)
    axes[0].plot([r["alpha_deg"] for r in forced], [r["cl"] for r in forced], "-",
                 color=C_XT, lw=1.6, label="XFOIL, forced xtr 0.05")
    axes[0].plot([r["alpha_deg"] for r in free], [r["cl"] for r in free], "--",
                 color=C_XF, lw=1.3, label="XFOIL, free Ncrit 9")
    axes[0].plot([p["alpha_deg"] for p in rans], [p["cl"] for p in rans], "o",
                 color=C_RANS, ms=7, label="RANS SST fine O-grid")
    axes[0].set_xlabel("alpha, deg", color=C_TEXT2, fontsize=9)
    axes[0].set_ylabel("cl", color=C_TEXT2, fontsize=9)
    axes[0].legend(fontsize=8.5, frameon=False, labelcolor=C_TEXT, loc="lower right")
    axes[1].plot([r["cd"] * 1e4 for r in forced], [r["cl"] for r in forced], "-",
                 color=C_XT, lw=1.6)
    axes[1].plot([r["cd"] * 1e4 for r in free], [r["cl"] for r in free], "--",
                 color=C_XF, lw=1.3)
    axes[1].plot([p["cd"] * 1e4 for p in rans], [p["cl"] for p in rans], "o",
                 color=C_RANS, ms=7)
    axes[1].set_xlabel("cd, counts", color=C_TEXT2, fontsize=9)
    axes[1].set_ylabel("cl", color=C_TEXT2, fontsize=9)
    axes[1].set_title("full range: the free-transition laminar bucket",
                      color=C_TEXT2, fontsize=8.5)

    # zoom: the comparison band, which the laminar bucket otherwise compresses
    rcd = [p["cd"] * 1e4 for p in rans]
    lo, hi = min(rcd) - 4, max(rcd) + 6
    axes[2].plot([r["cd"] * 1e4 for r in forced], [r["cl"] for r in forced], "-",
                 color=C_XT, lw=1.6)
    axes[2].plot([r["cd"] * 1e4 for r in free], [r["cl"] for r in free], "--",
                 color=C_XF, lw=1.3)
    axes[2].plot(rcd, [p["cl"] for p in rans], "o", color=C_RANS, ms=7)
    axes[2].set_xlim(lo, hi)
    axes[2].set_ylim(min(p["cl"] for p in rans) - 0.06, max(p["cl"] for p in rans) + 0.06)
    axes[2].set_xlabel("cd, counts", color=C_TEXT2, fontsize=9)
    axes[2].set_ylabel("cl", color=C_TEXT2, fontsize=9)
    axes[2].set_title("zoom on the RANS points", color=C_TEXT2, fontsize=8.5)
    fig.suptitle("EET faired baseline, Condition CR (Re 0.79e6): XFOIL vs RANS - "
                 "CODE-TO-CODE ONLY (no experimental reference)",
                 color=C_TEXT, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out / "xfoil_crosscheck.png", dpi=200, facecolor=C_SURFACE,
                bbox_inches="tight")
    print(f"RANS points: {len(rans)}; XFOIL forced pts: {len(forced)}; free pts: {len(free)}")
    print(f"lift slope /deg (fit -2..4): RANS {r_slope}, XFOIL forced {f_slope}, "
          f"free {n_slope}; alpha(cl=0): {r_a0} / {f_a0} / {n_a0}")
    print("at matched cl (XFOIL forced minus RANS):")
    for x in at_cl:
        print(f"  cl {x['cl']:.4f}: RANS {x['rans_cd']*1e4:6.2f} cts (a {x['rans_alpha_deg']:+.3f}) "
              f"| XFOIL {x['xfoil_cd']*1e4:6.2f} cts (a {x['xfoil_alpha_deg']:+.3f}) "
              f"| dcd {x['dcd_counts']:+.2f}")
    print(f"wrote {out}/")


if __name__ == "__main__":
    main()
