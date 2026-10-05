#!/usr/bin/env python3
"""collect_trim_deltas.py: task-3 deliverable - fixed-cl morphed-section deltas.

Pairs the trimmed baseline and morphed 2D cases per eta station (targets =
DSO .lod baseline sectional cls, D017), writes schema-valid run cards with
the trim block from trim_history.json, and produces the delta table:
delta cd with pressure/viscous split per eta at matched sectional cl on
matched (same-generator, same-level) meshes.

Outputs: results/eet2d_cr/trim_deltas.csv / .json, delta figure, run cards.
"""

import csv
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from collect_2d_grid import parse_case, run_card  # noqa: E402
sys.path.insert(0, str(REPO / "validation_data/naca0012_tmr"))
from collect_results import history  # noqa: E402

PAIRS = [
    (0.70, "cases/baseline/CR2d_trimb_eta0p70_fine_a0p92",
           "cases/mbr_c01_l011/CR2d_trimm_eta0p70_fine_am0p41", 0.4777),
    (0.80, "cases/baseline/CR2d_trimb_eta0p80_fine_a0p77",
           "cases/mbr_c01_l011/CR2d_trimm_eta0p80_fine_am1", 0.4609),
    (0.90, "cases/baseline/CR2d_trimb_eta0p90_fine_a0p42",
           "cases/mbr_c01_l011/CR2d_trimm_eta0p90_fine_am0p81", 0.4221),
]

C_SURFACE = "#fcfcfb"; C_TEXT = "#0b0b0b"; C_TEXT2 = "#52514e"; C_GRID = "#e8e7e3"
C_P = "#2a78d6"; C_V = "#eb6834"; C_TOT = "#1baf7a"


def true_axis(case, r):
    """Recompute cl/cd (with p/v split) in the TRUE wind axes, defined by the
    U freestream vector (authoritative record of the final block's flow angle).
    Guards against stale forceCoeffs liftDir/dragDir: the baseline trim batch
    ran its final block with axes one secant update behind the U rotation
    (trim2d controlDict regex no-op, fixed 2026-07-25), biasing forceCoeffs cd
    by ~cl x axis lag. Cross-checks the forcesAll vector route against the
    rotated forceCoeffs values; cm is rotation-invariant and kept as parsed."""
    m = re.search(r"freestreamValue uniform \(\s*([-0-9eE+.]+)\s+([-0-9eE+.]+)\s+([-0-9eE+.]+)\s*\)",
                  (case / "0.orig/U").read_text())
    alpha_u = math.degrees(math.atan2(float(m.group(3)), float(m.group(1))))
    dd = re.search(r"dragDir\s+\(\s*([-0-9eE+.]+)\s+[-0-9eE+.]+\s+([-0-9eE+.]+)\s*\)",
                   (case / "system/controlDict").read_text())
    alpha_dirs = math.degrees(math.atan2(float(dd.group(2)), float(dd.group(1))))

    meta = r["meta"]
    a = math.radians(alpha_u)
    drag = np.array([math.cos(a), 0.0, math.sin(a)])
    lift = np.array([-math.sin(a), 0.0, math.cos(a)])
    qref = 0.5 * 1.225 * meta["U_m_s"] ** 2 * meta["chord_m"]
    last = history(case, "forcesAll", "forces*.dat")[-1]
    t_forces = float(last.split()[0])
    vp, vv = [np.array([float(v) for v in g.split()])
              for g in re.findall(r"\(([-0-9.eE+ ]+)\)", last)[:2]]
    cdp, cdv = float(vp @ drag / qref), float(vv @ drag / qref)
    cl = float((vp + vv) @ lift / qref)

    # forcesAll writes on writeTime, forceCoeffs every 10 iterations, so the two
    # routes must be compared AT THE SAME TIME (a still-evolving tail otherwise
    # reads as a bogus axis disagreement). The p/v split necessarily comes from
    # the last forcesAll write, so the run must also be stationary between that
    # write and its final time, else the split would not describe the end state.
    fc = np.loadtxt(history(case, "forces", "forceCoeffs*.dat"))
    i = int(np.argmin(np.abs(fc[:, 0] - t_forces)))
    cd_fc, cl_fc, t_fc = float(fc[i, 2]), float(fc[i, 3]), float(fc[i, 0])
    delta = math.radians(alpha_u - alpha_dirs)
    cd_rot = cd_fc * math.cos(delta) + cl_fc * math.sin(delta)
    gap_cts = abs(cd_rot - (cdp + cdv)) * 1e4
    tail_cts = abs(float(fc[-1, 2]) - cd_fc) * 1e4
    # forceCoeffs samples on a fixed interval; a writeNow stop lands between
    # samples, so match the nearest one within one interval rather than exactly.
    step = float(np.median(np.diff(fc[:, 0]))) if len(fc) > 2 else 1.0
    assert abs(t_fc - t_forces) <= 1.5 * step, (
        f"{case.name}: nearest forceCoeffs sample is t={t_fc:g}, "
        f"{abs(t_fc - t_forces):g} iterations from the forcesAll write at "
        f"t={t_forces:g} (sample interval {step:g})")
    assert gap_cts < 0.05, (f"{case.name}: forcesAll-vs-rotated-forceCoeffs cd "
                            f"disagree by {gap_cts:.3f} counts at t={t_forces:g}")
    assert tail_cts < 0.05, (
        f"{case.name}: cd moved {tail_cts:.3f} counts between the last forcesAll "
        f"write (t={t_forces:g}) and the final time (t={fc[-1, 0]:g}); the p/v split "
        f"would not describe the end state. Let the case reach a write time.")
    return {"alpha_U_deg": round(alpha_u, 6), "alpha_forceCoeffs_dirs_deg": round(alpha_dirs, 6),
            "axis_lag_deg": round(alpha_u - alpha_dirs, 6),
            "cl": cl, "cd": cdp + cdv, "cdp": cdp, "cdv": cdv,
            "crosscheck_gap_counts": round(gap_cts, 4),
            "tail_drift_counts": round(tail_cts, 4),
            "time": t_forces, "final_time": float(fc[-1, 0])}


def aft_attachment(case, meta):
    """Aft-region (x/c 0.5-0.997) attachment from the D014 wallShearStress raw
    output at the last write time. org-7 wallShearStress is the stress ON THE
    FLUID: attached +x flow gives tau_x < 0; reversal shows as tau_x > 0.
    Faces are assigned to a surface by proximity to the meshed .dat surfaces
    (a z-sign test fails on drooped-TE morphed sections)."""
    pts = np.array([[float(v) for v in ln.split()] for ln in
                    (REPO / meta["geometry"]).read_text().splitlines()[1:]])
    ile = int(np.argmin(pts[:, 0]))
    up, lo = pts[:ile + 1][::-1], pts[ile:]

    raws = sorted((case / "postProcessing/surfaceFields").glob("*/wallShearStress_airfoil.raw"),
                  key=lambda p: float(p.parent.name))
    d = np.loadtxt(raws[-1], comments="#")
    x, z, tau_x = d[:, 0] / meta["chord_m"], d[:, 2] / meta["chord_m"], d[:, 3]
    zu = np.interp(np.clip(x, 0, 1), up[:, 0], up[:, 1])
    zl = np.interp(np.clip(x, 0, 1), lo[:, 0], lo[:, 1])
    is_up = np.abs(z - zu) < np.abs(z - zl)

    out = {}
    for name, mask in (("upper", is_up), ("lower", ~is_up)):
        aft = mask & (x > 0.5) & (x < 0.997)
        rev = tau_x[aft] > 1e-6
        out[name] = {"attached": not bool(rev.any()), "n_aft_faces": int(aft.sum()),
                     "n_reversed_faces": int(rev.sum()),
                     "max_taux_m2_s2": float(tau_x[aft].max()),
                     "x_c_first_reversal": float(np.sort(x[aft][rev])[0]) if rev.any() else None}
    return out


def card_with_trim(case, r, eta, target):
    ax = true_axis(case, r)
    r.update(cl=ax["cl"], cd=ax["cd"], cdp=ax["cdp"], cdv=ax["cdv"])
    assert abs(ax["cl"] - target) <= 1.5e-4, \
        f"{case.name}: true-axis cl {ax['cl']:.6f} misses target {target}"
    run_card(case, r)
    card = json.loads((case / "run_card.json").read_text())
    card["results"]["wind_axis_check"] = ax
    card["case"]["eta"] = eta
    card["case"]["description"] = (f"Task-3 fixed-cl trim at eta {eta} "
                                   f"(target sectional cl {target} from the DSO .lod, D017)")
    trim = json.loads((case / "trim_history.json").read_text())
    card["trim"] = {"target_CL": trim["target_CL"], "tolerance": trim["tolerance"],
                    "alpha_trim_deg": trim["alpha_trim_deg"],
                    "history": trim["history"], "converged": trim["converged"]}
    card["conditions"]["alpha_deg"] = trim["alpha_trim_deg"]
    card["results"]["aft_attachment"] = aft_attachment(case, r["meta"])
    card["notes"] = ("Task 3 (session-2): fixed-cl comparison per eta at the DSO baseline "
                    "sectional cl; matched D015 O-grid fine meshes; delta table in "
                    "results/eet2d_cr/trim_deltas.csv. Coefficients computed in TRUE wind "
                    "axes from the U freestream vector (results.wind_axis_check); the "
                    "baseline batch's forceCoeffs axes lagged one secant update "
                    "(D019 amendment 2026-07-25).")
    (case / "run_card.json").write_text(json.dumps(card, indent=2) + "\n")
    return card


def main():
    out = REPO / "results/eet2d_cr"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for eta, bpath, mpath, target in PAIRS:
        b, m = REPO / bpath, REPO / mpath
        rb, rm = parse_case(b), parse_case(m)
        cb = card_with_trim(b, rb, eta, target)
        cm = card_with_trim(m, rm, eta, target)
        att_b = cb["results"]["aft_attachment"]
        att_m = cm["results"]["aft_attachment"]
        tb = json.loads((b / "trim_history.json").read_text())
        tm = json.loads((m / "trim_history.json").read_text())
        rows.append({
            "eta": eta, "target_cl": target,
            "cl_baseline": rb["cl"], "cl_morphed": rm["cl"],
            "alpha_baseline_deg": tb["alpha_trim_deg"], "alpha_morphed_deg": tm["alpha_trim_deg"],
            "dalpha_deg": round(tm["alpha_trim_deg"] - tb["alpha_trim_deg"], 5),
            "cd_baseline_counts": round(rb["cd"] * 1e4, 3),
            "cd_morphed_counts": round(rm["cd"] * 1e4, 3),
            "dcd_counts": round((rm["cd"] - rb["cd"]) * 1e4, 3),
            "dcd_pressure_counts": round((rm["cdp"] - rb["cdp"]) * 1e4, 3),
            "dcd_viscous_counts": round((rm["cdv"] - rb["cdv"]) * 1e4, 3),
            "dcm": round(rm["cm"] - rb["cm"], 5),
            "trim_converged": tb["converged"] and tm["converged"],
            "aft_attached_baseline": att_b["upper"]["attached"] and att_b["lower"]["attached"],
            "aft_attached_morphed": att_m["upper"]["attached"] and att_m["lower"]["attached"],
        })
        print(f"eta {eta}: cd {rb['cd']*1e4:.2f} -> {rm['cd']*1e4:.2f} counts "
              f"(dcd {rows[-1]['dcd_counts']:+.2f}: p {rows[-1]['dcd_pressure_counts']:+.2f}, "
              f"v {rows[-1]['dcd_viscous_counts']:+.2f}); dalpha {rows[-1]['dalpha_deg']:+.3f} deg")

    with (out / "trim_deltas.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    (out / "trim_deltas.json").write_text(json.dumps(rows, indent=2) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.4), facecolor=C_SURFACE)
    ax.set_facecolor(C_SURFACE); ax.grid(True, color=C_GRID, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(C_TEXT2)
    ax.tick_params(colors=C_TEXT2, labelsize=9)
    etas = [r["eta"] for r in rows]
    ax.plot(etas, [r["dcd_counts"] for r in rows], "o-", color=C_TOT, ms=8, lw=1.6, label="total")
    ax.plot(etas, [r["dcd_pressure_counts"] for r in rows], "s--", color=C_P, ms=7, label="pressure")
    ax.plot(etas, [r["dcd_viscous_counts"] for r in rows], "^--", color=C_V, ms=7, label="viscous")
    ax.axhline(0, color=C_TEXT2, lw=0.8)
    ax.set_xlabel("eta station", color=C_TEXT2, fontsize=9)
    ax.set_ylabel("delta cd, counts (morphed - baseline)", color=C_TEXT2, fontsize=9)
    ax.legend(fontsize=9, frameon=False, labelcolor=C_TEXT)
    ax.set_title("2D profile-drag penalty of aft camber at fixed sectional cl\n"
                 "(Condition CR, DSO .lod targets, matched fine O-grids)",
                 color=C_TEXT, fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "trim_deltas.png", dpi=200, facecolor=C_SURFACE, bbox_inches="tight")
    print(f"wrote {out/'trim_deltas.csv'} / .json / .png")

    cards = [str(REPO / p / "run_card.json") for _, bp, mp, _ in PAIRS for p in (bp, mp)]
    subprocess.run([sys.executable, str(REPO / "scripts/validate_run_card.py")] + cards,
                   cwd=REPO, check=True)


if __name__ == "__main__":
    main()
