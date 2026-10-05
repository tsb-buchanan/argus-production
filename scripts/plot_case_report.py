#!/usr/bin/env python3
"""plot_case_report.py: the full figure set for one case -- Cp, Cf, y+, forces, U and p.

WHY (requirement, 2026-08-09): Cp and Cf data to plot, and the final time step of U and p,
both to validate the flow and for the report.

scripts/validation_figures.py already does Cp and Cf, but it REQUIRES a reference to compare
against. The ARGUS sections have none: they are a design geometry, not a validation case.
Demanding a reference would mean either no plots or a fabricated one, so this plots the case
on its own terms and adds the two things a design case needs that a validation case does not:
y+ along the surface, and the final field.

WHAT EACH PANEL IS FOR, since a figure with no reference has to justify itself:
    Cp        where the shock sits, and whether the recovery is smooth or separated
    Cf        SIGN is the diagnostic: Cf <= 0 on the upper surface IS separation, and at a
              transonic condition that is the failure mode worth looking for
    y+        the wall regime actually achieved, per station, not the run-average
    forces    cl and cd against iteration with the D023 window shaded, so "converged" is
              visible rather than asserted
    U, p      the final field: shock position, its extent off the surface, and any wake
              thickening that the surface plots cannot show

NORMALISATION COMES FROM scripts/case_normalisation.py, hash-guarded against the case files.
Nothing here recomputes a denominator: for a COMPRESSIBLE solver OpenFOAM writes
wallShearStress rho-weighted (Pa), and using the kinematic form made our Cf read 1.1664x high
for a day (ARG-148).
"""
import argparse
import json
import pathlib
import re
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]


def surf(case, field):
    g = sorted(pathlib.Path(case).glob("postProcessing/surfaces/*/%s_*.raw" % field))
    return np.loadtxt(g[-1]) if g else None


def split_surfaces(x, y, v, nbin=None):
    """Upper and lower, by walking the LOOP, not by any test on a single point.

    THE SAMPLES ARE NOT PAIRED IN x. A C-grid puts surface points around the loop, so upper
    and lower land at DIFFERENT stations: 0.86184 upper, 0.86594 lower, 0.87058 upper, and so
    on. Every per-point test therefore fails --

        sign of y        the aft camber puts the lower surface ABOVE y = 0 past x/c 0.79
        binning in x     no bin holds both surfaces, so the fallback has to guess
        y vs the median  the guess, and it guessed wrong

    -- and each produced the same saw-tooth in Cp, Cf and y+ over the aft chord: the two
    surfaces interleaved, which reads as a numerical instability and is a plotting artefact.

    THE LOOP IS THE ONLY THING THAT KNOWS. Order the points by a nearest-neighbour walk from
    the trailing edge, then cut at the leading edge: everything before the cut is one surface
    and everything after is the other. No assumption about camber, chord line or pairing.
    """
    P = np.column_stack([np.asarray(x), np.asarray(y)])
    val = np.asarray(v)
    n = len(P)
    if n < 4:
        return (P[:, 0], val), (P[:, 0][:0], val[:0])
    # scale y so the walk is not dominated by the chordwise spacing
    S = P / np.array([1.0, max(P[:, 1].ptp(), 1e-9) / max(P[:, 0].ptp(), 1e-9)])
    start = int(np.argmax(P[:, 0]))                 # trailing edge
    order, used = [start], np.zeros(n, bool)
    used[start] = True
    for _ in range(n - 1):
        d = np.hypot(*(S - S[order[-1]]).T)
        d[used] = np.inf
        k = int(np.argmin(d))
        order.append(k)
        used[k] = True
    order = np.array(order)
    Pl, vl_ = P[order], val[order]
    le = int(np.argmin(Pl[:, 0]))                   # leading edge splits the loop
    a_x, a_v = Pl[:le + 1, 0], vl_[:le + 1]
    b_x, b_v = Pl[le:, 0], vl_[le:]
    # the branch with the greater mean height is the upper surface
    if Pl[:le + 1, 1].mean() < Pl[le:, 1].mean():
        (a_x, a_v), (b_x, b_v) = (b_x, b_v), (a_x, a_v)
    o1, o2 = np.argsort(a_x), np.argsort(b_x)
    return (a_x[o1], a_v[o1]), (b_x[o2], b_v[o2])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", required=True)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--out-dir", default="results/argus2d")
    a = ap.parse_args()

    case = REPO / a.case if not pathlib.Path(a.case).is_absolute() else pathlib.Path(a.case)
    tag = a.tag or case.name
    out = REPO / a.out_dir
    out.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(REPO / "scripts"))
    import case_normalisation as cn
    if not (case / cn.NAME).exists():
        (case / cn.NAME).write_text(json.dumps(cn.derive(str(case)), indent=2) + "\n")
    nrm = cn.load(str(case), require_fresh=False)
    f, gg, n = nrm["flow"], nrm["geometry"], nrm["normalisation"]
    chord = gg["chord"]
    print("frame: M %.4f  Re %.3e  U %.2f  rho %.5f  Cp den %.2f  Cf den %.2f (%s)"
          % (f["Mach"], f["Re_chord"], f["U_inf"], f["rho_inf"],
             n["Cp_denominator"], n["Cf_denominator"],
             n["wallShearStress_units"].split()[0]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    made = {}

    # ---------- Cp, Cf, y+ : three stacked surface panels ----------
    p_raw, t_raw, y_raw = surf(case, "p"), surf(case, "wallShearStress"), surf(case, "yPlus")
    if p_raw is not None:
        fig, ax = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
        x, y = p_raw[:, 0] / chord, p_raw[:, 1]
        cp = (p_raw[:, 3] - f["p_inf"]) / n["Cp_denominator"]
        (xu, cu), (xl, cl_) = split_surfaces(x, y, cp)
        ax[0].plot(xu, -cu, "-", lw=1.8, color="#1f6feb", label="upper")
        ax[0].plot(xl, -cl_, "-", lw=1.8, color="#d1660f", label="lower")
        # sonic line: Cp* for the case's own Mach, so "supersonic" is marked not guessed
        M, g_ = f["Mach"], 1.4
        cps = (((2 + (g_ - 1) * M * M) / (g_ + 1)) ** (g_ / (g_ - 1)) - 1) * 2 / (g_ * M * M)
        ax[0].axhline(-cps, color="#888", ls="--", lw=1.1)
        ax[0].annotate("$C_p^*$ (sonic), M=%.3f" % M, xy=(0.02, -cps), fontsize=8, color="#666",
                       va="bottom")
        ax[0].set_ylabel(r"$-C_p$"); ax[0].legend(frameon=False, fontsize=9)
        ax[0].set_title("%s   M %.3f   Re %.2e" % (tag, f["Mach"], f["Re_chord"]),
                        loc="left", fontsize=11)
        made["cp"] = True

        if t_raw is not None:
            tau = np.hypot(t_raw[:, 3], t_raw[:, 4]) * np.sign(t_raw[:, 3])
            cf = tau / n["Cf_denominator"]
            (xu2, fu), (xl2, fl) = split_surfaces(t_raw[:, 0] / chord, t_raw[:, 1], cf)
            ax[1].plot(xu2, np.abs(fu), "-", lw=1.8, color="#1f6feb")
            ax[1].plot(xl2, np.abs(fl), "-", lw=1.8, color="#d1660f")
            # SIGN IS THE DIAGNOSTIC: reversed streamwise shear is separation
            sep = fu > 0
            if sep.any():
                ax[1].plot(xu2[sep], np.abs(fu[sep]), "o", ms=3, color="#c00",
                           label="reversed flow: %.1f%% of upper" % (100 * sep.mean()))
                ax[1].legend(frameon=False, fontsize=9)
            ax[1].set_ylabel(r"$|C_f|$"); ax[1].set_ylim(0, None)
            made["cf"] = True

        if y_raw is not None and y_raw.shape[1] > 3:
            (xu3, yu), _ = split_surfaces(y_raw[:, 0] / chord, y_raw[:, 1], y_raw[:, 3])
            ax[2].plot(xu3, yu, "-", lw=1.8, color="#158463")
            ax[2].axhspan(30, 300, color="#158463", alpha=.10)
            ax[2].annotate("wall-function band 30-300", xy=(0.02, 300), fontsize=8,
                           color="#158463", va="bottom")
            ax[2].set_ylabel(r"$y^+$"); ax[2].set_yscale("log")
            made["yplus"] = True
        ax[2].set_xlabel("x/c")
        for b in ax:
            b.grid(alpha=.25)
            for s in ("top", "right"):
                b.spines[s].set_visible(False)
        fig.tight_layout()
        fig.savefig(out / ("surface_%s.png" % tag), dpi=140)
        plt.close(fig)
        print("wrote %s" % (out / ("surface_%s.png" % tag)).relative_to(REPO))

    # ---------- force history, with the D023 window shaded ----------
    from watch_run import read_dat, d023
    hdr, rows = read_dat(str(case), "forces_coeffs", "forceCoeffs*.dat")
    if rows:
        rows = [[q for q in r if isinstance(q, float)] for r in rows]
        w = max(len(r) for r in rows)
        arr = np.array([r for r in rows if len(r) == w])
        c = {q: i for i, q in enumerate(hdr)} if hdr else {}
        it, cdv, clv = arr[:, 0], arr[:, c.get("Cd", 2)], arr[:, c.get("Cl", 3)]
        fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
        ax[0].plot(it, clv, lw=1.3, color="#1f6feb"); ax[0].set_ylabel(r"$c_l$")
        ax[1].plot(it, cdv, lw=1.3, color="#d1660f"); ax[1].set_ylabel(r"$c_d$")
        w0 = it[max(0, len(it) - 2000)]
        for b in ax:
            b.axvspan(w0, it[-1], color="#888", alpha=.12)
            b.grid(alpha=.25)
            for s in ("top", "right"):
                b.spines[s].set_visible(False)
        ax[0].set_title("%s   D023 %.3f counts over the shaded window   final cl %.5f  cd %.6f"
                        % (tag, d023(cdv.tolist()), clv[-1], cdv[-1]), loc="left", fontsize=10)
        ax[1].set_xlabel("iteration")
        fig.tight_layout(); fig.savefig(out / ("forces_%s.png" % tag), dpi=140); plt.close(fig)
        print("wrote %s   (cl %.5f  cd %.6f  D023 %.3f)"
              % ((out / ("forces_%s.png" % tag)).relative_to(REPO), clv[-1], cdv[-1],
                 d023(cdv.tolist())))
        made["forces"] = True

    # ---------- U and p fields at the final time ----------
    import subprocess
    r = subprocess.run([sys.executable, "scripts/plot_fields.py", "--case", str(case),
                        "--uinf", "%.6f" % f["U_inf"], "--pinf", "%.6f" % f["p_inf"],
                        "--rho", "%.8f" % f["rho_inf"], "--zoom", "1.5",
                        "--out", "%s/fields_%s.png" % (a.out_dir, tag)],
                       cwd=str(REPO), capture_output=True, text=True, timeout=1800)
    if r.returncode == 0:
        print("wrote %s/fields_%s.png" % (a.out_dir, tag))
        made["fields"] = True
    else:
        print("fields FAILED: %s" % (r.stdout + r.stderr)[-200:])

    (out / ("report_%s.json" % tag)).write_text(json.dumps(
        {"case": str(case.relative_to(REPO)), "tag": tag, "figures": sorted(made),
         "normalisation": n, "flow": f, "geometry": gg}, indent=2) + "\n")
    print("figures: %s" % ", ".join(sorted(made)))
    return 0 if made else 1


if __name__ == "__main__":
    raise SystemExit(main())
