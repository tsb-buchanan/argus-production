#!/usr/bin/env python3
"""validate_onera_m6.py: our ONERA M6 run against the AGARD measurements, station by station.

THE 3D VALIDATION ANCHOR. Everything 3D in this project so far is internal consistency: grid
convergence, delta cancellation, sign tests. None of it is anchored to measured data. This is.

WHY M6 AND WHY Cp. It is the canonical 3D swept transonic wing and it carries measured SURFACE
PRESSURE at seven span stations. The hinge moment we owe the DSO team is a pressure integral,
so the validation has to test the pressure field and not just the forces. A drag-only anchor
would leave the deliverable untested.

PER STATION, NEVER POOLED (GEO-087). A single pooled error hides which station failed, and on
this wing the stations are not equivalent: the lambda shock merges into one shock outboard of
y/b 0.8, so the outer stations test something the inner ones do not.

THE COMPARISON IS BINNED THE WAY THE AUTHOR BINNED IT. Cell centres do not land on the
measurement plane, so the surface data is binned in a band of +/- 0.005 b about each station.
That is Alletto's own choice in the tutorial's plot.py and it is kept so that our numbers and
his are formed the same way.

Cp = (p - p_inf) / (0.5 rho_inf U_inf^2). For a COMPRESSIBLE solver OpenFOAM writes p as
ABSOLUTE PRESSURE in Pa, so p_inf must be subtracted; forgetting that is the standing trap
this project records for the 2D compressible work.
"""
import argparse
import pathlib
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path(__file__).resolve().parents[1]
STATIONS = [0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99]
SPAN_M = 1.1963                     # ONERA M6 semispan, metres, published
BLUE, ORANGE, INK, INK2 = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
GRIDC, SURFC = "#e8e7e3", "#fcfcfb"


def experiment():
    """AGARD AR-138 run 308, decoded per ARG-139. Returns {station index: (x/c, Cp)}."""
    txt = (REPO / "validation_data/onera_m6/case_2308.dat").read_text()
    out = {}
    for i, z in enumerate(re.split(r"^ZONE ", txt, flags=re.M)[1:]):
        rows = [[float(v) for v in ln.split()] for ln in z.splitlines()[1:]
                if len(ln.split()) == 5]
        A = np.array(rows)
        out[i] = (A[:, 2], A[:, 4], A[:, 3])       # x/L, Cp, z/L (z sign gives the surface)
    return out


def shock_mask(xq, x, cp, halfwidth=0.05, thresh=5.0):
    """Which query points sit in ANY compression band, not just 'the' shock.

    THE M6 HAS A LAMBDA SHOCK: two shocks inboard of y/b 0.8 which merge outboard, as the
    reference page itself states. Locating a single shock and excluding one band around it
    therefore leaves the other shock inside the 'off-shock' error, which is what a
    single-shock model quietly did here. This instead flags every location where the CFD
    Cp gradient exceeds `thresh` and excludes a band about all of them, so it does not care
    how many shocks there are.
    """
    g = np.abs(np.gradient(cp, x))
    hot = x[g > thresh]
    if not len(hot):
        return np.zeros(len(xq), bool)
    return (np.abs(xq[:, None] - hot[None, :]) <= halfwidth).any(axis=1)


def shock_x(x, cp):
    """Upper-surface shock: the steepest ADVERSE dCp/dx over the mid-chord.

    Restricted to 0.05 < x/c < 0.85. The lower bound keeps the leading-edge suction peak,
    which is steeper than any shock, from winning. The UPPER bound keeps the trailing-edge
    pressure recovery out: at 0.95 the detector returned "shocks" at x/c 0.92, 0.94 and 0.87
    on three stations, which is the TE recovery, and because the shock location sets which
    points are excluded from the off-shock RMS, a misplaced shock corrupts the error metric
    as well as the reported position. No shock sits aft of 0.85 on this wing at this
    condition.
    """
    m = (x > 0.05) & (x < 0.85)
    if m.sum() < 10:
        return None
    xs, cs = x[m], cp[m]
    g = np.gradient(cs, xs)
    i = int(np.argmax(g))            # Cp RISES through a shock, so the max positive gradient
    return float(xs[i]) if g[i] > 2.0 else None


# A LAMBDA LEG IS 0.25 TO 0.30 CHORD FROM ITS PARTNER ON THIS WING; a position error is
# a few hundredths. So 0.10 separates "the two detectors found the same shock" from "they
# found different legs" with a factor of 2.5 clearance on either side, and it is not a
# tuned constant: no value between 0.08 and 0.20 changes which stations pair.
LAMBDA_PAIR_XC = 0.10


def pair_shocks(xsh_rans, xsh_agard):
    """The offset, or a refusal, and NEVER a number formed from two different shocks.

    A SINGLE-SHOCK DETECTOR APPLIED TO BOTH SIDES CAN SELECT DIFFERENT LEGS. Inboard of
    y/b 0.8 the M6 carries a lambda shock, two compressions, and the tap data resolves it
    with 23 points where the CFD has hundreds. When each detector picks the stronger of
    the two AT ITS OWN RESOLUTION, subtracting the positions measures the LEG SPACING and
    reports it as a position error. Measured here: y/b 0.20 gives 0.350 against 0.060 and
    y/b 0.65 gives 0.403 against 0.150, both about a quarter chord, which is the lambda
    structure and not a quarter-chord error in our shock.

    So the offset is quoted only where the two are close enough to be the same shock, and
    withheld with a reason otherwise. Withholding a number is a result; quoting the wrong
    one is not (D060 in spirit: the failure has a known signature, so test for it).
    """
    if xsh_rans is None or xsh_agard is None:
        return None, "-"
    d = xsh_rans - xsh_agard
    if abs(d) > LAMBDA_PAIR_XC:
        return None, "AMBIG"
    return float(d), ""


def raw(path):
    """OpenFOAM `raw` surface sample: x y z value."""
    A = np.loadtxt(path, comments="#")
    return A[:, :3], A[:, 3]


def _dictval(txt, key):
    """A scalar entry, following ONE level of OpenFOAM `$variable` indirection."""
    m = re.search(r"^\s*%s\s+([^;]+);" % re.escape(key), strip_comments(txt), re.M)
    if not m:
        return None
    v = m.group(1).strip().split()[-1]
    if v.startswith("$"):
        m2 = re.search(r"^\s*%s\s+([-\d.eE+]+)\s*;" % re.escape(v[1:]),
                       strip_comments(txt), re.M)
        return float(m2.group(1)) if m2 else None
    try:
        return float(v)
    except ValueError:
        return None


def strip_comments(t):
    t = re.sub(r"/\*.*?\*/", "", t, flags=re.S)
    return "\n".join(l.split("//")[0] for l in t.splitlines())


def freestream_from_case(case):
    """p_inf, T_inf, R and U READ OUT OF THE CASE. Never defaulted, or the frame is a guess.

    THIS IS D068 IN OUR OWN INSTRUMENT AND IT COST A REAL BIAS. The function used to take
    pinf = 101325 Pa and T = 293 K as DEFAULTS while this case's own `0.orig/p` sets
    pOut 1e5 and its `0.orig/T` sets Tinlet 298. Wrong on both, so q came out 50689 Pa
    against the case's 49187, 3.06% high, and Cp carried a uniform offset of +0.027 on top.
    A coefficient normalised on a number supplied at the command line is exactly the
    unauditable frame this project keeps being bitten by.

    IT WAS CAUGHT BY A KNOWN ANSWER, not by inspection. Isentropic stagnation Cp at this
    Mach is 1.1882, computed from gas dynamics and not from anything here. The maximum
    sampled wall Cp reads 1.1035 on the defaulted frame, 7.1% low, and 1.1641 on the
    case's own, 2.0% low. THE SIGN OF THE RESIDUAL IS FIXED BY PHYSICS: a face centre on
    a coarse mesh cannot sit exactly at the stagnation point, so the sampled maximum can
    only UNDER-read the true peak, never exceed it (D060). Both frames under-read; only
    one under-reads by an amount a 1 M-cell mesh can explain.
    """
    case = pathlib.Path(case)
    src = {}
    p = t = None
    for d in ("0.orig/p", "0/p"):
        if (case / d).exists():
            p = _dictval((case / d).read_text(), "internalField")
            if p is None:
                p = _dictval((case / d).read_text(), "pOut")
            if p is not None:
                src["p_inf"] = d
                break
    for d in ("0.orig/T", "0/T"):
        if (case / d).exists():
            t = _dictval((case / d).read_text(), "internalField")
            if t is None:
                t = _dictval((case / d).read_text(), "Tinlet")
            if t is not None:
                src["T_inf"] = d
                break
    W = 28.9
    tp = case / "constant/thermophysicalProperties"
    if tp.exists():
        w = _dictval(tp.read_text(), "molWeight")
        if w:
            W, src["molWeight"] = w, "constant/thermophysicalProperties"
    U = None
    for d in ("0.orig/U", "0/U"):
        f = case / d
        if f.exists():
            m = re.search(r"(?:Uinlet|internalField\s+uniform)\s*\(\s*([-\d.eE+]+)\s+"
                          r"([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", strip_comments(f.read_text()))
            if m:
                U = float(np.linalg.norm([float(x) for x in m.groups()]))
                src["U_inf"] = d
                break
    return dict(p_inf=p, T_inf=t, molWeight=W, U=U, R=8314.47 / W, sources=src)


def load_surface(case, time, pinf=None, tinf=None, uinf=None, rhoinf=None,
                 verbose=True):
    """Our wall sample, reduced to (points, {quantity: value}, freestream dict).

    EXTRACTED FROM main() SO A SECOND ENTRY POINT CANNOT DRIFT FROM IT. The TMR overlay
    script needs exactly this state, and a copy of these twenty lines would have two
    freestream definitions in the repo the moment either was touched. main() calls this,
    so anything that changes here changes both figures or neither.

    pinf/tinf/uinf/rhoinf are OVERRIDES, not defaults. Unset, every one comes from the
    case; see freestream_from_case for why that distinction is not cosmetic.
    """
    F = freestream_from_case(case)
    # HIS plot.py USES sqrt(290^2 + 13.73^2), i.e. alpha 2.712 deg, but the case's own 0/U
    # sets Uinlet (290 0 15.5), i.e. alpha 3.060 deg, which is the measured condition. The
    # case file wins; the plot script is stale. It matters not at all for q (0.03%) and a
    # great deal for what angle we claim to have run, so it is recorded rather than smoothed.
    U = uinf if uinf else F["U"]
    pinf = pinf if pinf else F["p_inf"]
    tinf = tinf if tinf else F["T_inf"]
    if U is None or pinf is None or tinf is None:
        raise SystemExit(
            "cannot read the freestream from %s: U %s, p_inf %s, T_inf %s. Cp and Cf are "
            "normalised on these, so a guess here is an unauditable frame (D068). Pass "
            "--pinf/--tinf/--uinf explicitly if the case genuinely does not carry them."
            % (case, U, pinf, tinf))
    rho = rhoinf if rhoinf else pinf / (F["R"] * tinf)
    q = 0.5 * rho * U * U
    a = float(np.sqrt(1.4 * F["R"] * tinf))
    mach = U / a
    cp_stag = float(((1 + 0.2 * mach ** 2) ** 3.5 - 1) / (0.7 * mach ** 2))
    if verbose:
        print("freestream FROM THE CASE: U %.3f m/s   p_inf %.1f Pa   T_inf %.1f K"
              % (U, pinf, tinf))
        print("  -> rho %.5f   q %.1f Pa   a %.2f   M %.4f   (%s)"
              % (rho, q, a, mach, ", ".join("%s from %s" % kv for kv in F["sources"].items())))

    S = pathlib.Path(case) / "postProcessing/samplePwall/surface" / str(time)
    P, p = raw(S / "p_patch_wing.raw")
    if verbose:
        print("wall sample: %d faces   p range %.1f .. %.1f Pa" % (len(p), p.min(), p.max()))
    QTY = {"Cp": (p - pinf) / q}
    # KNOWN-ANSWER GATE ON THE FRAME ITSELF, and it is a SIGN test rather than a tolerance
    # (D060). Isentropic stagnation Cp comes from gas dynamics, not from this code. A face
    # centre on a finite mesh cannot sit exactly at the stagnation point, so the sampled
    # maximum can only UNDER-read it. EXCEEDING IT IS THEREFORE IMPOSSIBLE and means the
    # normalisation is wrong, which is precisely how the 101325/293 default was caught.
    cp_max = float(QTY["Cp"].max())
    if verbose:
        print("  Cp_max %.4f against the isentropic stagnation value %.4f (%.1f%% low)"
              % (cp_max, cp_stag, 100.0 * (cp_stag - cp_max) / cp_stag))
    if cp_max > cp_stag * 1.02:
        raise SystemExit(
            "Cp_max %.4f EXCEEDS the isentropic stagnation Cp %.4f. A sampled face centre "
            "cannot out-read the stagnation point, so this is a normalisation error, not "
            "a flow feature. Check p_inf, T_inf and |U| against the case." % (cp_max, cp_stag))
    # Cf and y+ on the same faces. wallShearStress is a VECTOR in the raw file, so its
    # magnitude is taken; Cf = |tau_w| / q, which is his definition in plot.py.
    try:
        W = np.loadtxt(S / "wallShearStress_patch_wing.raw", comments="#")
        QTY["Cf"] = np.linalg.norm(W[:, 3:6], axis=1) / q
    except Exception as exc:
        print("  wallShearStress unavailable: %s" % exc)
    try:
        QTY["yPlus"] = np.loadtxt(S / "yPlus_patch_wing.raw", comments="#")[:, 3]
    except Exception as exc:
        print("  yPlus unavailable: %s" % exc)
    if verbose:
        print("  quantities: %s" % ", ".join(QTY))
    return P, QTY, {"U": U, "rho": rho, "pinf": pinf, "tinf": tinf, "q": q,
                    "a": a, "mach": mach, "cp_max": cp_max,
                    "cp_stagnation_isentropic": cp_stag,
                    "R": F["R"], "molWeight": F["molWeight"],
                    "sources": F["sources"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", default="cases/validation/onera_m6_A_alletto")
    ap.add_argument("--time", default="500")
    ap.add_argument("--pinf", type=float, default=None,
                    help="OVERRIDE; default is the case's own 0.orig/p")
    ap.add_argument("--rhoinf", type=float, default=None,
                    help="default: perfect gas from pinf and T")
    ap.add_argument("--tinf", type=float, default=None,
                    help="OVERRIDE; default is the case's own 0.orig/T")
    ap.add_argument("--uinf", type=float, default=None,
                    help="OVERRIDE; default is the case's own 0.orig/U")
    # 0.0025 IS HIS NUMBER, from plot.py: binwdith = 0.0025*b, filtered on |y - yPos|.
    # An earlier run used 0.005 and therefore averaged over twice his span band.
    ap.add_argument("--band", type=float, default=0.0025, help="+/- band in y/b, as the author")
    ap.add_argument("--out", default="results/figures/deck/onera_m6_validation.png")
    ap.add_argument("--json", default="results/onera_m6_validation.json",
                    help="machine-readable per-station summary; '' to skip")
    a = ap.parse_args()

    case = REPO / a.case
    P, QTY, _fs = load_surface(case, a.time, a.pinf, a.tinf, a.uinf, a.rhoinf)
    x, y, z = P[:, 0], P[:, 1], P[:, 2]
    eta = y / SPAN_M

    exp = experiment()
    # 3 ROWS (Cp, Cf, yPlus) x 7 STATIONS, which is his plot.py's full output laid out as one
    # figure so it can be set beside his published page panel for panel.
    order = [k for k in ("Cp", "Cf", "yPlus") if k in QTY]
    fig, ax = plt.subplots(len(order), len(STATIONS),
                           figsize=(3.05 * len(STATIONS), 3.15 * len(order)),
                           facecolor=SURFC, squeeze=False)
    rows = []
    for r, qty in enumerate(order):
        val = QTY[qty]
        for k, st in enumerate(STATIONS):
            A = ax[r][k]
            A.set_facecolor(SURFC)
            A.grid(True, color=GRIDC, lw=0.8)
            A.set_axisbelow(True)
            for sp in ("top", "right"):
                A.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                A.spines[sp].set_color(GRIDC)
            A.tick_params(colors=INK2, labelsize=8, length=0)
            m = np.abs(eta - st) <= a.band
            if m.sum() < 20:
                A.text(0.5, 0.5, "%d faces" % m.sum(), transform=A.transAxes, ha="center",
                       color="#c0392b", fontsize=9)
                continue
            xb, zb, vb = x[m], z[m], val[m]
            xc = (xb - xb.min()) / (xb.max() - xb.min())
            ex, ecp, ez = exp[k]
            for lab, selc, sele, col in (("upper", zb > 0, ez > 0, BLUE),
                                         ("lower", zb < 0, ez < 0, ORANGE)):
                xs, cs = xc[selc], vb[selc]
                o = np.argsort(xs)
                xs, cs = xs[o], cs[o]
                A.plot(xs, cs, "-", lw=1.5, color=col, label="RANS %s" % lab)
                if qty != "Cp":
                    continue
                exs, ecs = ex[sele], ecp[sele]
                A.plot(exs, ecs, "o", ms=4, mfc="none", mew=1.2, color=col,
                       label="AGARD %s" % lab)
                good = (exs >= xs.min()) & (exs <= xs.max())
                if not good.sum():
                    continue
                d = np.interp(exs[good], xs, cs) - ecs[good]
                # PLAIN RMS OVER EVERY TAP, NO EXCLUSION, AND THAT IS DELIBERATE.
                # Three attempts were made to quote an "off-shock" error and each acquired
                # its own failure mode: a single-shock model missed the second leg of the
                # lambda shock; widening the search caught the trailing-edge recovery and
                # called it a shock; a gradient-threshold mask excluded EVERY tap at three
                # stations because a 1 M-cell mesh has noisy gradients. Each fix needed a
                # tuned constant, and a metric that needs tuning to look right is not
                # measuring the thing any more. RMS over all taps needs no constant, cannot
                # misfire, and the shock's contribution is visible in the max beside it.
                xsh = shock_x(xs, cs) if lab == "upper" else None
                # THE SAME DETECTOR ON THE TAPS, so the position error is formed by one
                # operator applied to both sides rather than by comparing our detector
                # against a number read off a plot. ITS RESOLUTION IS BOUNDED BY THE TAP
                # SPACING and that bound is reported beside it: a shock located between
                # two taps cannot be placed more finely than the gap between them, so a
                # position error smaller than the local spacing is not a measurement.
                xsh_e, tapgap = None, None
                if lab == "upper" and good.sum() >= 10:
                    oe = np.argsort(exs[good])
                    xe, ce = exs[good][oe], ecs[good][oe]
                    xsh_e = shock_x(xe, ce)
                    if xsh_e is not None:
                        near = np.abs(xe - xsh_e) < 0.15
                        tapgap = float(np.median(np.diff(xe[near]))) if near.sum() > 2 else None
                rows.append((st, lab, float(np.sqrt(np.mean(d ** 2))),
                             float(np.abs(d).max()), xsh, xsh_e, tapgap, int(good.sum())))
            if qty == "Cp":
                A.invert_yaxis()
            if qty == "yPlus":
                A.set_yscale("log")
            if r == 0:
                A.set_title("y/b = %.2f" % st, fontsize=11, color=INK, fontweight="bold")
            if r == len(order) - 1:
                A.set_xlabel("x/c", color=INK2, fontsize=9)
            if k == 0:
                A.set_ylabel({"Cp": "$C_p$", "Cf": "$C_f$", "yPlus": "$y^+$"}[qty],
                             color=INK2, fontsize=12)
            if r == 0 and k == 0:
                A.legend(fontsize=7.5, frameon=False, labelcolor=INK2)

    fig.suptitle("ONERA M6 as a new-geometry test: AGARD AR-138 run 308, M 0.8395, "
                 "alpha 3.06 deg, Re 11.72e6", fontsize=15, color=INK, fontweight="bold",
                 x=0.006, ha="left")
    fig.text(0.006, 0.012, "Rows: surface pressure, skin friction, wall resolution. Binned "
             "+/- %.4f in y/b about each station, matching the reference tutorial's own "
             "plot.py so the panels are directly comparable to its published figures. Only "
             "Cp has measured data; Cf and y+ are ours against his figures." % a.band,
             fontsize=9.5, color=INK2)
    fig.tight_layout(rect=[0, 0.03, 1, 0.955])
    o = REPO / a.out
    o.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(o, dpi=130, facecolor=SURFC)
    plt.close(fig)
    print("  wrote %s" % o.relative_to(REPO))

    print("\nPER-STATION agreement against the taps (never pooled).")
    print("  RMS is over EVERY tap, with no shock exclusion (see the note in the source).")
    print("  max is dominated by the shock, where a small position offset gives a large dCp.")
    print("  %-6s %-7s %-10s %-10s %-11s %-11s %-9s %s"
          % ("y/b", "surf", "RMS all", "max", "shk RANS", "shk AGARD", "d x/c", "taps"))
    allrms = []
    for st, lab, rms, mx, xsh, xsh_e, gap, nall in rows:
        allrms.append(rms)
        dsh, why = pair_shocks(xsh, xsh_e)
        print("  %-6.2f %-7s %-10.4f %-10.4f %-11s %-11s %-9s %d"
              % (st, lab, rms, mx,
                 ("%.3f" % xsh) if xsh else "-", ("%.3f" % xsh_e) if xsh_e else "-",
                 ("%+.3f" % dsh) if dsh is not None else why, nall))
    if allrms:
        up = [r for (st, lab, r, *_) in rows if lab == "upper"]
        lo = [r for (st, lab, r, *_) in rows if lab == "lower"]
        print("\n  RMS in Cp over all taps: upper %.4f to %.4f, lower %.4f to %.4f"
              % (min(up), max(up), min(lo), max(lo)))
        print("  %d of %d station-surfaces compared" % (len(allrms), 2 * len(STATIONS)))
        print("  A shock position error below the local tap spacing is NOT a measurement;")
        print("  the spacing is carried per station in the JSON as agard_tap_spacing_xc.")

    # MACHINE-READABLE, because record-keeping rule 3 asks for it and because a number
    # that exists only in a terminal cannot be checked by anything. Written LAST, from the
    # same `rows` the table printed, so the file and the table cannot disagree.
    if a.json:
        import json as _json
        stations = []
        for st, lab, rms, mx, xsh, xsh_e, gap, nall in rows:
            d, why = pair_shocks(xsh, xsh_e)
            stations.append(dict(
                eta=st, surface=lab, cp_rms=rms, cp_max_abs=mx,
                shock_xc_rans=xsh, shock_xc_agard=xsh_e, shock_dxc=d,
                shock_pairing=("paired" if d is not None else
                               ("no shock detected" if why == "-" else
                                "AMBIGUOUS: the detectors selected different lambda legs")),
                agard_tap_spacing_xc=gap, n_taps=nall))
        paired = [s for s in stations if s["shock_dxc"] is not None]
        rec = dict(
            _what="ONERA M6 surface pressure against AGARD AR-138 run 308, per station.",
            _never_pooled=("A single pooled error hides which station failed, and the "
                           "stations are not equivalent: the lambda shock merges into one "
                           "shock outboard of y/b 0.8 (GEO-087)."),
            _rms_definition=("RMS of (RANS - tap) in Cp over EVERY tap in the station band, "
                             "with NO shock exclusion. Three exclusion schemes were tried "
                             "and each needed a tuned constant; the shock's contribution is "
                             "visible in cp_max_abs beside it."),
            _shock_caveat=("shock_xc_rans and shock_xc_agard come from ONE detector applied "
                           "to both sides. The AGARD value cannot be finer than the local "
                           "tap spacing, which is reported per station. shock_dxc is null "
                           "where the two detectors selected DIFFERENT LEGS of the lambda "
                           "shock: subtracting those would report the leg spacing, about a "
                           "quarter chord, as a position error. See pair_shocks()."),
            shock_position_paired=dict(
                n_paired=len(paired), n_upper_stations=len(STATIONS),
                max_abs_dxc=max((abs(s["shock_dxc"]) for s in paired), default=None),
                stations=[s["eta"] for s in paired]),
            case=a.case, time=a.time, band_eta=a.band,
            reference="AGARD AR-138 (Cook, McDonald & Firmin 1979) run 308, decoded ARG-139",
            condition=dict(mach=0.8395, alpha_deg=3.06, Re=11.72e6),
            freestream=_fs, semispan_m=SPAN_M,
            n_compared=len(rows), n_possible=2 * len(STATIONS),
            cp_rms_upper_range=[min(up), max(up)] if allrms else None,
            cp_rms_lower_range=[min(lo), max(lo)] if allrms else None,
            stations=stations,
            figure=str(pathlib.Path(a.out)),
            generated_by="scripts/validate_onera_m6.py",
        )
        jp = REPO / a.json
        jp.parent.mkdir(parents=True, exist_ok=True)
        jp.write_text(_json.dumps(rec, indent=2) + "\n")
        print("  wrote %s" % a.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
