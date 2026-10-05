#!/usr/bin/env python3
"""Hard gates on an exported wing OML, measured FROM THE SURFACE.

Runs on the HOST (numpy only). Everything here is measured from the triangle
soup, never read back from the .vsp3 parms: the question these gates answer is
whether the EXPORT carries the placement transform, so reading the design
values would beg the question.

Gates (docs/decisions.md D048):
  1. DIHEDRAL PRESENT  tip quarter-chord ~0.139 m above root quarter-chord.
  2. INCIDENCE PRESENT root section chord line ~1.873 deg nose-up.
  3. CHORDS            root 0.640081 m, tip 0.152400 m, measured in the
                       section plane and so dihedral-independent.
  4. SWEEP AND TWIST   LE sweep 28.86199 deg; twist 0 at root to -4.11495 deg
                       at tip, taken as the tip-minus-root difference.

Reported but NOT gated: projected y-extent and planform area, because VSP's
Bref convention for projected versus true span is ambiguous here.

NOTE ON WHY SPAN AND AREA ARE NOT GATES: the superseded in-house loft measures
3.657598 m span and 1.240926 m2 area, the exact reference values, while
missing the dihedral entirely. Those quantities pass on the defect.
"""
import argparse
import json
import math
import struct
import sys

import numpy as np

from section_conventions import le_te_base

TARGETS = {
    "dihedral_rise_m": 0.139,
    "root_incidence_deg": 1.873,
    "root_chord_m": 0.640081,
    "tip_chord_m": 0.152400,
    "le_sweep_deg": 28.86199,
    "tip_twist_deg": -4.11495,
}


def read_stl_vertices(path):
    with open(path, "rb") as fh:
        if fh.read(5) == b"solid":
            fh.seek(0)
            verts = []
            for line in fh:
                s = line.split()
                if s and s[0] == b"vertex":
                    verts.append((float(s[1]), float(s[2]), float(s[3])))
            return np.asarray(verts, dtype=float)
        fh.seek(80)
        n = struct.unpack("<I", fh.read(4))[0]
        buf = fh.read(n * 50)
        arr = np.frombuffer(
            np.frombuffer(buf, dtype=np.uint8).reshape(n, 50)[:, 12:48].tobytes(),
            dtype="<f4").reshape(n * 3, 3)
        return arr.astype(float)


def section_at(verts, axis, s_target, tol):
    """Vertices whose station along `axis` is within tol of s_target."""
    s = verts @ axis
    return verts[np.abs(s - s_target) < tol]


def le_te(sec, te_ref="midpoint", band_frac=0.005):
    """DELEGATES to scripts/section_conventions.py (GEO-084).

    The rule used to live here as a flag, and it did not travel: it was re-derived
    wrongly in check_synthetic_gates.py, which read 0.0602 deg of twist on a wing
    built with exactly zero. A convention that lives as a flag on one caller does
    not travel; one that lives in a shared module cannot fail to. The docstring
    that explained the trap stayed with the rule, in section_conventions.py.
    """
    return le_te_base(sec, te_ref=te_ref, band_frac=band_frac)


def analyse(path, half_tol=0.004, te_ref="midpoint"):
    v = read_stl_vertices(path)

    # Work on the right half only (y >= 0) when the export is symmetric, so the
    # tip is unambiguous.
    if v[:, 1].min() < -1e-6 and v[:, 1].max() > 1e-6:
        symmetric = True
        v = v[v[:, 1] >= -1e-9]
    else:
        symmetric = False

    # Spanwise axis: the wing is the unswept-in-y loft rotated about x by the
    # dihedral, so the true spanwise direction lies in the y-z plane. Recover it
    # from the vertex cloud rather than assuming it, since the dihedral is the
    # very thing under test. Take the direction in the y-z plane that maximises
    # extent, via the principal axis of the (y,z) cloud weighted to the outboard
    # half where the dihedral offset dominates the airfoil thickness.
    yz = v[:, 1:3]
    outer = yz[yz[:, 0] > 0.5 * yz[:, 0].max()]
    c = outer - outer.mean(axis=0)
    w, vecs = np.linalg.eigh(c.T @ c)
    d = vecs[:, np.argmax(w)]
    if d[0] < 0:
        d = -d
    axis = np.array([0.0, d[0], d[1]])
    axis /= np.linalg.norm(axis)

    s = v @ axis
    s_root, s_tip = s.min(), s.max()

    root = section_at(v, axis, s_root, half_tol)
    tip = section_at(v, axis, s_tip, half_tol)

    le_r, te_r, base_r = le_te(root, te_ref)
    le_t, te_t, base_t = le_te(tip, te_ref)

    qc_r = le_r + 0.25 * (te_r - le_r)
    qc_t = le_t + 0.25 * (te_t - le_t)

    # Gate 1: quarter-chord rise.
    rise = qc_t[2] - qc_r[2]

    # Gate 3: chord measured as the 3D distance in the section plane.
    chord_r = float(np.linalg.norm(te_r - le_r))
    chord_t = float(np.linalg.norm(te_t - le_t))

    # Gate 2 and 4: section incidence measured in the LOCAL section plane, so
    # the dihedral rotation does not contaminate the angle. e1 is the global x
    # axis made perpendicular to the spanwise axis; e3 completes the frame.
    e1 = np.array([1.0, 0.0, 0.0])
    e1 = e1 - (e1 @ axis) * axis
    e1 /= np.linalg.norm(e1)
    e3 = np.cross(axis, e1)
    if e3[2] < 0:
        e3 = -e3

    def incidence(le, te):
        ch = te - le
        return math.degrees(math.atan2(-(ch @ e3), ch @ e1))

    inc_r = incidence(le_r, te_r)
    inc_t = incidence(le_t, te_t)

    # Naive x-z angle, reported alongside to show the frame choice matters.
    def incidence_xz(le, te):
        ch = te - le
        return math.degrees(math.atan2(-(ch[2]), ch[0]))

    # Gate 4: LE sweep from root LE to tip LE, spanwise distance taken in the
    # plane normal to x.
    dle = le_t - le_r
    sweep = math.degrees(math.atan2(dle[0], math.hypot(dle[1], dle[2])))

    # ---- Rigorous recovery of the parent transform, surface-measured -------
    # The rise threshold alone is a weak gate. The root-to-tip quarter-chord
    # offset plus the independently measured root incidence over-determines the
    # parent rotation, so solve it exactly and check the recovered dihedral and
    # the recovered WING-FRAME half-span against 5.000 deg and 6.0000 ft.
    #
    # Wing frame: qc line flat (no dihedral), offset (dx0, dy0, 0).
    # Apply Y rotation (incidence) then X rotation (dihedral):
    #   dx'' = dx0 cos(ty)
    #   dy'' = dy0 cos(tx) + A sin(tx)      with A = dx0 sin(ty)
    #   dz'' = dy0 sin(tx) - A cos(tx)
    # Hence  dy'' sin(tx) - dz'' cos(tx) = A, solved in closed form below.
    dqc = qc_t - qc_r
    ty = math.radians(inc_r)
    dx0 = dqc[0] / math.cos(ty)
    A = dx0 * math.sin(ty)
    R = math.hypot(dqc[1], dqc[2])
    phi = math.atan2(dqc[2], dqc[1])
    # Raise rather than return an empty dict: a silently absent recovery block
    # reads as "not applicable" when it actually means the measurement failed.
    if not (R > 0 and abs(A / R) <= 1.0):
        raise RuntimeError(
            "transform recovery is not solvable for %s: R=%.6g A=%.6g. "
            "This is a measurement failure, not a geometry verdict." % (path, R, A))
    tx = phi + math.asin(A / R)
    dy0 = math.sqrt(max(R * R - A * A, 0.0))
    recovered = {
        "dihedral_deg": math.degrees(tx),
        "wing_frame_half_span_m": dy0,
        "wing_frame_half_span_ft": dy0 / 0.3048,
        "wing_frame_qc_x_offset_m": dx0,
        "residual_rise_check_m": (dy0 * math.sin(tx) - A * math.cos(tx)) - dqc[2],
    }

    # Reported, not gated.
    y_extent_half = float(v[:, 1].max() - v[:, 1].min())
    span_true_half = float(s_tip - s_root)

    res = {
        "file": path,
        "vertices": int(len(v)),
        "export_is_symmetric_full_span": symmetric,
        "spanwise_axis": axis.tolist(),
        "pca_axis_dihedral_deg_CRUDE": math.degrees(math.atan2(axis[2], axis[1])),
        "recovered_transform": recovered,
        "gates": {
            "1_dihedral_rise_m": rise,
            "2_root_incidence_deg": inc_r,
            "3_root_chord_m": chord_r,
            "3_tip_chord_m": chord_t,
            "4_le_sweep_deg": sweep,
            "4_tip_minus_root_twist_deg": inc_t - inc_r,
        },
        "te_reference": te_ref,
        "reported_not_gated": {
            "root_blunt_base_thickness_m": base_r,
            "root_blunt_base_pct_chord": 100.0 * base_r / chord_r,
            "tip_blunt_base_thickness_m": base_t,
            "tip_blunt_base_pct_chord": 100.0 * base_t / chord_t,
            "tip_incidence_deg": inc_t,
            "root_incidence_deg_naive_xz": incidence_xz(le_r, te_r),
            "tip_incidence_deg_naive_xz": incidence_xz(le_t, te_t),
            "projected_y_extent_half_m": y_extent_half,
            "projected_y_extent_full_m": 2.0 * y_extent_half,
            "true_span_along_axis_half_m": span_true_half,
            "true_span_along_axis_full_m": 2.0 * span_true_half,
            "root_qc_xyz": qc_r.tolist(),
            "tip_qc_xyz": qc_t.tolist(),
        },
    }
    return res


def verdict(res, tol):
    g = res["gates"]
    rows = [
        ("1 DIHEDRAL PRESENT", g["1_dihedral_rise_m"], TARGETS["dihedral_rise_m"], tol["rise"], "m"),
        ("2 INCIDENCE PRESENT", g["2_root_incidence_deg"], TARGETS["root_incidence_deg"], tol["ang"], "deg"),
        ("3 ROOT CHORD", g["3_root_chord_m"], TARGETS["root_chord_m"], tol["chord"], "m"),
        ("3 TIP CHORD", g["3_tip_chord_m"], TARGETS["tip_chord_m"], tol["chord"], "m"),
        ("4 LE SWEEP", g["4_le_sweep_deg"], TARGETS["le_sweep_deg"], tol["ang"], "deg"),
        ("4 TIP TWIST", g["4_tip_minus_root_twist_deg"], TARGETS["tip_twist_deg"], tol["ang"], "deg"),
    ]
    out, allpass = [], True
    for name, got, want, t, unit in rows:
        ok = abs(got - want) <= t
        allpass &= ok
        out.append((name, got, want, got - want, t, unit, "PASS" if ok else "FAIL"))
    return out, allpass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stl")
    ap.add_argument("--json-out")
    ap.add_argument("--tol-rise", type=float, default=0.010)
    ap.add_argument("--tol-ang", type=float, default=0.25)
    ap.add_argument("--tol-chord", type=float, default=0.002)
    ap.add_argument("--te-ref", default="midpoint", choices=("midpoint", "maxx"),
                    help="chord TE reference; 'midpoint' is blunt/sharp insensitive")
    args = ap.parse_args()

    res = analyse(args.stl, te_ref=args.te_ref)
    tol = {"rise": args.tol_rise, "ang": args.tol_ang, "chord": args.tol_chord}
    rows, allpass = verdict(res, tol)

    print("TE reference: %s" % res["te_reference"])
    print("%-22s %12s %12s %12s %8s  %s" % ("GATE", "MEASURED", "TARGET", "DELTA", "TOL", ""))
    for name, got, want, d, t, unit, status in rows:
        print("%-22s %12.6f %12.6f %12.6f %8.3f  %s  (%s)" % (name, got, want, d, t, status, unit))
    print()
    print("OVERALL:", "PASS" if allpass else "FAIL")
    globals()["_ALLPASS"] = allpass
    print()
    print("reported, not gated:")
    for k, val in res["reported_not_gated"].items():
        print("   %-32s %s" % (k, val))
    print()
    print("recovered parent transform (surface-measured, closed form):")
    for k, val in res["recovered_transform"].items():
        print("   %-32s %.6f" % (k, val))
    print("   %-32s %.4f deg  (crude estimator, diagnostic only)"
          % ("pca_axis_dihedral", res["pca_axis_dihedral_deg_CRUDE"]))

    res["verdict"] = {"all_gates_pass": bool(allpass),
                      "rows": [{"gate": r[0], "measured": r[1], "target": r[2],
                                "delta": r[3], "tol": r[4], "unit": r[5],
                                "status": r[6]} for r in rows]}
    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(res, fh, indent=2)


if __name__ == "__main__":
    # GEO-082: this printed OVERALL: FAIL and STILL EXITED 0, on every input, so
    # nothing but a human reading stdout could ever act on the verdict. A gate that
    # cannot fail is not a gate. Verified against a deliberately non-conforming
    # surface before and after. No caller depended on the old status: the only
    # programmatic caller in the repo is check_synthetic_gates.py --demo-argus-suite.
    main()
    sys.exit(0 if globals().get("_ALLPASS") else 1)
