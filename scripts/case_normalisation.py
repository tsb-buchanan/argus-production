#!/usr/bin/env python3
"""case_normalisation.py: ONE source of truth for every quantity that normalises a QoI.

WHY (requirement, 2026-08-08). One case-settings JSON holding the values Cf and Cp are
normalised with, so they cannot be mistyped, and updated whenever the case changes.

The motivating failure: our skin friction read 13% high for a day because Cf was divided by
0.5*U^2 instead of 0.5*rho*U^2. OpenFOAM writes wallShearStress rho-weighted (Pa) for a
COMPRESSIBLE solver and kinematic (m2/s2) for an incompressible one, and the plotting script
carried a docstring asserting the kinematic form unconditionally. Every downstream figure
inherited it. The error was not in the solver, the mesh or the forces; it was in ONE
NORMALISATION APPLIED IN ONE PLACE, and there was no single place to check it.

WHAT THIS FIXES, structurally:
  1. Every normalising quantity is DERIVED FROM THE CASE and written to
     <case>/qoi_normalisation.json, each with the file it came from. Nothing is typed in.
  2. THE Cp AND Cf DENOMINATORS ARE STORED EXPLICITLY, as numbers, beside the convention that
     produced them. A consumer divides by a stored denominator; it does not reconstruct one.
  3. THE wallShearStress CONVENTION IS DERIVED FROM THE SOLVER NAME, not assumed, with the
     OpenFOAM source line that establishes it recorded in the file.
  4. THE SOURCE FILES ARE HASHED. `--verify` re-derives from the case and fails if anything
     has moved, so "make sure these get updated if it changes" is CHECKED rather than
     remembered. A stale normalisation is exactly as dangerous as a wrong one and much harder
     to see.

USAGE
    python3 scripts/case_normalisation.py --case <case>            # derive and write
    python3 scripts/case_normalisation.py --case <case> --verify   # fail if stale or absent
    from case_normalisation import load; n = load(case)            # consumers read, not derive
"""
import argparse
import hashlib
import json
import math
import pathlib
import re
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
NAME = "qoi_normalisation.json"

# OpenFOAM-7 src/functionObjects/field/wallShearStress/wallShearStress.C
#   line 189  compressible    Reff = model.devRhoReff()   -> DYNAMIC   [Pa]
#   line 196  incompressible  Reff = model.devReff()      -> KINEMATIC [m2/s2]
COMPRESSIBLE = ("rhoSimpleFoam", "rhoPimpleFoam", "rhoCentralFoam", "sonicFoam",
                "buoyantSimpleFoam", "buoyantPimpleFoam")


def _get(text, key):
    m = re.search(r"^\s*%s\s+([^;]+);" % re.escape(key), text, re.M)
    return m.group(1).strip() if m else None


def _resolve(text, expr):
    """One-level OpenFOAM $macro substitution."""
    expr, n = expr.strip(), 0
    while expr.startswith("$") and n < 4:
        m = re.search(r"^\s*%s\s+([^;]+);" % re.escape(expr[1:].strip().strip("{}")), text, re.M)
        if not m:
            break
        expr, n = m.group(1).strip(), n + 1
    return expr


def _read(case, rel):
    p = pathlib.Path(case) / rel
    return p.read_text(errors="replace") if p.exists() else None


def _sha(case, rel):
    p = pathlib.Path(case) / rel
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16] if p.exists() else None


def derive(case):
    case = pathlib.Path(case)
    src = ["system/controlDict", "constant/thermophysicalProperties", "0/U", "0/p", "0/T",
           "constant/polyMesh/points", "constant/polyMesh/boundary",
           "constant/turbulenceProperties"]

    cd = _read(case, "system/controlDict")
    if cd is None:
        raise SystemExit("%s: no system/controlDict" % case)
    app = _get(cd, "application")

    # --- flow state, from the field files and the thermophysical model -------------
    u, p0, t0 = _read(case, "0/U"), _read(case, "0/p"), _read(case, "0/T")
    tp = _read(case, "constant/thermophysicalProperties")
    def _freestream(txt, what, n_comp):
        """The freestream value, from internalField if it is uniform and from the FREESTREAM
        BOUNDARY if it is not.

        A WARM-STARTED CASE HAS NO UNIFORM internalField. mapFields and any restart write
        `nonuniform List<vector> ...`, and the old code sliced `len("uniform")` characters off
        that string and regex-scraped whatever remained: on `nonuniform List<vector>` it kept
        "m List<vector>" and tried to float() the "e" of "vector". THE FAILURE WAS LOUD, WHICH
        IS THE ONLY REASON IT IS SAFE -- a parser that scrapes numbers out of a string it does
        not understand can equally well find SOME number and return it, and this one is the
        denominator of every Cp and Cf in the report.

        The freestream patch is the right source in both cases and is what the coefficient is
        normalised on. When internalField IS uniform the two agree by construction, and the
        caller's own self-check (magUInf vs the value derived here) still fires if they do not.
        """
        raw = _get(txt, "internalField") or ""
        if raw.startswith("uniform"):
            v = _resolve(txt, raw[len("uniform"):])
        else:
            m = re.search(r"\n\s{4}(?:freestream|inlet)\s*\n\s{4}\{"
                          r"(?:[^}]*?)(?:freestreamValue|inletValue|value)\s+uniform\s+"
                          r"([^;]+);", txt, re.S)
            if not m:
                raise SystemExit(
                    "%s: 0/%s internalField is %r and no freestream/inlet patch carries a "
                    "uniform value to fall back on. Cannot establish the normalising %s "
                    "without guessing." % (case, what, raw.split("\n")[0][:40], what))
            v = _resolve(txt, m.group(1))
        got = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", v)]
        if len(got) != n_comp:
            raise SystemExit("%s: 0/%s freestream resolves to %r, wanted %d component(s)"
                             % (case, what, v, n_comp))
        return got

    comp = _freestream(u, "U", 3)
    U = math.sqrt(sum(c * c for c in comp))
    pinf = _freestream(p0, "p", 1)[0]
    T = _freestream(t0, "T", 1)[0]
    mu = float(_get(tp, "mu"))
    R = 8314.47 / float(_get(tp, "molWeight"))
    gamma = 1.4
    rho = pinf / (R * T)
    nu = mu / rho
    a = math.sqrt(gamma * R * T)

    # --- geometry, measured from the mesh -------------------------------------------
    pts = np.array(re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)",
                              _read(case, "constant/polyMesh/points")), dtype=float)
    ext = pts.max(axis=0) - pts.min(axis=0)
    span = float(ext.min())                       # the collapsed direction of a 2D mesh
    b = _read(case, "constant/polyMesh/boundary")
    m = re.search(r"(\w+)\s*\{[^}]*?type\s+wall;[^}]*?nFaces\s+(\d+);[^}]*?startFace\s+(\d+);",
                  b, re.S)
    wall, nf, sf = m.group(1), int(m.group(2)), int(m.group(3))
    faces = [np.array(g.split(), dtype=int) for g in
             re.findall(r"\d+\(([\d ]+)\)", _read(case, "constant/polyMesh/faces"))]
    wp = pts[np.unique(np.concatenate([faces[sf + i] for i in range(nf)]))]
    wext = wp.max(axis=0) - wp.min(axis=0)
    sd = int(np.argmin(np.abs(wext - span)))      # the span direction, matched to the domain
    chord = float(max(e for k, e in enumerate(wext) if k != sd))
    thickness = float(min(e for k, e in enumerate(wext) if k != sd))

    # --- the conventions that actually bit us ---------------------------------------
    compressible = app in COMPRESSIBLE
    tau_units = "Pa (dynamic, rho-weighted)" if compressible else "m2/s2 (kinematic)"

    q_dyn = 0.5 * rho * U * U            # Cp and (compressible) Cf denominator
    cf_den = q_dyn if compressible else 0.5 * U * U

    d = {
        "case": case.name,
        "generated_by": "scripts/case_normalisation.py",
        "application": app,
        "compressible": compressible,

        "flow": {
            "U_inf": U, "p_inf": pinf, "T_inf": T, "rho_inf": rho,
            "mu": mu, "nu": nu, "R_specific": R, "gamma": gamma, "a_inf": a,
            "Mach": U / a, "Re_chord": U * chord / nu,
            "_source": "0/U, 0/p, 0/T, constant/thermophysicalProperties; "
                       "rho from p/(R T), no value entered by hand"},

        "geometry": {
            "chord": chord, "span": span, "thickness": thickness,
            "t_over_c": thickness / chord, "wall_patch": wall, "wall_faces": nf,
            "Aref_expected": chord * span,
            "_source": "constant/polyMesh points+boundary+faces, measured"},

        "normalisation": {
            "q_inf": q_dyn,
            "Cp_denominator": q_dyn,
            "Cp_formula": "(p - p_inf) / (0.5 rho_inf U_inf^2)",
            "Cp_note": "for a COMPRESSIBLE solver p is ABSOLUTE pressure in Pa; for an "
                       "incompressible one it is kinematic (p/rho) and rho must be dropped",
            "Cf_denominator": cf_den,
            "Cf_formula": ("|tau_w| / (0.5 rho_inf U_inf^2)" if compressible
                           else "|tau_w| / (0.5 U_inf^2)"),
            "wallShearStress_units": tau_units,
            "wallShearStress_source": (
                "OpenFOAM-7 src/functionObjects/field/wallShearStress/wallShearStress.C: "
                "line 189 compressible Reff=devRhoReff() -> Pa; "
                "line 196 incompressible Reff=devReff() -> m2/s2"),
            "force_coefficient_Aref": float(_get(cd, "Aref")) if _get(cd, "Aref") else None,
            "force_coefficient_lRef": float(_get(cd, "lRef")) if _get(cd, "lRef") else None,
            "force_coefficient_magUInf": float(_get(cd, "magUInf")) if _get(cd, "magUInf") else None,
            "force_coefficient_rhoInf": float(_get(cd, "rhoInf")) if _get(cd, "rhoInf") else None},

        "source_sha256": {s: _sha(case, s) for s in src},
    }

    # internal consistency, checked HERE so a consumer never has to
    n = d["normalisation"]
    problems = []
    if n["force_coefficient_Aref"] is not None:
        e = abs(n["force_coefficient_Aref"] - d["geometry"]["Aref_expected"])
        if e / d["geometry"]["Aref_expected"] > 1e-4:
            problems.append("Aref %.6f but chord x span = %.6f (factor %.4f)" % (
                n["force_coefficient_Aref"], d["geometry"]["Aref_expected"],
                d["geometry"]["Aref_expected"] / n["force_coefficient_Aref"]))
    if n["force_coefficient_rhoInf"] is not None:
        if abs(n["force_coefficient_rhoInf"] - rho) / rho > 2e-3:
            problems.append("forceCoeffs rhoInf %.6f but p/(RT) = %.6f"
                            % (n["force_coefficient_rhoInf"], rho))
    if n["force_coefficient_magUInf"] is not None:
        if abs(n["force_coefficient_magUInf"] - U) / U > 1e-4:
            problems.append("forceCoeffs magUInf %.6f but |0/U| = %.6f"
                            % (n["force_coefficient_magUInf"], U))
    d["consistency"] = {"ok": not problems, "problems": problems}
    return d


def load(case, require_fresh=True):
    """Read the stored normalisation. Consumers use THIS; they do not re-derive.

    Re-deriving in each consumer is how the conventions drifted apart in the first place: two
    scripts each computed 'the' Cf denominator and one of them was wrong.
    """
    p = pathlib.Path(case) / NAME
    if not p.exists():
        raise SystemExit("%s: no %s. Run scripts/case_normalisation.py --case %s"
                         % (case, NAME, case))
    d = json.loads(p.read_text())
    if require_fresh:
        stale = [s for s, h in d["source_sha256"].items() if h != _sha(case, s)]
        if stale:
            raise SystemExit("%s/%s is STALE: %s changed since it was written. Re-run "
                             "scripts/case_normalisation.py --case %s"
                             % (case, NAME, ", ".join(stale), case))
    if not d["consistency"]["ok"]:
        raise SystemExit("%s: normalisation is INCONSISTENT: %s"
                         % (case, "; ".join(d["consistency"]["problems"])))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", required=True, nargs="+")
    ap.add_argument("--verify", action="store_true",
                    help="do not write; fail if absent, stale or inconsistent")
    a = ap.parse_args()
    bad = 0
    for c in a.case:
        c = pathlib.Path(c)
        if a.verify:
            try:
                d = load(c)
                print("  %-24s FRESH   Re %.3e  M %.4f  q_inf %.2f  Cf den %.2f  (%s)"
                      % (c.name, d["flow"]["Re_chord"], d["flow"]["Mach"],
                         d["normalisation"]["q_inf"], d["normalisation"]["Cf_denominator"],
                         d["normalisation"]["wallShearStress_units"].split()[0]))
            except SystemExit as e:
                print("  %-24s %s" % (c.name, e)); bad += 1
            continue
        d = derive(c)
        (c / NAME).write_text(json.dumps(d, indent=2) + "\n")
        f, n, g = d["flow"], d["normalisation"], d["geometry"]
        print("  %s -> %s" % (c.name, NAME))
        print("    U %.4f  p %.1f  T %.2f  rho %.6f  nu %.4e" % (
            f["U_inf"], f["p_inf"], f["T_inf"], f["rho_inf"], f["nu"]))
        print("    Re_c %.4e   M %.4f   chord %.4f   span %.4f   t/c %.4f" % (
            f["Re_chord"], f["Mach"], g["chord"], g["span"], g["t_over_c"]))
        print("    Cp denominator %.4f    Cf denominator %.4f" % (
            n["Cp_denominator"], n["Cf_denominator"]))
        print("    wallShearStress: %s  (%s)" % (n["wallShearStress_units"], d["application"]))
        if not d["consistency"]["ok"]:
            print("    *** INCONSISTENT: %s" % "; ".join(d["consistency"]["problems"]))
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
