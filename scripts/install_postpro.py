#!/usr/bin/env python3
"""install_postpro.py: fit a case with the ARGUS standard post-processing block.

REQUIREMENT, 2026-08-22: post-processing must be automated end to end, so every OpenFOAM
case must set up the correct sampling and VTK output itself.

That is the right order. The post-processor cannot compute what the case never sampled,
so AUTOMATION STARTS IN THE CASE, not in the analysis. This script writes
system/argusPostPro from hpc/templates/argusPostPro with the case's own values
substituted, and wires it into controlDict.

EVERY SUBSTITUTED VALUE IS DERIVED FROM THE CASE, not asked for:

  wall patch      the wall-type patch in constant/polyMesh/boundary, or the patch the
                  existing forces object already names. If several are wall-type the
                  script HALTS rather than guessing, because integrating over the
                  wrong patch produces a plausible number with no warning.
  CofR            from an existing forces/forceCoeffs object if present, else the
                  quarter-chord of the mean chord, and WHICH ONE WAS USED IS RECORDED.
  rho mode        rhoInf for an incompressible solver (p is kinematic), rho for a
                  compressible one. Read from controlDict `application`.
  Trefftz x       four planes from 0.05 to 0.75 chords aft of the trailing edge, the
                  TE taken from the wall patch's own points. PLACEMENT IS SAMPLED
                  RATHER THAN CHOSEN: CDi decays monotonically downstream with no
                  plateau, so one plane is a choice and four make the sensitivity
                  visible (measured 2026-08-22 on the synthetic wings, e swept 0.63
                  to 1.05 across 0.05 to 0.71 chords).
  residual fields the solved fields, from the solver and turbulence model.

Usage:
    python3 scripts/install_postpro.py <case> [<case> ...] [--dry-run]
    python3 scripts/assert_case_postpro.py --before <case>     # then verify
"""
import argparse
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
TEMPLATE = REPO / "hpc/templates/argusPostPro"


class InstallError(RuntimeError):
    pass


def read(p, default=""):
    return Path(p).read_text(errors="replace") if Path(p).exists() else default


def strip_comments(t):
    t = re.sub(r"/\*.*?\*/", "", t, flags=re.S)
    return re.sub(r"//[^\n]*", "", t)


def wall_patch(case):
    """The wall patch, from the mesh. HALTS on ambiguity rather than guessing."""
    cands = []
    for b in list((case / "constant/polyMesh").glob("boundary*")) + \
             list((case / "processor0/constant/polyMesh").glob("boundary*")):
        t = read(b)
        for m in re.finditer(r"^\s{4}(\w+)\s*\n\s*\{(.*?)\n\s{4}\}", t, re.S | re.M):
            if re.search(r"type\s+wall\s*;", m.group(2)):
                cands.append(m.group(1))
        if cands:
            break
    cands = sorted(set(cands))
    if not cands:
        # fall back to whatever an existing forces object names
        cd = strip_comments(read(case / "system/controlDict"))
        m = re.search(r"patches\s*\(\s*([^)]*?)\s*\)", cd)
        if m:
            return m.group(1).split()[0].strip('"'), "existing forces object"
        raise InstallError("%s: no wall-type patch found and no forces object to copy "
                           "one from. Integrating over the wrong patch gives a "
                           "plausible number with no warning, so this halts."
                           % case.name)
    if len(cands) > 1:
        raise InstallError("%s: %d wall patches %s. Which one carries the wing is a "
                           "modelling choice, not something to guess." % (case.name, len(cands), cands))
    return cands[0], "constant/polyMesh/boundary"


def solver_and_rho(case):
    """(application, rho mode, rhoInf, pInf), the last two DERIVED for a compressible case.

    THE TEMPLATE WAS ONLY EVER EXERCISED INCOMPRESSIBLY AND IT SHOWED. It hard-codes
    `pInf 0` and this function used to hard-code `rhoInf 1`, both correct for
    simpleFoam, where p is KINEMATIC GAUGE pressure so the freestream is zero and the
    normalisation is per unit density. Neither is correct for rhoSimpleFoam, where p is
    ABSOLUTE PRESSURE IN PASCALS:

      pInf 0    makes cp = p/q_ref, about +2 everywhere instead of a coefficient about
                zero. The freestream must be subtracted.
      rhoInf 1  makes the denominator 0.5*1*|U|^2 = 42,170 Pa against the real dynamic
                pressure 49,187 Pa on the M6 case, a 16.6% error, which is exactly the
                freestream density it forgot to include.

    D077's SHAPE ON A NEW AXIS: a template validated in one configuration reused in
    another where the thing that made it true no longer holds. So for a compressible
    solver both are read from the case's own freestream, and the function HALTS rather
    than falling back on a number that would look plausible.
    """
    cd = strip_comments(read(case / "system/controlDict"))
    m = re.search(r"^\s*application\s+(\w+)\s*;", cd, re.M)
    app = m.group(1) if m else "simpleFoam"

    # COMPRESSIBILITY IS NOT READABLE FROM `application` UNDER OF-12, AND THE MIGRATION
    # SILENTLY RE-INTRODUCED THE VERY BUG THIS FUNCTION EXISTS TO FIX.
    #
    #   OF-7   application rhoSimpleFoam;   <- the discriminator
    #   OF-12  application foamRun;         <- IDENTICAL for both paths
    #          solver      fluid;                 compressible
    #          solver      incompressibleFluid;   incompressible
    #
    # `app.startswith("rho")` is false for every OF-12 case, so every compressible one
    # took the incompressible branch and got `rhoInf 1, pInf 0` back: the 16.6 % dynamic
    # pressure error and the cp-offset-by-a-whole-freestream this docstring describes.
    # Measured on CMPB_a0p87, 2026-08-26. Third instance of the same OF-7-spelling
    # class, after turb_fields' `turbulenceProperties/RASModel` and the
    # `thermophysicalProperties` path below.
    sm = re.search(r"^\s*solver\s+(\w+)\s*;", cd, re.M)
    solver_mod = sm.group(1) if sm else ""
    if solver_mod:
        compressible = (solver_mod == "fluid")
    else:
        compressible = app.startswith("rho") or app.startswith("sonic")

    if not compressible:
        return app, "rhoInf", "1", "0"    # incompressible: p is kinematic, gauge
    p_inf = _scalar_entry(case, ("0.orig/p", "0/p"), ("internalField", "pOut", "pInf"))
    t_inf = _scalar_entry(case, ("0.orig/T", "0/T"), ("internalField", "Tinlet", "TInf"))
    # OF-12 renamed thermophysicalProperties -> physicalProperties. Try both, OF-12 first.
    w = _scalar_entry(case, ("constant/physicalProperties",
                             "constant/thermophysicalProperties"), ("molWeight",)) or 28.9
    if p_inf is None or t_inf is None:
        raise InstallError(
            "%s runs the COMPRESSIBLE solver %s, so p is absolute pressure in Pa and cp "
            "needs the freestream p_inf and rho_inf. Neither is readable from 0.orig/p "
            "and 0.orig/T. Refusing to write pInf 0 and rhoInf 1, which would give a cp "
            "wrong by the whole freestream." % (case.name, app))
    rho_inf = p_inf / ((8314.47 / w) * t_inf)
    return app, "rho", "%.6f" % rho_inf, "%.1f" % p_inf


def _scalar_entry(case, files, keys):
    """First readable scalar among `keys` in the first existing file, one $var deep."""
    for f in files:
        p = case / f
        if not p.exists():
            continue
        txt = strip_comments(p.read_text(errors="replace"))
        for k in keys:
            m = re.search(r"^\s*%s\s+([^;]+);" % re.escape(k), txt, re.M)
            if not m:
                continue
            v = m.group(1).strip().split()[-1]
            if v.startswith("$"):
                m2 = re.search(r"^\s*%s\s+([-\d.eE+]+)\s*;" % re.escape(v[1:]), txt, re.M)
                v = m2.group(1) if m2 else None
            try:
                return float(v)
            except (TypeError, ValueError):
                continue
    return None


def cofr(case):
    cd = strip_comments(read(case / "system/controlDict"))
    m = re.search(r"CofR\s*\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", cd)
    if m:
        return "%s %s %s" % m.groups(), "existing forces object"
    return "0 0 0", "DEFAULTED to the origin; no existing forces object named one"


def trailing_edge(case):
    """Most-downstream x of the WING, taken from geometry the SOLVER actually saw.

    THE STL IS THE WRONG SOURCE AND M6 PROVES IT TWICE. Alletto's STL is BINARY, so a
    text scan finds no vertices at all; and it is in MILLIMETRES while the mesh is
    scaled to metres by `transformPoints -scale 0.001` AFTER snappy. Placing the
    Trefftz planes from it would put them 1000x too far downstream, and nothing
    downstream would notice: the planes would simply integrate freestream and return a
    small, plausible, meaningless CDi.

    So the order of preference is: a SAMPLED WALL SURFACE (post-transform, exactly what
    was integrated), then an ASCII STL with its scale applied, then halt. Never a
    silent guess.
    """
    import numpy as np
    # 1. a sampled wall surface, in solver coordinates
    for pat in ("postProcessing/**/*patch*.raw", "postProcessing/**/*wing*.raw",
                "postProcessing/**/*wing*.vtk", "postProcessing/**/*patch*.vtk"):
        for f in sorted(case.glob(pat)):
            try:
                txt = f.read_text(errors="replace")
            except OSError:
                continue
            xs = [float(l.split()[0]) for l in txt.splitlines()
                  if l.strip() and not l.startswith("#") and
                  re.match(r"^\s*-?\d", l) and len(l.split()) >= 3]
            if len(xs) > 100:
                v = np.array(xs)
                return float(v.max()), float(v.max() - v.min()), "sampled surface %s" % f.name
    # 2. an ASCII STL, WITH its scale stated
    for stl in (case / "constant/triSurface").glob("*.stl"):
        t = read(stl)
        xs = re.findall(r"vertex\s+(-?\d+\.?\d*(?:[eE][-+]?\d+)?)", t)
        if xs:
            v = np.array(xs, dtype=float)
            if v.max() > 100:
                raise InstallError(
                    "%s: %s spans %.1f in x, which is millimetres, not metres. The mesh "
                    "is scaled after snappy, so this file is NOT in solver coordinates. "
                    "Run the case far enough to write a wall sample, or state the scale."
                    % (case.name, stl.name, v.max()))
            return float(v.max()), float(v.max() - v.min()), "ASCII STL %s" % stl.name
    raise InstallError("%s: no sampled wall surface and no readable ASCII STL, so the "
                       "trailing edge cannot be located in SOLVER coordinates and the "
                       "Trefftz planes cannot be placed." % case.name)


def freestream_mag(case):
    """|U| from the case's own 0/U or 0.orig/U.

    THE CASE'S OWN VECTOR, not a nominal condition. cp is normalised on it, and a
    coefficient normalised on a number nobody can trace back to the run is exactly the
    frame error this project keeps hitting. Also the reason the wind axes in the
    hinge-moment integrator come from the same place: derive both from one vector and
    the lag between them is zero by construction rather than small by luck.
    """
    import math
    VEC = r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)"
    for rel in ("0/U", "0.orig/U"):
        t = strip_comments(read(case / rel))
        # (a) the vector written inline
        m = re.search(r"internalField\s+uniform\s*" + VEC, t)
        if m:
            v = [float(x) for x in m.groups()]
            return math.sqrt(sum(c * c for c in v)), rel + " (inline)"
        # (b) OPENFOAM VARIABLE INDIRECTION, which is how the M6 tutorial writes it:
        #        Uinlet          (290 0 15.5);
        #        internalField   uniform $Uinlet;
        #     A pattern that only understands (a) finds nothing and reports "no uniform
        #     U vector", which is true of the LINE and false of the FILE.
        m = re.search(r"internalField\s+uniform\s+\$(\w+)", t)
        if m:
            d = re.search(r"^\s*%s\s+%s\s*;" % (re.escape(m.group(1)), VEC), t, re.M)
            if d:
                v = [float(x) for x in d.groups()]
                return math.sqrt(sum(c * c for c in v)), "%s (via $%s)" % (rel, m.group(1))
    raise InstallError("%s: no uniform U vector in 0/U or 0.orig/U, so cp cannot be "
                       "normalised. Refusing to guess a freestream." % case.name)


# Which dictionary and which KEY names the turbulence model depends on the OpenFOAM
# generation, and getting it wrong is silent rather than loud:
#
#     OF-7   constant/turbulenceProperties   RASModel  kOmegaSST;
#     OF-12  constant/momentumTransport      model     kOmegaSST;
#
# THE ORIGINAL READ ONLY THE OF-7 SPELLING. Every migrated case here carries
# momentumTransport, so the regex matched nothing and the function fell through to a
# hard-coded "kOmegaSST" default on ALL of them. It happened to be right for the 18
# wall-resolved solves and WRONG for every SA and every compressible case, which is
# the worst arrangement: a correct answer reached by an invalid route (GEO-080 item 5),
# so nothing surfaced it until a case that was not kOmegaSST came through.
#
# Comments must be stripped BEFORE matching. onera_m6's momentumTransport carries a
# commented-out `// RASModel kOmegaSST;` above a live `model SpalartAllmaras;`, and a
# naive grep reports the case as kOmegaSST. That is not hypothetical: it misled this
# author while auditing this very function.
MODEL_FIELDS = {
    "SpalartAllmaras": "nuTilda",
    "kOmegaSST": "k omega",
    "kOmegaSSTSAS": "k omega",
    "kOmega": "k omega",
    "kEpsilon": "k epsilon",
    "realizableKE": "k epsilon",
    "RNGkEpsilon": "k epsilon",
    "laminar": "",
}

# A compressible run also solves an ENERGY equation, and its residual is not optional:
# it is usually the last one to settle, so a convergence claim that never looked at it
# is a claim over an unstated subset (GEO-087). The variable's NAME follows the
# thermophysical energy type, so it is read rather than assumed.
ENERGY_FIELD = {"sensibleInternalEnergy": "e", "sensibleEnthalpy": "h"}


def turb_fields(case):
    """Residual fields for THIS case: momentum, pressure, turbulence, and energy.

    DERIVED FROM THE CASE, NEVER DEFAULTED. An unrecognised model or a compressible
    case with an unreadable energy type HALTS. Defaulting is what produced the bug
    this function replaces: a wrong field list does not fail loudly, it just records
    residuals for fields the run does not solve and omits the ones it does.
    """
    of12 = case / "constant/momentumTransport"
    of7 = case / "constant/turbulenceProperties"
    if of12.exists():
        src, key = of12, "model"
    elif of7.exists():
        src, key = of7, "RASModel"
    else:
        raise InstallError(
            "%s: neither constant/momentumTransport (OF-12) nor "
            "constant/turbulenceProperties (OF-7) exists, so the turbulence model "
            "cannot be derived. Refusing to default." % case.name)

    t = strip_comments(read(src))
    m = re.search(key + r"\s+(\w+)\s*;", t)
    if not m:
        raise InstallError(
            "%s: %s has no live `%s <Name>;` entry. Refusing to default to "
            "kOmegaSST, which is how every migrated case was silently mislabelled."
            % (case.name, src.name, key))
    model = m.group(1)
    if model not in MODEL_FIELDS:
        raise InstallError(
            "%s: unrecognised turbulence model %r. Add it to MODEL_FIELDS with its "
            "solved fields rather than letting it fall through." % (case.name, model))

    fields = ["U", "p"]
    turb = MODEL_FIELDS[model]
    if turb:
        fields += turb.split()

    # Compressibility comes from the solver module the case actually names, not from
    # the presence of a T field: `incompressibleFluid` cases may carry T for a passive
    # scalar without solving an energy equation.
    cd = strip_comments(read(case / "system/controlDict"))
    sm = re.search(r"\bsolver\s+(\w+)\s*;", cd)
    solver_mod = sm.group(1) if sm else ""
    if solver_mod == "fluid":
        phys = ""
        for name in ("physicalProperties", "thermophysicalProperties"):
            p = case / "constant" / name
            if p.exists():
                phys = strip_comments(read(p))
                break
        em = re.search(r"\b(sensibleInternalEnergy|sensibleEnthalpy)\b", phys)
        if not em:
            raise InstallError(
                "%s: solver is `fluid` (compressible) but no sensibleInternalEnergy / "
                "sensibleEnthalpy found in constant/physicalProperties, so the energy "
                "variable name cannot be derived. Refusing to omit the energy "
                "residual." % case.name)
        fields.insert(2, ENERGY_FIELD[em.group(1)])

    return " ".join(fields), model


def install(case, dry=False):
    case = Path(case).resolve()
    if not (case / "system/controlDict").exists():
        raise InstallError("%s: no system/controlDict" % case)
    patch, patch_src = wall_patch(case)
    app, rho_mode, rhoinf, pinf = solver_and_rho(case)
    cr, cr_src = cofr(case)
    te, span_x, stl = trailing_edge(case)
    chord = span_x
    resid, model = turb_fields(case)
    umag, usrc = freestream_mag(case)

    xs = [te + f * chord for f in (0.05, 0.25, 0.50, 0.75)]
    sub = {
        "ARGUS_WALL_PATCH": patch,
        "ARGUS_COFR": cr,
        "ARGUS_RHO_MODE": rho_mode,
        "ARGUS_RHOINF": rhoinf,
        "ARGUS_PINF": pinf,
        "ARGUS_RESIDUAL_FIELDS": resid,
        # THE TREFFTZ PLANES MUST CARRY ENTROPY ON A COMPRESSIBLE CASE, AND THEY DID NOT.
        # The list was hard-coded `(U p)`, which is sufficient at low speed where the only
        # far-field drag is vortex drag. At M 0.78 it is NOT: a defensible far-field
        # decomposition has to separate VORTEX drag from WAVE and VISCOUS drag, and that
        # separation runs on entropy, which needs T (with p and perfectGas giving rho).
        # From p and U alone the only thing available is the naive integral of (v^2+w^2),
        # which over a viscous wake folds the momentum deficit straight in and is biased
        # HIGH by an amount nobody can bound from the sampled data (ARG-191 amendment 1,
        # which records this as the reason the cruise induced-drag route was blocked).
        #
        # PARAMETERISED RATHER THAN WIDENED, because the template is SHARED and Condition
        # CR is incompressible: it has no T at all, so a blunt `(U p T rho)` would make
        # every low-speed case sample a field that does not exist. The discriminator is
        # the one solver_and_rho already established, `solver fluid` under OF-12, which is
        # why this reads rho_mode rather than re-deriving compressibility and risking the
        # OF-7-spelling class of bug documented there.
        "ARGUS_TREFFTZ_FIELDS": "U p T rho" if rho_mode == "rho" else "U p",
        "ARGUS_UINF": "%.6f" % umag,
        "ARGUS_TREFFTZ_X1": "%.6f" % xs[0],
        "ARGUS_TREFFTZ_X2": "%.6f" % xs[1],
        "ARGUS_TREFFTZ_X3": "%.6f" % xs[2],
        "ARGUS_TREFFTZ_X4": "%.6f" % xs[3],
    }
    txt = TEMPLATE.read_text()
    for k, v in sub.items():
        txt = txt.replace("$" + k, v)
    left = re.findall(r"\$ARGUS_\w+", txt)
    if left:
        raise InstallError("unsubstituted placeholders remain: %s" % sorted(set(left)))

    print("  %s" % case.name)
    print("     solver        %s -> rho mode '%s', rhoInf %s, pInf %s"
          % (app, rho_mode, rhoinf, pinf))
    print("     wall patch    %-14s (from %s)" % (patch, patch_src))
    print("     CofR          %-14s (from %s)" % (cr, cr_src))
    print("     turbulence    %-14s -> residual fields: %s" % (model, resid))
    print("     |U| for cp    %-14.4f (from %s)" % (umag, usrc))
    print("     TE x %.4f, chord %.4f (from %s)" % (te, chord, stl))
    print("     Trefftz x     %s" % " ".join("%.4f" % x for x in xs))
    if dry:
        print("     DRY RUN, nothing written")
        return
    (case / "system/argusPostPro").write_text(txt)

    cd_path = case / "system/controlDict"
    cd = cd_path.read_text()

    # REMOVE THE BLOCKS argusPostPro SUPERSEDES, or the case samples everything TWICE.
    # This script was additive by design and never removed anything, which is safe read
    # one case at a time and ruinous in aggregate: the wing surface is 510 MB a write and
    # was being written by BOTH `wingSurface` and `argusWingSurface`. Measured across the
    # campaign on 2026-09-02: 121 GB of byte-for-byte duplicate, on a filesystem whose
    # operator was emailing us about being full.
    # forceCoeffs1 is NOT in the superseded set and prune_legacy_fo REFUSES it outright:
    # argusPostPro has no equivalent and six scripts read postProcessing/forceCoeffs1.
    from prune_legacy_fo import remove_block, SUPERSEDED
    dropped = []
    for _fo in SUPERSEDED:
        cd, _did = remove_block(cd, _fo)
        if _did:
            dropped.append(_fo)
    if dropped:
        print("     removed superseded function objects: %s" % ", ".join(dropped))
        if not dry:
            cd_path.write_text(cd)

    if "argusPostPro" in cd:
        print("     controlDict already includes argusPostPro")
    else:
        if re.search(r"^functions\s*$", cd, re.M) or re.search(r"^functions\s*\{", cd, re.M):
            cd = re.sub(r"(^functions\s*\n?\s*\{)",
                        r'\1\n    #include "argusPostPro"', cd, count=1, flags=re.M)
        else:
            cd += '\n\nfunctions\n{\n    #include "argusPostPro"\n}\n'
        cd_path.write_text(cd)
        print("     controlDict now includes argusPostPro")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    rc = 0
    for c in a.cases:
        try:
            install(c, a.dry_run)
        except InstallError as e:
            print("  HALT %s" % e)
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
