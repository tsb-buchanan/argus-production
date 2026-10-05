#!/usr/bin/env python3
"""Derive run-card fields FROM CASE STATE instead of authoring them beside the action.

GEO-080. THE RULE: DERIVE ASSERTIONS FROM STATE; DO NOT AUTHOR THEM ALONGSIDE THE
ACTION. A claim about a file's contents comes from the file; a claim about a run's
outcome comes from the log; nothing is asserted by the step that was supposed to
cause it.

Run cards previously hard-coded `"solver": "simpleFoam"`, `"turbulence_model":
"kOmegaSST"`, `"wall_treatment": "low-Re resolved"`, `"openfoam_version":
"openfoam-org-7"` and the residualControl targets as LITERALS in the collector.
Every one of those is a claim ABOUT THE CASE written by a hand that never read the
case. Edit constant/turbulenceProperties to kEpsilon and the card still says
kOmegaSST, with nothing anywhere to signal it. That is the same defect as a
Dockerfile `|| true`: the step that was supposed to establish the fact is not the
step that reports it.

Every field here is READ BACK from the case:
    solver              system/controlDict `application`, CROSS-CHECKED against the
                        solver log's own `Exec :` line (two independent sources)
    turbulence_model    constant/turbulenceProperties RAS/RASModel
    wall_treatment      the nut boundary condition type on the wall patch
    openfoam_version    the solver log's `Build :` line
    convergence_criteria system/fvSolution residualControl

and each carries a PIN assertion, so a case that silently drifts off the project's
frozen recipe HALTS instead of being relabelled. Wall treatment is the one that
matters most: wall functions are off-limits without a decision entry, and the old
literal would have reported "low-Re resolved" for a wall-function case.

NOT DERIVABLE, and said so rather than quietly derived-looking: cost.machine is a
property of the executing host, not of the case, and stays an input. The log's
`Host :` line is recorded next to it so the two can be compared by a reader.
"""
import argparse
import gzip
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Project pins. A derived value outside these is a HALT, not a relabel.
PIN_OPENFOAM_MAJOR = "7"                 # D009, org-7 pinned project-wide
# THE OF-12 MIGRATION ADDS A SECOND SANCTIONED MAJOR (project ratification,
# 2026-08-26). D009 pinned org-7 project-wide and that pin is NOT withdrawn: OF-7
# cases remain valid and all 22 committed run cards are OF-7. Both are admitted
# because the migration is a TRANSITION and the repo holds cases from both sides of
# it. A major outside this set is still a HALT.
PIN_OPENFOAM_MAJORS = ("7", "12")

# OF-12 REPLACED TWO APPLICATIONS WITH ONE APPLICATION PLUS A SOLVER MODULE, so the
# recipe key is the PAIR, never the application alone:
#
#     OF-7   application simpleFoam;      OF-12  application foamRun; solver incompressibleFluid;
#     OF-7   application rhoSimpleFoam;   OF-12  application foamRun; solver fluid;
#
# KEYING ON `foamRun` ALONE WOULD SILENTLY MERGE THE TWO RECIPES and let a
# wall-function incompressible case pass the very check that exists to catch it,
# because the compressible recipe sanctions wall functions and the incompressible one
# deliberately does not.
#
# THIS IS A MAPPING ONTO THE EXISTING SANCTIONED RECIPES, NOT A NEW RECIPE. The
# equivalence is exactly what the migration's own acceptance tests established:
#   phase 1  argus2d      OF-12 delta +0.716 counts vs OF-7 +0.710; the gate was the
#                         DELTA, not the absolute, and it held.
#   phase 2  cr2d         incompressibleFluid reproduces simpleFoam to max abs dCd
#                         0.40 counts across three grids.
#   phase 4  ONERA M6     fluid reproduces rhoSimpleFoam to -0.06 % on CD and -0.13 %
#                         on CL on the same mesh, surface Cp RMS 1.2e-3 over 100,805
#                         faces, and the OF-7 result was validated against AGARD
#                         AR-138 run 308.
OF12_SOLVER_TO_RECIPE = {
    "incompressibleFluid": "simpleFoam",
    "fluid": "rhoSimpleFoam",
}

# THERE ARE NOW TWO RECIPES, AND ONE PIN CANNOT EXPRESS BOTH (2026-08-21, ARG-170).
# This module pinned simpleFoam + kOmegaSST as THE recipe, which was true when the
# project had one. It now has two: the incompressible wall-resolved recipe (ARG-097)
# and the compressible recipe settled in 2D and validated in 3D on ONERA M6 (ARG-161,
# ARG-165). A single pin does not merely fail to cover the second, it ACTIVELY REFUSES
# it: derive() raised "solver rhoSimpleFoam is off the pinned recipe" on a case that is
# exactly on its own recipe, so no compressible case could ever get a run card and
# record-keeping rule 8 then blocked quoting any of them.
#
# THE PIN IS STILL A HALT, NOT A RELABEL. The recipe is SELECTED by the solver the case
# actually ran (cross-checked against the log's Exec line), and then every other field
# must match THAT recipe. A case mixing rhoSimpleFoam with kOmegaSST still halts, which
# is the property worth keeping: the point was never "one recipe", it was "no case gets
# relabelled to whatever it happens to contain".
RECIPES = {
    "simpleFoam": {
        "name": "incompressible wall-resolved",
        "turbulence": {"kOmegaSST"},
        "record": "ARG-097 (wall-resolved is the standard), D009",
    },
    "rhoSimpleFoam": {
        "name": "compressible transonic (flin2)",
        "turbulence": {"SpalartAllmaras", "kOmegaSST"},
        "record": "ARG-161 (flin2 on SA, kOmegaSST as the control), ARG-165",
        # WALL FUNCTIONS ARE SANCTIONED ON THIS RECIPE AND ONLY THIS ONE. ARG-151
        # sanctioned them for 2D transonic and said the 3D recipe was untouched;
        # ARG-170 extends that to 3D compressible, on the D077 argument that the
        # production intent is wall-modelled at y+ ~80 and validating a recipe we do
        # not use is not validation. The incompressible entry above deliberately has
        # no such key, so a wall-function case there still halts.
        "wall_functions_sanctioned_by": "ARG-151 (2D transonic), ARG-170 (3D compressible)",
    },
}
PIN_SOLVER = "simpleFoam"                # retained: the DEFAULT recipe for 3D production
PIN_TURBULENCE = "kOmegaSST"             # retained: back-compat for callers importing it
# nut BC types that mean WALL FUNCTIONS. Anything not listed is not evidence of a
# wall treatment either way: `calculated` is what a farfield patch carries and says
# nothing about the wall, so it must not be read as "low-Re resolved". Classify by
# the PRESENCE of a wall-function type, never by mapping every patch to a verdict.
NUT_WALL_FUNCTION_TYPES = {
    "nutUWallFunction", "nutkWallFunction", "nutkRoughWallFunction",
    "nutURoughWallFunction", "nutUSpaldingWallFunction", "nutUBlendedWallFunction",
    "nutkAtmRoughWallFunction",
}
# nutLowReWallFunction is NOT a wall function despite the name: it sets nut = 0 at
# the wall and is the resolved-sublayer BC. It must never be added to the set above.
NUT_LOW_RE_TYPES = {"nutLowReWallFunction"}


class DerivationError(RuntimeError):
    """A field could not be read from state, or contradicts a project pin."""


def strip_comments(text):
    """Remove OpenFOAM comments, /* */ then // to end of line.

    A COMMENTED-OUT DIRECTIVE IS NOT STATE, AND READING ONE AS LIVE IS THIS MODULE'S
    OWN FAILURE MODE TURNED INWARD (2026-08-21, ARG-170). Found on ONERA M6, whose
    constant/turbulenceProperties reads:

        //    RASModel            kOmegaSST;
           RASModel SpalartAllmaras;

    `re.search(r"RASModel\\s+(\\w+)\\s*;")` returns the EARLIEST match, so it read the
    DISABLED model and would have stamped both M6 run cards `kOmegaSST` while the runs
    were Spalart-Allmaras. Nothing downstream could have caught it: the card would have
    been internally consistent, schema-valid, and wrong.

    This module exists because run cards once carried `"turbulence_model": "kOmegaSST"`
    as a LITERAL authored by a hand that never read the case. Reading the case and then
    parsing a commented-out line lands in exactly the same place by a longer route, and
    is worse for looking derived. A CORRECT ANSWER REACHED BY AN INVALID ROUTE IS A
    DEFECT, NOT A PASS (GEO-080 item 5); here the answer was not even correct.

    Applied in _read so every parse in this module is covered at once, rather than at
    the one call site that happened to get bitten.
    """
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def _read(case, rel, raw=False):
    p = case / rel
    if not p.exists():
        raise DerivationError("%s: missing %s, so the field it backs cannot be "
                              "derived and must not be asserted" % (case.name, rel))
    t = p.read_text()
    return t if raw else strip_comments(t)


def _one(pattern, text, what, case):
    # re.M is REQUIRED: every pattern here is ^-anchored to a line, and without it
    # ^ binds to the start of the whole file and each lookup silently misses.
    m = re.search(pattern, text, re.M)
    if not m:
        raise DerivationError("%s: could not read %s from state" % (case.name, what))
    return m.group(1)


def solver_log(case):
    """The solver log, under any of the layouts this project actually produces.

    TWO CONVENTIONS EXIST AND ONLY ONE WAS LISTED (2026-08-21, ARG-170). Cluster
    cases write `log.simpleFoam` in the case root; the locally-run validation cases
    (ONERA M6, RAE 2822, NACA 0012) write `logs/<solver>.log` or `log/<solver>.log`,
    because they were built from the upstream tutorials' own Allrun layout. A case
    whose log this function cannot find is reported as having NO RUN OUTCOME, which
    is a true statement about this function and a false one about the case.

    rhoSimpleFoam was also absent from the list entirely, so the compressible
    branch could not be read even where the layout matched.
    """
    # THE SOLVER THE CASE NAMES COMES FIRST, because a fixed preference order picks the
    # WRONG FILE whenever a case holds more than one solver's log, and this project's
    # cases routinely do. Twice now the same failure has been diagnosed backwards:
    #
    #   2026-08-21  log.potentialFoam beat log.simpleFoam.gz -> "controlDict says
    #               simpleFoam but the log was produced by potentialFoam", 49 times,
    #               on 49 perfectly good cases.
    #   2026-08-26  log.simpleFoam.gz (an OF-7 run, 10 Aug) beat log.foamRun (the
    #               OF-12 run, 25 Aug) in cases/of12/cr2d/* -> "controlDict says
    #               foamRun but the log was produced by simpleFoam".
    #
    # Both times the mismatch was REAL and the conclusion WRONG: the files disagreed
    # because the wrong file was read. Reading `application` first removes the guess
    # rather than lengthening the list, so a THIRD solver pairing cannot reproduce it.
    # The full list remains as a fallback so a case with no readable controlDict still
    # gets a diagnosis instead of a crash.
    names = ("foamRun", "simpleFoam", "rhoSimpleFoam", "pimpleFoam", "potentialFoam")
    try:
        app = _one(r"^\s*application\s+(\w+)\s*;", _read(case, "system/controlDict"),
                   "controlDict application", case)
        names = (app,) + tuple(n for n in names if n != app)
    except Exception:
        pass  # no readable controlDict: fall back to the fixed order, and say so below
    # SOLVE LOGS ARE GZIPPED IN THIS REPO, and missing that cost 49 false verdicts
    # (2026-08-21, ARG-172). 50 of 53 cases store `log.simpleFoam.gz`, not
    # `log.simpleFoam`. The bare-name search fell through to `log.potentialFoam`, which
    # every case also has because potentialFoam initialises the field, and the audit
    # then reported "controlDict says simpleFoam but the log was produced by
    # potentialFoam" on 49 perfectly good cases. The mismatch was REAL and the
    # conclusion was WRONG: the two files disagreed because the wrong file was read.
    # Suffixes are tried in preference order per solver, so the ACTIVE solver's log
    # always beats a different solver's, whatever the compression.
    # THREE LAYOUTS, and the v5 cluster one is NESTED PER LEG. `solve_v5.slurm` writes
    # log/leg1/simpleFoam.log, log/leg2/... because a solve is resumed across several
    # jobs. A flat `log/<solver>.log` search misses all of them, so no 3D v5 case could
    # ever have been carded. Found 2026-08-26 by pointing this script at an OF-12 case.
    def candidates(n):
        out = []
        for suf in ("", ".gz"):
            p = case / ("log." + n + suf)                       # cluster flat
            if p.exists():
                out.append(p)
            # A CONTINUATION LOG IS A DIFFERENT FILE, NOT A DIFFERENT SOLVER.
            # cases/of12/cr2d/* hold log.foamRun AND log.foamRun.cont; the second is
            # the later leg, so returning the first derives convergence from a run
            # that was subsequently continued.
            p = case / ("log." + n + ".cont" + suf)
            if p.exists():
                out.append(p)
            for d in ("logs", "log"):                           # tutorial-derived
                p = case / d / (n + ".log" + suf)
                if p.exists():
                    out.append(p)
                for leg in sorted((case / d).glob("leg*")):      # v5 nested per leg
                    p = leg / (n + ".log" + suf)
                    if p.exists():
                        out.append(p)
        return out

    for n in names:
        found = candidates(n)
        if found:
            # NEWEST WINS, derived from the filesystem rather than from an ordering
            # this function asserts. Several legs of ONE solver are the normal case;
            # the last one written is the one whose tail carries the final state.
            best = max(found, key=lambda p: p.stat().st_mtime)
            return str(best.relative_to(case))
    raise DerivationError(
        "%s: no solver log in any known layout (log.<solver>[.cont][.gz] in the case "
        "root, logs/<solver>.log, log/<solver>.log, or log/leg<N>/<solver>.log), so no "
        "run outcome is derivable. Solvers searched: %s"
        % (case.name, ", ".join(names)))


def derive(case):
    case = Path(case)
    log_name = solver_log(case)
    # Only the header is needed and these logs run to hundreds of MB, compressed or not.
    _open = gzip.open if log_name.endswith(".gz") else open
    with _open(case / log_name, "rt", errors="replace") as fh:
        head = "".join(next(fh, "") for _ in range(40))

    ctrl = _read(case, "system/controlDict")
    app = _one(r"^\s*application\s+(\w+)\s*;", ctrl, "controlDict application", case)
    exec_ = _one(r"^Exec\s+:\s+(\S+)", head, "the log's Exec line", case)
    if app != exec_:
        raise DerivationError(
            "%s: controlDict says application %s but the log was produced by %s. "
            "The dictionary and the run disagree; the card may not assert either."
            % (case.name, app, exec_))
    # SELECT the recipe from the solver the case actually ran, then hold it to THAT
    # recipe. Selection is not permission: an unknown solver still halts.
    # Resolve the recipe from the PAIR under OF-12, from the application under OF-7.
    solver_mod = ""          # bound unconditionally: OF-7 cases never enter the branch
    if app == "foamRun":
        solver_mod = re.search(r"^\s*solver\s+(\w+)\s*;", ctrl, re.M)
        solver_mod = solver_mod.group(1) if solver_mod else ""
        if not solver_mod:
            raise DerivationError(
                "%s: application is foamRun but controlDict names no `solver` module. "
                "foamRun ALONE IS NOT A RECIPE: it is the same executable for the "
                "incompressible and compressible paths, which have different "
                "sanctioned turbulence models and opposite wall-function rules. "
                "Refusing to guess which one this is." % case.name)
        if solver_mod not in OF12_SOLVER_TO_RECIPE:
            raise DerivationError(
                "%s: solver module %s maps to no sanctioned recipe (%s). Adding one "
                "needs a decision entry, not an entry in this table."
                % (case.name, solver_mod, ", ".join(sorted(OF12_SOLVER_TO_RECIPE))))
        recipe_key = OF12_SOLVER_TO_RECIPE[solver_mod]
    else:
        recipe_key = app
    if recipe_key not in RECIPES:
        raise DerivationError(
            "%s: solver %s matches no sanctioned recipe (%s). Adding one needs a "
            "decision entry, not an entry in this table."
            % (case.name, app, ", ".join(sorted(RECIPES))))
    recipe = RECIPES[recipe_key]

    build = _one(r"^Build\s+:\s+(\S+)", head, "the log's Build line", case)
    major = build.split("-")[0].split(".")[0]
    if major not in PIN_OPENFOAM_MAJORS:
        raise DerivationError("%s: ran on OpenFOAM %s, sanctioned are org-%s (D009 and "
                              "the OF-12 migration)"
                              % (case.name, build, "/org-".join(PIN_OPENFOAM_MAJORS)))

    # OF-7 named this constant/turbulenceProperties with a `RASModel` key; OF-12 renamed
    # the file to constant/momentumTransport and the key to `model`. Reading only the
    # OF-7 spelling is the defect that made install_postpro.py default every migrated
    # case to kOmegaSST, and it would do the same here. Try OF-12 first, fall back.
    if (case / "constant/momentumTransport").exists():
        turb_txt, turb_key, turb_src = (_read(case, "constant/momentumTransport"),
                                        "model", "momentumTransport")
    else:
        turb_txt, turb_key, turb_src = (_read(case, "constant/turbulenceProperties"),
                                        "RASModel", "turbulenceProperties")
    sim = _one(r"simulationType\s+(\w+)\s*;", turb_txt, "simulationType", case)
    if sim != "RAS":
        raise DerivationError("%s: simulationType %s, expected RAS" % (case.name, sim))
    model = _one(turb_key + r"\s+(\w+)\s*;", turb_txt, turb_key + " in " + turb_src, case)
    if model not in recipe["turbulence"]:
        raise DerivationError(
            "%s: turbulence model %s is off the %s recipe, which admits %s (%s). "
            "Changing it needs a decision entry."
            % (case.name, model, recipe["name"],
               " or ".join(sorted(recipe["turbulence"])), recipe["record"]))

    nut = _read(case, "0/nut")
    nut_bf = re.search(r"boundaryField\s*\{(.*)\}", nut, re.S)
    if not nut_bf:
        raise DerivationError("%s: 0/nut has no boundaryField block" % case.name)
    # Scope to INSIDE boundaryField, or the enclosing block's own name is parsed as
    # a patch and a wall-function case is reported as "mixed treatment": the right
    # verdict on the wrong evidence, which is the defect this module exists to stop.
    patch_nut = dict(re.findall(r"^\s*(\w+)\s*\{[^}]*?type\s+(\w+)\s*;",
                                nut_bf.group(1), re.M | re.S))
    patch_nut.pop("boundaryField", None)
    wf = sorted(p for p, t in patch_nut.items() if t in NUT_WALL_FUNCTION_TYPES)
    lowre = sorted(p for p, t in patch_nut.items() if t in NUT_LOW_RE_TYPES)
    # WALL FUNCTIONS ARE OFF-LIMITS *WITHOUT A DECISION ENTRY*, and the entry is what
    # this branch checks for (2026-08-21, ARG-170). The rule was implemented as a flat
    # prohibition, which is stricter than the project rules actually say and which therefore
    # refused the compressible branch outright: ONERA M6 is wall-modelled BY DESIGN at
    # y+ ~80, because validating wall-resolved and then running wall-modelled would
    # validate a recipe we do not use (D077). The sanction is RECIPE-SCOPED and cites
    # the entry that granted it, so a wall-function case on the incompressible recipe
    # still halts exactly as before.
    if wf and not recipe.get("wall_functions_sanctioned_by"):
        raise DerivationError(
            "%s: nut BC on patch(es) %s is %s, i.e. WALL FUNCTIONS, and the %s recipe "
            "does not sanction them. Off-limits without a decision entry; halting "
            "rather than recording it."
            % (case.name, wf, sorted(set(patch_nut[p] for p in wf)), recipe["name"]))
    if wf and lowre:
        raise DerivationError(
            "%s: MIXED wall treatment. Wall functions on %s, resolved on %s. A card "
            "may not assert one treatment for a case that carries both."
            % (case.name, wf, lowre))
    if not wf and not lowre:
        raise DerivationError(
            "%s: no nut wall BC found in 0/nut, so the wall treatment is not "
            "derivable and must not be asserted. Patches seen: %s"
            % (case.name, sorted(patch_nut)))
    if wf:
        # The schema enum is the vocabulary; the SANCTION travels in its own field
        # rather than being smuggled into the enum value, which would have made the
        # card unvalidatable and the vocabulary meaningless.
        treatment = "wall functions"
        wall_types = {p: patch_nut[p] for p in wf}
    else:
        treatment = "low-Re resolved"
        wall_types = {p: patch_nut[p] for p in lowre}

    # Top-level wall_treatment block (schema d1abb50). Its own commit message
    # asserted "Every card now names nut, k and omega" while no generator was
    # changed and no card ever carried it: the schema and the claim were authored
    # together, and neither was ever checked against a card. Derived here instead.
    # THE WALL PATCH SET AND THE TURBULENCE FIELDS BOTH DEPEND ON THE CASE, and this
    # block assumed one answer to each (2026-08-21, ARG-170). It read `lowre`, which is
    # EMPTY on a sanctioned wall-modelled case, and it demanded `k` and `omega`, which
    # SPALART-ALLMARAS DOES NOT HAVE: its transported variable is nuTilda. On M6 the
    # first produced an IndexError on an empty list, i.e. a crash rather than a
    # diagnosis. Both now follow from what the case actually is.
    wall_patches = wf if wf else lowre
    wall_bc = {"nut": sorted(set(patch_nut[p] for p in wall_patches))[0]}
    TURB_FIELDS = {"kOmegaSST": ("k", "omega"), "SpalartAllmaras": ("nuTilda",)}
    if model not in TURB_FIELDS:
        raise DerivationError(
            "%s: no turbulence-field list known for model %s, so the wall_treatment "
            "block cannot be derived. Add it deliberately rather than omitting the "
            "block, which would let a card claim a wall treatment it never read."
            % (case.name, model))
    for field in TURB_FIELDS[model]:
        txt = _read(case, "0/%s" % field)
        bf = re.search(r"boundaryField\s*\{(.*)\}", txt, re.S)
        if not bf:
            raise DerivationError("%s: 0/%s has no boundaryField block" % (case.name, field))
        types = dict(re.findall(r"^\s*(\w+)\s*\{[^}]*?type\s+(\w+)\s*;", bf.group(1),
                                re.M | re.S))
        hit = [types[p] for p in wall_patches if p in types]
        if not hit:
            raise DerivationError("%s: no %s BC on the wall patch(es) %s"
                                  % (case.name, field, wall_patches))
        wall_bc[field] = sorted(set(hit))[0]

    fvsol = _read(case, "system/fvSolution")
    rc = re.search(r"residualControl\s*\{(.*?)\}", fvsol, re.S)
    criteria = {}
    if rc:
        for key, val in re.findall(r'"?([\w|()]+)"?\s+([0-9.eE+-]+)\s*;', rc.group(1)):
            criteria[key.strip('"')] = float(val)

    return {
        # `foamRun` ALONE DOES NOT IDENTIFY THE RECIPE, so an OF-12 card carries the
        # module too. Emitted only when there is one, because OF-7 has no such concept
        # and the schema forbids unknown keys (model.additionalProperties is false).
        "model": dict({"solver": app, "turbulence_model": model,
                       "wall_treatment": treatment},
                      **({"solver_module": solver_mod}
                         if app == "foamRun" and solver_mod else {})),
        "provenance": {"openfoam_version": "openfoam-org-%s" % major,
                       "openfoam_build": build},
        "numerics": {"convergence_criteria": criteria},
        "host_in_log": (re.search(r'^Host\s+:\s+"?([^"\n]+)"?', head, re.M)
                        or [None, None])[1],
        "wall_nut_bc": {k: v for k, v in wall_types.items()},
        # The patches the treatment was READ FROM, whichever kind they are. This
        # reported `lowre` unconditionally, so a sanctioned wall-modelled case
        # published an EMPTY patch list beside a real wall_treatment string: a
        # verdict over a set the field says was empty (GEO-087).
        "wall_patches": wall_patches,
        "wall_functions_sanctioned_by": recipe.get("wall_functions_sanctioned_by"),
        "wall_treatment_block": wall_bc,
    }


# --- audit: the committed card is a stored ground truth, the derivation is a
# --- DIFFERENT CODE PATH reading the case itself. That is a null with an
# --- independently known answer (D056/D057), not a self-comparison.
def audit(card_path):
    card = json.loads(Path(card_path).read_text())
    case = Path(card_path).parent
    got = derive(case)
    bad = []
    for field in ("solver", "turbulence_model", "wall_treatment"):
        if card.get("model", {}).get(field) != got["model"][field]:
            bad.append(("model.%s" % field, card.get("model", {}).get(field),
                        got["model"][field]))
    if card.get("provenance", {}).get("openfoam_version") != got["provenance"]["openfoam_version"]:
        bad.append(("provenance.openfoam_version",
                    card.get("provenance", {}).get("openfoam_version"),
                    got["provenance"]["openfoam_version"]))
    wt = card.get("wall_treatment")
    if wt is None:
        bad.append(("wall_treatment", "ABSENT",
                    "schema has REQUIRED it since d1abb50; derivable as %s"
                    % got["wall_treatment_block"]))
    else:
        for field, val in got["wall_treatment_block"].items():
            if wt.get(field) != val:
                bad.append(("wall_treatment.%s" % field, wt.get(field), val))

    claimed = card.get("numerics", {}).get("convergence_criteria")
    if claimed is not None and got["numerics"]["convergence_criteria"]:
        norm = {k.replace("|", "|"): v for k, v in got["numerics"]["convergence_criteria"].items()}
        for k, v in claimed.items():
            # The block also carries free-text annotations (e.g. "note"). Only
            # numeric entries are residual targets and comparable to fvSolution.
            if not isinstance(v, (int, float)) or isinstance(v, bool):
                continue
            key = k if k in norm else k.replace("k|omega", "(k|omega)")
            if key not in norm:
                bad.append(("numerics.convergence_criteria[%s]" % k, v, "ABSENT from fvSolution"))
            elif float(norm[key]) != float(v):
                bad.append(("numerics.convergence_criteria[%s]" % k, v, norm[key]))
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--audit-all", action="store_true",
                    help="audit every tracked run_card.json against its case state")
    ap.add_argument("--case", help="derive and print the fields for one case dir")
    args = ap.parse_args()

    if args.case:
        print(json.dumps(derive(args.case), indent=2))
        return 0

    if not args.audit_all:
        ap.error("pass --audit-all or --case")

    cards = sorted(REPO.glob("cases/**/run_card.json"))
    # CARDS WHOSE CASE LIVES ON THE CLUSTER, NOT IN THIS REPO. This audit re-derives a
    # card's claims from its case directory, so a card for an HPC12 case cannot be audited
    # here at all. That is a real limit, but OMITTING them is not an option: before this
    # line existed the audit globbed cases/** only, found 57 cards, and reconciled
    # "57 = 56 + 0 + 1" perfectly while FOUR 3D r2 cards sat outside its glob entirely.
    # A reconciliation that balances over a set it silently narrowed is the
    # narrowed-input-set defect (GEO-087) in the auditor itself. They are named and
    # counted here, with how they ARE verified, rather than left invisible.
    remote_cards = sorted(REPO.glob("results/run_cards_r2/*.json"))
    checked, skipped, failures = 0, [], []
    underivable = []          # counted separately: NOT audited, so not in `checked`
    for c in cards:
        if not (c.parent / "system").exists():
            skipped.append(c)
            continue
        try:
            bad = audit(c)
        except DerivationError as e:
            failures.append((c, [("<derivation>", "", str(e))]))
            underivable.append(c)
            continue
        checked += 1
        if bad:
            failures.append((c, bad))

    print("RUN-CARD ASSERTION AUDIT (GEO-080): are the card's claims true of the case?")
    print("  cards found      %d" % len(cards))
    print("  audited          %d" % checked)
    print("  skipped (no case state on disk) %d" % len(skipped))
    print("  NOT LOCALLY AUDITABLE (case state on HPC12) %d: %s"
          % (len(remote_cards), ", ".join(c.name for c in remote_cards) or "none"))
    print("      these are verified instead by scripts/validate_run_card.py against the "
          "schema, and are built by scripts/build_r2_run_cards.py from facts collected "
          "on the cluster by collect_run_facts.sh, which refuses to emit an invalid card.")
    # COVERAGE ASSERTION (GEO-087): every card found must land in exactly one
    # bucket, or this reports a verdict over an unstated subset. The first version
    # of this check omitted `underivable` and, worse, returned BEFORE the failure
    # list was printed, so it announced 4 unaccounted cards while suppressing the
    # very lines that said which they were. A coverage check that hides the
    # diagnosis is its own defect.
    print("  underivable (halted before audit) %d" % len(underivable))
    unaccounted = len(cards) - (checked + len(skipped) + len(underivable))
    print("  reconciled       %d found = %d audited + %d skipped + %d underivable%s"
          % (len(cards), checked, len(skipped), len(underivable),
             "" if unaccounted == 0 else "  *** %d UNACCOUNTED ***" % unaccounted))
    for c, bad in failures:
        print("  MISMATCH %s" % c.relative_to(REPO))
        for field, was, now in bad:
            print("      %-42s card %-22r state %r" % (field, was, now))
    print()
    if unaccounted:
        print("COVERAGE FAILURE: %d card(s) fell out of every bucket." % unaccounted)
    if failures or unaccounted:
        print("FAIL: %d card(s) assert something the case does not support%s."
              % (len(failures),
                 "" if not unaccounted else ", and %d are unaccounted" % unaccounted))
        return 1
    print("PASS: every audited card's model, OpenFOAM version and convergence")
    print("criteria are reproduced by reading the case, so the previously authored")
    print("literals were CORRECT as well as now being DERIVED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
