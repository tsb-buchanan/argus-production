#!/usr/bin/env python3
"""solver_recipe.py: the compressible solver recipe as ONE NAMED CHOICE, applied and verified.

WHY THIS EXISTS (project plan, 2026-08-10). A recipe is a SET of settings across three
dictionaries, and every previous attempt applied it as a scatter of independent regex patches.
That is how the "damped" runs came back BIT-IDENTICAL to the undamped ones (ARG-156 item 11a):
one anchor missed, nothing said so, and the run reported a configuration it was not in. A
half-applied recipe is indistinguishable from an applied one in every artefact anyone reads.

TWO STRUCTURAL DEFENCES, not diligence:

1. THE FILES ARE WRITTEN WHOLE, NOT PATCHED. system/fvSchemes and system/fvSolution are
   generated from the spec below, with NO commented alternatives in them. The wolfdynamics
   dictionaries ship the tutorial recipe as a COMMENTED BLOCK sitting directly beneath the live
   one, and `$urf1` variables above that. In such a file `x in text` is never a valid test, and
   this project has now been bitten by that four separate times (`urf_rho` swallowing a `rho`
   check, a commented `turbulence` line, a `pInf` read out of a /* */ block, a `rho` entry that
   was never written). Generated files have nothing to match by accident.

2. THE VERIFICATION IS A DIFFERENT PARSER FROM THE WRITER (GEO-080). verify() does not re-read
   its own regexes: it shells out to `foamDictionary -expand`, i.e. OpenFOAM's own dictionary
   reader, which strips comments and expands $variables, and compares the expanded result
   against the spec. Nothing here vouches for a step it did not observe.

THE THREE RECIPES, and what separates them:

    wolf      what every case in this repo currently runs. rhoSimpleFoam, SIMPLEC,
              `transonic yes`, an EMPTY fields{} relaxation block.
    tutorial  OpenFOAM-7's own tutorials/compressible/rhoSimpleFoam/aerofoilNACA0012, at
              U 250 m/s. Plain SIMPLE, `transonic` off, fields{p 0.7; rho 0.01;}.
    lts       rhoPimpleFoam with localEuler local time stepping. No relaxation factors exist
              at all; robustness comes from the local pseudo-timestep instead.

WHY `transonic yes` IS THE SUSPECT, read off the source rather than inferred. OpenFOAM-7
rhoSimpleFoam pEqn.H (plain SIMPLE) and pcEqn.H (SIMPLEC) BOTH end:

    rho = thermo.rho();
    if (!simple.transonic()) { rho.relax(); }

so under `transonic yes` THE DENSITY IS NEVER RELAXED, on either algorithm. `p.relax()` is
called unconditionally in both branches -- but it is guarded inside GeometricField::relax() by
mesh().relaxField("p"), which is FALSE when the fields{} block has no p entry. The wolf recipe
has an empty fields{} block. So it runs with no density damping AND no pressure-field damping,
and the only damping left is pEqn.relax() plus the equation factors. That is survivable at
Re 6.5e6 and is not at Re 1.7-2.1e7.

(CORRECTION, recorded rather than silently fixed. The plan text said p.relax() runs only in the
non-transonic branch. It does not: it is unconditional. The load-bearing claim -- that rho is
never relaxed under transonic yes, and that an empty fields{} block leaves p unrelaxed too --
is unaffected, and is what the experiment tests. Only the mechanism description moved.)

THE PRESSURE CLAMPS ARE A PHYSICAL INPUT, NOT A TUNING KNOB, and that is why
pressure_clamp_admits() is in here. OpenFOAM-7 pressureControl.C takes pMin/pMax from the
VALUE-FIXING patches and multiplies by pMinFactor/pMaxFactor, so on our aerofoil (outlet fixed
at 1e5 Pa) pMinFactor 0.5 clamps static pressure at 50,000 Pa. THE AGARD CASE 9 SUCTION PEAK IS
51,144 Pa. The clamp is 2.3% below the measured answer, and the existing case09 run shows it
firing 19 times. A limiter transplanted from a tutorial has no idea what pressure range the new
case needs, and a clamped run converges beautifully to the wrong Cp (ARG-096 item 10: a gate on
the output cannot see a wrong input). The gate here is DERIVED from the case's own q_inf and a
declared Cp floor, not chosen.
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
FOAM = ("export PATH=/usr/bin:$PATH; source /opt/openfoam7/etc/bashrc >/dev/null 2>&1; "
        "export PATH=/usr/bin:$PATH; ")

# Cp floor the pressure clamp must admit. AGARD measures -1.31 on RAE case 9 and -1.22 on
# case 6; a converging solution overshoots below its final peak, so the clamp is required to
# leave headroom rather than merely clear the answer. -2.0 is the declared floor: it is well
# below anything these cases can physically reach and still leaves p positive.
CP_FLOOR = -2.0

# ---------------------------------------------------------------------------
# The gradient schemes referenced BY NAME from divSchemes. `linearUpwind limitedT` looks up a
# gradScheme called limitedT, so these entries are load-bearing and cannot be flattened away.
GRAD_SCHEMES = """gradSchemes
{
    limitedU        cellMDLimited Gauss linear 1;
    limitedT        cellLimited Gauss linear 1;
    limitedE        cellLimited Gauss linear 1;
    default         Gauss linear;
    grad(U)         $limitedU;
    grad(k)         $limitedT;
    grad(omega)     $limitedT;
}"""

# Common to all three arms: only the entries named in the plan's comparison table move.
DIV_COMMON = {
    "default": "none",
    "div(phi,U)": "bounded Gauss Minmod",
    "div(phi,k)": "bounded Gauss linearUpwind limitedT",
    "div(phi,omega)": "bounded Gauss linearUpwind limitedT",
    "div(phi,h)": "bounded Gauss limitedLinear 1",
    "div(phi,K)": "bounded Gauss limitedLinear 1",
    "div(((rho*nuEff)*dev2(T(grad(U)))))": "Gauss linear",
}

# ---- `bounded` IS A STEADY-STATE DEVICE AND MUST COME OFF FOR THE TRANSIENT ARM.
# boundedConvectionScheme subtracts fvm::Sp(fvc::div(phi), psi). At steady convergence
# div(phi) is zero, so the term vanishes and all it does is help the path there. IN A
# TRANSIENT COMPRESSIBLE RUN div(phi) IS NOT ZERO -- it is the physical -drho/dt -- so the
# term subtracts a real part of the transport equation on every step.
#
# MEASURED, NOT ARGUED. A 2x2 probe on this exact mesh and case, 1,000 LTS steps per cell:
#     bounded   + transonic yes   CRASHED at step 559, h residual at 1.0 on 244 steps, |Ux| 1.4e13
#     bounded   + transonic no    CRASHED at step 577, h residual at 1.0 on 451 steps, |Ux| 1.5e13
#     unbounded + transonic yes   ran clean, h never reached 1, |Ux| 283 m/s
#     unbounded + transonic no    ran clean, h never reached 1, |Ux| 302 m/s
# It is `bounded`, at both levels of `transonic`, and it is not marginal. A further 1,500-step
# probe confirmed the SINGLE change is enough: unbounded while keeping A and B's limitedLinear
# energy scheme also runs clean (|Ux| 269 m/s), so arm C does not have to take the LTS
# tutorial's first-order upwind energy and be penalised for accuracy it never needed to lose.
#
# CORROBORATION FROM THE DISTRIBUTION: of the 39 OpenFOAM-7 tutorials whose fvSchemes contains
# `bounded Gauss`, 38 are steadyState. The single exception is an INCOMPRESSIBLE pisoFoam case,
# where div(phi) is identically zero and the term cannot bite.
#
# THIS IS THE SAME DEFECT AS THE PRESSURE CLAMP, one level down: a setting inherited from a
# configuration where it was correct, carried into one where it is not, and invisible to every
# output gate because the run simply dies with no statement of why.
DIV_TRANSIENT = {k: v.replace("bounded ", "") for k, v in DIV_COMMON.items()}

SOLVERS_STEADY = """solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-6;
        relTol          0.001;
        minIter         2;
    }
    U
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0;
        minIter         2;
    }
    "(e|h)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0;
        minIter         2;
    }
    "(k|omega)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0;
        minIter         2;
    }
}"""

# rhoPimpleFoam asks for FOUR things rhoSimpleFoam never does (ARG-145 amendment 1): a `rho`
# entry, Final-corrector entries for every solved field, div(phiv,p), and PIMPLE controls as
# single values. All four are here rather than patched in later.
SOLVERS_LTS = """solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
        minIter         2;
    }
    pFinal
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0;
        minIter         2;
    }
    "(rho|U|e|h|k|omega)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0.01;
        minIter         2;
    }
    "(rho|U|e|h|k|omega)Final"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0;
        minIter         2;
    }
}"""

RECIPES = {
    # ---- ARM A. What every case in this repo runs today, WITH ONE DELIBERATE DEPARTURE.
    # The historical cases carry pMinFactor 0.5, i.e. a static-pressure floor of 50,000 Pa
    # against an AGARD case 9 suction peak of 51,144 Pa, and the existing case09 log shows
    # pressureControl firing 19 times. Left at 0.5, arm A would be handicapped by a defect
    # that has nothing to do with the question this experiment asks, and B and C would get
    # credit for fixing it. D071's test decides it: the clamp is a defect COMMON TO BOTH SIDES
    # of the comparison only if all three arms carry the same one, so it is corrected in all
    # three and cancels. Arm A therefore differs from the historical runs in this one respect,
    # which is stated rather than absorbed.
    "wolf": dict(
        application="rhoSimpleFoam",
        ddt="steadyState",
        algo_block="SIMPLE",
        div=dict(DIV_COMMON, **{"div(phid,p)": "Gauss limitedLinear 1"}),
        algo={"consistent": "yes", "transonic": "yes",
              "nNonOrthogonalCorrectors": "1", "pMinFactor": "0.1", "pMaxFactor": "2"},
        fields_relax={},
        equations_relax={"p": "0.3", "U": "0.3", "e": "0.3", "h": "0.3",
                         '"(k|omega)"': "0.3"},
        residual_control={"p": "1e-6", "U": "1e-6",
                          '"(k|omega|e|h)"': "1e-6"},
        solvers=SOLVERS_STEADY,
    ),
    # ---- ARM B. OpenFOAM-7 tutorials/compressible/rhoSimpleFoam/aerofoilNACA0012.
    # NOTE THE ENERGY KEY. The shipped tutorial relaxes `e`, because it selects
    # sensibleInternalEnergy. Our thermophysicalProperties selects sensibleEnthalpy, so the
    # solved field is `h` and an `e` entry RELAXES NOTHING -- the lookup is by field name.
    # Transplanting the tutorial block verbatim would have left the energy equation
    # completely unrelaxed while the dictionary looked right. Both keys are written.
    # `p` is deliberately ABSENT from equations{}: the non-transonic branch of pEqn.H/pcEqn.H
    # never calls pEqn.relax(), so a p equation factor there is dead, and a dead entry in a
    # generated file is a lie about what is damping the solve.
    "tutorial": dict(
        application="rhoSimpleFoam",
        ddt="steadyState",
        algo_block="SIMPLE",
        div=dict(DIV_COMMON, **{"div(phid,p)": "Gauss upwind",
                                "div((phi|interpolate(rho)),p)": "bounded Gauss upwind"}),
        algo={"consistent": "no", "transonic": "no",
              "nNonOrthogonalCorrectors": "1", "pMinFactor": "0.1", "pMaxFactor": "2"},
        fields_relax={"p": "0.7", "rho": "0.01"},
        equations_relax={"U": "0.3", "e": "0.7", "h": "0.7", '"(k|omega)"': "0.7"},
        residual_control={"p": "1e-6", "U": "1e-6",
                          '"(k|omega|e|h)"': "1e-6"},
        solvers=SOLVERS_STEADY,
    ),
    # ---- ARM C. LTS, from tutorials/compressible/rhoPimpleFoam/RAS/angledDuctLTS.
    # ITS PRESSURE CLAMPS DO NOT COME WITH IT, and this is the one deliberate departure from
    # the plan's pre-registered table. angledDuctLTS runs pMinFactor 0.9 / pMaxFactor 1.5,
    # which is fine in a duct where static pressure barely moves and is FLATLY WRONG on a
    # transonic aerofoil: 0.9 clamps at 90,000 Pa against a measured suction peak of 51,144 Pa
    # and a stagnation pressure of 142,540 Pa, so it would clip BOTH ends of the real field.
    # The arm would have returned a beautifully converged wrong answer. Clamps are set from
    # the physics instead and checked by pressure_clamp_admits().
    "lts": dict(
        application="rhoPimpleFoam",
        ddt="localEuler",
        algo_block="PIMPLE",
        div=dict(DIV_TRANSIENT, **{"div(phid,p)": "Gauss upwind",
                                   "div(phiv,p)": "Gauss upwind"}),
        algo={"momentumPredictor": "yes", "transonic": "yes",
              "nOuterCorrectors": "1", "nCorrectors": "2",
              "nNonOrthogonalCorrectors": "1",
              "pMinFactor": "0.1", "pMaxFactor": "2",
              "maxCo": "0.2", "rDeltaTSmoothingCoeff": "0.1",
              "rDeltaTDampingCoeff": "1", "maxDeltaT": "1"},
        fields_relax=None,        # None means: no relaxationFactors block at all
        equations_relax=None,
        solvers=SOLVERS_LTS,
    ),
}


# ---- ARM C2. LTS with A's SECOND-ORDER PRESSURE FLUX, added after arm C ran.
# Arm C as pre-registered converged, cleanly and monotonically, TO THE WRONG ANSWER: no
# supersonic pocket at all (min T 271 K against wolf's 238 K, i.e. peak local M 0.94 against
# the measured ~1.2), the leading-edge suction peak absent, C_l 0.443 against AGARD's 0.803,
# Cp RMS 0.492 against wolf's 0.076. The LOWER surface matched; only the upper was wrong.
# Eddy viscosity is not the cause -- max k is 551 against wolf's 1137, the same order.
#
# THE ONE ACCURACY-RELEVANT SETTING THAT DIFFERS FROM ARM A IS `div(phid,p) Gauss upwind`,
# FIRST ORDER on the dominant compressible pressure-flux term, taken from the plan's table
# via the shipped angledDuctLTS. Reporting "LTS is inaccurate" while that is what actually
# changed would blame the time integration for the discretisation, so this arm exists to
# separate them: identical to C except div(phid,p), which is A's `Gauss limitedLinear 1`.
# ARM C IS KEPT AND REPORTED. This does not replace it.
#
# ---- ARM D. ALLETTO'S ONERA M6, the only VALIDATED transonic external-aero rhoSimpleFoam
# setup we have found, and it is at Re ~1e7, between RAE's 6.5e6 and ARGUS's 1.68-2.11e7.
# Source: wiki.openfoam.com/OneraM6_by_Michael_Alletto (local copy in dso_reference/) and
# gitlab.com/mAlletto/openfoamtutorials/-/tree/master/OneraM6Wing. M 0.84, alpha 3.06 deg,
# validated against the NASA TMR experimental Cp at seven span stations.
#
# IT REFUTES THE PHASE 0 PREMISE, which is why it matters more than a fourth arm normally
# would. The plan's hypothesis was that ARGUS dies because `transonic yes` makes rho.relax()
# unreachable and the empty fields{} block leaves p undamped. Alletto runs `transonic yes`
# AND `consistent yes` AND `fields { p 1; }` AND `equations { p 1; ... }` -- that is p
# relaxation of ONE at both levels, i.e. LESS pressure damping than we have, no density
# damping either, at Re 1e7, converging to excellent agreement. Missing damping is therefore
# not sufficient to explain the ARGUS failure.
#
# WHAT ACTUALLY DIFFERS FROM OUR wolf RECIPE, all of it now in this arm:
#   relaxation   wolf runs p/U/e/h/k/omega ALL at 0.5 (0.3 in the ARGUS attempts) under
#                SIMPLEC. Alletto notes explicitly that SIMPLEC takes MUCH HIGHER factors
#                than SIMPLE, and uses p 1, U 0.9, e 0.8, k 0.9. We are 3x OVER-damped
#                relative to the reference, which is not a safety margin: it is off-recipe
#                for the algorithm and it makes convergence glacial without buying stability.
#   `bounded`    Alletto uses NONE, on a steadyState case. Ours is on all five divSchemes.
#   convection   linearUpwind limited throughout, not Minmod / limitedLinear. The wiki is
#                explicit that a LINEAR scheme on every convective term is what captures the
#                shock accurately; the detail is in his NACA0012 schemes-variation tutorial.
#   laplacian    `Gauss linear corrected` and snGrad `corrected`, against our `limited 0.5`.
#   nNonOrth     0, against our 1 (and 2 in the ARGUS builder).
#   solvers      GAMG for EVERYTHING including U and the turbulence, relTol 0.1,
#                nCellsInCoarsestLevel 20, against our GAMG-for-p plus PBiCGStab.
#   clamps       pMinFactor 0.1 / pMaxFactor 2, which is what pressure_clamp_admits()
#                independently derived from the physics. Our historical 0.5 was wrong.
#
# TWO DELIBERATE DEPARTURES FROM THE REFERENCE, both recorded rather than absorbed:
#   1. TURBULENCE MODEL. Alletto uses Spalart-Allmaras, chosen in his own model-variation
#      study for the best shock position and strength. Changing the turbulence default is
#      OFF-LIMITS without a decision entry (project rules), so this arm keeps kOmegaSST. That
#      makes it NOT a reproduction of ONERA M6, and any disagreement has this as a live
#      candidate cause.
#   2. ENERGY VARIABLE. Alletto solves sensibleInternalEnergy (e); our thermophysical model
#      selects sensibleEnthalpy (h). Both keys are written, since the relaxation and the
#      divSchemes are looked up by the SOLVED field's name and an `e` entry on an h case
#      silently relaxes nothing.
#
# RISK TO WATCH, stated in advance so it is not rationalised afterwards: `corrected` with
# nNonOrthogonalCorrectors 0 is an unlimited non-orthogonal correction with no corrector
# sweep, and our meshes are noMax 45.5 (RAE) and 38.3 (ARGUS) against a snappyHexMesh block
# mesh that is likely much better. If this arm goes unstable, that pairing is the first
# suspect, not the relaxation.
DIV_ALLETTO = {
    "default": "none",
    "div(phi,U)": "Gauss linearUpwind limited",
    "div(phi,k)": "Gauss linearUpwind limited",
    "div(phi,omega)": "Gauss linearUpwind limited",
    "div(phi,h)": "Gauss linearUpwind limited",
    "div(phi,e)": "Gauss linearUpwind limited",
    "div(phi,K)": "Gauss linearUpwind limited",
    "div(phi,Ekp)": "Gauss linearUpwind limited",
    "div(phid,p)": "Gauss limitedLinear 1",
    "div(phiv,p)": "Gauss upwind",
    "div((phi|interpolate(rho)),p)": "Gauss upwind",
    "div(((rho*nuEff)*dev2(T(grad(U)))))": "Gauss linear",
}

SOLVERS_ALLETTO = """solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-08;
        relTol          0.1;
        smoother        GaussSeidel;
        nCellsInCoarsestLevel 20;
    }
    "(U|e|h|k|omega|epsilon|nuTilda)"
    {
        solver          GAMG;
        tolerance       1e-08;
        relTol          0.1;
        smoother        GaussSeidel;
        nCellsInCoarsestLevel 20;
    }
}"""

RECIPES["alletto"] = dict(
    application="rhoSimpleFoam",
    ddt="steadyState",
    algo_block="SIMPLE",
    div=dict(DIV_ALLETTO),
    algo={"consistent": "yes", "transonic": "yes",
          "nNonOrthogonalCorrectors": "0", "pMinFactor": "0.1", "pMaxFactor": "2"},
    fields_relax={"p": "1"},
    equations_relax={"p": "1", "U": "0.9", "e": "0.8", "h": "0.8",
                     '"(k|omega)"': "0.9"},
    residual_control={"p": "1e-5", "U": "1e-6", '"(k|omega|e|h)"': "1e-5"},
    solvers=SOLVERS_ALLETTO,
    grad_schemes="""gradSchemes
{
    default         Gauss linear;
    limited         cellLimited Gauss linear 1;
    grad(U)         $limited;
    grad(k)         $limited;
    grad(omega)     $limited;
}""",
    laplacian="Gauss linear corrected",
    sngrad="corrected",
)
RECIPES["lts_ll"] = dict(
    RECIPES["lts"],
    div=dict(DIV_TRANSIENT, **{"div(phid,p)": "Gauss limitedLinear 1",
                               "div(phiv,p)": "Gauss upwind"}),
)


# ---- ARM E. WOLF'S SCHEMES WITH ALLETTO'S ITERATION (objective, 2026-08-11: the
# accuracy of wolf with the stability of alletto). The two are separable, and the evidence says which
# side of the recipe owns which property.
#
# WHAT MAKES ALLETTO STABLE is the ITERATION, not the discretisation: SIMPLEC-appropriate
# relaxation (p 1, U 0.9, e 0.8, k 0.9 against wolf's uniform 0.3), GAMG on every variable,
# and a residual stop. Measured: alletto reached C_d 0.0175 / C_l 0.803 by iteration 600
# where wolf needed ~7,000, and its tail drift is 0.07 counts against wolf's 1.12.
#
# WHAT MAKES ALLETTO UNPHYSICAL is the DISCRETISATION, and the evidence is threefold:
#   1. T_max reaches 400 K against an adiabatic stagnation ceiling of 318.9 K -- 81 K of
#      ENERGY CREATED, sustained, above T0 on 51.8% of iterations and pinned at the limiter
#      on 14.0%. An adiabatic flow cannot do that; only a scheme can.
#   2. A standing wave train on the upper surface at x/c 0.28 / 0.39 / 0.52, about 13 cells
#      per wavelength, with matching spikes in C_f and finger striations in the field. RAE
#      2822 at M 0.730 has ONE shock; three compressions are numerical.
#   3. `Gauss linearUpwind limited` is a second-order upwind RECONSTRUCTION whose GRADIENT is
#      limited; the face interpolation itself is not TVD, so it overshoots across a shock.
#      wolf uses Minmod on momentum and limitedLinear on energy, both of which are.
#      Alletto's own mesh hides this: snappyHexMesh surface refinement level 8-9 around the
#      wing, against our uniform ds/c 0.00875 aft of the shoulder.
# The third point is INFERENCE from the scheme definitions; the first two are measurements.
#
# SO ARM E TAKES THE ITERATION AND LEAVES THE DISCRETISATION: every scheme is wolf's, every
# relaxation factor, linear solver and stopping criterion is Alletto's. If the diagnosis is
# right this is stable AND accurate. If it is stable and still wiggles, the cause is the
# `corrected` laplacian/snGrad or nNonOrthogonalCorrectors 0 rather than the convection, and
# arm F separates that.
RECIPES["hybrid"] = dict(
    application="rhoSimpleFoam",
    ddt="steadyState",
    algo_block="SIMPLE",
    div=dict(DIV_COMMON, **{"div(phid,p)": "Gauss limitedLinear 1"}),
    algo={"consistent": "yes", "transonic": "yes",
          "nNonOrthogonalCorrectors": "1", "pMinFactor": "0.1", "pMaxFactor": "2"},
    fields_relax={"p": "1"},
    equations_relax={"p": "1", "U": "0.9", "e": "0.8", "h": "0.8", '"(k|omega)"': "0.9"},
    residual_control={"p": "1e-5", "U": "1e-6", '"(k|omega|e|h)"': "1e-5"},
    solvers=SOLVERS_ALLETTO,
    laplacian="Gauss linear limited 0.5",
    sngrad="limited 0.5",
)

# ---- ARM F. ALLETTO WITH TVD CONVECTION, to attribute the wiggle rather than merely avoid
# it. Identical to `alletto` except the convective schemes become wolf's limited ones. If the
# wave train and the 81 K of created energy disappear, the convection owned them and the
# `corrected` laplacian on a noMax 45.5 mesh is exonerated; if they survive, it did not.
# Worth one run because the same choice has to be made again for the 3D unstructured mesh,
# where the reference case's own refinement will not be there to hide it either.
RECIPES["alletto_tvd"] = dict(
    RECIPES["alletto"],
    div=dict(DIV_ALLETTO, **{"div(phi,U)": "bounded Gauss Minmod",
                             "div(phi,h)": "bounded Gauss limitedLinear 1",
                             "div(phi,e)": "bounded Gauss limitedLinear 1",
                             "div(phi,K)": "bounded Gauss limitedLinear 1",
                             "div(phi,Ekp)": "bounded Gauss limitedLinear 1",
                             "div(phi,k)": "bounded Gauss linearUpwind limitedT",
                             "div(phi,omega)": "bounded Gauss linearUpwind limitedT"}),
)


# ---- ARM G. THE COMPLEMENT OF ARM F, and it exists because F DIED rather than answered.
# F took Alletto and swapped in wolf's `bounded` TVD convection, holding `corrected`
# laplacian/snGrad and nNonOrthogonalCorrectors 0. It crashed at ITERATION 110 with a
# floating point exception inside GAMG's coarsest-level solve, p pinned at both clamps and T
# at both limiter stops. So the two halves of Alletto's recipe are NOT independently
# swappable, and F cannot attribute the wiggle because it never produced a solution.
#
# G swaps the OTHER half instead: Alletto's convection is kept exactly as published, and only
# the laplacian, snGrad and non-orthogonal corrector become wolf's. The attribution is then
# clean in the surviving direction:
#   G still wiggles  -> the unlimited `linearUpwind` CONVECTION owns it, as inferred
#   G is clean       -> `corrected` with nNonOrth 0 on a noMax 45.5 mesh owned it, and the
#                       inference from the scheme definitions was wrong
# Either outcome is worth having before the 3D unstructured mesh, where non-orthogonality
# will be worse than 45.5 and the reference case's surface refinement will not be there.
RECIPES["alletto_lim"] = dict(
    RECIPES["alletto"],
    algo=dict(RECIPES["alletto"]["algo"], nNonOrthogonalCorrectors="1"),
    laplacian="Gauss linear limited 0.5",
    sngrad="limited 0.5",
)


# ---- ARM H. THE HYBRID WITH KRYLOV MOMENTUM, after arm E died on its LINEAR SOLVER.
#
# WHAT KILLED ARM E, read off its own log rather than inferred:
#     GAMG:  Solving for Ux, Initial residual = 5.9e-04, Final residual = 7.0e+101,
#            No Iterations 1000
# The GAMG solve for the MOMENTUM equation diverged by a hundred orders of magnitude and hit
# its iteration cap. The flow was healthy at the time -- min p 56,791 Pa, T_max 318.68 against
# T0 318.9, zero energy-residual events at iteration 19,213 -- so this is not the physics
# failing, it is the LINEAR ALGEBRA failing on a perfectly good matrix-free-of-trouble state.
# The same thing killed arm E on ARGUS at iteration 14.
#
# WHY, and this is the transferable part. GAMG IS A MULTIGRID METHOD AND MULTIGRID ASSUMES THE
# OPERATOR IS ELLIPTIC-DOMINATED: its coarse-grid correction is only a correction if smooth
# error components on the fine grid remain smooth on the coarse one. The PRESSURE equation is
# elliptic (or hyperbolic-elliptic under `transonic`) and GAMG is the right tool for it. THE
# MOMENTUM EQUATION AT HIGH CELL REYNOLDS NUMBER IS CONVECTION-DOMINATED and strongly
# non-symmetric, and on this mesh it is also viciously anisotropic: the near-wall cells run an
# aspect ratio of order 1e4-1e5. Agglomeration across cells that anisotropic produces coarse
# operators that no longer approximate the fine one, and the V-cycle amplifies instead of
# damping. PBiCGStab with DILU makes no such assumption: it is a Krylov method with an
# incomplete-LU preconditioner, which is what non-symmetric convection-dominated systems want.
#
# WHY ALLETTO GETS AWAY WITH GAMG-ON-EVERYTHING AND WE DO NOT. Two differences, both ours:
# his mesh is a snappyHexMesh 3D block-structured-ish grid whose cells are near-cubic away
# from the five prism layers, and his convection is unbounded linearUpwind, which produces a
# different and more diagonally dominant momentum matrix than wolf's `bounded Gauss Minmod`.
# Arms `alletto` and `alletto_lim` both ran 60,000 iterations with GAMG-everything on this
# same mesh; arm E, which differs from them in the convective schemes, did not. So it is the
# PAIRING of bounded TVD convection with GAMG momentum that fails, not either alone.
#
# SO ARM H TAKES ALLETTO'S RELAXATION AND LEAVES HIS LINEAR SOLVERS. That is the correct
# separation anyway: the relaxation fix is about SIMPLEC being under-relaxed by a factor of
# three, and the linear-solver choice is orthogonal to it. GAMG for p, Krylov for everything
# else, which is both the OpenFOAM convention and what the wolf recipe already had right.
# AND THE LINEAR TOLERANCES COME FROM ALLETTO TOO, which the first attempt got wrong.
# Built with wolf's solvers block verbatim (p relTol 0.001, everything else relTol 0 to
# tolerance 1e-8) it ran 6.5x SLOWER than the other arms on the same mesh: 216 GAMG iterations
# for the pressure solve against alletto_lim's 2, i.e. 1.0 outer iterations per second against
# 6.7.
#
# THE REASON IS THE PAIRING, NOT EITHER SETTING. In a SIMPLE-type scheme the linear system is
# ITSELF AN APPROXIMATION -- the operator is lagged and the right-hand side is about to change
# -- so solving it far tighter than the outer iteration's own error is work thrown away. wolf's
# relTol 0.001 was calibrated against relaxation 0.3, where the pressure correction per outer
# iteration is small. Alletto's p relaxation of 1 makes that correction LARGE, and demanding
# three orders of magnitude on a large correction is the expensive way to compute something
# that is about to be superseded. Tight linear tolerances and aggressive relaxation are the
# worst combination: you pay for exactness in a system whose input you are about to replace.
# relTol 0.1 per outer iteration is the standard discipline and is what Alletto uses.
SOLVERS_HYBRID = """solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-08;
        relTol          0.1;
        smoother        GaussSeidel;
        nCellsInCoarsestLevel 20;
    }
    "(U|e|h|k|omega)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-08;
        relTol          0.1;
    }
}"""

RECIPES["hybrid_krylov"] = dict(
    RECIPES["hybrid"],
    solvers=SOLVERS_HYBRID,
)


# ---- ARMS I AND J. SEPARATING SPEED FROM CONVERGENCE LEVEL.
# hybrid_krylov runs everywhere -- RAE, ARGUS normal, ARGUS freestream -- and is the only arm
# that does. But on RAE it lands at shock x/c 0.481 against a measured 0.550 and C_l -8.6%,
# where wolf, WITH IDENTICAL DISCRETISATION, gets 0.551 and -1.0%. A converged steady solution
# cannot depend on the iteration path, so one of them is not converged, and the logs say which:
#     wolf           outer p residual  1.2e-05
#     hybrid_krylov  outer p residual  2.2e-03      -- 190x higher
#
# TWO CANDIDATE CAUSES, both introduced by me and both testable in one run each:
#   I  relTol 0.1 on the linear solves. It bought a 6.2x speedup, but a linear solve that only
#      reduces its residual tenfold per outer iteration limits how far the OUTER loop can
#      settle. `hybrid_tight` sets relTol 0.01, between wolf's 0.001 and Alletto's 0.1.
#   J  p relaxation of 1. SIMPLEC needs no under-relaxation to be STABLE, but the outer
#      iteration is still a fixed-point map and a factor below 1 contracts it. wolf runs 0.3.
#      `hybrid_relax` keeps relTol 0.1 and drops p to 0.7, U to 0.7.
# One variable each, so whichever recovers the shock position identifies the mechanism. If
# BOTH matter the study says so rather than picking the first that works.
SOLVERS_HYBRID_TIGHT = SOLVERS_HYBRID.replace("relTol          0.1;", "relTol          0.01;")

RECIPES["hybrid_tight"] = dict(RECIPES["hybrid_krylov"], solvers=SOLVERS_HYBRID_TIGHT)
RECIPES["hybrid_relax"] = dict(
    RECIPES["hybrid_krylov"],
    fields_relax={"p": "0.7"},
    equations_relax={"p": "0.7", "U": "0.7", "e": "0.8", "h": "0.8", '"(k|omega)"': "0.9"},
)


# ---- ARM K. TVD LIMITING WITHOUT THE BOUNDED SOURCE TERM. What the theory actually points at.
#
# THE THREE-POINT ATTRIBUTION, all on RAE case 9, all at the same mesh and cold start:
#     bounded Minmod   + relax 0.3       shock error 0.001 c   (wolf)
#     bounded Minmod   + relax 1.0/0.9   shock error 0.069 c   (hybrid_krylov)
#     unbnd linearUpwd + relax 1.0/0.9   shock error 0.008 c   (alletto_lim)
# The first and third are fine and the middle one is not, so it is neither the relaxation
# alone nor the convection alone: IT IS THE PAIRING. And it is not the linear tolerance --
# hybrid_krylov and alletto_lim run the same relTol 0.1 and reach the same outer residual,
# 5.3e-4 against 5.4e-4, and differ by 9x in shock error.
#
# MECHANISM. `bounded Gauss X` subtracts fvm::Sp(fvc::div(phi), psi). That term is zero only
# at convergence; during the solve it is a SOLUTION-DEPENDENT SOURCE proportional to the
# continuity error. Relaxation of 0.3 damps the per-iteration field change, so div(phi) stays
# small and the source is negligible. At p relaxation 1 the field moves hard every iteration,
# div(phi) stays large, and the source persists -- biasing the shock forward by 0.07 c.
#
# SO TAKE THE TVD LIMITING AND LEAVE THE SOURCE. Minmod and limitedLinear are flux limiters
# and are what capture a shock without the wiggle train that unlimited linearUpwind produced
# on `alletto` (81 K of created energy, three spurious compressions). `bounded` is a
# convergence-acceleration device for steady solves and is separable from them. This arm is
# hybrid_krylov with the `bounded` prefix removed and NOTHING else changed, so if the shock
# comes back the attribution is complete.
#
# IT ALSO POINTS THE RIGHT WAY FOR 3D. `bounded` is worth least exactly where continuity is
# hardest to satisfy, which is an unstructured mesh with high non-orthogonality.
RECIPES["hybrid_unbounded"] = dict(
    RECIPES["hybrid_krylov"],
    div=dict(DIV_TRANSIENT, **{"div(phid,p)": "Gauss limitedLinear 1"}),
)


# ---- ARM L. WOLF'S ACCURACY WITH THE ONE THING THAT MIGHT BE KILLING IT CHANGED.
# Neither relTol (hybrid_tight, C_l 0.741) nor `bounded` (hybrid_unbounded, 0.741) recovers
# the lift that hybrid_krylov loses, so the culprit in that family is the RELAXATION acting on
# wolf's LIMITERS: Minmod and cellMDLimited are nonlinear switches whose activation depends on
# the local state, and at p relaxation 1 the field moves so far per iteration that the limiter
# chatters and locks the shock into a smeared forward position. wolf's 0.3 damps that; 0.7
# (hybrid_relax) diverged outright, so the window is narrow.
#
# THAT MEANS WOLF'S RELAXATION MUST BE KEPT, and the question becomes what ELSE kills wolf at
# Re 1.68e7. Its energy LINEAR solves are healthy there -- final residual 1e-9 to 1e-11 in
# 2-15 iterations -- while the INITIAL residual pegs at 1.0, so the algebra is fine and the
# outer loop is diverging. hybrid_krylov, with the SAME SCHEMES, survives ARGUS. The only
# things left between them are the linear-solver setup and residualControl.
# So: wolf, unchanged, except GAMG-for-p / Krylov-for-the-rest at relTol 0.1 and a residual
# stop. If it survives ARGUS this is the toolbox recipe -- wolf accuracy, hybrid stability.
RECIPES["wolf_gamg"] = dict(
    RECIPES["wolf"],
    solvers=SOLVERS_HYBRID,
    residual_control={"p": "1e-6", "U": "1e-6", '"(k|omega|e|h)"': "1e-6"},
)


# ---- ARM M. TVD WITH A *SMOOTH* LIMITER, WHICH IS WHAT THE EVIDENCE NOW POINTS AT.
# Three failed hypotheses narrowed it to one. relTol did not fix hybrid (0.741 at 0.01 vs
# 0.734 at 0.1). `bounded` did not fix it (0.741 unbounded). Relaxation 0.7 diverged. What
# survives is LIMITER CHATTER: Minmod and cellMDLimited are nonlinear SWITCHES whose
# activation depends on the local state, and at p relaxation 1 the field moves far enough per
# iteration that they toggle, stalling the shock in a smeared forward position. wolf's
# relaxation of 0.3 damps the field change enough that they settle; Alletto avoids the problem
# by using no flux limiter at all, and pays for it with a wiggle train and 81 K of created
# energy on the arm that used `corrected`.
#
# MINMOD IS THE MOST AGGRESSIVE OF THE COMMON LIMITERS -- its switching function is the
# steepest, so it chatters first. limitedLinear and vanLeer are smooth over most of their
# range while still being TVD, so they bound the shock WITHOUT the on/off behaviour that
# interacts with the relaxation. That is the combination nothing has tested yet: TVD (so no
# wiggles, unlike alletto_lim) with a smooth limiter (so no chatter, unlike hybrid_*), driven
# by Alletto's iteration (so it survives Re 1.68e7, unlike wolf).
#
# PRECEDENT: OpenFOAM-7's own transonic `rhoPimpleFoam/RAS/nacaAirfoil` uses
# `div(phi,U) Gauss limitedLinearV 1` with `transonic yes` -- the V form being the
# vector-aware limiter, which limits on the direction of steepest change rather than
# component-wise and is the right choice for momentum.
RECIPES["alletto_smooth"] = dict(
    RECIPES["alletto_lim"],
    div=dict(RECIPES["alletto_lim"]["div"], **{
        "div(phi,U)": "Gauss limitedLinearV 1",
        "div(phi,h)": "Gauss limitedLinear 1",
        "div(phi,e)": "Gauss limitedLinear 1",
        "div(phi,K)": "Gauss limitedLinear 1",
        "div(phi,Ekp)": "Gauss limitedLinear 1",
        "div(phi,k)": "Gauss limitedLinear 1",
        "div(phi,omega)": "Gauss limitedLinear 1"}),
)


# ---- ARMS N AND O. HOW MUCH UNDER-RELAXATION DOES THE ACCURATE PACKAGE NEED AT HIGH Re?
#
# The evidence now says a FLUX LIMITER AND NEAR-UNITY RELAXATION ARE INCOMPATIBLE. Four arms,
# three different limiters -- bounded Minmod, unbounded Minmod, limitedLinearV -- all land at
# C_l 0.741 +/- 0.002 against a measured 0.803, while the two arms that are accurate are the
# ones that either under-relax (wolf, 0.794 at relaxation 0.3) or use no flux limiter at all
# (alletto_lim, 0.796 at relaxation 1.0). A limiter is a NON-DIFFERENTIABLE operator: the
# SIMPLE iteration lands on a limit cycle of the limiter switching, and the average of that
# cycle solves neither branch. The standard cures are under-relaxation, freezing the limiter,
# or a smooth limiter, and the third has now been tried and crashed.
#
# SO THE QUESTION IS THE AMOUNT. wolf needs 0.3 at Re 6.5e6 and diverges at 1.68e7, where the
# physical viscous damping is 2.6x weaker; hybrid_relax at 0.7 diverged even at 6.5e6. If the
# required factor simply scales with the weakening of the physical damping, more of it should
# recover the accurate package at ARGUS conditions, and nothing else has to change.
#
# N: everything down to 0.15, uniformly, which is the direct test of that scaling.
# O: only the ENERGY equation damped, to 0.1, because h is the equation that actually
#    diverges -- its INITIAL residual pegs at 1.0 on ARGUS while its linear solves stay
#    healthy at 1e-9. If O works and N is unnecessary, the cost is one equation's convergence
#    rate rather than the whole system's.
RECIPES["wolf_damped"] = dict(
    RECIPES["wolf"],
    # 0.2, NOT THE 0.15 I FIRST WROTE. The pre-run audit refused 0.15: its SIMPLEC band for a
    # uniform factor is 0.2 to 0.7. That floor is a heuristic rather than a derived bound, so
    # I could have widened it -- and widening a gate to admit my own experiment is precisely
    # the antipattern this project has a standing rule against. 0.2 sits ON the floor and is
    # still 1.5x more damping than wolf's 0.3, which is enough to test the scaling. If 0.2 is
    # insufficient AND the trend points lower, THAT is the evidence that justifies moving the
    # floor, with the reason recorded.
    equations_relax={"p": "0.2", "U": "0.2", "e": "0.2", "h": "0.2",
                     '"(k|omega)"': "0.2"},
)
RECIPES["wolf_hdamp"] = dict(
    RECIPES["wolf"],
    equations_relax={"p": "0.3", "U": "0.3", "e": "0.1", "h": "0.1",
                     '"(k|omega)"': "0.3"},
)


# ---- ARM P. PLAIN SIMPLE WITH FLUX LIMITERS: the internally consistent pairing.
#
# THE CONFLICT THE STUDY HAS UNCOVERED. Two requirements that cannot both be met under SIMPLEC:
#   1. SIMPLEC AT HIGH Re WANTS p RELAXATION NEAR 1. Its velocity correction already carries
#      the neighbour contribution (rAtU = 1/(1/rAU - UEqn.H1())), so the derivation assumes the
#      full correction is applied; under-relaxing p hard puts momentum and continuity out of
#      step. Measured: wolf (SIMPLEC, p 0.3) DIES at Re 1.68e7 at iteration 1,958, while
#      hybrid_krylov -- IDENTICAL SCHEMES, p 1.0 -- survives 60,000 at both ARGUS conditions.
#      Less damping is more stable, which only makes sense if the damping is off-algorithm.
#   2. FLUX LIMITERS WANT LOW RELAXATION. Four arms with three limiters at p 1.0 all land at
#      C_l 0.741 +/- 0.002 against 0.803: a non-differentiable operator driven hard enough to
#      chatter, whose limit-cycle average solves neither branch.
#
# EVERY HYBRID FAILED BECAUSE IT TRIED TO SATISFY BOTH. The resolution is to change ALGORITHM
# rather than to keep tuning the factor: PLAIN SIMPLE genuinely REQUIRES under-relaxation --
# its correction neglects the neighbour terms, so p ~ 0.3 / U ~ 0.7 is not a safety margin but
# the algorithm's actual operating point. That is the same low relaxation the limiters need,
# so the two requirements become one instead of fighting.
#
# WHY THIS IS NOT THE `tutorial` ARM, which was also plain SIMPLE and failed. That arm ran
# `transonic OFF` with fields{p 0.7, rho 0.01} and had C_d,pressure negative on 46% of its
# iterations. This one keeps `transonic yes` (the pressure equation stays hyperbolic, which is
# what a shock needs), keeps wolf's validated TVD schemes, and uses SIMPLE's own textbook
# factors. Only `consistent` and the relaxation move.
RECIPES["simple_tvd"] = dict(
    RECIPES["wolf"],
    algo=dict(RECIPES["wolf"]["algo"], consistent="no"),
    # CORRECTED. My first version put p 0.3 in equations{}, which under plain SIMPLE relaxes
    # the pressure MATRIX (diagonal dominance) and leaves the pressure FIELD taking the full
    # correction every iteration -- i.e. under-damped SIMPLE, not the textbook pairing at all.
    # It diverged within 500 iterations on BOTH cases, including RAE, which every other arm
    # survives. The classic SIMPLE pairing is a FIELD relaxation on p of ~0.3 with an EQUATION
    # relaxation on U of ~0.7, and p.relax() is what applies it.
    fields_relax={"p": "0.3"},
    equations_relax={"U": "0.7", "e": "0.7", "h": "0.7", '"(k|omega)"': "0.7"},
    solvers=SOLVERS_HYBRID,
    residual_control={"p": "1e-6", "U": "1e-6", '"(k|omega|e|h)"': "1e-6"},
)

# ---- ARM Q. ALLETTO'S PACKAGE WITH THE LIMITER ONLY WHERE THE DEFECT IS.
# alletto_lim is the accuracy leader that survives -- shock 0.008 c, C_l -0.8%, clean at RAE
# and at ARGUS freestream. Its two faults are the wiggle train and T above the adiabatic
# ceiling, and BOTH ARE ENERGY SYMPTOMS: an unlimited convective scheme overshooting across
# the shock puts heat where there is none. alletto_smooth limited EVERYTHING and fell into the
# C_l 0.741 chatter basin with the rest.
# So limit ONLY the energy transport and leave momentum unlimited. If the chatter comes from
# the MOMENTUM limiter interacting with the relaxation -- which is what the 0.741 cluster
# suggests, since every member of it limited div(phi,U) -- then bounding h and K alone should
# remove the created energy without moving the lift.
RECIPES["alletto_elim"] = dict(
    RECIPES["alletto_lim"],
    div=dict(RECIPES["alletto_lim"]["div"], **{
        "div(phi,h)": "Gauss limitedLinear 1",
        "div(phi,e)": "Gauss limitedLinear 1",
        "div(phi,K)": "Gauss limitedLinear 1",
        "div(phi,Ekp)": "Gauss limitedLinear 1"}),
)


def to_spalart_allmaras(case, nu_inf=None, verbose=True):
    """Convert a kOmegaSST case to Spalart-Allmaras. Project decision, 2026-08-11.

    WALL FUNCTIONS AND TURBULENCE-MODEL DEFAULTS ARE BOTH OFF-LIMITS WITHOUT A DECISION ENTRY
    (project rules). It was authorised explicitly by project decision, accepting that every
    geometry is then re-run, and it needs ARG-161 to carry it before anything from an SA run is
    quotable.

    THE MOTIVATION IS MEASURED, NOT PREFERENCE. Every arm on RAE clips k or omega tens of
    thousands of times: wolf 45,001, alletto_elim 60,001, alletto_lim 62,527 bounding events in
    60,000 iterations. omega goes like 6*nu/(beta*y^2) at a wall, so the omega equation is
    STIFF wherever the first cell is small, and it is the equation being clipped. SA HAS NO
    OMEGA EQUATION: one transport equation for nuTilda, no near-wall singularity, and it is the
    standard choice for external aerodynamics at high Reynolds number for exactly this reason.
    Alletto's own model-variation study chose it over kOmegaSST for shock position and strength.

    IT IS A DIFFERENT MODEL, NOT A TUNING KNOB. Everything computed with it must be re-run --
    RAE for the validation claim and every ARGUS geometry -- because a delta between a kOmega
    baseline and an SA morph would be a turbulence-model difference wearing a morphing label.

    Written from the shipped incompressible/simpleFoam/airFoil2D case, which is the direct
    analogue: 2D aerofoil, wall-modelled, nutUSpaldingWallFunction, nuTilda fixedValue 0 at the
    wall and freestream elsewhere.
    """
    case = pathlib.Path(case)
    tp = case / "constant/turbulenceProperties"
    s = tp.read_text()
    s = re.sub(r"RASModel\s+\w+;", "RASModel        SpalartAllmaras;", s)
    tp.write_text(s)
    got = re.search(r"RASModel\s+(\w+);", _expand(case, "../constant/turbulenceProperties"))
    if not got or got.group(1) != "SpalartAllmaras":
        raise SystemExit("RASModel did not land: %s" % (got.group(1) if got else None))

    # nuTilda freestream value. THE CONVENTIONAL SEED IS 3-5x THE LAMINAR nu, not the 0.14 of
    # the shipped tutorial -- that case is at nu 1e-5 and a different scale entirely, and a
    # dimensional constant copied across scales is the defect this project keeps finding.
    # Derived from THIS case's own molecular viscosity.
    if nu_inf is None:
        sys.path.insert(0, str(REPO / "scripts"))
        import case_normalisation
        nu_inf = case_normalisation.derive(case)["flow"]["nu"]
    nt = 4.0 * nu_inf

    b = (case / "constant/polyMesh/boundary").read_text(errors="replace")
    patches = re.findall(r"\n    (\w+)\s*\n    \{[^}]*?type\s+(\w+);", b, re.S)
    bf = ""
    for name, ptype in patches:
        if ptype == "empty":
            bf += "    %s\n    {\n        type            empty;\n    }\n" % name
        elif ptype == "wall":
            bf += ("    %s\n    {\n        type            fixedValue;\n"
                   "        value           uniform 0;\n    }\n" % name)
        else:
            bf += ("    %s\n    {\n        type            freestream;\n"
                   "        freestreamValue uniform %.6e;\n    }\n" % (name, nt))
    (case / "0/nuTilda").write_text(
        _header("nuTilda").replace('location    "system";', 'location    "0";')
        .replace("class       dictionary;", "class       volScalarField;")
        + "dimensions      [0 2 -1 0 0 0 0];\ninternalField   uniform %.6e;\n\n"
          "boundaryField\n{\n%s}\n" % (nt, bf))
    for f in ("0/k", "0/omega"):
        if (case / f).exists():
            (case / f).unlink()
    # ---- THE RESIDUALS FUNCTION OBJECT MUST FOLLOW THE MODEL (project decision, 2026-08-11).
    # system/residuals still asked for `fields (p U h k omega)`. k and omega no longer exist,
    # so the function object SKIPS THEM SILENTLY and nuTilda -- the only turbulence equation
    # there is now -- was not being logged at all. That is exactly the failure this repo already
    # documents for the energy residual: `e` requested on a case solving `h`, silently dropped,
    # and every residual plot missing the one equation that mattered. Same trap, my turn.
    rq = case / "system/residuals"
    if rq.exists():
        s2 = re.sub(r"fields\s*\([^)]*\);", "fields (p U h nuTilda);", rq.read_text())
        rq.write_text(s2)
        if "fields (p U h nuTilda);" not in rq.read_text():
            raise SystemExit("residuals field list did not land for SA")
    if verbose:
        print("  SA: nuTilda freestream %.4e (4x nu = %.4e), k and omega removed" % (nt, nu_inf))
    return nt


# ---- SA VARIANTS. Same numerics as their kOmegaSST parents; only the turbulence field
# changes, so any difference measured is the MODEL and not the discretisation. div(phi,nuTilda)
# takes the scheme its parent used for k and omega, and nuTilda joins the solver and relaxation
# groups in the same place they did.
def _sa_of(name, parent):
    d = dict(RECIPES[parent])
    dv = {k: v for k, v in d["div"].items() if "phi,k" not in k and "phi,omega" not in k}
    turb = RECIPES[parent]["div"].get("div(phi,k)", "Gauss linearUpwind limited")
    dv["div(phi,nuTilda)"] = turb
    d["div"] = dv
    if d.get("equations_relax"):
        er = {k: v for k, v in d["equations_relax"].items() if "k|omega" not in k}
        er["nuTilda"] = RECIPES[parent]["equations_relax"].get('"(k|omega)"', "0.9")
        d["equations_relax"] = er
    if d.get("residual_control"):
        d["residual_control"] = {"p": "1e-5", "U": "1e-6", '"(nuTilda|e|h)"': "1e-5"}
    d["solvers"] = d["solvers"].replace("k|omega|epsilon|nuTilda", "nuTilda").replace(
        "(U|e|h|k|omega)", "(U|e|h|nuTilda)").replace(
        '"(k|omega)"', '"(nuTilda)"')
    RECIPES[name] = d

_sa_of("sa_elim", "alletto_elim")
_sa_of("sa_lim", "alletto_lim")
_sa_of("sa_wolf", "wolf")
# THE PURE ONERA COMBINATION (review question, 2026-08-11: the ONERA settings as published). Alletto
# tuned his numerics WITH Spalart-Allmaras -- corrected laplacian, nNonOrth 0, linearUpwind
# throughout, p relaxation 1 -- and every arm so far has broken that pairing in one direction
# or the other: `alletto` ran his numerics with kOmegaSST, `sa_elim` runs SA with MY laplacian
# and an energy limiter he does not have. Neither is his case. This one is.
_sa_of("sa_alletto", "alletto")


# ---- WHAT EACH RECIPE IS, IN THE AXES THAT ACTUALLY VARY (project decision, 2026-08-11:
# names must say what a recipe does). The IDs encode ANCESTRY -- alletto, wolf,
# hybrid -- so reading one tells you where it came from and not what it does, and the worst
# pair is `alletto_lim` vs `alletto_elim`: ONE LETTER apart, where `lim` is a limited LAPLACIAN
# and `elim` is an ENERGY FLUX LIMITER. That single entry is worth 0.8 points of C_l and the
# difference between a clean field and a wiggle train, and it is the basis of the whole
# two-recipe split. A distinction that load-bearing must not hang on one letter.
# The IDs are frozen (renaming would orphan every case and figure produced today); THIS is the
# authority on what they mean, and docs/naming_recipe_trial.md carries the same table.
RECIPE_AXES = {
    "wolf":            dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="low 0.3", solvers="krylov"),
    "alletto":         dict(algo="SIMPLEC", mom="lup", energy="lup",
                            lap="CORRECTED", nno=0, relax="high", solvers="gamg-all"),
    "alletto_lim":     dict(algo="SIMPLEC", mom="lup", energy="LUP",
                            lap="lim05", nno=1, relax="high", solvers="gamg-all"),
    "alletto_elim":    dict(algo="SIMPLEC", mom="lup", energy="TVD",
                            lap="lim05", nno=1, relax="high", solvers="gamg-all"),
    "alletto_smooth":  dict(algo="SIMPLEC", mom="tvd", energy="tvd",
                            lap="lim05", nno=1, relax="high", solvers="gamg-all"),
    "alletto_tvd":     dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="CORRECTED", nno=0, relax="high", solvers="gamg-all"),
    "hybrid":          dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="high", solvers="gamg-all"),
    "hybrid_krylov":   dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="high", solvers="krylov"),
    "hybrid_unbounded":dict(algo="SIMPLEC", mom="tvd", energy="tvd",
                            lap="lim05", nno=1, relax="high", solvers="krylov"),
    "hybrid_tight":    dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="high", solvers="krylov relTol 0.01"),
    "hybrid_relax":    dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="0.7", solvers="krylov"),
    "wolf_damped":     dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="0.2", solvers="krylov"),
    "wolf_hdamp":      dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="low, h 0.1", solvers="krylov"),
    "wolf_gamg":       dict(algo="SIMPLEC", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="low 0.3", solvers="krylov relTol 0.1"),
    "simple_tvd":      dict(algo="SIMPLE", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="p0.3 field", solvers="krylov"),
    "tutorial":        dict(algo="SIMPLE", mom="tvd+bounded", energy="tvd+bounded",
                            lap="lim05", nno=1, relax="p0.7 rho0.01", solvers="krylov"),
    "lts":             dict(algo="LTS", mom="lup", energy="lup, 1st-order p",
                            lap="lim05", nno=1, relax="none", solvers="krylov"),
    "lts_ll":          dict(algo="LTS", mom="lup", energy="lup, 2nd-order p",
                            lap="lim05", nno=1, relax="none", solvers="krylov"),
}
for _n in list(RECIPE_AXES):
    RECIPE_AXES["sa_" + _n] = dict(RECIPE_AXES[_n], model="SpalartAllmaras")
# AND THE TWO I NAMED INCONSISTENTLY. The derived keys are `sa_alletto_elim` and
# `sa_alletto_lim`, but I registered the recipes as `sa_elim` and `sa_lim` -- dropping the
# parent, so they no longer say which recipe they are the SA version of. Caught by this table's
# own completeness check reporting them as undocumented, which is the check earning its keep.
# Aliased rather than renamed: sa_elim already has two case directories and figures.
RECIPE_AXES["sa_elim"] = dict(RECIPE_AXES["sa_alletto_elim"])
RECIPE_AXES["sa_lim"] = dict(RECIPE_AXES["sa_alletto_lim"])


def describe(name):
    """One line saying what a recipe IS, for printing beside any result."""
    a = RECIPE_AXES.get(name)
    if not a:
        return "%s (undocumented -- add it to RECIPE_AXES)" % name
    return ("%-17s %-8s mom %-12s energy %-16s lap %-10s nNonOrth %s  relax %-12s %s"
            % (name, a["algo"], a["mom"], a["energy"], a["lap"], a["nno"], a["relax"],
               a["solvers"]))



# ---- THE ENERGY SCHEME IS A DIAL, NOT A SWITCH (review question, 2026-08-11: is
# there a scheme between the two that counters this known problem).
#
# `alletto_lim` and `alletto_elim` have been treated all day as two options. They are the
# ENDPOINTS OF ONE FAMILY. OpenFOAM-7's limitedLinear limiter is
#     max(min((2/k)*r, 1), 0)
# so k -> 0 saturates the limiter at 1 and gives PURE LINEAR (= alletto_lim), and k = 1 gives
# full TVD limiting (= alletto_elim). Anything between limits only where the gradient ratio is
# genuinely adverse. If the middle works, the shock rule becomes a SETTING chosen by shock
# strength rather than a choice between two recipes -- one pipeline instead of two, which
# removes the "validated recipe is not the production recipe" objection entirely.
#
# AND ONE SCHEME BUILT FOR EXACTLY THIS DEFECT. filteredLinear2's own header: "remove
# high-frequency modes with staggering characteristics by comparing the face gradient with
# both neighbouring cell gradients and introduce small amounts of upwind in order to damp
# these modes." That is a literal description of the finger structure in the supersonic
# pocket. Unlike a TVD limiter it does not clamp wherever r looks bad; it DETECTS the
# staggering and adds only the dissipation needed. Its two coefficients are k (gradient-ratio
# scaling, 0 linear / 1 fully limited) and l (permitted overshoot relative to the face
# difference), so 0.2 0 is light detection with no overshoot allowed.
#
# THE RISK, STATED BEFORE THE RUNS. alletto_elim dies at ARGUS Re. If the mechanism is ANY
# energy limiting rather than the AMOUNT, k=0.3 dies too and only filteredLinear2 is
# interesting, because it is the only one of the three that adds dissipation selectively.
for _nm, _sch in (("alletto_k03", "Gauss limitedLinear 0.3"),
                  ("alletto_flin2", "Gauss filteredLinear2 0.2 0"),
                  ("alletto_vanleer", "Gauss vanLeer")):
    RECIPES[_nm] = dict(
        RECIPES["alletto_lim"],
        div=dict(RECIPES["alletto_lim"]["div"], **{"div(phi,h)": _sch,
                                                   "div(phi,e)": _sch}),
    )
    RECIPE_AXES[_nm] = dict(RECIPE_AXES["alletto_lim"], energy=_sch)
    _sa = dict(RECIPES[_nm])
    _dv = {k: v for k, v in _sa["div"].items() if "phi,k" not in k and "phi,omega" not in k}
    _dv["div(phi,nuTilda)"] = RECIPES["alletto_lim"]["div"].get("div(phi,k)",
                                                                "Gauss linearUpwind limited")
    _sa["div"] = _dv
    _er = {k: v for k, v in _sa["equations_relax"].items() if "k|omega" not in k}
    _er["nuTilda"] = "0.9"
    _sa["equations_relax"] = _er
    _sa["residual_control"] = {"p": "1e-5", "U": "1e-6", '"(nuTilda|e|h)"': "1e-5"}
    _sa["solvers"] = _sa["solvers"].replace("k|omega|epsilon|nuTilda", "nuTilda")
    RECIPES["sa_" + _nm.replace("alletto_", "")] = _sa
    RECIPE_AXES["sa_" + _nm.replace("alletto_", "")] = dict(RECIPE_AXES[_nm],
                                                            model="SpalartAllmaras")


# ---- SA VARIANTS THAT FIX THE EQUATION I LEFT OUT (project decision, 2026-08-12: improve
# the flin2 SA recipe rather than change turbulence model).
#
# sa_flin2 applies filteredLinear2 to the ENERGY equation and leaves div(phi,nuTilda) on
# unlimited linearUpwind. In Spalart-Allmaras the transported variable IS the eddy viscosity,
# so an overshoot in nuTilda across the shock puts the wrong viscosity exactly where the shock
# sits -- and shock position is what the lift error is made of. The same defect I fixed in the
# energy equation is still live in the turbulence equation.
#
# THE ERROR PATTERN SAYS THIS IS THE RIGHT PLACE TO LOOK. sa_flin2 is +2.2% on case 07,
# -1.6% on 06 and -2.8% on 09: it CHANGES SIGN with alpha, which is a lift-SLOPE error, not an
# offset. A shock sitting progressively too far forward as the shock strengthens produces
# exactly that, and nothing in the energy equation would.
#
# sa_flin2_t   filteredLinear2 on nuTilda as well as energy -- the single missing entry.
# sa_flin2_k05 k = 0.5 instead of 0.2 on both, i.e. earlier detection of staggering. Tests
#              whether the amount of selective dissipation matters once it is applied
#              everywhere it should be.
for _nm, _e, _tu in (("sa_flin2_t", "Gauss filteredLinear2 0.2 0", "Gauss filteredLinear2 0.2 0"),
                     ("sa_flin2_k05", "Gauss filteredLinear2 0.5 0", "Gauss filteredLinear2 0.5 0")):
    _d = dict(RECIPES["sa_flin2"])
    _d["div"] = dict(_d["div"], **{"div(phi,h)": _e, "div(phi,e)": _e,
                                   "div(phi,nuTilda)": _tu})
    RECIPES[_nm] = _d
    RECIPE_AXES[_nm] = dict(RECIPE_AXES["sa_flin2"], energy=_e, mom=RECIPE_AXES["sa_flin2"]["mom"])


# ---------------------------------------------------------------------------
# sa_m6  THE ONLY COMPRESSIBLE SETUP IN THIS PROJECT WITH A VALIDATED ANSWER
#        BEHIND IT. Transcribed from cases/validation/onera_m6_B_flin2, which
#        reproduces AGARD AR-138 run 308 and, re-run cold on 2026-08-24, returns
#        its own committed forces to better than 0.02 counts on CL and 0.001 on CD.
#        Its run card names it for ARG-170, the 3D compressible decision.
#
# WHY THIS ENTRY HAD TO EXIST. Until now the M6 setup lived ONLY as a case, so no
# script could reach it, while the recipe that scripts COULD reach (sa_flin2) is
# the one with nothing validated behind it. Any 3D transonic ARGUS case built
# through this machinery would have inherited the unvalidated numerics by default.
#
# IT DIFFERS FROM sa_flin2 IN EXACTLY THREE PLACES, all measured from the case:
#   1. laplacian/snGrad `corrected` rather than `limited 0.5`. The RAE2822
#      validation cases also use `corrected`, so the two families with experiment
#      behind them agree and sa_flin2 is the outlier.
#   2. nNonOrthogonalCorrectors 0 rather than 1. sa_flin2's correctors are
#      justified by the 2D mesh running 38-45 degrees non-orthogonality; THE M6
#      RUNS ZERO CORRECTORS AT 64.9 DEGREES and converges in 500 iterations, which
#      undercuts that justification rather than supporting it.
#   3. nuTilda relaxation 0.8 rather than 0.9.
# The solver block also widens to "(U|e|k|epsilon|nuTilda)" so the entry works
# under kOmegaSST as well as Spalart-Allmaras.
#
# NOT A CLAIM THAT sa_m6 IS BETTER ON THE 2D SECTION. It is a claim that this is
# the setup that has been checked against measured data, and that anything
# transonic should start here and depart from it deliberately.
_m6 = dict(RECIPES["sa_flin2"])
_m6["algo"] = dict(_m6["algo"], nNonOrthogonalCorrectors="0")
_m6["equations_relax"] = dict(_m6["equations_relax"], nuTilda="0.8")
_m6["laplacian"] = "Gauss linear corrected"
_m6["sngrad"] = "corrected"
_m6["solvers"] = _m6["solvers"].replace('"(U|e|h|nuTilda)"', '"(U|e|h|k|epsilon|nuTilda)"')
# THE M6 CASE IS MODEL-AGNOSTIC AND THE RECIPE MUST BE TOO. It declares the
# `turbulence` and `energy` div ALIASES and carries k, epsilon and omega entries so
# the same dictionaries run under kOmegaSST as well as Spalart-Allmaras. Dropping
# them would make sa_m6 a Spalart-only recipe that no longer reproduces the case it
# was transcribed from, which is the one property this entry exists to have.
_m6["div"] = dict(_m6["div"], **{
    "turbulence": "Gauss linearUpwind limited",
    "energy": "Gauss linearUpwind limited",
    "div(phi,k)": "Gauss linearUpwind limited",
    "div(phi,omega)": "Gauss linearUpwind limited",
})
# BOTH ENERGY FIELDS ARE NAMED, and this is a portability fix, not a departure from the
# M6. The M6 case uses `sensibleInternalEnergy`, so it solves `e` and never sees `h`.
# The ARGUS 2D template uses `sensibleEnthalpy`, so it solves `h`. A recipe that names
# only `e` leaves `h` with NO relaxation entry, OpenFOAM defaults it to 1.0, and the run
# FPEs inside the first GaussSeidel smooth. MEASURED 2026-08-25: sa_m6 applied to the
# ARGUS 2D baseline died on iteration 1 for exactly this reason. sa_flin2 was portable
# because it always carried both. Naming both changes NOTHING about what the M6 case
# runs, because `h` is not constructed there; it only stops the recipe silently
# depending on which energy form the target case happens to use.
_m6["equations_relax"] = {"p": "1", "U": "0.9", "e": "0.8", "h": "0.8",
                          "k": "0.9", "epsilon": "0.9", "nuTilda": "0.8"}
_m6["residual_control"] = {"p": "1e-05", "U": "1e-06",
                           '"(k|epsilon|nuTilda|e|h)"': "1e-05"}
RECIPES["sa_m6"] = _m6
RECIPE_AXES["sa_m6"] = dict(RECIPE_AXES["sa_flin2"])


# ---- sa_m6_nolim: sa_m6 with the CELL-LIMITED GRADIENT REMOVED ---------------------------
#
# WHY (measured 2026-08-24, 67 controlled runs on the ARGUS 2D transonic set). The entry
# `limited  cellLimited Gauss linear 1` is a DISCRETE SWITCH: it clips or does not clip per
# cell, so an arbitrarily small field change flips it somewhere and that is a finite change
# in the discretisation. Fed back through SIMPLE it produces an exact period-2 limit cycle.
# It surfaced in the VISCOUS force, which alternated 10.656 / 9.258 N on consecutive
# iterations while the pressure force stayed converged to five digits, so Cl and Cm looked
# clean and Cd sat in a 7.55-count bracket whose value depended on the PARITY of the stopping
# iteration.
#
# TEN VARIANTS, ONE CHANGE EACH, WITH A CONTROL THAT REPRODUCED ITS PARENT TO SIX DIGITS.
# Upwinding div(phi,nuTilda), div(phi,h) or div(phid,p) individually each left the bracket at
# 7.5 counts, so no div scheme is responsible. Removing the gradient limiter collapsed it.
# Across all 19 cases: bracket 5.37 -> 0.04 counts mean, ZERO cases above 0.1. At M 0.695,
# 15 of 15 converged to 1.4e-6 and a cold start agreed with a restart to 0.013 counts, i.e.
# the fixed recipe has a genuine unique fixed point. At M 0.78, 0 of 4 converged, from a
# SECOND cause that is not this one and is not yet identified.
#
# THIS IS NOT A CLAIM THAT THE LIMITER IS WRONG. The validated M6 case KEEPS it and
# reproduces AGARD AR-138 run 308 exactly, in 3D at M 0.8395. The chatter is triggered by the
# ARGUS 2D mesh (max aspect ratio 2893, max non-orthogonality 45.4 deg), not by the recipe.
# So sa_m6 remains the default and this is the documented departure from it.
#
# HOW TO TELL WHICH YOU NEED, rather than guessing: split Cd and the residuals by iteration
# parity. Two separated levels means the limiter is biting. One line means leave it alone.
_m6n = dict(RECIPES["sa_m6"])
_m6n["grad_schemes"] = _m6n["grad_schemes"].replace(
    "cellLimited Gauss linear 1", "Gauss linear")
assert _m6n["grad_schemes"] != RECIPES["sa_m6"]["grad_schemes"], \
    "sa_m6_nolim: the gradient-limiter anchor did not match; the edit did not land"
RECIPES["sa_m6_nolim"] = _m6n
RECIPE_AXES["sa_m6_nolim"] = dict(RECIPE_AXES["sa_m6"])


def _header(obj):
    return ("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
            "    class       dictionary;\n    location    \"system\";\n"
            "    object      %s;\n}\n\n" % obj)


def _fvschemes(spec):
    div = "\n".join("    %-38s %s;" % (k, v) for k, v in spec["div"].items())
    return (_header("fvSchemes")
            + "ddtSchemes\n{\n    default         %s;\n}\n\n" % spec["ddt"]
            + spec.get("grad_schemes", GRAD_SCHEMES) + "\n\n"
            + "divSchemes\n{\n%s\n}\n\n" % div
            + "laplacianSchemes\n{\n    default         %s;\n}\n\n"
              % spec.get("laplacian", "Gauss linear limited 0.5")
            + "interpolationSchemes\n{\n    default         linear;\n}\n\n"
            + "snGradSchemes\n{\n    default         %s;\n}\n\n"
              % spec.get("sngrad", "limited 0.5")
            + "wallDist\n{\n    method          meshWave;\n}\n")


def _fvsolution(spec):
    algo = "\n".join("    %-28s %s;" % (k, v) for k, v in spec["algo"].items())
    rc = spec.get("residual_control")
    if rc:
        algo += ("\n    residualControl\n    {\n"
                 + "".join("        %-24s %s;\n" % (k, v) for k, v in rc.items())
                 + "    }")
    t = (_header("fvSolution") + spec["solvers"] + "\n\n"
         + "%s\n{\n%s\n}\n" % (spec["algo_block"], algo))
    if spec["fields_relax"] is not None:
        f = "".join("        %-24s %s;\n" % (k, v)
                    for k, v in spec["fields_relax"].items())
        e = "".join("        %-24s %s;\n" % (k, v)
                    for k, v in spec["equations_relax"].items())
        t += ("\nrelaxationFactors\n{\n    fields\n    {\n%s    }\n"
              "    equations\n    {\n%s    }\n}\n" % (f, e))
    return t


def _expand(case, rel):
    """Read a dictionary back through OPENFOAM'S OWN PARSER, not through my regexes.

    -expand strips every comment and resolves every $variable, so what comes out is what the
    solver will act on. A check that re-reads the writer's own text with the writer's own
    pattern proves only that the writer is self-consistent.
    """
    r = subprocess.run(["bash", "-lc", FOAM + "foamDictionary -expand system/%s" % rel],
                       cwd=str(case), capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        raise SystemExit("foamDictionary could not read system/%s in %s:\n%s"
                         % (rel, case, (r.stderr or r.stdout)[-800:]))
    return r.stdout


def _block(text, name):
    """Body of a named brace block in already-expanded (comment-free) dictionary text."""
    m = re.search(r"(?m)^\s*%s\s*\n\s*\{" % re.escape(name), text)
    if not m:
        return None
    i = text.index("{", m.start())
    depth, j = 0, i
    while j < len(text):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[i + 1:j]
        j += 1
    return None


def _entries(body):
    """Flat `key value;` pairs at the TOP level of a block body, nested blocks skipped."""
    out, depth, buf = {}, 0, ""
    for ch in body:
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            buf = ""
            continue
        if depth == 0:
            if ch == ";":
                parts = buf.strip().split(None, 1)
                if len(parts) == 2:
                    out[parts[0]] = " ".join(parts[1].split())
                buf = ""
            else:
                buf += ch
    return out


def _same(a, b):
    """Compare two dictionary values, NUMERICALLY when both are numbers.

    foamDictionary re-emits what it parsed, so `1e-6` comes back as `1e-06`. A pure string
    comparison then reports a mismatch on a file that is exactly right, which is a gate that
    ALWAYS FIRES -- as useless as one that cannot, and quicker to get switched off. Numbers
    compare as numbers; everything else, including scheme names, stays exact.
    """
    if a is None or b is None:
        return a == b
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return " ".join(str(a).split()) == " ".join(str(b).split())


def _same_map(got, want):
    return set(got) == set(want) and all(_same(got[k], want[k]) for k in want)


def verify(case, name, verbose=True):
    """Assert the case IS in recipe `name`, reading every claim out of the case (GEO-080).

    Returns a list of problems. Empty list means verified. Nothing is asserted from the fact
    that apply() ran; the same function is used to audit a case nobody here built.
    """
    case, spec, bad = pathlib.Path(case), RECIPES[name], []
    sch, sol = _expand(case, "fvSchemes"), _expand(case, "fvSolution")
    cd = _expand(case, "controlDict")

    app = _entries(cd).get("application")
    if app != spec["application"]:
        bad.append("controlDict application is %s, recipe %s wants %s"
                   % (app, name, spec["application"]))

    ddt = _entries(_block(sch, "ddtSchemes") or "").get("default")
    if ddt != spec["ddt"]:
        bad.append("ddtSchemes default is %s, want %s" % (ddt, spec["ddt"]))

    got_div = _entries(_block(sch, "divSchemes") or "")
    for k, v in spec["div"].items():
        if not _same(got_div.get(k), v):
            bad.append("divSchemes %s is %r, want %r" % (k, got_div.get(k), v))
    # EXTRA ENTRIES ARE A FAILURE TOO, not a skip: a stale div(phid,p) left behind next to a
    # new one is exactly the half-applied state this module exists to make impossible.
    for k in set(got_div) - set(spec["div"]):
        bad.append("divSchemes has UNEXPECTED entry %s %s" % (k, got_div[k]))

    # THE LAPLACIAN AND snGrad ARE PART OF THE RECIPE, NOT SCENERY. Alletto runs `corrected`
    # where the wolf set runs `limited 0.5`; on a noMax 45 mesh with nNonOrthogonalCorrectors 0
    # that is a live stability difference, so it is verified rather than left to the writer.
    for blk_name, key, dflt in (("laplacianSchemes", "laplacian", "Gauss linear limited 0.5"),
                                ("snGradSchemes", "sngrad", "limited 0.5")):
        want = spec.get(key, dflt)
        got = _entries(_block(sch, blk_name) or "").get("default")
        if not _same(got, want):
            bad.append("%s default is %r, want %r" % (blk_name, got, want))

    got_algo = _entries(_block(sol, spec["algo_block"]) or "")
    if _block(sol, spec["algo_block"]) is None:
        bad.append("fvSolution has no %s block" % spec["algo_block"])
    for k, v in spec["algo"].items():
        if not _same(got_algo.get(k), v):
            bad.append("%s/%s is %r, want %r" % (spec["algo_block"], k, got_algo.get(k), v))
    rc_want = spec.get("residual_control")
    rc_body = _block(_block(sol, spec["algo_block"]) or "", "residualControl")
    if rc_want:
        got_rc = _entries(rc_body or "")
        if not _same_map(got_rc, rc_want):
            bad.append("%s/residualControl is %s, want %s"
                       % (spec["algo_block"], got_rc, rc_want))
    elif rc_body is not None:
        bad.append("recipe %s takes no residualControl but the block is present" % name)
    other = "PIMPLE" if spec["algo_block"] == "SIMPLE" else "SIMPLE"
    if _block(sol, other) is not None:
        bad.append("fvSolution still carries a %s block beside the %s one"
                   % (other, spec["algo_block"]))

    rf = _block(sol, "relaxationFactors")
    if spec["fields_relax"] is None:
        if rf is not None:
            bad.append("recipe %s takes NO relaxation factors but the block is present" % name)
    elif rf is None:
        bad.append("fvSolution has no relaxationFactors block")
    else:
        for sub, want in (("fields", spec["fields_relax"]),
                          ("equations", spec["equations_relax"])):
            got = _entries(_block(rf, sub) or "")
            if not _same_map(got, want):
                bad.append("relaxationFactors/%s is %s, want %s" % (sub, got, want))

    if verbose:
        tag = "VERIFIED" if not bad else "MISMATCH"
        print("  recipe %-8s %-8s %s" % (name, tag, case))
        for b in bad:
            print("      %s" % b)
    return bad


def pressure_clamp_admits(case, cp_floor=CP_FLOOR, verbose=True):
    """The pressure clamp must admit the pressure the case will physically reach.

    A ONE-DIRECTIONAL GATE (D060). pressureControl can only ever pull the field TOWARD the
    clamp, so a clamp inside the physical range does not add error in an unknown direction: it
    removes suction peak and it removes stagnation pressure, biasing Cp, cl and the shock
    position all one way. There is no tolerance to argue about.

    Both bounds are DERIVED, not chosen:
        pMin must sit below  p_inf + cp_floor * q_inf     (the suction side)
        pMax must sit above  p_inf * (1 + 0.2 M^2)^3.5    (isentropic stagnation; a shock can
                                                           only LOWER downstream p0, never
                                                           raise it, so this is a hard ceiling)
    """
    sys.path.insert(0, str(REPO / "scripts"))
    import case_normalisation
    n = case_normalisation.derive(pathlib.Path(case))
    p_inf, q_inf, M = n["flow"]["p_inf"], n["normalisation"]["q_inf"], n["flow"]["Mach"]
    sol = _expand(pathlib.Path(case), "fvSolution")
    blk = _block(sol, "SIMPLE") or _block(sol, "PIMPLE") or ""
    e = _entries(blk)
    if "pMinFactor" not in e and "pMaxFactor" not in e:
        if verbose:
            print("  clamp     NONE SET  (pressureControl inactive)")
        return []
    # pressureControl takes its base from the VALUE-FIXING patches. On these cases that is the
    # farfield/outlet at p_inf, which is what n["flow"]["p_inf"] reads out of 0/p.
    p_min = float(e.get("pMinFactor", 0)) * p_inf
    p_max = float(e.get("pMaxFactor", 1e30)) * p_inf
    need_min = p_inf + cp_floor * q_inf
    need_max = p_inf * (1.0 + 0.2 * M * M) ** 3.5
    bad = []
    if p_min > need_min:
        bad.append("pMinFactor %s clamps p at %.0f Pa; Cp %.1f needs %.0f Pa (Cp at the clamp "
                   "is only %.3f)" % (e.get("pMinFactor"), p_min, cp_floor, need_min,
                                      (p_min - p_inf) / q_inf))
    if p_max < need_max:
        bad.append("pMaxFactor %s clamps p at %.0f Pa; isentropic stagnation at M %.3f is "
                   "%.0f Pa" % (e.get("pMaxFactor"), p_max, M, need_max))
    if verbose:
        print("  clamp     %-8s p in [%.0f, %.0f] Pa; physical range needs [%.0f, %.0f] "
              "(Cp_floor %.1f, p0 at M %.3f)"
              % ("OK" if not bad else "TOO TIGHT", p_min, p_max, need_min, need_max,
                 cp_floor, M))
        for b in bad:
            print("      %s" % b)
    return bad


def apply(case, name, iters=None, verbose=True):
    """Write the recipe, then VERIFY it from the case. Never returns having only written."""
    case = pathlib.Path(case)
    if name not in RECIPES:
        raise SystemExit("unknown recipe %r; have %s" % (name, ", ".join(RECIPES)))
    spec = RECIPES[name]
    (case / "system/fvSchemes").write_text(_fvschemes(spec))
    (case / "system/fvSolution").write_text(_fvsolution(spec))
    cd = case / "system/controlDict"
    t = cd.read_text()
    t = re.sub(r"^application\s+\S+;", "application     %s;" % spec["application"],
               t, flags=re.M)
    if iters:
        t = re.sub(r"^endTime\s+\S+;", "endTime         %d;" % iters, t, flags=re.M)
    # LTS advances one PSEUDO-step per iteration and sets its own local dt from maxCo, so
    # adjustTimeStep must be off and deltaT is the iteration counter.
    t = re.sub(r"^adjustTimeStep\s+\S+;", "adjustTimeStep  no;", t, flags=re.M)
    t = re.sub(r"^deltaT\s+\S+;", "deltaT          1;", t, flags=re.M)
    t = re.sub(r"^writeControl\s+\S+;", "writeControl    timeStep;", t, flags=re.M)
    cd.write_text(t)
    bad = verify(case, name, verbose=verbose)
    bad += pressure_clamp_admits(case, verbose=verbose)
    if bad:
        raise SystemExit("recipe %s did not land in %s:\n  %s"
                         % (name, case, "\n  ".join(bad)))
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("case", nargs="+")
    ap.add_argument("--recipe", choices=sorted(RECIPES), required=True)
    ap.add_argument("--iters", type=int, default=None)
    ap.add_argument("--verify-only", action="store_true",
                    help="audit an existing case against the recipe; change nothing")
    ap.add_argument("--json", default=None, help="write the verdict per case")
    a = ap.parse_args()
    out, rc = {}, 0
    for c in a.case:
        print("%s" % c)
        if a.verify_only:
            bad = verify(c, a.recipe) + pressure_clamp_admits(c)
        else:
            try:
                apply(c, a.recipe, a.iters)
                bad = []
            except SystemExit as e:
                bad = [str(e)]
        out[c] = {"recipe": a.recipe, "problems": bad}
        rc |= 1 if bad else 0
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=2) + "\n")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
