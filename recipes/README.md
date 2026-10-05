# recipes: the setup files, one folder per (recipe, condition)

**This is the only place a solver setup is defined.** A recipe folder holds `0.orig/`,
`constant/`, `system/` and nothing else: no geometry, no mesh, no results. Build a case by
copying a folder and adding the geometry.

Every folder carries `MANIFEST.sha256` (integrity) and `RECIPE_FILES.json` (provenance).
Verify one with `cd <folder> && sha256sum -c MANIFEST.sha256`.

## The four recipes

| folder | regime | model | wall treatment | condition |
|---|---|---|---|---|
| `L11_wallResolved_WT` | incompressible | kOmegaSST | wall-RESOLVED, `nutLowReWallFunction` | M 0.12, U 40.8, Re_cref 0.95e6 |
| `L11_wallResolved_CR` | incompressible | kOmegaSST | wall-RESOLVED, `nutLowReWallFunction` | M 0.10, U 34.0, Re_cref 0.79e6 |
| `M6_wallModelled_CRUISE_EARLY` | compressible | Spalart-Allmaras | wall-MODELLED, `nutUSpaldingWallFunction` | M 0.78, CL 0.52929708745, 74,500 kg at 10 km |
| `M6_wallModelled_CRUISE_LATE` | compressible | Spalart-Allmaras | wall-MODELLED, `nutUSpaldingWallFunction` | M 0.78, CL 0.54411480021, 56,000 kg at 12 km |

That is the whole picture. There is one incompressible recipe at two conditions, and one
compressible recipe at two conditions.

## WAKE REFINEMENT, added 2026-09-23. Read this before running anything.

`L11_wallResolved_CR`, `M6_wallModelled_CRUISE_EARLY` and `M6_wallModelled_CRUISE_LATE` had
their `system/snappyHexMeshDict` changed in three places, and in three places only:

```
maxGlobalCells      40000000;   ->   maxGlobalCells      130000000;

geometry{}          + box_tref_a  min ( 1.80 0.00 -0.60) max ( 3.40 2.20  0.60)
                    + box_tref_b  min ( 1.80 0.00 -2.00) max ( 3.40 2.60  2.00)

refinementRegions{} + box_tref_a  { mode inside; levels ((1e15 7)); }
                    + box_tref_b  { mode inside; levels ((1e15 6)); }
```

Those are the only wake-refinement changes. The same three folders were later aligned with the
r2 cases in two more files: `system/decomposeParDict` (all three) and `system/controlDict`
(`L11_wallResolved_CR` only). See "Decomposition and run length" below.

**No turbulence model, wall treatment or layer setting was touched.** The compressible
recipes remain wall-modelled at y+ ~80, which is the production intent and unchanged.

**Why the cap was the defect.** 40M was binding on **all 18 meshes**. W, H and F exhausted it,
so *shell* refinement received nothing at all: final level-7 count 94,543 against B/C/M's
5.87M. **Corrected 2026-09-24:** an earlier version of this paragraph said W, H and F used up
the cap in level-11 *feature* refinement, because W's `wing.eMesh` carries 2067 points against
B's 1273. The meshing logs show otherwise. The F, H and W mesh cases held surfaces with an
unwelded root seam (796 open edges, which is also where the extra `wing.eMesh` points come
from). At surface level 11 that seam let snappyHexMesh's first cell-removal step keep the wing
interior, and refining the interior used up the budget: 43.5M to 44.3M cells after the first
shell iteration against 29.6M to 33.8M for B/C/M, while feature refinement added only 1k to 25k
cells. A later removal step deleted the interior, so the final meshes are not flooded. The
conclusions below (raise the cap; 130M; the wake boxes) are unchanged. The starved published
W, H and F meshes' Trefftz induced drag self-certified as corrupt, returning span efficiency
`e > 1`, which is physically impossible for a planar wake. Build F, H and W from the
registered, welded surfaces. A first raise to 70M was **also**
bound, at shell iteration 9 on SWW_r2 (69,599,457 of 70,000,000), and a truncation there would
have landed on whichever region refines last, which no output gate can see. 130M is the value
all six wake-refined meshes carry. Meshing memory is not the constraint: 68.6M cells measured
at 115 GB of the 1004 GB allocated.

**Level 7 is 7.8125 mm**, the wake-core cell size the healthy meshes already had. The starved
meshes carried 62.5 mm there, three levels coarser. That is the whole story of the `e > 1`
artefact.

**`L11_wallResolved_WT` was deliberately NOT changed** and still carries 40M. The wake boxes
exist to resolve the Trefftz integration planes for induced drag, which is not what the
wind-tunnel condition is for. If a WT run ever needs a wake integral, apply the same three
changes and say so in a decision entry. It did not receive the decomposition or run-length
changes either.

### THE FRAME WARNING, and it is the important part

**Every delivered result is on this recipe's mesh settings at its trim.** Condition CR was
re-run on it first (the wake-refined generation). The early- and late-cruise trims were re-run on
it from the registered, welded surfaces between 25 September and 2 October 2026 (the welded
generation): in the B, F and W mesh cases at both cruise points, 13 of the 14 files under
`system/` and `constant/` of `M6_wallModelled_CRUISE_EARLY` and `_LATE` are byte-identical to this
folder, and `controlDict` differs only by the staging steps (`purgeWrite` 0, the post-processing
include, the case's own force axes). So:

1. A new early- or late-cruise run built from this folder reproduces the delivered cruise set-up.
   No case has been rebuilt from a fresh clone and compared number for number, so do not treat
   the delivered numbers as a regression test until that has been done.
2. The earlier cruise results, on the old recipe, are superseded; `data/rans_forces.json` keeps
   them under `superseded_entry` only.
3. Any table placing a Condition CR number beside a cruise number is crossing an operating point
   (and a surface generation for F, H and W, see the README, Known gaps, item 8). State both.

`MANIFEST.sha256` was regenerated for the three changed folders and verifies clean.

**The two wall treatments are both deliberate and both sanctioned.** Wall-resolved is the
incompressible programme's basis and wall functions there are off-limits without a decision
entry. Wall-modelled at y+ ~80 is the compressible production intent (ARG-152, and ARG.md:3596:
*"validating wall-resolved and then running wall-modelled would validate a recipe we do not
use"*). The sanction is recipe-scoped in code, so applying one regime's rule to the other is
an error. It has been made.

## Decomposition and run length, aligned with the r2 cases

The r2 cases are the six wake-refined Condition CR cases (B, C, F, H, M, W), meshed and solved
on HPC12 with OpenFOAM-12. `L11_wallResolved_CR` and the two cruise folders now carry the
decomposition those cases ran, and `L11_wallResolved_CR` carries their cold-leg run length.

**Two decompositions, one per stage.**

| stage | ranks | method | set by |
|---|---|---|---|
| mesh (parallel snappyHexMesh) | 192 | `scotch` | `system/decomposeParDict` in this folder |
| solve (foamRun) | 128 by default | `hierarchical`, 8 by 2 by 8 | `hpc/stage_solve_case.sh`, which re-decomposes the solve case |

1. **The mesh decomposition is part of the mesh recipe.** Parallel snappyHexMesh is
   decomposition-dependent: a different rank count or method builds a different mesh. All six
   r2 mesh cases carry 192 subdomains with `method scotch`, and their snappyHexMesh logs report
   `nProcs : 192`. `hpc/mesh_hpc12.pbs` reads the rank count from this file and refuses an
   allocation with fewer slots.
2. **The solve case is re-decomposed.** All six r2 solves ran on 128 ranks, hierarchical; the
   8 by 2 by 8 split matches the 128-rank case definitions. The 128-rank legs ran on four 32-core
   typej nodes (fpt-medium) or on one typen node at 128 cores (fpt-large). Re-decomposing a
   finished mesh does not change it. `hpc/solve_hpc12.pbs` reads the rank count from the case,
   so a solve case that kept this folder's `decomposeParDict` would run on 192 ranks.
3. **What the method costs.** Measured on the baseline at 4 levels below production, scotch
   against hierarchical: wing wetted area 1.28019 against 1.28026 m2 (0.005 percent), total
   cells 243,379 against 241,911 (0.60 percent). Small, but it is a mesh change, which is why
   the method is not a scheduling setting. DelftBlue's openfoam-org/12 module has no scotch,
   so meshing there means editing this file, and therefore changing the mesh.
4. `L11_wallResolved_WT` is unchanged and keeps its older setting: 192 ranks, hierarchical,
   8 by 3 by 8.

Whether the meshes behind the shipped EC and LC numbers were decomposed this way is not recorded
in this folder. The frame warning above already says those numbers are not a regression test
for a run built from here.

**`L11_wallResolved_CR/system/controlDict` carries the r2 cold leg.**

| entry | was | now |
|---|---|---|
| `endTime` | 4500 | 4000 |
| `writeInterval` | 150 | 500 |

Every other byte is unchanged, including `liftDir` and `dragDir` at full precision.
`hpc/solve_hpc12.pbs` refuses an `endTime` that is not a multiple of `writeInterval`, and when a
leg reaches `endTime` unconverged it asks for the next leg with `endTime` raised by 1500
(passed as `ARGUS_ENDTIME`): 4000, then 5500, then 7000, all multiples of 500. All six r2 cases
finished at `endTime` 6000. Their controlDicts also carry `#include "argusPostPro"`; the recipe
does not, because `scripts/install_postpro.py` adds it to the case.

**The compressible controlDicts were not changed**, because their run lengths were never part of
r2. Their `endTime` 10000 with `writeInterval` 2000 does **not** survive the +1500 raise:
11500 is not a multiple of 2000, so that next leg would be refused. Submit it with the next
multiple of 2000 instead, and settle this before running a compressible case that may need a
second leg.

## Where the compressible numerics come from

`system/fvSchemes` and `system/fvSolution` are **byte-copied** from
`cases/validation/onera_m6_B_flin2`, the only compressible setup in this repository with a
validated answer behind it (AGARD AR-138 run 308). Verify with:

```
sha256sum cases/validation/onera_m6_B_flin2/system/fvSchemes \
          hpc/recipes/M6_wallModelled_CRUISE_EARLY/system/fvSchemes
```

`constant/turbulenceProperties` is **derived, not byte-copied**, and the reason is worth
knowing. The M6's file carries a commented-out `RASModel kOmegaSST` directly above the live
`RASModel SpalartAllmaras`. OpenFOAM ignores it, so the M6 ran Spalart-Allmaras correctly, but
a commented alternative in a recipe file is exactly the hazard `solver_recipe.py` was written
against, and this repository has been bitten by that shape four times. The shipped file has one
`RASModel` line and no alternatives. Equivalence is verified with OpenFOAM's own reader rather
than a text match:

```
diff <(foamDictionary -expand cases/validation/onera_m6_B_flin2/constant/turbulenceProperties) \
     <(foamDictionary -expand hpc/recipes/M6_wallModelled_CRUISE_EARLY/constant/turbulenceProperties)
```

which differs only in the file path header. The same rule is applied to the recipes' own files:
no dictionary here carries a settable value inside a comment.

Boundary-condition **types** are the M6's, re-expressed on the ARGUS patch names
(`symm`→`symmetry`, `Back`→`outboard`, `topAndBottom`→`topBottom`). Reference quantities are
the DSO basis at model scale and match `L11_*` exactly, so incompressible and compressible
coefficients are directly comparable.

### Relationship to `scripts/solver_recipe.py`

They are complementary, not alternatives. `solver_recipe.py` writes the three **solver**
dictionaries (`fvSchemes`, `fvSolution`, and the algorithm block) and verifies them through
`foamDictionary` rather than its own regexes. It does not carry meshes, boundary conditions or
operating points. This folder carries those.

`RECIPES["sa_m6"]` **is** a faithful transcription of the validated M6 setup: 0 non-orthogonal
correctors, nuTilda 0.8, laplacian and snGrad `corrected`, the M6's div aliases, relaxation and
residual control, and it is deliberately model-agnostic so the same dictionaries run under
kOmegaSST or Spalart-Allmaras. It was registered in commit `8da89d5`. Using
`--recipe sa_m6` on a case built from this folder is consistent and is the intended route.

Two things not to misread. The 31 entries are not 31 production options: three are real
(`wolf`, `tutorial`, `lts`), `sa_m6` is the validated compressible one, and the rest are trial
variants from a numerics study driven by `overnight_recipe_study.py`. And ARG-187 records that
no registered recipe was the M6 setup; **that gap has since been closed** by `8da89d5`, so the
entry is history rather than current state.

**`sa_m6` had an energy-form portability bug, fixed 2026-08-25.** It named relaxation for `e`
only. The M6 case uses `sensibleInternalEnergy` and solves `e`, so it was fine there. The ARGUS
2D template uses `sensibleEnthalpy` and solves `h`, which was then left with no relaxation
entry, defaulted to 1.0, and the run died inside the first GaussSeidel smooth. Measured: applying
`sa_m6` to the ARGUS 2D baseline FPE'd on iteration 1. `sa_flin2` was portable because it always
named both. Both fields are now named, which changes nothing about what the M6 case runs because
`h` is never constructed there. **The general point: a recipe that names a field by its energy
form silently depends on the target case's thermo, and nothing in the recipe says so.**

**`sa_m6_nolim`** is registered beside it: `sa_m6` with the cell-limited gradient removed, which
is the documented departure for a mesh where the limiter chatters. See the section below.

## The gradient limiter, and the check that decides it

Both compressible folders ship `system/fvSchemes` with `cellLimited Gauss linear 1`, exactly as
the M6 validates it, and a variant `variants/fvSchemes.nolimiter` with the gradient limiter
removed.

**Default to the shipped one.** On the ARGUS 2D transonic mesh that limiter chatters: it is a
discrete switch, it flips state on alternate iterations, and it put Cd in a 7.55-count period-2
limit cycle across all 19 cases while Cl and Cm looked clean. Removing it collapsed the cycle to
0.04 counts and dropped the Ux residual from 2.6e-4 to 1.4e-6. But the M6 keeps the same limiter
in 3D at M 0.8395 and reproduces exactly, so this is a defect of that 2D mesh (max aspect ratio
2,893, max non-orthogonality 45.4 deg), not of the recipe.

So detect rather than pre-emptively change. Split the drag and the residuals by iteration parity:

```
awk '/^Time = /{t=$3} /Solving for Ux/{split($0,a,"residual = ");
     printf "%d %s %s\n", t%2, t, a[2]}' log.solver | tail -40
```

Two separated levels means the limiter is biting and you switch to the variant. One line means
leave it alone. The same split on `forceCoeffs.dat` odd versus even iterations is more direct.

## What is NOT settled

1. **Scale.** The cruise condition is defined at aircraft scale (chord 3.196 m, b 38.36 m).
   `geometry/derived/wing3d/*.stl` is model scale (half-span 1.8254 m, Cref 0.393957 m). The
   recipes resolve this by **Re-matching**: model-scale geometry with `mu` adjusted so Re on
   Cref equals the aircraft cruise Re. That is the same device the 2D transonic work already
   uses, but **it is not ratified by a decision entry**. It is the one modelling choice in these
   folders that has not been signed off. See `CONDITION.json` → `scale_resolution`, which also
   states the two alternatives. Change it there, not in `thermophysicalProperties`.
2. **Alpha.** `0.orig/U` is set from the VSPAERO/VLM trim alpha (0.8655 deg early, 0.9967 deg
   late). That is a **starting guess, not a RANS trim.** Each geometry still has to be trimmed
   to its own alpha at the target CL to |dCL| < 1e-4.
3. **M 0.78 convergence.** In 2D, the limiter fix converged 15 of 15 cases at M 0.695 but 0 of 4
   at M 0.78, from a second cause not yet identified. The M6 runs clean at M 0.8395, so this may
   be specific to the 2D setup, but it is unproven at M 0.78 in 3D.

## Correction to `docs/ARGUS_reference_data.md` section 9

Section 9 item 6 says of the transonic condition: *"Both CLs must be pinned before any such
phase is scoped, and neither is currently in this file."* **That is stale.** Both were pinned on
2026-08-09 (ARG-152) and are confirmed independently in the DSO delivery at
`dso_reference/ARGUS_CDIW_audit_delivery_2026-08-19/01_cruise_cdiw_audit/rigid_baseline/state_summary.csv`:

```
early_cruise  target_CL 0.5292970874539268   alpha_trim 0.8654766754810572 deg
late_cruise   target_CL 0.5441148002120662   alpha_trim 0.9966996379212477 deg
```

They were never back-filled into section 9, which is exactly the failure that section warns
about: an operating point that lives only in a decision log. Section 9 is owned by the geometry
workstream and is not edited here; this note records the discrepancy until it is.

## The two Mach conventions are one condition, not two

The 2D transonic work runs `normal` and `freestream` as a pair (project decision). They are the
**same** M 0.78 flight condition expressed two ways: the wing has 27 deg of quarter-chord sweep,
so a section sees M_n = 0.78 cos 27 = 0.695. The `normal` convention takes both the Mach **and**
the chord normal to the quarter-chord, `Re_n = rho (V cos L)(c cos L) / mu`, because taking one
without the other is half a transformation.

**This distinction is 2D only.** A 3D wing is at freestream M 0.78, so these recipes carry that
and no convention flag.
