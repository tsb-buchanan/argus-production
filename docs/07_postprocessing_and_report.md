# 7. Post-processing, run cards and the report

> **At a glance**
> 1. **What it does:** installs the post-processing block in a case before it runs and checks it;
>    after the run, reads the forces and the `.converged` verdict, checks what the case wrote,
>    samples y+ and spanwise cuts where the mesh is, derives and validates run cards, regenerates
>    the report figures and rebuilds the report.
> 2. **What you need:** a case built as in [03](./03_build_a_case.md) or [04](./04_send_to_hpc.md)
>    (for section 4), run as in [05](./05_submit_and_run.md) and [06](./06_trim_to_target_cl.md),
>    copied back as in [04](./04_send_to_hpc.md); OpenFOAM-org 12, Python and LaTeX as in
>    [01](./01_setup.md).
> 3. **What you have at the end:** window-mean C_L and C_D; for a case that carried the block,
>    the pressure and viscous drag, the wing-surface file and the four Trefftz planes; the
>    case-derived run-card fields; the six delivered run cards validated; 7 of the 23 report
>    figures regenerated; the 34-page report rebuilt; the current Condition CR table
>    (section 9.2). A complete run card for a new run, and the induced drag of a new case's
>    Trefftz planes, need tools that are not shipped (section 11).

This guide covers what happens after a solve has finished: reading the forces, getting the
wake and surface data, writing and checking a run card, and turning the data into the figures
and the all-geometry report (`report/argus_rans_validation_report.pdf`). Its section 4 is the one step that must
happen **before** the solve: a case without the post-processing block never writes the wake
and surface data, and they cannot be produced afterwards.

It follows [06_trim_to_target_cl.md](./06_trim_to_target_cl.md). The toolchain is
OpenFOAM-org 12 on your own machine and on HPC12. Setting it up is in
[01_setup.md](./01_setup.md); what the three operating points and the two mesh generations
mean is in the [README](../README.md).

> **READ THIS FIRST.** Most of the post-processing tools ship and run from this repository. A
> few scripts still point at directories of the source repository that are not here, or need
> inputs that are not shipped. This guide says, script by script, what runs from this
> repository, what fails and why. Where a step cannot be done from this repository, the guide
> says so instead of giving you a command that will fail later.

## Contents

1. [Status at a glance](#1-status-at-a-glance)
2. [Before you start](#2-before-you-start)
3. [What a case writes](#3-what-a-case-writes)
4. [The post-processing block and the solve-job gate](#4-the-post-processing-block-and-the-solve-job-gate)
5. [Forces and the convergence verdict](#5-forces-and-the-convergence-verdict)
6. [Wake, surface and y+ data](#6-wake-surface-and-y-data)
7. [Run cards](#7-run-cards)
8. [From data to figures and the report](#8-from-data-to-figures-and-the-report)
9. [Quoting results](#9-quoting-results)
10. [Script reference](#10-script-reference)
11. [Post-processing pieces that are not in this repository](#11-post-processing-pieces-that-are-not-in-this-repository)

## 1. Status at a glance

| Task | Works from this repository? | Section |
|---|---|---|
| Install the post-processing block and check it before a run | Yes, on your own machine (the script needs numpy) | 4.1 |
| Make the HPC12 solve job check the block before every leg | Yes, with `ARGUS_GATE` in every leg's `qsub -v` list (verified on a real HPC12 job) | 4.5 |
| Read forces and C_L, C_D history (`forceCoeffs1`) | Yes | 5 |
| Read the convergence verdict (`.converged`) | Yes, for cases run with `hpc/solve_hpc12.pbs` | 5 |
| Check what a finished case wrote (`--after`) | Yes, but three of its eight lines are wrong on every real mesh: check those with `grep` | 4.6 |
| Spanwise cuts of p and U (`planeSample`) | Only where the mesh and the solved fields are: on HPC12, in the decomposed case (parallel form, NOT VERIFIED in an HPC12 job). Not on a copy pulled as in 04 (6.1) | 6.1 |
| y+ and wall shear stress | Written during the solve by a case with the block. After the run, only where the mesh and the solved fields are, as for planeSample | 6.2 |
| Trefftz planes | Sampled during the solve by a case with the block | 3, 6.3 |
| Induced drag and span efficiency from the planes of a new case | **No**: the routine that integrates the planes reads directories of the source repository. The delivered values re-select byte for byte from a shipped file | 6.3 |
| Spanwise loading from a wing-surface file | Yes, for a case that carried the block | 6.4 |
| Derive the model and wall-treatment fields of a run card | Yes | 7.2 |
| Assemble a complete run card | **No**: the card builder is not shipped | 7.1 |
| Validate a run card | Yes, on your own machine (the HPC12 login node has no jsonschema) | 7.5 |
| Regenerate figures | 7 of the 23 report figures, in a separate build directory | 8.3 |
| Rebuild `report/argus_rans_validation_report.pdf` | Yes, in the build directory: 36 pages, the same text as the shipped PDF | 8.4 |

## 2. Before you start

1. Set `$REPO` to your clone, as in [01_setup.md](./01_setup.md), and work from it:

   ```bash
   export REPO=$HOME/argus-production   # or wherever you cloned it
   cd "$REPO"
   ```

   `$CASE` below means the case directory you are post-processing. Set it as an **absolute**
   path, because several steps run from `$REPO` and a relative path would then point at the
   wrong place. For a case brought back as in [04_send_to_hpc.md](./04_send_to_hpc.md), for
   example:

   ```bash
   export CASE=$HOME/argus_results/SWB_trim   # absolute path of the case you copied back
   ```

   Section 4 is the exception: there `$CASE` is the case you are about to run, as that
   section says.

2. Load OpenFOAM-org 12 and check that its post-processing utility is on your path:

   ```bash
   source /opt/openfoam12/etc/bashrc
   which foamPostProcess
   ```

   Expected: one path ending in `/bin/foamPostProcess` (on the packaged install,
   `/opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/foamPostProcess`). If nothing is printed,
   your shell has sourced a different OpenFOAM version. The comment at the top of
   `recipes/*/system/planeSample` uses the older name `postProcess`. In OpenFOAM-org 12 that
   name is only a wrapper script: it prints a notice that `postProcess` has been superseded and
   then runs `foamPostProcess` with the same arguments. Call `foamPostProcess` directly.

   On HPC12, load OpenFOAM-org 12 with the three lines of [01_setup.md](./01_setup.md),
   section B4, typed at the login prompt one at a time in the same shell (never through a
   pipe), then check:

   ```bash
   module load devtoolset/11
   module load mpi/openmpi-4.1.2
   source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc
   which foamPostProcess
   ```

   The `source` line prints two harmless lines, `dirname: missing operand` and
   `Try 'dirname --help' for more information.` Then `which`
   prints `$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/foamPostProcess`, with
   your home directory in place of `$HOME`.

3. Check your Python. The run-card validator needs `yaml` and a `jsonschema` new enough to
   provide `Draft202012Validator` (JSON Schema draft 2020-12); `scripts/install_postpro.py` and
   the figure scripts need `numpy`, and the figure scripts also need `matplotlib`.

   ```bash
   python3 -c 'import jsonschema, yaml; jsonschema.Draft202012Validator; print("ok")'
   python3 -c 'import numpy, matplotlib; print("ok")'
   ```

   Expected: `ok` twice. An `AttributeError: module 'jsonschema' has no attribute
   'Draft202012Validator'` means your `jsonschema` is too old (version 3.2.0 fails this way;
   4.19.2 works). A `ModuleNotFoundError` means the package is missing. Use a Python that passes
   both lines. This guide was checked with Python 3.11.8, jsonschema 4.19.2, PyYAML 6.0.1,
   numpy 1.26.4 and matplotlib 3.8.0.

   On HPC12 the login node's `python3` is 3.6.8, with numpy 1.12.1 and no jsonschema. It runs
   the checks of section 4, which need only the Python standard library, but validate run
   cards on your own machine.

4. If that Python is an Anaconda installation that comes first on your `PATH`, its `mpirun`
   (MPICH) also comes first, and parallel OpenFOAM runs then fail. This matters only for the
   parallel commands in sections 6.1 and 6.2. Use the two-terminal arrangement of
   [01_setup.md](./01_setup.md), section A5: Python scripts in a terminal where the Python
   above comes first, parallel OpenFOAM commands in a terminal where `mpirun` is the system's
   Open MPI.

## 3. What a case writes

A finished OpenFOAM case keeps its sampled data in `$CASE/postProcessing/`, one directory per
function object. Which function objects exist is set by `system/controlDict`. This is what the
shipped recipes set up on their own:

| Recipe | Function objects in `system/controlDict` | Separate sampling file |
|---|---|---|
| `L11_wallResolved_CR` (Condition CR) | `forceCoeffs1` only | `system/planeSample` |
| `L11_wallResolved_WT` (tunnel-matched, not one of the three delivered conditions) | `forceCoeffs1` only | `system/planeSample` |
| `M6_wallModelled_CRUISE_EARLY`, `M6_wallModelled_CRUISE_LATE` | `forceCoeffs1`, `yPlus`, `#includeFunc MachNo`, `#includeFunc residuals` | `system/planeSample` |

`L11_wallResolved_WT` has the same layout as the Condition CR recipe but its own free stream:
`magUInf 40.8` m/s and a `dragDir` for 2.10 deg. So its force-file header (section 5.2) differs
from the Condition CR one, and that is not a fault.

1. **`forceCoeffs1`** writes the force and moment coefficients of the `wing` patch to
   `postProcessing/forceCoeffs1/<start time>/forceCoeffs.dat`. Every analysis in this
   repository reads it, including the convergence gate. See section 5.
2. **`yPlus`** (cruise recipes only) writes `postProcessing/yPlus/<time>/yPlus.dat`. Installing
   the post-processing block removes it, because the block's `argusYPlus` replaces it.
3. **`#includeFunc residuals`** as written in the cruise recipes (with no field list) writes a
   `residuals.dat` holding only a time column. This was checked under OpenFOAM-org 12. The form
   that does write residuals names the fields, for example
   `#includeFunc residuals(U, p, k, omega)` for a kOmegaSST case, which writes
   `postProcessing/residuals(U,p,k,omega)/0/residuals.dat` with one column per field component.
   The block's `argusResiduals` names its fields.
4. **`system/planeSample`** is **not** run during the solve. It defines seven cuts at constant
   y (spanwise stations eta 0.30, 0.50, 0.60, 0.70, 0.80, 0.90 and 0.95) and is run afterwards
   (section 6.1). These are spanwise sections through the flow. They are **not** Trefftz planes.

### What the post-processing block adds

A case that carries the post-processing block (section 4) also writes the directories below,
all under `$CASE/postProcessing/`. The six delivered Condition CR trims carried it. A case built
from `recipes/` without the block writes `forceCoeffs1` (plus the cruise extras above) and
nothing else: **no Trefftz planes, no wing-surface file and no pressure and viscous split.**

| Directory in `postProcessing/` | Contents | Used for |
|---|---|---|
| `argusForces/<start>/forces.dat` | pressure and viscous parts of the force and moment on `wing`, every iteration. For Condition CR they are per unit density (the block sets `rhoInf 1`, because p is kinematic); a cruise block uses the density field | the C_D pressure and viscous split in the run card (see below) |
| `argusResiduals/<start>/residuals.dat` | initial residuals of Ux, Uy, Uz, p, k and omega (e and nuTilda instead of k and omega for cruise), every iteration. Each leg starts a new folder named after its first iteration | the final residuals of the run card |
| `argusYPlus/<start>/yPlus.dat` | y+ minimum, maximum and average on `wing`, at every saved iteration | the y+ fields of the run card |
| `argusWallShearStress/<start>/wallShearStress.dat` | wall shear stress summary on `wing` | |
| `argusWingSurface/<iteration>/wingSurface.vtk` | the wing surface with p, cp, the wall shear stress vector and y+ on every face, at iteration 0 and at every saved iteration; about 0.6 GB each (628,068,096 bytes for the delivered B at iteration 6000) | spanwise loading (section 6.4), pressure and friction sections, surface maps, hinge moments |
| `argusTrefftz/<iteration>/trefftz_x1.vtk` to `trefftz_x4.vtk` | four planes of constant x behind the wing carrying U and p (a cruise block adds T and rho); for Condition CR at x = 2.394, 2.626, 2.916 and 3.206 m | induced drag C_Di and span efficiency e (section 6.3) |

`purgeWrite` removes none of this. For the delivered B solve, `argusWingSurface` holds 13 files
(7,698 MiB) and `argusTrefftz` 329 MiB. Allow for that on scratch and when you copy results
back (section 5.1).

The pressure and viscous split is not in `forceCoeffs1`: its `Cl(f)` and `Cl(r)` are front and
rear, not pressure and viscous (section 5.2). From `argusForces`, the pressure part of C_D is the
pressure force vector dotted with `dragDir` (section 5.2), divided by 0.5 |U|² A_ref; the viscous
part is the same with the viscous force. For Condition CR, with forces per unit density, that
is 0.5 × 34.0² m²/s² × 0.620462 m² (half model, DSO basis). This was checked on the baseline's
first leg: at iteration 4000 the two parts added up to `forceCoeffs1`'s `Cd` at the same
iteration to eight figures.

Cutting the wake planes after the run would need the whole volume field on disk, and the volume
field is the first thing lost when cluster scratch is cleaned, so the planes are sampled during
the solve.

## 4. The post-processing block and the solve-job gate

The post-processing block is one file, `system/argusPostPro`, plus one line in
`system/controlDict` that includes it. It makes the solver write what this guide needs:

1. the pressure and viscous parts of the forces;
2. the wing-surface fields p, cp, y+ and the wall-shear-stress vector;
3. four Trefftz planes behind the wing;
4. the residuals.

`scripts/install_postpro.py` installs it. `scripts/assert_case_postpro.py` checks a case against
eight requirements: with `--before` it reads `system/controlDict` (and what it includes) and
says what the case **will** write; with `--after` it reads `postProcessing/` and says what the
case **did** write. The `--before` check is an input check: it tells you before a multi-day
solve that the data you want will not be written.

Install the block once, **on your own machine**, in the case you build from the recipe,
**before** the case goes to HPC12. The mesh job ignores the block, the solve case inherits it,
and the solve job checks it before every leg. In this section `$CASE` is the case you are about
to run: on your own machine in 4.1, and the solve case on HPC12 in 4.3 to 4.5, where `$REPO` is
the copy of this repository you put on the cluster ([04](./04_send_to_hpc.md), Step 4).

### 4.1 Install the block (your machine, before meshing)

Run this in your Python terminal ([01](./01_setup.md), section A5), because the script needs
numpy. Do it after `wing.stl` is in `constant/triSurface` ([04](./04_send_to_hpc.md), Step 5c;
[03](./03_build_a_case.md), Step 5 in the tutorial) and before you send the case
([04](./04_send_to_hpc.md), Step 6). You can set the angle of
attack ([06](./06_trim_to_target_cl.md)) before or after this step: the block depends on the
flow speed, not on the angle.

```bash
python3 $REPO/scripts/install_postpro.py "$CASE"
python3 $REPO/scripts/assert_case_postpro.py --before "$CASE"
```

**What you should see** (baseline B at Condition CR; the first line is your case folder's name):

```
  B_CR
     solver        foamRun -> rho mode 'rhoInf', rhoInf 1, pInf 0
     wall patch    wing           (from existing forces object)
     CofR          1.34 0 0       (from existing forces object)
     turbulence    kOmegaSST      -> residual fields: U p k omega
     |U| for cp    34.0000        (from 0.orig/U (via $Uinf))
     TE x 2.3363, chord 1.1597 (from ASCII STL wing.stl)
     Trefftz x     2.3942 2.6262 2.9161 3.2061
     controlDict now includes argusPostPro
```

Then the check prints a line `=== <your case path> (before) ===`, eight `ok` lines, and
ends with the line below.

```
  ---- 8 of 8 requirements satisfied ----
```

Before the install, the same check on a case built from `recipes/L11_wallResolved_CR` prints
this, with two explanation lines (`consumer:` and `why:`) under each `MISSING`, and exits with
status 1. The cruise recipes give the same result.

```
  ok      compressible_cp_frame
  MISSING forces_split
  MISSING wing_surface_p
  MISSING wing_surface_tau_VECTOR
  MISSING wing_surface_yplus
  MISSING surface_format_carries_topology
  MISSING trefftz_planes
  MISSING residuals
  ---- 1 of 8 requirements satisfied ----
```

How to read the numbers:

1. `TE x` and `Trefftz x` are in metres, in the frame of `wing.stl`, with x pointing downstream.
2. `chord` is the streamwise length of `wing.stl` (1.1597 m for B). It is not a reference chord.
   The four planes sit 0.05, 0.25, 0.50 and 0.75 of that length behind the trailing edge.
3. `|U| for cp` is the Condition CR speed in m/s.
4. `rhoInf 1, pInf 0` is right for an incompressible case, where p is kinematic. The forces in
   `argusForces` are therefore per unit density.
5. For B, the file this writes is byte-for-byte the one the delivered B solve used.
6. A cruise recipe prints `rho mode 'rho'` with its own freestream values: early cruise
   `rhoInf 0.412710, pInf 26495.9`, late cruise `rhoInf 0.310830, pInf 19374.0`. A cruise block
   also samples `U p T rho` on the Trefftz planes, and the install prints
   `removed superseded function objects: yPlus` (section 3, item 2).

Running the script a second time changes nothing. It prints
`controlDict already includes argusPostPro`.

After the install, `sha256sum -c MANIFEST.sha256` prints `system/controlDict: FAILED`. That is
expected. The mesh job does not compare `system/controlDict` with the manifest; it checks the
angle in it instead. `system/argusPostPro` is not in the manifest at all.

Two more points about the check:

1. The script has its own test, which should end with `self-test: PASS, both directions`:

   ```bash
   python3 $REPO/scripts/assert_case_postpro.py --self-test
   ```

2. Do not read the `ok compressible_cp_frame` line as a pass. That check only recognises the
   OpenFOAM-org 7 compressible solver names, and every OpenFOAM-org 12 case says
   `application foamRun`, so it reports `ok` without having checked anything. The cruise blocks
   are nonetheless written with real freestream values (item 6 above); the check simply cannot
   catch a wrong cp frame on an OpenFOAM-org 12 compressible case.

Quote only the `ok`, `MISSING` and `---- N of 8 ----` lines of this script. Some `consumer:`
lines under a `MISSING` name scripts that are not in this repository
(`hpc/make_run_card.py`, `scripts/make_validation_run_card.py` and
`scripts/export_wing_vtk.py`), and some `consumer:` and `why:` lines refer to the source
project's own records, which are not shipped.

**If it does not:**

1. `ModuleNotFoundError: No module named 'numpy'`: this terminal's `python3` has no numpy. Use
   your Python terminal. Nothing was written.
2. `HALT <case>: no sampled wall surface and no readable ASCII STL ...`:
   `constant/triSurface/wing.stl` is missing. Do [04](./04_send_to_hpc.md), Step 5c (or
   [03](./03_build_a_case.md), Step 5), first. Nothing was written.

### 4.2 Send the case and mesh it as usual

Send the case as in [04](./04_send_to_hpc.md) and submit the mesh job as in
[05](./05_submit_and_run.md). The mesh job's two input gates accept a case with the block
installed. In the mesh job's output file you should see these lines (the first line shows your
angle):

```
    alpha 1.299999 deg, |U| 34.0000 m/s, axis lag -2.04e-08 deg
  GATE A: 19 manifest entries = 17 intact + 0 altered + 0 missing + 2 alpha-bearing
  recipe intact
  GATE B constant/triSurface/wing.stl: tris=192588 open=0 over=0 y_min=-20.000mm OK
  surface closed
```

For a cruise recipe the Gate A line reads
`22 manifest entries = 20 intact + 0 altered + 0 missing + 2 alpha-bearing`.

The block does not change the mesh. The tutorial mesh of [03](./03_build_a_case.md) was built
twice on 8 cores, once with the block and once without. Both runs used the mesh job's order of
programs. Both gave 1,193,347 cells and the same four failed checks. The logs differed only in
times and process numbers, and no meshing program wrote a `postProcessing` folder. A 192-rank
production mesh job with the block installed has run on HPC12 (2026-09-24, an early-cruise mesh
of W, 73.7 million cells, the same four failed checks). There was no mesh of the same surface
without the block to compare it with, so at production size "does not change the mesh" is
NOT VERIFIED.

**If it does not:** an `ALTERED` line names a recipe file that no longer matches the manifest;
it is never `system/argusPostPro`. Start again from a fresh copy of the recipe, then repeat
section 4.1.

### 4.3 Stage the solve case and check it (HPC12 login node)

```bash
$REPO/hpc/stage_solve_case.sh $MESH $CASE
python3 $REPO/scripts/assert_case_postpro.py --before $CASE
```

`$MESH` is the finished mesh case and `$CASE` the new solve case, as in
[05](./05_submit_and_run.md), section 0. The staging script copies `system/` whole, so the block
goes with it. On HPC12 only its input checks have been run; its copying, verifying and writing
steps have run only on a local Linux machine: NOT VERIFIED on HPC12
([05](./05_submit_and_run.md), section 6). Under `NEXT STEPS` it prints:

```
1. Post-processing block.
   system/controlDict already includes argusPostPro.
```

The check should end `---- 8 of 8 requirements satisfied ----`. It needs only the Python
standard library. The login node's `python3` (3.6.8) runs it: it gave 8 of 8 on each of the six
delivered Condition CR solve cases.

**If it does not:** `system/controlDict does NOT include argusPostPro` means the mesh case was
built without the block. Do section 4.4 before you submit.

### 4.4 Only if the block is missing: install it on the solve case

Do this before the first leg, at the HPC12 login prompt:

```bash
python3 $REPO/scripts/install_postpro.py $CASE
python3 $REPO/scripts/assert_case_postpro.py --before $CASE
```

What was tested on HPC12:

1. The login node has numpy 1.12.1.
2. The script's derivations ran there with `--dry-run` under Python 3.6.8, on a delivered mesh
   case, and printed the same values as in section 4.1. There the wall patch `wing` is read from
   `constant/polyMesh/boundary`.
3. The writing step has not been run on HPC12: NOT VERIFIED.

The `--before` check afterwards tells you whether the install worked; it must end `8 of 8`. On
your own machine, the same install on a staged solve case wrote a byte-identical block and gave
8 of 8.

**If it does not:** run section 4.1 on your local copy of the case. Copy only its
`system/argusPostPro` into `$CASE/system/` (from your machine, for example with `scp` as in
[04](./04_send_to_hpc.md), Step 6b, item 4). In `$CASE/system/controlDict`, add the line
`    #include "argusPostPro"` directly under the `{` that follows `functions`. Then repeat the
check (NOT VERIFIED on HPC12).

### 4.5 Switch the gate on for every leg

`hpc/solve_hpc12.pbs` runs the same check at the start of each leg, but only if it finds the
checker file. It looks at the path in `ARGUS_GATE`. Without `ARGUS_GATE`, it looks at
`$HOME/argus_hpc12/assert_case_postpro.py`, which does not exist on a new account.

So put the shipped copy in the `-v` list of **every** leg you submit. On HPC12 you submit every
leg by hand: a job's own `qsub` for its next leg fails there with
`could not connect to trqauthd`. That job-made `qsub` would not pass `ARGUS_GATE` on anyway.

```bash
ls $REPO/scripts/assert_case_postpro.py
cd $CASE
qsub -N argus-<name> -q fpt-medium -l nodes=4:ppn=32:typej,walltime=72:00:00 \
     -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py \
     $REPO/hpc/solve_hpc12.pbs
```

1. The `ls` line must print the file's path.
2. `$REPO` must be set in the login shell and must be an absolute path. It is expanded when you
   press Enter, and the job looks for the file only after it has moved into the case folder.
3. For later legs, add `ARGUS_ENDTIME` as in [05](./05_submit_and_run.md), section 10a, and keep
   `ARGUS_GATE`.

Once the job has started, read the verdict. While the job runs, HPC12 keeps its output in
`$HOME/<job number>.hpc12.hpc.OU`. The output file `argus-<name>.o<job number>` appears in the
folder you ran `qsub` from, here `$CASE`, only when the job ends
([05](./05_submit_and_run.md), section 11b):

```bash
grep -E "requirements satisfied|gate NOT run|gate refused" $HOME/<job number>.hpc12.hpc.OU
grep -E "requirements satisfied|gate NOT run|gate refused" $CASE/argus-<name>.o<job number>
```

Use the first line while the job runs and the second once it has ended.

**What you should see:**

```
    ---- 8 of 8 requirements satisfied ----
```

The check runs after the job has started and loaded OpenFOAM. It runs before the write probe,
`decomposePar`, `potentialFoam` and `foamRun`, using the compute node's `python3` (3.6). The
delivered legs ran the same `--before` check there. A refusal still costs the queue wait, which
is why section 4.3 runs the check before `qsub`.

**If it does not:**

1. `NOTE: <path> absent, post-processing gate NOT run`: nothing exists at the path in
   `ARGUS_GATE`. It may be unset, misspelled or relative, or `$REPO` was not set. The leg carries
   on without the check. Fix the `-v` list for the next leg.
2. `---- N of 8 requirements satisfied ----` (N below 8), then
   `REFUSED: post-processing gate refused`: the job stopped before it touched the case. It leaves
   only an empty `log/leg<N>`; see [05](./05_submit_and_run.md), section 10a, point 5. Do
   section 4.4, then resubmit.

How this was tested: the gate lines of the job script were run unchanged on staged cases.

1. With the block: 8 of 8, and the job carries on.
2. Without the block: `1 of 8`, then `REFUSED: post-processing gate refused`, exit status 2.
3. With `ARGUS_GATE` unset or misspelled: the `NOTE` line, and the job carries on.

It has been: on 2026-09-24 a solve leg on HPC12 (an early-cruise case of W, job output file of leg 1) received `ARGUS_GATE` through `-v` (it is in the job's `qstat -f` variable list) and ran the check before its first iteration, ending `8 of 8 requirements satisfied`.

In `log/leg<N>/foamRun.log`, a case with the block shows `wallShearStress argusWallShearStress:`,
then `processing wall patches:` and `wing`, then two `Reading surface description:` blocks. The
first lists `wingSurface`; the second lists `trefftz_x1` to `trefftz_x4`. The block adds no
warnings.

### 4.6 After the run: what `--after` says

```bash
python3 $REPO/scripts/assert_case_postpro.py --after $CASE
```

On a correctly run case, expect **5 of 8**. Each `MISSING` line is followed by two lines of
explanation.

```
  ok      forces_split
  ok      surface_format_carries_topology
  ok      compressible_cp_frame
  ok      trefftz_planes
  ok      residuals
  MISSING wing_surface_p
  MISSING wing_surface_tau_VECTOR
  MISSING wing_surface_yplus
  ---- 5 of 8 requirements satisfied ----
```

**The three `MISSING` lines are wrong.** The check looks for each field name only in the first
20,000 characters of a wing-surface file. OpenFOAM-12 writes the field names after the list of
points and faces, so on any real mesh (the tutorial mesh of [03](./03_build_a_case.md) or
larger) they are beyond that limit. In the delivered B file (628,068,096 bytes), none of the
three names is in the first 20,000 characters. Before it gives up, the check also reads every
wing-surface file in full: 13 files and 7,698 MiB for the delivered B. Check these fields
directly instead:

```bash
cd $CASE
T=$(ls postProcessing/argusWingSurface | sort -n | tail -1)
grep -a -m4 -E "^(p|cp|yPlus|wallShearStress) [0-9]+ [0-9]+ float" postProcessing/argusWingSurface/$T/wingSurface.vtk
ls postProcessing/argusTrefftz/$T
```

**What you should see** (the delivered B solve, where `T` is 6000):

```
p 1 5276772 float
cp 1 5276772 float
yPlus 1 5276772 float
wallShearStress 3 5276772 float
trefftz_x1.vtk  trefftz_x2.vtk  trefftz_x3.vtk  trefftz_x4.vtk
```

The `3` after `wallShearStress` shows that the vector was written. The next number is the number
of wing faces (82,922 on the tutorial mesh).

Other `--after` results:

1. `1 of 8`, with only `compressible_cp_frame` ok: the case ran without the block. A Condition
   CR case then wrote only `forceCoeffs1`; a cruise case also wrote `yPlus` and a `residuals`
   file holding only a time column (section 3, items 2 and 3), and neither of those counts. The
   missing data cannot be produced afterwards.
2. `<case>: no postProcessing/ at all` and `0 of 8`: the case has not run, or you gave the wrong
   folder.
3. `compressible_cp_frame` is always `ok` after a run, because nothing a run writes can show it.
4. `--time` appears in the script's usage text but does not exist. It stops with
   `error: unrecognized arguments`.

The check reads only `postProcessing/`, so it also runs on a copy brought back as in
[04](./04_send_to_hpc.md).

## 5. Forces and the convergence verdict

### 5.1 What to bring back from HPC12

For each finished case, copy these back to your machine (how to copy is in
[04_send_to_hpc.md](./04_send_to_hpc.md)):

1. `.converged` (present only if the case finished and passed the gate, see 5.3);
2. `postProcessing/` (for a case with the post-processing block this is large: for the
   delivered B, 7,698 MiB of wing-surface files and 329 MiB of Trefftz planes, see section 3);
3. `log/` (the solver logs, one subdirectory per leg: `log/leg1/foamRun.log` and so on);
4. `system/`, `constant/momentumTransport` and `0/` (the run-card fields in section 7 are read
   from these);
5. the PBS output file of each leg (`<job name>.o<job id>`), which carries the time-stamped
   progress of the leg. PBS writes it into the directory `qsub` was run from, and it appears
   there when the job ends. [05_submit_and_run.md](./05_submit_and_run.md) runs every `qsub`
   from inside the case (`cd $CASE` first), and on HPC12 you submit every leg yourself (05,
   section 10a), so these files are in `$CASE` on the cluster. A leg you submitted from
   somewhere else, for example your cluster home, has its file there.

Do not copy `processor*` directories. HPC12 scratch (`/home/scratch/$USER/...`) deletes files
older than 50 days with no recovery (stated in `hpc/solve_hpc12.pbs`), so bring results back
promptly. (`hpc/README_HPC.md` describes an earlier setup on another cluster; it is not the
HPC12 route.)

### 5.2 The force file

Each solve leg starts a new subdirectory named after the iteration it started from, for
example `postProcessing/forceCoeffs1/0/`, then `.../4000/` after a restart. Look at the header
of the first one:

```bash
head -9 "$CASE/postProcessing/forceCoeffs1/0/forceCoeffs.dat"
```

For a case built from `recipes/L11_wallResolved_CR`, OpenFOAM-org 12 writes:

```text
# Force coefficients
# liftDir       : (-2.26873230e-02 0.00000000e+00 9.99742610e-01)
# dragDir       : (9.99742610e-01 0.00000000e+00 2.26873230e-02)
# pitchAxis     : (0.00000000e+00 1.00000000e+00 0.00000000e+00)
# magUInf       : 3.40000000e+01
# lRef          : 3.93957000e-01
# Aref          : 6.20462000e-01
# CofR          : (1.34000000e+00 0.00000000e+00 0.00000000e+00)
# Time          	Cm              	Cd              	Cl              	Cl(f)           	Cl(r)
```

What the header tells you:

1. **The frame.** `Aref` 0.620462 m² is the half-model reference area on the DSO basis (full
   wing S_ref 1.24092 m²), so a half-model coefficient here equals the full-wing coefficient on
   the DSO basis. `lRef` 0.393957 m is the reference chord, `magUInf` 34.0 m/s the free-stream
   speed, and `Cm` is taken about `CofR` (x = 1.34 m) around the y axis. Drag counts are
   `Cd` × 10⁴.
2. **The angle of attack the solver actually used**, from `dragDir`:
   alpha = atan2(dragDir_z, dragDir_x). The recipe's own vectors give 1.30 deg; a trimmed case
   carries its own (see [06_trim_to_target_cl.md](./06_trim_to_target_cl.md)). The scripts
   read alpha from this header, not from `controlDict`. On a trimmed case each leg's file has
   its own header, and the surface scripts of section 6.4 read the one of the latest leg.
3. **`Cl(f)` and `Cl(r)` are the front and rear parts of the lift**, not a pressure and viscous
   split. The split needs the `forces` function object (`argusForces` in section 3).
4. In the Condition CR recipe `forceCoeffs1` writes every 5 iterations (`writeInterval 5` in
   the function object), while the fields are written every 500 iterations (`writeInterval 500`
   in `controlDict`).

If the header is missing lines, or the columns differ, the case was not run with
OpenFOAM-org 12 or its `forceCoeffs1` block was changed; check `system/controlDict`.

### 5.3 The convergence verdict: `.converged`

At the end of every leg, `hpc/solve_hpc12.pbs` reads all `forceCoeffs1` files, keeps the later
value where two legs overlap, and applies the convergence gate to the last 200 samples:

1. mean drift of C_D, the difference between the mean of the second 100 samples and the mean
   of the first 100, at most 0.05 counts;
2. mean drift of C_L, the same way, at most 1.0 count;
3. C_L span, the difference between the last and the first sample of the window, at most
   1.0 count.

With samples every 5 iterations, the window is the last 1000 iterations. If the case has
reached its `endTime` and passes all three, the job writes `$CASE/.converged`:

```bash
cat "$CASE/.converged"
```

For the delivered baseline trim (`SWB_trim`) it reads (reconstructed from the `notes` field of
`data/run_cards/SWB_r2_trim.json`, with the line breaks the job writes):

```text
case SWB_trim
finalised 2026-09-21T18:16:17Z on n12-142, leg 3
  CONVERGENCE GATE over 200 samples: MEAN drift Cd 0.0003 ct (<=0.05), Cl 0.0399 ct (<=1.0), Cl SPAN 0.0777 ct (<=1.0) -> CONVERGED
  REPORT THESE: Cd 0.0188031  Cl 0.4282682  (window means)
```

The `REPORT THESE` values are the 200-sample window means, and they are the numbers the report
quotes. **No `.converged` file means the case is not a result**: it has not finished, or it
finished without passing the gate. Do not quote the last line of `forceCoeffs.dat` instead.

### 5.4 One force history per case (optional)

To plot C_D and C_L across all restart legs, `scripts/splice_postpro.py` joins the
`forceCoeffs1` pieces of every case in a directory into one file per case. Run it on **copies**
of your cases, because it writes into them. Give the directory as an absolute path, because
the command runs from `$REPO`:

```bash
cd "$REPO"
python3 scripts/splice_postpro.py <absolute path of a directory holding one directory per case>
```

Expected: a table with one line per case (`case samples iter Cd (ct) Cl dCd dCl gate`), a file
`<case>/forceCoeffs_spliced.dat` in each case, and `INDEX.json` in the directory you named. Its
own verdict column is **not** the project's gate: it has no C_L span test, it calls a case
quotable only if it reached iteration 4500 (`ENDTIME = 4500`, a fixed number in the script), and
its `Cd (ct)` column is the last iteration, not the window mean. The Condition CR recipe now ends
its first leg at 4000, so a case that converged there is listed as not quotable. The verdict
that counts is `.converged`.

### 5.5 Where the delivered numbers are

`data/rans_forces.json` holds every delivered converged case under `cases`, and the trimmed
result of each geometry at each condition under `trims`. Of the `trims` rows, only M's
Condition CR row carries a `mesh_generation` label. For every row, take the label from the
`cases` record that the row names in `trim_case`, as this snippet does. It prints the
Condition CR rows with their labels:

```bash
cd "$REPO"
python3 - <<'EOF'
import json
d = json.load(open("data/rans_forces.json"))
t = d["trims"]["condition_CR"]
print("target C_L", t["target_CL"], " Mach", t["mach"])
for g, r in sorted(t["geometries"].items()):
    c = d["cases"][r["trim_case"]]
    print(g, r["trim_case"], c["mesh_generation"],
          "C_D %.3f ct" % r["Cd_counts"],
          "dC_D %+.3f ct" % r["delta_Cd_counts_vs_baseline"])
EOF
```

**What you should see:**

```text
target C_L 0.428277635108  Mach 0.1
B SWB_trim wake_refined C_D 188.031 ct dC_D +0.000 ct
C SWC_trim wake_refined C_D 187.840 ct dC_D -0.191 ct
F SWF_trim wake_refined C_D 187.216 ct dC_D -0.815 ct
H SWH_trim wake_refined C_D 188.022 ct dC_D -0.009 ct
M SWM_trim wake_refined C_D 188.263 ct dC_D +0.232 ct
W SWW_trim wake_refined C_D 187.163 ct dC_D -0.868 ct
```

These are the values of section 9.2 (Condition CR, M 0.10, target C_L 0.428277635108, DSO
basis). **If it does not** show six `wake_refined` lines with these values, your `data/` is not
the delivered one.

Take results from the `trims` rows only:

1. Each row keeps the value it replaced under `superseded_trim_entry`, and M's row also keeps
   its earlier published-mesh result under `superseded_published_row` (C_D 191.246 ct). These
   are superseded: never quote them.
2. Inside `cases`, the six `SW?_trim` records labelled `wake_refined` take their values from
   the wake-refined case: `Cd`, `Cd_counts`, `Cl`, `alpha_deg`, `iterations` and `endTime`
   (6000), `marker_text`, the convergence fields (`drift_Cd_ct`, `drift_Cl_ct`, `Cl_span_ct`,
   equal to that marker's), and the averaging window (`window_t_first` 5005 to
   `window_t_last` 6000, 200 samples). `alpha_deg_resid_vs_provenance` is null on purpose, with
   a note. In each `sources` map three entries, `geometry_sha256`, `provenance_alpha_deg` and
   `turbulence_model`, still cite the published case archive, as the map's own
   `_not_rederived_on_r2` key says. `fields_replaced_at_r2_repair` holds the values each record
   carried just before its 2026-09-24 repair and `superseded_published_sources` the earlier
   source map; neither is the published result, so never quote them.
3. `data/rans_forces_current.json` has the shape of the raw harvest, carries no
   `mesh_generation` labels at all, and holds earlier Condition CR values for C, F, H and M
   (for example M `Cd` 0.0191246). Do not quote from it.

## 6. Wake, surface and y+ data

### 6.1 Spanwise cuts: `planeSample`

> **These commands need the mesh and the solved fields, and a copy pulled as in
> [04_send_to_hpc.md](./04_send_to_hpc.md) has neither.** On HPC12 the solve job leaves the
> solved fields only inside `processor*/`, one piece per rank, and the pull in 04 leaves
> `processor*/` and `constant/polyMesh/` behind. On such a copy the commands stop at once with
> `--> FOAM FATAL ERROR:` and
> `Cannot find file "points" in directory "polyMesh" in times "0" down to constant` (tested for
> this guide on a copy laid out as the 04 pull). On a case that has its mesh but only the `0/`
> time directory, they do not warn you. `-latestTime` then means `0`, and `planeSample` samples
> the starting fields and finishes normally (seen on a recipe case: it ran at time 0 and wrote
> seven files). Before trusting any output, check that the time directory it was written under
> (`postProcessing/planeSample/<time>/`) is the final iteration of the run, for example `6000`,
> and not `0`.

In a case directory that holds its mesh and the solved fields as a normal time directory, with
OpenFOAM-org 12 loaded:

```bash
cd "$CASE"
foamPostProcess -func planeSample -latestTime
```

Expected: the log lists the seven surfaces `y_eta0p30` to `y_eta0p95`, reads `p` and `U`, and
ends with `End`. The output is `postProcessing/planeSample/<time>/y_eta0p30.vtk` and so on,
seven files, which open in ParaView. A warning
`Removing patchGroup 'symmetry' which clashes with patch 0 of the same name` is harmless.

If it stops instead, the message says which of three things is wrong:

1. `FOAM FATAL ERROR: cannot find file ".../system/controlDict"`: you are not in a case
   directory. `cd "$CASE"` and run it again.
2. `FOAM FATAL IO ERROR: Cannot find configuration file planeSample`: you are in the case, but
   `system/planeSample` is missing. Copy it from the recipe the case was built from.
3. `FOAM FATAL ERROR:` followed by `Cannot find file "points" in directory "polyMesh" ...`: the
   case has no mesh, for example a copy pulled as in 04. Run the step where the mesh and the
   solved fields are: on HPC12, the parallel form below.

Where the fields are still decomposed (in `processor*/`, as they are on HPC12), add `-parallel`
and start it with `mpirun` using the case's own rank count (`numberOfSubdomains` in
`system/decomposeParDict`). On HPC12, load OpenFOAM-org 12 first, as in section 2, item 2.
Then check that `mpirun` is the Open MPI that OpenFOAM uses:

```bash
mpirun --version | head -1
```

Expected on a workstation set up as in [01_setup.md](./01_setup.md): `mpirun (Open MPI) 4.0.3`.
On HPC12, after the three lines of section 2, item 2: `mpirun (Open MPI) 4.1.2`.
If it prints `HYDRA build details:`, an Anaconda `mpirun` comes first on your `PATH`, and the
parallel command below fails with every rank printing
`attempt to run parallel on 1 processor`; fix it as in [01_setup.md](./01_setup.md),
section A5. Then:

```bash
mpirun -np <ranks> foamPostProcess -parallel -func planeSample -latestTime
```

This was checked on a 2-rank decomposition on a workstation, with the system's Open MPI.
Running it inside an HPC12 job is **NOT VERIFIED**, and this repository has no job script for
it.

For Condition CR, `p` in these files is kinematic pressure (pressure divided by density,
m²/s²), because the solver is incompressible.

### 6.2 y+ and wall shear stress after the run

A case that carried the post-processing block already has both: `argusYPlus` holds the y+
minimum, maximum and average on `wing` at every saved iteration, and the wing-surface file
carries y+ and the wall shear stress vector on every face (section 3). The rest of this section
is for a case without the block. The Condition CR recipe alone does not compute y+ during the
solve.

You can compute y+ and the wall shear stress from the last written time. The warning at the
top of section 6.1 applies here too: this needs the mesh and the solved fields, and the time
directory the output is written under must be the final iteration, not `0`.

```bash
cd "$CASE"
foamPostProcess -solver incompressibleFluid -func yPlus -latestTime
foamPostProcess -solver incompressibleFluid -func wallShearStress -latestTime
```

Expected: each ends with `End`. The yPlus run prints a line
`patch wing y+ : min = ..., max = ..., average = ...` and writes
`postProcessing/yPlus/<time>/yPlus.dat` (columns `Time patch min max average`); both write
their field into the time directory, for viewing in ParaView.

The `-solver incompressibleFluid` option is required. Without it the yPlus run stops with
`FOAM FATAL ERROR: Unable to find turbulence model in the database`. The cruise recipes already
compute y+ during the solve.

Where the fields are still decomposed (in `processor*/`, as they are on HPC12), run the same two
commands in parallel, with the case's own rank count and the Open MPI `mpirun` checked in
section 6.1 (on HPC12, load OpenFOAM-org 12 first, as in section 2, item 2):

```bash
cd "$CASE"
mpirun -np <ranks> foamPostProcess -parallel -solver incompressibleFluid -func yPlus -latestTime
mpirun -np <ranks> foamPostProcess -parallel -solver incompressibleFluid -func wallShearStress -latestTime
```

This form was checked on a 2-rank decomposition of an OpenFOAM-12 tutorial case on a
workstation. Both ended with `End`, the yPlus run printed one
`patch <name> y+ : min = ..., max = ..., average = ...` line per wall patch, and the files were
written under `postProcessing/yPlus/<time>/` and `postProcessing/wallShearStress/<time>/`.
Running it inside an HPC12 job is **NOT VERIFIED**, and this repository has no job script for
it.

What to expect for a wall-resolved Condition CR result: the delivered baseline trim reports y+
minimum 0.0936, mean 1.24 and maximum 81.9 on the wing (`mesh.yplus` in
`data/run_cards/SWB_r2_trim.json`). The mean near 1 is the design target; the high maximum is a
small patch in the log layer and is a caveat on surface friction, not on the far-field results.

### 6.3 Trefftz-plane induced drag and span efficiency

The delivered induced drag C_Di and span efficiency e come from the four `argusTrefftz` planes
(section 3):

1. C_Di = ∫(v² + w²) dS / (U² A_ref), where v is the spanwise velocity and w is the vertical
   velocity minus the free-stream part U sin(alpha). A_ref is the half-model 0.620462 m²,
   matching `forceCoeffs1`, because the planes cover the half model.
2. e = C_L² / (π AR C_Di), with AR 10.7807 on the DSO basis (b_ref 3.6576 m,
   S_ref 1.24092 m²) and C_L the case's own.
3. All four planes are reported, never one. C_Di falls with distance behind the wing and has
   no plateau: for the baseline Condition CR trim on the wake-refined mesh it drops by 9.46 %
   between the first plane (x = 2.394 m) and the last (x = 3.206 m), and by 9.46 to 10.22 %
   for the six wake-refined geometries in `data/span_efficiency.json`. The decay is numerical
   (the grid dissipating the wake), so it depends on the mesh. That spread is the
   uncertainty. So the absolute value of e depends on which plane, and which mesh, it was taken
   on: compare e between geometries only **at the same plane on the same mesh generation**, and
   never compare an e from one plane, mesh or code with an e from another.

What can be done from this repository:

1. **Sampling the planes: yes.** A case with the post-processing block writes them during the
   solve (section 4).
2. **Integrating the planes of a new case: no shipped route.** The routine that does it,
   `span_efficiency()` in `scripts/case_derived_quantities.py`, only looks in
   `results/postpro_latest/<case>/` (a layout made by a result-pulling script that is not
   shipped), reads the case's `controlDict` from `cases/of12/v5/solve/<case>/`, and writes to
   `results/derived/`. None of these directories exists here. Run with `--all` on a clone of this
   repository it finds no case and prints `error: name a case or pass --all`. The two formulas
   above are what it computes.
3. **Rebuilding `data/span_efficiency.json` from the delivered plane integrals: yes.** The
   routine's output for the six delivered trims ships as `data/span_eff_r2trim_raw.json`, one
   record per plane. `report/scripts/build_span_efficiency.py` selects from it, and reproduces
   the shipped file byte for byte, if you tell Python where `case_derived_quantities.py` is:

   ```bash
   cd "$REPO"
   PYTHONPATH="$REPO/scripts" python3 report/scripts/build_span_efficiency.py \
       --raw data/span_eff_r2trim_raw.json --out "$HOME/span_efficiency_check.json"
   cmp data/span_efficiency.json "$HOME/span_efficiency_check.json" && echo identical
   ```

   **What you should see:** `geometries: 6 of 6 present under tag 'r2trim'`, one line per
   geometry, a `wrote` line, a table (for example
   `M  CDi  62.035 ->  55.696 ct  (decay 10.22%)   e_last 0.9724`), and last `identical`.

   **If it does not:**
   1. `ModuleNotFoundError: No module named 'case_derived_quantities'`: `PYTHONPATH` was not
      set. The script looks for `scripts/` three directory levels above `report/`.
   2. `FileNotFoundError` naming `report/data/span_efficiency.json`: `--out` was left out. Its
      default is a directory that does not exist in your clone. Always give `--out`, and never
      point it at `data/span_efficiency.json`.

Section 10 lists the other wake scripts and why they do not help here.

### 6.4 Spanwise loading from a wing-surface file

`scripts/spanwise_dense_report.py` is the script behind the spanwise-loading figure. It needs a
wing-surface VTK that carries `p` and `cp` (an `argusWingSurface` file, section 3), and the
case directory, from which it reads `dragDir` in `forceCoeffs1`. Only a case that carried the
post-processing block has such a file.

Use the wing-surface file of the final iteration. The `argusWingSurface` folder also holds one
for iteration 0, with the starting fields, which the script cannot use. The `T=` line picks the
highest iteration number, and `echo` should print the last iteration of the run (6000 for the
delivered trims), not `0`:

```bash
cd "$REPO"
T=$(ls "$CASE/postProcessing/argusWingSurface" | sort -n | tail -1)
echo "$T"
python3 scripts/spanwise_dense_report.py \
    --vtk "$CASE/postProcessing/argusWingSurface/$T/wingSurface.vtk" \
    --case <case name> --case-dir "$CASE" --time "$T" --out <output directory>
```

The script creates the output directory if it does not exist. This is what it printed for the delivered M trim at iteration 6000 (shipped as
`data/spanwise/SWM_trim.log`; that run asked for 78 stations and gave `--cl-ref`, and the path
in the first line is shortened here):

```text
reading <case>/postProcessing/argusWingSurface/6000/wingSurface.vtk
  5276994 faces, fields: cp, p, wallShearStress, yPlus
  DERIVED FROM THIS FILE: q = 578, p_inf = 9.37052e-08 (resid max 5.78e-06, p99 1.21e-06)
  dragDir [ 0.999759193  0.           0.021944404] -> alpha 1.257423 deg (from forceCoeffs.dat)
  semispan from surface 1.825395 m (registered 1.828800)
  stations: 78 asked = 78 integrated + 0 skipped
  CL from strips (2x half-wing / Sref 1.24092) = 0.408063
  CL from forceCoeffs = 0.428278, closure gap -4.720%
  wrote SWM_trim_spanwise_dense.csv and SWM_trim_spanwise_dense.json
```

1. `q = 578` is the kinematic dynamic pressure 0.5 × 34.0² m²/s² of Condition CR, derived from
   the file's own p against cp. A different value means the file is not from a Condition CR
   case.
2. The angle is read from the `dragDir` of the **latest** leg's `forceCoeffs1` directory (the
   directories are sorted as numbers), so on a trimmed case it is the trimmed angle: 1.257423
   deg for M, its `trim.alpha_trim_deg` 1.25742268 to the printed precision. Compare the printed
   angle with your case's trimmed angle. For a case directory holding the baseline's two
   `dragDir` vectors, in `0/` (the cold leg, 1.742369 deg) and in `4000/` and `5000/` (trimmed),
   the routine returned 1.747609 deg.
3. The strip C_L carries the half-wing factor of two and the full-wing S_ref 1.24092 m² (DSO
   basis). Add `--cl-ref <window-mean C_L>` to have the script compare it with the case's own
   C_L. The strips cover only part of the span (eta 0.025 to 0.988 in that run; 0.02 to 0.99
   with the defaults of item 4) and carry pressure only, so a strip C_L about 5 % below the
   case's C_L is expected; the gap is an upper bound on the sampling error.
4. `--n` (default 63), `--eta-min` (0.02), `--eta-max` (0.99) and `--band` (0.004 of the
   semispan) set the stations.
5. The report's figure script reads these outputs from `data/spanwise/`, which ships M's files
   only (section 8.3).

**If it does not:**

1. `FATAL: cp is constant over the surface; cannot derive q`: the file is the iteration-0
   surface, which holds the starting fields. Use the final iteration, as above.
2. `FATAL: too few finite (p, cp) pairs to derive q`: the file has fewer than 1,000 faces
   with finite p and cp, so it is not the wing surface of a real mesh (the tutorial mesh of
   [03](./03_build_a_case.md) has 82,922 wing faces).

`scripts/section_slices_from_surface.py` reads the angle the same way; `spanwise_dense_report.py`
imports that routine from it.

## 7. Run cards

A run card is one JSON file per completed run recording what was run, on what, and what came
out. The six delivered cards, one per Condition CR trim, are in `data/run_cards/`. **A card that
does not validate may not be quoted.**

### 7.1 What writes a card

The six cards were assembled by a card builder in the source repository from facts collected on
HPC12 by a cluster-side script. The shipped tooling names them `scripts/build_r2_run_cards.py`
and `collect_run_facts.sh`. **Neither is in this repository**, so you cannot assemble a
complete card for a new run from the shipped tools alone. Ask the project for them (see
section 11).

What is shipped:

1. `scripts/derive_run_card_fields.py`: reads the fields that must come from the case itself
   (section 7.2);
2. `scripts/run_card.schema.json`: the schema every card must satisfy;
3. `scripts/validate_run_card.py`: checks a card against the schema, the rules between fields
   and the registry (section 7.5).

### 7.2 Derive the case-dependent fields

```bash
cd "$REPO"
python3 scripts/derive_run_card_fields.py --case "$CASE"
```

It reads, from the case: `application` and `solver` in `system/controlDict`, cross-checked
against the `Exec` line of the solver log; the OpenFOAM version from the log's `Build` line;
the turbulence model from `constant/momentumTransport`; the wall treatment from the `wing`
boundary conditions in `0/nut`, `0/k` and `0/omega`; and the `residualControl` targets from
`system/fvSolution`. The solver log is looked for as `log.foamRun`, `log/foamRun.log`,
`logs/foamRun.log` or `log/leg<N>/foamRun.log` (the HPC12 layout); the newest one wins.

For a case built from `recipes/L11_wallResolved_CR` and run with OpenFOAM-org 12 the output is:

```json
{
  "model": {
    "solver": "foamRun",
    "turbulence_model": "kOmegaSST",
    "wall_treatment": "low-Re resolved",
    "solver_module": "incompressibleFluid"
  },
  "provenance": {
    "openfoam_version": "openfoam-org-12",
    "openfoam_build": "12-86e126a7bc4d"
  },
  "numerics": {
    "convergence_criteria": {
      "p": 1e-09,
      "U": 1e-10,
      "(k|omega)": 1e-10
    }
  },
  "host_in_log": "<your-host>",
  "wall_nut_bc": {
    "wing": "nutLowReWallFunction"
  },
  "wall_patches": [
    "wing"
  ],
  "wall_functions_sanctioned_by": null,
  "wall_treatment_block": {
    "nut": "nutLowReWallFunction",
    "k": "kLowReWallFunction",
    "omega": "omegaWallFunction"
  }
}
```

The build string after `12-` depends on your OpenFOAM-org 12 build. `nutLowReWallFunction` is,
despite its name, the boundary condition for a **resolved** wall (it sets the turbulent
viscosity to zero at the wall), which is why the treatment reads `low-Re resolved`.

**When it stops instead.** The script refuses a case that is off its recipe rather than
recording it. Two examples, produced deliberately on a copy of the test case:

```text
DerivationError: <case>: turbulence model kEpsilon is off the incompressible wall-resolved recipe, which admits kOmegaSST (...). Changing it needs a decision entry.
DerivationError: <case>: nut BC on patch(es) ['wing'] is ['nutUSpaldingWallFunction'], i.e. WALL FUNCTIONS, and the incompressible wall-resolved recipe does not sanction them. ...
```

"Needs a decision entry" refers to the project's own record of approved changes to the setup,
which is not shipped. For you it means only this: the case does not match the recipe, and the
script will not describe it as if it did.

If you see one of these on a Condition CR case, the case is not the delivered setup: compare its
`constant/momentumTransport` and `0/nut` with the recipe. For the cruise cases (solver `fluid`)
wall functions are part of the recipe and are accepted. A message saying no solver log was
found means the log is not in one of the places listed above.

**Do not use `--audit-all` in this repository.** It audits cards found under `cases/**` and
`results/run_cards_r2/`, neither of which exists here, so it finds 0 cards and still prints
`PASS`. That is not a check.

### 7.3 What is in a card: `data/run_cards/SWB_r2_trim.json`, field by field

The "class" column comes from the schema the cards are validated against,
`scripts/run_card.schema.json` (section 7.5). **derived** fields are recomputed from the case
whenever a card is regenerated; **recorded at run time** fields describe the run itself and must
never be recomputed, because regenerating them would overwrite the record of the run with the
circumstances of the regeneration. A blank class means the schema does not class that field. In
that schema, 21 fields are classed derived and 9 recorded at run time.

| Field | Value in `SWB_r2_trim.json` | Meaning | Class |
|---|---|---|---|
| `schema_version` | 1 | card format version | |
| `candidate` | `baseline` | registry key in `registry/candidates.yaml` (letter B) | |
| `case.path` | `/home/scratch/$USER/argus/solve_r2/SWB_trim` | where the case ran on HPC12 | |
| `case.dimensionality`, `case.condition`, `case.eta` | `3D`, `CR`, null | a 3D wing at Condition CR (eta is for 2D sections) | |
| `case.description` | "Wake-refined re-run of the SWB CR trim, ..." | free text; this is where the card says it is on the wake-refined mesh | |
| `geometry.file` | `geometry/derived/wing3d/baseline_oml_placed_laddercap.stl` | meshing surface, a path in the source repository (you build the same surface with [02](./02_geometry_from_openvsp.md)) | |
| `geometry.sha256` | `20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641` | checksum of that surface; the validator checks it against the registry | |
| `provenance.git_commit`, `provenance.git_dirty` | `e1648ab3549643cb383df377b13a37ed2facf27c`, true | commit of the **source** repository when the card was written, not a commit of this repository (the notes say so) | recorded at run time |
| `provenance.openfoam_version` | `openfoam-org-12` | | recorded at run time |
| `provenance.date` | `2026-09-21` | | recorded at run time |
| `mesh.generator` | `snappyHexMesh` | | |
| `mesh.cells` | 99111506 | cell count | derived |
| `mesh.max_skewness` | 15.844881 (dimensionless) | maximum skewness, as `checkMesh` reports it | derived |
| `mesh.max_non_orthogonality` | 70.957803 (degrees) | maximum non-orthogonality, as `checkMesh` reports it. The six cards give 70.957803 (B) to 71.967297 (W) | derived |
| `mesh.yplus.min`, `.mean`, `.max` | 0.0935698388, 1.24115357, 81.880779 | y+ on `wing` from the yPlus function object at the final iteration | |
| `mesh.checkMesh_pass` | false | `checkMesh` flags 4 checks (maximum skewness, face tets, cells with small determinant, concave cells), the same 4 as on the published mesh (notes) | derived |
| `mesh.mesh_generation` | ranks 192, method `scotch`, inferred false, machine `hpc12` | **how the mesh was built** (rank count and decomposition method). This is **not** the `wake_refined`/`published` label of section 9 | |
| `model.solver`, `model.turbulence_model`, `model.wall_treatment` | `foamRun`, `kOmegaSST`, `low-Re resolved` | as derived in 7.2 (the shipped cards do not carry `solver_module`) | derived |
| `bc.farfield`, `bc.airfoil` | inletOutlet on the outer patches with the free-stream turned to the trim alpha; noSlip on `wing` | free-text summary | |
| `numerics.schemes_summary` | steady SIMPLE, second-order upwind, GAMG on p | free text | |
| `numerics.convergence_criteria` | p 1e-09, U 1e-10, (k\|omega) 1e-10 | the `residualControl` targets in `system/fvSolution`. The final residuals did not reach them; the run ended at its `endTime` and the force gate decided convergence | derived |
| `convergence.iterations` | 6000 | | derived |
| `convergence.final_residuals` | p 9.16180477e-08, U 2.70476012e-08, k 1.11472875e-07, omega 9.14216084e-10 | | derived |
| `convergence.converged`, `convergence.criterion` | true; the 200-sample gate of section 5.3 | | derived (converged) |
| `conditions.U_mag_m_s`, `alpha_deg` | 34.0, 1.747608659 | free-stream speed and the trimmed angle of attack | |
| `conditions.nu_m2_s`, `Re_ref`, `mach_nominal`, `rho_ref_kg_m3` | 1.46e-05, 917434.1 (rounded), 0.1, 1.225 | Reynolds number on `lRef` 0.393957 m | |
| `conditions.reference` | Aref 0.620462, lRef 0.393957, CofR (1.34, 0, 0), convention text | the frame of every coefficient in the card: DSO basis, half model | |
| `results.CL`, `results.CD`, `results.Cm` | 0.428270799, 0.0188030265, -0.459200543 | force coefficients on the frame above. **Not** the window means: compare `Cd 0.0188031` in the notes | derived |
| `results.CD_pressure`, `results.CD_viscous` | 0.00868908, 0.01011394 (rounded) | pressure and viscous drag, from `argusForces` (section 3); they add up to `results.CD` | derived |
| `results.CD_counts` | 188.030265 | `results.CD` × 10⁴ | |
| `trim.target_CL`, `trim.tolerance` | 0.428277635108, 0.0001 | the Condition CR target and the trim tolerance in C_L (1 count) | |
| `trim.alpha_trim_deg`, `trim.history`, `trim.converged` | 1.747608659; (1.742369 deg, C_L 0.427826554), (1.747608659 deg, C_L 0.428270799); true | the two converged points of the trim | |
| `cost.wall_clock_s`, `n_cores`, `core_hours`, `machine` | 70443, 128, 2504.6, `n12-142` | summed over the 3 legs; core hours = wall clock × cores / 3600. `machine` is an HPC12 compute node on all six cards; the cluster is in `hpc.machine` | recorded at run time |
| `decomposition.method`, `decomposition.n_subdomains` | `hierarchical`, 128 | how the **solve** was split. The six solve cases on HPC12 carry hierarchical (8 2 8), which is what `hpc/stage_solve_case.sh` writes by default. The mesh was built on 192 ranks with scotch (`mesh.mesh_generation`) | recorded at run time (n_subdomains) |
| `hpc` | machine `hpc12.tudelft.net`, solve_ranks 128, `hierarchical`, walltime_requested `72:00:00` | cluster record; must agree with `decomposition` | |
| `wall_treatment` | nut `nutLowReWallFunction`, k `kLowReWallFunction`, omega `omegaWallFunction`, y+ as above, convention `cell-centre` | the wall boundary conditions actually used | derived (nut, k, omega) |
| `notes` | free text | job output file names and times of the three legs, and the `.converged` text quoted in 5.3. It also names a document of the source repository that is not shipped | |

The frame of every coefficient in the card is stated in `conditions.reference`: DSO basis,
half-model A_ref 0.620462 m², lRef 0.393957 m, moments about x = 1.34 m, Condition CR
(M 0.10, target C_L 0.428277635108).

### 7.4 The six delivered cards

| Letter | Card | `candidate` |
|---|---|---|
| B | `data/run_cards/SWB_r2_trim.json` | `baseline` |
| C | `data/run_cards/SWC_r2_trim.json` | `cte_i002_c04` |
| F | `data/run_cards/SWF_r2_trim.json` | `cfft_b02_c01` |
| H | `data/run_cards/SWH_r2_trim.json` | `chc_g02_c06` |
| M | `data/run_cards/SWM_r2_trim.json` | `mcv2_i002_c01` |
| W | `data/run_cards/SWW_r2_trim.json` | `cffw_b01_c01` |

All six are Condition CR trims on the wake-refined mesh. There are no cards for the early- and
late-cruise results.

This mapping is taken from each card's `candidate` and `geometry.sha256` fields (the validator
checks the checksum against the registry, section 7.5), and it agrees with the geometry table
of the report. **Do not take it from the free-text `STATUS` notes in
`registry/candidates.yaml`**: for three candidates they name the wrong case and card. The note
under `chc_g02_c06` (H) names `SWW_trim` and `SWW_r2_trim.json`, the note under `cfft_b02_c01`
(F) names `SWH_trim` and `SWH_r2_trim.json`, and the note under `cffw_b01_c01` (W) names
`SWF_trim` and `SWF_r2_trim.json`. The checksums in the same entries are right and match the
cards; only the notes are wrong.

### 7.5 Validate a card

Run this on your own machine, with the Python of section 2, item 3. The HPC12 login node has no
jsonschema.

```bash
cd "$REPO"
sha256sum scripts/run_card.schema.json
python3 scripts/validate_run_card.py --self-test
python3 scripts/validate_run_card.py data/run_cards/*.json
```

**What you should see:**

```text
1f198e7f3ffea91e3c1cbb9a8f6a51a22e239fb5a0bc7907334fb453bcfcd9e2  scripts/run_card.schema.json
self-test PASS (example accepted; 4 corrupted variants rejected)
PASS data/run_cards/SWB_r2_trim.json
PASS data/run_cards/SWC_r2_trim.json
PASS data/run_cards/SWF_r2_trim.json
PASS data/run_cards/SWH_r2_trim.json
PASS data/run_cards/SWM_r2_trim.json
PASS data/run_cards/SWW_r2_trim.json
```

Each command exits with status 0. The card values this guide quotes, in sections 6.2, 7.3 and
9.2, rest on this check.

**If it does not:**

1. A different checksum on the first line: you have a different version of the schema, and
   the results above may not hold for it. Ask the project which version you have before relying
   on a `PASS`.
2. `AttributeError` or `ModuleNotFoundError`: this Python lacks `jsonschema` 4 or `yaml`
   (section 2, item 3).
3. A `FAIL` line: it is followed by the reason. Fix the card's source, not the card.

What the validator checks:

1. the card against the schema: every required block present, no unknown fields, the right
   types;
2. rules between fields: `hpc.solve_ranks` must equal `decomposition.n_subdomains` and the two
   decomposition methods must agree; a mesh built serially must say so;
3. against `registry/candidates.yaml` (skipped with `--no-registry-check`): the `candidate` must
   be registered, and `geometry.sha256` must be registered for that candidate or for the
   baseline.

Deliberately broken copies of `SWB_r2_trim.json` were rejected as expected, one per rule, each
with exit status 1: a zeroed `geometry.sha256`
(`geometry sha256 000000000000... not registered for 'baseline' (or baseline)`), a removed
`results.CD_viscous` (`results: 'CD_viscous' is a required property`), `hpc.solve_ranks` set to
192 (`hpc/solve_ranks 192 disagrees with decomposition/n_subdomains 128; ...`),
`hpc.solve_decomposition_method` set to `scotch`
(`hpc/solve_decomposition_method 'scotch' disagrees with decomposition/method 'hierarchical'; ...`),
and a removed `mesh.mesh_generation` (`mesh: 'mesh_generation' is a required property`).

The rule still holds for a new run: do not quote a card you have not been able to validate.

## 8. From data to figures and the report

### 8.1 How the delivered files were produced

In order, with what each step needs and whether it runs from this repository:

1. **During the solve, on HPC12**: `forceCoeffs1` and the post-processing block wrote
   `postProcessing/` (section 3). Both ship: the block's template is
   `hpc/templates/argusPostPro`, installed by `scripts/install_postpro.py` (section 4).
2. **End of each leg, on HPC12**: `hpc/solve_hpc12.pbs` applied the gate and wrote `.converged`
   (section 5.3). Shipped.
3. **On HPC12**: `scripts/harvest_rans_forces.py` read every case's `.converged` and
   `forceCoeffs1` and printed one JSON document (the shape of `data/rans_forces_current.json`).
   It takes no arguments. It reads two roots: `$ARGUS_SOLVE` if that is set, else
   `/home/scratch/$USER/argus/solve` (with your user name filled in), and `~/argus_archive`.
   If either is missing it stops with a `FileNotFoundError` naming it. So as it stands it does
   not run on your own machine, or on an HPC12 account without an `~/argus_archive`. With
   `ARGUS_SOLVE` set to a directory holding one directory per case, and an existing
   `~/argus_archive` (an empty one will do), it runs and prints the document. Unless
   `ARGUS_SOLVE` points there, it never looks in `solve_r2`, where the wake-refined Condition CR
   trims ran.
4. **Locally**: `report/scripts/harvest_case_provenance.py [CASE_ROOT] [-o OUT.json]` read each
   archived case's provenance and wrote `data/case_provenance_current.json`. It runs if you give
   both a case root (a directory holding one directory per case) and `-o`; its defaults point
   at a drive of the original workstation and at `report/data/`, which do not exist here.
5. **Locally**: `report/scripts/build_rans_forces.py <raw harvest> <provenance>` merged steps 3
   and 4 into `rans_forces.json`. **Do not run it to "refresh" `data/`**. Rebuilding from the
   shipped inputs was tried: the result has the same 75 cases, but the C, F, H and M Condition
   CR trims come back with superseded values (C 190.857, F 187.257, H 188.065, M 191.246 ct),
   M's ΔC_D (+3.215 ct) is then taken against the wake-refined B, which is a difference across
   mesh generations, and no case carries a `mesh_generation` label. It writes to `data/` next
   to the `scripts/` directory it sits in. In your clone that is `report/data/`, which does not
   exist, so it stops with a `FileNotFoundError` naming `report/data/rans_forces.json.tmp`. **In
   the build directory of section 8.2 it overwrites `$BUILD/data/rans_forces.json`**, and every
   figure you make there afterwards is drawn from the regressed file.
6. **Locally**: `report/scripts/update_cr_from_r2.py` patched the Condition CR trims with
   wake-refined values that are written into the script. They predate the final trims (its H
   entry is 188.065 ct against the current 188.022 ct). It writes to the same place as item 5,
   so in the build directory it also overwrites `$BUILD/data/rans_forces.json` (and
   `rans_forces_current.json`). **Do not run it.** The `mesh_generation` labels were added
   afterwards by a step whose script is not in this repository. If you ran either script in the
   build directory by mistake, delete `$BUILD/data` and copy `$REPO/data` there again.
7. **Locally**: `report/scripts/reconcile_cr_derived.py` recomputes fields derived from others
   in the same record (C_L minus target, the trim tolerance flag, L/D). It only reports unless
   given `--apply`.
8. **Wake data**: the planes were pulled from HPC12 by a script that is not shipped and
   integrated by `span_efficiency()` in `scripts/case_derived_quantities.py`; its output for the
   six trims is shipped as `data/span_eff_r2trim_raw.json`, and
   `report/scripts/build_span_efficiency.py` wrote `data/span_efficiency.json` from it. The last
   step reproduces byte for byte here (section 6.3); the integration of a new case's planes has
   no shipped route.
9. **VLM data**: `report/scripts/extract_vlm.py` wrote `data/vlm.json` from the design study's
   VSPAERO files. Its root is `$ARGUS_ROOT` if that is set, else your clone, and it needs the
   design study's files there in a folder `dso_reference/`, which is not shipped. Without it the
   script stops with `extract_vlm.py: <root>/dso_reference is not a directory. Set ARGUS_ROOT
   to the directory that holds dso_reference/.` and exit status 1. Given that folder, it writes
   `$ARGUS_ROOT/docs/report/all_geometry_2026-09-15/data/vlm.json` (not `data/vlm.json`), and
   that file was byte-identical to the shipped `data/vlm.json` (23 runs).
10. **Figures**: `report/fig/make_*.py` (report) and `report/fig/deck_*.py` (presentation
    slides).
11. **Tables and generated text**: `report/scripts/make_cr_table.py`, `make_rank_table.py`,
    `spaneff_prose_values.py` and `spanwise_prose_values.py` print or write LaTeX from the
    data; `make_q_table.py` needs inputs that are not shipped.
12. **Checks**: `report/scripts/check_cr_tables.py`, `check_no_mesh_narration.py`,
    `check_figure_currency.py`.
13. **LaTeX**: `report/src/argus_rans_validation_report.tex` builds `argus_rans_validation_report.pdf`.

### 8.2 The one path fix-up: a separate build directory

The scripts in `report/scripts/` find `data/` themselves: `report/scripts/report_paths.py` looks
for a `data/rans_forces.json` in the script's folder and in every folder above it, and stops with
`HALT: cannot locate the report data directory` if there is none. So the scripts that read only
`data/` run straight from your clone, for example:

```bash
cd "$REPO/report"
python3 scripts/make_cr_table.py
python3 scripts/reconcile_cr_derived.py
```

The first prints the six Condition CR table rows (the W row reads
`W & SWW\_trim & 1.5205 & 187.163 & $-0.868$ & ...`); the second ends
`DRY RUN. 0 field(s) would change. Re-run with --apply.`

Everything else still expects the layout of the source repository, with `fig/`, `sections/` and
the report's `.tex` file **next to** `data/`: the scripts in `report/fig/`, the report scripts
that also read `fig/` or `sections/` (`make_rank_table.py`, `check_cr_tables.py`,
`spanwise_prose_values.py`, `spaneff_prose_values.py` and `check_generated_tables.py`), and the
LaTeX build, because the report looks for its figures in `fig/` beside
`argus_rans_validation_report.tex`. In this repository `fig/` and `sections/` sit under `report/`
and `report/src/` while `data/` is at the top level, so rebuild that layout in a new directory
outside your clone; nothing in the clone is changed:

```bash
export BUILD=$HOME/argus-report-build
mkdir "$BUILD"
cp -r "$REPO/report/src/argus_rans_validation_report.tex" "$REPO/report/src/sections" \
      "$REPO/report/fig" "$REPO/report/scripts" "$REPO/data" "$BUILD/"
cd "$BUILD"
ls
```

Expected: `argus_rans_validation_report.tex  data  fig  scripts  sections`. If `mkdir` says the directory exists,
choose a new name rather than copying into an old build.

### 8.3 Regenerate the figures that can be regenerated

Run `make_vlmerr.py` first: three of the slide scripts read the file it writes.

```bash
cd "$BUILD"
python3 fig/make_vlmerr.py
python3 fig/make_bars.py
python3 fig/make_polars.py
python3 fig/make_spaneff.py
```

Expected: each exits with status 0 and prints `WROTE` or `wrote` lines naming its files in
`fig/`. Between them they rewrite 7 of the 23 figures the report includes:

| Script | Report figures it writes | Also writes |
|---|---|---|
| `make_vlmerr.py` | `fig_vlmerr_ranking.pdf` | `fig_vlmerr_correlation.pdf`, `fig_vlmerr_counts.pdf`, `fig_vlmerr_values.json` |
| `make_bars.py` | `fig_trimmed_performance.pdf` | `fig_trimmed_performance_values.csv`, `fig_trimmed_performance_notes.txt` |
| `make_polars.py` | `cl_alpha_condition_CR.pdf`, `cl_alpha_early_cruise.pdf`, `drag_decomposition_rans_vs_vlm_condition_CR.pdf`, `trim_alpha_rans_vs_vlm_condition_CR.pdf` | `polar_cl_cd_condition_CR.pdf`, `polar_cl_cd_early_cruise.pdf` |
| `make_spaneff.py` | `fig_spaneff.pdf` | `fig_spaneff_values.csv`, `make_spaneff_LOG.txt` |

`make_vlmerr.py` reads `data/rans_forces.json`, `data/vlm.json` and `data/dso_vlm_claims.json`
(the design study's own statements about its VLM results, used for the noise-floor band).
`make_spaneff.py` prints one line per geometry, for example
`M  CDi  62.035 ->  55.696 ct over 0.812 m (decay 10.22%),  e_last 0.9724,  CL +0.001 ct of target`.

**Do this before you build the report.** Two of the shipped figure files in `report/fig/`,
`fig_trimmed_performance.pdf` and `cl_alpha_condition_CR.pdf`, still show M's earlier Condition
CR result (for example C_D 191.25 ct and alpha 1.252 deg, on the published mesh); the copies of
these two figures inside the shipped `report/argus_rans_validation_report.pdf` do not. Regenerating them here replaces
them with M's current values
(188.26 ct, 1.257 deg). The other five report figures these scripts write came out with the
same text as the shipped files. Three further files they write, which the report does not use
(`fig_vlmerr_counts.pdf`, `fig_vlmerr_correlation.pdf`, `polar_cl_cd_condition_CR.pdf`), also
change slightly on regeneration, because their shipped copies carry M's earlier result (for
example `fig_vlmerr_counts.pdf` prints `RANS signal span -0.87 to 0.24 ct`).

The figures are drawn from the copy of `data/` you made, so they show whatever that copy holds.
None of them prints a mesh generation on the figure. **`fig_trimmed_performance.pdf` (from
`make_bars.py`) mixes mesh generations without labelling them**: its Condition CR panels are
wake-refined, and its early- and late-cruise panels are welded. Its values file,
`fig_trimmed_performance_values.csv`, has no `mesh_generation` column either. Never read a
difference across its panels, and never copy its numbers into one table without showing both
generations (section 9.1). `fig_spaneff.pdf` shows all six geometries on the wake-refined mesh,
as the `mesh_generation` column of `fig_spaneff_values.csv` says.

`fig_vlmerr_values.json` and `make_spaneff_LOG.txt` record the full path of your build
directory. Check them before you pass regenerated files to anyone.

Five other report figure scripts cannot redraw their figures here. The last line each prints is
shown, so you can recognise it:

| Script | Last line printed | Stops because |
|---|---|---|
| `make_cp.py` | `AssertionError: ('CR', [])` | the pressure section files (`data/sections_nflag/`, `data/sections_r2/`) are not shipped, so it finds no Condition CR sections |
| `make_cf.py` | `TypeError: must be real number, not NoneType` | the same section files are missing |
| `make_hinge.py` | `HALT: cannot locate results/derived. Looked for 'results/derived/SWB_trim.json' in <path>/fig and every parent. ...` (exit status 1) | it reads `results/derived/`, which is not shipped |
| `make_rootbending.py` | `HALT: no VLM baseline, so no percentage can be formed` (exit status 1) | it reads the design-study files and `results/derived/`, neither of which is shipped; it stops at the first |
| `make_spanwise.py` | exits with status 0 (see below) | `data/spanwise/` ships M's files only |

**Do not run `make_spanwise.py`.** It finishes normally, but finds spanwise data for M alone
(`no_spanwise_data` for SWB_trim, SWC_trim, SWF_trim, SWH_trim and SWW_trim), so it overwrites
`fig/spanwise_loading_condition_CR.pdf` and `fig/spanwise_caption_values.json` in `$BUILD` with
a version that has one RANS curve instead of six. If you ran it, put the shipped files back:

```bash
cp "$REPO/report/fig/spanwise_loading_condition_CR.pdf" \
   "$REPO/report/fig/spanwise_caption_values.json" "$BUILD/fig/"
```

The remaining 16 report figures (the four delta maps, the 3D concept rendering, the lineage
staircase, the five pressure and friction section figures, the two hinge figures, the
root-bending figure and the two spanwise-loading figures) are shipped as built in `report/fig/`
and cannot be regenerated from this repository. In five of them, `fig_cp_sections_CR.pdf`,
`fig_cf_sections_upper.pdf`, `fig_cf_sections_lower.pdf`, `delta_map_CR_cp.png` and
`delta_map_CR_cf_s.png`, M's Condition CR panels come from its published-mesh case, while B, C,
F, H and W are wake-refined
([README](../README.md#which-mesh-each-result-is-on-read-this-before-using-any-number), rule 7).
The two hinge figures and the root-bending figure are wake-refined at Condition CR for all six
geometries, at iteration 6000; the early- and late-cruise hinge rows are welded
(README, same section, rule 8). Moving Condition CR onto the wake-refined trims changed no
verdict: C and M still exceed the 6.8 per cent screen in RANS, and the two methods now differ by
at most 0.40 percentage points (H). It did move F, H and W more than B, C and M. The hinge
moment here is `M_h_total_Nm` in `report/fig/hinge_moment_values.csv`: pressure plus viscous,
half wing, about the morph line, at q 708.05 Pa. It is negative, so a positive change means a
smaller magnitude. From the published to the wake-refined trims it changed by +0.0159 to
+0.0180 N m for F, H and W and by +0.0036 to +0.0045 N m for B, C and M, about four times as
far; the trim angle changed by +0.029 to +0.030 deg against +0.005 deg, about six times. F, H
and W are the three meshes built from the unwelded surfaces (README, Known gaps, item 8).

Slide figures: 8 of the 18 `deck_*.py` scripts run and write PNGs to `fig/deck/`
(`deck_dcd_cr`, `deck_dcd_cruise`, `deck_decomp`, `deck_magnitude`, `deck_pipeline`,
`deck_pipeline_simple`, `deck_ranking`, `deck_trim_alpha`), for example:

```bash
python3 fig/deck_dcd_cr.py
```

Expected: `wrote fig/deck/deck_dcd_cr.png  1340x1030 px  (6.70 x 5.15 in)`. The other ten stop
with exit status 1. Nine of them stop on an input that is not shipped (a `FileNotFoundError`,
or a `HALT: no binned map ...` line), among them `deck_spanwise.py`, with a `FileNotFoundError`
naming `data/spanwise/SWB_trim_spanwise_dense.csv` (or `KeyError: 'B'` after
`make_spanwise.py`). `deck_spaneff.py` stops with `KeyError: 'subset'`, because
`data/span_efficiency.json` no longer has the entry it looks for. `deck_common.py` and
`deck_maps_common.py` are shared modules, not figure scripts.

### 8.4 The report PDF

First regenerate the figures (section 8.3). Then rewrite the span-efficiency paragraph from
`data/span_efficiency.json`, and build:

```bash
cd "$BUILD"
python3 scripts/spaneff_prose_values.py
latexmk -pdf -interaction=nonstopmode -halt-on-error argus_rans_validation_report.tex
grep "Output written" argus_rans_validation_report.log
grep -c "undefined" argus_rans_validation_report.log
```

**What you should see:** `spaneff_prose_values.py` prints `wrote
sections/12_spanefficiency_result.tex (3 paragraphs)` and the three paragraphs. `latexmk` exits
with status 0 and ends `Latexmk: All targets (argus_rans_validation_report.pdf) are up-to-date`. The first `grep`
prints `Output written on argus_rans_validation_report.pdf (36 pages, ...)`, and the second prints `0` (no undefined
references). Built this way, the text of `argus_rans_validation_report.pdf` was identical to the text of the shipped
`report/argus_rans_validation_report.pdf` (both converted to text and compared; re-checked for the release of
9 October 2026).

Why the paragraph step: `spaneff_prose_values.py` rewrites
`sections/12_spanefficiency_result.tex` from `data/span_efficiency.json`, so the paragraph cannot
drift from the data. The shipped source already carries M's current 10.22 %, and the rewrite
reproduced it byte for byte (checked), so this step is a check rather than a correction. It
writes in `$BUILD/sections/`, not in your clone. Built without section 8.3, the PDF still has 36
pages and no undefined references, but it carries M's earlier values in the two figures of
section 8.3.

**If it does not:** `! LaTeX Error: File ... not found` names a file that did not reach `$BUILD`:
redo section 8.2 in a new directory. A missing package or command means your LaTeX installation
is incomplete: see [01_setup.md](./01_setup.md), section A9.


### 8.5 The report checks

Run in `$BUILD`:

```bash
python3 scripts/check_cr_tables.py
python3 scripts/check_no_mesh_narration.py
python3 scripts/make_cr_table.py
python3 scripts/make_rank_table.py
python3 scripts/reconcile_cr_derived.py
```

Expected:

1. `check_cr_tables.py` compares every Condition CR row in the report's tables with
   `data/rans_forces.json`. It ends `rows compared: 6; disagreeing: 0; unreadable: 0`. Any
   disagreement means the report text and the data have drifted apart. It deliberately never
   compares the cruise tables with Condition CR data. A clean result shows only that the
   report and `data/` agree with each other, not that either is current: it does not look at
   figures or at the span-efficiency paragraph (sections 8.3 and 8.4).
2. `check_no_mesh_narration.py` ends `scanned 11 files, 20 forbidden terms, 0 hit(s), 6 pinned
   exception(s)`.
3. `make_cr_table.py` prints two comment lines first, the rank order
   (`% rank order by dCd, cheapest first: W, F, C, H, M`) and
   `% 6 of 6 rows within 1.0 count of target CL`, then the Condition CR table rows as LaTeX,
   one per geometry. The M row reads `M & SWM\_trim & 1.2574 & 188.263 & $+0.232$ & ...`.
4. `make_rank_table.py` (run `fig/make_vlmerr.py` first) prints the selection-test table and
   the rank statistics the text quotes. Its last line is
   `% RANS increment range at Condition CR: 1.100 ct (was 6.729 before the re-trim)`. Among its
   rows are the trim-angle rows of the report's Table 16 (section 5.3 of the report, "Trim
   angle"), in the block headed `%% table body for tab:qv:lift`, including the current M row,
   which reads `M & 1.257423 & 1.24972 & $-0.0077$ & ...`. The shipped `report/argus_rans_validation_report.pdf`
   and `sections/t5_results.tex` carry the same row.
5. `reconcile_cr_derived.py` ends `DRY RUN. 0 field(s) would change. Re-run with --apply.` If it
   reports fields that would change, the data file contradicts itself; report it rather than
   applying.

`report/scripts/spanwise_prose_values.py` prints (and does not write) the numbers the spanwise
paragraphs quote, from the shipped `fig/spanwise_caption_values.json`. It runs in `$BUILD` and
exits with status 0, but only while that file is the shipped one (not after `make_spanwise.py`,
section 8.3).

`check_figure_currency.py` compares file dates and says which figures are older than the data.
In a fresh copy every file has the date of the copy, so it reports every figure as current and
tells you nothing. Use it only in a directory where the dates are real.

## 9. Quoting results

### 9.1 The rules

1. **Never difference, fit or rank across operating points.** Condition CR, early cruise and
   late cruise each have their own baseline trim. A drag change is always a candidate minus B
   **at the same condition**.
2. **State the frame with every coefficient**: the reference basis (DSO: S_ref 1.24092 m²,
   C_ref 0.39396 m; half-model A_ref 0.620462 m²), the condition (Mach and target C_L), and the
   mesh generation. Every RANS coefficient delivered here is on the DSO basis. The early- and
   late-cruise VLM results in `data/vlm.json` keep the TP-1580 basis of their input decks, and
   whether the cruise target C_L values are on the DSO basis is unresolved
   ([README](../README.md), "Reference quantities and units", items 6 and 7). The project also
   reports a second basis alongside, the TP-1580 wind-tunnel basis (S_ref 1.1148 m²,
   c_ref 0.3404 m, the same b_ref 3.6576 m). The two areas differ by a factor of 1.1131
   (11.31 %), which on a C_D of 188 counts is about 21 counts, against a morphing effect of
   1.100 counts (the Condition CR ΔC_D range, section 9.2). Never compare a coefficient on one
   basis with one on the other without converting it.
3. **Mesh generations.** Condition CR trims are on the wake-refined mesh (label
   `wake_refined` in `data/rans_forces.json`, spelled `wake-refined` with a hyphen in
   `data/span_efficiency.json`; search for both). Early- and late-cruise trims, and the off-trim
   legs that reached them, are welded (label `welded`). The published cruise alpha-bracket cases
   and the Condition CR angle-of-attack sweep cases that bracket each trim (names like
   `SWB_a1p75`) are on the earlier, published mesh (label `published`); no cruise result uses the
   former. Any Condition CR trim value kept from the published mesh is superseded (section 5.5),
   and so is every published cruise trim value
   ([README](../README.md#which-mesh-each-result-is-on-read-this-before-using-any-number),
   rule 9). **Never put a Condition CR trim and a published number in one table without
   showing both generations**,
   and never difference across them. The report text itself does not say which mesh each number
   is on (its own check, `check_no_mesh_narration.py`, keeps mesh-history terms out of it), so
   take the generation from this section and the [README](../README.md), not from the report.
4. **The published Condition CR span efficiencies W 1.2861, H 1.2974 and F 1.3028 are the
   unresolved-wake artefact** of the earlier mesh, wherever they appear. A span efficiency above
   1 is impossible for a planar wake. Do not quote them.
5. **Compare span efficiency only at the same plane and on the same mesh generation**
   (section 6.3). Its absolute value depends on the plane it was taken at.
6. **Quote window means from `.converged`**, not a single iteration (section 5.3).
7. **Quote a run card only after it validates** (section 7.5).

### 9.2 Current Condition CR results

Condition CR: M 0.10, U 34.0 m/s, every geometry at its own trimmed alpha, target
C_L 0.428277635108; all six within 1.0 count of target C_L. C_D in counts (1 count = 0.0001) on
the DSO basis, half-model A_ref 0.620462 m²; ΔC_D is candidate minus B at the same condition,
negative is a drag saving. e is the span efficiency at the last of the four Trefftz planes
(x = 3.206 m), AR 10.7807: `e_last` in `data/span_efficiency.json`, rounded to four decimals.
`data/rans_forces.json` carries a separate `span_efficiency_x4`, for B, F, H and W only, that
differs in the fourth decimal for F (0.9858) and H (0.9771). Mesh generation: **wake_refined**
for all six.

| Letter | `candidate` | C_D (ct) | ΔC_D (ct) | e, last plane |
|---|---|---|---|---|
| B | `baseline` | 188.031 | 0 (reference) | 0.9673 |
| W | `cffw_b01_c01` | 187.163 | -0.868 | 0.9837 |
| F | `cfft_b02_c01` | 187.216 | -0.815 | 0.9859 |
| H | `chc_g02_c06` | 188.022 | -0.009 | 0.9773 |
| C | `cte_i002_c04` | 187.840 | -0.191 | 0.9773 |
| M | `mcv2_i002_c01` | 188.263 | +0.232 | 0.9724 |

The ΔC_D range is -0.868 to +0.232 ct (1.100 ct), and the order by C_D is W < F < C < H < M.
Every e in the table is at the same plane (x = 3.206 m) on the same mesh generation, so the
column can be compared down the table: what the result says is each candidate's e relative to
B's. The absolute values are not a property of the wing alone, because C_Di, and with it e,
still changes by 9.46 to 10.22 % between the first and the last plane (section 6.3). An e
taken at another plane, on another mesh or by another method is not comparable with these.

For B, C_D 188.031 ct is the window mean `Cd 0.0188031` in its `.converged` marker; the run
card's `results.CD_counts` (188.030265) is not the window mean and differs in the third
decimal.

The `trims` rows of `data/rans_forces.json` (section 5.5), `data/span_efficiency.json` and the
shipped `report/argus_rans_validation_report.pdf` carry exactly these values. Five shipped files in
`report/fig/` still carry M's earlier result until you rebuild as in section 8.3: the report
figures `fig_trimmed_performance.pdf` and `cl_alpha_condition_CR.pdf`, and three files the report
does not use (`fig_vlmerr_counts.pdf`, `fig_vlmerr_correlation.pdf`,
`polar_cl_cd_condition_CR.pdf`).

Two places in the report are not wholly wake-refined at Condition CR even after that rebuild:
M's panels in the section and delta-map figures (section 8.3), which are on the published mesh;
and the lift-curve slopes (Table 16's dC_L/dalpha and implied-dC_L columns, the section 5.3
"Slope" sentence, the slope in `cl_alpha_condition_CR.pdf`), which are fitted over published
sweep cases plus the wake-refined trim (README rule 2). Do not quote either as wake-refined. The hinge-moment and root-bending results at Condition CR are wake-refined for
all six geometries
([README](../README.md#which-mesh-each-result-is-on-read-this-before-using-any-number),
rule 8).

## 10. Script reference

Every script in this table exists at the path given. (Three scripts named elsewhere in this
guide are **not** shipped and are named only so you can ask for them:
`scripts/build_r2_run_cards.py`, `collect_run_facts.sh` and `scripts/wing_maps_of12.py`.
Three more, not shipped either, appear only in the `consumer:` lines of
`assert_case_postpro.py`, section 4.1.)
"Here" means run from a clone of this repository (for `report/` scripts, from the build
directory of section 8.2).

| Script | Arguments | Reads | Runs here? |
|---|---|---|---|
| `scripts/install_postpro.py` | case directories, `--dry-run` | the case's `system/controlDict`; `constant/polyMesh/boundary` (or `processor0`'s) for the wall patch, else the patch an existing forces object names; `0/U` if it exists, else `0.orig/U`; `constant/momentumTransport`; for a cruise case also `0.orig/p`, `0.orig/T` (or `0/p`, `0/T`) and `constant/physicalProperties` (or `thermophysicalProperties`); a sampled wall surface under `postProcessing/`, else the ASCII STL in `constant/triSurface`; and `hpc/templates/argusPostPro` | **Yes**, needs numpy (section 4.1); stops with `HALT` if the case has no ASCII STL |
| `scripts/assert_case_postpro.py` | `--before` or `--after`, case directories; `--self-test` | `system/controlDict` and its includes, or `postProcessing/` | **Yes** (section 4); `--after` gives three wrong `MISSING` lines on a real mesh (4.6) |
| `scripts/prune_legacy_fo.py` | `controlDicts`, `--dry-run`, `--self-test` | a `controlDict` | not needed for this guide |
| `scripts/splice_postpro.py` | one directory holding case directories | `postProcessing/forceCoeffs1/` | **Yes**, on copies (section 5.4) |
| `scripts/harvest_rans_forces.py` | none | `$ARGUS_SOLVE` or `/home/scratch/$USER/argus/solve`, and `~/argus_archive` | **No** as it stands: its default roots are those of the account that ran the campaign. It runs with `ARGUS_SOLVE` set to a solve root and an existing `~/argus_archive`, even an empty one (section 8.1 item 3) |
| `scripts/spanwise_dense_report.py` | `--vtk --case --out`, optional `--case-dir --time --cl-ref --eta-min --eta-max --n --band` | a wing-surface VTK with p and cp; the `dragDir` of the latest `forceCoeffs1` directory | **Yes**, given a wing-surface VTK (section 6.4) |
| `scripts/spanwise_dense.py` | `--vtk --tag --alpha`, optional `--uinf` (default 40.8, the tunnel speed, not Condition CR's 34.0), `--out` (default `results/`) | a wing-surface VTK | not used for the report; forms Cp as p / (0.5 U²), which is wrong for the compressible cruise cases |
| `scripts/section_slices_from_surface.py` | `--vtk --case --out` (all required), optional `--eta --band --semi --case-dir --surface-flag` | a wing-surface VTK; the `dragDir` of the latest `forceCoeffs1` directory, as in section 6.4 | needs a wing-surface VTK; its outputs feed `make_cp.py` and `make_cf.py`, whose other inputs are not shipped |
| `scripts/wing_surface_maps.py` | `--dir --tag --alpha` (all required), optional `--case --time --uinf` (default 40.8, the tunnel speed, not Condition CR's 34.0), `--n-eta --n-xc --out` (default `results/wing_maps`) | a wing-surface time directory | needs a wing-surface file, which only a case with the block writes |
| `scripts/render_delta_maps_of12.py` | `--quantity --condition`, or `--all`; `--null`, `--outdir` | binned maps in `results/wing_maps_of12/`, made by `scripts/wing_maps_of12.py` | **No**: neither the maps nor `wing_maps_of12.py` is shipped. With `--all` it refuses every figure (`PARTITION: 6 requested = 0 drawn + 6 refused`). The delta-map figures in the report are shipped as built |
| `scripts/render_concept_3d.py` | `--outdir --size --no-caption` | `geometry/stl/<name>.stl.gz` | **No**: stops with `FileNotFoundError` on `geometry/stl/baseline.stl.gz`, because no geometry is shipped. Used by `report/fig/deck_concept.py` |
| `scripts/hinge_moment_integrate.py` | `--vtk --tau-vtk --no-viscous --alpha --q --semispan --out`, `--self-test` | a wing-surface VTK | `--self-test` runs and ends `self-test: PASS, 9/9 ...`; a real run needs a wing-surface VTK |
| `report/fig/section_source.py` | none: a shared module, not a script | `data/sections_r2/`, `data/sections_nflag/` | imported by `make_cp.py`, `make_cf.py`, `deck_common.py` and `deck_sections_eta080.py`; its inputs are not shipped |
| `scripts/case_derived_quantities.py` | case names or `--all`, `--quick` | `results/postpro_latest/`, `cases/of12/v5/solve/` | **No** (section 6.3); imported by `build_span_efficiency.py` |
| `scripts/trefftz_trims.py` | `--json` | `docs/report/all_geometry_2026-09-15/data/rans_forces.json` | **No**: that path does not exist here |
| `scripts/trefftz.py` | `--case --alpha`, optional `--uinf` (default 40.8), `--x --reconstructed --te --chord`; run with `pvbatch` | the case's volume field | **No**: not the route of the delivered numbers, and under ParaView 5.10.1 `pvbatch` it exits 0 having printed nothing, because its final `os._exit` discards the buffered output |
| `scripts/derive_run_card_fields.py` | `--case <dir>` or `--audit-all` | controlDict, solver log, momentumTransport, `0/nut`, `0/k`, `0/omega`, fvSolution | **Yes** with `--case` (section 7.2); `--audit-all` finds nothing here |
| `scripts/validate_run_card.py` | card files, `--self-test`, `--no-registry-check` | `scripts/run_card.schema.json`, `registry/candidates.yaml` | **Yes**, with jsonschema 4 and yaml (section 7.5); not on the HPC12 login node |
| `scripts/watch_run.py` | cases, `--once --plot --interval --ref-cd --ref-cl` | an older case layout | reported 0 iterations on every case this guide was checked with, including a delivered Condition CR trim; not used here |
| `report/scripts/harvest_case_provenance.py` | `[CASE_ROOT] [-o OUT.json]` | case directories | **Yes** with both given |
| `report/scripts/build_rans_forces.py` | raw harvest, provenance | the two JSON files | **do not use**: fails in the clone, and in the build directory it overwrites `data/rans_forces.json` with a regressed file (section 8.1 item 5) |
| `report/scripts/update_cr_from_r2.py` | none | values written into the script | **do not run** (section 8.1 item 6) |
| `report/scripts/reconcile_cr_derived.py` | `--apply` | `data/rans_forces.json` | **Yes**, from `$REPO/report` or the build directory (section 8.2) |
| `report/scripts/build_span_efficiency.py` | `--raw` (required), `--tag`, `--out`, `--allow-subset --reason` | `data/span_eff_r2trim_raw.json` | **Yes**, with `PYTHONPATH="$REPO/scripts"` and `--out` (section 6.3) |
| `report/scripts/extract_vlm.py` | none | `$ARGUS_ROOT/dso_reference/` (else the clone's `dso_reference/`) | **No**: the design study's files are not shipped (section 8.1 item 9) |
| `report/scripts/make_cr_table.py` | none | `data/rans_forces.json` | **Yes**, from `$REPO/report` or the build directory (section 8.2) |
| `report/scripts/report_paths.py` | (a module, not run directly) | finds the folder whose `data/rans_forces.json` exists, starting from the script's own folder | imported by the report scripts (section 8.2) |
| `report/scripts/make_rank_table.py` | none | `data/rans_forces.json`, `data/vlm.json`, `data/span_efficiency.json`, `fig/fig_vlmerr_values.json` | **Yes**, in the build directory, after `fig/make_vlmerr.py` |
| `report/scripts/check_generated_tables.py` | none | the output of the table generators and `sections/*.tex` | **Yes**, in the build directory; it checks that every row the generators print is in the report source, and ends `VERDICT: OK, every generated row is in the document` (12 rows checked, 0 stale) |
| `report/scripts/make_q_table.py` | `--sections` (required), `--time CASE=TIME`, `--check` | a sections directory with logs | **No**: sections not shipped |
| `report/scripts/spaneff_prose_values.py` | none | `data/span_efficiency.json` | **Yes**, in the build directory; writes `sections/12_spanefficiency_result.tex` (section 8.4) |
| `report/scripts/spanwise_prose_values.py` | none | `fig/spanwise_caption_values.json` | **Yes**, in the build directory, with the shipped file; prints only (section 8.5) |
| `report/scripts/check_cr_tables.py`, `check_no_mesh_narration.py` | none | `sections/*.tex`, `argus_rans_validation_report.tex`, `data/` | **Yes**, in the build directory |
| `report/scripts/check_figure_currency.py` | `--since`, `--fail-stale` | file dates | runs, but see section 8.5 |
| `report/fig/make_vlmerr.py`, `make_bars.py`, `make_polars.py`, `make_spaneff.py` | none | `data/` | **Yes**, in the build directory |
| `report/fig/make_spanwise.py` | none | `data/spanwise/`, `data/vlm.json`, `data/rans_forces.json` | runs, but draws M alone: **do not run it** (section 8.3) |
| `report/fig/make_cp.py`, `make_cf.py`, `make_hinge.py`, `make_rootbending.py` | none | unshipped inputs | **No** (section 8.3) |
| `report/fig/deck_*.py` | none | `data/`, `fig/` | 8 of 18 (section 8.3) |

## 11. Post-processing pieces that are not in this repository

1. The run-card builder and the cluster-side script that collects its facts
   (`scripts/build_r2_run_cards.py` and `collect_run_facts.sh`, as the shipped tooling names
   them). Without them you cannot assemble a complete card for a new run (section 7.1).
2. The script that pulled `postProcessing/` from HPC12 into the layout
   `scripts/case_derived_quantities.py` reads (`results/postpro_latest/`), and that layout
   itself. Without them the planes of a new case cannot be integrated by a shipped route
   (section 6.3). The integrals of the six delivered trims are shipped
   (`data/span_eff_r2trim_raw.json`).
3. The step that added the `mesh_generation` labels to `data/rans_forces.json`.
4. The inputs of the non-regenerable figures: `data/sections_nflag/`, `data/sections_r2/`,
   the spanwise files of B, C, F, H and W (`data/spanwise/` holds M's only), `results/derived/`,
   `results/wing_maps_of12/` (and `scripts/wing_maps_of12.py`, which made it), `geometry/stl/`,
   and the design study's files (the `dso_reference/` folder of section 8.1 item 9).
5. Raw fields, logs and `postProcessing/` of the delivered cases. What was read from them is in
   `data/`.

Ask the project for item 1 before you need a complete run card for a new run. This repository
gives no contact address; ask whoever gave you access to it.
