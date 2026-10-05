# 03. Build a case: baseline B at Condition CR, on your own machine

> **At a glance.** This guide builds and checks a tutorial-resolution mesh of the baseline
> wing B at Condition CR on your own computer, then (optionally) checks that the solver starts.
> **You need:** [01, setup](./01_setup.md) done, and the baseline wing surface that you build
> from B's `.vsp3` file with [02](./02_geometry_from_openvsp.md) (see the box below). **You
> end with:** a mesh of about 1.19 million cells whose `checkMesh` shows only the four
> expected failures (Step 13), the same four that every delivered production mesh shows, for
> learning the pipeline, never for results. **Time:** building the surface with 02 took about
> 30 s, `snappyHexMesh` 94 to 112 s on 8 cores, and the optional solver check 1.5 to 4
> minutes; allow about 2.2 GB of extra memory and 440 MB of disk.

This is the worked example. You take one wing (B, the baseline), one flow condition
(Condition CR), and go from an empty folder to a checked mesh on your own computer, before
anything touches the cluster. At the end there is an optional check that the solver starts.

It comes after [01, setup](./01_setup.md) and [02, geometry](./02_geometry_from_openvsp.md),
and before [04, sending a case to HPC12](./04_send_to_hpc.md). Start at the
[README](../README.md) if you have not read it.

Every command below was run for this guide under OpenFOAM-org 12 (the Ubuntu `openfoam12`
package, installed in `/opt/openfoam12`) with Open MPI 4.0.3, on a machine with 12 physical
cores and 23 GB of memory, and then run a second time exactly as printed here, from a fresh
folder. Every "What you should see" block is copied from those runs (in a terminal, `ls`
spreads its output over columns rather than lines). Where something could not be checked,
the text says `NOT VERIFIED`.

> [!WARNING]
> **TUTORIAL RESOLUTION: FOR LEARNING THE PIPELINE, NOT FOR RESULTS.**
> The mesh built here is deliberately coarse so that it fits on an ordinary computer:
> about 1.19 million cells against 99,111,506 for the delivered baseline mesh
> (`data/run_cards/SWB_r2_trim.json`). Its forces are **not comparable to anything in
> `data/`**, and must never be reported, plotted beside project results, or used to judge a
> geometry. The production mesh is built from the **unedited** recipe on the cluster
> ([04](./04_send_to_hpc.md), [05](./05_submit_and_run.md)).

> [!IMPORTANT]
> **NO GEOMETRY IS SHIPPED. YOU BUILD THE WING SURFACE YOURSELF, BEFORE THIS GUIDE.**
> The wings come as OpenVSP `.vsp3` files (see [02](./02_geometry_from_openvsp.md), Step 1,
> for where each one is). For this guide you need one of them, the baseline B, `baseline_wing_only_refined.vsp3`, with sha256
> `55cbdc7c48a15d847da602b540b95712b38e893490823aa7b5b997302100ba02` (the `delivered_vsp3`
> entry for `baseline` in `registry/candidates.yaml`). Build the meshing surface from it with
> [02](./02_geometry_from_openvsp.md):
>
> 1. Section 4 (set-up), then Section 5, Step 1: download the `.vsp3` and check it against
>    that checksum.
> 2. Section 5, Steps 2 to 7: build the surface. Section 11, variant A, gives the same build
>    as one block of commands.
> 3. Section 5, Steps 8 and 9: check it with `surfaceCheck` and compare its checksum with the
>    registry.
>
> The result is `$HOME/argus-geometry/B/baseline_oml_placed.stl` if you followed Section 5,
> or `$HOME/argus-geometry/B/baseline_wing_only_refined_oml_placed.stl` if you used
> Section 11. Either way it is 52,623,882 bytes, and its sha256 must be
> `20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641`, the
> checksum of the registered baseline surface (`baseline_oml_placed_laddercap.stl` in
> `registry/candidates.yaml`). Rebuilt that way for this guide, it matched exactly. Nothing
> below works without it.

---

## What you will build

| Item | Value | Where it comes from |
|---|---|---|
| Geometry | B, registry candidate `baseline` | `candidate` field of `data/run_cards/SWB_r2_trim.json`; `registry/candidates.yaml` |
| Surface | the file you build with [02](./02_geometry_from_openvsp.md), `baseline_oml_placed.stl` (the registry calls the same surface `baseline_oml_placed_laddercap.stl`): half wing, in metres, root extended 20 mm through the symmetry plane and capped | `registry/candidates.yaml`, entry `baseline`, `geometry_3d_laddercap` |
| Condition | CR: Mach 0.10, U 34.0 m/s, kinematic viscosity 1.46e-05 m2/s | `recipes/README.md`, `0.orig/U`, `constant/physicalProperties` |
| Reynolds number | about 0.79 million on the TP-1580 chord (0.3404 m); about 0.92 million on the DSO reference chord (0.39396 m). Same flow, two chord bases. | U x chord / nu |
| Recipe | `recipes/L11_wallResolved_CR`: kOmegaSST, steady, wall-resolved (the first cell next to the wing sits at y+ about 1, with low-Reynolds-number wall boundary conditions; see Step 3) | `recipes/README.md` |
| Angle of attack | 1.30 deg, the recipe's **starting value**, not a trim | `0.orig/U` and `system/controlDict` |

This guide does not trim the wing to its target lift coefficient (C_L 0.428277635108 at
Condition CR). That is [06](./06_trim_to_target_cl.md). It also uses only Condition CR. The
project has three operating points (Condition CR, early cruise, late cruise), each with its
own baseline trim, and **results from different operating points are never differenced
against each other**.

For OpenFOAM-12 the recipe folder **is** the case skeleton. You copy the recipe.

---

## Before you start

1. [01_setup.md](./01_setup.md) is done: OpenFOAM-12 loads in your terminal (section A2),
   `mpirun` is Open MPI (section A5), and ParaView is installed (section A7).
2. You have built the baseline surface with [02](./02_geometry_from_openvsp.md), and its
   sha256 is the one in the box above.
3. What this costs, measured on the machine described above:

   | Step | Time | Peak memory |
   |---|---|---|
   | Building the surface ([02](./02_geometry_from_openvsp.md), Section 11, variant A) | about 30 s | not measured |
   | `snappyHexMesh` on 8 cores | 94 to 112 s (nine runs) | about 2.2 GB more than the machine was already using |
   | `snappyHexMesh` on 4 cores | 127 s and 302 s (two runs) | about 1.4 GB more |
   | `checkMesh` | 28 s | small |
   | Optional solver check (Step 15) | 1.5 to 4 minutes (six runs) | small |

   Times depend on what else the machine is doing; the ranges above are what was seen.

   The finished case takes 440 MB of disk (450 MB with the optional Step 15a). Keep it
   **outside** your clone of this repository.

---

## Step 1. Open a terminal with OpenFOAM-12 and Open MPI

Every command in this guide runs in this terminal (the two exceptions, optional Python steps
in Step 6 and Step 15a, say so). If you open a new terminal later, repeat this step (including
the fix in item 2 below, if you needed it), the `export` lines from Steps 2, 3 and 5, and
then `cd "$CASE"`: from Step 4 on, every command uses paths relative to the case folder.

```bash
source /opt/openfoam12/etc/bashrc
foamVersion
which foamRun
mpirun --version | head -1
export REPO=$HOME/argus-production   # or wherever you cloned it
ls $REPO/recipes
```

**What you should see:**

```
OpenFOAM-12
/opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/foamRun
mpirun (Open MPI) 4.0.3
L11_wallResolved_CR  L11_wallResolved_WT  M6_wallModelled_CRUISE_EARLY  M6_wallModelled_CRUISE_LATE  README.md
```

**If it does not:**

1. `foamVersion` prints something other than `OpenFOAM-12` (for example `OpenFOAM-7`, or
   `bash: foamVersion: command not found`), or `which foamRun` prints nothing, or `foamRun`
   gives `bash: foamRun: command not found`. Look at what the `source` line printed:
   1. `bash: /opt/openfoam12/etc/bashrc: No such file or directory`: OpenFOAM-12 is not
      installed there. Do [01_setup.md](./01_setup.md), section A1, first.
   2. No such message: the `source` line was most likely not run in this terminal, for
      example because you opened a new one. On the project's machine the login shell
      loaded OpenFOAM-org 7, which has no `foamRun`. Run the `source` line again (loading
      OpenFOAM-12 on top of version 7 works: tested), then fix your `~/.bashrc` as in
      [01_setup.md](./01_setup.md), sections A2 and A4.
2. The `mpirun` line prints `HYDRA build details:` instead of `mpirun (Open MPI) 4.0.3`:
   another MPI (usually Anaconda's) is ahead of the system Open MPI on your `PATH`. If you
   carry on, every parallel step below fails: every rank prints
   `attempt to run parallel on 1 processor` (reproduced for this guide). In the tests for
   this guide the command then returned within a few seconds, but
   [01_setup.md](./01_setup.md), section A5, records a run that hung instead: if a parallel
   command has not returned when it should have, press Ctrl+C and look in its log. Fix it
   as in [01_setup.md](./01_setup.md), section A5:

   ```bash
   export PATH=/usr/bin:$PATH
   source /opt/openfoam12/etc/bashrc
   mpirun --version | head -1
   ```

3. Sourcing prints `x86_64-conda-linux-gnu-cc: error: unrecognized command-line option
   '--showme:link'`: same cause as item 2, same fix.

---

## Step 2. Decide how many cores to use

```bash
lscpu | grep -E "^Core\(s\) per socket|^Socket\(s\)"
```

**What you should see** (the numbers are your machine's; these are the test machine's):

```
Core(s) per socket:                      12
Socket(s):                               1
```

Your number of physical cores is the product of the two. Pick a number **no larger** than
that, and leave a couple free if you want to use the computer meanwhile. This guide uses 8:

```bash
export NP=8
```

Do not use the number that `nproc` prints. On the test machine `nproc` printed 24 (it
counts hyper-threads), and Open MPI then refused a 16-rank run with
`There are not enough slots available in the system to satisfy the 16 slots that were
requested by the application`.

---

## Step 3. Copy the recipe into a new case folder

```bash
export CASE=$HOME/argus_cases/B_CR_tutorial
mkdir -p $HOME/argus_cases
mkdir "$CASE"
cp -r "$REPO/recipes/L11_wallResolved_CR/." "$CASE/"
cd "$CASE"
ls -A
```

**What you should see** (the order may differ):

```
0.orig  MANIFEST.sha256  RECIPE_FILES.json  constant  system
```

**If it does not:** if `mkdir` says `File exists`, that folder is already in use. Pick
another name rather than copying on top of old files.

What you copied. You do not need to understand every file to finish this guide, but it helps
to know where things are set:

| File | What it sets |
|---|---|
| `0.orig/U`, `k`, `omega`, `nut`, `p` | Freestream and wall values. `U` holds the inflow vector `(33.991249 0 0.771369)`: 34.0 m/s at 1.30 deg. On the wing, `nut` uses `nutLowReWallFunction` (it sets the turbulent viscosity to 0 at the wall), `k` uses `kLowReWallFunction`, and `omega` uses `omegaWallFunction` (in OpenFOAM-12 this takes its viscous-sublayer value wherever y+ is below about 11.5). See "Wall-resolved" below. |
| `constant/momentumTransport` | The turbulence model, `kOmegaSST` |
| `constant/physicalProperties` | Kinematic viscosity, `1.46e-05` |
| `system/blockMeshDict` | The background box: 60 x 30 x 60 m of 1 m cubes, and the far-field patches `symmetry`, `inlet`, `outlet`, `outboard`, `topBottom` |
| `system/snappyHexMeshDict` | The mesh recipe. It pulls in three more files: |
| `system/meshLevels` | ...the refinement levels |
| `system/layerSettings` | ...the prism layers on the wing: 11 layers, first layer 4.5 micrometres, growth ratio 1.45 |
| `system/meshQualityDict` | ...the quality limits snappyHexMesh must respect |
| `system/surfaceFeaturesDict` | Which sharp edges of `wing.stl` to extract |
| `system/decomposeParDict` | How many pieces the mesh is split into for the parallel mesh build: 192, method `scotch`, as all six delivered Condition CR meshes were built. The solve on HPC12 is split differently (128 pieces): `hpc/stage_solve_case.sh` writes the solve case's own copy of this file ([05](./05_submit_and_run.md)) |
| `system/controlDict` | `application foamRun`, `solver incompressibleFluid`, a first run of 4000 iterations (`endTime 4000`), fields saved every 500 iterations (`writeInterval 500`), and the force output (`forceCoeffs1`) |
| `system/fvSchemes`, `system/fvSolution` | Numerical schemes and linear solvers |
| `system/planeSample` | Cutting planes for post-processing ([07](./07_postprocessing_and_report.md)) |
| `MANIFEST.sha256`, `RECIPE_FILES.json` | Checksums and provenance of the recipe files |

One file is not in the recipe: `system/argusPostPro`, the post-processing block, which makes
the solver write what the report needs. A production case gets it before it goes to HPC12
([04](./04_send_to_hpc.md), Step 5e). This tutorial does not need it; the optional Step 15a
shows what it does.

**Wall-resolved.** y+ is the height of the first cell above the wall, measured in the
wall's own viscous length scale. About 1 means the mesh resolves the thin viscous layer
next to the wall instead of bridging it with a formula (a wall function). Being
wall-resolved needs two things together: the mesh (the production layer settings put the
first cell at y+ about 1; the delivered baseline solve reports a mean y+ of 1.24 in
`data/run_cards/SWB_r2_trim.json`) and wall boundary conditions that suit such a mesh (the
three above). Neither alone is enough. The tutorial mesh keeps these boundary conditions
but its first cell is 8 times taller (Step 7), so it is **not** wall-resolved.

Four things in the comments will confuse you if nobody warns you. Trust the values, not
the comments:

1. Several comments say "OpenFOAM-org 7" or "OpenFOAM-7". That is history: the recipe was
   first written for version 7. The files themselves are OpenFOAM-12 files
   (`application foamRun`, `constant/momentumTransport`, `constant/physicalProperties`).
2. The comment at the top of `0.orig/U` talks about 40.8 m/s and 2.1028 deg, and the comment
   in `0.orig/k` also works with 40.8 m/s. Those comments are left over from the
   wind-tunnel condition. The value that counts is the `Uinf` line,
   `(33.991249 0 0.771369)`, which is 34.0 m/s at 1.30 deg, and it agrees with `liftDir`,
   `dragDir` and `magUInf 34.0` in `system/controlDict`.
3. The size comments in `system/meshLevels` assume an older 0.25 m background cell, and the
   header of `system/layerSettings` describes an older 36-layer stack at 40.8 m/s. The live
   values are what count: the background cell is 1 m, so level n is 1/2^n m (level 11 is
   0.488 mm, level 8 is 3.906 mm), and the layer entries are 11 layers, first layer
   4.5e-06 m, ratio 1.45 (a 0.586 mm stack).
4. The comment above the `potentialFlow` block in `system/fvSolution` says that without the
   block `potentialFoam` aborts on every rank. Under OpenFOAM-12 it does not (tested): it
   runs to the end and takes the `SIMPLE` block's `nNonOrthogonalCorrectors`, which is 1 in
   this recipe, instead of the block's 10. Keep the block; Step 15b shows how to check it.

---

## Step 4. Check that the recipe arrived intact

```bash
sha256sum -c MANIFEST.sha256
```

**What you should see:** 19 lines, every one ending in `OK`:

```
0.orig/U: OK
0.orig/k: OK
...
system/snappyHexMeshDict: OK
system/surfaceFeaturesDict: OK
```

**If it does not:** a `FAILED` line means that file in your clone differs from the recipe
as shipped. Do not carry on. Delete the case folder, restore the clone, and start Step 3
again. To see what changed in your clone:

```bash
git -C $REPO status
```

---

## Step 5. Put the wing surface in place

> [!IMPORTANT]
> **The surface is not shipped** (see the box at the top). You built it from B's `.vsp3`
> with [02](./02_geometry_from_openvsp.md). It is
> `$HOME/argus-geometry/B/baseline_oml_placed.stl` if you followed 02, Section 5, Steps 4
> to 7, or `$HOME/argus-geometry/B/baseline_wing_only_refined_oml_placed.stl` if you used
> 02, Section 11, unless you chose another folder for `GEOM` in 02, Section 4. Set `STL` to
> your file. The checksum below decides whether it is the registered surface, not the name.

```bash
export STL=$HOME/argus-geometry/B/baseline_oml_placed.stl
mkdir -p constant/triSurface
cp "$STL" constant/triSurface/wing.stl
sha256sum constant/triSurface/wing.stl
grep -o "baseline_oml_placed_laddercap.stl, sha256: [0-9a-f]*" "$REPO/registry/candidates.yaml"
grep "^solid" constant/triSurface/wing.stl
```

**What you should see:** the same checksum twice, then one solid:

```
20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641  constant/triSurface/wing.stl
baseline_oml_placed_laddercap.stl, sha256: 20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641
solid wing_closed
```

**If it does not:**

1. `cp: cannot stat '...': No such file or directory`: `STL` does not point at your
   surface. List what you have, then set `STL` again:

   ```bash
   ls $HOME/argus-geometry/B/*.stl
   ```

   Use the file whose name ends in `_oml_placed.stl`. `baseline_oml.stl` is an unplaced
   intermediate from 02, Step 4, and must not be meshed.
2. The checksums differ: you have a different file. It may still be a usable wing, but it
   is not the registered baseline. Check the `.vsp3` (02, Section 5, Step 1) and the build
   commands (02, Section 11) and build it again. For this tutorial you may carry on to learn
   the steps; the surface check in Step 6 then matters even more.
3. More than one `solid` line: the file holds more than one surface (for example parts that
   are not the wing). The 02 route writes only the wing, so this comes from some other
   source, such as an OpenVSP STL export. Remove the other parts by name as described in
   [02](./02_geometry_from_openvsp.md), Section 7, before carrying on.

Why `constant/triSurface/wing.stl` and nothing else:

1. The recipe refers to the surface as `wing.stl`. The name matters; the folder it came from
   does not.
2. The cluster mesh job (`hpc/mesh_hpc12.pbs`) checks only `constant/triSurface/*.stl`.
3. **Never create `constant/geometry/`.** OpenFOAM-12 chooses the **folder**, not the file:
   if a `constant/geometry/` folder exists, `surfaceFeatures` and `snappyHexMesh` read
   surfaces from it and ignore `constant/triSurface/` completely. Both cases were tested.
   With an empty decoy `wing.stl` in `constant/geometry/` and the real surface in
   `constant/triSurface/`, `surfaceFeatures` read the decoy, reported
   `Triangles    : 0`, and still finished without an error. With an empty
   `constant/geometry/` folder and the real surface in `constant/triSurface/`, it stopped
   with `Cannot read ".../constant/geometry/wing.stl"`.

---

## Step 6. Check the surface with `surfaceCheck` (the input gate)

This is the most important check in the whole pipeline, and it is a check on the **input**.
A mesh built from a surface that is open, in the wrong units or the wrong wing can still pass
every mesh-quality check afterwards: `checkMesh` judges the mesh it is given, not whether it
is the domain you wanted. [02](./02_geometry_from_openvsp.md), Section 1, explains the
history, and its Step 8 shows a good and a bad result.
Run it on the copy inside the case, because that is the file that will be meshed.

```bash
surfaceCheck constant/triSurface/wing.stl > log.surfaceCheck 2>&1
grep -E "^Triangles|^Bounding Box|illegal triangles|is closed|unconnected parts|Number of zones" log.surfaceCheck
```

**What you should see:**

```
Triangles    : 192588
Bounding Box : (1.17653 -0.02 -0.175402) (2.33626 1.82541 0.0136577)
Surface has no illegal triangles.
Surface is closed. All edges connected to two faces.
Number of unconnected parts : 1
Number of zones (connected area with consistent normal) : 1
```

How to read it:

1. **Units.** The bounding box is (x, y, z) minimum, then (x, y, z) maximum. The span runs
   along y. A maximum y of **1.82541** means metres: the wing tip reaches y = 1.8254 m.
   (That is not the reference semispan, which is half of bref 3.6576 m, so 1.8288 m. The
   surface is placed with the wing's 5 deg dihedral, so y is not exactly the distance
   along the span.) If you see about **5.99** instead, the file is in feet. Go to Step 6a.
2. **Root extension.** A minimum y of **-0.02** means the surface reaches 20 mm through the
   symmetry plane (y = 0), as intended. The mesh cuts it off at y = 0. The cluster mesh job
   refuses a surface whose minimum y is not below -0.010 m.
3. **Closed, one part, one zone.** An open surface lets the mesh leak into the wing. Stop and
   go back to [02](./02_geometry_from_openvsp.md) if any of these three lines is different.

The log also says `Dumping bad quality faces to "badFaces"` and writes a file called
`badFaces` into the folder you ran from. The registered baseline surface does this too; it
is harmless. The minimum triangle quality it reports (`min 6.23599e-13`) is also normal for
this surface.

Optional, if you have the Python environment from [01_setup.md](./01_setup.md) (numpy
needed): the project's own case-level surface gate, described in
[02](./02_geometry_from_openvsp.md), section 10a. Run this one in your Python terminal
([01_setup.md](./01_setup.md), section A5), **not** in the Step 1 terminal: after the MPI fix
in Step 1, item 2, `python3` there is the system Python, which has no numpy (it stops with
`ModuleNotFoundError: No module named 'numpy'`). Set `REPO` and `CASE` in that terminal too.
The first line below checks for numpy and must print nothing:

```bash
python3 -c 'import numpy'
python3 $REPO/scripts/assert_surface_closed.py "$CASE" --require-overhang
```

The second should print `<your case path>: all surfaces closed, outward, and crossing y=0`.

### Step 6a. Only if the surface is in feet

OpenVSP works in feet. A surface that came straight from OpenVSP is in feet and must be
scaled by 0.3048. Scale from your original file into the case, then check again:

```bash
surfaceTransformPoints "scale=(0.3048 0.3048 0.3048)" "$STL" constant/triSurface/wing.stl
surfaceCheck constant/triSurface/wing.stl > log.surfaceCheck 2>&1
grep -E "^Bounding Box|is closed" log.surfaceCheck
```

Tested by first converting the baseline to feet (its bounding box then read
`(3.86 -0.0656168 -0.575466) (7.66489 5.98888 0.0448087)`) and scaling it back: the result
was the metre bounding box shown in Step 6, and still closed. A scaled file no longer has the
registered checksum, because the numbers are written out again. That is expected, but it
also tells you that you do not have the registered surface; for production, build the
surface with [02](./02_geometry_from_openvsp.md), whose scripts write metres directly and
reproduce the registered checksum. A scaled file is also slightly less
precise: `surfaceTransformPoints` writes coordinates to 6 significant figures (the
registered file carries 10). In the round trip tested here, points moved by up to
6.5 micrometres, the same order as the 4.5 micrometre production first layer. One more
reason to use the registered surface for production.

**If it does not:** if the maximum y now reads about **0.556** (tested:
`(0.358606 -0.006096 -0.0534625) (0.712092 0.556386 0.00416286)`), the file was already in
metres and you scaled it anyway. Copy the original again as in Step 5 and skip Step 6a.

---

## Step 7. Make the four tutorial-resolution edits

The production recipe cannot be meshed on a laptop. The delivered baseline mesh has
99,111,506 cells, and `recipes/README.md` records 68.6 million cells measured at 115 GB of
memory. So for this tutorial you change four files. **Nothing else changes.**

| File | Entry | Production (as shipped) | Tutorial |
|---|---|---|---|
| `system/meshLevels` | every refinement level | as shipped | **3 levels lower** (cells 8 times larger in each direction) |
| `system/snappyHexMeshDict` | `box_tref_a`, `box_tref_b` levels (the two wake boxes; their levels are written in this file, not in `meshLevels`) | 7 and 6 | 4 and 3 |
| `system/layerSettings` | `firstLayerThickness`, `minThickness` | 4.5e-06 m, 4.05e-06 m | 3.6e-05 m, 3.24e-05 m (8 times thicker, like the cells) |
| `system/decomposeParDict` | `numberOfSubdomains` | 192 (method `scotch`) | your `$NP` (method stays `scotch`) |

The surface cell becomes 3.906 mm instead of 0.488 mm (level 8 instead of 11, on a 1 m
background cell).

### 7a. Refinement levels

Replace `system/meshLevels` completely:

```bash
cat > system/meshLevels <<'EOF'
// TUTORIAL RESOLUTION: for learning the pipeline, not for results.
// Every level is the production value minus 3 (cells 8x larger in each direction).
surfaceLevelMin  8;   // production 11
surfaceLevelMax  8;   // production 11
levelOuter       2;   // production 5
levelMid         3;   // production 6
levelNear        4;   // production 7
levelWing        5;   // production 8
levelWake        4;   // production 7
levelCone        4;   // production 7
featureLevel     8;   // production 11
levelTEouter     8;   // production 11
levelTEinner     8;   // production 11
levelMorph       4;   // production 7
EOF
```

### 7b. The two wake boxes

```bash
sed -i 's/box_tref_a { mode inside; levels ((1e15 7)); }/box_tref_a { mode inside; levels ((1e15 4)); }/' system/snappyHexMeshDict
sed -i 's/box_tref_b { mode inside; levels ((1e15 6)); }/box_tref_b { mode inside; levels ((1e15 3)); }/' system/snappyHexMeshDict
```

You can make the same change in a text editor instead: it is the number after `1e15` on
the two lines starting `box_tref_a { mode inside;` and `box_tref_b { mode inside;`.

### 7c. The first prism layer

```bash
sed -i 's/^    firstLayerThickness 4.5e-06;/    firstLayerThickness 3.6e-05;/' system/layerSettings
sed -i 's/^    minThickness        4.05e-06;/    minThickness        3.24e-05;/' system/layerSettings
```

### 7d. The number of pieces for the parallel run

The shipped file already says `method scotch`, which can split a mesh into any number of
pieces, so only the number changes:

```bash
foamDictionary -writePrecision 12 -entry numberOfSubdomains -set $NP system/decomposeParDict
```

**What you should see:** `New entry numberOfSubdomains 8;` (with your `$NP`).

This is the first `foamDictionary` command in this guide, and three rules apply to every
one of them:

1. **Run it from inside the case folder, with a relative path** such as
   `system/decomposeParDict`. `foamDictionary` resolves every file path against the folder
   you are in, even a path that starts with `/`. So
   `foamDictionary -entry numberOfSubdomains -value $CASE/system/decomposeParDict` fails,
   from inside the case folder or anywhere else, with
   `file "<current folder>//<your case path>/system/decomposeParDict" does not exist`
   (tested). If you see that, run `cd "$CASE"` and give the path from there.
2. `foamDictionary -set` rewrites the whole file and removes its comments. That is fine in
   a case folder. Never run it on the files in `recipes/` themselves.
3. Always pass `-writePrecision 12` to `-set`: without it, `foamDictionary` rewrites every
   number in the file at 6 significant figures. (Tested: changing only `endTime` in
   `system/controlDict` without it shortened `liftDir` from `( -0.022687323 0 0.999742610 )`
   to `( -0.0226873 0 0.999743 )`. See also [06](./06_trim_to_target_cl.md), section 6.4.1.)

### 7e. Check that all four edits took

`sed` does nothing, silently, if the text it looks for is not there. So read the values
back:

```bash
foamDictionary -entry surfaceLevelMin -value system/meshLevels
foamDictionary -entry castellatedMeshControls/refinementRegions/box_tref_a/levels -value system/snappyHexMeshDict
foamDictionary -entry castellatedMeshControls/refinementRegions/box_tref_b/levels -value system/snappyHexMeshDict
foamDictionary -entry addLayersControls/firstLayerThickness -value system/layerSettings
foamDictionary -entry addLayersControls/minThickness -value system/layerSettings
foamDictionary -entry numberOfSubdomains -value system/decomposeParDict
foamDictionary -entry method -value system/decomposeParDict
sha256sum -c MANIFEST.sha256 | grep -v ": OK$"
```

**What you should see** (with `NP=8`):

```
8
( ( 1e+15 4 ) )
( ( 1e+15 3 ) )
3.6e-05
3.24e-05
8
scotch
sha256sum: WARNING: 4 computed checksums did NOT match
system/decomposeParDict: FAILED
system/layerSettings: FAILED
system/meshLevels: FAILED
system/snappyHexMeshDict: FAILED
```

(The `WARNING` line can also come after the four `FAILED` lines.)

**If it does not:** a value that still shows the production number means that edit did not
take. Open the file in a text editor and change it by hand.

The four `FAILED` lines are this case's fingerprint. **This case must never be sent to the
cluster.** Before it meshes, the cluster mesh job compares the recipe files with
`MANIFEST.sha256`. Two files are exempt from that comparison, `system/controlDict` and
`0.orig/U`, because they carry the angle of attack; for those two it checks instead that
the angle is consistent. It refuses a case where any other file differs, so it would refuse
this one. For the production mesh you start again from an untouched copy of the recipe
([04](./04_send_to_hpc.md), Step 5).

### Why these four edits and not fewer

1. **Levels.** Each level down halves the cell size in all three directions. Three levels
   down gave about 1.19 million cells, meshed in under two minutes on 8 cores.
2. **Wake boxes.** Their levels are written in `snappyHexMeshDict` itself. Left at 7 and 6,
   they would refine two large blocks of the wake to production cell size and undo most of
   the saving from 7a.
3. **First layer.** This is the edit people leave out, so it was tested both ways. With the
   production first layer (4.5 micrometres) on tutorial-size cells, snappyHexMesh **finished
   without complaint and added no layers at all**:
   1. three levels down: 0 of 11 layers, 0% of wall faces extruded (107 s);
   2. two levels down: an average of 0.0479 of 11 layers, 1.42% of the requested thickness,
      after 625 s and about 6.7 GB of extra memory.

   In the three-levels-down run, the first layer-addition pass rejected 1,069,929 faces as
   `faces on cells with determinant < 0.001`: a layer cell 4.5 micrometres thick under a
   3.9 mm cell is too flat to pass the quality limits. With the first layer scaled with the
   cells (36 micrometres), the last layer to surface cell ratio is the production one, about
   0.38, and the layers go in (Step 11). The price: the first cell is 8 times taller, so its
   y+ is about 8 times the production value at the same wall shear, and the tutorial mesh is
   not wall-resolved (Step 3). Another reason these forces mean nothing.
4. **Pieces.** The number of pieces must equal the number of processes you start. The
   shipped method, `scotch`, can split into any number of pieces, so only the number
   changes.

   **Only if your OpenFOAM-12 has no scotch.** `decomposePar` in Step 10 then stops with
   `You are trying to use scotch but do not have the scotch library loaded.` Use the
   `hierarchical` method instead. It needs three factors whose product is exactly `$NP`
   (here 2 x 2 x 2 = 8), and the shipped file has no `hierarchicalCoeffs` block to hold
   them. So set the method **and** write a complete block, both from inside the case folder:

   ```bash
   foamDictionary -writePrecision 12 -entry method -set hierarchical system/decomposeParDict
   foamDictionary -writePrecision 12 -entry hierarchicalCoeffs -set "{ n (2 2 2); delta 0.001; order xyz; }" system/decomposeParDict
   foamDictionary -entry hierarchicalCoeffs/n -value system/decomposeParDict
   ```

   **What you should see:** `New entry method          hierarchical;`, then the new block,

   ```
   New entry hierarchicalCoeffs
   {
       n               ( 2 2 2 );
       delta           0.001;
       order           xyz;
   }
   ```

   and last `( 2 2 2 )`. Then run Step 10 again. It prints
   `Selecting decomposer hierarchical`.

   Half of this is not enough. Setting only `hierarchicalCoeffs/n` stops with
   `keyword hierarchicalCoeffs is undefined`, and setting only the method makes
   `decomposePar` stop with `keyword n is undefined` (both tested).

   How this was tested: the Ubuntu package used for this guide has scotch, so the failure
   was reproduced by putting OpenFOAM-12's stand-in ("dummy") scotch library ahead of the
   real one. `decomposePar` then stopped with the message above. After the three commands
   it printed `Selecting decomposer hierarchical`. The mesh built from that split (with the
   stand-in library removed again, because the stand-in set also replaces the parallel
   library) had 1,186,226 cells, the `hierarchical` row of the table under "Why the rank
   count is part of the mesh recipe". So this changes the tutorial mesh a little, which is
   harmless here. The comment in the shipped `decomposeParDict` records the same failure on
   another cluster's build of OpenFOAM-12. The HPC12 build has scotch: all six delivered
   Condition CR meshes were split with it.

---

## Step 8. Build the background mesh: `blockMesh`

```bash
blockMesh > log.blockMesh 2>&1
grep -E "nCells|name:" log.blockMesh
```

`> log.blockMesh 2>&1` sends everything the program prints into a file called
`log.blockMesh`. Every step below does the same, so you can read the logs afterwards.

**What you should see** (less than a second):

```
  nCells: 108000
  patch 0 (start: 316800 size: 3600) name: symmetry
  patch 1 (start: 320400 size: 1800) name: inlet
  patch 2 (start: 322200 size: 1800) name: outlet
  patch 3 (start: 324000 size: 3600) name: outboard
  patch 4 (start: 327600 size: 3600) name: topBottom
```

That is the 60 x 30 x 60 m box of 1 m cubes. The wing is not in it yet.

**If it does not:**

1. `cannot find file ".../system/controlDict"`: you are not in the case folder. Run
   `cd "$CASE"` and try again.
2. `blockMesh: command not found`, or any other `FOAM FATAL`: go back to Step 1.

---

## Step 9. Extract the sharp edges: `surfaceFeatures`

The wing has a blunt trailing edge. snappyHexMesh can only keep its corners sharp if it is
told where they are.

```bash
surfaceFeatures > log.surfaceFeatures 2>&1
grep -E "^points|^edges|open edges|Writing featureEdgeMesh" log.surfaceFeatures
```

**What you should see:**

```
points      : 1273
edges       : 1275
    open edges                     :        0
Writing featureEdgeMesh to "constant/triSurface/wing.eMesh"
```

1273 points on the feature-edge mesh (`wing.eMesh`) is the value the project recorded for
B. The log also has a separate line `feature points : 4`; that is normal, do not confuse
the two. `open edges : 0` is one more confirmation that the surface is closed.

**If it does not:**

1. `Cannot read ".../constant/triSurface/wing.stl"`: the surface is missing or misnamed.
   Check Step 5: the file must be exactly `constant/triSurface/wing.stl`.
2. `Cannot read ".../constant/geometry/wing.stl"`: either a `constant/geometry/` folder
   exists (delete it; see Step 5, item 3), or there is no `constant/triSurface/` folder at
   all (Step 5 was not done). Both were tested and give this message.
3. `Triangles    : 0` near the top of the log: OpenFOAM read an empty `wing.stl`, most
   likely from a `constant/geometry/` folder. Delete that folder.
4. `open edges` is not 0: the surface is not closed. Stop; go back to Step 6.

---

## Step 10. Split the background mesh into pieces: `decomposePar`

```bash
decomposePar > log.decomposePar 2>&1
grep -E "Selecting decomposer|Max number of cells" log.decomposePar
ls -d processor*
```

**What you should see** (with `NP=8`):

```
Selecting decomposer scotch
Max number of cells = 13500 (0% above average 13500)
processor0  processor1  processor2  processor3  processor4  processor5  processor6  processor7
```

**If it does not:**

1. `Processor meshes exist but have no addressing.`: `processor*` folders from an earlier
   `snappyHexMesh` run are still there. Remove what the earlier attempt made (see
   "Starting again" near the end of this guide) and repeat from Step 8.
2. `Case is already decomposed with 8 domains, use the -force option` (with your
   numbers): the case was split before into a different number of pieces, for example
   after you changed `$NP`. If `snappyHexMesh` has not run yet, repeat the split with
   `-force`:

   ```bash
   decomposePar -force > log.decomposePar 2>&1
   ```

   If it has, use "Starting again".

`decomposePar` does not always stop on old `processor*` folders. If you repeat Step 8
(`blockMesh`) after a finished `snappyHexMesh` run without cleaning up, `decomposePar`
splits the new box without complaint but leaves files from the old mesh in
`processor*/constant/polyMesh`. `snappyHexMesh` in Step 11 then stops within seconds with
`Number of cells in mesh:13500 does not equal size of cellLevel:...` (tested). So before
re-meshing, always run the "Starting again" command.

---

## Step 11. Build the mesh: `snappyHexMesh` in parallel

This is the long step. It refines the box towards the wing, cuts the wing out, snaps the cell
faces onto the surface, then grows prism layers on the wing.

```bash
mpirun -np $NP snappyHexMesh -parallel -overwrite > log.snappyHexMesh 2>&1
```

If you want to watch it, open a second terminal and run this (use your own path if you chose
another case name in Step 3; Ctrl+C there stops the watching, not the mesher):

```bash
tail -f $HOME/argus_cases/B_CR_tutorial/log.snappyHexMesh
```

When it has finished, read the summary lines:

```bash
grep -E "^Refined mesh :|^Snapped mesh :|^Layer mesh :|^Finished meshing" log.snappyHexMesh
grep -A4 "^patch faces    layers   overall thickness" log.snappyHexMesh
grep "^Extruding" log.snappyHexMesh | tail -1
grep -c "reached limit" log.snappyHexMesh
```

**What you should see** (107 s on 8 cores):

```
Refined mesh : cells:570083  faces:1841888  points:702209
Snapped mesh : cells:570083  faces:1824444  points:690852
Layer mesh : cells:1193347  faces:3746782  points:1366326
Finished meshing with 54 illegal faces (concave, zero area or negative cell pyramid volume)
Finished meshing in = 107.19933 s.
patch faces    layers   overall thickness
                       [m]       [%]
----- -----    ------   ---       ---
wing 82922    7.52     0.00419   89.7    

Extruding 78615 out of 82922 faces (94.805962%). Removed extrusion at 0 faces.
0
```

Your counts may differ a little. This case, with the same files, on 8 cores with `scotch`,
was meshed six times: five runs gave exactly the output above, and one gave 1,195,711
cells, 82,924 wing faces, 7.54 layers and 51 illegal faces. Parallel snappyHexMesh does
not always reproduce itself to the last cell. A difference of a few thousand cells is that;
a difference of hundreds of thousands is not.

How to read it:

1. `Finished meshing` must be there. If it is missing, the mesher did not finish (see below).
2. **Layers.** The wing has 82,922 faces. On average 7.52 of the 11 requested layers were
   built, the achieved stack is 0.00419 m thick on average, and 94.8% of the wing faces
   carry layers. The 89.7% is the average, over the wing faces, of achieved thickness
   divided by the thickness wanted at that face. (The request table printed earlier in the
   log gives the wanted overall thickness as 0.00458 m on average; the nominal 11-layer
   stack at tutorial settings is 4.69 mm.) These **achieved** values can only be at or
   below the request. Make sure you are reading the right table: the request table printed
   earlier starts with the same `patch faces    layers` words, but it always shows exactly 11
   layers and has no `[%]` column. The `grep` above picks the result table.
3. **54 illegal faces**, out of 3.75 million, are faces that did not meet every quality
   limit and were left in place. A few tens is normal for this recipe at this resolution.
4. `reached limit` must be **0**. A non-zero count means the cell cap (`maxGlobalCells`)
   stopped refinement part-way, and whatever region refines last silently gets nothing.
   That happened on the project's earlier production meshes; see `recipes/README.md`. Do not
   lower `maxGlobalCells`.

**If it does not:**

1. **Every rank prints `attempt to run parallel on 1 processor`:** the wrong `mpirun`.
   Because the output goes into the log, you see nothing on the screen. In the tests for
   this guide the command returned in about 1 second, but [01_setup.md](./01_setup.md),
   section A5, records a run that hung instead. Either way, check the log (from a second
   terminal in the case folder if the command has not returned):

   ```bash
   grep -c "attempt to run parallel" log.snappyHexMesh
   ```

   Anything other than 0 means this problem (it printed 8, one per rank, in the test):
   press Ctrl+C if the command is still running, then go back to Step 1, item 2.
2. **`There are not enough slots available in the system`:** `$NP` is larger than your
   number of physical cores. Go back to Step 2, choose a smaller number, redo Step 7d, and
   then split again with `-force` instead of the plain Step 10 command (plain
   `decomposePar` refuses a case that is already split into a different number of pieces):

   ```bash
   decomposePar -force > log.decomposePar 2>&1
   ```
3. **`... specifies 8 processors but job was started with 4 processors`**, or
   **`number of processor directories = 8 is not equal to the number of processors = 10`**
   (with your numbers): `-np` does not match `numberOfSubdomains`. The first message means
   `-np` was too small, the second that it was too large. Use the same `$NP` for both.
4. **No `Finished meshing` line.** Look at the end of the log:

   ```bash
   tail -30 log.snappyHexMesh
   ```

   The usual cause when meshing too many cells on too small a machine is that the operating
   system kills the processes for lack of memory: the log simply stops, and `mpirun`
   reports that a process was killed (signal 9). The recipe's own comments record exactly
   that for an early version of the production mesh ("OOM-killed the run at 21.4 M cells
   (signal 9)"). Not reproduced for this guide (NOT VERIFIED). Check that you made the edits
   in Step 7, and that this shows enough free memory (in GB, `available` column):

   ```bash
   free -g
   ```
5. **Layers near zero** (`layers` close to 0, `Extruding` a few per cent or less): the first
   layer is too thin for the cell size. Check Step 7c. You will also see, during layer
   addition, lines such as
   `faces on cells with determinant < 0.001 : 1069929`.
6. **The cell count is far from 1.19 million.** Differences of up to about 1% are normal
   if you used a different `$NP` or splitting method (see "Why the rank count is part of the
   mesh recipe" below). A mesh of about 730,000 cells whose far-field patches later vanish
   in `checkMesh` means `locationInMesh` is inside the wing (Troubleshooting, at the end).

---

## Step 12. Collect the mesh: `reconstructPar -constant`

After `snappyHexMesh -parallel -overwrite`, the new mesh exists only inside the
`processor*` folders. The mesh in `constant/polyMesh` is still the empty box from Step 8.
This step joins the pieces.

```bash
reconstructPar -constant > log.reconstructPar 2>&1
grep "Reconstructing meshes" log.reconstructPar
grep -A4 "^    wing" constant/polyMesh/boundary
```

**What you should see** (6 s):

```
Reconstructing meshes
    wing
    {
        type            wall;
        inGroups        List<word> 1(wall);
        nFaces          82922;
```

The box from Step 8 has no `wing` patch, so seeing it here proves you have the snapped
mesh, not the box.

**If it does not:**

1. No `wing` block: the mesh was not collected. Run the `reconstructPar` line again and read
   `log.reconstructPar`.
2. Do not use `reconstructParMesh`. In OpenFOAM-12 it only prints a message saying it has
   been replaced by `reconstructPar`, and builds nothing.

---

## Step 13. Check the mesh: `checkMesh`

```bash
checkMesh -allGeometry -allTopology -constant > log.checkMesh 2>&1
echo "checkMesh exit status: $?"
grep -E "^ +cells:|Number of regions|\*\*\*|^Failed|^Mesh OK" log.checkMesh
```

Use both `-allGeometry` and `-allTopology`: plain `checkMesh` leaves some checks out. This
is the same command the cluster mesh job runs.

**What you should see** (28 s):

```
checkMesh exit status: 0
    cells:            1193347
    Number of regions: 1 (OK).
 ***Max skewness = 5.2082764, 3 highly skew faces detected which may impair the quality of the results
 ***Error in face tets: 43 faces with low quality or negative volume decomposition tets.
 ***Cells with small determinant (< 0.001) found, number of cells: 114109
 ***Concave cells (using face planes) found, number of cells: 64725
Failed 4 mesh checks.
```

(The one run that gave 1,195,711 cells failed the same four checks, with slightly different
counts.)

Note the exit status: **`checkMesh` returns 0 even when it reports failed checks.** A script
that only looks at the exit status will never see a failure. Always read the log.

### Which failures are expected

The four above are expected. They come from how snappyHexMesh builds cells around the prism
layers and the snapped surface:

1. `Max skewness` (here 3 faces, maximum 5.21 against the limit of 4);
2. `Error in face tets`;
3. `Cells with small determinant`;
4. `Concave cells`.

**The production meshes fail the same four, and only those.** The `checkMesh` logs of all
six delivered Condition CR meshes (HPC12 job logs; they are not in this repository) each
end `Failed 4 mesh checks.`, and the four are exactly the ones above:

| Failed check | This tutorial mesh | Delivered baseline B mesh | All six delivered meshes |
|---|---|---|---|
| `Max skewness` | 5.21, 3 faces | 15.84, 7 faces | 8.86 to 17.35, 4 to 29 faces |
| `Error in face tets` | 43 faces | 3,162 faces | 3,162 to 3,572 faces |
| `Cells with small determinant` | 114,109 cells | 13,020,715 cells | 13.00 to 13.05 million cells |
| `Concave cells` | 64,725 cells | 2,106,988 cells | 2.10 to 2.12 million cells |

The delivered meshes have 99,111,506 (B) to 99,208,470 cells. Their one-star warnings are
the same four kinds as this tutorial's (listed below). The run cards
(`data/run_cards/SW*_r2_trim.json`) record the verdict as `checkMesh_pass: false`.

Two more numbers from those logs appear in the run cards, and both are fine:

1. `mesh.max_non_orthogonality` is the largest face non-orthogonality, in degrees, from the
   log line `Mesh non-orthogonality Max:`: 70.957803 for B, and 70.96 to 71.97 over the six
   (70.35 on this tutorial mesh). Faces above 70 deg are only counted in a one-star warning,
   `*Number of severely non-orthogonal (> 70 degrees) faces` (1 to 4 faces on the delivered
   meshes, 1 here). A real non-orthogonality of 90 deg or more would have added a
   `***Number of non-orthogonality errors` failure.
2. `mesh.max_skewness` is the `Max skewness` value in the table: 15.844881 for B.

The maximum aspect ratio, which the cards do not carry, is printed as `OK`: 112.81923 on
five of the delivered meshes, 113.77236 on W, and 110.51766 on this tutorial mesh.

### Which failures are not expected (stop and find the cause)

1. `Number of regions` other than 1 (it may be printed with one star, `*Number of regions:
   13`, and then it is not counted in `Failed n mesh checks`: read it anyway).
2. Any other line starting `***`, for example
   `Zero or negative cell volume`, `incorrectly oriented` face pyramids,
   `Number of non-orthogonality errors`, or open cells.
3. `cells: 108000`: you checked the empty box. Step 12 did not run.
4. The wing patch missing, or the far-field patches (`inlet`, `outlet`, `outboard`,
   `topBottom`) missing: see "locationInMesh inside the wing" under Troubleshooting.

Fewer failed checks is not automatically good news. A mesh of the wrong domain can fail
fewer checks than the right one (Troubleshooting, "locationInMesh inside the wing"). Check the
cell count and the number of regions first, then the list of failures.

Four more lines in the log are warnings, not failures (one star, not three), and are normal
here. To see them:

```bash
grep -E "^ +\*[A-Za-z]" log.checkMesh
```

```
   *Number of severely non-orthogonal (> 70 degrees) faces: 1.
   *Edges too small, min/max edge length = 3.5458763e-05 1.0000043, number too small: 178562
   *There are 1621 faces with concave angles between consecutive edges. Max concave angle = 79.484861 degrees.
   *There are 14 faces with ratio between projected and actual area < 0.8
```

(The maximum non-orthogonality on this mesh is 70.35 deg.) If this `grep` also prints a
`*Number of regions` line, see the list above.

A useful last look in the same log: the six patches, and the wing patch's bounding box.

```bash
grep -E "^ +(symmetry|inlet|outlet|outboard|topBottom|wing) +[0-9]+ +[0-9]+ +ok" log.checkMesh
grep -E "^ +wing \(" log.checkMesh
```

```
                symmetry     9908    10586  ok (non-closed singly connected)
                   inlet     1800     1891  ok (non-closed singly connected)
                  outlet     1800     1891  ok (non-closed singly connected)
                outboard     3600     3721  ok (non-closed singly connected)
               topBottom     3600     3782  ok (non-closed singly connected)
                    wing    82922    89224  ok (non-closed singly connected)
                        wing (1.1768103 0 -0.17539476) (2.3362596 1.8254118 0.013652058)
```

All six patches must be there. The minimum y of the wing patch is 0: the 20 mm root extension
was cut off at the symmetry plane, as intended.

---

## Step 14. Inspect the mesh

### 14a. From the logs

How many cells ended up at each refinement level:

```bash
grep -A10 "^Layer mesh :" log.snappyHexMesh
```

```
Layer mesh : cells:1193347  faces:3746782  points:1366326
Cells per refinement level:
    0	107727
    1	1545
    2	3091
    3	11477
    4	35770
    5	8679
    6	23518
    7	74602
    8	926938
```

Level 0 is the 1 m background; level 8 (3.906 mm) is the wing surface and its layers.

### 14b. In ParaView

Do this **before** the optional Step 15: Step 15b rebuilds the `processor*` folders, which
hold the layer information shown below.

```bash
paraFoam -builtin -touch
```

This prints `Created 'B_CR_tutorial.foam'`. Open that file in ParaView, as described in
[01_setup.md](./01_setup.md), section A7. Then, in the reader's properties panel:

NOT VERIFIED in ParaView's graphical interface (it was not opened for this guide). The same
reader, with the same settings, was loaded from ParaView's Python interface (`pvpython`) on
this case, and gave the result described in item 1.

1. **Layer coverage on the wing.** Set `Case Type` to `Decomposed Case`, untick
   `Skip Zero Time` (it may be under the advanced properties, the gear icon), tick only
   `patch/wing` under `Mesh Regions`, tick `thicknessFraction` under `Cell Arrays`, and press
   Apply. Colour by `thicknessFraction`: 1 is a full layer stack, 0 is no layers. On this
   case its range was 0 to 1 over the 82,922 wing faces.
2. **The volume mesh.** Set `Case Type` back to `Reconstructed Case`, tick `internalMesh`,
   press Apply, add a `Slice` filter with normal (0, 1, 0) and origin y = 0.9127 m (half way
   out along the wing, in the frame of Step 6), and show it as `Surface With Edges`. You
   should see the cells shrink towards the wing in steps, the wake refinement behind it,
   and the thin layers on the wing surface.

---

## Step 15 (optional). Does it start? A smoke test, not a result

This runs the solver for 50 iterations to prove that the case is complete: fields, boundary
conditions, numerics and force output all work. **It is a smoke test, not a result.** The
mesh is the tutorial mesh, the angle is the recipe's starting value (1.30 deg), not a trim,
and 50 iterations is nowhere near converged (a production solve runs thousands).

### 15a (optional). Add the post-processing block

A production case carries one file more than the recipe, `system/argusPostPro`, and one
more line in `system/controlDict` that includes it. Together they make the solver write what
the report needs ([07](./07_postprocessing_and_report.md)): the pressure and viscous parts
of the forces; the pressure, pressure coefficient (cp), y+ and wall shear stress on the wing
surface; four Trefftz planes behind the wing; and the residuals. For production you install
it on your own machine before the case goes to HPC12 ([04](./04_send_to_hpc.md), Step 5e).
This tutorial case never goes there, so here it is optional: do it if you want to see what
those files look like. It changes neither the mesh nor the forces of Step 15b.

Run these two commands in your Python terminal, as in Step 6 (numpy needed; set `REPO` and
`CASE` there too):

```bash
python3 $REPO/scripts/install_postpro.py "$CASE"
python3 $REPO/scripts/assert_case_postpro.py --before "$CASE"
```

**What you should see** (the first line is your case folder's name):

```
  B_CR_tutorial
     solver        foamRun -> rho mode 'rhoInf', rhoInf 1, pInf 0
     wall patch    wing           (from constant/polyMesh/boundary)
     CofR          1.34 0 0       (from existing forces object)
     turbulence    kOmegaSST      -> residual fields: U p k omega
     |U| for cp    34.0000        (from 0.orig/U (via $Uinf))
     TE x 2.3363, chord 1.1597 (from ASCII STL wing.stl)
     Trefftz x     2.3942 2.6262 2.9161 3.2061
     controlDict now includes argusPostPro
```

Then the check prints a line with your case path and `(before)`, eight `ok` lines, and
ends with:

```
  ---- 8 of 8 requirements satisfied ----
```

How to read the numbers:

1. `TE x` and `Trefftz x` are in metres, in the frame of `wing.stl` (the frame of Step 6),
   with x pointing downstream.
2. `chord` 1.1597 m is the streamwise length of `wing.stl`, not a reference chord. The four
   Trefftz planes sit 0.05, 0.25, 0.50 and 0.75 of that length behind the trailing edge.
3. `|U| for cp` is the Condition CR speed in m/s.
4. `rhoInf 1, pInf 0` is right for this incompressible case, where p is kinematic (pressure
   divided by density).
5. `wall patch ... (from constant/polyMesh/boundary)`: the script found the wing patch in
   the mesh of Step 12. On a production case, installed before meshing, this line reads
   `(from existing forces object)` instead.

Running the install a second time changes nothing; it prints
`controlDict already includes argusPostPro`. `sha256sum -c MANIFEST.sha256` now reports
five changed files, `system/controlDict` beside the four of Step 7e
(`sha256sum: WARNING: 5 computed checksums did NOT match`).

**If it does not:**

1. `ModuleNotFoundError: No module named 'numpy'`: this terminal's `python3` has no numpy.
   Use your Python terminal. Nothing was written.
2. `HALT <case>: no sampled wall surface and no readable ASCII STL ...`:
   `constant/triSurface/wing.stl` is missing. Do Step 5 first. Nothing was written.

### 15b. Run 50 iterations

Only in this tutorial case, lower the iteration count from the recipe's 4000 to 50:

```bash
foamDictionary -writePrecision 12 -entry endTime -set 50 system/controlDict
grep -E "^endTime|liftDir" system/controlDict
```

**What you should see:** after the `New entry endTime` message,

```
endTime         50;
        liftDir         ( -0.022687323 0 0.99974261 );
```

The lift direction kept all its digits (see Step 7d, rule 3). Check digits with `grep`, or give
`foamDictionary -value` the same `-writePrecision 12`: without it, `-value` prints 6
significant figures (`( -0.0226873 0 0.999743 )`) whatever the file holds. With it,

```bash
foamDictionary -writePrecision 12 -entry functions/forceCoeffs1/liftDir -value system/controlDict
```

prints `( -0.022687323 0 0.99974261 )`.

Then create the starting fields, split the **finished** mesh and the fields into pieces,
initialise the flow, and run:

```bash
cp -r 0.orig 0
decomposePar -force > log.decomposePar.solve 2>&1
mpirun -np $NP potentialFoam -parallel > log.potentialFoam 2>&1
mpirun -np $NP foamRun -parallel > log.foamRun 2>&1
```

`potentialFoam` computes a smooth starting flow. It is not optional for this recipe: a
comment in its `system/fvSolution` records a full-size (124.5 million cell) wall-resolved
case that started from a uniform flow and died with a floating-point exception at
iteration 7. On the cluster, the solve job runs it for you
([05](./05_submit_and_run.md)).

**What you should see.** `potentialFoam` (11 s to 32 s in four runs):

```bash
grep -E "Continuity error|Interpolated velocity error|^End" log.potentialFoam
```

```
Continuity error = 6.0371439e-07
Interpolated velocity error = 8.6248765e-05
End
```

`foamRun` (80 s to 175 s for the 50 iterations in five runs on 8 cores):

```bash
grep -E "^Selecting solver|^Selecting RAS" log.foamRun
tail -4 log.foamRun
```

```
Selecting solver incompressibleFluid
Selecting RAS turbulence model kOmegaSST

End

Finalising parallel run
```

Each iteration in `log.foamRun` starts with a line like `Time = 50s` followed by one
residual line per equation, for example (from an earlier run of this case)
`smoothSolver:  Solving for Ux, Initial residual = 8.9388352e-05, ...` and
`GAMG:  Solving for p, Initial residual = 8.6244869e-05, ...`.

And the force output, written every 5 iterations:

```bash
head -9 postProcessing/forceCoeffs1/0/forceCoeffs.dat
tail -3 postProcessing/forceCoeffs1/0/forceCoeffs.dat
```

```
# Force coefficients
# liftDir       : (-2.26873230e-02 0.00000000e+00 9.99742610e-01)
# dragDir       : (9.99742610e-01 0.00000000e+00 2.26873230e-02)
# pitchAxis     : (0.00000000e+00 1.00000000e+00 0.00000000e+00)
# magUInf       : 3.40000000e+01
# lRef          : 3.93957000e-01
# Aref          : 6.20462000e-01
# CofR          : (1.34000000e+00 0.00000000e+00 0.00000000e+00)
# Time          	Cm              	Cd              	Cl              	Cl(f)           	Cl(r)           
40              	-3.65244492e-01	1.79029725e-02	3.25407248e-01	-2.02540868e-01	5.27948116e-01
45              	-3.75833774e-01	1.80293169e-02	3.37318488e-01	-2.07174530e-01	5.44493017e-01
50              	-3.84468891e-01	1.80905865e-02	3.47077567e-01	-2.10930108e-01	5.58007675e-01
```

The frame of these numbers is in the header: coefficients on the **half-model** reference
area Aref 0.620462 m2 (half of the DSO-basis Sref 1.24092 m2), reference length 0.393957 m,
speed 34.0 m/s, moments about (1.34, 0, 0) m, at 1.30 deg, on the **tutorial mesh**, after
50 iterations. Cl is still climbing by about 0.01 every 5 iterations. That the file exists
and fills in is the whole point of this step. **The numbers are not a result and are not
comparable to anything in `data/`.**

**If you did Step 15a**, the force coefficients are the same to the last digit (tested), and
the block's own output is there too:

```bash
ls postProcessing
grep -a -m4 -E "^(p|cp|yPlus|wallShearStress) [0-9]+ [0-9]+ float" postProcessing/argusWingSurface/0/wingSurface.vtk
ls postProcessing/argusTrefftz/0
```

```
argusForces  argusResiduals  argusTrefftz  argusWallShearStress  argusWingSurface  argusYPlus  forceCoeffs1
p 1 82922 float
cp 1 82922 float
yPlus 1 82922 float
wallShearStress 3 82922 float
trefftz_x1.vtk  trefftz_x2.vtk  trefftz_x3.vtk  trefftz_x4.vtk
```

1. 82922 is the number of wing faces (Step 12). The `3` after `wallShearStress` shows that
   the shear stress was written as a vector.
2. The wing-surface and Trefftz files are written at iteration 0 and at every saved
   iteration. This run saves nothing ("Also expected", item 3, below), so there is only the
   `0` folder.
3. `log.foamRun` shows `wallShearStress argusWallShearStress:`, then
   `processing wall patches:` and `wing`, then two `Reading surface description:` blocks: the
   first lists `wingSurface`, the second `trefftz_x1` to `trefftz_x4`. The block adds no
   warnings.
4. Do not judge this case with `python3 $REPO/scripts/assert_case_postpro.py --after "$CASE"`.
   It reports `---- 5 of 8 requirements satisfied ----`, with `wing_surface_p`,
   `wing_surface_tau_VECTOR` and `wing_surface_yplus` marked `MISSING` (tested). Those three
   lines are wrong: the check looks for each field name only in the first 20,000
   characters of the file, and OpenFOAM-12 writes the field names after the list of points
   and faces. The `grep` above is the check to use.

Also expected:

1. Both logs print a warning `Removing patchGroup 'symmetry' which clashes with patch 0 of
   the same name.` Harmless: the patch called `symmetry` is also of type `symmetry`.
2. `log.foamRun` prints `bounding k, min: ... max: ...` lines, here from iteration 14 on,
   in 35 of the 50 iterations (from iteration 13, in 37, in the earlier run). It means the
   turbulence energy went negative in some cells and was clipped. On this coarse mesh in the
   first iterations it does not stop the run.
3. No new time folder (such as `processor0/50`) is written: fields are saved every 500
   iterations (`writeInterval 500`), and the run stopped at 50. `ls processor0` shows only
   `0` and `constant`.

**If it does not:**

1. `potentialFoam` stops with `keyword Phi is undefined in dictionary "IOstream/solvers"`:
   your `system/fvSolution` is not the recipe's. Copy it again from the recipe. Note that a
   missing `potentialFlow` block does **not** stop `potentialFoam` (tested, serial and on 4
   ranks): it runs silently with the `SIMPLE` block's non-orthogonal corrections (1 in this
   recipe) instead of the `potentialFlow` block's 10, whatever the comment in the file says
   (Step 3, the list of confusing comments, item 4). So check the file against the
   manifest:

   ```bash
   sha256sum -c MANIFEST.sha256 | grep fvSolution
   ```

   It must print `system/fvSolution: OK`.
2. `potentialFoam` stops with `cannot find file ".../processor0/0/U"`: the `0` folder was
   missing when you ran `decomposePar`. `decomposePar` does not complain about that; it
   prints `(no FV fields)` and splits only the mesh (tested). Create the fields and split
   again, then carry on with the two `mpirun` lines:

   ```bash
   cp -r 0.orig 0
   decomposePar -force > log.decomposePar.solve 2>&1
   ```
3. `foamRun` stops with `Floating point exception` in the first iterations: check that
   `potentialFoam` ran and ended with `End`.

---

## What changes for the production mesh

| | This tutorial | Production |
|---|---|---|
| Where | your machine | HPC12, through [04](./04_send_to_hpc.md) and [05](./05_submit_and_run.md) |
| Recipe files | four edited (Step 7) | **unedited**: `sha256sum -c MANIFEST.sha256` says `OK` for every file except `system/controlDict`, which changes when the post-processing block is installed, and `0.orig/U`, which changes (with `system/controlDict`) when you set the angle of attack ([06](./06_trim_to_target_cl.md)) |
| Surface cell | level 8, 3.906 mm | level 11, 0.488 mm |
| Wake boxes `box_tref_a`, `box_tref_b` | levels 4 and 3 | levels 7 and 6 |
| First prism layer | 36 micrometres | 4.5 micrometres (11 layers, ratio 1.45) |
| Pieces for `snappyHexMesh` | your `$NP` (8 here), `scotch` | 192, `scotch`: the shipped `decomposeParDict`, unchanged. All six delivered Condition CR meshes were built that way: their run cards record ranks 192 and method `scotch`, and their `snappyHexMesh` logs report `nProcs : 192` (HPC12 job logs) |
| Cells | about 1.19 million | 99,111,506 for the delivered baseline mesh (`data/run_cards/SWB_r2_trim.json`), 99,111,506 to 99,208,470 over the six. A new build from the same recipe has not been compared with these counts; Step 11 shows that parallel `snappyHexMesh` does not always reproduce itself to the last cell |
| Layers on the wing | 7.52 of 11 on average, over 82,922 wing faces; 89.7% of the wanted thickness | 9.79 of 11 on average, over 5,276,772 wing faces; 95.2% of the wanted thickness, a 0.000557 m stack on average (delivered baseline mesh, HPC12 job logs) |
| `checkMesh` | the four expected failures (Step 13) | the same four on all six delivered meshes (Step 13) |
| Memory | about 2.2 GB | 68.6 million cells measured at 115 GB (`recipes/README.md`); the header of `hpc/mesh_hpc12.pbs` asks for six 32-core typej nodes |
| Mesh job | the commands in Steps 8 to 13 | `hpc/mesh_hpc12.pbs`, which runs the same programs |
| Post-processing block | optional (Step 15a) | required: installed on your machine before the case goes to HPC12 ([04](./04_send_to_hpc.md), Step 5e) |
| Solve | 50 iterations on `$NP` ranks (Step 15b) | `hpc/stage_solve_case.sh` makes a solve case from the finished mesh, split into 128 pieces (`hierarchical`, 8 x 2 x 8), and `hpc/solve_hpc12.pbs` runs it. The first leg runs to `endTime` 4000 with fields saved every 500 iterations; the six delivered solves ran on to iteration 6000 ([05](./05_submit_and_run.md)) |

The cluster mesh job runs the same programs, with `surfaceFeatures` before `blockMesh` (the
reverse of Steps 8 and 9; both orders work) and `decomposePar -force` in place of Step 10's
plain `decomposePar`: `surfaceFeatures`, `blockMesh`, `decomposePar -force`,
`snappyHexMesh -parallel -overwrite` on 192 ranks, `reconstructPar -constant`, and
`checkMesh -allGeometry -allTopology -constant`. Before it meshes, it refuses a case whose
recipe files do not match `MANIFEST.sha256` (apart from `system/controlDict` and
`0.orig/U`, which carry the angle of attack and are checked for a consistent angle instead;
see [04](./04_send_to_hpc.md), Step 5d), and a surface that is open or does not reach 10 mm
below y = 0. `system/argusPostPro` is not in the manifest, so the block does not trouble
either check. Read the job's header comments: they explain each check.
(`hpc/README_HPC.md` and `hpc/submit.sh` describe an earlier setup on another cluster; they
are not the HPC12 route.)

### Why the rank count is part of the mesh recipe

A *rank* is one of the parallel processes a run is split into; there is one rank per piece
(`mpirun -np $NP` starts `$NP` ranks). snappyHexMesh in parallel refines, balances and snaps
each piece of the domain, so the final mesh depends on how the domain was split. It is not
only a question of speed. The header of `hpc/mesh_hpc12.pbs` says it plainly: "192 RANKS IS
PART OF THE RECIPE, NOT A SCHEDULING CHOICE".

This guide measured it on the tutorial case: the same files, meshed with different splits.
"Repeats" says how often the same split reproduced itself.

| Split | Cells | Wing faces | Average layers | Repeats |
|---|---|---|---|---|
| `hierarchical` (2 2 2), 8 pieces | 1,186,226 | 82,921 | 7.52 of 11 | 3 of 3 runs identical |
| `hierarchical` (2 2 1), 4 pieces | 1,194,174 | 82,923 | 7.51 of 11 | 1 run |
| `scotch`, 8 pieces | 1,193,347 | 82,922 | 7.52 of 11 | 5 of 6 runs identical; the sixth gave 1,195,711 |
| `scotch`, 4 pieces | 1,193,097 | 82,923 | 7.51 of 11 | 1 run |

What this shows:

1. With `hierarchical`, which reproduced itself exactly, going from 8 to 4 pieces changed
   the mesh by 7,948 cells (0.67%). The number of pieces changes the mesh.
2. At 8 pieces, changing only the method changed the mesh by 7,121 cells (0.60%). The
   method changes the mesh too.
3. With `scotch`, the 4-piece and 8-piece meshes differ by only 250 cells, less than the
   2,364-cell difference between two `scotch` runs on 8 pieces. So with `scotch` the effect
   of the piece count was not resolved here.

Small, but not zero, and larger than the repeat-to-repeat spread. Two geometries meshed with
different splits carry a meshing difference inside their difference. So for any mesh whose
numbers you will compare with the delivered ones, use the production split: 192 pieces with
`scotch`. That is exactly what the shipped `decomposeParDict` says, and how all six delivered
Condition CR meshes were built (`mesh.mesh_generation` in `data/run_cards/SW*_r2_trim.json`:
ranks 192, method `scotch`, machine hpc12; their `snappyHexMesh` logs report `nProcs : 192`).
So send the case with that file unedited ([04](./04_send_to_hpc.md)). Solving is different
again: re-splitting a finished mesh does not change it, so the solver can run on another
number of pieces. The six delivered Condition CR solves ran on 128 pieces, `hierarchical`
(8 x 2 x 8), and `hpc/stage_solve_case.sh` writes that split for you
([05](./05_submit_and_run.md)).

---

## Starting again

To re-mesh this case from the beginning (for example after changing `system/meshLevels`),
remove what the mesh steps made, then repeat from Step 8:

```bash
cd "$CASE"
rm -rf processor* constant/polyMesh constant/extendedFeatureEdgeMesh constant/triSurface/wing.eMesh 0 postProcessing
```

Or delete the whole case folder (`rm -rf "$CASE"`) and start from Step 3. Deleting the case
removes about 440 MB (450 MB if you did Step 15a).

---

## Troubleshooting

| What you see | Most likely cause | What to do |
|---|---|---|
| `bash: foamRun: command not found`, or `foamVersion` does not print `OpenFOAM-12` | OpenFOAM-12 not loaded in this terminal (or not installed) | Step 1, item 1 |
| `attempt to run parallel on 1 processor` from every rank | Anaconda's (or another) `mpirun` is used | Step 1, item 2 |
| `There are not enough slots available in the system` | more ranks than physical cores | Step 2 |
| `... specifies N processors but job was started with M processors` | `-np` smaller than `numberOfSubdomains` | use `$NP` for both (Steps 7d and 11) |
| `number of processor directories = N is not equal to the number of processors = M` | `-np` larger than `numberOfSubdomains` | use `$NP` for both (Steps 7d and 11) |
| `cannot find file ".../system/controlDict"` | you are not in the case folder | `cd "$CASE"` |
| surfaceCheck bounding box reaches y = 5.99 | the surface is in feet | Step 6a |
| surfaceCheck bounding box reaches only y = 0.556 | a surface already in metres was scaled again | Step 6a, "If it does not" |
| surfaceCheck says the surface is not closed, or more than one part | broken or wrong surface | stop; [02](./02_geometry_from_openvsp.md) |
| `Cannot read ".../constant/triSurface/wing.stl"` | surface missing or misnamed | Step 5 |
| `Cannot read ".../constant/geometry/wing.stl"` | a `constant/geometry/` folder exists, or there is no `constant/triSurface/` folder | Step 9, "If it does not", item 2 |
| `Triangles    : 0` in `log.surfaceFeatures` | an empty `wing.stl` in `constant/geometry/` was read first | delete `constant/geometry/` |
| no `Finished meshing` line | the mesher died, often for lack of memory | Step 11, item 4 |
| layers about 0 | first layer too thin for the cell size | Step 7c |
| `checkMesh` says `cells: 108000` | you checked the empty box | Step 12 |
| about 730,000 cells, `Number of regions: 13`, no `inlet`/`outlet`/`outboard`/`topBottom` patches, while `checkMesh` fails **fewer** checks than expected (or plain `checkMesh` says `Mesh OK.`) | `locationInMesh` is inside the wing | see below |
| any `***` line in `log.checkMesh` other than the four in Step 13 | a real mesh problem | stop; find the cause before solving |
| `Processor meshes exist but have no addressing.` from `decomposePar` | `processor*` folders left from an earlier attempt | "Starting again", then Step 8 |
| `Case is already decomposed with N domains, use the -force option` | the case was split before into a different number of pieces | Step 10, "If it does not", item 2 |
| `Number of cells in mesh:... does not equal size of cellLevel:...` from `snappyHexMesh` | files from an earlier mesh left in `processor*` | "Starting again", then Step 8 |
| `You are trying to use scotch but do not have the scotch library loaded.` | OpenFOAM built without scotch | Step 7, "Why these four edits", item 4 |
| `keyword hierarchicalCoeffs is undefined` from `foamDictionary`, or `keyword n is undefined` from `decomposePar` | `hierarchical` chosen without a complete `hierarchicalCoeffs` block | Step 7, "Why these four edits", item 4 |
| `file "<current folder>//<your path>" does not exist` from `foamDictionary` | a path starting with `/` (such as `$CASE/system/...`) was given | Step 7d, rule 1: `cd "$CASE"`, then use `system/...` |
| `ModuleNotFoundError: No module named 'numpy'` | a Python script was run in the OpenFOAM terminal | use your Python terminal (Step 6, Step 15a) |
| `potentialFoam`: `cannot find file ".../processor0/0/U"` | `0` was missing when the case was split | Step 15b, "If it does not", item 2 |
| `assert_case_postpro.py --after` reports `5 of 8`, with the wing-surface fields `MISSING` | a known false result of that check on OpenFOAM-12 output | Step 15b, "If you did Step 15a", item 4 |

### `locationInMesh` inside the wing

`snappyHexMeshDict` contains `locationInMesh (-1.0 2.5 1.5);`: a point in the air, ahead of
and outboard of the wing. snappyHexMesh keeps the region that contains this point. If the
point is inside the wing, it keeps the **inside** of the wing and throws the air away.

This was tested by moving the point to `(1.60 0.50 -0.09)`, inside the wing. snappyHexMesh
finished normally ("Finished meshing"). The mesh had 732,446 cells, its only patches were
`symmetry` and `wing`, and its overall bounding box was the size of the wing:
`(1.1767168 0 -0.17539476) (2.2979176 1.8254119 0.013652361)`. Plain `checkMesh` reported
**`Mesh OK.`** The full `checkMesh` of Step 13 reported `Failed 3 mesh checks.`: face tets,
small determinant and concave cells, all three from the "expected" list, and **fewer**
failures than the correct mesh. Only its `*Number of regions: 13` line (one star, so not
counted as a failure) gave it away. A valid-looking mesh of the wrong domain. That is why the
checks in Steps 6 and 13 look at **what** was meshed (cell count, number of regions, patch
list, bounding box), not only at the list of failed checks. The recipe's point is correct for
every wing in this project. Do not change it.

---

## Where to go next

1. [04, sending a case to HPC12](./04_send_to_hpc.md): build the production case fresh from
   the unedited recipe, install the post-processing block, and send it to the cluster.
2. [05, submitting and running](./05_submit_and_run.md): the mesh job, turning the finished
   mesh into a solve case with `hpc/stage_solve_case.sh`, and the solve job on HPC12.
3. [06, trimming to the target lift](./06_trim_to_target_cl.md): how the angle of attack is
   set (by rotating the inflow, never the wing) and how each geometry is trimmed.
4. [07, post-processing and the report](./07_postprocessing_and_report.md): from a finished
   solve to a run card and the report.
