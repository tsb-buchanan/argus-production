# 5. Submit and run on HPC12

> **At a glance**
> 1. **What it gets done:** one PBS job builds the mesh at 192 ranks, `hpc/stage_solve_case.sh` turns the finished mesh into a solve case at 128 ranks, then one or more PBS solve jobs ("legs") run it with OpenFOAM-org 12.
> 2. **Before you start:** OpenFOAM-12 built on HPC12 ([01](./01_setup.md)), a closed wing surface ([02](./02_geometry_from_openvsp.md)), and on the cluster a copy of this repository plus the mesh case, with the post-processing block installed ([04](./04_send_to_hpc.md)).
> 3. **What you have at the end:** a solve case on HPC12 scratch carrying a `.converged` marker and the post-processing output, ready for [06](./06_trim_to_target_cl.md) and [07](./07_postprocessing_and_report.md).
> 4. **How long:** the two delivered mesh jobs whose output was kept took about 2.2 and 2.4 hours (section 5). A 4000-iteration first solve leg needs about 20.2 to 22.1 hours of solver time on 4 typej hosts, or about 11.3 hours on one typen host (section 8), plus time waiting in the queue. The six delivered cases took 17.0 to 34.2 hours of wall time each in total, over 2 or 3 legs (section 10c).

This guide takes a case that is already on HPC12 and runs it: one PBS job builds the
mesh, a script on the login node turns the mesh into a solve case, then one or more PBS
jobs solve it. Everything here is OpenFOAM-org 12 (`foamRun`) and PBS (`qsub`, `qstat`).

Read [01 setup](./01_setup.md) (OpenFOAM-12 on HPC12) and
[04 send to HPC](./04_send_to_hpc.md) (getting the case onto the cluster) first. The
trim procedure that follows a first solve is in
[06 trim to target C_L](./06_trim_to_target_cl.md).

A "leg" in this guide is one PBS solve job. A case that needs more iterations than fit
in one job is continued in a second leg, a third, and so on, each restarting from the
last saved state.

---

## 0. Before you start

Open a terminal on HPC12 from your own machine:

```bash
ssh <your-netid>@hpc12.tudelft.net
```

Replace `<your-netid>` with your own login name, without the angle brackets. This is the
address [04 send to HPC](./04_send_to_hpc.md) uses in Step 2. SSH access and keys for
HPC12 are your own set-up. Every command in this guide, and every HPC12 command in
[06](./06_trim_to_target_cl.md), runs in that terminal.

You need all of the following on HPC12.

1. **OpenFOAM-org 12 built in your home directory** at
   `$HOME/OpenFOAM/OpenFOAM-12`. Both job scripts look there by default. See
   [01 setup](./01_setup.md).
2. **A copy of this repository on the cluster**, made in
   [04 send to HPC](./04_send_to_hpc.md), Step 4. `$REPO` is that copy:

   ```bash
   export REPO=$HOME/argus-production   # or wherever you put the copy on HPC12
   ```

   `$REPO` must be an absolute path, because it goes into the `qsub` lines (section 7b).
   A path starting with `$HOME` is absolute.
3. **A mesh case**, sent as described in [04 send to HPC](./04_send_to_hpc.md): a recipe
   folder from `recipes/` copied **whole** (it must keep its `MANIFEST.sha256`), with the
   wing surface at `constant/triSurface/wing.stl` and the post-processing block installed
   ([04](./04_send_to_hpc.md), Step 5; section 7a of this guide says what the block is).

This guide uses two case directories, laid out as the delivered runs were (the run cards
record `mesh/<name>` and `solve_r2/<name>`; any names work):

```bash
export MESH=/home/scratch/$USER/argus/mesh/<name>       # the mesh case from 04
export CASE=/home/scratch/$USER/argus/solve_r2/<name>   # the solve case, made in step 2
```

`$CASE` must **not** exist yet: step 2 (section 6) creates it. The full paths of `$CASE`
and `$REPO` must contain no space and no comma, because both go into the comma-separated
`qsub -v` list.

For the baseline wing meshed in [04](./04_send_to_hpc.md) as `B_CR`, a worked example is:

```bash
export MESH=/home/scratch/$USER/argus/mesh/B_CR
export CASE=/home/scratch/$USER/argus/solve_r2/B_CR
```

The mesh case and the solve case may share a name: they live in different folders. Never
point either at a delivered case (`mesh/SWB_r2`, `solve_r2/SWB_trim` and the like) on the
account that holds them: the mesh job deletes the mesh in the case it is given, and
`hpc/stage_solve_case.sh` refuses a solve folder that already exists. The delivered run
used `mesh/SWB_r2` and `solve_r2/SWB_trim`.

The delivered `SWB_trim` started its first leg at 1.742369 deg, B's trim angle on the
published mesh (`trim.history` in `data/run_cards/SWB_r2_trim.json`). A solve case gets
its angle in section 6 ("Then set the angle of attack").

**`export` lasts only as long as the terminal it was typed in.** Legs run for days and you
submit each one by hand, so you will log in to HPC12 again. In every new terminal, run the
`export` lines for `REPO`, `MESH` and `CASE` again before any command in this guide. If you
forget, `$REPO` expands to nothing: `$REPO/hpc/stage_solve_case.sh $MESH $CASE` then fails
with `bash: /hpc/stage_solve_case.sh: No such file or directory`. With only `MESH` or
`CASE` missing, the staging script refuses (section 6, "If it does not", item 1).

**Anything in angle brackets in this guide, such as `<name>`, `<N>` or `<job number>`, is a
placeholder.** Replace it, brackets included, before you run the command. Pasted as it
stands, bash reads `<` and `>` as redirections and the command fails with an error.

The mesh job builds the mesh in `$MESH`. Step 2 copies the finished mesh into a new
`$CASE`, and every solve leg runs there. Use absolute paths everywhere: the jobs run on
other machines and do not know which directory you were in.

> **GEOMETRY IS NOT IN THIS REPOSITORY.** No `wing.stl` is shipped. You build it yourself
> from the design study's `.vsp3` files (Step 1 of
> [02 geometry from OpenVSP](./02_geometry_from_openvsp.md) says where each one is) as described in
> [02 geometry from OpenVSP](./02_geometry_from_openvsp.md). The mesh job refuses to
> start without it.

Two things in this repository are **not** part of this route:

1. `hpc/README_HPC.md`, `hpc/submit.sh` and `hpc/push_to_cluster.sh` describe an earlier
   DelftBlue/SLURM setup and are not the HPC12 route.
2. The long comments at the top of `hpc/solve_hpc12.pbs` describe an older set-up: 192
   ranks split 8 x 3 x 8 on a single typen node, and an endTime of 4500. Ignore them. What
   counts are the `#PBS` lines at the very top of the file (queue `fpt-medium`, 4 nodes x
   32 cores = 128) and the case's own files, which the script reads.

---

## 1. The job scripts at a glance

All three are in `hpc/`. The two `.pbs` files are PBS jobs; their `#PBS` header lines set
these defaults. The staging script runs directly at the login prompt.

| | `hpc/mesh_hpc12.pbs` | `hpc/stage_solve_case.sh` | `hpc/solve_hpc12.pbs` |
|---|---|---|---|
| Where it runs | a PBS job on the compute nodes | the login node, typed at the prompt | a PBS job on the compute nodes |
| What it runs | surfaceFeatures, blockMesh, decomposePar, snappyHexMesh (parallel), reconstructPar, checkMesh | copies the finished mesh case into a new solve case and writes its 128-rank `decomposeParDict` | decomposePar and potentialFoam (first leg only), foamRun, then a convergence check |
| Queue (`-q`) | `fpt-large` | none | `fpt-medium` |
| Resources (`-l`) | `nodes=6:ppn=32:typej` | none | `nodes=4:ppn=32:typej` |
| Slots (cores) | 6 x 32 = 192 | none | 4 x 32 = 128 |
| Walltime | `12:00:00` | none | `72:00:00` |
| What you must give it | `ARGUS_CASE` | the mesh case and the new solve case folder | `ARGUS_CASE`, `ARGUS_SCRIPT` and `ARGUS_GATE` (section 7) |

Each job script loads the compiler and MPI modules and sources OpenFOAM-12 itself. You do
not need OpenFOAM loaded in your login shell to submit, and the staging script needs
neither OpenFOAM nor Python.

---

## 2. What the `-l` fields mean

Take the solve script's request as the example:

```
-l nodes=4:ppn=32:typej,walltime=72:00:00
```

1. `nodes=4`: four separate machines (compute nodes).
2. `ppn=32`: 32 cores ("processors per node") on each of them.
3. `typej`: a node property. It restricts the job to nodes of type typej. The node
   types are listed in section 3.
4. `walltime=72:00:00`: the longest the job may run, in hours:minutes:seconds of real
   time. When it is reached PBS kills the job, whatever it is doing.

`nodes x ppn` is the number of **slots** the job gets. PBS writes one line per slot into
a file named by `$PBS_NODEFILE`. Both job scripts count those lines and compare the count
with the case's rank count (section 4).

The other header lines are `-N` (job name), `-q` (queue) and `-j oe` (write error
messages into the same output file as normal output).

**Options on the `qsub` command line override the `#PBS` lines in the script.** The
project used this to send the same solve script to different queues and node types
without editing it.

---

## 3. Queues and node types on HPC12

The project measured the node inventory with `pbsnodes` in September 2026 and recorded
it in the header of `hpc/solve_hpc12.pbs`:

| Node type | Partition | Nodes | Cores per node | Memory per node | Queues |
|---|---|---|---|---|---|
| typei | p8/p9 | 32 | 20 | 125 GB | asm-medium, fpt-small |
| typej | p12 | 36 | 32 | 93 GB | fpt-medium, fpt-large |
| typel | p14 | 8 | 32 | 187 GB | asm-small |
| typem | p15 | 12 | 96 | 251 GB | fpt-large |
| typen | p16 | 6 | 192 | 377 GB | fpt-large |
| typeo | p17 | 2 | 192 | 755 GB | asm-large |

The queue limits below were read on HPC12 in September 2026. To see a queue's current
limits yourself:

```bash
qstat -Qf fpt-medium | grep -E "max_user_queuable|resources_max|resources_default|acl_groups"
```

For fpt-medium this printed `max_user_queuable = 3`, `resources_max.nodes = 4`,
`resources_max.walltime = 72:00:00`, `resources_default.walltime = 72:00:00` and an
`acl_groups` line naming `lr-fpt`.

1. **fpt-medium** gives typej nodes only, at most 4 nodes per job, so **at most 128
   cores**. Each user may have at most 3 jobs in it; a fourth `qsub` is refused.
2. **fpt-large** gives typej, typem and typen nodes and accepts **one job per user,
   queued or running** (`max_user_queuable = 1`). A second `qsub` to fpt-large is refused
   when you submit it; it does not wait in the queue. It is the only queue that can give
   192 slots in one job (6 x typej, or 2 typem nodes at 96 cores, or 1 typen node).
3. **fpt-small** takes at most 1 node per job and up to 10 jobs per user. The project used
   it only for small maintenance work.
4. All three fpt queues allow at most `72:00:00` of walltime, and give `72:00:00` when a
   job asks for none. They admit members of the group `lr-fpt`.
5. A mesh job asks for 6 nodes, so it can only go to fpt-large. **While a mesh job of
   yours is queued or running, a second fpt-large `qsub` (another mesh job, or a solve on
   typen or typem nodes) is refused.** Wait until the first job has left
   `qstat -u $USER` (section 11a), then submit.
6. The asm-* queues (typel, typeo) were not open to the project's account. Your account's
   access may differ.

---

## 4. The rank-count rule: mesh at 192, solve at 128

The **rank count** is the number of MPI processes a parallel OpenFOAM run uses. It is
set by `numberOfSubdomains` in `system/decomposeParDict`, and the case is split into
that many pieces, one per rank.

1. Both job scripts read `numberOfSubdomains` from the case. They never assume it.
2. The solve script **refuses to start when the job has fewer slots than ranks**. With
   more slots than ranks it prints a `NOTE` and leaves the extra cores idle. The
   delivered baseline's first leg got a 192-slot host for its 128 ranks and printed
   `NOTE: 192 slots for 128 ranks; the surplus 64 will idle` (HPC12 job logs).
3. Every recipe ships `numberOfSubdomains 192`, the mesh-stage decomposition. The CR and
   the two cruise recipes use method `scotch`; `L11_wallResolved_WT` uses `hierarchical`
   (8 x 3 x 8). The mesh job's 6 x 32 = 192 slots match that. The solve job's default
   4 x 32 = 128 does not, so a case submitted to the solve script with the recipe's own
   dictionary stops at once with:

   ```
   REFUSED: allocation gives 128 slots but the case decomposes to 192; fix the -l line
   ```

4. **The delivered runs were meshed at 192 ranks and solved at 128 ranks.** All six
   Condition CR run cards in `data/run_cards/` record the mesh built at 192 ranks with
   method `scotch` (`mesh.mesh_generation`, read from each mesh case) and the solve at 128
   ranks with method `hierarchical` (`decomposition`, and `hpc.solve_decomposition_method`).
   Read back on HPC12, all six solve cases carry `hierarchicalCoeffs` with `n (8 2 8)`,
   `delta 0.001`, `order xyz`. The reason for 128 is the queues in section 3: at 128 ranks
   a solve fits fpt-medium, where you can have three jobs, instead of competing for the
   single fpt-large place.
5. **The mesh must be built with the shipped 192-rank dictionary.** The mesh job checks
   every recipe file against `MANIFEST.sha256`, including `system/decomposeParDict`, and
   refuses a case where any of them was changed. So the switch to 128 happens **after**
   the mesh job, in a separate solve case made by `hpc/stage_solve_case.sh` (step 2,
   section 6). The mesh case keeps its shipped dictionary.
6. How the switch works: the mesh job finishes by collecting the whole mesh into
   `$MESH/constant/polyMesh`. The staging script copies that collected mesh into `$CASE`,
   without any `processor*` directories, and writes a new `system/decomposeParDict` there:
   128 subdomains, `hierarchical`, `n (8 2 8)`. The first solve leg then runs
   `decomposePar -force`, which splits the mesh into as many pieces as `decomposeParDict`
   asks for. Re-splitting a finished mesh does not change the mesh.

A note on the mesh-stage split method. The CR and cruise recipes split the mesh with
`scotch` at 192 ranks, which is how the six delivered wake-refined meshes were built
(their mesh cases carry 192 subdomains, method `scotch`, and their snappyHexMesh logs
report `nProcs : 192`; HPC12 job logs). The method is part of the mesh recipe, not a
scheduling setting: parallel snappyHexMesh builds a slightly different mesh on a
different split. Measured on the baseline wing at 4 refinement levels below production,
`scotch` against `hierarchical` changed the wetted area by 0.005% (1.28019 against
1.28026 m2) and the total cell count by 0.60% (243,379 against 241,911)
([recipes/README.md](../recipes/README.md)). So do not edit the recipe's
`decomposeParDict`; the mesh job would refuse it anyway (section 5, Gate A).

---

## 5. Step 1: submit the mesh job

From the login node:

```bash
cd $MESH
qsub -N argus-mesh-<name> -v ARGUS_CASE=$MESH $REPO/hpc/mesh_hpc12.pbs
```

`qsub` prints the job identifier (a number followed by the server name). Keep the
number. Job names are up to you; the project used `argus-mesh-<case>` for mesh jobs.

Where to look while it runs and after it ends:

1. `qstat -u $USER` shows whether the job is waiting (`Q`) or running (`R`); once it has
   gone from the list, the job has ended (section 11a).
2. While it runs, its output so far is in `$HOME/<job number>.hpc12.hpc.OU` (section 11b).
3. When it ends, PBS writes the job's output file, `argus-mesh-<name>.o<job number>`, into
   `$MESH`, the folder you ran `qsub` from (section 11b). "The job's output file" in the
   rest of this section means this file, or the live copy while the job runs.

Environment variables the mesh script reads:

| Variable | Required | Default in the script | What to set |
|---|---|---|---|
| `ARGUS_CASE` | yes | none; the job stops with `ARGUS_CASE not set` | absolute path of the mesh case, `$MESH` |
| `ARGUS_OF_BASHRC` | no | `$HOME/OpenFOAM/OpenFOAM-12/etc/bashrc` | only if your OpenFOAM-12 is somewhere else |

**About the mesh job's resources.** The command above uses the script's header defaults
(fpt-large, 6 x 32 typej, 12 hours). The project submitted its late-cruise meshes exactly
this way. The delivered wake-refined Condition CR meshes ran in fpt-large with 192 slots
across 1, 3 or 4 hosts, and the exact `-l` lines were not recorded (HPC12 job logs). The
two of those mesh jobs whose output files were kept took about 2.2 and 2.4 hours. Whether
the defaults are enough for a 99-million-cell wake-refined CR mesh is **NOT VERIFIED**.
If `snappyHexMesh` is killed for lack of memory, the script's own header suggests two
typem nodes instead (2 x 96 = 192 cores, 502 GB in total) rather than assuming the mesh is
at fault:

```bash
cd $MESH
qsub -N argus-mesh-<name> -q fpt-large -l nodes=2:ppn=96:typem,walltime=12:00:00 \
     -v ARGUS_CASE=$MESH $REPO/hpc/mesh_hpc12.pbs
```

**NOT VERIFIED**: this request does not appear in any recorded submission.

### What the mesh job checks before it meshes

1. **OpenFOAM-12 is present.** Otherwise `REFUSED: OpenFOAM-12 not found at ...` or
   `REFUSED: snappyHexMesh not on PATH after sourcing ...`. If you see either, OpenFOAM-12
   is not built where the script looks, or `ARGUS_OF_BASHRC` points at a different
   OpenFOAM version. Go back to [01 setup](./01_setup.md).
2. **Slots against ranks**, as in section 4.
3. **Gate A, the recipe is intact.** Every file listed in `MANIFEST.sha256` must be
   byte-for-byte what the recipe shipped. Two files are allowed to differ,
   `system/controlDict` and `0.orig/U`, because they carry the angle of attack (and
   `system/controlDict` also carries the post-processing block's include line). For those
   two the job instead checks that the inflow direction, `dragDir` and `liftDir` all
   imply the same angle. `system/argusPostPro` is not in the manifest, so the block never
   trips this gate. For the CR recipe with the block installed, a pass looks like:

   ```
       alpha 1.299999 deg, |U| 34.0000 m/s, axis lag -2.04e-08 deg
     GATE A: 19 manifest entries = 17 intact + 0 altered + 0 missing + 2 alpha-bearing
     recipe intact
   ```

   The first line shows your angle. For a cruise recipe the Gate A line reads
   `22 manifest entries = 20 intact + 0 altered + 0 missing + 2 alpha-bearing`.

   A changed file is named on an `ALTERED` line (for example
   `    ALTERED system/decomposeParDict`); it is never `system/argusPostPro`. The job then
   ends with two lines:

   ```
   REFUSED: the recipe dicts in this case do not match the manifest they
            shipped with. Meshing would build an unrecorded recipe.
   ```

   Start again from a fresh copy of the recipe ([04](./04_send_to_hpc.md), Step 5).

   **The WT recipe as shipped is refused here.** A fresh copy of
   `recipes/L11_wallResolved_WT` stops at Gate A with the lines below (each `ALPHA` line
   ends with a bracketed reference, left out here):

   ```
       alpha 2.102800 deg, |U| 40.8000 m/s, axis lag -3.66e-06 deg
       ALPHA dragDir implies alpha 2.102797, Uinf implies 2.102800
       ALPHA liftDir implies alpha 2.102797, Uinf implies 2.102800
     GATE A: 19 manifest entries = 17 intact + 0 altered + 0 missing + 2 alpha-bearing
   REFUSED: the recipe dicts in this case do not match the manifest they
            shipped with. Meshing would build an unrecorded recipe.
   ```

   Its `liftDir` and `dragDir` are written to 6 significant figures, so they lag the
   inflow by more than the gate's 1e-6 deg. WT is not one of the three delivered
   conditions. If you mesh it, set its angle in the mesh case before you submit (best on
   your own machine, before sending it), as in [06](./06_trim_to_target_cl.md), section
   6.4.2, with `UMAG=40.8`, even if you keep its 2.1028 deg. Tested on a fresh copy with
   06's `sed` route and `ALPHA=2.1028`: the gate then printed `axis lag -2.14e-12 deg` and
   `recipe intact`. The CR and cruise recipes pass as shipped.
4. **Gate B, the surface is closed and crosses the symmetry plane.** Every STL in
   `constant/triSurface` must have no open edges and no over-connected edges (an edge
   shared by more than two triangles), and must reach at least 10 mm below `y = 0`. For
   the baseline meshing surface (sha256 `20bd6acd...`, the one bound to the baseline in
   `registry/candidates.yaml`) a pass looks like:

   ```
     GATE B constant/triSurface/wing.stl: tris=192588 open=0 over=0 y_min=-20.000mm OK
     surface closed
   ```

   **Do not work around this gate.** A mesh built from an open surface fills the inside
   of the wing, and `checkMesh` passes it, because it is a perfectly valid mesh of the
   wrong domain. See [02 geometry from OpenVSP](./02_geometry_from_openvsp.md).
   If the STL is simply missing you get `GATE B: no STL in constant/triSurface` followed
   by the open-shell `REFUSED` message; the first line is the real reason.

   **The surface must be an ASCII (text) STL.** The gate reads only the text `vertex`
   lines, so a binary STL gives
   `GATE B COULD NOT RUN: ValueError: min() arg is an empty sequence` followed by
   `REFUSED: gate B could not run. NOT a verdict about the geometry.` The registered
   surfaces and the ones [02](./02_geometry_from_openvsp.md) builds are ASCII.
   `head -c 5 $MESH/constant/triSurface/wing.stl; echo` must print `solid` (02, Section 7),
   and for the registered baseline `grep "^solid" $MESH/constant/triSurface/wing.stl`
   prints `solid wing_closed` ([03](./03_build_a_case.md), Step 5). If yours is binary,
   rebuild it with 02, then resubmit.

**The post-processing block does not change the mesh.** The tutorial mesh of
[03](./03_build_a_case.md) was built twice on 8 cores, once with the block and once
without, using the mesh job's order of programs. Both gave 1,193,347 cells and the same
four failed checks, the logs differed only in times and process numbers, and no meshing
program wrote a `postProcessing` folder. A 192-rank production mesh job with the block
installed has run on HPC12 (2026-09-24, an early-cruise mesh of W, 73.7 million cells, the same
four failed checks), but with no mesh of the same surface without the block to compare, "does
not change the mesh" is **NOT VERIFIED** at production size.

### What a finished mesh job prints

At the end of the job's output file there is a block read out of the logs:

```
--- VERDICTS, read from the logs ---
  checkMesh : <Mesh OK, or Failed N mesh checks>
  cells     : <cell count>
  non-orth  : <maximum non-orthogonality line>
  aspect    : <maximum aspect ratio line>
  layers ...
  extrusion : Extruding <n> out of <m> faces (<p>%)
```

What the six delivered Condition CR meshes looked like, so you can compare (run cards and
HPC12 job logs):

1. **Cells:** 99,111,506 to 99,208,470.
2. **checkMesh:** every one ends `Failed 4 mesh checks.`, so do not expect `Mesh OK` from
   this recipe. The four are `***Max skewness` (8.86 to 17.35, on 4 to 29 faces),
   `***Error in face tets` (3162 to 3572 faces), `***Cells with small determinant`
   (13.00 to 13.05 million cells) and `***Concave cells` (2.10 to 2.12 million cells).
   The maximum non-orthogonality is 70.96 to 71.97 deg, with one-star warnings only; the
   run cards carry it as `mesh.max_non_orthogonality`, so you can compare your
   `non-orth` line with it. The maximum aspect ratio is 112.81923 (113.77236 for W). To
   see which checks your mesh failed, list the lines checkMesh marks with `***`:

   ```bash
   grep -F '***' $MESH/log/checkMesh.log
   ```

   Each failed check prints one such line, for example
   ` ***Max skewness = 5.2082764, 3 highly skew faces detected which may impair the quality of the results`
   (from the tutorial mesh in 03). Four failures are expected, the ones listed above and
   in 03, Step 13. Any other `***` line means the mesh is broken and you should not
   solve it. That includes `***Number of non-orthogonality errors`,
   `***Zero or negative cell volume detected` and `***Open cells found`.
3. **Layer extrusion:** 97.40% to 97.51% of wing faces. A value near 0% means the
   boundary layers failed.
4. **Layers:** the output prints a requested row and an achieved row. The achieved
   thickness can only be at or below the requested one. On the delivered meshes the
   achieved row read about 9.8 layers of 11, 95.2 to 95.3% of the requested thickness.

Check that OpenFOAM-12 built the mesh (a mesh from another version is not comparable).
The job's own output file prints, near the top, a line
`OpenFOAM: <path>/snappyHexMesh  WM_PROJECT_VERSION=12`. The snappyHexMesh log says the
same in its header:

```bash
grep -m1 -E '^Build +: 12' $MESH/log/snappyHexMesh.log
```

It must print a line starting `Build  : 12`. On HPC12 the project's OpenFOAM-12 build
prints `Build  : 12` alone (read from the delivered mesh logs); a packaged install
elsewhere prints `Build  : 12-<code>`. If the command prints nothing, read the first
`Build` line of the log yourself: any other number means another OpenFOAM version built
the mesh.

The full logs are in `$MESH/log/`: `surfaceFeatures.log`, `blockMesh.log`,
`decomposePar.log`, `snappyHexMesh.log`, `reconstructParMesh.log` (despite its name it
holds the output of `reconstructPar -constant`, which collects the mesh under
OpenFOAM-12) and `checkMesh.log`. Each new mesh job in the same directory overwrites
them, so copy `log/` somewhere else first if you want to keep a failed run's logs.

**When the mesh job stops early**, the last lines of its output file say why:

| You see | What to do |
|---|---|
| `REFUSED: gate A could not run` or `REFUSED: gate B could not run` | the check itself failed to run, which says nothing about the case; read the `COULD NOT RUN` line above it (for gate B, usually a binary STL) |
| `surfaceFeatures failed` or `blockMesh failed` | read the matching file in `$MESH/log/` |
| `decomposePar failed` | read `$MESH/log/decomposePar.log`. If it says `You are trying to use scotch but do not have the scotch library loaded`, your OpenFOAM-12 build has no `scotch`, which the CR and cruise recipes need. Build OpenFOAM-12 with it ([01 setup](./01_setup.md), Part B); do not edit the recipe (section 4) |
| `REFUSED: decomposePar produced N dirs, expected M` | read `$MESH/log/decomposePar.log` |
| `snappyHexMesh failed; see log/snappyHexMesh.log` | the job prints the last 25 lines of that log below this line. If a process was killed for lack of memory, see the resources note above |
| the output simply stops after `snappyHexMesh on 192 ranks` | PBS killed the job, normally at its walltime; the script is killed with it and prints nothing more. See the resources note above |
| `REFUSED: the reconstructed mesh has no wall patch` | the mesh was not collected; read `$MESH/log/reconstructParMesh.log` and `$MESH/log/snappyHexMesh.log` |

Resubmitting to the same mesh case is safe once the cause is fixed: after its checks pass,
the job deletes the previous attempt's `processor*`, `constant/polyMesh` and
`constant/extendedFeatureEdgeMesh` and builds again.

**Never point the mesh job at a solve case.** Once its checks pass, it deletes
`processor*` and `constant/polyMesh`. A solve case made in step 2 is refused by Gate A
before that happens (`ALTERED system/decomposeParDict`, because the staging script
replaced that file); do not undo that protection.

---

## 6. Step 2: make the solve case at 128 ranks with `hpc/stage_solve_case.sh`

Only after the mesh job has finished and its verdicts look right. Run this at the HPC12
login prompt. It needs only bash, coreutils, sed and grep; it does not need OpenFOAM or
Python.

```bash
$REPO/hpc/stage_solve_case.sh $MESH $CASE
```

Usage (`--help` prints this line, with your own path to the script):

```
usage: $REPO/hpc/stage_solve_case.sh <mesh-case-dir> <solve-case-dir> [--ranks N]
```

What each argument and option does:

1. `<mesh-case-dir>` (`$MESH`): a case meshed by `hpc/mesh_hpc12.pbs` from a recipe
   folder. It must hold the reconstructed `constant/polyMesh`, `0.orig/`, `system/`, the
   other `constant/` files (including `constant/triSurface/`), `MANIFEST.sha256` and
   `RECIPE_FILES.json`, plus `CONDITION.json` for a cruise recipe. The script only reads
   this folder.
2. `<solve-case-dir>` (`$CASE`): the new solve case. It must not exist yet. The script
   never overwrites and never deletes anything. Its full path, and the path of your copy
   of this repository, must contain no space and no comma, because all of them go into
   the `qsub -v` list.
3. `--ranks N` (or `--ranks=N`): the number of ranks for the solve. The default is 128.
   `ARGUS_RANKS=N` in the environment also sets it. If you give both, the flag wins and
   the script says the environment value was ignored.
   1. 128: `hierarchical`, n (8 2 8), delta 0.001, order xyz. All six delivered Condition
      CR solve cases on HPC12 carry exactly this.
   2. 192: `hierarchical`, n (8 3 8), delta 0.001, order xyz. These are the case
      definitions of the published 192-rank solves.
   3. Any other N: `scotch`, with a printed note that this departs from the delivered
      configuration.
4. `-h` or `--help`: prints the usage line and stops.

What it does, in order:

1. It checks the inputs: the mesh case, the new solve case, and every path the `qsub` line
   will carry. If a check fails, it prints `REFUSED: ...` and `Nothing was created.`, and
   exits with status 2.
2. It copies `0.orig/`, `system/` (including `system/argusPostPro` if present),
   `constant/`, `MANIFEST.sha256`, `RECIPE_FILES.json` and `CONDITION.json` into
   `$CASE.incomplete_<UTC time>`. It never copies `0/`, `processor*` or `log/`, and it
   prints every top-level entry as copied or not copied.
3. It writes `system/decomposeParDict` for the rank count. It then reads everything back
   from the files: each copied file is compared with its original, the rank count is read
   the way `solve_hpc12.pbs` reads it, and the cold-start conditions are checked.
4. It writes `SOLVE_SETUP.json`, renames the folder to `$CASE`, and prints the next steps
   with ready-to-run `qsub` lines.

Exit status: 0 means the case is staged. 2 means refused (see "If it does not" below for
the two refusals that can leave a folder behind). 1 means a failure after copying began
(`FAILED: ...`). In that case the partial folder `$CASE.incomplete_<UTC time>` stays in
place for you to read and then remove yourself.

**What you should see.** This is the tutorial mesh of baseline B at Condition CR
([03](./03_build_a_case.md)), built from the current CR recipe with the post-processing
block installed. The script prints full paths where this example shows `$MESH`, `$CASE`
and `$REPO`. A production mesh has 99 million cells, so its copy takes longer and the
case is far larger.

```
=== stage_solve_case.sh 2026-09-23T23:23:08Z ===
  ranks: 128 (the default, as the r2 solves; neither --ranks nor ARGUS_RANKS was given)
  rule:  128 ranks: hierarchical (8 2 8), the r2 solves' configuration
  mesh case:  $MESH
  solve case: $CASE
  mesh:       nCells:1193347, 12 files in constant/polyMesh, 215811108 bytes to copy in all
  top level of the mesh case:
    copied        0.orig
    copied        MANIFEST.sha256
    copied        RECIPE_FILES.json
    copied        constant/ (polyMesh verified file by file; the rest by sha256)
    NOT copied    log: the mesh job's logs stay with the mesh case
    copied        system
    NOT copied    processor* (8 directories): the mesh job's decomposition;
                  the first solve leg runs decomposePar at 128 ranks
    14 entries = 5 copied + 9 excluded
  staging in $CASE.incomplete_20260923T232308Z
  copying constant/polyMesh (12 files; the slow part)
  mesh copied in 0 s
  verifying
    decomposeParDict: 128 subdomains, hierarchical (8 2 8), product asserted
    mesh: 12 files = 11 identical by sha256 + 1 identical by size (over 64 MiB)
    other files: 25 = 24 identical by sha256 + 1 decomposeParDict replaced
    staged file set equals the source set; no symlinks
    cold start: 0.orig holds 5 fields (U k nut omega p); no 0/, processor*, log/ or .converged
    nCells 1193347 (owner header)
  STAGED: $CASE (206M)

NEXT STEPS
1. Post-processing block.
   system/controlDict already includes argusPostPro.
   Check it here before you submit. The check needs only the Python standard
   library, and it must end "---- 8 of 8 requirements satisfied ----":
       python3 $REPO/scripts/assert_case_postpro.py --before $CASE
2. Angle of attack. The case carries the angle the mesh case was given:
       Uinf            (33.991249 0 0.771369);
       liftDir         ( -0.022687323 0 0.999742610 );
       dragDir         ( 0.999742610 0 0.022687323 );
   To run another angle, set it as in docs/06_trim_to_target_cl.md, section 6.4.2
   (a cold start), and check it with section 6.4.4.
3. Iterations. system/controlDict: endTime 4000, writeInterval 500. The r2
   solves ran their cold leg to endTime 4000 with writeInterval 500. This script
   changes neither.
4. Submit (128 ranks) from the HPC12 login node, from inside the solve case: PBS
   copies the job's output file, argus-B_CR.o<job id>, into the folder qsub runs from when
   the job ends. While it runs, HPC12 keeps that text in $HOME/<job id>.hpc12.hpc.OU.
       cd $CASE
       qsub -N argus-B_CR -q fpt-medium -l nodes=4:ppn=32:typej,walltime=72:00:00 \
            -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py \
            $REPO/hpc/solve_hpc12.pbs
   Faster, if you have no other fpt-large job (fpt-large holds one job per user),
   in place of that qsub line:
       qsub -N argus-B_CR -q fpt-large -l nodes=1:ppn=128:typen,walltime=72:00:00 \
            -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py \
            $REPO/hpc/solve_hpc12.pbs
   The r2 legs ran at about 18.2 to 19.9 s per iteration on 4 typej hosts and
   about 10.1 s per iteration on one typen host (solver ClockTime).
   ARGUS_GATE makes solve_hpc12.pbs run the step 1 check at the start of every leg,
   before the write probe, decomposePar and foamRun, and stop a leg that fails it
   with "REFUSED: post-processing gate refused". Without ARGUS_GATE the job looks
   for $HOME/argus_hpc12/assert_case_postpro.py, which does not exist on a new
   account, and runs the leg unchecked. Once the job has started:
       grep -E "requirements satisfied|gate NOT run|gate refused" $HOME/<job id>.hpc12.hpc.OU
   should print "---- 8 of 8 requirements satisfied ----". After the job ends, grep
   argus-B_CR.o<job id> in the solve case instead.
   Add ARGUS_OF_BASHRC=<path> to the -v list if OpenFOAM-12 is not in
   $HOME/OpenFOAM/OpenFOAM-12.
5. Later legs. Submit every leg yourself, from the same folder, with the same -v
   list, ARGUS_GATE included. A job's own qsub for its next leg fails on HPC12
   (could not connect to trqauthd), and would not pass ARGUS_GATE on anyway.
   Add ARGUS_ENDTIME where docs/05_submit_and_run.md, section 10a, says so. It
   must be a multiple of writeInterval (500), or solve_hpc12.pbs refuses the leg.
```

"r2" in the output is the project's name for the six delivered wake-refined Condition CR
cases. The job name comes from the solve case folder's name (`argus-<name>`). Sections 7a
and 7b go through the next steps one by one. The script writes `<job id>` for the job
number that `qsub` prints; the rest of this guide calls it `<job number>`.

With `--ranks` other than 128, item 4 prints `-q <queue> -l <resources giving at least N
slots>` for you to fill in, because the delivered runs used no such line (section 7b has
the one the project used for 192).

If the mesh case was built without the post-processing block, item 1 starts instead with
the lines below, followed by the install commands and the fallback of section 7a:

```
1. Post-processing block.
   WARNING: system/controlDict does NOT include argusPostPro: the mesh case
   was built without the block.
```

If `system/controlDict` has the include line but `system/argusPostPro` is missing, the
first lines read `WARNING: system/controlDict includes argusPostPro, but the file` /
`system/argusPostPro is missing.`, and the fallback says
`The include line is already in system/controlDict.` A mention of argusPostPro inside a
comment does not count as an include.

**If it does not.** Every `REFUSED` exits with status 2. Almost every one also prints
`Nothing was created.` and creates nothing. The exceptions are two refusals that can only
come when the script makes the staging folder, just before the copy:
`REFUSED: cannot create the parent of '...'` and `REFUSED: cannot create '...'`. They do
not print `Nothing was created.`, and they may leave behind parent folders of `$CASE` that
the script had just made. Remove those yourself if you do not want them. The refusals you
are most likely to meet:

1. `REFUSED: expected a mesh case and a solve case, got 1 path(s). usage: ...` (or
   `got 0 path(s)`): `$MESH` or `$CASE` is not set in this terminal, so the script was
   given fewer than two paths. Run the `export` lines of section 0 again.
2. `REFUSED: mesh case '...' does not exist or is not a directory`: `$MESH` points at the
   wrong folder. Check the path, and set it again as in section 0.
3. `REFUSED: solve case '...' already exists.`: choose a new name, or move the old case
   aside yourself.
4. `REFUSED: the solve case's full path '...' contains a space or a comma.`: the check
   runs on the full path, so a plain name such as `B_CR` is refused if you are inside a
   folder whose path has a space or a comma. Choose a solve case folder whose full path
   has neither.
5. `REFUSED: '.../hpc/solve_hpc12.pbs' contains a space or a comma.` (or
   `'.../scripts/assert_case_postpro.py'`): your copy of the repository sits in such a
   folder. Put it in a folder whose path has no space and no comma, and run the script
   from there.
6. `REFUSED: mesh case '...' has N processor* directories but no reconstructed constant/polyMesh`:
   the mesh lives only in `processor*`. The mesh job did not collect it; read its output
   file.
7. `REFUSED: constant/polyMesh/boundary has no wall patch.`: that is the background
   blockMesh, so the mesh job did not finish.
8. `REFUSED: rank count '...' (--ranks) is not a positive integer` (or `is below 1`):
   give a whole number of ranks.
9. `REFUSED: the copy needs N KiB and '...' has M KiB free`: make room on scratch first.
   On HPC12 the `constant/polyMesh` of the delivered mesh case SWB_r2 occupies
   12,162,360 KiB (the `ls -l` total), far more than the tutorial mesh above.
10. `FAILED: mesh files did not copy identically: ...` (status 1): the copy was damaged.
    Read the message, remove the `.incomplete_` folder yourself, and run the script again.

What has been run where:

1. On HPC12 (bash 4.2.46, GNU coreutils 8.22), read-only: `bash -n` passes. Against the
   delivered mesh case SWB_r2, every mesh input check passed, and the script refused at
   the free-space and permission check because the destination was deliberately
   unwritable. A destination with a space was refused. Nothing was created.
2. The whole script, copying, verifying, writing the dictionary and printing the next
   steps, has run on HPC12 once: on 2026-09-24 it staged a 192-rank early-cruise mesh case
   of W (73.7 million cells) at the default 128 ranks. It ended `STAGED: <solve case> (16G)`,
   wrote `SOLVE_SETUP.json`, and the solve leg submitted with its printed `qsub` line passed
   the post-processing gate (`8 of 8 requirements satisfied`). It had run on a local Linux
   machine before that (bash 5.0.17, coreutils 8.30).

### Then set the angle of attack

The CR recipe ships at 1.30 deg. That is a starting value, not a trim. Before you submit,
give the solve case the angle you mean to run.
[06](./06_trim_to_target_cl.md), sections 6.5.1 and 6.5.2, says how to choose it.

1. The staging script copies the angle the mesh case carries, and prints it under
   `NEXT STEPS`, item 2. If you set the angle in the mesh case on your own machine before
   sending it ([04](./04_send_to_hpc.md), Step 5), the solve case already carries it.
2. Otherwise, set it now in `$CASE` with [06](./06_trim_to_target_cl.md), section 6.4.2.
   The staged case has no `0/` and no `processor*`, so it starts cold, as 6.4.2 needs.
3. `foamDictionary` runs at the HPC12 login prompt once OpenFOAM-12 is loaded in the same
   shell. Type these lines one at a time; do not send the `source` line through a pipe,
   which runs it in a separate shell and throws the settings away:

   ```bash
   module load devtoolset/11
   module load mpi/openmpi-4.1.2
   source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc
   foamDictionary -help | grep -E "Using|writePrecision"
   ```

   **What you should see:** the `source` line prints two harmless lines,
   `dirname: missing operand` and `Try 'dirname --help' for more information.` (every
   job log on HPC12 shows them too). Then:

   ```
     -writePrecision <label>
   Using: OpenFOAM-12 (see https://openfoam.org)
   ```

   Tested on HPC12. Before these lines `foamDictionary` is not on the PATH. If your
   OpenFOAM-12 is somewhere else, source that one's `etc/bashrc` instead.

   > **`foamDictionary` reads every file path relative to the folder you are in, even a
   > path that starts with `/`.** `foamDictionary -entry numberOfSubdomains -value $CASE/system/decomposeParDict`
   > fails with `file "<current folder>/<the path you gave>" does not exist`. Always
   > `cd $CASE` first and give paths such as `system/controlDict`.

   For example, these two lines read the staged split back from inside the solve case:

   ```bash
   cd $CASE
   foamDictionary -entry numberOfSubdomains -value system/decomposeParDict
   foamDictionary -entry method -value system/decomposeParDict
   ```

   They print `128` and `hierarchical` (read this way on HPC12 from a delivered solve
   case, and locally from a freshly staged one).
4. Either way, run the check in 06, section 6.4.4, in `$CASE`. It needs only `python3`.
   All three angles must agree and must be the angle you intended.

---

## 7. Step 3: submit the solve job

### 7a. Check the post-processing block first

The post-processing block is one file, `system/argusPostPro`, plus one line in
`system/controlDict` that includes it. It makes the solver write what
[07 post-processing and report](./07_postprocessing_and_report.md) needs:

1. the pressure and viscous parts of the forces;
2. the wing-surface fields p, cp, y+ and the wall-shear-stress vector;
3. four Trefftz planes behind the wing;
4. the residuals.

You install it once, on your own machine, before the case is sent
([04](./04_send_to_hpc.md), Step 5, with `scripts/install_postpro.py`). The mesh job
ignores it, the staging script copies it into the solve case with `system/`, and the
solve job checks it before every leg (section 7b). Without the block the solve still
runs, and the convergence check at the end of each leg still works, because it needs
only `forceCoeffs1`. But none of the four is written, and the missing data cannot be
produced afterwards.

Check the solve case at the HPC12 login prompt before you submit:

```bash
python3 $REPO/scripts/assert_case_postpro.py --before $CASE
```

It needs only the Python standard library. The login node's `python3` (3.6.8) runs it: it
gave 8 of 8 on each of the six delivered Condition CR solve cases.

**What you should see** (the first line names your solve case):

```
=== /home/scratch/<your-netid>/argus/solve_r2/<name> (before) ===
  ok      forces_split
  ok      wing_surface_p
  ok      wing_surface_tau_VECTOR
  ok      wing_surface_yplus
  ok      surface_format_carries_topology
  ok      compressible_cp_frame
  ok      trefftz_planes
  ok      residuals
  ---- 8 of 8 requirements satisfied ----
```

`compressible_cp_frame` is always `ok` under OpenFOAM-12: it only recognises the
OpenFOAM-org 7 compressible solver names, so it checks nothing there
([07](./07_postprocessing_and_report.md), section 4). The `forceCoeffs1` force
coefficients that every recipe writes are not one of the eight.

**If it does not** end `8 of 8`, some lines read `MISSING` instead of `ok`, each followed
by two lines of explanation, and the last line is `---- N of 8 requirements satisfied ----`
with N below 8. `1 of 8`, with only `compressible_cp_frame` ok, means the mesh case was
built without the block, and the staging script said
`system/controlDict does NOT include argusPostPro`. Install the block on the solve case
before the first leg, at the HPC12 login prompt:

```bash
python3 $REPO/scripts/install_postpro.py $CASE
python3 $REPO/scripts/assert_case_postpro.py --before $CASE
```

The install prints a short summary of what it derived from the case. On a solve case
staged from the tutorial mesh of baseline B at Condition CR, built without the block, it
read (the first line is your case folder's name):

```
  <name>
     solver        foamRun -> rho mode 'rhoInf', rhoInf 1, pInf 0
     wall patch    wing           (from constant/polyMesh/boundary)
     CofR          1.34 0 0       (from existing forces object)
     turbulence    kOmegaSST      -> residual fields: U p k omega
     |U| for cp    34.0000        (from 0.orig/U (via $Uinf))
     TE x 2.3363, chord 1.1597 (from ASCII STL wing.stl)
     Trefftz x     2.3942 2.6262 2.9161 3.2061
     controlDict now includes argusPostPro
```

These are the same values as the install on your own machine
([04](./04_send_to_hpc.md), Step 5): the `TE x` and `Trefftz x` values are in metres, in
the frame of `wing.stl`, and `|U| for cp` is the Condition CR speed in m/s. Running the
script a second time changes nothing and prints `controlDict already includes argusPostPro`.
The `--before` check afterwards tells you whether the install worked; it must end
`8 of 8`.

What was tested on HPC12:

1. The login node has numpy 1.12.1, which the install script needs.
2. The script's derivations ran there with `--dry-run` under Python 3.6.8, on a delivered
   mesh case, and printed the same values as on your own machine. There the wall patch
   `wing` is read from `constant/polyMesh/boundary`.
3. The writing step has not been run on HPC12: **NOT VERIFIED**. On a local machine, the
   same install on a staged solve case wrote a byte-identical block and gave 8 of 8.

If the install fails on HPC12:

1. `ModuleNotFoundError: No module named 'numpy'`: the `python3` you ran has no numpy.
   Nothing was written. Use the login node's own `python3`, which has it.
2. Anything else: run the install on your own machine's copy of the case
   ([04](./04_send_to_hpc.md), Step 5). Copy only its `system/argusPostPro` into
   `$CASE/system/` (from your machine, for example with `scp` as in
   [04](./04_send_to_hpc.md), Step 6b, item 4, with the solve case's `system/` folder on
   HPC12 as the destination). In `$CASE/system/controlDict`, add the line
   `    #include "argusPostPro"` directly under the `{` that follows `functions`. Then
   repeat the check (**NOT VERIFIED** on HPC12).

### 7b. Submit, with the gate switched on

`hpc/solve_hpc12.pbs` runs the same check at the start of each leg, but only if it finds
the checker file. It looks at the path in `ARGUS_GATE`. Without `ARGUS_GATE`, it looks at
`$HOME/argus_hpc12/assert_case_postpro.py`, which does not exist on a new account. So put
the shipped copy in the `-v` list of **every** leg you submit:

```bash
ls $REPO/scripts/assert_case_postpro.py
cd $CASE
qsub -N argus-<name> -q fpt-medium -l nodes=4:ppn=32:typej,walltime=72:00:00 \
     -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py \
     $REPO/hpc/solve_hpc12.pbs
```

This is the line the staging script printed, with your paths already filled in.

1. The `ls` line must print the file's path.
2. `$REPO` and `$CASE` must be set in the login shell and must be absolute paths. They are
   expanded when you press Enter, and the job looks for the files only after it has moved
   into the case folder.
3. `-v` takes a comma-separated list with no spaces.
4. The `-q` and `-l` values are the same as the script's header. They are written out so
   you can see what you are asking for.
5. `qsub` prints the job identifier. Keep the number. PBS writes the job's output file,
   `argus-<name>.o<job id>`, into the folder you ran `qsub` from, here `$CASE`, when the
   job ends (section 11b).
6. On HPC12 you submit every leg by hand (section 10a). For later legs, add
   `ARGUS_ENDTIME` as in section 10a, and keep `ARGUS_GATE`.

**Faster, if you have no other fpt-large job.** One typen host ran the delivered legs at
about 10.1 s per iteration, against 18.2 to 19.9 s on 4 typej hosts (section 8). Use this
in place of the `qsub` line above:

```bash
cd $CASE
qsub -N argus-<name> -q fpt-large -l nodes=1:ppn=128:typen,walltime=72:00:00 \
     -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py \
     $REPO/hpc/solve_hpc12.pbs
```

fpt-large holds one job per user, so while a mesh job of yours is queued or running this
`qsub` is refused (section 3). The project's tooling also submitted 128-rank solves as
`-q fpt-large -l nodes=2:ppn=64:typem,walltime=72:00:00`; its speed on the wake-refined
meshes was not recorded.

**A 192-rank solve** needs 192 slots in fpt-large. Stage the solve case with
`--ranks 192` (section 6), which writes `hierarchical` (8 3 8):

```bash
$REPO/hpc/stage_solve_case.sh $MESH $CASE --ranks 192
```

The project ran a 192-rank solve on HPC12 with one typen node:

```bash
cd $CASE
qsub -N argus-<name> -q fpt-large -l nodes=1:ppn=192:typen,walltime=72:00:00 \
     -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py \
     $REPO/hpc/solve_hpc12.pbs
```

That request ran 4500 iterations on an earlier mesh, not on a wake-refined one. The
solve script's header comment also suggests 6 typej nodes
(`-l nodes=6:ppn=32:typej,walltime=72:00:00`). **NOT VERIFIED**: no recorded solve used
that request.

**Read the gate's verdict** once the job has started. While the job runs, HPC12 keeps its
output in `$HOME/<job number>.hpc12.hpc.OU` (section 11b); after it ends, the same text is
in `argus-<name>.o<job number>` in `$CASE`:

```bash
grep -E "requirements satisfied|gate NOT run|gate refused" $HOME/<job number>.hpc12.hpc.OU
grep -E "requirements satisfied|gate NOT run|gate refused" $CASE/argus-<name>.o<job number>
```

Use the first line while the job runs and the second once it has ended.

**What you should see:**

```
    ---- 8 of 8 requirements satisfied ----
```

The check runs after the job has started and loaded OpenFOAM, before the write probe,
`decomposePar`, `potentialFoam` and `foamRun`, using the compute node's `python3` (3.6).
The delivered legs ran a `--before` check there, from the project's own copy of the
checker. A refusal still costs the queue wait, which is why section 7a runs the check
before `qsub`.

**If it does not:**

1. `NOTE: <path> absent, post-processing gate NOT run`: nothing exists at the path in
   `ARGUS_GATE`. It may be unset, misspelled or relative, or `$REPO` was not set. The leg
   carries on without the check. Fix the `-v` list for the next leg.
2. `---- N of 8 requirements satisfied ----` (N below 8), then
   `REFUSED: post-processing gate refused`: the job stopped before it touched the case.
   It leaves only an empty `log/leg<N>` (section 10a, point 5). Install the block as in
   section 7a, then resubmit.

How this was tested: the gate lines of the job script were run unchanged on staged cases.

1. With the block: 8 of 8, and the job carries on.
2. Without the block: `1 of 8`, then `REFUSED: post-processing gate refused`, exit status 2.
3. With `ARGUS_GATE` unset or misspelled: the `NOTE` line, and the job carries on.

It has been: on 2026-09-24 a solve leg on HPC12 (an early-cruise case of W, job output file of leg 1) received `ARGUS_GATE` through `-v` (it is in the job's `qstat -f` variable list) and ran the check before its first iteration, ending `8 of 8 requirements satisfied`.

In `log/leg<N>/foamRun.log`, a case with the block shows `wallShearStress argusWallShearStress:`,
then `processing wall patches:` and `wing`, then two `Reading surface description:`
blocks. The first lists `wingSurface`; the second lists `trefftz_x1` to `trefftz_x4`.
The block adds no warnings.

### 7c. Environment variables the solve script reads

| Variable | Required | Default in the script | What to set |
|---|---|---|---|
| `ARGUS_CASE` | yes | none; the job stops with `ARGUS_CASE not set` | absolute path of the case |
| `ARGUS_SCRIPT` | yes, in practice | none | absolute path of `solve_hpc12.pbs` on the cluster, `$REPO/hpc/solve_hpc12.pbs`. Without it a leg that needs a follow-on ends with `ARGUS_SCRIPT: leg N needs to requeue but ARGUS_SCRIPT is unset; refusing to guess` and does not print the endTime you need (section 10). |
| `ARGUS_GATE` | yes, on every leg | `$HOME/argus_hpc12/assert_case_postpro.py`, which does not exist on a new account | `$REPO/scripts/assert_case_postpro.py` (section 7b) |
| `ARGUS_ENDTIME` | no | the case's own `endTime`: 4000 in the CR recipe, 4500 in the WT recipe, 10000 in the two cruise recipes | only to override. It **must be a multiple of the case's `writeInterval`** (500 in CR, 150 in WT, 2000 in cruise) or the leg stops with `REFUSED: endTime ... is not a multiple of writeInterval ...`. For example 4250 is refused for the CR recipe, and 5500 is accepted. |
| `ARGUS_MAXLEG` | no | 4 | the leg number after which the script stops asking for another leg |
| `ARGUS_LEG` | no | worked out from the `log/leg<N>` directories already in the case | **do not set it by hand.** A leg number that is already used overwrites that leg's solver log. |
| `ARGUS_OF_BASHRC` | no | `$HOME/OpenFOAM/OpenFOAM-12/etc/bashrc` | only if your OpenFOAM-12 is somewhere else |

**About `ARGUS_GATE`.** With the shipped checker in the `-v` list, every leg runs the
check of section 7a first and a case that fails it stops at once with
`REFUSED: post-processing gate refused`. Without `ARGUS_GATE`, the job looks at the
default path. On a new account nothing is there, so every leg prints (with your own home
directory written out)
`NOTE: /home/<your-netid>/argus_hpc12/assert_case_postpro.py absent, post-processing gate NOT run`
and carries on unchecked. Keep the gate on for every leg.

Every variable you add to the `-v` list of the first leg (`ARGUS_GATE`,
`ARGUS_OF_BASHRC`, `ARGUS_MAXLEG`) must be added again to every later leg you submit by
hand (section 10a). A later leg without `ARGUS_GATE` runs without the gate. A later leg
without `ARGUS_OF_BASHRC` looks for OpenFOAM-12 in the default place.

---

## 8. Walltime: always `72:00:00` for a solve

Measured on the delivered wake-refined Condition CR meshes (99.1 to 99.2 million cells)
at 128 ranks, from the solver's ClockTime (HPC12 job logs):

1. **4 typej hosts** (`-q fpt-medium -l nodes=4:ppn=32:typej`): 18.2 to 19.9 s per
   iteration on the four 4000-iteration first legs that ran there. One later leg, F's
   second (2000 iterations), averaged about 24 s per iteration.
2. **One typen host** (`-q fpt-large -l nodes=1:ppn=128:typen`): about 10.1 s per
   iteration, so a 4000-iteration leg takes about 11.3 h.

| Iterations in the leg | Solver time at 18.2 s/it | at 19.9 s/it |
|---|---|---|
| 4000 (the CR recipe's first leg, as the delivered first legs) | 72,800 s (20.2 h) | 79,600 s (22.1 h) |
| 2000 (a typical later leg) | 36,400 s (10.1 h) | 39,800 s (11.1 h) |
| 1500 (one automatic raise of endTime, section 10a) | 27,300 s (7.6 h) | 29,850 s (8.3 h) |

On top of that the first leg spends time in `decomposePar` (8 to 13 minutes on these
meshes) and `potentialFoam`. The whole first-leg jobs on 4 typej hosts took 73,823 s to
80,684 s (20.5 to 22.4 h) of wall time (run cards' notes). So a 24-hour request is
marginal on typej: it would have left between about 1.6 and 3.5 hours to spare. 72 hours
holds any leg the delivered runs needed.

1. Always submit solves with `walltime=72:00:00`. It is also the most the fpt queues
   allow (section 3).
2. On a busy machine a 72-hour request on fpt-medium can wait longer to start than a
   short one, because the scheduler has fewer gaps it fits into. In September 2026 the
   scheduler estimated start times about 1.6 to 2.5 days away for three 72-hour
   fpt-medium jobs. Two 72-hour fpt-large jobs, submitted on different occasions (the
   queue holds one job per user, section 3), each started almost as soon as it was
   submitted. A long wait in the `Q` state is not a fault.
3. You cannot rescue a job whose walltime is too short: `qalter` cannot raise the
   walltime of a job that is already running (**NOT VERIFIED**: recorded by the project,
   not reproduced in any script or job log in this repository). Get it right at
   submission.
4. If a leg is killed at its walltime anyway, see section 9: it restarts from its last
   checkpoint.

---

## 9. What one solve leg does, and how a restart works

### 9a. The sequence inside one leg

1. **Refuses a finished case.** If `$CASE/.converged` exists the job prints the marker
   and stops (section 10c).
2. **Works out its leg number** from the `log/leg<N>` directories already present, and
   creates `log/leg<N>/`.
3. **Checks slots against ranks** (section 4).
4. **Loads OpenFOAM-12** and checks that `foamRun` is on the path.
5. **Runs the post-processing gate** (section 7b): the checker at `ARGUS_GATE`, or at the
   default path. If no file is there it prints a `NOTE` and carries on.
6. **Write test:** writes, syncs and reads back 8 small files in the case. Fewer than 7
   good ones stops the job, because a filesystem that cannot write will lose a checkpoint
   silently.
7. **Decides cold start or restart by reading the case** (9b).
8. **Resolves endTime** from the case, or from `ARGUS_ENDTIME`, checks it is a multiple of
   `writeInterval`, writes it into `system/controlDict`, and checks that `controlDict`
   says `startFrom latestTime`.
9. **Cold start only:** copies `0.orig` to `0`, runs `decomposePar -force`, checks it made
   exactly one `processor*` directory per rank, then runs `potentialFoam` to give the
   solver a smooth starting flow field.
10. **Runs `foamRun`** on all ranks. A launch that produces no iterations is tried up to
    3 times in total, 60 s apart. A `FOAM FATAL` error is not retried, because it will
    fail the same way again.
11. **Reads the outcome from disk**: the newest time directory, and the convergence check
    (section 10b).
12. **Finalises** (writes `.converged`) or **asks for another leg** (section 10a).

A healthy first leg's output file looks like this. The lines are the script's own
messages; the numbers are what the CR recipe at 128 ranks gives, which are the settings
the delivered first legs ran with. Times, node names and paths differ:

```
[hh:mm:ss] === <name> leg 1/4 on <node> ===
[hh:mm:ss] NRANK=128 (read from system/decomposeParDict)
[hh:mm:ss] allocation: 128 slots across 4 host(s)
dirname: missing operand
Try 'dirname --help' for more information.
[hh:mm:ss] foamRun: <your OpenFOAM-12>/platforms/<options>/bin/foamRun
[hh:mm:ss] WM_OPTIONS=<options>

  === /home/scratch/<your-netid>/argus/solve_r2/<name> (before) ===
    ok      forces_split
    ok      wing_surface_p
    ok      wing_surface_tau_VECTOR
    ok      wing_surface_yplus
    ok      surface_format_carries_topology
    ok      compressible_cp_frame
    ok      trefftz_planes
    ok      residuals
    ---- 8 of 8 requirements satisfied ----
[hh:mm:ss] write probe: 8 ok, 0 bad of 8
[hh:mm:ss] required fields (from 0.orig): U k nut omega p
[hh:mm:ss] COLD leg: no decomposed field to continue from
[hh:mm:ss] endTime: 4000 (from the case)
[hh:mm:ss] checkpoints: 8 writes at every 500 iterations
[hh:mm:ss] 0.orig -> 0 (5 fields)
[hh:mm:ss] all 5 required fields present in 0/
[hh:mm:ss] decomposePar
[hh:mm:ss] decomposed onto 128 ranks
[hh:mm:ss] potentialFoam
[hh:mm:ss] foamRun returned 0
[hh:mm:ss] state: latestTime=4000 reached=1
  CONVERGENCE GATE over 200 samples: MEAN drift Cd ... ct (<=0.05), Cl ... ct (<=1.0), Cl SPAN ... ct (<=1.0) -> CONVERGED
  REPORT THESE: Cd ...  Cl ...  (window means)
[hh:mm:ss] -> FINALISE: reached endTime AND converged
[hh:mm:ss] wrote .converged marker
```

The `dirname` pair comes from sourcing OpenFOAM-12 on HPC12 and is harmless. If
`potentialFoam returned non-zero (continuing; it is an initialiser)` appears, read
`log/leg1/potentialFoam.log`. The job carries on, but a wall-resolved case started
without a smooth initial field can crash in its first few iterations.

### 9b. Restart from checkpoint

1. **Checkpoints.** The solver saves the full flow state every `writeInterval`
   iterations (500 in the CR recipe, 150 in WT, 2000 in the cruise recipes) into
   `processor*/<iteration>/`. The CR and WT recipes set `purgeWrite 3` (the cruise
   recipes 2), so only the newest three are kept. `purgeWrite` does not touch
   `postProcessing/`: the post-processing block's files from every saved iteration stay
   (section 12).
2. **`startFrom latestTime`** in `system/controlDict` makes `foamRun` continue from the
   newest checkpoint instead of from zero. The script refuses a case without it.
3. **How the script classifies a leg.** It never trusts a leg number; it reads the case:
   1. No `processor0` directory: **cold start**.
   2. `processor*` directories with saved iterations: **restart** from the newest
      checkpoint that is complete on **every** rank. It lists checkpoints across all
      ranks, not only rank 0, because a checkpoint missing from some ranks stops
      `foamRun` with `Start time is not the same for all processors`.
   3. A checkpoint that is incomplete on some ranks (for example one being written when
      the job was killed) is deleted from every rank that holds it, but only when an older
      complete one exists. If no complete checkpoint exists the job stops with
      `REFUSED: no complete checkpoint exists; every time directory is torn. Nothing deleted.`
      That case needs manual recovery; the tools for it are not in this repository.
   4. A restart prints lines like:

      ```
      [hh:mm:ss]   checkpoint t=3000: whole=128 torn=0 of 128
      [hh:mm:ss]   latestTime agreement: all 128 ranks resolve to t=3000
      [hh:mm:ss] RESTART leg: resuming from t=3000
      ```

4. **What a walltime kill costs**: the iterations since the last complete checkpoint, at
   most one `writeInterval` (500 iterations in the CR recipe, about 2.5 to 2.8 hours at
   18.2 to 19.9 s/it), plus the wait in the queue for the next leg.
5. **Never delete the `processor*` directories of a case you mean to continue.** Without
   them the next leg is treated as a cold start and `potentialFoam` overwrites the flow
   field. Rebuilding them from a collected time is not covered by this repository.
6. **Do not change `endTime` in a case that is running.** The job reads it once, at the
   start of the leg. If you lower it while the solver runs, the solver stops at the new
   value but the job then believes the run was cut short. See
   [06 trim to target C_L](./06_trim_to_target_cl.md). To stop earlier than planned,
   follow 06, section 6.7: stop the job (`qdel`, section 11a) and resubmit with
   `ARGUS_ENDTIME` set to the first multiple of `writeInterval` at or above the iteration
   it reached. To run longer, let the leg finish and submit the next one with a higher
   `ARGUS_ENDTIME`.

---

## 10. Legs, resubmitting, and the `.converged` marker

### 10a. How a leg ends, and why you resubmit by hand on HPC12

Every leg ends in one of four ways. The output file says which:

| Last lines of the output file | Meaning | What you do |
|---|---|---|
| `-> FINALISE: reached endTime AND converged` then `wrote .converged marker` | done | nothing more to solve. Go on to [06](./06_trim_to_target_cl.md) or [07](./07_postprocessing_and_report.md). |
| `-> reached endTime but NOT converged. Raising endTime and continuing.` | ran to endTime but the forces are still moving | submit the next leg with `ARGUS_ENDTIME` set to the value on the `requeueing as leg ... with endTime ...` line. If that line is missing, use the endTime printed near the start of the leg plus 1500. For a cruise case see point 2. |
| `-> did not reach endTime (walltime). Continuing at the same endTime.` | foamRun stopped before endTime while the job was still alive, for example after an error (point 3 below) or because endTime was changed during the run (section 9b, point 6) | read the end of `$CASE/log/leg<N>/foamRun.log` to find out why. After an error: fix the case, then submit the next leg with no `ARGUS_ENDTIME`. If endTime was changed during the run: submit the next leg with `ARGUS_ENDTIME` set to the first multiple of `writeInterval` at or above the iteration the leg reached (section 9b, point 6). Without it the leg takes the edited endTime from the case, and refuses it with `REFUSED: endTime ... is not a multiple of writeInterval ...` unless it happens to be a multiple |
| none of these: the output simply stops, usually after the `potentialFoam` line (first leg) or the `checkpoints:` line (later legs) | PBS killed the job, normally at its walltime. The script is killed with it, so it prints nothing more. | submit the next leg with no `ARGUS_ENDTIME`; it restarts from the last complete checkpoint (section 9b) |

The script then tries to submit the next leg itself with `qsub`, from inside the job.
**On HPC12 that does not work**: compute nodes cannot reach the PBS server, and the
`qsub` inside the job fails with
`qsub: cannot connect to server (null) (errno=15137) could not connect to trqauthd`
(HPC12 job logs). That job-made `qsub` would not pass `ARGUS_GATE` on anyway. So **you
submit every leg yourself**, with the same `qsub` line as in section 7b plus
`ARGUS_ENDTIME` where the table says so. Example for the second row, after a CR leg that
ended unconverged at endTime 4000:

```bash
cd $CASE
qsub -N argus-<name> -q fpt-medium -l nodes=4:ppn=32:typej,walltime=72:00:00 \
     -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py,ARGUS_ENDTIME=5500 \
     $REPO/hpc/solve_hpc12.pbs
```

If the first leg's `-v` list also carried `ARGUS_OF_BASHRC` or `ARGUS_MAXLEG`, add them to
this line too (section 7c).

Points to get right:

1. **Why `ARGUS_ENDTIME` matters in the second row.** `system/controlDict` still holds the
   old endTime. Resubmitting without it gives a leg that starts at endTime, runs 0
   iterations (`0 iterations: already at endTime ...`), fails the same check and ends
   the same way. For the CR recipe the raise by 1500 gives 5500, then 7000, all multiples
   of its `writeInterval` of 500.
2. **Cruise recipes.** Their `writeInterval` is 2000, so the old endTime + 1500 is not a
   multiple of it: after a first leg to 10000, `ARGUS_ENDTIME=11500` is refused with
   `REFUSED: endTime 11500 is not a multiple of writeInterval 2000; the final state would never be written`.
   Use the next multiple of 2000 instead, 12000 in that example.
3. **The third row is also what a crash looks like.** After a `FOAM FATAL` error the output
   shows `---- DETERMINISTIC FAULT, NOT RETRYING ----` with the error text, and then
   still says `did not reach endTime (walltime)`. Read the error and fix the case before
   you resubmit; a resubmission fails the same way. The solver log is
   `$CASE/log/leg<N>/foamRun.log`.
4. **Never run two jobs on one case at the same time.** Two solvers writing the same
   `processor*` directories destroy the case, and the solve script does not check for
   this. Before resubmitting, make sure no job for this case is still queued or running
   (section 11a). In particular, if a job named `argus-<case directory name>-L<N>` has
   appeared, the automatic follow-on did get submitted; do not add another.
5. **Leg numbers only count up.** A job that stops at a `REFUSED` check still leaves an
   empty `log/leg<N>` directory behind, so the next job becomes leg N+1. The script stops
   asking for more legs once the leg number reaches `ARGUS_MAXLEG` (default 4); it then
   prints `MAXLEG 4 reached, not requeueing`. You can still submit further legs by hand.
   An empty directory left by a refused job can be removed; `rmdir` only removes a
   directory that is empty, so it cannot delete a real leg's logs:

   ```bash
   rmdir $CASE/log/leg<N>
   ```

### 10b. The convergence check

At the end of every leg the script reads `postProcessing/forceCoeffs1/*/forceCoeffs.dat`,
sorted by iteration, and takes the last 200 samples. `forceCoeffs1` writes a sample every
5 iterations in the CR recipe, so that window is the last 1000 iterations. A drag
count (ct) is 0.0001 in a force coefficient. The check passes when all three hold:

1. the mean of `Cd` over the second half of the window differs from the mean over the
   first half by at most 0.05 ct;
2. the same for `Cl`, at most 1.0 ct;
3. `Cl` at the end of the window differs from `Cl` at its start by at most 1.0 ct.

It checks how much the **mean** moves, not the size of the wiggle: a converged steady
solve keeps oscillating slightly around a fixed mean forever. The window means it prints
on the `REPORT THESE:` line are the values the project reports. A leg that begins with a
deliberate change of angle (a trim leg) needs the window to fall after that change; see
[06](./06_trim_to_target_cl.md).

### 10c. The `.converged` marker

`$CASE/.converged` is written **only** when a leg has reached endTime **and** passed the
check. It holds the evidence: the case name, when and where it was finalised, the leg,
and the check's own output. Example, from the baseline B run card
(`data/run_cards/SWB_r2_trim.json`, Condition CR, wake-refined mesh, U 34.0 m/s;
coefficients on the half-model `Aref 0.620462 m2`, which is half of the DSO-basis
`Sref 1.24092 m2`, and `lRef 0.393957 m`):

```
case SWB_trim
finalised 2026-09-21T18:16:17Z on n12-142, leg 3
  CONVERGENCE GATE over 200 samples: MEAN drift Cd 0.0003 ct (<=0.05), Cl 0.0399 ct (<=1.0), Cl SPAN 0.0777 ct (<=1.0) -> CONVERGED
  REPORT THESE: Cd 0.0188031  Cl 0.4282682  (window means)
```

1. **A case carrying `.converged` is never restarted.** Submitting it again prints
   `REFUSING to restart <name>: .converged present`, then the marker, then
   `Remove the marker deliberately to extend this case.` and ends without touching
   anything.
2. **To extend a finished case on purpose**, rename the marker rather than deleting it,
   so the evidence is kept, then submit with a higher `ARGUS_ENDTIME`:

   ```bash
   mv $CASE/.converged $CASE/.converged.leg<N>
   ```

3. **Only the file counts.** A case that looks converged part-way through a leg (for
   example by reading its force file yourself) is still running, and its window means can
   still move. Wait for `.converged`.

What to expect overall: the six delivered Condition CR cases each ran 6000 iterations in
2 or 3 legs, for 61,125 s to 123,274 s (17.0 to 34.2 hours) of wall time in total at
128 ranks (`cost.wall_clock_s` in `data/run_cards/`). Each ran a first leg to endTime
4000 with `writeInterval` 500, the CR recipe's own settings, then trim legs up to 6000
(see [06](./06_trim_to_target_cl.md)). In iterations per leg: B 4000 + 1000 + 1000; C, F, H and W 4000 + 2000; M
4000 + 1500 + 500, after its endTime was changed during a leg (section 9b, point 6).
Their legs mixed 4 typej hosts and the roughly twice as fast single typen host
(section 8), which is the main reason the shortest total is about half the longest.

---

## 11. Monitoring a job

### 11a. `qstat`, and the column that lies

List your jobs:

```bash
qstat -u $USER
```

The `S` column is the state: `Q` waiting, `R` running. Two traps:

1. **The column headed "Elap Time" is not elapsed time.** On HPC12 it shows CPU time
   summed over every core of the job, so on a 128-core solve it races ahead of the clock
   by up to about the number of cores. Measured example: a 128-rank solve showed
   `72:12:06` there, next to a requested time of `72:00:00`, when it had been running for
   `00:58:06`; the displayed figure was about 75 times the real running time. It looks
   like a job about to be killed and is not.
2. **The job name column is cut to 16 characters.** Longer names are truncated and
   different jobs can look identical.

For one job, get the real figures from the full listing:

```bash
qstat -f <job number> | grep -E 'Job_Name|job_state|resources_used.walltime'
```

`resources_used.walltime` is the real time the job has been running. Compare it with
the 72:00:00 you asked for.

To cancel a job of yours (the project's own tooling used this on HPC12):

```bash
qdel <job number>
```

### 11b. Where the output goes

1. **The job's own output file** (the `[hh:mm:ss]` lines above) is written as
   `<job name>.o<job number>` into the directory you ran `qsub` from, which is why the
   commands above `cd $MESH` or `cd $CASE` first. It appears there when the job ends.
2. **While the job runs**, HPC12 keeps the same text in `$HOME/<job number>.hpc12.hpc.OU`,
   which you can read:

   ```bash
   tail -n 20 $HOME/<job number>.hpc12.hpc.OU
   ```

   The project's own monitoring on HPC12 read running mesh and solve jobs this way. The
   file name carries the job number, not the job name. The file is gone from `$HOME` once
   the job has ended; the text is then in the output file of point 1. Once the solver
   starts it goes quiet on purpose, because solver output goes to its own log.
3. **The solver's logs**, per leg: `$CASE/log/leg<N>/decomposePar.log`,
   `potentialFoam.log` (cold legs only) and `foamRun.log`. A launch retry keeps the
   earlier attempt as `foamRun.attempt<k>.log`.
4. **The force history**: `$CASE/postProcessing/forceCoeffs1/<start iteration>/forceCoeffs.dat`.
   Each restart starts a new directory named after the iteration it resumed from.

### 11c. What a healthy solve looks like

Where the solver is now:

```bash
grep "^Time = " $CASE/log/leg<N>/foamRun.log | tail -n 1
```

OpenFOAM-12 prints the iteration with a trailing `s`, for example `Time = 1500s`. It
should keep increasing. Each iteration in `foamRun.log` looks like:

```
Time = 1500s

smoothSolver:  Solving for Ux, Initial residual = ..., Final residual = ..., No Iterations ...
smoothSolver:  Solving for Uy, ...
smoothSolver:  Solving for Uz, ...
GAMG:  Solving for p, Initial residual = ..., Final residual = ..., No Iterations ...
GAMG:  Solving for p, ...
time step continuity errors : sum local = ..., global = ..., cumulative = ...
smoothSolver:  Solving for omega, ...
smoothSolver:  Solving for k, ...
ExecutionTime = ... s  ClockTime = ... s
```

Seconds per iteration: compare `ClockTime` on the last two iterations. On the
wake-refined CR meshes at 128 ranks expect about 18 to 20 on 4 typej hosts (about 24 was
seen once), and about 10 on one typen host (section 8).

```bash
grep "^ExecutionTime" $CASE/log/leg<N>/foamRun.log | tail -n 2
```

The latest force coefficients (columns: `Time Cm Cd Cl Cl(f) Cl(r)`, normalised on the
half-model `Aref 0.620462 m2` and `lRef 0.393957 m` set in the recipe's `forceCoeffs1`):

```bash
tail -n 2 $CASE/postProcessing/forceCoeffs1/*/forceCoeffs.dat
```

OpenFOAM writes its log in bursts, so a few minutes without new lines is normal.

### 11d. When to worry

| You see | Meaning | What to do |
|---|---|---|
| `REFUSED: allocation gives 128 slots but the case decomposes to 192; fix the -l line` | the case still carries the mesh recipe's 192-rank `decomposeParDict` (it was not made with the staging script), or you staged it with `--ranks 192` and asked for 128 slots | stage a solve case as in section 6, or ask for 192 slots (section 7b), then resubmit |
| `REFUSED: OpenFOAM-12 not found at ...` or `REFUSED: foamRun not on PATH after sourcing ...` | OpenFOAM-12 is not where the script looks, or `ARGUS_OF_BASHRC` points at a different OpenFOAM version | [01 setup](./01_setup.md) |
| `NOTE: ... absent, post-processing gate NOT run` | `ARGUS_GATE` is missing from the `-v` list, misspelled or relative | the leg runs unchecked; fix the `-v` list for the next leg (section 7b) |
| `REFUSED: post-processing gate refused` | the gate ran and the case fails it | install the block as in section 7a, check 8 of 8, then resubmit |
| `REFUSED: filesystem failed N of 8 write probes` | the scratch filesystem is misbehaving | resubmit later |
| `REFUSED: 0.orig is empty or missing` | the solve case was not made with the staging script, which refuses such a mesh case | stage a new solve case as in section 6 |
| `REFUSED: endTime ... is not a multiple of writeInterval ...` | bad `ARGUS_ENDTIME` | pick a multiple (section 7c) |
| `REFUSED: controlDict is not startFrom latestTime` | the `startFrom` line in `system/controlDict` was changed | set it back to `startFrom       latestTime;` as in the recipe, then resubmit |
| the last lines of `decomposePar.log`, then `REFUSED: decomposePar failed` | the first leg could not split the mesh | read `$CASE/log/leg<N>/decomposePar.log`. If it says `You are trying to use scotch but do not have the scotch library loaded`, the case was staged with a rank count other than 128 or 192, and your OpenFOAM-12 build has no `scotch`: use the fallback below this table |
| `REFUSED: ARGUS_SCRIPT=... does not exist` | the path in the `-v` list is wrong | check `ls $REPO/hpc/solve_hpc12.pbs` and resubmit |
| `attempt N produced NO iterations ... retrying` | the MPI launch failed; the job retries by itself | only act if all 3 attempts fail |
| `---- DETERMINISTIC FAULT, NOT RETRYING ----` | a `FOAM FATAL` error; the text follows | fix the case, then resubmit (section 10a, point 3) |
| `DIAGNOSIS: ranks disagree on latestTime` | the checkpoints differ between ranks | manual recovery; not covered by this repository |
| `REFUSED: the convergence gate could not run` | no `forceCoeffs1` output | check that `controlDict` still has `forceCoeffs1` on patch `wing` |
| residuals growing, `nan`, or a floating point exception in `foamRun.log` | the solution is diverging | stop and inspect the mesh and set-up; do not simply resubmit |

**Fallback when your OpenFOAM-12 has no `scotch`.** This only concerns a solve case staged
with `--ranks` other than 128 or 192: those two already use `hierarchical`. The simplest
fix is to stage a new solve case at 128 or 192 (section 6). To keep your rank count,
switch the case to `hierarchical` and add a complete `hierarchicalCoeffs` block whose three
numbers multiply to that count. The example is for 64 ranks, 8 x 2 x 4:

```bash
cd $CASE
sed -i 's/^method .*/method          hierarchical;/' system/decomposeParDict
cat >> system/decomposeParDict <<'EOF'

hierarchicalCoeffs
{
    n               (8 2 4);
    delta           0.001;
    order           xyz;
}
EOF
grep -E '^(numberOfSubdomains|method)|^ +n ' system/decomposeParDict
```

The `grep` must print these three lines, with your own count and split:

```
numberOfSubdomains 64;
method          hierarchical;
    n               (8 2 4);
```

Tested on a solve case staged locally with `--ranks 64` from the tutorial mesh: under
OpenFOAM-12, `decomposePar -force` then printed `Selecting decomposer hierarchical` and
made 64 `processor*` directories. Changing only the `method` line is not enough:
`decomposePar` then stops with `keyword n is undefined` (tested the same way).

Then resubmit as in section 7b. The failed leg made no `processor*` time directories to
continue from, so the new leg starts cold again and runs `decomposePar` with the new
split.

---

## 12. After the run

1. **Scratch is not permanent.** HPC12 deletes files under scratch that are older than
   50 days, with no recovery. Copy back what you need (see
   [04 send to HPC](./04_send_to_hpc.md)) well before that.
2. **The post-processing output is large.** `purgeWrite` removes none of it. For the
   delivered B solve, `postProcessing/argusWingSurface` holds 13 files (7,698 MiB) and
   `postProcessing/argusTrefftz` 329 MiB. Allow for that on scratch and when you copy
   results back.
3. **The `--after` check reports 5 of 8 on a correctly run case.**
   `python3 $REPO/scripts/assert_case_postpro.py --after $CASE` marks
   `wing_surface_p`, `wing_surface_tau_VECTOR` and `wing_surface_yplus` as `MISSING`
   because it looks for the field names only near the start of each wing-surface file,
   and OpenFOAM-12 writes them after the points and faces. It also reads every
   wing-surface file in full before it gives up, which takes a while. Check those fields
   directly instead:

   ```bash
   cd $CASE
   T=$(ls postProcessing/argusWingSurface | sort -n | tail -1)
   grep -a -m4 -E "^(p|cp|yPlus|wallShearStress) [0-9]+ [0-9]+ float" postProcessing/argusWingSurface/$T/wingSurface.vtk
   ls postProcessing/argusTrefftz/$T
   ```

   For the delivered B solve (where `T` is 6000) this printed:

   ```
   p 1 5276772 float
   cp 1 5276772 float
   yPlus 1 5276772 float
   wallShearStress 3 5276772 float
   trefftz_x1.vtk  trefftz_x2.vtk  trefftz_x3.vtk  trefftz_x4.vtk
   ```

   The `3` after `wallShearStress` shows that the vector was written. The same check on
   a short local smoke run of the tutorial case printed the same four lines with 82922
   wing faces, while `--after` reported `5 of 8`. `--time` appears in the checker's usage
   text but does not exist: it stops with `error: unrecognized arguments`. What the block
   writes, file by file, is in
   [07 post-processing and report](./07_postprocessing_and_report.md), section 3 ("What the
   post-processing block adds"); the other `--after` results are in 07, section 4.6.
4. A first solve is usually followed by a trim to the target lift coefficient:
   [06 trim to target C_L](./06_trim_to_target_cl.md).
5. Run cards, figures and the report: [07 post-processing and report](./07_postprocessing_and_report.md).

---

## Quick reference

```bash
# on HPC12, with REPO, MESH and CASE set (section 0)

# 1. mesh at 192 ranks (fpt-large, 6 x 32 typej, 12 h)
cd $MESH
qsub -N argus-mesh-<name> -v ARGUS_CASE=$MESH $REPO/hpc/mesh_hpc12.pbs

# 2. after the mesh job: make the solve case at 128 ranks, hierarchical (8 2 8)
$REPO/hpc/stage_solve_case.sh $MESH $CASE
python3 $REPO/scripts/assert_case_postpro.py --before $CASE    # must end 8 of 8
# 2b. set the angle of attack (06, section 6.4.2) and check it (06, section 6.4.4)

# 3. solve at 128 ranks (fpt-medium, 4 x 32 typej, always 72 h), gate on
cd $CASE
qsub -N argus-<name> -q fpt-medium -l nodes=4:ppn=32:typej,walltime=72:00:00 \
     -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py \
     $REPO/hpc/solve_hpc12.pbs

# watch
qstat -u $USER
qstat -f <job number> | grep -E 'Job_Name|job_state|resources_used.walltime'
grep -E "requirements satisfied|gate NOT run|gate refused" $HOME/<job number>.hpc12.hpc.OU
grep "^Time = " $CASE/log/leg<N>/foamRun.log | tail -n 1

# finished?
cat $CASE/.converged
```

Back to the [README](../README.md).
