# 4. Sending a case to HPC12 and bringing results back

> **At a glance**
> 1. **What this does:** copies this repository to HPC12, builds one mesh case fresh from a recipe on your machine (the set-up files, the wing surface and the post-processing block), sends it, and later brings the results back.
> 2. **Before you start:** an HPC12 account, the set-up from [01](./01_setup.md) (including a Python with numpy), the wing surface you built with [02](./02_geometry_from_openvsp.md), and the case-building walkthrough in [03](./03_build_a_case.md).
> 3. **At the end:** a mesh case in `/home/scratch/<your-netid>/argus/mesh/<name>` that passes the post-processing check (8 of 8), ready for [05](./05_submit_and_run.md); after the run, the results and the mesh logs in `$HOME/argus_results/`.
> 4. **How long:** this guide gives no step timings. The pull (Step 8) usually comes days after the push and can be several GiB per case, and scratch deletes files older than 50 days (Step 9).

This guide moves a case from your machine to the HPC12 cluster, and later brings the
results back. It comes after [03, building a case](./03_build_a_case.md) and before
[05, submitting and running](./05_submit_and_run.md). Everything here is OpenFOAM-org 12
on HPC12 with PBS.

> **GEOMETRY IS NOT IN THIS REPOSITORY.**
> There is no STL and no `.vsp3` file anywhere in it. The designs come as `.vsp3` files
> ([02, geometry from OpenVSP](./02_geometry_from_openvsp.md), Step 1, says where each one is).
> You build the wing surface you mesh from the `.vsp3`
> yourself, on your own machine, with [02, geometry from OpenVSP](./02_geometry_from_openvsp.md).
> Step 5c lists the checksum each finished surface must have. Nothing in this guide works
> until you have that surface on your own disk.

Every command below runs **in a bash terminal on your own machine** (Linux or WSL, set up
as in [01](./01_setup.md)). The ones that start with `ssh` run their quoted part on the
cluster.

---

## Before you start

1. An HPC12 account. Your login name is your netID; this guide writes it as
   `<your-netid>`.
2. A clone of this repository. If you have not set it up yet, see
   [01, setup](./01_setup.md), which introduces:

   ```bash
   export REPO=$HOME/argus-production   # or wherever you cloned it
   ```

3. `ssh` and `rsync` on your machine. Check both:

   ```bash
   ssh -V
   rsync --version | head -1
   ```

   *Expected:* a line starting `OpenSSH_` and a line starting `rsync  version`. If
   either says `command not found`, install it with your system's package manager first.
   HPC12 has `rsync` too (version 3.1.2), which the copies in this guide need at both ends.
4. The wing surface, built with [02](./02_geometry_from_openvsp.md) (see the box above).
5. A terminal whose `python3` has numpy, for Step 5e. [01](./01_setup.md), section A5,
   calls it your Python terminal.

---

## The route, and why it is this one

**You send the case set-up and the surface. The mesh is built on HPC12.** You never send
a mesh, and you never send `processor*` directories.

Three reasons, each read from files in this repository or from HPC12 job logs and
listings:

1. **The rank count and the split method are part of the mesh recipe.** A *rank* is one
   of the parallel processes a run is split into. `snappyHexMesh` does not produce the
   same mesh at different numbers of ranks (header of `hpc/mesh_hpc12.pbs`), and a
   different split method changes the mesh too. The recipe's own `decomposeParDict`
   comment gives the measured effect of the method on the baseline wing, meshed 4
   refinement levels below production: wetted area 1.28019 m2 with `scotch` against
   1.28026 m2 with `hierarchical` (0.005%), and 0.60% in cell count. The Condition CR and
   cruise recipes fix both in `system/decomposeParDict`: 192 ranks, method `scotch`. That
   is what all six delivered wake-refined Condition CR meshes carry (HPC12 job logs: every
   one of their `snappyHexMesh` logs reports `nProcs : 192`; also `mesh.mesh_generation`
   in `data/run_cards/*.json`). The mesh job refuses any edit to that file. So a
   Condition CR mesh built by this route is split the same way as the six delivered
   Condition CR meshes. Whether it reproduces their cell counts exactly has not been
   tested. The delivered early- and late-cruise meshes are different: they were built
   with the older recipe, before its wake refinement was added, split `hierarchical`
   (8 3 8) at 192 ranks. A cruise mesh built now therefore will not reproduce the
   delivered cruise numbers exactly, and that is expected (`recipes/README.md`). The solve
   is split differently (128 ranks by default), which does not change the mesh;
   `hpc/stage_solve_case.sh` sets that up when it makes the solve case
   ([05](./05_submit_and_run.md), section 6).
2. **A production mesh does not fit on a workstation.** The delivered wake-refined CR
   meshes have 99,111,506 to 99,208,470 cells (`mesh.cells` in the same run cards), and
   meshing 68.6 million cells used 115 GB of memory (`recipes/README.md`).
3. **The set-up is small and the mesh is not.** The six registered wing surfaces are 52.6
   to 53.3 MB each, and in a local test `rsync -z` sent the baseline case, surface
   included, as 8.33 MB. On HPC12 the `constant/polyMesh` of the delivered wake-refined
   baseline mesh case alone occupies 12,162,360 KiB (about 11.6 GiB, the `ls -l` total),
   before any `processor*` directory is counted.

What goes and what stays:

| In your case directory | Send it? | Why |
|---|---|---|
| `0.orig/`, `system/` (including `system/argusPostPro`, Step 5e), `constant/momentumTransport`, `constant/physicalProperties` | yes | the recipe itself, plus the post-processing block |
| `MANIFEST.sha256`, `RECIPE_FILES.json` (and `CONDITION.json`, `variants/` in the two cruise recipes) | yes | the mesh job refuses a case without its manifest |
| `constant/triSurface/wing.stl` | yes | the only geometry input |
| `processor*/` | **never** | created on the cluster; by far the largest part of a case |
| `constant/polyMesh/` | no | the mesh is rebuilt on the cluster at 192 ranks |
| `constant/extendedFeatureEdgeMesh/`, `constant/triSurface/wing.eMesh` | no | the mesh job rebuilds them with `surfaceFeatures` |
| time directories (`0/`, `500/`, `4000/` ...) | no | the solve job creates `0/` from `0.orig/` itself |
| `postProcessing/`, `log/`, `log.*`, `*.foam`, `VTK/`, `dynamicCode/` | no | output of a local run, meaningless on the cluster |

---

## Step 1. Set the cluster variables

Replace `<your-netid>` with your own login name (no angle brackets), then paste:

```bash
export NETID=<your-netid>
export HPC=$NETID@hpc12.tudelft.net
export SCRATCH=/home/scratch/$NETID/argus
echo "$HPC"  "$SCRATCH"
```

*Expected:* your login followed by `@hpc12.tudelft.net`, then
`/home/scratch/<your-netid>/argus` with your login filled in.

**Every variable in this guide lasts only until you close the terminal.** Sending a case
and bringing its results back are often days apart, so before any step, check that the
variables it uses are set in the terminal you are typing in:

| Variable | Set in | Holds |
|---|---|---|
| `REPO` | [01](./01_setup.md) (also shown under "Before you start") | your clone of this repository |
| `NETID`, `HPC`, `SCRATCH` | Step 1 | your login, the cluster address, your scratch folder |
| `REPO_FILTER` | Step 4a | what to leave out when copying this repository |
| `NAME`, `CASE` | Step 5b | the mesh case name, and its folder on your machine |
| `PUSH_FILTER` | Step 6a | the list of files to send |
| `RNAME`, `RCASE` | Step 8a | the solve case name, and its folder on the cluster |
| `PULL_FILTER` | Step 8c | the list of files to bring back |

Step 5e runs in your Python terminal, which is a different terminal: set `REPO` and `CASE`
there too.

Other guides reuse the name `CASE` for their own case directory. In
[05](./05_submit_and_run.md) and [06](./06_trim_to_target_cl.md), `CASE` is the solve case
on the cluster, the directory this guide calls `RCASE`. In
[07](./07_postprocessing_and_report.md), `CASE` is the local copy you pull in Step 8d
(`$HOME/argus_results/<name>`). Set it again whenever you move from one guide to another.

`echo "$NAME"` (for example) prints an empty line when a variable is not set. Several
commands below start with a guard line such as `: "${NAME:?set NAME first (Step 5b)}" && \`
(some run over two lines, the first ending in `\`; paste the whole block). A guard line
does nothing when its variables are set, and stops the command with the message when one
is not, because an empty variable would point `rsync` at the wrong folder, or, for `REPO`,
at the whole of your machine.

---

## Step 2. Check that you can log in

```bash
ssh $HPC 'hostname; echo $USER'
```

*Expected:* the name of the machine you landed on (the login node prints `hpc12`), then
your netID. The single quotes matter: `$USER` is expanded on the cluster, not on your
machine.

Two prompts you may see first:

1. **The first time you connect**, `ssh` prints
   `The authenticity of host 'hpc12.tudelft.net (...)' can't be established` and asks
   `Are you sure you want to continue connecting (yes/no/[fingerprint])?`. This is `ssh`
   meeting the cluster for the first time. Type `yes` only if the name is exactly
   `hpc12.tudelft.net`; if you are unsure, ask the cluster administrators before
   answering. After `yes` it prints `Permanently added ...` and does not ask again.
2. **Unless you have set up SSH keys**, every `ssh`, `scp` and `rsync` command in this
   guide asks for your password. That is normal.

If it does not work, read the first line of the error:

| What you see | What it means | What to do |
|---|---|---|
| `Could not resolve hostname` | the name `hpc12.tudelft.net` is unknown from where you are | check the spelling, then check the name on its own (below) |
| nothing, then `Connection timed out` | the name resolves but the cluster is not reachable from where you are | ask the cluster administrators how you are meant to reach it |
| `Permission denied` | you reached the cluster and it refused your login | check your netID and your password or key |

If `HPC` is not set in this terminal, the error does not name `hpc12.tudelft.net`. `ssh`
then prints `hostname contains invalid characters`, or `Could not resolve hostname` with a
word from the command as the name, and `rsync` prints
`ssh: Could not resolve hostname : Name or service not known`, with nothing before the
colon. Repeat Step 1 in this terminal.

For more detail, run the login with `-v`:

```bash
ssh -v $HPC hostname
```

*Expected:* many lines starting `debug1:`. On a working login, a line starting
`Authenticated to hpc12.tudelft.net` comes near the end, the machine name follows on a line
of its own, and the last line is `debug1: Exit status 0`. On a failing login, read the
few lines just above the error: they show how far the connection got.

On Linux and WSL you can check the name on its own with:

```bash
getent hosts hpc12.tudelft.net
```

*Expected:* one line with an IP address and the name. No output means the name does not
resolve from your machine.

SSH access and keys (logging in without typing a password) are your own set-up, and a
private key, a key fingerprint or a password never goes into this repository, a script,
or a command line.

---

## Step 3. Create the directories on the cluster

HPC12 gives you two storage areas, and the difference decides where everything goes:

| Area | Path | Put here | Rules |
|---|---|---|---|
| home | `$HOME` on the cluster, `/home/<your-netid>` | the copy of this repository, the OpenFOAM-12 installation, results you want to keep | permanent, but has a quota |
| scratch | `/home/scratch/<your-netid>` | every case directory: meshes, `processor*`, time directories | large, but **files older than 50 days are deleted automatically and cannot be recovered** (recorded in `hpc/solve_hpc12.pbs`) |

Create the case folders the run cards record:

```bash
ssh $HPC "mkdir -p $SCRATCH/mesh $SCRATCH/solve_r2"
ssh $HPC "ls -ld $SCRATCH/mesh $SCRATCH/solve_r2"
```

*Expected:* the first command prints nothing. The second prints two lines starting with
`d`, one per directory.

Where each thing lives, and why:

| On HPC12 | Holds | Where the path comes from |
|---|---|---|
| `$HOME/argus-production/` | a copy of this repository, made in Step 4 | [05](./05_submit_and_run.md) sets `REPO=$HOME/argus-production` on the cluster, submits `$REPO/hpc/mesh_hpc12.pbs`, runs `$REPO/hpc/stage_solve_case.sh`, and submits `$REPO/hpc/solve_hpc12.pbs` with the check `$REPO/scripts/assert_case_postpro.py` |
| `$HOME/OpenFOAM/OpenFOAM-12/` | OpenFOAM-org 12 | both PBS scripts source `$HOME/OpenFOAM/OpenFOAM-12/etc/bashrc` by default |
| `/home/scratch/<your-netid>/argus/mesh/<name>/` | a case being meshed | the run cards record the meshes as `mesh/SWB_r2` and so on |
| `/home/scratch/<your-netid>/argus/solve_r2/<name>/` | a case being solved, made from a finished mesh case by `hpc/stage_solve_case.sh` ([05](./05_submit_and_run.md)) | the run cards record `case.path` as `/home/scratch/$USER/argus/solve_r2/SWB_trim` and so on |
| `$HOME/argus_hpc12/` | nothing: this route does not use it, and it does not exist on a new account (Step 4d) | `solve_hpc12.pbs` looks there for the post-processing check only when a leg is submitted without `ARGUS_GATE` |

Keep every one of these paths free of spaces and commas. `hpc/stage_solve_case.sh` refuses
a solve case path, or a clone path, that contains either, because they all go into the
`qsub -v` list.

Check that OpenFOAM-12 is where the scripts expect it:

```bash
ssh $HPC 'ls $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc'
```

*Expected:* the path printed back. `No such file or directory` means OpenFOAM-12 is not
installed there yet, and both PBS scripts will stop at once with
`REFUSED: OpenFOAM-12 not found at ...`. Building it is covered in
[01, setup](./01_setup.md), Part B, section B3.

To see how much of your quotas is used:

```bash
ssh $HPC 'quota -s; df -h $HOME /home/scratch/$USER'
```

*Expected:* `Disk quotas for user <your-netid> ...`, then one line per filesystem with the
space used, the soft quota and the hard limit. A `*` after the space used means you are
over the soft quota. `df` then names the filesystem behind each area: the one
`Mounted on /home/<your-netid>` is home, the one `Mounted on /home/scratch` is scratch.

---

## Step 4. Put a copy of this repository on the cluster

The jobs in [05](./05_submit_and_run.md) are submitted on the cluster from a copy of this
repository in your cluster home, `$HOME/argus-production`, which 05 calls `$REPO`. This
step makes that copy from your clone. The whole repository is small (about 14 MB without
its `.git/` folder, which the cluster does not need).

### 4a. Define the repository filter once

```bash
REPO_FILTER=(
  --exclude='/.git/'
  --exclude='__pycache__/'
)
```

### 4b. Dry run, then the real copy

Tested against HPC12 on 2026-09-24 (see "How the filters in this guide were checked").

Both commands below start with a guard line (Step 1). Here it matters most: with `REPO`
unset, `"$REPO/"` would be `/`, the whole of your machine.

```bash
: "${REPO:?set REPO first (Before you start)}" "${HPC:?set HPC first (Step 1)}" && \
rsync -avhz --dry-run "${REPO_FILTER[@]}" "$REPO/" "$HPC:argus-production/"
```

*Expected:* `sending incremental file list`, the first time also
`created directory argus-production`, then the files that **would** be sent
(`README.md`, `data/...`, `docs/...`, `hpc/...`, `recipes/...` and so on), ending in
`(DRY RUN)`. No path may start with `.git/`. Now copy it:

```bash
: "${REPO:?set REPO first (Before you start)}" "${HPC:?set HPC first (Step 1)}" && \
rsync -avhz "${REPO_FILTER[@]}" "$REPO/" "$HPC:argus-production/"
```

*Expected:* the same list without `(DRY RUN)`, then a `sent ... bytes` line. Keep the
trailing slash after `$REPO`, and write the destination exactly as `argus-production/`:
it is relative, so it lands in your cluster home, which is where 05 looks.

### 4c. Check the copy

Compare the files the cluster actually runs: the two job scripts, the script that makes a
solve case, and the post-processing check. Run both lines and compare:

```bash
: "${REPO:?set REPO first (Before you start)}" "${HPC:?set HPC first (Step 1)}" && \
sha256sum "$REPO/hpc/mesh_hpc12.pbs" "$REPO/hpc/solve_hpc12.pbs" "$REPO/hpc/stage_solve_case.sh" "$REPO/scripts/assert_case_postpro.py" && \
ssh $HPC 'cd $HOME/argus-production && sha256sum hpc/mesh_hpc12.pbs hpc/solve_hpc12.pbs hpc/stage_solve_case.sh scripts/assert_case_postpro.py && test -x hpc/stage_solve_case.sh && echo "stage_solve_case.sh is executable"'
```

*Expected:* four 64-character hashes from your clone, then the same four, in the same
order, from the cluster, then `stage_solve_case.sh is executable`. If that last line is
missing, 05 cannot start the script. Check it in your clone:

```bash
ls -l "$REPO/hpc/stage_solve_case.sh"
```

*Expected:* a line starting `-rwx`. If it does not, your copy of the repository lost the
permission (a download as a zip file can do that): clone it again with git, then repeat 4b
and 4c. Then check every other file by content:

```bash
: "${REPO:?set REPO first (Before you start)}" "${HPC:?set HPC first (Step 1)}" && \
rsync -avh --dry-run --checksum "${REPO_FILTER[@]}" "$REPO/" "$HPC:argus-production/"
```

*Expected:* `sending incremental file list`, then no file names, then only the
`sent ...` and `total size is ...` lines. Any file listed differs between your clone and
the cluster copy: run the copy in 4b again. Directory names ending in `/` on their own are
harmless. (Tested locally: a job script changed on the receiving side without changing its
size or date was listed by this check and missed without `--checksum`.)

Whenever you update your clone, repeat 4b and 4c: the jobs run the copy on the cluster, not
the one on your machine. **Never add `--delete`** to these commands: it removes anything
in the cluster copy that is not in your clone.

### 4d. The post-processing check: pass the shipped copy with `ARGUS_GATE`

A *leg* is one PBS solve job; a solve that needs more iterations than fit in one job
continues in a second leg, and so on. At the start of every leg, `hpc/solve_hpc12.pbs` can
run the post-processing check, `scripts/assert_case_postpro.py --before`, on the case. The
check asks whether the case will write what [07](./07_postprocessing_and_report.md) needs.
If the case meets fewer than all 8 of its requirements, the job stops the leg with
`REFUSED: post-processing gate refused` before it touches the case. A case reaches 8 of 8
once the post-processing block is installed, which you do on your machine in Step 5e.

Where the job looks for the check:

1. At the path in the variable `ARGUS_GATE`, when you give it in the `qsub -v` list. The
   route in these guides gives the shipped copy, `ARGUS_GATE=$REPO/scripts/assert_case_postpro.py`,
   in the `-v` list of **every** leg ([05](./05_submit_and_run.md), section 7b). On the
   cluster that file is `$HOME/argus-production/scripts/assert_case_postpro.py`, copied in
   4b. It needs only the Python standard library, and the HPC12 login node's `python3`
   (3.6.8) runs it.
2. Without `ARGUS_GATE`, at `$HOME/argus_hpc12/assert_case_postpro.py`, which does not
   exist on a new account: the job prints
   `NOTE: <path> absent, post-processing gate NOT run` and runs the leg unchecked. Only a
   leg submitted without `ARGUS_GATE` looks there, so you do not need to create anything
   in `$HOME/argus_hpc12/`.

The job script's gate lines were run unchanged on staged cases: with the block they gave
`---- 8 of 8 requirements satisfied ----` and the job carried on; without it they gave
`1 of 8`, then `REFUSED: post-processing gate refused`. It has been: on 2026-09-24 a solve leg on HPC12 (an early-cruise case of W, job output file of leg 1) received `ARGUS_GATE` through `-v` (it is in the job's `qstat -f` variable list) and ran the check before its first iteration, ending `8 of 8 requirements satisfied`.

### 4e. Files in `hpc/` that are not part of this route

`hpc/README_HPC.md`, `hpc/submit.sh` and `hpc/push_to_cluster.sh` are from an earlier DelftBlue/SLURM setup, depend on files not in this repository, and are not the HPC12 route. The HPC12 route uses `hpc/mesh_hpc12.pbs`, `hpc/stage_solve_case.sh`, `hpc/solve_hpc12.pbs`, and the template `hpc/templates/argusPostPro`, which `scripts/install_postpro.py` reads in Step 5e.

---

## Step 5. Build the case you will send (on your machine)

Build it **fresh from the recipe folder**. Do not send the case you made in
[03](./03_build_a_case.md) if you changed any dictionary to make it run on your machine: the
mesh job compares every dictionary with the recipe's `MANIFEST.sha256` and refuses the case
if one differs, apart from the two files that carry the angle of attack (Step 5d). (Tested
locally with the mesh job's own check: changing only the cell cap in
`system/snappyHexMeshDict` gave `ALTERED system/snappyHexMeshDict` and a refusal.)

### 5a. Pick the recipe

| Condition | Recipe folder |
|---|---|
| CR (low speed, M 0.10) | `recipes/L11_wallResolved_CR` |
| Early cruise (M 0.78) | `recipes/M6_wallModelled_CRUISE_EARLY` |
| Late cruise (M 0.78) | `recipes/M6_wallModelled_CRUISE_LATE` |

The three conditions are separate operating points, each with its own baseline. Never
difference a result from one against a result from another. (`recipes/L11_wallResolved_WT`
is not one of the three conditions, and the mesh job refuses even a fresh copy of it: its
angle check prints `ALPHA dragDir implies alpha 2.102797, Uinf implies 2.102800`, because
that recipe writes `liftDir` and `dragDir` to 6 significant figures;
[05](./05_submit_and_run.md), section 5, Gate A, shows how to mesh it anyway.)

### 5b. Copy the recipe into a new, empty directory

This example is the baseline wing, B, at Condition CR, named `B_CR`. Any name works, but
never reuse a name the project's own cases already carry on the cluster (`SWB_r2`,
`SWB_trim` and the like: SW for Condition CR, CMP for early cruise, LC for late cruise,
then the geometry letter). If you work on the account that holds the delivered cases, a
push into one of their folders overwrites its dictionaries and surface, and a mesh job
started there deletes that case's mesh before it meshes.

```bash
export NAME=B_CR
export CASE=$HOME/argus_cases/$NAME
mkdir -p $HOME/argus_cases
mkdir "$CASE"
cp -r "$REPO/recipes/L11_wallResolved_CR/." "$CASE/"
ls -A "$CASE"
```

*Expected:* these five names, in an order that depends on your system's language
settings: `0.orig`, `MANIFEST.sha256`, `RECIPE_FILES.json`, `constant`, `system`. The
cruise recipes also show `CONDITION.json` and `variants`. If `mkdir` says `File exists`,
that directory is already in use: pick another `NAME` rather than copying on top of old
files.

`NAME` and `CASE` are needed again in Steps 5e, 6 and 8. In a new terminal, set them again
with the first two lines of the block above (and the Step 1 variables).

### 5c. Put the surface in place

The mesh job looks for the surface **only** in `constant/triSurface/`, and the recipe calls
it `wing.stl`. OpenFOAM-12 itself also accepts `constant/geometry/`, but the mesh job's
surface check does not, and refuses the case (tested). The source path below is where
[02](./02_geometry_from_openvsp.md) puts B's surface when you follow its Section 5, Steps 4
to 7, with the `GEOM` folder its Section 4 sets. For another letter, for 02's Section 11
commands (which name the file `baseline_wing_only_refined_oml_placed.stl`), or for another
folder, use the path of the surface you built (its file name does not matter; its checksum
does):

```bash
mkdir -p "$CASE/constant/triSurface"
cp "$HOME/argus-geometry/B/baseline_oml_placed.stl" "$CASE/constant/triSurface/wing.stl"
sha256sum "$CASE/constant/triSurface/wing.stl"
```

*Expected:* the hash for your geometry from this table. The letter-to-candidate mapping is
read from the `candidate` and `geometry` fields of `data/run_cards/SW<letter>_r2_trim.json`
and agrees with `registry/candidates.yaml`. The last column is the file name the registry
gives the surface; your own file may be called something else.

| Letter | Candidate | sha256 | Registry file name |
|---|---|---|---|
| B | `baseline` | `20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641` | `baseline_oml_placed_laddercap.stl` |
| C | `cte_i002_c04` | `790c3ded3568fa877bf8daf1a7296d7f66e7e865e42f4a1e9b90a14501c45d3b` | `cte_i002_c04_oml_placed_laddercap.stl` |
| F | `cfft_b02_c01` | `82bf2e1e570187ab5b3dfa987f122142ab01415995ef5c47fd27838901d0037b` | `cfft_b02_c01_oml_placed_laddercap.stl` |
| H | `chc_g02_c06` | `18329aa33c0c1d971b85b40bffe588cd45eb25cba83683715816e5be29186fc3` | `chc_g02_c06_oml_placed_laddercap.stl` |
| M | `mcv2_i002_c01` | `dae22d0c4bf3cba849e39a3b59d84d1412ffed3430aa4e10c0ebdf92556941dc` | `mcv2_i002_c01_oml_placed_laddercap.stl` |
| W | `cffw_b01_c01` | `200fa64ca9d6eecc54bbecaa62c384edeab95c2560ce5f514c4a9ffc830de692` | `cffw_b01_c01_oml_placed_laddercap.stl` |

These are the files [02](./02_geometry_from_openvsp.md) reproduces, so your surface must
match this table. For B, C and M, the table's file is also, byte for byte, the file the
delivered meshes on HPC12 were built from. For F, H and W, the delivered meshes (all three
conditions) were built from an earlier copy of each surface: the closed but unwelded
version, as [02](./02_geometry_from_openvsp.md), Step 6, leaves it, before the weld of
Step 7 (read in the mesh cases on HPC12). It has the same 195,266 triangles and differs
from the table's file only along the root seam, by at most 7.8e-07 m (0.78 micrometre, in
the metre frame of `wing.stl`). Exact point matching finds 796 open edges there:
`surfaceCheck` reports `connected to one face : 796` on these copies, as it does for the
closed but unwelded baseline in [02](./02_geometry_from_openvsp.md), Step 8. The trim
solve cases (names ending `_trim`) and the run cards carry the table's welded file; the
angle-sweep solve cases carry the unwelded copy. A new F, H or W mesh built from the
table's file therefore starts from the same shape as the delivered one, but not from
identical bytes, and its root seam is welded.

The same surface serves all three conditions: for each letter, the project meshed early
and late cruise from the same file as Condition CR.

These surfaces are already in **metres**, are a half model, and extend 20 mm through the
symmetry plane (`registry/candidates.yaml`; the mesh job's surface check, run locally on
the baseline, reports `y_min=-20.000mm`). If your hash differs, you have a different file:
do not send it until you have been through [02](./02_geometry_from_openvsp.md), including
its `surfaceCheck` step. The mesh job's surface check only confirms the surface is closed
and crosses the symmetry plane by more than 10 mm; it cannot tell you the file is in the
wrong units or is the wrong wing.

### 5d. Check the case before it leaves your machine

```bash
cd "$CASE" && sha256sum -c MANIFEST.sha256
grep -E "^(numberOfSubdomains|method)" "$CASE/system/decomposeParDict"
```

*Expected:* every line of the first command ends in `OK` (19 lines for the CR recipe, 22
for either cruise recipe), and the second prints `numberOfSubdomains 192;` and
`method          scotch;`.

**The mesh does not depend on the angle of attack** (header of `hpc/mesh_hpc12.pbs`: one
mesh serves every angle). So you can do either of these. Both set the angle with
`foamDictionary`, which reads every file path relative to the folder you are in, even a
path that starts with `/`. So `cd` into the case first and give it paths such as
`system/controlDict`, never `$CASE/system/controlDict`.

1. Send the mesh case with the recipe's angle (1.30 deg for Condition CR), and set the
   angle later in each solve case that [05](./05_submit_and_run.md) makes from the
   finished mesh ([06](./06_trim_to_target_cl.md), section 6.4.2), from inside that solve
   case. At the HPC12 login prompt, `foamDictionary` works once these three lines have run
   in the same shell (checked on HPC12: afterwards `foamDictionary -help` reports
   `Using: OpenFOAM-12` and lists `-writePrecision`; before them, `foamDictionary` is not
   found):

   ```bash
   module load devtoolset/11
   module load mpi/openmpi-4.1.2
   source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc
   ```

   The `source` line prints `dirname: missing operand` and
   `Try 'dirname --help' for more information.` on HPC12. That is harmless; every job log
   shows the same pair. Type the three lines as they are: run inside a pipe, they would
   load OpenFOAM into a throwaway shell and leave yours unchanged.
2. Set the angle for the first solve now, on your machine, with 06, section 6.4.2, and
   check it with 06, section 6.4.4. Run them inside this guide's `$CASE` (`cd "$CASE"`
   first), with OpenFOAM-12 loaded as in [01](./01_setup.md), section A2. The mesh job
   accepts the new angle. Tested locally on a copy of the Condition CR recipe
   with the post-processing block installed and the angle set to 1.747609 deg, its recipe
   check printed
   `GATE A: 19 manifest entries = 17 intact + 0 altered + 0 missing + 2 alpha-bearing`.
   Every solve case copied from this mesh case starts at that angle. Any further angle (a
   sweep point or a trim step) is set in its own solve case.

If a line says `FAILED`, that file differs from the recipe. Copy the recipe again into a
new directory. The one exception: if you have already set the angle of attack in this case,
exactly two lines say `FAILED`, for `system/controlDict` and `0.orig/U`. The mesh job
accepts that, because it checks those two files for consistency with each other instead of
by hash. Any other `FAILED` line means the mesh job will refuse the case.

### 5e. Install the post-processing block

The post-processing block is one file, `system/argusPostPro`, plus one line in
`system/controlDict` that includes it. It makes the solver write what
[07](./07_postprocessing_and_report.md) needs:

1. the pressure and viscous parts of the forces;
2. the wing-surface fields p, cp, y+ and the wall-shear-stress vector;
3. four Trefftz planes behind the wing;
4. the residuals.

Install it here, on your machine, before the case goes to HPC12. The mesh job ignores the
block, the solve case made from this mesh case inherits it, and the solve job checks it
before every leg (Step 4d). A case run without it cannot produce these data afterwards.

Run this in your Python terminal ([01](./01_setup.md), section A5), because the script
needs numpy. The surface must already be in place (5c). You can set the angle of attack
(5d) before or after this step: the block depends on the flow speed, not on the angle.

```bash
: "${REPO:?set REPO first (Before you start)}" "${CASE:?set CASE first (Step 5b)}" && \
python3 "$REPO/scripts/install_postpro.py" "$CASE" && \
python3 "$REPO/scripts/assert_case_postpro.py" --before "$CASE"
```

*Expected* (baseline B at Condition CR; the first line is your case folder's name, and the
`===` line shows its full path):

```
  SWB_r2
     solver        foamRun -> rho mode 'rhoInf', rhoInf 1, pInf 0
     wall patch    wing           (from existing forces object)
     CofR          1.34 0 0       (from existing forces object)
     turbulence    kOmegaSST      -> residual fields: U p k omega
     |U| for cp    34.0000        (from 0.orig/U (via $Uinf))
     TE x 2.3363, chord 1.1597 (from ASCII STL wing.stl)
     Trefftz x     2.3942 2.6262 2.9161 3.2061
     controlDict now includes argusPostPro

=== <your case folder> (before) ===
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

Before the install, the same check ended `---- 1 of 8 requirements satisfied ----`.

How to read the numbers:

1. `TE x` and `Trefftz x` are in metres, in the frame of `wing.stl`, with x pointing
   downstream.
2. `chord` is the streamwise length of `wing.stl` (1.1597 m for B). It is not a reference
   chord. The four planes sit 0.05, 0.25, 0.50 and 0.75 of that length behind the trailing
   edge.
3. `|U| for cp` is the Condition CR speed in m/s.
4. `rhoInf 1, pInf 0` is right for an incompressible case, where p is kinematic. The
   forces the block writes are therefore per unit density.
5. For B, the file this writes is byte-for-byte the one the delivered B solve used.
6. A cruise recipe prints `rho mode 'rho'` with its own freestream density and pressure
   (kg/m3 and Pa): early cruise `rhoInf 0.412710, pInf 26495.9`, late cruise
   `rhoInf 0.310830, pInf 19374.0`. A cruise block also samples `U p T rho` on the Trefftz
   planes, where Condition CR samples `U p`.

Running the script a second time changes nothing. It prints
`controlDict already includes argusPostPro`.

Now check the manifest again:

```bash
cd "$CASE" && sha256sum -c MANIFEST.sha256 | grep -v ': OK$'
```

*Expected:* `system/controlDict: FAILED` and
`sha256sum: WARNING: 1 computed checksum did NOT match`. If you have also set the angle
(5d), `0.orig/U: FAILED` is listed too, and the warning says `2 computed checksums`. That
is expected. The mesh job does not compare `system/controlDict` with the manifest; it
checks the angle in it instead. `system/argusPostPro` is not in the manifest at all. Any
other `FAILED` line means the mesh job will refuse the case.

**If it does not:**

1. `ModuleNotFoundError: No module named 'numpy'`: this terminal's `python3` has no numpy.
   Use your Python terminal. Nothing was written.
2. `HALT <case>: no sampled wall surface and no readable ASCII STL ...`:
   `constant/triSurface/wing.stl` is missing. Do 5c first. Nothing was written.

The block does not change the mesh. The tutorial mesh of [03](./03_build_a_case.md) was
built twice on 8 cores, once with the block and once without, both times with the mesh
job's order of programs. Both gave 1,193,347 cells and the same four failed checks, the
logs differed only in times and process numbers, and no meshing program wrote a
`postProcessing` folder. A 192-rank production mesh job with the block installed has run on
HPC12 (2026-09-24, an early-cruise mesh of W, 73.7 million cells, the same four failed checks),
but with no mesh of the same surface without the block to compare, "does not change the mesh"
is **NOT VERIFIED** at production size.

The same install on copies of the early- and late-cruise recipes also gave 8 of 8; their
manifest check then lists `system/controlDict: FAILED` in the same way.

---

## Step 6. Send the case

### 6a. Define the push filter once

This is a list of rules that `rsync` reads top to bottom; the first rule that matches a
file decides it. Everything listed is sent, and the last rule drops everything else. Paste
it as one block:

```bash
PUSH_FILTER=(
  --prune-empty-dirs
  --include='/0.orig/***'
  --include='/system/***'
  --include='/variants/***'
  --include='/constant/'
  --include='/constant/momentumTransport'
  --include='/constant/physicalProperties'
  --include='/constant/triSurface/'
  --include='/constant/triSurface/*.stl'
  --include='/MANIFEST.sha256'
  --include='/RECIPE_FILES.json'
  --include='/CONDITION.json'
  --exclude='*'
)
```

Like the variables in Step 1, this lasts until the terminal closes.

### 6b. Dry run, then the real thing

Tested against HPC12 on 2026-09-24 (see "How the filters in this guide were checked").

```bash
: "${CASE:?set CASE first (Step 5b)}" "${NAME:?set NAME first (Step 5b)}" \
  "${HPC:?set HPC first (Step 1)}" "${SCRATCH:?set SCRATCH first (Step 1)}" && \
rsync -avhz --dry-run "${PUSH_FILTER[@]}" "$CASE/" "$HPC:$SCRATCH/mesh/$NAME/"
```

*Expected:* `building file list ... done`, the first time also a `created directory` line,
then a list of what **would** be sent, ending in `(DRY RUN)`. For the CR recipe it
lists `MANIFEST.sha256`, `RECIPE_FILES.json`, five files in `0.orig/`,
`constant/momentumTransport`, `constant/physicalProperties`,
`constant/triSurface/wing.stl` and twelve files in `system/`, `system/argusPostPro`
among them: 22 files. Nothing named `processor`, `polyMesh`, `postProcessing` or `log`
may appear. If one does, stop and check that you pasted the whole filter. If
`system/argusPostPro` is missing, do Step 5e first.

Now send it:

```bash
: "${CASE:?set CASE first (Step 5b)}" "${NAME:?set NAME first (Step 5b)}" \
  "${HPC:?set HPC first (Step 1)}" "${SCRATCH:?set SCRATCH first (Step 1)}" && \
rsync -avhz "${PUSH_FILTER[@]}" "$CASE/" "$HPC:$SCRATCH/mesh/$NAME/"
```

*Expected:* the same list without `(DRY RUN)`, then a `sent ... bytes` line.

Details that matter:

1. **Keep both trailing slashes.** `"$CASE/"` means "the contents of the case", and the
   destination ending in `$NAME/` is the directory they go into. `rsync` creates that last
   directory itself; its parent must exist (Step 3). Without the slash after `$CASE`, the
   filter matches nothing and the dry run ends in `total size is 0` (tested).
2. **Never add `--delete`.** Nothing here needs it, and on a case that has already been
   meshed or run it can remove cluster files.
3. **Do not push into a mesh case that has a mesh job queued or running.** The job reads
   the case's dictionaries and surface when it starts and while it runs, so files changed
   under it give a mesh of a recipe nobody recorded. Once its two input checks have
   passed, and before it meshes, it also deletes `processor*`, `constant/polyMesh` and
   `constant/extendedFeatureEdgeMesh` in the case. Check first:

   ```bash
   ssh $HPC 'qstat -u $USER'
   ```

   *Expected:* nothing at all when you have no jobs, or a table with one line per job
   of yours. Reading that table is covered in [05](./05_submit_and_run.md), section 11a.
4. To replace a single file later, `scp` works too, for example the surface:

   ```bash
   : "${CASE:?set CASE first (Step 5b)}" "${NAME:?set NAME first (Step 5b)}" \
     "${HPC:?set HPC first (Step 1)}" "${SCRATCH:?set SCRATCH first (Step 1)}" && \
   scp "$CASE/constant/triSurface/wing.stl" "$HPC:$SCRATCH/mesh/$NAME/constant/triSurface/"
   ```

   *Expected:* a progress line that starts with `wing.stl` and reaches `100%`.

### 6c. Confirm what arrived

```bash
: "${NAME:?set NAME first (Step 5b)}" \
  "${HPC:?set HPC first (Step 1)}" "${SCRATCH:?set SCRATCH first (Step 1)}" && \
ssh $HPC "cd $SCRATCH/mesh/$NAME && sha256sum -c MANIFEST.sha256 | grep -c ': OK' && sha256sum constant/triSurface/wing.stl"
```

*Expected:* the number of `OK` lines your machine gave after Step 5e: `18` for the CR recipe
(`17` if you also set the angle), `21` for a cruise recipe (`20` with the angle). Somewhere
in the output is also the `sha256sum: WARNING: ... did NOT match` line from 5e, and then
the same surface hash as in Step 5c. A smaller count than on your machine means a file is
missing or changed on the way: run the `rsync` in 6b again.

### How the filters in this guide were checked

**Against HPC12, 2026-09-24.** Every command in Steps 4, 6 and 8 was run from a WSL laptop
against HPC12, into a test directory under `/home/scratch/<your-netid>/` that was removed
afterwards:

1. Step 4 pushed the repository (14.37 MB, no `.git/`). The four hashes of 4c matched on
   both sides, `stage_solve_case.sh is executable` was printed, and the `--checksum` dry
   run listed no file.
2. Step 6 pushed a Condition CR case built as in Step 5 with the block installed, plus
   planted debris (`processor0/`, `constant/polyMesh/`, `postProcessing/`, `log/`,
   `log.blockMesh`). The dry run listed exactly the 22 files named in 6b and no debris;
   6c then gave `18` OK lines, the `WARNING: 1 computed checksum did NOT match` line and
   the surface hash `20bd6acd...`. The single-file `scp` of item 4 arrived intact.
3. Step 8e pulled a delivered mesh case's `log/` and `system/` (50 files, 1.08 MB) and no
   `constant/polyMesh/`. Step 8d pulled the Step 6 test case back byte-identical, and the
   8f `--checksum` dry run listed no file. On a delivered solve case the 8d dry run listed
   no `processor*`, no `constant/polyMesh/` and no top-level time directory other than `0/`.

**Locally, earlier.** Every filter in this guide was also run locally, between two directories on one machine. The push filter was pointed at a Condition CR case with the
post-processing block installed that also held every kind of debris a local run leaves
(`processor0/`, `processor191/`, `processors192/`, `constant/polyMesh/`,
`constant/extendedFeatureEdgeMesh/`, `constant/triSurface/wing.eMesh`, `0/`, `500/`,
`4000/`, `postProcessing/`, `log/`, `log.blockMesh`, `case.foam`, `VTK/`, `dynamicCode/`).

1. **Arrived:** exactly the 22 set-up files listed in 6b. `sha256sum -c MANIFEST.sha256`
   on the receiving side gave 18 `OK` and one `FAILED`, for `system/controlDict`.
2. **Did not arrive:** every one of the debris items above.
3. With the real baseline surface in place, the two checks the mesh job runs first passed
   on the received copy:
   `GATE A: 19 manifest entries = 17 intact + 0 altered + 0 missing + 2 alpha-bearing` and
   `GATE B constant/triSurface/wing.stl: tris=192588 open=0 over=0 y_min=-20.000mm OK`.
4. The same filter on the early-cruise recipe, block installed, delivered `CONDITION.json`
   and `variants/` as well (26 files), gave 21 `OK`, and its recipe check printed
   `GATE A: 22 manifest entries = 20 intact + 0 altered + 0 missing + 2 alpha-bearing`.
5. The repository filter in Step 4 sent everything except `.git/` and `__pycache__/`
   folders; the four files compared in 4c arrived with the same sha256 as in the clone, and
   `hpc/stage_solve_case.sh` arrived executable.
6. The pull filter in 8e, run with `NAME` empty against a stand-in `mesh/` folder holding
   two mesh cases, listed both cases' `constant/polyMesh/`; with `NAME` set it listed only
   that case's set-up files, `log/` and the job's output file. That is why 8e sets `NAME`
   and dry-runs first.
7. The pull filter in 8c, run against a stand-in finished solve case laid out as the solve
   job and the post-processing block leave it, brought back the items in the table in 8d
   and left behind every item in its right-hand column.

---

## Step 7. What happens next on the cluster

The case is now ready for the mesh job, which you submit as described in
[05](./05_submit_and_run.md), pointing it at the absolute path
`/home/scratch/<your-netid>/argus/mesh/<name>`. Three things about that job affect what
you just sent:

1. It checks the case before meshing, with two input checks that accept a case with the
   block installed. In the mesh job's output file you should see these lines (the first
   shows your angle):

   ```
       alpha 1.299999 deg, |U| 34.0000 m/s, axis lag -2.04e-08 deg
     GATE A: 19 manifest entries = 17 intact + 0 altered + 0 missing + 2 alpha-bearing
     recipe intact
     GATE B constant/triSurface/wing.stl: tris=192588 open=0 over=0 y_min=-20.000mm OK
     surface closed
   ```

   For a cruise recipe the Gate A line reads
   `22 manifest entries = 20 intact + 0 altered + 0 missing + 2 alpha-bearing`.
2. It refuses if a dictionary differs from the manifest (a line starting
   `REFUSED: the recipe dicts in this case do not match the manifest`) or if the surface
   is missing, open, or does not cross the symmetry plane (a line starting
   `REFUSED: the surface is an OPEN SHELL or does not cross y=0.`). The lines above each
   refusal name the file or the count that failed. An `ALTERED` line names a recipe file
   that no longer matches the manifest; it is never `system/argusPostPro`. Start again
   from a fresh copy of the recipe (Step 5), then repeat 5e.
3. Once its two input checks have passed, and before it meshes, it deletes `processor*`,
   `constant/polyMesh` and `constant/extendedFeatureEdgeMesh` in the case. Keep nothing of
   your own there.

Once the mesh job has finished, the mesh case becomes a solve case under `solve_r2/`
through `hpc/stage_solve_case.sh`, run at the HPC12 login prompt. It needs only bash,
coreutils, sed and grep, not OpenFOAM or Python. It copies the finished mesh case into a
new folder (it never overwrites and never deletes anything), including `system/argusPostPro`,
and writes the solve decomposition: 128 ranks by default, method `hierarchical` (8 2 8),
which is what all six delivered Condition CR solve cases carry. It then prints whether
the block is there and the `qsub` lines for the first leg, with `ARGUS_GATE` in them. The
solve case keeps the recipe's run length: for Condition CR, `endTime 4000` and
`writeInterval 500`, the first leg of the delivered Condition CR solves. On HPC12 its input
checks have been run against a delivered mesh case, and they passed. Copying, verifying,
writing the decomposition and printing the next steps have run only on a local Linux
machine: **NOT VERIFIED** on HPC12. The full procedure, its output and the `qsub` lines are
in [05](./05_submit_and_run.md), sections 6 and 7.

The job scripts you submit are the ones in your cluster copy of this repository
(`$HOME/argus-production/hpc/`, Step 4).

---

## Step 8. Bring results back

This step is usually days after Step 6, in a new terminal. Set the Step 1 variables again
first (`NETID`, `HPC`, `SCRATCH`), and `NAME` as in Step 5b if you will do 8e.

### 8a. Is it finished?

Set the case you want: the solve case name you gave `hpc/stage_solve_case.sh` in
[05](./05_submit_and_run.md). The project's Condition CR baseline solve was called
`SWB_trim`:

```bash
export RNAME=SWB_trim
: "${HPC:?set HPC first (Step 1)}" "${SCRATCH:?set SCRATCH first (Step 1)}" && \
export RCASE=$SCRATCH/solve_r2/$RNAME && \
ssh $HPC "cat $RCASE/.converged"
```

The guard line (Step 1) matters here: with `SCRATCH` unset, `RCASE` would point at a folder
that does not exist, and the `No such file or directory` below would then say nothing
about the run.

*Expected:* the marker the solve job writes when, and only when, the run reached its
`endTime` **and** passed its convergence gate (`leg <N>` is the solve job it finished in;
see Step 4d). It looks like:

```text
case SWB_trim
finalised <date and time, UTC> on <node>, leg <N>
  CONVERGENCE GATE over 200 samples: MEAN drift Cd <x> ct (<=0.05), Cl <y> ct (<=1.0), Cl SPAN <z> ct (<=1.0) -> CONVERGED
  REPORT THESE: Cd <mean>  Cl <mean>  (window means)
```

`ct` is a count, 0.0001 in a coefficient. These `Cd` and `Cl` come from the
`forceCoeffs1` function object: force divided by 0.5 x `rhoInf` x `magUInf`^2 x `Aref`,
with `Aref 0.620462` m2 (the half-model area, half the DSO-basis Sref of 1.24092 m2, so the
values equal full-wing DSO-basis coefficients), at the case's own freestream (Condition
CR: `magUInf 34.0` m/s), and resolved along the `liftDir` and `dragDir` set for the case's
angle of attack (wind axes). The `lRef 0.393957` m in the same dictionary enters only the
moment coefficients, not `Cd` or `Cl`. Values from different conditions are never
compared.

`No such file or directory` means the run has **not** finished and passed: it may still be
queued or running, it may need another leg, or it stopped early. Do not pull it as a final
result; see [05](./05_submit_and_run.md) and [06](./06_trim_to_target_cl.md). The file
name starts with a dot, so a plain `ls` does not show it; use `ls -A`.

### 8b. See how big it is

```bash
: "${RCASE:?set RCASE first (8a)}" "${HPC:?set HPC first (Step 1)}" && \
ssh $HPC "du -sh $RCASE/postProcessing $RCASE/log"
```

Most of `postProcessing/` is the post-processing block's output. It writes a wing-surface
file at iteration 0 and at every saved iteration, and `purgeWrite` removes none of them.
For the delivered baseline solve (B at Condition CR, run to iteration 6000),
`postProcessing/argusWingSurface/` holds 13 files, 7,698 MiB in all (one of them is
628,068,096 bytes), and `postProcessing/argusTrefftz/` holds 329 MiB. The logs are small
by comparison: that solve's three `log/leg<N>/foamRun.log` files are 3,785,217, 947,767
and 948,538 bytes, and the largest log file in any of the project's case folders on HPC12
scratch is 4,262,031 bytes. Check you have the space for the wing-surface files.

The project's own trim cases on HPC12 also hold an `alpha0_results/` folder: the first
(cold) leg's `log/` and `postProcessing/` at the starting angle, saved before the angle was
changed because `purgeWrite 3` would otherwise delete them (its `README.txt` says so). It is
5.5 GB in each of the six delivered Condition CR trims, and the pull filter brings it back
too: the 8d dry run of the delivered `SWB_trim` gave `Total transferred file size: 14.23G
bytes`. Size the whole case except its mesh and processor folders with:

```bash
: "${RCASE:?set RCASE first (8a)}" "${HPC:?set HPC first (Step 1)}" && \
ssh $HPC "du -sh --exclude='processor*' --exclude=polyMesh $RCASE"
```

*Expected* for the delivered `SWB_trim`: `14G`, matching the dry run.

### 8c. Define the pull filter once

The opposite of the push filter: take everything **except** the heavy parts that stay on
the cluster.

```bash
PULL_FILTER=(
  --prune-empty-dirs
  --exclude='processor*'
  --exclude='/constant/polyMesh/'
  --exclude='/constant/extendedFeatureEdgeMesh/'
  --exclude='/constant/triSurface/'
  --exclude='/dynamicCode/'
  --include='/0.orig/'
  --include='/0/'
  --exclude='/[0-9]*/'
)
```

The order matters: `0.orig/` and `0/` are let through before the last rule drops every
other directory whose name starts with a digit, which is every time directory. Rules that
start with `/` are *anchored*: they match only at the top level of the case, which is why
the trailing slashes in the commands below matter.

### 8d. Pull the solve case

Tested against HPC12 on 2026-09-24 (see "How the filters in this guide were checked").

```bash
mkdir -p $HOME/argus_results
: "${RNAME:?set RNAME first (8a)}" "${RCASE:?set RCASE first (8a)}" \
  "${HPC:?set HPC first (Step 1)}" && \
rsync -avh --dry-run --stats "${PULL_FILTER[@]}" "$HPC:$RCASE/" "$HOME/argus_results/$RNAME/"
```

*Expected:* the file list, then statistics including `Total transferred file size:`, which
is how much would come down. Nothing named `processor` and no `constant/polyMesh/` may
appear, and no time directory at the top level of the case other than `0/`. Numbered
folders **inside** `postProcessing/` (for example `postProcessing/forceCoeffs1/0/` or
`postProcessing/argusWingSurface/4000/`) are expected: each function object writes its
output under the time it started from, or under the iteration it sampled.

**Keep the trailing slash after `$RCASE`.** Without it every path gains a leading
`SWB_trim/`, the anchored rules above stop matching, and the dry run then lists
`constant/polyMesh/` and the time directories for download (tested). If the list starts
with the case name, stop and add the slash. If the list is what you want:

```bash
: "${RNAME:?set RNAME first (8a)}" "${RCASE:?set RCASE first (8a)}" \
  "${HPC:?set HPC first (Step 1)}" && \
rsync -avh --partial "${PULL_FILTER[@]}" "$HPC:$RCASE/" "$HOME/argus_results/$RNAME/"
```

`--partial` keeps a half-transferred file if the connection drops, so running the same
command again resumes rather than starting over.

On WSL, keep the destination inside your Linux home (`$HOME`), as written, not on a
Windows drive under `/mnt/`. On the project's machine, `rsync` into a Windows drive under
`/mnt/` failed with `failed to set times ... Operation not permitted` and exit code 23,
because that drive does not accept the file dates `-a` preserves.

What comes back and what stays (from the local test described after Step 6):

| Comes back | Stays on the cluster |
|---|---|
| `.converged` | `processor*/` (which is where the solve job leaves every solved field) |
| `log/` (`leg1/foamRun.log`, `leg1/decomposePar.log`, `leg1/potentialFoam.log`, later legs) | any time directory other than `0/` at the top level of the case |
| `postProcessing/` (all of it: `forceCoeffs1/` and the block's `argusForces/`, `argusResiduals/`, `argusYPlus/`, `argusWallShearStress/`, `argusWingSurface/` and `argusTrefftz/`) | `constant/polyMesh/` |
| `system/` (with `system/argusPostPro`), `constant/momentumTransport`, `constant/physicalProperties` | `constant/triSurface/` (you already have the surface) |
| `0/`, `0.orig/` | `constant/extendedFeatureEdgeMesh/` |
| `MANIFEST.sha256`, `RECIPE_FILES.json`, `SOLVE_SETUP.json` (written by `hpc/stage_solve_case.sh`), and `CONDITION.json` for a cruise case | `dynamicCode/` |
| the PBS output file of every leg, `<job name>.o<job number>`, if it was written into the case (see 8g) | |

What each of the block's folders holds is described in
[07](./07_postprocessing_and_report.md), section 3.

### 8e. Pull the mesh case's logs

Tested against HPC12 on 2026-09-24 (see "How the filters in this guide were checked").

The mesh statistics a run card reports (cell count, skewness and maximum
non-orthogonality, from `checkMesh`) and the rank count and split method the mesh was
built with (192, `scotch`) live in the **mesh** case, not the solve case. The same filter
brings back its `log/` and `system/` and leaves the mesh itself behind.

`NAME` must be the **mesh** case name you used in Step 5b. In a new terminal it is empty,
and an empty `NAME` makes the source the whole `mesh/` folder, where the anchored rules no
longer match and every mesh case's `constant/polyMesh/` would come down (tested). So set
it, then dry-run first:

```bash
export NAME=B_CR        # the mesh case name from Step 5b
: "${NAME:?set NAME first (Step 5b)}" \
  "${HPC:?set HPC first (Step 1)}" "${SCRATCH:?set SCRATCH first (Step 1)}" && \
rsync -avh --dry-run --stats "${PULL_FILTER[@]}" "$HPC:$SCRATCH/mesh/$NAME/" "$HOME/argus_results/mesh_$NAME/"
```

*Expected:* the set-up files you sent in Step 6 except the surface (`0.orig/`, `system/`,
`MANIFEST.sha256` and so on), `log/` with its six logs, and the mesh job's output file if
it was written there. No `constant/polyMesh/`, no `processor` and no path starting with a
case name may appear. Then the real pull:

```bash
: "${NAME:?set NAME first (Step 5b)}" \
  "${HPC:?set HPC first (Step 1)}" "${SCRATCH:?set SCRATCH first (Step 1)}" && \
rsync -avh --partial "${PULL_FILTER[@]}" "$HPC:$SCRATCH/mesh/$NAME/" "$HOME/argus_results/mesh_$NAME/"
ls "$HOME/argus_results/mesh_$NAME/log"
```

*Expected:* `blockMesh.log  checkMesh.log  decomposePar.log  reconstructParMesh.log
snappyHexMesh.log  surfaceFeatures.log`, the six logs the mesh job writes.

### 8f. Check the copy is complete

Run the pull again as a dry run that compares file **contents**:

```bash
: "${RNAME:?set RNAME first (8a)}" "${RCASE:?set RCASE first (8a)}" \
  "${HPC:?set HPC first (Step 1)}" && \
rsync -avh --dry-run --checksum "${PULL_FILTER[@]}" "$HPC:$RCASE/" "$HOME/argus_results/$RNAME/"
```

*Expected:* `receiving file list ... done`, then no file names at all, then only the
`sent ... received ...` and `total size ...` lines. Any file listed differs from the cluster
copy; run the pull in 8d again. Directory names ending in `/` on their own are harmless. In
the local test, this check caught a file that had been altered without changing its size
or date, which a plain re-run without `--checksum` did not. With the wing-surface files it
reads several GiB at each end, so it takes a while.

### 8g. What the run card needs

Items 1 to 4 are what `scripts/derive_run_card_fields.py` reads from a case. Items 5 to 8
are what the project fetched from the cluster to build its run cards: the mesh rank count
and `checkMesh` statistics, the force, residual and y+ histories, and the cost of each leg.
The two pulls above bring back every one of them that exists, except item 8 when the job
was submitted from somewhere other than the case. The residuals file in item 4, the
`argusForces/` history in item 6 and item 7 exist only in a case that carried the
post-processing block (Step 5e). Keep the `-a` in every pull: it preserves file dates,
and the project's run cards take the date of a run from its solver log's modification
time.

1. `system/controlDict` and `system/fvSolution`.
2. `constant/momentumTransport`.
3. `0/nut`, plus `0/k` and `0/omega` (Condition CR) or `0/nuTilda` (the cruise
   conditions). These are the `0/` the solve job made from `0.orig/`, which is why `0/` is
   let through the filter.
4. The solver log, `log/leg<N>/foamRun.log`. `scripts/derive_run_card_fields.py` reads
   its first 40 lines (the `Exec` and `Build` lines). The project's run cards took the
   date of the run from the log's modification time. Every iteration's residuals are in
   the log, but the delivered cards' final residuals came from the block's residuals file,
   `postProcessing/argusResiduals/<start>/residuals.dat`: the initial residuals of Ux, Uy,
   Uz, p, k and omega (for cruise, e and nuTilda in place of k and omega) at every
   iteration, in a new folder for each leg, named after the iteration the leg started from
   (for the delivered B, whose legs ran 1 to 4000, 4001 to 5000 and 5001 to 6000: `0`,
   `4000` and `5000`). No shipped script copies the final residuals into a card.
5. From the mesh case: `log/checkMesh.log` and `system/decomposeParDict`.
6. The force histories under `postProcessing/`: `forceCoeffs1/` (the coefficients the
   convergence gate reads) and the block's `argusForces/` (the pressure and viscous parts).
7. The y+ history, `postProcessing/argusYPlus/<start>/yPlus.dat`: the minimum, maximum and
   average y+ on the wing at every saved iteration, which fills the card's y+ fields. (The
   install in 5e replaces the cruise recipes' own `yPlus` function object, which wrote
   `postProcessing/yPlus/`.) A Condition CR case run without the block has no y+ history.
   Its y+ can then only be computed where the mesh and the solved fields are, which is
   inside the decomposed case on the cluster; [07](./07_postprocessing_and_report.md),
   section 6.2, gives the parallel command. Running it on HPC12 is **NOT VERIFIED**, and no job script for it ships. The
   pulls in this step do not make it possible, because they leave the mesh and
   `processor*/` behind.
8. The PBS output file of every leg and of the mesh job, `<job name>.o<job number>`,
   which the project's cards took their wall-clock and core-hour figures from. The
   project read the first and last time-stamped lines of each of these files. PBS writes
   each one into the directory `qsub` was run from. [05](./05_submit_and_run.md) runs every
   `qsub` from inside the case (`cd $MESH` or `cd $CASE` first), so the pulls above bring
   these files back. If you ran `qsub` from somewhere else, for example your cluster home,
   fetch them from there.

Generating and validating the card is in
[07, post-processing and report](./07_postprocessing_and_report.md).

### 8h. Time directories, `processor*` and meshes: only if you mean to

1. **Never pull `processor*`.** It is the decomposed copy of everything, the largest part
   of the case, and nothing on your machine reads it.
2. **Solved fields** are only worth bringing back if you intend to open the flow field
   itself. The solve job does not reassemble them: they exist only inside `processor*/`,
   one piece per rank. Opening them on your machine would first need OpenFOAM's
   `reconstructPar` run on the cluster to write a whole time directory (`6000/` and so on),
   and then a pull of that time directory together with `constant/polyMesh/`. That step is
   not covered by these guides. Expect many gigabytes. Fields are **unreadable without the
   exact mesh they were computed on**: they are stored cell by cell, and a rebuilt mesh
   numbers its cells differently, so fields laid on a re-meshed case are garbage, not an
   approximation.
3. The rule "never transfer a mesh you can rebuild" is about **running** new cases. For
   **reading** old fields, the mesh cannot be rebuilt.

---

## Step 9. Before the 50 days are up

Scratch deletes files older than 50 days with no recovery. Two consequences:

1. **Copy off what you need as soon as a run finishes.** The pull in Step 8 is the minimum.
   For its earlier campaign the project also kept copies of finished cases in cluster home
   (fields and post-processing, without `processor*`), checked file by file with sha256;
   the scripts that did that are not in this repository. Home has a quota, so check it
   with `quota -s` (Step 3) before copying anything large there.
2. **`rsync -a` keeps each file's original date.** A file copied today can therefore
   arrive already looking weeks old: the project's meshes copied from another cluster
   arrived carrying dates 31 days old. Assume the 50 days may count from the date the file
   carries, not from the day you copied it, and keep your own master copy of everything
   you send.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ssh` prints `Connection timed out` | the cluster is not reachable from where you are | ask the cluster administrators how you are meant to reach it |
| `rsync: mkdir "..." failed: No such file or directory (2)` | the parent directory on the cluster does not exist | run Step 3 again |
| The dry run in 6b lists `processor...` or `polyMesh` | the filter was not pasted in full, or `PUSH_FILTER` is empty in this terminal (a new terminal forgets it) | paste the block in 6a again, then repeat the dry run |
| The dry run in 6b lists no files and ends in `total size is 0` | the trailing slash on `"$CASE/"` was left off | keep both trailing slashes exactly as written |
| The dry run in 6b lists no `system/argusPostPro` | the post-processing block was not installed | do Step 5e, then repeat the dry run |
| `bash: NAME: set NAME first (Step 5b)` (or the same for `CASE`, `REPO`, `HPC`, `SCRATCH`, `RNAME`, `RCASE`) | that variable is not set in this terminal | set it again (table in Step 1), then repeat the command |
| `ssh: Could not resolve hostname : Name or service not known` (nothing before the colon), `hostname contains invalid characters`, or `Could not resolve hostname` naming a word from the command | `HPC` is not set in this terminal, in a command without a guard line | repeat Step 1 in this terminal |
| Step 5e prints `ModuleNotFoundError: No module named 'numpy'` | this terminal's `python3` has no numpy | run it in your Python terminal ([01](./01_setup.md), section A5); nothing was written |
| Step 5e prints `HALT <case>: no sampled wall surface and no readable ASCII STL ...` | `constant/triSurface/wing.stl` is missing | do Step 5c first; nothing was written |
| After Step 5e, `sha256sum -c MANIFEST.sha256` prints `system/controlDict: FAILED` | the include line the block adds | expected; see Step 5e |
| The dry run in 8d lists paths starting with the case name, including `constant/polyMesh/` or `6000/` | the trailing slash after `$RCASE` was left off | add it back before running the real pull |
| The dry run in 8e lists `constant/polyMesh/` under several case names | `NAME` is empty or wrong, so the source is the whole `mesh/` folder | `export NAME=<mesh case name>` and dry-run again |
| A pull ends with `failed to set times ... Operation not permitted` and exit code 23 | the destination is on a Windows drive under `/mnt/` (WSL) | pull into `$HOME/argus_results` inside your Linux home |
| A solve leg stops with `REFUSED: post-processing gate refused` | the check ran and the case meets fewer than 8 requirements: the block was not installed before sending | install it on the solve case as [05](./05_submit_and_run.md), section 7a, describes, check it reaches 8 of 8, then resubmit; for the next case, do Step 5e before sending |
| A solve leg's output shows `NOTE: <path> absent, post-processing gate NOT run` | `ARGUS_GATE` was missing from the `-v` list, misspelled or relative, or `REPO` was not set when you typed `qsub` | fix the `-v` list for the next leg (Step 4d, [05](./05_submit_and_run.md), section 7b) |
| In 05, `qsub` or `ls` cannot find `$REPO/hpc/mesh_hpc12.pbs`, `$REPO/hpc/stage_solve_case.sh` or `$REPO/hpc/solve_hpc12.pbs` | `REPO` is not set in your terminal on the cluster, or the copy in Step 4 was not made | on the cluster, `export REPO=$HOME/argus-production`; if `ls $REPO/hpc` then fails, repeat Step 4 |
| `$REPO/hpc/stage_solve_case.sh` says `Permission denied` | the script lost its permission to run on the way | repeat Step 4 from a clone, then check it as in 4c |
| Mesh job output shows `ALTERED <file>` or `MISSING <file>` | a dictionary was edited, or a file was not sent | rebuild the case from the recipe (Step 5) and push again |
| Mesh job output shows `GATE A: no MANIFEST.sha256` | the case did not come from a recipe folder | build from `recipes/` (Step 5) |
| Mesh job output shows `GATE B: no STL in constant/triSurface` | the surface is missing, misnamed, or in `constant/geometry/` | put it at `constant/triSurface/wing.stl` (Step 5c) |
| Mesh job output shows `GATE B ... open=` or `over=` above 0, or `FAIL` | the surface is not a closed shell, or does not cross the symmetry plane | fix the surface first; see [02](./02_geometry_from_openvsp.md) |
| `ssh $HPC "cat $RCASE/.converged"` says `No such file` | the run has not finished and passed its gate | see [05](./05_submit_and_run.md) and [06](./06_trim_to_target_cl.md) |
| A write fails in cluster home | home quota full | check with `quota -s`; cases belong in scratch, not home |
