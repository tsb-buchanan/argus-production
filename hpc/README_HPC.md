# ARGUS on DelftBlue: the transfer procedure

Written to be followed at a terminal. Literal commands in order, with the
expected output after each so you can tell whether it worked.

Two halves: **laptop** (sections 4 to 6) and **DelftBlue** (sections 7 to 10).

---

## 1. These scripts have never been executed against SLURM

Everything under `hpc/` was written without cluster access. **No script here has
been run on DelftBlue, no `sbatch` line has been submitted, no module line has
been resolved on the machine, and no `srun` has been issued.** The cluster facts
in `hpc/env.sh` (partition `compute-p1`, account `research-ae-fpt`, QOS `normal`,
48 cores and 186 GiB per node, the `2024r1 openmpi openfoam-org/7` module line,
the 5 TiB `/scratch` quota, the `00:01:00` default walltime) came from the
machine and are treated as verified; **everything built on top of them is
untested.** Job 0 exists for exactly this reason: it is a one-task, five-minute
allocation that proves the account, the partition, the QOS, the module stack and
the helper functions before a single multi-rank job is committed. It is the cheapest
proof available and it should not be skipped.

What **has** been tested, offline, on the laptop:

```
bash ~/argus-validation/hpc/selftest.sh
```

50 known-answer gates: the mesh-identity assertion with one positive and five
negatives, checksum verification against a one-byte corruption, the
`pull_results` exclusion logic against `processor0/`, `processors48/`, a
root-level time directory and a field file hidden under a non-numeric directory
name, the stale-serial-mesh detector, the OpenFOAM-absent gate, `bash -n` on
every file, a sweep for ESI-only utility names, and a check that each `--ntasks`
header matches the constant `unpack_case.sh` decomposes at. Every gate is
exercised twice: once where it must pass and once where it must fail. All 50 pass. **A check that
cannot fail is not a check.**

`package_case.sh --dry-run` has also been run against all six real case
directories, read-only.

---

## 2. The storage rule

A rule you can check yourself, and the reason it is not negotiable:

| location | what goes there | why |
|---|---|---|
| `/home/$USER/argus-validation/` | **the git repo only.** Scripts, run cards, results CSVs, the decision log. | 30.00 GiB hard quota, **25.02 GiB already used**, 51,742 files. Under 5 GiB free. Small, backed up. |
| `/scratch/$USER/argus/cases/<case>/` | **everything runnable.** Meshes, fields, time directories, `processor*`. | 5.00 TiB, 1,000,000 chunk files, currently empty. Large, fast, **wiped after 6 months, not backed up.** |

Three consequences, each of which is enforced somewhere in `hpc/`:

1. **Case bundles transfer directly to `/scratch`, never staged through
   `/home`.** A single 61 MB bundle would fit; six would not leave room to work.
   The `rsync` line in section 6 targets `/scratch`.
2. **Results come back to `/home` selectively:** `postProcessing/`, logs, the run
   card. **Never a time directory, never a `processor*` directory.**
   `pull_results.sh` refuses both by explicit rule and then verifies what
   actually landed.
3. **`/scratch` is not backed up and is wiped after six months.** Anything that
   matters comes back to the repo.

Check it yourself at any time:

```
# on the cluster
quota -s                     # or the DHPC-specific quota command
du -sh /scratch/$USER/argus
du -sh $HOME/argus-validation
```

---

## 3. What is being transferred, and the distinction that governs it

**The six solves' meshes already exist and must not be rebuilt.**

There are two different statements about decomposition and they are not the same
statement:

1. **Mesh generation in parallel IS decomposition-dependent.** `snappyHexMesh` at
   24 ranks and at 48 ranks does not produce the same mesh: castellation,
   snapping and especially layer addition are handled differently across
   processor boundaries. **Rank count is part of the mesh recipe**, on the same
   footing as the refinement levels and the layer specification, and it is now a
   mandatory run-card field alongside the decomposition method.
2. **Solving on an existing mesh is NOT.** `decomposePar` only partitions a mesh
   that already exists. The geometry is identical; only the linear-algebra
   partitioning changes. A converged steady solution agrees to well inside
   iterative tolerance.

So the route is: **reconstructed `constant/polyMesh` transferred from the laptop,
`decomposePar` at 16 on the cluster.** Not remeshed. The coverage-ladder
calibration and the measured noise floor of **0.106 to 0.199 counts at fixed CL**
(D058) were established on those exact meshes; rebuilding on the cluster would
silently invalidate the thing the ladder is being compared against.

**Rank count is a run-card field because parallel `snappyHexMesh` is
decomposition-dependent. Re-decomposing an existing mesh for a solve is not.**
Those are two separate numbers and the run card carries both under separate
names. For these six: `mesh.mesh_generation.ranks` is **3** and `hpc.solve_ranks`
is **16**. Only the first is part of the recipe; the second is recorded for cost
and reproducibility and could be changed tomorrow without touching a result.

**The solve runs on 16 ranks, not 48** (ruled 2026-07-29). These bundles are 1.20
to 1.34 M cells; at 48 ranks that is about 27 k cells/rank, where communication
overhead dominates and wall-clock can be *worse* than at fewer ranks while
consuming four times the allocation. 16 ranks puts them at 75 to 84 k cells/rank,
inside OpenFOAM's efficiency band. The earlier 48-rank sizing assumed 10 M cells
and does not apply. **48 is for the production and wall-resolved meshes**, when
they exist, and `hpc/slurm/mesh.slurm` keeps it.

A consequence worth stating, because it is a deviation from the original brief:
`solve.slurm` uses `--mem-per-cpu=4G` rather than `--exclusive --mem=0`.
`--exclusive` charges the whole node whatever `--ntasks` says, and SLURM has no
`--no-exclusive`, so leaving it in the header would have delivered 16 ranks at
48-rank cost and could not have been overridden at submit time. The DHPC warning
that motivated `--exclusive --mem=0` is about serial steps being capped inside a
`--mem-per-cpu` job, and `solve.slurm` has no serial memory-heavy step:
`decomposePar` runs in `unpack_case.sh`, outside the allocation. 16 x 4G = 64 GB
against a requirement near 3 GB. `mesh.slurm` keeps `--exclusive --mem=0`, where
the serial phases *are* the peak-memory steps and the warning applies directly.

The mesh is unchanged by that round trip, and `unpack_case.sh` **asserts** it:
cell count, face count, point count, internal-face count and a content digest of
`constant/polyMesh` must match the laptop values **exactly**, and `checkMesh` max
skewness and max non-orthogonality to a relative tolerance of 1e-6. If they do
not, it stops and reports.

**Good news you do not have to act on:** all six cases already carry a current
serial `constant/polyMesh`. The serial cell count equals the sum over the
processor directories for every one of them, so **no `reconstructParMesh` is
needed** and `package_case.sh` will say so:

```
package_case: serial mesh is CURRENT (counts match). No reconstructParMesh needed.
```

### Known future item, not part of this deliverable

When production meshes at 20 to 50 M cells are eventually built **on the
cluster**, **the noise floor must be re-measured there at the fixed rank count**,
because it will be a different mesher configuration from the laptop's. That is
two solves and a self-comparison. `hpc/slurm/mesh.slurm` exists for that day and
is clearly marked as unused now. Do not do it yet.

---

# LAPTOP HALF

## 4. One-time setup: network access, then SSH keys

### 4.0 EduVPN first, or nothing below works

The DHPC documentation is unambiguous: **"A direct SSH to DelftBlue from outside
of the university network is impossible!"** Off campus, connect **EduVPN
(Institute Access)** before anything else. The alternative is to route through
the TU Delft bastion, which every command here supports:

```
ssh -J <netid>@linux-bastion.tudelft.nl <netid>@login.delftblue.tudelft.nl
scp -J <netid>@linux-bastion.tudelft.nl ...
hpc/push_to_cluster.sh --bastion <netid>@linux-bastion.tudelft.nl
```

Students use `student-linux.tudelft.nl` instead of `linux-bastion.tudelft.nl`.

Quick check that you have a path at all, before debugging anything else:

```
getent hosts login.delftblue.tudelft.nl
```

*Expected:* three A records (`131.180.67.13/14/15`). If this resolves but `ssh`
hangs, that is the VPN, not your key.

### 4.1 SSH keys, so you stop typing your password

Password auth means every `rsync` asks, and a password typed repeatedly is a
password that ends up pasted somewhere it should not be. Fix it once.

**One caveat the DHPC docs raise, recorded because it is not obvious:** with SSH
keys "you won't have access to the `/tudelft.net` network drives on DelftBlue
until you issue the `kinit` command on the login node." Password login creates
the Kerberos ticket implicitly; key login does not. **This does not affect the
ARGUS workflow**, which uses only `/home` and `/scratch`, but if you ever reach
for the project drive after a key login, run `kinit` first.

```
ls ~/.ssh/id_ed25519.pub
```

*Expected:* the file, or `No such file`. If it is missing:

```
ssh-keygen -t ed25519 -C "argus-delftblue"
```

*Expected:* prompts for a file location (accept the default) and a passphrase
(use one), then `Your identification has been saved in $HOME/.ssh/id_ed25519`.

Copy it up. **This is the one time you type the password:**

```
ssh-copy-id -i ~/.ssh/id_ed25519.pub <netid>@login.delftblue.tudelft.nl
```

*Expected:* `Number of key(s) added: 1`.

Now add the host alias so every later command is short. Append to `~/.ssh/config`:

```
Host delftblue
    HostName login.delftblue.tudelft.nl
    User <netid>
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 60
    ServerAliveCountMax 10
```

`ServerAliveInterval` matters: a long `rsync` over a quiet link gets dropped
without it.

```
chmod 600 ~/.ssh/config
ssh delftblue 'hostname; echo $USER; df -h /scratch | tail -1'
```

*Expected:* a login-node hostname, your netid, and a `/scratch` line showing
multi-TB capacity, **with no password prompt.** If it still asks for a password,
the key did not land: re-run `ssh-copy-id`.

Add a passphrase-less session so you type the passphrase once per boot:

```
eval "$(ssh-agent -s)" && ssh-add ~/.ssh/id_ed25519
```

## 5. Build the bundle

Always dry-run first. It writes nothing and runs no OpenFOAM utility:

```
cd ~/argus-validation
bash hpc/package_case.sh cases/noise_floor/A2_base_a A2_base_a --dry-run
```

*Expected, and this is the real output from this case:*

```
package_case: case is decomposed over 3 processor directories
package_case: sum over processors: 1341934 cells; serial constant/polyMesh: 1341934
package_case: serial mesh is CURRENT (counts match). No reconstructParMesh needed.

--- included ---
    0                        24K
    0.orig                   24K
    constant                 104M
    system                   32K
--- excluded ---
    2500  (time directory, not requested)
    postProcessing  (excluded by rule)
    processor0  (excluded by rule)
    ...
package_case: mesh sha256 0e42d6c476bc6026a97f92c115eb4a10b8c87c157de8aa6e1a517486d95dc077
package_case: mesh counts nPoints 1551792  nCells 1341934  nFaces 4223796  nInternalFaces 4047927
package_case: mesh-generation ranks 3, method scotch (inferred: true)
package_case: WARNING: the mesh-generation rank count was INFERRED ...
```

**Read that warning.** The rank count is inferred from the number of `processor*`
directories present, not recorded at meshing time. A case can be re-decomposed
after meshing and the evidence is then gone. If you know what `snappyHexMesh` was
actually run at, pass it and the manifest records it as fact rather than
inference:

```
bash hpc/package_case.sh cases/noise_floor/A2_base_a A2_base_a \
     --mesh-ranks 3 --mesh-method scotch
```

**This runs `checkMesh` against the case directory.** Only do it when no solve is
in flight on that case.

*Expected tail:*

```
package_case: DONE
  bundle      $HOME/argus-bundles/A2_base_a.tar.gz
  compressed  61MB  (63477790 bytes)
  cells       1341934
  mesh sha    0e42d6c4...
```

Repeat for the other four ladder cases. Sizes measured on this repo:

| bundle | source | cells | compressed |
|---|---|---|---|
| `A2_base_a` | `cases/noise_floor/A2_base_a` | 1,341,934 | 61 MB |
| `sw_nlayer4` | `cases/layersweep/sw_nlayer4` | 1,294,446 | 59 MB |
| `lad_morph_n4` | `cases/layersweep/lad_morph_n4` | 1,293,270 | 59 MB |
| `sw_nlayer3` | `cases/layersweep/sw_nlayer3` | 1,196,259 | 55 MB |
| `lad_morph_n3` | `cases/layersweep/lad_morph_n3` | 1,195,144 | 55 MB |

### The `--with-fields` case, and a gate that fires today

Jobs 6 and 7 restart from the converged field, so they need
`--with-fields 2500`. Try it now and it stops:

```
bash hpc/package_case.sh cases/noise_floor/A2_base_a A2_base_a_wf \
     --with-fields 2500 --dry-run
```

*Expected:*

```
package_case: FATAL: time '2500' exists but is NOT a usable restart: missing U p k omega nut.
This is the real state of the A2 cases in this repo: 2500/ holds only uniform/
and yPlus.gz. Reconstruct the solved fields first:
    cd .../A2_base_a && reconstructPar -time 2500
```

That is correct behaviour, not a bug. The reconstructed `2500/` holds only
`uniform/` and `yPlus.gz`, because only the `yPlus` function object's output was
ever reconstructed. It looks like a converged time directory and would restart
from nothing.

### Jobs 6 and 7: two paths, and the choice is not the HPC workstream's

`reconstructPar -time 2500` belongs to the ARGUS workstream, which owns those case
directories, and was deliberately not run from here. Both paths below reach the
same converged answer, because this is a **steady** solve and the initial
condition does not survive convergence. Whoever runs the jobs picks one; the
literal `sbatch` lines for both are in `hpc/jobs/job_index.md`.

| | path (a): restart from the converged field | path (b): cold, from `potentialFoam` |
|---|---|---|
| prerequisite | ARGUS workstream runs `reconstructPar -time 2500` on `A2_base_a` and `A2_morph_a`, when no solve is in flight | none |
| packaging | `--with-fields 2500` | no extra flag |
| `nut` BC edit | in the `2500/` fields | in `0/nut` |
| submit with | `ARGUS_ENDTIME=5000` | `ARGUS_POTENTIALFOAM=1` |
| bundle size | about 126 MB each | 61 MB each, same as job 1 |
| iterations | fewer | more |
| dependency | on another workstream's action | none |

`solve.slurm` supports both. Under path (b) it echoes the `potentialFoam`
initialisation loudly and tells you to record it in the run card's notes:
**which path was taken is part of the run's description** even though it does not
change the converged answer. If the two members of the pair ever end up on
*different* paths, say so explicitly when the delta is reported: they still
difference correctly, but an unstated difference in initialisation between two
halves of a pair is precisely the shape of unstated frame this project keeps
getting caught by.

## 6. Transfer

Two separate transfers with different destinations, and only the first is needed
before job 0.

### 6.1 The scripts, to `/home`. This is the one to run now.

```
bash ~/argus-validation/hpc/push_to_cluster.sh --dry-run
bash ~/argus-validation/hpc/push_to_cluster.sh
```

Off campus, add `--bastion <netid>@linux-bastion.tudelft.nl`, or connect EduVPN
first (section 4.0). Without a `~/.ssh/config` block, add
`--remote <netid>@login.delftblue.tudelft.nl`.

*Expected tail:*

```
push_to_cluster: commit a21193a..., repo dirty false, hpc/ dirty false
push_to_cluster: stamped .../hpc/COMMIT
push_to_cluster: ssh OK
--- verifying what landed ---
    ok      env.sh
    ok      slurm/smoke.slurm
    ...
    all present, smoke.slurm parses
```

It transfers **`hpc/` only, 192 KB**, stamps the git commit into `hpc/COMMIT`,
creates `/scratch/$USER/argus/{bundles,cases}`, and verifies on the far side that
every file arrived and `smoke.slurm` parses under the cluster's own `bash`.

**Why only `hpc/`.** Every cluster-side script sources `${REPO}/hpc/env.sh` and
nothing else outside `hpc/`. The analysis scripts run on the laptop after
`pull_results.sh`.

**Two corrections to earlier versions of this section, since both would have
misled you at the terminal:**

1. It said the repo minus `cases/` and `.git/` was "a few MB". **It is 700 MB**,
   dominated by 519 MB of untracked geometry and 131 MB of results, none of which
   the cluster reads. Against under 5 GiB free on `/home` that is a seventh of
   the remaining space for no benefit. The claim that `cases/` is 1.5 GB was also
   wrong: **it is 24 GB.**
2. It offered `git clone` on the cluster as the preferred route. **That route is
   currently dead**: nothing is pushed (every commit including `hpc/` is
   local-only) and `github.com/<owner>/argus-validation` returns HTTP 404
   unauthenticated, so a clone would need a token *and* would not contain `hpc/`.
   If the repo is later pushed and made reachable, `git clone` becomes the better
   route and `hpc/COMMIT` stops being necessary: `argus_provenance` prefers a
   live git repository whenever it finds one.

**Provenance, and why a stamp is not a live reading.** There is no git repository
on the cluster, but `git_commit` is a mandatory run-card field.
`push_to_cluster.sh` writes `hpc/COMMIT` at transfer time and `argus_provenance`
reads it back, labelled `STAMPED at transfer time, not read live`. It records the
laptop's HEAD at the moment of transfer and cannot notice anything after it, so
**re-push `hpc/` whenever the repo moves.** This is the inferred-versus-recorded
distinction again, in a second place.

### 6.2 The case bundles, to `/scratch`. Not yet.

**These do not exist.** `package_case.sh` has only been dry-run; building them
for real runs `checkMesh` against the case directories and must wait until no
solve is in flight. Do this after job 0 passes.

```
ssh delftblue 'mkdir -p /scratch/$USER/argus/bundles'
rsync -avz --progress ~/argus-bundles/A2_base_a.tar.gz \
      delftblue:/scratch/$USER/argus/bundles/
```

*Expected:* a progress line reaching `100%`, then `sent ... bytes  received ...
bytes` and `total size is 63477790  speedup is 1.00`.

**Note the destination.** `/scratch`, not `/home`. A 61 MB bundle has no business
passing through a quota with under 5 GiB free, and there are five of them.

All five ladder bundles at once:

```
rsync -avz --progress ~/argus-bundles/*.tar.gz delftblue:/scratch/$USER/argus/bundles/
```

`rsync -av` is the DHPC-documented form; `-z` is added because these compress
further in flight. `scp -p` (and `scp -pr` for directories) is the documented
alternative if `rsync` is unavailable on either end.

---

# DELFTBLUE HALF

## 7. Log in and unpack

```
ssh delftblue
```

*Expected:* the DHPC login banner and a login-node prompt.

```
bash $HOME/argus-validation/hpc/unpack_case.sh \
     /scratch/$USER/argus/bundles/A2_base_a.tar.gz
```

*Expected, in this order:*

```
unpack_case: verifying SHA256SUMS over 24 files
unpack_case: SHA256SUMS OK
unpack_case: installed at /scratch/<netid>/argus/cases/A2_base_a
    OK: case is on /scratch
argus: module load 2024r1 openmpi openfoam-org/7
unpack_case: running checkMesh -allGeometry -allTopology (serial) -> .../log/checkMesh.cluster.serial.log

--- MESH ROUND-TRIP ASSERTION ---
  mesh assertion, manifest -> this host
    nPoints        1551792  ->  1551792
    nCells         1341934  ->  1341934
    nFaces         4223796  ->  4223796
    nInternalFaces 4047927  ->  4047927
    mesh_sha256    0e42d6c476bc6026...  ->  0e42d6c476bc6026...
    nonOrtho max   ...  ->  ...
    skewness max   ...  ->  ...
  MESH ASSERTION PASSED: cluster mesh is identical to the packaged mesh.
unpack_case: preserved the mesh-generation dictionary at system/decomposeParDict.meshgen
unpack_case: setting numberOfSubdomains 16, method scotch for the SOLVE
unpack_case: decomposed COLLATED into processors16/
```

**The `checkMesh` here is deliberately serial**, matching the serial run that
produced the manifest values on the laptop. Serial against serial removes the
serial-versus-parallel confound from the assertion. The mandatory
`checkMesh -parallel -allGeometry -allTopology` at the solve rank count is a separate,
additional record and runs inside `solve.slurm`.

**If the assertion fails,** it stops with `MESH ASSERTION FAILED` and names every
field that disagreed. Do not solve. On the two floating-point fields: if it fails
at the default 1e-6 but passes at `--tol 1e-4`, that is a floating-point finding
about the toolchain, **not** a mesh change. Report it. Do not widen the tolerance
to make it green.

**If `SHA256SUMS` fails,** the bundle is corrupt or was modified after packaging.
Re-transfer; do not install it.

## 8. Submit

Every literal `sbatch` line, in order, with what each proves and what to do if it
fails, is in **`hpc/jobs/job_index.md`**. Submit them **by hand, one at a time**.
There is no driver script, no job array and no dependency chain, deliberately:
every one of jobs 2 to 7 has a pre-registered branch attached, and each result is
meant to be read before the next is committed.

**Job 0, the whole line, nothing else needed:**

```
sbatch $HOME/argus-validation/hpc/slurm/smoke.slurm
```

*Expected:* `Submitted batch job <id>`, then within a few minutes a file
`argus-smoke-<id>.out` in the directory you submitted from, ending
`############ SMOKE TEST PASSED ############`.

**Copy the `--- module stack, resolved ---` block out of that file and into
`docs/environment.md`**, which still carries the module line as UNCONFIRMED.

Then job 1, the migration gate, alone. Then jobs 2 to 7, which are mutually
independent and can go in one sitting.

## 9. Monitor

```
squeue -u $USER
```

*Expected:* one line per job. `ST` is `PD` (pending), `R` (running) or absent
once finished. `NODELIST(REASON)` explains a pend: `Priority` and `Resources` are
normal, `AssocMaxWallDurationPerJobLimit` means the walltime exceeds what the QOS
allows.

```
sacct -j <id> --format=JobID,JobName,State,Elapsed,MaxRSS
```

*Expected after a good run:* `State` `COMPLETED`, `Elapsed` well under the
requested `--time`, `MaxRSS` a few GB per node.

`sacct` is also how you confirm a walltime kill: `State` reads `TIMEOUT` and
`Elapsed` equals the requested `--time` exactly.

```
scancel <id>
```

**When something fails, read `log_slurm/*.err` FIRST.** It carries the
shell-level error: a missing module, a failed pre-flight gate, a `set -e` abort.
`log/simpleFoam.log` carries the solver's own output and is the **second** place
to look, not the first. A solver that never started leaves nothing useful in
`log/`, and reading it first wastes the trip.

```
tail -40 /scratch/$USER/argus/cases/<case>/log_slurm/*.err
tail -60 /scratch/$USER/argus/cases/<case>/log/simpleFoam.log
```

## 10. Retrieve

On the **laptop**:

```
cd ~/argus-validation
bash hpc/pull_results.sh A2_base_a --dry-run
bash hpc/pull_results.sh A2_base_a
```

*Expected tail:*

```
--- verifying what landed ---
    no processor* directories, no time directories, no field files. OK.
pull_results: DONE
  destination $ARGUS_ROOT/results/hpc/A2_base_a
```

It brings back `postProcessing/`, `log/`, `log_slurm/`, `MANIFEST.json`,
`MESH_RECIPE.json` and the run card. It refuses `processor*`, `processors*` and
every time directory, first by rsync exclusion and then by inspecting what
actually arrived, testing by content as well as by name: a directory holding
`U`, `p`, `k`, `omega` or `nut` is a time directory whatever it is called.

`log_slurm/` is included because section 9 tells you to read it first, so it has
to be there.

---

## 11. Troubleshooting

| symptom | cause | fix |
|---|---|---|
| job killed after ~60 seconds, no useful error | **`--time` was omitted.** The `compute-p1` default is `00:01:00`. One minute. | every header in `hpc/slurm/` sets `--time` explicitly; do not remove it. Override on the `sbatch` line, not by editing the header |
| `sbatch: error: invalid account` or `Invalid account or account/partition combination` | the `--account` line | `research-ae-fpt`, or `innovation`. Check with `sacctmgr show assoc user=$USER format=account,qos` |
| `sbatch: error: ... doesn't specify the amount of memory per-cpu` | a **non-exclusive** job using `--mem=` | use `--mem-per-cpu=`. Measured 2026-07-29 |
| `sbatch: error: ... doesn't specify the amount of memory per-node. When requesting complete nodes exclusively` | an **`--exclusive`** job using `--mem-per-cpu=` | use `--mem=0`. The rule is conditional on `--exclusive`, in both directions |
| `sbatch: error: ... memory per CPU (4096 MB) exceeds the available memory per CPU (3968 MB)` | over compute-p1's per-CPU ceiling | 3968 MB = RealMemory 190464 / 48 cores. It binds however few tasks you ask for; the node total is not the constraint. `--mem-per-cpu=4G` is over it |
| a **solve** job dies with an out-of-memory kill, or `MaxRSS` near the request | `--mem-per-cpu` too low for the cell count | `solve.slurm` requests 4G/rank, i.e. 64 GB at 16 ranks against roughly 3 GB needed at 1.34 M cells. If a much larger case is ever run through it, raise it on the `sbatch` line: `--mem-per-cpu=8G` |
| a **mesh** job dies with an out-of-memory kill | **the `--exclusive` and `--mem=0` pair is missing.** Without it SLURM caps the job at a default per-CPU share, and the serial phases hit that cap first | `--exclusive --mem=0` **together** is the documented DelftBlue idiom for taking the whole node's 186 GiB. Both, not one. `mesh.slurm` has both; do not remove either |
| `simpleFoam: command not found`, or any foam utility not found | **`hpc/env.sh` was not sourced**, or the module line did not resolve | every `.slurm` script sources it and then calls `argus_require_foam`, which fails loudly. If you are running a utility by hand: `source $HOME/argus-validation/hpc/env.sh && argus_load_modules` |
| `Lmod has detected an error` on the module line | wrong stack year or module name | `module spider openfoam` on the login node, then fix `ARGUS_MODULE_LINE` in `hpc/env.sh` and record the working line in `docs/environment.md` |
| `attempt to run parallel on 1 processor` | MPI mismatch, or `numberOfSubdomains` disagrees with `--ntasks` | `solve.slurm` checks the second before submitting work. For the first, resubmit job 0 with `--ntasks=8`: it decomposes the cavity and runs `-parallel`, which is the MPI proof |
| `MESH ASSERTION FAILED` | the mesh on the cluster is not the mesh that was packaged | stop. Do not solve. Re-transfer and re-unpack. If only the two floating-point fields disagree, see section 7 |
| checksum mismatch at unpack | truncated or corrupted transfer | re-run the `rsync`; it resumes |
| `startFrom latestTime, latest time 2500 >= endTime 2500` | a restart with nothing to do | add `ARGUS_ENDTIME=<t>` to the `--export` list |
| solve finishes suspiciously fast and `postProcessing/` is empty | function objects missing | `solve.slurm` refuses to start without both `yPlus` and `forceCoeffs`; if you see this, the case was run some other way |

---

## 12. The line to add to every run card

**Decomposition method and rank count. Two of each, and they are different
numbers:**

```json
"mesh": {
  "mesh_generation": {
    "ranks": 3, "method": "scotch", "inferred": true,
    "inferred_from": "count of processor* directories at package time; method from system/decomposeParDict",
    "machine": "laptop-WSL2"
  }
},
"hpc": {
  "machine": "DelftBlue", "partition": "compute-p1", "account": "research-ae-fpt",
  "nodes": 1, "solve_ranks": 16, "solve_decomposition_method": "scotch",
  "file_handler": "collated", "walltime_requested": "02:00:00",
  "bundle_mesh_sha256": "0e42d6c4..."
}
```

`mesh.mesh_generation.ranks` **is part of the mesh recipe.** `hpc.solve_ranks`
**is not.** For all six of these solves those numbers are 3 and 16 respectively:
three ranks built the mesh on the laptop, sixteen solved it on the cluster.
A single field called `ranks` would have to hold one or the other and would be
wrong either way: a correct number with an unstated frame.

Both values are printed by `solve.slurm` in its `RUN-CARD FIELDS` block and
carried in `MESH_RECIPE.json` inside every bundle.

**The schema change is drafted but NOT applied.** `scripts/run_card.schema.json`
belongs to another workstream and is not edited from here. The exact fragment and the
three edits needed to merge it are in **`hpc/run_card_hpc_fields.json`**.

---

## 13. Two settled points, recorded so they are not re-litigated

**(a) `surfaceFeatures` is the org utility. RULED 2026-07-29, and the HPC workstream's
charter had it inverted.**

`surfaceFeatureExtract` was **replaced by `surfaceFeatures` in OpenFOAM-org at
version 6**; ESI kept the old name. The charter under which `hpc/` was written
stated the opposite (that `surfaceFeatureExtract` was the org-7 utility and
`surfaceFeatures` the ESI one) and instructed that the ESI name never be
written. **That instruction was wrong and is superseded.** It is recorded here,
in the header of `hpc/slurm/mesh.slurm`, and in D078, so that nobody later
"corrects" this repo back to the wrong name on the strength of the charter text.

The evidence, measured on the pinned local install:

1. **Both binaries ship in org-7**, in
   `/opt/openfoam7/platforms/linux64GccDPInt32Opt/bin`, and they read different
   dictionaries: `surfaceFeatures` reads `system/surfaceFeaturesDict`,
   `surfaceFeatureExtract` reads `system/surfaceFeatureExtractDict`.
2. **All 16 org-7 tutorials that extract features use `surfaceFeatures`.** Zero
   use `surfaceFeatureExtract`.
3. **This repo already follows the tutorials.**
   The OpenFOAM-org 7 wing template's `Allmesh` (no longer shipped) called `surfaceFeatures`, with the
   comment `# org-7 utility (ESI: surfaceFeatureExtract - do not use)`, and every
   existing 3D case carries `system/surfaceFeaturesDict` and
   `constant/extendedFeatureEdgeMesh`.

**The dictionary-selection design in `mesh.slurm` is kept anyway**, and that was
also ruled. `mesh.slurm` hardcodes neither name: it selects the utility from the
dictionary the case actually carries and refuses to guess if the case carries
both or neither. Hardcoding even the correct name assumes every case will always
carry the matching dictionary. Selecting from what is actually there cannot go
wrong, and it is what surfaced the discrepancy in the first place.

Nothing in this deliverable depended on the answer either way, because none of
the six meshes is rebuilt.

**(b) `cluster/delftblue/` is superseded and NEEDS A DEPRECATION POINTER.**

The repo contains `cluster/delftblue/00_smoke_test.sbatch`,
`01_mesh_serial.sbatch` and `02_solve_parallel.sbatch`, written from the DHPC
docs on 2026-07-26 with `partition`, `account` and `qos` left as marked
`REPLACE_ME` placeholders. Those placeholders are now answered
(`compute-p1`, `research-ae-fpt`, `normal`), and the sizing guidance is
superseded twice over: it recommends a **2-node default** for the fine level,
whereas the measured capacity of about 2.3 GB per million cells puts all three
planned grid levels on **one node**, and its `--ntasks-per-node=48` is the rank
count that ruling 2 moved off for cases this size.

Those three files were not edited, moved or deleted: `hpc/` is the only directory
written from here, and **`cluster/` is the ARGUS workstream's path.** **Two
directories for the same job is how someone uses the wrong one in three months.**
A pointer at the top of each of the three files, or a `cluster/delftblue/README`
saying "superseded by `hpc/`, see `hpc/README_HPC.md`", is the minimum. That
edit is the ARGUS workstream's to make.

---

## 14. File map

| file | side | what it does |
|---|---|---|
| `README_HPC.md` | both | this document |
| `env.sh` | both | cluster constants, module loading, and the mesh-identity helpers. Sourced, never executed |
| `slurm/smoke.slurm` | cluster | job 0. 1 task, 5 minutes, proves the environment |
| `slurm/solve.slurm` | cluster | jobs 1 to 7. Generic solve on an existing mesh. Never meshes |
| `slurm/mesh.slurm` | cluster | **NOT USED FOR THIS DELIVERABLE.** Production meshes only, later. Refuses to run on a transferred bundle |
| `package_case.sh` | laptop | builds a self-verifying bundle |
| `unpack_case.sh` | cluster | verifies, installs, decomposes, and asserts the mesh round trip |
| `pull_results.sh` | laptop | selective retrieval, with a post-transfer verifier |
| `push_to_cluster.sh` | laptop | puts `hpc/` on the cluster, stamps the commit, verifies the far side. **the user runs it; it takes no credentials** |
| `selftest.sh` | laptop | 50 known-answer gates, offline |
| `run_card_hpc_fields.json` | n/a | proposed run-card schema addition. **Not applied** |
| `jobs/job_index.md` | n/a | the eight literal `sbatch` lines, in order |
