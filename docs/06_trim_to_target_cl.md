# 6. Trim every geometry to its target lift coefficient

> **At a glance**
>
> 1. **What this guide does:** brings one wing to its target lift coefficient (C_L) at one operating point, by changing the angle of attack in an OpenFOAM-12 case, and tells you when a run counts as finished.
> 2. **Before you start:** OpenFOAM-12 ([./01_setup.md](./01_setup.md)), a wing surface built from the design study's `.vsp3` file ([./02_geometry_from_openvsp.md](./02_geometry_from_openvsp.md)), a case that meshes and runs ([./03_build_a_case.md](./03_build_a_case.md)), sent to HPC12 ([./04_send_to_hpc.md](./04_send_to_hpc.md)), meshed there and made into a solve case with `hpc/stage_solve_case.sh` ([./05_submit_and_run.md](./05_submit_and_run.md)).
> 3. **At the end you have:** a solved case whose `.converged` marker reports a C_L within 1 count (`|dC_L| < 1e-4`) of the target, and the record listed in section 6.9.
> 4. **How long:** every angle you try costs at least one solve leg on HPC12. The delivered 4000-iteration Condition CR legs took 20.2 to 22.1 hours on four typej nodes and about 11.3 hours on one typen node; [./05_submit_and_run.md](./05_submit_and_run.md), section 8, gives the details.

This page explains how each wing is brought to its target lift coefficient (C_L) before
its drag is compared, how the angle of attack is actually changed in an OpenFOAM-12 case,
how to decide that a run is finished, and one trap in the HPC12 job script that can leave
a finished leg without its `.converged` marker.

Before you start you need:

1. A case that meshes and runs (see [./03_build_a_case.md](./03_build_a_case.md)) and that
   you can submit on HPC12 (see [./05_submit_and_run.md](./05_submit_and_run.md)).
2. OpenFOAM-12 loaded in your shell (see [./01_setup.md](./01_setup.md)). On HPC12, load
   it at the login prompt with the three lines of section 6.5.4, step 1
   ([./01_setup.md](./01_setup.md), section B4, has the full check).
3. **GEOMETRY NOT INCLUDED.** No wing surface is shipped in this repository. You build it
   yourself from the design study's `.vsp3` file (see
   [./02_geometry_from_openvsp.md](./02_geometry_from_openvsp.md)).

`$REPO` is your clone of this repository. On HPC12 it is the copy of the repository in your
cluster home, for example `export REPO=$HOME/argus-production` (or wherever you put the copy on HPC12),
set as in section 0 of [./05_submit_and_run.md](./05_submit_and_run.md). `$MESH` is a
finished mesh case and `$CASE` is the solve case you are working on. On HPC12 both are
absolute paths, laid out as in 05, for example `/home/scratch/$USER/argus/mesh/<name>` and
`/home/scratch/$USER/argus/solve_r2/<case name>`.
Commands marked **NOT VERIFIED** could not be checked when this page was written; treat
them with care. Replace anything in angle brackets, such as `<job id>`, before running a
command.

**Frame for every coefficient on this page.** C_L and C_D are the values OpenFOAM's
`forceCoeffs1` function object prints: half model, `Aref 0.620462 m2` (half of the
reference area of the design study, the design-space optimisation or "DSO",
`Sref 1.24092 m2`), `lRef 0.393957 m`, so they equal full-wing coefficients on the DSO
basis. The other basis this page mentions, the TP-1580 basis, uses the reference area and
chord of NASA TP-1580, the report of this wing's low-speed tests. Both bases are listed in
[../README.md](../README.md), "Reference quantities and units". One **count** is 0.0001 of
a coefficient ("ct"). Each number below also names its operating point and its mesh
generation. Every angle of attack is measured from the model's x axis, not from the wing
chord (section 6.4).

---

## 6.1 Why a trim is needed

Every geometry is compared **at the same C_L, not at the same angle of attack**. A morphed
wing produces the target lift at a different incidence from the baseline. If you compared
two wings at one fixed angle, part of the drag difference would simply be a lift
difference.

How large the incidence difference is, read from the delivered run cards
(`data/run_cards/`, Condition CR, wake-refined mesh):

| Geometry | Run card | Trim angle of attack | C_L reached |
|---|---|---|---|
| B (candidate `baseline`) | `SWB_r2_trim.json` | 1.747608659 deg | 0.428270799 |
| C (candidate `cte_i002_c04`) | `SWC_r2_trim.json` | 1.268311286 deg | 0.428278932 |

Both sit on the same target (0.428277635108) to better than 1 count, yet C needs 0.479 deg
less incidence than B to get there.

**Target tolerance.** A trim is accepted when the converged C_L is within
`|dC_L| < 1e-4` (1 count) of the target.

---

## 6.2 The three operating points

| | Condition CR (low speed) | Early cruise | Late cruise |
|---|---|---|---|
| Case-name prefix in `data/` | `SW` | `CMP` | `LC` |
| Mach | 0.10 | 0.78 | 0.78 |
| Freestream speed \|U\| (`magUInf`) | 34.0 m/s | 233.934028 m/s | 230.501789 m/s |
| **Target C_L (applied on the DSO basis)** | **0.428277635108** | **0.529297087450** | **0.544114800210** |
| Turbulence model and wall | kOmegaSST, wall resolved, y+ ~1 | Spalart-Allmaras, wall modelled | Spalart-Allmaras, wall modelled |
| Recipe folder | `recipes/L11_wallResolved_CR` | `recipes/M6_wallModelled_CRUISE_EARLY` | `recipes/M6_wallModelled_CRUISE_LATE` |
| Angle the recipe ships with | 1.3000 deg | 0.8654766754810572 deg | 0.9966996379212477 deg |
| Baseline trim case (delivered) | `SWB_trim` | `CMPB_trim` | `LCB_trim` |
| Mesh generation of the delivered numbers | wake-refined | published (earlier) | published (earlier) |

Where the numbers come from: the target C_L values are the project's pinned values (the
cruise recipes carry the same two numbers as `0.52929708745` and `0.54411480021` in
`CONDITION.json`, which is the same value written with one fewer trailing zero). The
speeds are `magUInf` in each recipe's `system/controlDict`. The CR angle is derived from
the `Uinf` vector in `recipes/L11_wallResolved_CR/0.orig/U`; the cruise angles are
`alpha.value_deg` in each `CONDITION.json`. The baseline case names are `baseline_case` in
`data/rans_forces.json`.

The Condition CR target is defined on the DSO basis. For the two cruise targets that is an
open question. `data/vlm.json` places the design study's cruise targets on the TP-1580
basis, and `data/rans_forces.json` records the disagreement as unresolved
(`_open_cross_file_frame_conflict`; [../README.md](../README.md), "Reference quantities and
units", item 7). The delivered RANS cruise trims were trimmed to these numbers on the DSO
basis.

> [!WARNING]
> **THE THREE OPERATING POINTS ARE NEVER DIFFERENCED AGAINST EACH OTHER.**
>
> Each operating point carries **its own baseline trim** (B at that condition), and every
> drag difference is "geometry minus B, **same condition**, same mesh generation". Never
> subtract an early-cruise number from a Condition CR number, never subtract a late-cruise
> number from an early-cruise number, and never compare a percentage from one condition
> with a percentage from another. The conditions differ in Mach number, target C_L,
> turbulence model and wall treatment, and (for the delivered numbers) in mesh generation.
> A difference across them measures all of that at once and means nothing.
>
> This is the easiest mistake a new user can make with this data. If a table you are
> building has columns from two conditions, it must never have a column that subtracts
> one from the other.

`recipes/L11_wallResolved_WT` is a fourth recipe (wind-tunnel condition, |U| 40.8 m/s).
It has no delivered trims and this page does not trim it. One thing on this page
still applies to it: as shipped, its force axes are written to 6 significant figures, and
the mesh job's recipe check (Gate A) refuses it. If you mesh it, first set its angle as in
section 6.4.2 with `UMAG=40.8`, even if you keep its 2.1028 deg
([./05_submit_and_run.md](./05_submit_and_run.md), section 5). On a fresh copy, the
`sed` route of section 6.4.2 with `ALPHA=2.1028` made Gate A print
`axis lag -2.14e-12 deg` instead of refusing.

---

## 6.3 Which mesh each delivered trim is on

1. **Condition CR results are on the wake-refined mesh** (the six `SW*_trim` cases whose
   run cards are in `data/run_cards/`).
2. **Early- and late-cruise results are on the published (earlier) mesh.** They were not
   re-run.
3. **The recipes in `recipes/` build the wake-refined mesh for all three conditions**
   (see `recipes/README.md`).

A trim on the wake-refined mesh and a trim on the published mesh are **not comparable**,
even for the same geometry at the same condition. For B at Condition CR
(`data/rans_forces.json`, `trims.condition_CR.geometries.B`):

| B, Condition CR | Trim angle | C_D |
|---|---|---|
| published mesh (superseded) | 1.742369 deg | 191.006 ct |
| wake-refined mesh (current) | 1.747608659 deg | 188.031 ct |

The mesh change alone moved B's trimmed drag by -2.975 counts. The whole Condition CR
field of six geometries spans 1.100 counts, so the mesh effect is 2.7 times larger than
the differences being measured.

Consequences for you:

1. If you trim a geometry at early or late cruise with the shipped recipe, your mesh is
   wake-refined and **your result cannot be differenced with the delivered early- or
   late-cruise numbers**. Trim B yourself on the same mesh generation and difference
   against that.
2. At Condition CR, the bracket cases in `data/rans_forces.json` (`SW*_a1p*`) are on the
   published mesh. The wake-refined Condition CR trims are the six whose run cards are in
   `data/run_cards/`. The published bracket slope was used only to size the first angle
   step on the new mesh (section 6.5.3).
3. The per-case rows in `data/rans_forces.json` carry a `mesh_generation` field
   (`wake_refined` or `published`). Read it before you put two numbers side by side. At
   Condition CR, the six `SW*_trim` rows read `wake_refined` and are the current results;
   every Condition CR row that reads `published` belongs to the superseded mesh and must
   not be set beside them. In the run cards, `mesh.mesh_generation` means something else:
   how the mesh was built (ranks, decomposition method, machine).

---

## 6.4 How the angle of attack is set: rotate the freestream, not the wing

The wing geometry is **never rotated**. The angle of attack alpha is imposed by tilting the
incoming flow in the x-z plane at a fixed speed |U|. The axes are: x downstream, y along
the span (the symmetry plane is y = 0), z up.

**Frame of every angle on this page.** Alpha is the angle between the freestream and the x
axis of the model's coordinate system. The wing surface already includes its built-in
1.873 deg incidence (a rotation about y, applied when the geometry is placed; see
[./02_geometry_from_openvsp.md](./02_geometry_from_openvsp.md)). So alpha is **not** the
angle between the flow and the wing root chord, and it is not a wind-tunnel angle.

Alpha is written into **three entries that must always agree** (`UMAG` below is the
freestream speed |U| of the operating point):

| File | Entry | Value for angle alpha |
|---|---|---|
| `0.orig/U` | `Uinf` | `(UMAG*cos(alpha)  0  UMAG*sin(alpha))` |
| `system/controlDict` | `functions/forceCoeffs1/liftDir` | `(-sin(alpha)  0  cos(alpha))` |
| `system/controlDict` | `functions/forceCoeffs1/dragDir` | `(cos(alpha)  0  sin(alpha))` |

`magUInf` in the same `forceCoeffs1` block stays equal to |U| and is **not** changed. In
`0.orig/U` every farfield value refers to `$Uinf`, so changing `Uinf` changes them all.

Why the force axes matter: `forceCoeffs1` projects the wing force onto `liftDir` and
`dragDir`. If they lag the freestream by an angle, C_D is biased by C_L times that angle
(in radians). An earlier campaign's axes lagged by up to 0.0083 deg, which at C_L 0.43 is
about 0.6 counts of drag, larger than several of the differences this project reports,
and invisible in residuals and force histories. The delivered post-processing also reads
each case's angle back from the `dragDir` line in the header of `forceCoeffs.dat`
(`scripts/harvest_rans_forces.py`), so a wrong `dragDir` is also a wrong recorded angle.

**No script in this repository sets alpha on, or trims, a 3D wing case built from
`recipes/`.** `scripts/build_2d_case.py` sets alpha, but only for 2D sections, and does not
apply to the cases on this page. `scripts/trim2d_transonic.py` is an OpenFOAM-org 7 tool
that must not be run ([./01_setup.md](./01_setup.md), section A8). The project's
own 3D trim tools (`make_trim_case.py`, `trim_bracket.py`, and two helper scripts that
lived only on the cluster) are **not in this repository**. `hpc/stage_solve_case.sh`
copies the angle of the mesh case into the solve case and prints it, but never changes it.
The manual procedure below does the edits with `foamDictionary`. It was checked under
OpenFOAM-12 on copies of `recipes/L11_wallResolved_CR` and
`recipes/M6_wallModelled_CRUISE_EARLY` with the post-processing block installed, and on
solve cases made from them with `hpc/stage_solve_case.sh`.

### 6.4.1 Check your `foamDictionary` first

```bash
foamDictionary -help | grep -E "Using|writePrecision"
```

You should see both of these lines:

```
  -writePrecision <label>
Using: OpenFOAM-12 (see https://openfoam.org)
```

At the HPC12 login prompt, after the three loading lines of section 6.5.4, step 1, the same
command prints the same two lines.

If `Using:` names another version, or the `-writePrecision` line is missing, your shell
has sourced a different OpenFOAM version. Stop and fix that first
([./01_setup.md](./01_setup.md)); every command below depends on `-writePrecision`.

> [!WARNING]
> **Run `foamDictionary` from inside the case, and give it a path relative to the case.**
> `foamDictionary` joins every file path you give it to the current directory, even a
> path that starts with `/`. Tested under OpenFOAM-12, from inside a case directory:
>
> 1. `foamDictionary system/controlDict -entry endTime -value` printed `4000`.
> 2. `foamDictionary $CASE/system/controlDict -entry endTime -value` stopped with
>    `FOAM FATAL IO ERROR` and `file "<current directory>//<the path you typed>" does not exist`.
>
> The same failure was seen on HPC12. So every command on this page first goes into the
> case with `cd $CASE`, and then uses paths such as `system/controlDict`, `0.orig/U` and
> `processor0/4000/U`.

> [!WARNING]
> **Always pass `-writePrecision 12` when `foamDictionary` changes a file.**
> `foamDictionary -set` rewrites the **whole** file, and without `-writePrecision` it
> rewrites **every** number in it at 6 significant figures. Tested under OpenFOAM-12:
>
> 1. Changing only `endTime` in `system/controlDict` silently shortened `liftDir` from
>    `( -0.022687323 0 0.999742610 )` to `( -0.0226873 0 0.999743 )`.
> 2. Changing one boundary value in an ASCII field file (the cruise recipes write ASCII)
>    rewrote the whole solved velocity field at 6 figures, for example
>    `(233.90734 1.9820989e-14 3.533537)` became `(233.907 1.9821e-14 3.53354)`.
>
> With `-writePrecision 12` both files kept every digit. `foamDictionary` also removes
> the comments from the file it rewrites. That is harmless in a case directory, but
> **never run it on the files in `recipes/`**: always work on a copy
> (`recipes/*/MANIFEST.sha256` would no longer verify).

### 6.4.2 Setting alpha on a case that has not run yet (a cold start)

Run these in the case directory, with OpenFOAM-12 loaded. The example sets B's delivered
Condition CR trim angle.

Where to run them:

1. **The first angle of a new mesh:** either on your own machine, in the mesh case before
   you send it ([./04_send_to_hpc.md](./04_send_to_hpc.md), Step 5), or on HPC12 in the
   solve case. `hpc/stage_solve_case.sh` copies `0.orig/` and `system/` from the mesh case,
   so every solve case staged from it starts at the mesh case's angle, and the script
   prints that angle under `NEXT STEPS`, item 2.
2. **Every further cold-start angle** (a sweep point or a fitted angle, section 6.5.2): in
   a solve case of its own, staged from the same mesh case, at the HPC12 login prompt with
   OpenFOAM-12 loaded as in section 6.5.4, step 1.
3. **A new angle for a case that has already run** is not a cold start. Use section 6.5.4
   instead: it continues the case from its converged flow.

How this was tested: on a local OpenFOAM-12 installation, on a copy of the Condition CR
recipe with the post-processing block installed, and on solve cases staged with
`hpc/stage_solve_case.sh`. On HPC12, `foamDictionary` runs at the login prompt once
OpenFOAM-12 is loaded, and the Python parts of step 1 and of section 6.4.4 ran unchanged
there, inside a delivered solve case, under the login node's `python3` (3.6.8). Section
6.4.4 shows what its check printed there.

1. Choose the angle and the speed of the operating point (section 6.2 table), then compute
   the three vectors:

   ```bash
   cd $CASE
   ALPHA=1.747608659    # angle of attack in degrees
   UMAG=34.0            # |U| in m/s: 34.0 (CR), 233.934028 (early cruise), 230.501789 (late cruise), 40.8 (WT)
   read UX UZ LX LZ DX DZ <<< "$(python3 -c "
   import math
   a = math.radians($ALPHA); u = $UMAG
   print('%.12g %.12g %.12g %.12g %.12g %.12g' % (u*math.cos(a), u*math.sin(a), -math.sin(a), math.cos(a), math.cos(a), math.sin(a)))
   ")"
   echo "Uinf ($UX 0 $UZ)   liftDir ($LX 0 $LZ)   dragDir ($DX 0 $DZ)"
   ```

   Expected output for this example:

   ```
   Uinf (33.9841853945 0 1.03689105977)   liftDir (-0.0304967958757 0 0.999534864545)   dragDir (0.999534864545 0 0.0304967958757)
   ```

   This is a known-answer check: the delivered `SWB_trim` case printed
   `dragDir (0.999534865 0 0.030496796)` in its force file header
   (`data/rans_forces.json`, case `SWB_trim`), which is the same vector rounded to nine
   figures. If Python prints an error, or the line shows `( 0 )` or blank vectors instead
   of three numbers each, then `python3` is missing or `ALPHA`/`UMAG` are not set. Do not
   go on to step 2 until the line looks like the one above.

2. Write the three entries. The paths are relative to the case, as section 6.4.1 requires:

   ```bash
   foamDictionary 0.orig/U -entry Uinf -set "($UX 0 $UZ)" -writePrecision 12
   foamDictionary system/controlDict -entry functions/forceCoeffs1/liftDir -set "($LX 0 $LZ)" -writePrecision 12
   foamDictionary system/controlDict -entry functions/forceCoeffs1/dragDir -set "($DX 0 $DZ)" -writePrecision 12
   ```

   Each command prints one line starting `New entry` followed by the value it wrote. After
   the three commands the files contain (Condition CR recipe copy):

   ```
   0.orig/U:
   Uinf            ( 33.9841853945 0 1.03689105977 );

   internalField   uniform $Uinf;

   system/controlDict, inside functions { forceCoeffs1 { ... } }:
           magUInf         34;
           liftDir         ( -0.0304967958757 0 0.999534864545 );
           dragDir         ( 0.999534864545 0 0.0304967958757 );
   ```

   `magUInf 34.0` is rewritten as `34`; that is the same number. If the post-processing
   block is installed, `foamDictionary` also rewrites its include line as
   `#include        "argusPostPro"` (with more spaces). That is still the include line:
   afterwards `python3 $REPO/scripts/assert_case_postpro.py --before $CASE` still ended
   `---- 8 of 8 requirements satisfied ----`, on the Condition CR and on the early-cruise
   test case. On a mesh case edited this way, the mesh job's recipe check (Gate A, run
   unchanged on a local copy) printed
   `alpha 1.747609 deg, |U| 34.0000 m/s, axis lag +4.86e-12 deg` and
   `GATE A: 19 manifest entries = 17 intact + 0 altered + 0 missing + 2 alpha-bearing`.

   If you would rather not load OpenFOAM for this, make the same three edits with `sed`
   instead, in the same terminal as step 1. `sed` needs no OpenFOAM:

   ```bash
   sed -i "s/^Uinf .*/Uinf            ($UX 0 $UZ);/" 0.orig/U
   sed -i "s/^\( *liftDir \+\)(.*/\1($LX 0 $LZ);/; s/^\( *dragDir \+\)(.*/\1($DX 0 $DZ);/" system/controlDict
   grep -E '^Uinf' 0.orig/U
   grep -E '^ *(liftDir|dragDir)' system/controlDict
   ```

   The two `grep` lines must print the three vectors step 1 printed, one line each. Tested
   on local copies of the Condition CR and early-cruise recipes with no OpenFOAM loaded,
   and again on a Condition CR copy with the post-processing block installed: section
   6.4.4 then printed one angle three times with `|Uinf|` equal to `magUInf`, and
   `sha256sum -c MANIFEST.sha256` reported only `0.orig/U` and `system/controlDict` as
   changed, the two files the mesh job allows to differ. On HPC12 (GNU sed 4.2.2) the two
   patterns were run without `-i` inside a delivered solve case and printed the new lines;
   the in-place edit itself was tested locally. These lines cover a cold start only; a
   case that has already run is changed as in section 6.5.4.

3. Check that the three entries agree (section 6.4.4).

4. Make sure the job will really start cold. The HPC12 job script copies `0.orig` to `0`
   **only if `0/` does not exist** (`hpc/solve_hpc12.pbs`, lines 255 to 260), and it
   **restarts** instead of starting cold whenever `processor*` directories hold a complete
   time step (lines 179 to 218). A leftover `0/` or `processor*/` from an earlier angle
   would silently run the old angle.

   ```bash
   ls -d 0 processor* 2>/dev/null
   ```

   This must print nothing for a cold start. A solve case fresh from
   `hpc/stage_solve_case.sh` passes: the script never copies `0/`, `processor*` or `log/`,
   and says so in its output
   (`cold start: 0.orig holds 5 fields (U k nut omega p); no 0/, processor*, log/ or .converged`
   for a Condition CR case). If it prints `0`, delete that directory; the job script
   regenerates it from `0.orig`. Make sure you are in the case directory first, and never
   delete `0.orig`:

   ```bash
   cd $CASE
   rm -rf 0
   ```

   If it prints `processor...` directories, this is not a new case: either you meant the
   warm restart of section 6.5.4, or you copied a case that has already run.

   In the job's output, a correct cold start shows `COLD leg: no decomposed field to
   continue from` and `0.orig -> 0 (N fields)`. While the job runs, that output is in
   `$HOME/<job id>.hpc12.hpc.OU`, not yet in the case (section 6.5.4, step 11, shows how
   to read it). If you see `0 already exists, left alone` on a leg you meant to be cold,
   stop the job and check `0/U`.

### 6.4.3 Where the angle lives once a case has run

After a leg has run, the solver no longer reads `0.orig/U`. It restarts from the
**latest time directory of every `processor*` directory**, and there the farfield
velocity is written out explicitly in each boundary patch. Section 6.5.4 changes those
entries. The patches and entries that carry the freestream are:

| Condition | Velocity boundary type on the farfield | Patches | Entries to change |
|---|---|---|---|
| CR | `inletOutlet` | `inlet`, `topBottom`, `outboard` | `inletValue`, `value` |
| Early and late cruise | `freestreamVelocity` | `inlet`, `outlet`, `outboard`, `topBottom` | `freestreamValue`, `value` |

At Condition CR the `outlet` patch has `inletValue uniform (0 0 0)` and is left alone. No
other field file (p, k, omega, nut, T, nuTilda, alphat) needs editing: their farfield
values are scalars that do not change with the angle.

### 6.4.4 Check that the three entries agree

Run this in the case directory. It needs only a standard Python 3; on HPC12 the login
node's `python3` (3.6.8) runs it.

```bash
python3 - <<'EOF'
import math, re
def vec(path, key):
    txt = re.sub(r"//.*", "", open(path).read())
    m = re.search(r"\b%s\s*\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)" % key, txt)
    return [float(v) for v in m.groups()]
ux, uy, uz = vec("0.orig/U", "Uinf")
lx, ly, lz = vec("system/controlDict", "liftDir")
dx, dy, dz = vec("system/controlDict", "dragDir")
cd = re.sub(r"//.*", "", open("system/controlDict").read())
mag = float(re.search(r"\bmagUInf\s+(\S+?)\s*;", cd).group(1))
print("alpha from Uinf     %.9f deg" % math.degrees(math.atan2(uz, ux)))
print("alpha from dragDir  %.9f deg" % math.degrees(math.atan2(dz, dx)))
print("alpha from liftDir  %.9f deg" % math.degrees(math.atan2(-lx, lz)))
print("|Uinf| %.6f m/s, magUInf %.6f m/s" % (math.sqrt(ux*ux + uy*uy + uz*uz), mag))
EOF
```

Expected output after the example above:

```
alpha from Uinf     1.747608659 deg
alpha from dragDir  1.747608659 deg
alpha from liftDir  1.747608659 deg
|Uinf| 34.000000 m/s, magUInf 34.000000 m/s
```

The pass rule, the same tolerances the mesh job's Gate A applies
(`hpc/mesh_hpc12.pbs`):

1. The three angles agree to within 0.000001 deg (1e-6 deg).
2. `|Uinf|` is within 0.001 m/s of `magUInf`.

A case set with section 6.4.2 agrees exactly, as above. A case whose entries were written
with fewer digits shows small differences that are rounding, not a mistake. Run inside the
delivered B solve case on HPC12, whose entries were written by earlier tools with fewer
digits (its force axes read `liftDir ( -0.0304968 0 0.999535 )`), the check printed:

```
alpha from Uinf     1.747607675 deg
alpha from dragDir  1.747608659 deg
alpha from liftDir  1.747608659 deg
|Uinf| 34.000020 m/s, magUInf 34.000000 m/s
```

The angles differ by 0.000000984 deg and the speeds by 0.00002 m/s, both inside the rule.

This is what a mistake looks like. On a fresh copy of the Condition CR recipe, only
`Uinf` was changed (the two `foamDictionary` lines for the force axes were forgotten), and
the check printed:

```
alpha from Uinf     1.747608659 deg
alpha from dragDir  1.299999393 deg
alpha from liftDir  1.299999393 deg
|Uinf| 34.000000 m/s, magUInf 34.000000 m/s
```

The force axes still carry the recipe's 1.3 deg. The speeds agree, so the speed line does
not catch this mistake; only the angle lines do. A speed line where `|Uinf|` differs from
`magUInf` by more than 0.001 m/s means something else: the `Uinf` vector was computed with
the wrong `UMAG` (for example a cruise speed on a Condition CR case). In either case, fix
the entry that disagrees and run the check again.

---

## 6.5 The trim procedure

### 6.5.1 Starting angles

1. **Condition CR.** The recipe ships at 1.3000 deg. That is only a starting value, not the
   trim of any geometry: the delivered CR trims of B and C are 1.747608659 deg and
   1.268311286 deg (section 6.1).
2. **Early and late cruise.** The recipes ship at the VLM (vortex-lattice) trim angle,
   0.8654766754810572 deg (early) and 0.9966996379212477 deg (late). `CONDITION.json`
   labels this `STARTING GUESS from the VSPAERO/VLM trim, NOT a RANS trim`. It is low: the
   delivered RANS trims of B are 1.102080 deg (early cruise) and 1.263549 deg (late
   cruise) (`data/rans_forces.json`, welded mesh). Part of that gap may come from the
   unresolved frame. If the design study's cruise target is on the TP-1580 basis, its VLM
   trim angle is for less lift than the RANS target on the DSO basis, which alone would
   make it lower. Do not read the gap as a measure of VLM error.

Note that the recipe comment at the top of `recipes/L11_wallResolved_CR/0.orig/U` still
reads `|U| = 40.8 m/s` and `alpha = 2.1028 deg`. That comment belongs to the wind-tunnel
recipe. The values in the file are 34.0 m/s and 1.3000 deg; trust the values, not the
comment.

### 6.5.2 Three-point sweep and straight-line fit

This is how the delivered trims were found.

1. Make three solve cases from the same finished mesh, one per angle, and set a different
   angle in each. Pick three angles that you expect to straddle the trim. The delivered
   sweeps spaced their angles 0.10 to 0.33 deg apart, for example 1.60, 1.75 and 1.90 deg
   for B at Condition CR, and 0.8655, 1.20 and 1.35 deg for B at both cruise points.

   Each angle gets **its own solve case**, staged from the one mesh case `$MESH`. For
   **each** angle, in this order, at the HPC12 login prompt:

   1. Point `$CASE` at a new directory of its own that does not exist yet, for example:

      ```bash
      export CASE=/home/scratch/$USER/argus/solve_r2/<name>_a1p60
      ```

   2. Make the solve case there from the mesh case:

      ```bash
      $REPO/hpc/stage_solve_case.sh $MESH $CASE
      ```

      It prints `STAGED: <your solve case> (<size>)` and then a `NEXT STEPS` list; item 2
      shows the angle the new case carries, which is the mesh case's.
      [./05_submit_and_run.md](./05_submit_and_run.md), section 6, explains the rest of
      its output. On HPC12 only the script's input checks have been run; its copying,
      verifying and writing steps have run only on a local Linux machine: **NOT VERIFIED**
      on HPC12 (same section). If it prints `REFUSED: solve case '...' already exists.`,
      you forgot to change `$CASE`. Every call copies the whole mesh: the
      `constant/polyMesh` of the delivered mesh case SWB_r2 alone is 12,162,360 KiB (about
      11.6 GiB), so check your scratch space before you stage several angles.
   3. Set its angle (section 6.4.2) and check it (section 6.4.4). Check the post-processing
      block, which the staging script copied from the mesh case:

      ```bash
      python3 $REPO/scripts/assert_case_postpro.py --before $CASE
      ```

      It must end `---- 8 of 8 requirements satisfied ----`. Then submit the case with the
      `qsub` line the staging script printed (`NEXT STEPS`, item 4; see also section 6.5.4,
      step 10).

      **If it does not** end `8 of 8`: the usual cause is a mesh case built without the
      post-processing block, and then the staging script has already printed
      `WARNING: system/controlDict does NOT include argusPostPro`. Install the block on
      this solve case as in [./05_submit_and_run.md](./05_submit_and_run.md), section 7a,
      and repeat the check before you submit.

   This was tested locally with the shipped scripts: two solve cases staged from one test
   mesh case, set to 1.60 and 1.75 deg, each printed its own angle three times in section
   6.4.4 and each ended `8 of 8`.

   Export `$CASE` again before you work on the next angle, and check it with `echo $CASE`
   before every `qsub`. If `$CASE` is empty, `cd $CASE` quietly takes you to your home
   directory, and the job only refuses (`ARGUS_CASE not set`) after it has waited in the
   queue.
2. Run each case to its `.converged` marker ([./05_submit_and_run.md](./05_submit_and_run.md)
   and section 6.6).
3. Read each case's converged C_L from its marker. Use this line and nothing else; it is
   the mean over the convergence window, not a single noisy sample:

   ```bash
   grep "REPORT THESE" $CASE/.converged
   ```

   It prints a line of this form (B, Condition CR, wake-refined mesh):

   ```
     REPORT THESE: Cd 0.0188031  Cl 0.4282682  (window means)
   ```

   A marker written by this repository's `hpc/solve_hpc12.pbs` always has this line
   (lines 386 and 415 to 419). If `grep` prints nothing, the file is not a finished marker
   from this job script; do not use it.

   You may meet an older one-line marker format in the delivered data, for example
   `case SWB_a1p60: reached 4500 of endTime 4500; drift ... -> CONVERGED` in the
   `marker_text` field of some rows of `data/rans_forces.json` (12 rows have
   `marker_format` `legacy_single_line`). That format carries no window mean. For those
   rows the file's `Cl` and `Cd` values were recomputed from the force history as the mean
   of the last 200 samples (the row's `sources` field says so), so you can use them as they
   are.

   One exception: `SWB_trim` is one of the 12, but its `Cl` (0.4282682) and `Cd`
   (0.0188031) are the `REPORT THESE` window means of B's final wake-refined marker
   (section 6.5.4, step 3), and its `alpha_deg` (1.747608659) is that case's trim angle.
   Its `marker_text`, `endTime` (4000), window (`window_t_first` 3005 to `window_t_last`
   4000) and `sources` still describe an earlier solve of B on the published mesh. The
   same holds for the window, `marker_text` and `sources` fields of all six `SW*_trim`
   rows: their windows read 3005 to 4000, although the wake-refined cases ran to 6000
   (section 6.6.3). Use the `Cl`, `Cd` and `alpha_deg` values of these rows, and do not
   read their marker or window fields as a description of the run that produced them.
4. Fit a straight line C_L = m alpha + c through all three points by least squares, and
   invert it at the target. With B's published-mesh sweep at Condition CR (values from
   `data/rans_forces.json`, cases `SWB_a1p60`, `SWB_a1p75`, `SWB_a1p90`):

   ```bash
   python3 - <<'EOF'
   target = 0.428277635108                    # target C_L of this operating point (DSO basis)
   alpha  = [1.60, 1.75, 1.90]                # deg, one per bracket case
   cl     = [0.4160231, 0.4289662, 0.4418115] # window means, one per bracket case
   n  = len(alpha)
   ma = sum(alpha) / n
   mc = sum(cl) / n
   m  = sum((a - ma) * (c - mc) for a, c in zip(alpha, cl)) / sum((a - ma) ** 2 for a in alpha)
   c0 = mc - m * ma
   a_trim = (target - c0) / m
   print("slope dCL/dalpha  %.7f per deg" % m)
   print("fit residuals     " + "  ".join("%+.1e" % (c - (m * a + c0)) for a, c in zip(alpha, cl)))
   print("trim alpha        %.6f deg, %s the bracket" % (a_trim, "INSIDE" if min(alpha) <= a_trim <= max(alpha) else "OUTSIDE"))
   EOF
   ```

   Output:

   ```
   slope dCL/dalpha  0.0859613 per deg
   fit residuals     -1.6e-05  +3.3e-05  -1.6e-05
   trim alpha        1.742369 deg, INSIDE the bracket
   ```

   By hand: the mean angle is 1.75 deg and the mean C_L 0.4289336. The slope is
   0.0038683 / 0.045 = 0.0859613 per deg. The trim angle is
   1.75 + (0.428277635 - 0.4289336) / 0.0859613 = 1.75 - 0.0076309 = 1.742369 deg, which
   is B's published trim angle to six decimals. The residuals (at most 0.33 counts) say
   the three points lie on a straight line to well inside the tolerance. Large residuals
   mean the lift curve is bending and you should add a point.
5. **The trim angle must lie inside the bracket.** If the script prints `OUTSIDE`, the
   answer is an extrapolation: run a fourth angle beyond it and fit again. Three delivered
   geometries needed this: C at Condition CR (an extra case at 1.20 deg), and C and M at
   early cruise (an extra case at 0.65 deg).
6. Stage one more solve case from the same mesh case (step 1), set the fitted angle in it
   (section 6.4.2), run it to its marker, and read its C_L (step 3). That solved case is
   the result you quote, not the interpolation.
7. If `|target - C_L| < 1e-4`, the trim is done. If not, continue from this case with the
   two-leg pattern below rather than starting yet another cold case.

For the record: the delivered early- and late-cruise trims were not reached this way. They
were continued from one case with the two-leg pattern below, and all twelve finished inside
the tolerance, the widest `LCW_trim` at 0.895 counts. `data/rans_forces.json` lists every
residual under `_trim_quality`, and every leg under each cruise trim's `trim_history`.

### 6.5.3 The two-leg pattern

This is how the delivered trims were done: Condition CR on the wake-refined mesh, and both
cruise points, where three of the twelve cases needed a second nudge and so a third leg.
Instead of
building a new cold case for every angle, one case is continued from its own converged
flow:

1. **Leg 1 (cold):** run the case at a first angle alpha0 until it writes `.converged`.
2. Read its C_L (section 6.5.2 step 3). If `|target - C_L| < 1e-4`, stop: it is trimmed.
3. Otherwise compute a new angle, alpha1 = alpha0 + (target - C_L) / slope, using the slope
   from a sweep of the same geometry **at the same operating point**, preferably on the
   same mesh generation. The lift slope is not the same at Condition CR and at cruise, so
   a slope from one condition mis-sizes the step at another. (The delivered Condition CR
   trims on the wake-refined mesh sized their first step with the published-mesh slope;
   see the worked example below.)
4. **Leg 2 (warm):** change the angle in the finished case and continue from its converged
   state (section 6.5.4).
5. When leg 2 has its marker, check C_L again. If it is still outside the tolerance, repeat
   with the geometry's **own** slope, the secant through its two converged points.

**Worked example, from `data/run_cards/SWB_r2_trim.json`, field `trim.history`**
(B, Condition CR, wake-refined mesh):

| Step | alpha (deg) | C_L | target - C_L |
|---|---|---|---|
| leg 1, cold, at the published trim angle | 1.742369 | 0.427826554 | +4.511e-04 (4.5 counts, outside) |
| leg 2, warm, after one step | 1.747608659 | 0.428270799 | +6.836e-06 (0.07 counts, inside) |

The step after leg 1:

```bash
python3 - <<'EOF'
target = 0.428277635108             # target C_L (DSO basis)
alpha1, cl1 = 1.742369, 0.427826554 # leg 1: angle (deg) and converged C_L
slope = 0.0859613                   # dCL/dalpha per deg, from the sweep fit in 6.5.2
print("dCL         %+.3e" % (target - cl1))
print("step        %+.7f deg" % ((target - cl1) / slope))
print("next alpha  %.9f deg" % (alpha1 + (target - cl1) / slope))
EOF
```

```
dCL         +4.511e-04
step        +0.0052475 deg
next alpha  1.747616491 deg
```

The project's step was +0.0052397 deg (to 1.747608659 deg), which corresponds to a slope
of 0.08609 per deg, 0.15 percent from this fit. The difference is 0.0000078 deg, worth
under 0.01 counts of C_L. Leg 2 then landed 0.07 counts from the target, so no third step
was needed. Had one been needed, it would use B's own secant slope on the new mesh:

```bash
python3 - <<'EOF'
target = 0.428277635108
a1, c1 = 1.742369,    0.427826554   # leg 1
a2, c2 = 1.747608659, 0.428270799   # leg 2
s = (c2 - c1) / (a2 - a1)
print("dCL after leg 2   %+.3e" % (target - c2))
print("own secant slope  %.7f per deg" % s)
print("next alpha        %.9f deg" % (a2 + (target - c2) / s))
EOF
```

```
dCL after leg 2   +6.836e-06
own secant slope  0.0847851 per deg
next alpha        1.747689288 deg
```

Four things to know about this example and the delivered cases:

1. The C_L values in `trim.history` are not identical to the marker's window mean for the
   same leg: the marker of `SWB_trim` reports `Cl 0.4282682` against 0.428270799 in the
   card. Both are within the tolerance. For your own trims use the marker's
   `REPORT THESE` value.
2. Do not measure the slope from the first few hundred iterations after a step. Right
   after the angle changes, the lift overshoots and then relaxes as the wake re-forms. On
   B, C_L rose from 0.4278248 at iteration 4000 to a peak of 0.4283555 at iteration 4165
   and then fell at every later sample. At iteration 169 of that restart the apparent
   slope was 0.1008 per deg, 17 percent above the sweep slope. Use converged window means
   only.
3. An **angle step** is not the same thing as a **job leg**. B took two angles but three
   jobs, because its warm continuation was first submitted too short (section 6.6.3).
4. **The delivered Condition CR cases were set up the way this repository now sets up a
   case.** Their meshes used the recipe's 192-rank `scotch` decomposition. Their solve
   cases ran on 128 ranks, `hierarchical` (8 2 8), which is what `hpc/stage_solve_case.sh`
   writes by default. Their first leg ended at `endTime 4000` and they saved fields every
   500 iterations, as `recipes/L11_wallResolved_CR/system/controlDict` now does
   (`endTime 4000`, `writeInterval 500`); the staging script changes neither. So the
   leg boundary and save points quoted from the delivered cases on this page (a first leg
   ending at 4000, fields saved every 500 iterations) are the ones your own case will
   show. Two things differ. The timing of the restart transient is not set by the case:
   B's peak at iteration 4165 (point 2) depends on the geometry and on the size of the
   angle step. And the warm legs differ: all six delivered cases finished at 6000
   (section 6.6.3), while the rule in section 6.5.4, step 9, gives 5500 for a first warm
   leg.

### 6.5.4 Changing the angle of a finished case on HPC12, step by step

Do this on HPC12, in the case directory on scratch, when leg 1 has finished. The case must
still have its `processor*` directories. HPC12 scratch deletes files older than 50 days
(noted in `hpc/solve_hpc12.pbs`), so do not leave a case waiting that long between legs.

The commands are light: each `foamDictionary` call reads and rewrites one processor's
velocity file (reading that entry from ten processor files of a delivered 128-rank case
took under 2 seconds at the HPC12 login prompt). Run steps 1 to 10 in **one** terminal:
step 1 loads OpenFOAM-12, step 3 sets `CASE`, and the later steps share the shell
variables `CASE`, `UX`, `UZ`, `LX`, `LZ`, `DX`, `DZ`, `LT`, `WI` and `NEW`. Every
`foamDictionary` call below runs inside `$CASE` with a path relative to it (section
6.4.1).

What was checked where:

1. **On HPC12, at the login prompt:** the loading lines of step 1, and every read in this
   procedure (the marker in step 3, the time-step listing in step 6, the `foamDictionary`
   reads of steps 8 and 9) inside a delivered 128-rank solve case.
2. **On a local OpenFOAM-12 installation:** the whole procedure, including the steps that
   change files (5, 7 and 9). The test case was made from the Condition CR recipe with the
   post-processing block installed, staged with `hpc/stage_solve_case.sh`, and its legs
   were run through the unchanged `hpc/solve_hpc12.pbs`. It had 4 ranks; the expected
   output below is written for a 128-rank case.

The procedure:

1. Load OpenFOAM-12 at the login prompt. Type the three lines as they are, one per line,
   not through a pipe (a pipe runs them in a separate shell, and OpenFOAM is not loaded
   afterwards):

   ```bash
   module load devtoolset/11
   module load mpi/openmpi-4.1.2
   source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc
   command -v foamDictionary
   foamDictionary -help | grep -E "Using|writePrecision"
   ```

   **What you should see** (HPC12; your real home folder where this shows `$HOME`):

   ```
   dirname: missing operand
   Try 'dirname --help' for more information.
   $HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/foamDictionary
     -writePrecision <label>
   Using: OpenFOAM-12 (see https://openfoam.org)
   ```

   The `dirname` pair comes from the `source` line and is harmless: every HPC12 job log
   has it too. Before these lines, `command -v foamDictionary` prints nothing.

   **If it does not:** `command -v` prints nothing, or `Using:` names another version.
   See [./01_setup.md](./01_setup.md), section B4. If your OpenFOAM-12 is somewhere else,
   source that installation's `etc/bashrc` instead (the file you give the jobs as
   `ARGUS_OF_BASHRC`).

2. Make sure no job is running or queued for this case:

   ```bash
   qstat -u $USER
   qstat -f <job id> | grep Job_Name
   ```

   `qstat -u` shortens job names to 16 characters, so if you have several jobs, confirm the
   full name with `qstat -f` for each job id. Nothing may run in this case while you edit
   it.

3. Point `CASE` at the solve case (this is a new terminal, so it is not set yet), check
   it, and read the result of leg 1:

   ```bash
   export CASE=/home/scratch/$USER/argus/solve_r2/<name>
   ls $CASE/system/controlDict
   cd $CASE
   cat .converged
   ```

   `ls` must print the path of the case's `system/controlDict`. If it prints
   `No such file or directory`, `CASE` is wrong: fix the `export` line and run the four
   lines again. Do not go on while `CASE` is wrong or empty: an empty `CASE` makes
   `cd $CASE` take you quietly to your home directory.

   Every finished leg writes the same layout (`hpc/solve_hpc12.pbs`, lines 415 to 419).
   This example is the final marker of the delivered B trim (Condition CR, wake-refined
   mesh), which finished in its third job; yours will name your case and your leg:

   ```
   case SWB_trim
   finalised 2026-09-21T18:16:17Z on n12-142, leg 3
     CONVERGENCE GATE over 200 samples: MEAN drift Cd 0.0003 ct (<=0.05), Cl 0.0399 ct (<=1.0), Cl SPAN 0.0777 ct (<=1.0) -> CONVERGED
     REPORT THESE: Cd 0.0188031  Cl 0.4282682  (window means)
   ```

   If there is no `.converged` file, leg 1 is not finished. Do not continue (section 6.6.2).

4. Compute the new angle from the `Cl` value (section 6.5.3), then set `ALPHA`, `UMAG` and the
   three vectors exactly as in section 6.4.2 step 1.

5. **Rename the marker. Never delete it.** The job script refuses to run a case that has a
   `.converged` file (`hpc/solve_hpc12.pbs`, lines 84 to 89). It prints
   `REFUSING to restart ...: .converged present` and ends the job with a success status,
   so a forgotten marker costs you a queue wait and nothing else tells you.

   ```bash
   mv .converged .converged.alpha0_untrimmed
   ```

   The delivered B case kept its first-leg marker exactly this way, as
   `.converged.alpha0_untrimmed`.

6. Find the latest time step and check that every processor agrees on it:

   ```bash
   for p in processor*; do ls -d $p/[0-9]* | sed 's#.*/##' | sort -g | tail -1; done | sort | uniq -c
   LT=$(ls -d processor0/[0-9]* | sed 's#.*/##' | sort -g | tail -1); echo "latest time: $LT"
   ```

   The first command must print **one** line: the number of processors and one time. For a
   128-rank Condition CR case built from the recipe, after a first leg that ran to the
   recipe's `endTime` it prints `128 4000` (a cruise case: `128 10000`). On HPC12 the
   delivered B case, after its last leg, prints `128 6000`. More than one line means the
   ranks disagree; stop and read [./05_submit_and_run.md](./05_submit_and_run.md) on
   restarts before touching anything.

7. Write the new freestream into every processor's latest velocity file. Run the loop from
   `$CASE`: `$p/$LT/U` is then a path relative to the case, such as `processor0/4000/U`,
   as `foamDictionary` needs (section 6.4.1).

   **Condition CR:**

   ```bash
   V="uniform ($UX 0 $UZ)"
   for p in processor*; do
     foamDictionary $p/$LT/U -writePrecision 12 -set \
       "boundaryField/inlet/inletValue=$V, boundaryField/inlet/value=$V, boundaryField/topBottom/inletValue=$V, boundaryField/topBottom/value=$V, boundaryField/outboard/inletValue=$V, boundaryField/outboard/value=$V" \
       || echo "FAILED on $p"
   done
   ```

   **Early or late cruise:**

   ```bash
   V="uniform ($UX 0 $UZ)"
   for p in processor*; do
     foamDictionary $p/$LT/U -writePrecision 12 -set \
       "boundaryField/inlet/freestreamValue=$V, boundaryField/inlet/value=$V, boundaryField/outlet/freestreamValue=$V, boundaryField/outlet/value=$V, boundaryField/outboard/freestreamValue=$V, boundaryField/outboard/value=$V, boundaryField/topBottom/freestreamValue=$V, boundaryField/topBottom/value=$V" \
       || echo "FAILED on $p"
   done
   ```

   These loops print nothing when they succeed. Any `FAILED on processorN` line means that
   file was not changed; do not submit until it is fixed. A processor whose part of a
   patch has no faces carries `nonuniform List<vector> 0()` there; setting a uniform value
   on it is harmless (tested on a 12-rank split where 4 processors had no `topBottom`
   faces: the solver then ran without any new warning). The solved interior flow is not
   touched: in the tests the internal field read back byte for byte identical before and
   after (checked on all 4 processor files of the Condition CR test case, and earlier on
   one processor file of an early-cruise test case).

8. Confirm every processor now carries the new value (use `freestreamValue` instead of
   `inletValue` at cruise):

   ```bash
   for p in processor*; do foamDictionary $p/$LT/U -entry boundaryField/topBottom/inletValue -value -writePrecision 12; done | sort | uniq -c
   ```

   Expected: one line, the processor count and the new vector, for example
   `128 uniform ( 33.9841853945 0 1.03689105977 )`. Two lines means some files still carry
   the old angle. (Run before step 7, the same command lists the old vector, and on a
   case whose split leaves some processors without `topBottom` faces it also lists
   `nonuniform List<vector> 0()` for those. After step 7 they all carry the new vector.)

9. Update the record copy and the force axes, and raise `endTime`.

   The new `endTime` must satisfy three rules:

   1. **It must be above the latest time `$LT`.** This is the rule that matters most. If
      it is not, the leg runs zero iterations, the job script counts the leg as having
      reached `endTime`, and it runs the convergence gate on the force history of the
      **old** angle. If that history passes, the script writes a `.converged` marker whose
      `REPORT THESE` C_L belongs to the old angle, while `0.orig/U` and
      `system/controlDict` show the new one. Nothing in the output warns you except the
      line `0 iterations: already at endTime`. The job script itself only checks rule 3.
   2. **It must leave room for the flow to settle after the angle step.** At Condition CR,
      at least 1500 iterations beyond `$LT` (section 6.6.3 explains why). At early or late
      cruise, at least 2000 beyond `$LT`, which is one `writeInterval` of the cruise
      recipes and the smallest step they allow; no warm cruise leg has been run to
      measure how much is really needed (section 6.6.3).
   3. **It must be a whole multiple of `writeInterval`.** The job script refuses anything
      else with `endTime X is not a multiple of writeInterval Y` (`hpc/solve_hpc12.pbs`,
      lines 237 to 243).

   First write the new angle into the record copy and the force axes, and compute a new
   `endTime` from `$LT` (step 6) and the case's own `writeInterval`. Set `MIN` for your
   operating point:

   ```bash
   MIN=1500    # Condition CR; use MIN=2000 at early or late cruise
   foamDictionary 0.orig/U -entry Uinf -set "($UX 0 $UZ)" -writePrecision 12
   foamDictionary system/controlDict -entry functions/forceCoeffs1/liftDir -set "($LX 0 $LZ)" -writePrecision 12
   foamDictionary system/controlDict -entry functions/forceCoeffs1/dragDir -set "($DX 0 $DZ)" -writePrecision 12
   WI=$(foamDictionary system/controlDict -entry writeInterval -value)
   NEW=$(( (LT + MIN + WI - 1) / WI * WI ))
   echo "latest time $LT, endTime now $(foamDictionary system/controlDict -entry endTime -value), writeInterval $WI, new endTime $NEW"
   ```

   For a Condition CR case built from the recipe, after its first leg, the last line reads
   `latest time 4000, endTime now 4000, writeInterval 500, new endTime 5500`. For a cruise
   case it reads `latest time 10000, endTime now 10000, writeInterval 2000, new endTime 12000`.
   **Check it before going on.** If `latest time` is empty, or `new endTime` is not larger
   than `latest time`, then `LT` or `MIN` is not set in this terminal: repeat step 6 and
   this step. If `writeInterval` is empty or you see `division by 0`, `foamDictionary`
   could not read `system/controlDict`: check that you are in `$CASE` and that step 1
   worked.

   Only when the line is right, write the new `endTime`:

   ```bash
   foamDictionary system/controlDict -entry endTime -set $NEW -writePrecision 12
   ```

   It prints `New entry endTime` followed by your value, for example
   `New entry endTime         5500;`.

   Then run the check of section 6.4.4. All three angles must equal your new `ALPHA`. The
   post-processing block survives these edits: in the test,
   `python3 $REPO/scripts/assert_case_postpro.py --before $CASE` still ended
   `---- 8 of 8 requirements satisfied ----` afterwards, which matters because the next
   leg runs that check before it starts (step 10).

10. Resubmit the case from inside `$CASE` with the same `qsub` line you used for leg 1,
    which `hpc/stage_solve_case.sh` printed (`NEXT STEPS`, item 4). For a 128-rank case on
    `fpt-medium`:

    ```bash
    cd $CASE
    qsub -N argus-<name> -q fpt-medium -l nodes=4:ppn=32:typej,walltime=72:00:00 \
         -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py \
         $REPO/hpc/solve_hpc12.pbs
    ```

    The faster alternative, if you have no other `fpt-large` job (that queue holds one job
    per user), replaces the `-q` and `-l` part with
    `-q fpt-large -l nodes=1:ppn=128:typen,walltime=72:00:00`. `qsub` prints the job id, a
    number followed by the server name; keep the number, it is the `<job id>` below. PBS
    writes the job's output file, `argus-<name>.o<job id>`, in the folder `qsub` runs
    from, but only when the job ends. While the job runs, the same text is in
    `$HOME/<job id>.hpc12.hpc.OU` ([./05_submit_and_run.md](./05_submit_and_run.md),
    section 11b).

    1. Use the `-q` and `-l` that match your case's rank count; the script refuses a job
       whose allocation has fewer cores than the case has ranks (`hpc/solve_hpc12.pbs`,
       lines 121 and 122).
    2. Do not pass `ARGUS_ENDTIME` here: the new `endTime` is already in
       `system/controlDict` from step 9.
    3. Always pass `ARGUS_SCRIPT`: without it a leg that needs to continue ends with
       `leg N needs to requeue but ARGUS_SCRIPT is unset; refusing to guess`, which is what
       happened to B's second job.
    4. Always pass `ARGUS_GATE`: the job then runs the post-processing check at the start
       of the leg, before the write probe, `decomposePar` and `foamRun`, and stops a leg
       that fails it with `REFUSED: post-processing gate refused`. Such a leg leaves only
       an empty `log/leg<N>` directory in the case. Without `ARGUS_GATE` the job looks for
       `$HOME/argus_hpc12/assert_case_postpro.py`, which does not exist on a new account,
       and runs the leg unchecked. `$REPO` must be set, as an absolute path, in the shell
       you submit from. A real PBS job has received `ARGUS_GATE` through `-v` on HPC12 and
       ran the check (2026-09-24, `8 of 8 requirements satisfied`).
    5. If you passed `ARGUS_OF_BASHRC` on leg 1, add it here as well. Every leg is
       submitted by hand on HPC12, and nothing carries these settings over
       ([./05_submit_and_run.md](./05_submit_and_run.md), section 10a).

    `hpc/README_HPC.md` and `hpc/submit.sh` describe an earlier DelftBlue/SLURM setup and
    are not the HPC12 route.

11. Once the job has started, check its output. While the job runs (`qstat -u $USER`
    still lists it), the output is only in the live copy in your home folder. The output
    file in `$CASE` appears when the job ends, and the live copy is then removed:

    ```bash
    grep -E "requirements satisfied|gate NOT run|REFUS|RESTART leg|COLD leg|left alone|endTime:|already at endTime" $HOME/<job id>.hpc12.hpc.OU
    grep -E "requirements satisfied|gate NOT run|REFUS|RESTART leg|COLD leg|left alone|endTime:|already at endTime" $CASE/argus-<name>.o<job id>
    ```

    Use the first line while the job runs and the second once it has ended. Check early:
    the lines below are printed before the solver starts (21 seconds after the start of
    the delivered B case's second leg, HPC12 job logs), so a wrong start can be stopped
    before it costs a leg.

    **What you should see** (a Condition CR case built from the recipe, after step 9; the
    lines the job script writes itself start with the time of day):

    ```
        ---- 8 of 8 requirements satisfied ----
    [01:38:09] RESTART leg: resuming from t=4000
    [01:38:09] endTime: 5500 (from the case)
    ```

    The delivered B case's second job printed the same two lines with
    `resuming from t=4000` and `endTime: 5000 (from the case)` (HPC12 job logs). The
    `resuming from t=` value must be **smaller** than the `endTime:` value. If it is not,
    or if the output later says `0 iterations: already at endTime`, `endTime` was not
    raised (step 9): any `.converged` marker this leg writes reports the **old** angle.
    Stop the job if it is still running (section 6.7). If it has already ended and
    written a marker, rename that marker so that it is kept but no longer counts. Then
    redo step 9 and resubmit:

    ```bash
    cd $CASE
    mv .converged .converged.not_the_new_angle
    ```

    **If it does not:**

    1. `COLD leg: no decomposed field to continue from`, then `0 already exists, left alone`:
       the processor directories were not found. The job does **not** go back to
       `0.orig`. It starts cold from the case's own `0/`, which the first leg created and
       which still holds the first leg's angle, while the force axes in
       `system/controlDict` carry the new one. (The delivered B case's `0/U` still reads
       1.742369 deg, its first-leg angle.) Stop the job (section 6.7), then stage a new
       solve case from the mesh case and set the new angle there as a cold start (section
       6.5.2, step 1).
    2. `REFUSING to restart <case>: .converged present`: you did not rename the marker
       (step 5).
    3. `---- N of 8 requirements satisfied ----` with N below 8, then
       `REFUSED: post-processing gate refused`: the job stopped before the write probe,
       `decomposePar` and `foamRun`. It leaves only an empty `log/leg<N>` directory, so
       the next job's leg number goes up by one
       ([./05_submit_and_run.md](./05_submit_and_run.md), section 10a, point 5). See
       which requirement is missing:

       ```bash
       python3 $REPO/scripts/assert_case_postpro.py --before $CASE
       ```

       Fix it ([./05_submit_and_run.md](./05_submit_and_run.md), section 7a), and
       resubmit.
    4. `NOTE: <path> absent, post-processing gate NOT run`: nothing exists at the path in
       `ARGUS_GATE` (unset, misspelled or relative, or `$REPO` was not set). The leg runs
       unchecked. Fix the `-v` list for the next leg.

    All of these messages were produced locally by the unchanged job script, and the
    `grep` lines above print each of them.

12. When the new marker appears, go back to section 6.5.3 step 5.

**The warm route at the cruise conditions has not been used for any delivered result.**
The delivered early- and late-cruise trims were cold cases only (section 6.5.2). The cruise
commands in step 7 were checked on a small test case built from the early-cruise recipe,
which restarted cleanly and printed the new `dragDir` in its force file, but not on a
production cruise mesh. The step 7 cruise loop and the step 8 check (with
`freestreamValue`) were run again on a decomposed early-cruise test case with the
post-processing block installed: one line, `4 uniform ( 233.891216083 0 4.47532073066 )`
for 4 ranks at 1.096175 deg.

---

## 6.6 When is a run finished: the convergence gate

### 6.6.1 What the gate tests

The gate is implemented in `hpc/solve_hpc12.pbs`, in the Python block on lines 335 to 391.
It runs **once, at the end of each leg**, on the force history in
`postProcessing/forceCoeffs1/*/forceCoeffs.dat` (all restarts merged and sorted by time).
Thresholds, as written in the code:

1. **Window:** the last 200 force samples (`W = 200`, line 360).
2. **Mean drift of C_D:** the size (absolute value) of the difference between the mean of
   the second 100 samples and the mean of the first 100, times 10,000, must be **at most
   0.05 counts** (lines 369 to 371 and 382). A drift of -0.08 counts fails just as +0.08
   does.
3. **Mean drift of C_L:** the same absolute difference for C_L, **at most 1.0 count**
   (lines 372 and 382).
4. **C_L span:** the size (absolute value) of the difference between the last C_L sample
   and the first in the window, times 10,000, **at most 1.0 count** (lines 381 and 382).
   This catches a lift curve that turns round inside the window, which the half-means
   cannot see.

All three must pass. The gate then prints the window means (`REPORT THESE`). Those are the
numbers to quote, never the last sample: a converged steady solve keeps oscillating about
a steady mean.

**200 samples is not the same number of iterations at every condition.** The Condition CR
recipe writes forces every 5th iteration, so its window is the last **1000 iterations**.
The cruise recipes write every iteration, so their window is the last **200 iterations**.
(A delivered 4000-iteration CR case has 801 force samples and a cruise case 4001, in
`data/rans_forces.json`.)

The `.converged` marker is written **only** if the leg reached its `endTime` **and** the
gate passed (`hpc/solve_hpc12.pbs`, lines 410 to 420). If the leg reached `endTime` but
the gate failed, the script tries to continue the case for another 1500 iterations
(lines 421 to 424) by submitting a follow-on job from inside the job. On HPC12 that
submission does not get through (section 6.7, point 4), so you resubmit the next leg
yourself with `ARGUS_ENDTIME` raised, as described in
[./05_submit_and_run.md](./05_submit_and_run.md), section 10a. The job's output names the
value it asked for, in the line `requeueing as leg N with endTime X`. The value you pass
must be a multiple of `writeInterval`, or the job script refuses the leg (lines 237 to
243):

1. **Condition CR:** the value the script asks for is valid. 4000 + 1500 = 5500, and
   later 7000, are multiples of the recipe's `writeInterval` 500.
2. **Early and late cruise:** it is not. 10000 + 1500 = 11500 is not a multiple of the
   recipes' `writeInterval` 2000, and that leg would be refused. Pass the next multiple of
   2000 instead, for example `ARGUS_ENDTIME=12000`.

### 6.6.2 A live CONVERGED reading is not the `.converged` marker

Monitoring tools that apply the same test to the force history **while a job is still
running** (the campaign had one; it is not in this repository) can read `CONVERGED`
hundreds of iterations before the job ends. Measured in the campaign: one late-cruise case
read `CONVERGED` at iteration 2088 of 3000 and then ran 912 more iterations, during which
its window means were still free to move.

Rules:

1. **Only the `.converged` file means finished.** A `CONVERGED` reading on a running case
   is a forecast.
2. Never compute a trim, fit a slope or change an angle from a case that has no
   `.converged` file.
3. Count finished cases by their markers, for example:

   ```bash
   find /home/scratch/$USER/argus -maxdepth 3 -name .converged
   ```

   Adjust the path to where your cases live. Each line is one finished case.

### 6.6.3 A leg that starts with an angle step

The gate looks at the last 1000 iterations at Condition CR. If a warm leg (section 6.5.4)
is only 1000 iterations long, the window contains the deliberate angle step itself, and the
C_L span test fails even though the case is fine. This happened on B: its first warm job
ran 4000 to 5000 and the gate reported `Cl SPAN 4.3254 ct (<=1.0)`, `NOT CONVERGED`. The
windows that were tested starting 200, 400 and 700 iterations after the step all passed.
The fix is **not** to pick a window that passes; it is to run the warm leg long enough
that the gate's own window starts after the transient: **at least 1500 iterations at
Condition CR**. A warm leg set up by section 6.5.4, step 9, runs from 4000 to 5500, so its
window starts 500 iterations after the step. The delivered CR trims all ran to 6000, from
first legs that ended at 4000 (HPC12 job logs): C, F, H and W in one warm leg of 2000
iterations, B in two (to 5000, then 6000), and M to 5500 and then 6000, because its
`endTime` was lowered while it ran (section 6.7).

At the cruise conditions the window is only 200 iterations, so it is less likely to catch
the step, but the flow must still have settled after the step before the window starts. No
warm cruise leg has been run to measure how long that takes.

---

## 6.7 Warning: never change `endTime` in a running case

> [!WARNING]
> **If you must stop a leg earlier than planned, do not lower `endTime` in the running
> case. Stop the job and resubmit it with the new `endTime` instead.**

The recipes set `runTimeModifiable true`, so the solver re-reads `system/controlDict`
while it runs. Lowering `endTime` in the file therefore does make the solver stop at the
new value, cleanly (checked under OpenFOAM-12: a test run started with `endTime 3000` and
lowered to 95 while running stopped at 95 and exited normally). The problem is the job
script around it (`hpc/solve_hpc12.pbs`):

1. **It reads `endTime` once, when the job starts.** Lines 222 to 224 read `endTime`
   from `system/controlDict` (or take `ARGUS_ENDTIME` if you passed it) into the shell
   variable `ENDTIME`, and line 245 writes that value back into `system/controlDict`.
   After line 246 the script never reads `endTime` from `system/controlDict` again, so a
   later edit of `endTime` in the file does not change `ENDTIME`. `ENDTIME` is what the
   script compares against after the solve (line 330) and what it passes to its requeue
   (line 427).
2. **After the solver stops, it compares against the old value.** Lines 328 to 331 take
   the latest time directory of `processor0` and set `reached=1` only if it is at or beyond
   the remembered `ENDTIME`. A run you stopped early at a lower `endTime` is below it, so
   `reached=0`.
3. **So it treats your deliberate stop as a walltime kill.** With `reached=0` the script
   skips the branch that writes `.converged` (lines 410 to 420), **even if the gate
   passed**, prints `-> did not reach endTime (walltime). Continuing at the same endTime.`,
   and calls its requeue with the **old** `endTime` (lines 425 to 428).
4. **The requeue carries no queue or resource options, and the old `endTime`.** Its
   `qsub` (lines 404 to 406) passes only `-N` and `-v`, so if it were accepted the next
   leg would get the script header's `fpt-medium`, `nodes=4:ppn=32:typej`, 72 hours,
   whatever you had submitted with, and it would write the old, higher `endTime` back
   into the case and run past the point where you wanted to stop. In practice the
   campaign recorded that HPC12 compute nodes cannot reach the PBS server
   (`could not connect to trqauthd`), so this `qsub` from inside the job fails and nothing
   is queued (see [./05_submit_and_run.md](./05_submit_and_run.md)). If `ARGUS_SCRIPT` was
   not passed, the job ends with the `ARGUS_SCRIPT is unset` message instead. In every one
   of these outcomes the case is left **without** a `.converged` marker, looking
   unfinished.
5. **You may also lose iterations.** OpenFOAM-12 writes a time step only at multiples of
   `writeInterval`; it does not force a write when it stops at `endTime`. If the lowered
   value is not a multiple of `writeInterval`, everything after the last write is lost
   (checked: a test run stopped at `endTime 17` with `writeInterval 5` left 15 as its
   latest time step).

The project found this trap on 2026-09-23. The job script shipped in `hpc/` still behaves
this way.

**What to do instead:**

1. Do not edit `system/controlDict` of a case while its job is running, for any reason.
   The `endTime` a leg runs to is fixed when the job starts: `ARGUS_ENDTIME` if you passed
   it, otherwise whatever `system/controlDict` says at that moment.
2. To stop earlier than planned, stop the job and resubmit it with a new `endTime`:

   1. Stop the job, then list your jobs until it no longer appears:

      ```bash
      qdel <job id>
      qstat -u $USER
      ```

      `<job id>` is the job's number as `qstat -u $USER` lists it. `qdel` is how the
      project's own HPC12 tooling stopped jobs.

   2. Find the last iteration the stopped job reached. `<N>` is the leg number of the job
      you stopped (the `log/leg<N>` directory it wrote):

      ```bash
      grep "^Time = " $CASE/log/leg<N>/foamRun.log | tail -1
      ```

      OpenFOAM-12 prints this line as `Time = ` followed by the iteration and the letter
      `s`, for example `Time = 6000s` (the last line of this kind in the delivered B
      case's third leg, HPC12 job logs).

   3. Choose the new `endTime`: the first multiple of `writeInterval` **at or above** that
      iteration. Example: a Condition CR case built from the recipe (`writeInterval 500`)
      whose leg was submitted to 6000 and stopped at iteration 5437 gets 5500
      (5437 / 500 = 10.9, and 11 x 500 = 5500). The new leg restarts from the last saved
      time step (5000 here), runs to 5500, and its force samples replace the ones the
      stopped job wrote after 5000, so the gate's window ends at a state that is saved on
      disk. Do **not** choose the last saved time step or anything below it (5000 in the
      example): that leg runs zero iterations (`0 iterations: already at endTime`), and the
      gate window then still includes the stopped job's force samples from after the last
      saved time step. Both behaviours were checked in a local OpenFOAM-12 test. At early
      or late cruise the multiples are of 2000.

   4. Resubmit with that value passed explicitly:

      ```bash
      cd $CASE
      qsub -N argus-<name> -q fpt-medium -l nodes=4:ppn=32:typej,walltime=72:00:00 \
           -v ARGUS_CASE=$CASE,ARGUS_SCRIPT=$REPO/hpc/solve_hpc12.pbs,ARGUS_GATE=$REPO/scripts/assert_case_postpro.py,ARGUS_ENDTIME=<new endTime> \
           $REPO/hpc/solve_hpc12.pbs
      ```

      This `qsub` form with `ARGUS_ENDTIME` is the one the job script itself uses
      (line 405). The job's `endTime:` line must show your new value
      (`endTime: case says <old>, ARGUS_ENDTIME overrides to <new>`). A value that is not
      a multiple of `writeInterval` is refused (lines 237 to 243). Keep `ARGUS_GATE` as in
      section 6.5.4, step 10, and add `ARGUS_OF_BASHRC` if you used it on leg 1: every leg
      is submitted by hand on HPC12, and nothing carries these settings over
      ([./05_submit_and_run.md](./05_submit_and_run.md), section 10a).
3. To run longer: let the leg finish. If it has a marker, rename the marker (section
   6.5.4 step 5) and resubmit with a higher `endTime`.

---

## 6.8 Troubleshooting

| What you see | What it means | What to do |
|---|---|---|
| `REFUSING to restart <case>: .converged present` | The marker from the last leg is still there | Rename the marker (section 6.5.4 step 5), then resubmit |
| `COLD leg: no decomposed field to continue from` on a warm leg | No `processor*` directories. The job starts from the case's old `0/`, at the first leg's angle (`0 already exists, left alone`) | The warm route needs them. Stage a new solve case from the mesh case (section 6.5.2, step 1) and set the new angle there as a cold start |
| `0 already exists, left alone` on a cold leg | A stale `0/` may carry an old angle | Stop the job, check `0/U`, delete `0/` if in doubt |
| `endTime X is not a multiple of writeInterval Y` | New `endTime` does not fit the write cadence. At cruise this is what the job's own raise to 11500 gives (section 6.6.1) | Pick a multiple of `writeInterval` (500 at Condition CR, 2000 at cruise) that is also above the latest time (section 6.5.4 step 9) |
| `REFUSED: solve case '...' already exists.` from `hpc/stage_solve_case.sh` | `$CASE` still points at an earlier angle's solve case | Export a new `$CASE` for this angle (section 6.5.2, step 1); the script never overwrites a case |
| `---- N of 8 requirements satisfied ----`, then `REFUSED: post-processing gate refused` | The case does not carry the post-processing block, or its include line in `system/controlDict` is gone | Run `python3 $REPO/scripts/assert_case_postpro.py --before $CASE`, fix what it lists as `MISSING` ([./05_submit_and_run.md](./05_submit_and_run.md), section 7a), resubmit |
| `NOTE: <path> absent, post-processing gate NOT run` | `ARGUS_GATE` is missing from the `-v` list, misspelled or relative | The leg ran unchecked. Put `ARGUS_GATE=$REPO/scripts/assert_case_postpro.py` in the next leg's `-v` list (section 6.5.4 step 10) |
| `FOAM FATAL IO ERROR` with `file "<current directory>//<path>" does not exist` from `foamDictionary` | The file path was absolute | `cd $CASE` and give the path relative to the case (section 6.4.1) |
| The check in 6.4.4 prints different angles | One of the three entries was not updated | Rewrite the entry that disagrees |
| `leg N needs to requeue but ARGUS_SCRIPT is unset` | `ARGUS_SCRIPT` missing from `qsub -v` | Resubmit with `ARGUS_SCRIPT` |
| `-> did not reach endTime (walltime). Continuing at the same endTime.` | The solver stopped before the `endTime` the job started with, while the job was still alive. A crash does this (`FOAM FATAL` or `---- DETERMINISTIC FAULT, NOT RETRYING ----` earlier in the output), and so does lowering `endTime` in the running case. A real walltime kill prints nothing more (see [./05_submit_and_run.md](./05_submit_and_run.md), section 10a) | First read the end of `$CASE/log/leg<N>/foamRun.log` and the job output for `FOAM FATAL`; fix the cause before resubmitting. Only if there is no error, check whether `endTime` was changed during the run (section 6.7) |
| `0 iterations: already at endTime` on a warm leg | `endTime` was not raised above the latest time (section 6.5.4 step 9) | Any `.converged` this leg wrote reports the **old** angle: rename it (section 6.5.4 step 11), raise `endTime`, resubmit |
| `MAXLEG 4 reached, not requeueing` | The script's default limit of 4 legs was reached | Check why so many legs were needed, then resubmit |
| C_L jumps after an angle step, then falls | The normal restart transient | Wait for the marker; do not read the slope from it |
| `Removing patchGroup 'symmetry' which clashes with patch 0 of the same name` | Warning caused by the recipe's patch name | Harmless; it appears on cold and warm legs alike |

---

## 6.9 What to record for each trim

For every geometry and operating point, keep:

1. The target C_L and the operating point it belongs to.
2. The mesh generation (`wake_refined` or `published`).
3. Every converged (alpha, C_L) pair in order: the sweep points, the solved trim case and
   every warm leg. This is the `trim.history` list of a run card.
4. The angle read back from the case (section 6.4.4), not the angle you intended.
5. The marker text of the final leg, including the gate line and `REPORT THESE`.

Run cards and how they are generated and validated are covered in
[./07_postprocessing_and_report.md](./07_postprocessing_and_report.md). Back to the
overview: [../README.md](../README.md).
