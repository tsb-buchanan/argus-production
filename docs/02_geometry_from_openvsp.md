# 02. Geometry: from an OpenVSP `.vsp3` file to a meshable wing surface

> **At a glance**
> 1. **What it does:** turns one wing's OpenVSP `.vsp3` file into a closed, placed, welded half-wing surface in metres, and proves it is right before you mesh it.
> 2. **Before you start:** finish [01, setup](./01_setup.md) (OpenFOAM-12, Python with numpy and scipy). The wings come as OpenVSP `.vsp3` files: F, H and W from the design study's public GitHub repository, https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow (folder `results/geometries/low_speed/`); B, C and M, which that repository does not hold, from the ARGUS project drive, `\\tudelft.net\staff-umbrella\ARGUS Morphing wing CFD\argus_hpc12_backup\geometry_vsp3\`, which holds all six with a `SHA256SUMS` file. Step 1 shows how to get them and check them. No geometry is shipped in this repository.
> 3. **At the end:** a surface that passes `surfaceCheck` and whose sha256 checksum is listed in `registry/candidates.yaml`, ready to become `wing.stl` in [03, building a case](./03_build_a_case.md).
> 4. **Time:** everything runs on your own machine. The timed steps are short: Steps 4 and 7 take under 20 seconds each, Step 8 about two seconds.

This guide turns an OpenVSP model of one wing into the file the mesh recipes need:
`wing.stl`, a closed half-wing surface in metres. It then shows you how to prove the file
is right **before** you mesh it.

It comes after [01, setup](./01_setup.md) and before [03, building a case](./03_build_a_case.md).
Everything here runs on your own machine. Nothing here touches the cluster.

> **NO GEOMETRY IS SHIPPED. THE `.vsp3` AND STL FILES ARE NOT IN THIS REPOSITORY.**
> The wings come as OpenVSP `.vsp3` files: F, H and W from the design study's public GitHub repository, https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow (folder `results/geometries/low_speed/`); B, C and M, which that repository does not hold, from the ARGUS project drive, `\\tudelft.net\staff-umbrella\ARGUS Morphing wing CFD\argus_hpc12_backup\geometry_vsp3\`, which holds all six with a `SHA256SUMS` file. You download
> the `.vsp3` files, check each one against the checksums in Section 3 (Step 1), and build
> every meshing surface yourself with this guide. The route reproduces all six registered
> surfaces bit for bit (Section 11), so the checksum of your result tells you whether it is
> the right one.

Every command below has been run under OpenFOAM-org 12 and Python 3.11.8, and the output
shown is the output that was observed. The few exceptions are marked **NOT VERIFIED** where
they occur. One script (Step 7) starts `surfaceCheck` by itself; Step 7
says which OpenFOAM it uses and what it prints when it cannot run it.

---

## 1. Read this first: check the input, not only the output

A mesh can be perfectly valid and still be a mesh of the wrong thing. That is the most
expensive mistake this project has made, and it is why this guide exists.

What happened:

1. The wing surfaces were **open at the root**: the root section outline was a free edge with
   no cap (398 open edges on the project's wings).
2. After the wing was tilted to its 5 degree dihedral (Step 5 below), the root outline no
   longer lay in the symmetry plane `y = 0`. Part of it sat above the plane, leaving a slot up
   to **4.679 mm** wide into the inside of the wing.
3. `snappyHexMesh` meshed straight through that slot and filled the **inside** of the wing
   with fluid. The wing became a zero-thickness sheet with fluid on both sides. On one
   7.51 million cell mesh, **76.8 % of the wing patch faces** were back-to-back coincident
   pairs (74.5 % of the patch area), and the wing patch area was 1.905 times the surface
   area, i.e. both sides of a sheet.
4. `checkMesh` **passed every one of these meshes**. It checks the mesh it is given, and a
   flooded mesh is a perfectly valid mesh of the wrong domain.
5. Making the mesh finer made it **worse**, not better: the slot is fixed by the geometry. On
   the project's earlier meshes, whose surface cell at refinement level 9 was 3.891 mm, it was
   1.20 cells wide at level 9 and 4.79 cells wide at level 11. In the shipped recipe (1 m
   background cell, so level 11 is 0.488 mm; see [03](./03_build_a_case.md), Step 3) the same
   slot is about 9.6 cells wide at the production surface level. Finer meshes leak along more
   of the root.

This cost the project an entire 3D campaign: every 3D result built on those meshes was
suspended and had to be redone on repaired surfaces. The check that would have caught it in
one command, `surfaceCheck`, had never been run.

**`surfaceCheck` is an INPUT gate. Run it on every surface, every time, before meshing.**
`checkMesh` is an output gate and cannot see this. Step 8 shows exactly what a good and a bad
result look like.

(The numbers above are the project's own measurements, recorded in
`scripts/assert_surface_closed.py` and `scripts/close_wing_surface.py`. The one exception is
the shipped-recipe figure in item 5, which is 4.679 mm divided by the 0.488 mm level-11 cell.)

---

## 2. What a finished surface looks like

A surface is ready for [03](./03_build_a_case.md) when all of these hold. Each one is checked
somewhere below.

1. **Units are metres**, not feet. OpenVSP models for this wing are in feet (Step 4).
2. **It is placed**: 5 degrees of dihedral and 1.873 degrees of incidence are applied, with
   the root leading edge at x = 1.176528 m, z = -0.12192 m (Step 5).
3. **It is closed**: every edge belongs to exactly two triangles, one piece, one consistent
   normal orientation (Step 8).
4. **The root goes through the symmetry plane**: the root is extended 20 mm to `y = -0.020 m`
   and capped there, so the mesher cuts through solid wing at `y = 0` instead of finding a
   slot (Step 6).
5. **It is welded**: OpenFOAM, which matches points exactly, also sees it as closed (Step 7).
6. **Its sha256 checksum is listed in `registry/candidates.yaml`** (Step 9). The registry's own
   rule is that no case runs on an unregistered geometry.
7. It is **one solid** (one `solid ... endsolid` block), named `wing_closed`.

The meshing surface is a meshing construct, not the wing's true outer surface. Its own area
is 1.3497 m², which is 5.49 % more than the half-wing's wetted area of 1.2794 m² (upper plus
lower surface of one half-wing), because of the 20 mm root extension and its cap (both numbers
from `registry/candidates.yaml`). **Never take an area, span or volume from this file.** Take
wetted area from the solver's wing patch.

---

## 3. Which geometry is which

The six geometries in the results are known by one letter. The letter comes from the run card
file name (`data/run_cards/SW<letter>_r2_trim.json`) and the candidate id from that card's
`"candidate"` field. You can print the mapping yourself (`$REPO` is your clone of this
repository, as in [01_setup.md](./01_setup.md)):

```bash
export REPO=$HOME/argus-production   # or wherever you cloned it
cd $REPO
grep -H '"candidate"' data/run_cards/*.json
```

Expected output:

```text
data/run_cards/SWB_r2_trim.json:  "candidate": "baseline",
data/run_cards/SWC_r2_trim.json:  "candidate": "cte_i002_c04",
data/run_cards/SWF_r2_trim.json:  "candidate": "cfft_b02_c01",
data/run_cards/SWH_r2_trim.json:  "candidate": "chc_g02_c06",
data/run_cards/SWM_r2_trim.json:  "candidate": "mcv2_i002_c01",
data/run_cards/SWW_r2_trim.json:  "candidate": "cffw_b01_c01",
```

| Letter | Candidate id | What it is (from the registry) | Source file (from the design study repository) | Route variant |
|---|---|---|---|---|
| B | `baseline` | the baseline wing | `baseline_wing_only_refined.vsp3` | A |
| C | `cte_i002_c04` | re-optimised morphing candidate | `cte_i002_c04.vsp3` | A |
| M | `mcv2_i002_c01` | morphing candidate, kept as a geometry-validation case | `mcv2_i002_c01.vsp3` | A |
| H | `chc_g02_c06` | conventional hinged concept (discrete flap deflection) | `chc_g02_c06.vsp3` | B |
| F | `cfft_b02_c01` | continuous trailing-edge camber concept | `cfft_b02_c01.vsp3` | B |
| W | `cffw_b01_c01` | distributed twist concept (+1.80 to -1.00 deg of twist about the quarter chord, from 68 % of the half-span out to the tip) | `cffw_b01_c01.vsp3` | B |

Trust the run cards and the checksums for this mapping. If you read the `STATUS` notes under
`chc_g02_c06`, `cfft_b02_c01` and `cffw_b01_c01` in `registry/candidates.yaml`, you will find
them naming the wrong case letters (`SWW_trim`, `SWH_trim` and `SWF_trim` respectively). The
`sha256` values in those same entries agree with the run cards and with the table above.

The two route variants differ in two settings only. Both were verified to reproduce the
registered files exactly (Section 11):

| Variant | Tip cap setting (Step 4) | Weld tolerance (Step 7) |
|---|---|---|
| A (B, C, M) | `--cap delaunay --cap-edge 0.0030` | default (`1e-5` m) |
| B (H, F, W) | `--cap delaunay` (no `--cap-edge`) | `--tol 1e-6` |

For the baseline, use `baseline_wing_only_refined.vsp3` and nothing else. The registry records
that the registered baseline surface is reproduced only from that file, and that a surface
built from the similarly named `baseline_corrected.vsp3` is a different surface, registered as
unused.

If you open the registry you will see two kinds of note on the baseline surface. Comments and
a `MUST_NOT` field mark it as a different aerofoil-interpolation lineage from C and M, which the
registry calls corrected-pair candidates. Its `status` (`SANCTIONED_PRODUCTION_BASELINE`) and
`status_note` record the project's ruling that it is nonetheless the baseline, and that no case
is re-run on lineage grounds. This guide only builds the surface; it does not change which
surface is the baseline.

**Checksums of the source files.** Step 1 saves these lines as `vsp3.sha256` and checks every
`.vsp3` you download against them. The B, H, F and W values are the ones recorded in
`registry/candidates.yaml` (the `delivered_vsp3` entries). The registry does **not** record a
checksum for the C and M source files; the two values below were measured on the files the
project built C and M from.

```text
55cbdc7c48a15d847da602b540b95712b38e893490823aa7b5b997302100ba02  baseline_wing_only_refined.vsp3
aae4a07e4735df8a0b7598fdd445a03eea15aefd871f7575b16bc6d5e17c5bb6  cte_i002_c04.vsp3
54fe540dd4cf685f2ed6d1e4e4b3513b0c81c75f4d423b6d1f02dbc39c8ffb96  mcv2_i002_c01.vsp3
fd591b0c3064aa4fa8610c15d2988c664dce140b03aeb7571c56a3ea7d4b4bee  chc_g02_c06.vsp3
b3712ac50df554d17b37b851c02072cdc53c4794aab4e703cdbaf1ea1bf75a6f  cfft_b02_c01.vsp3
69d2be9b515283a5980a55bf8ede9a7d01a420efddbb27f667683425fbb02652  cffw_b01_c01.vsp3
```

**Checksums of the finished meshing surfaces** (all six are in `registry/candidates.yaml` under
`geometry_3d_laddercap`, and in each run card's `"geometry"` block). The file names are the
registry's. A surface you build yourself gets a different file name (Steps 4 to 7 and
Section 11) but, if it is right, the same checksum: always compare by checksum, not by name.

```text
20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641  baseline_oml_placed_laddercap.stl
790c3ded3568fa877bf8daf1a7296d7f66e7e865e42f4a1e9b90a14501c45d3b  cte_i002_c04_oml_placed_laddercap.stl
dae22d0c4bf3cba849e39a3b59d84d1412ffed3430aa4e10c0ebdf92556941dc  mcv2_i002_c01_oml_placed_laddercap.stl
18329aa33c0c1d971b85b40bffe588cd45eb25cba83683715816e5be29186fc3  chc_g02_c06_oml_placed_laddercap.stl
82bf2e1e570187ab5b3dfa987f122142ab01415995ef5c47fd27838901d0037b  cfft_b02_c01_oml_placed_laddercap.stl
200fa64ca9d6eecc54bbecaa62c384edeab95c2560ce5f514c4a9ffc830de692  cffw_b01_c01_oml_placed_laddercap.stl
```

The project does not hand out finished surfaces: you build them from the `.vsp3` files with
Steps 1 to 7, and Step 9 compares your result with this list. If a finished surface reaches you
some other way (for example a copy from a colleague, or from a mesh folder on the cluster),
check its checksum as in Step 9 before you use it, then run the gates (Step 8 and Section 10).

---

## 4. Before you start

1. Load OpenFOAM-12 in your terminal, as described in [01_setup.md](./01_setup.md), section A2:

   ```bash
   source /opt/openfoam12/etc/bashrc
   foamVersion
   ```

   Expected: `OpenFOAM-12`. If `foamVersion` prints anything else (for example `OpenFOAM-7`),
   or `foamRun` is not found, your shell has loaded a different OpenFOAM version; source the
   line above again. Note that `surfaceCheck` exists in other versions too, so finding it proves
   nothing about the version; if `surfaceCheck` is not found at all, no OpenFOAM is loaded.
   If sourcing prints `unrecognized command-line option '--showme:link'`, an Anaconda
   installation is ahead of the system tools on your `PATH`. That is harmless for this guide,
   which runs nothing in parallel; see [01_setup.md](./01_setup.md), section A5, before any
   parallel run.

   Step 8 and Sections 6 and 7 run OpenFOAM programs directly and need this terminal. The
   script in Step 7 finds OpenFOAM-12 by itself (Step 7, item 1).

2. Check Python (see [01_setup.md](./01_setup.md), section A8). The geometry scripts need
   numpy, and `weld_surface.py` (Step 7) also needs scipy:

   ```bash
   python3 --version
   python3 -c "import numpy, scipy; print(numpy.__version__, scipy.__version__)"
   ```

   The route was verified with `Python 3.11.8` and `1.26.4 1.11.4`. The bit-for-bit
   reproduction in Section 11 was only tested with these versions. If you see
   `ModuleNotFoundError`, this terminal's `python3` is not the Python from
   [01_setup.md](./01_setup.md), section A8. One way this happens: in your OpenFOAM terminal
   ([01_setup.md](./01_setup.md), section A5, `/usr/bin` first on your `PATH`), `python3` is
   Ubuntu's own, which has no numpy. Two set-ups work, and both have been run:

   1. **Your Python terminal** from [01_setup.md](./01_setup.md), section A5: your Python first
      on the `PATH` and OpenFOAM-12 loaded (it prints the harmless `--showme:link` message
      above). Every command in this guide runs there. This is the simplest choice.
   2. A terminal with your Python first on the `PATH` and **no** OpenFOAM loaded, for every
      `python3` command (Step 7 finds OpenFOAM-12 by itself, Step 7, item 1), plus **your
      OpenFOAM terminal** from [01_setup.md](./01_setup.md), section A5, for Step 8 and
      Sections 6 and 7. Set `$REPO` and `$GEOM` (item 3) in both terminals, and `cd $GEOM/B`
      in your OpenFOAM terminal before Step 8.

3. Set two variables. `$REPO` is your clone of this repository (as in
   [01_setup.md](./01_setup.md)). `$GEOM` is a working folder for geometry, **outside** the
   repository, because the files are large (the finished surfaces are 52.6 to 53.3 MB each,
   about 120 MB per geometry including intermediates) and because `surfaceCheck` leaves
   extra files in whatever folder you run it from.

   ```bash
   export REPO=$HOME/argus-production   # or wherever you cloned it
   export GEOM=$HOME/argus-geometry
   mkdir -p $GEOM/vsp3 $GEOM/B
   ```

4. Step 1 puts the `.vsp3` files from the design study repository into `$GEOM/vsp3/`.

The rest of this guide builds **B, the baseline**, as the worked example. For another letter,
change the file names and use that letter's route variant from Section 3 (Section 11 lists
the exact commands for all six).

---

## 5. The steps

### Step 1. Get the `.vsp3` files and check their checksums (Section 3)

F, H and W come from the design study's public GitHub repository, https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow,
folder `results/geometries/low_speed/`. B, C and M are not in that repository; they, and copies
of F, H and W, are on the ARGUS project drive, `\\tudelft.net\staff-umbrella\ARGUS Morphing wing CFD\argus_hpc12_backup\geometry_vsp3\`, with a `SHA256SUMS` file. The
repository holds more `.vsp3` files than the six in Section 3, some with similar names (Section 3
explains why `baseline_corrected.vsp3` is not the baseline). **The checksum decides which file
is right, not the name.** Do the four parts below in order.

**1a. Save the checksum list.** Save the source checksum list from Section 3 as
`$GEOM/vsp3/vsp3.sha256`:

```bash
cd $GEOM/vsp3
cat > vsp3.sha256 <<'EOF'
55cbdc7c48a15d847da602b540b95712b38e893490823aa7b5b997302100ba02  baseline_wing_only_refined.vsp3
aae4a07e4735df8a0b7598fdd445a03eea15aefd871f7575b16bc6d5e17c5bb6  cte_i002_c04.vsp3
54fe540dd4cf685f2ed6d1e4e4b3513b0c81c75f4d423b6d1f02dbc39c8ffb96  mcv2_i002_c01.vsp3
fd591b0c3064aa4fa8610c15d2988c664dce140b03aeb7571c56a3ea7d4b4bee  chc_g02_c06.vsp3
b3712ac50df554d17b37b851c02072cdc53c4794aab4e703cdbaf1ea1bf75a6f  cfft_b02_c01.vsp3
69d2be9b515283a5980a55bf8ede9a7d01a420efddbb27f667683425fbb02652  cffw_b01_c01.vsp3
EOF
```

**1b. Get the files.** Copy B, C and M (or all six) from the project drive folder above
into `$GEOM/vsp3/`. For F, H and W you can instead clone the design study's repository:

```bash
git clone --depth 1 https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow.git $HOME/design-study
```

Checked on 2026-09-25: the clone holds F, H and W under `results/geometries/low_speed/` with the
registered checksums, and 1c below found and copied exactly those three. It does not hold B, C
or M, so 1c prints no line for them from the clone; take those from the project drive.

**1c. If you cloned it, copy the six files out of the clone,** wherever they sit in it and
whatever they are called there:

```bash
cd $GEOM/vsp3
find $HOME/design-study -name '*.vsp3' -exec sha256sum {} + \
    | while read -r sum path; do
          name=$(awk -v s="$sum" '$1 == s {print $2}' vsp3.sha256)
          if [ -n "$name" ]; then cp -v "$path" "$name"; fi
      done
```

The command lists every `.vsp3` in the clone with its checksum, looks each checksum up in your
list, and copies a file only if its checksum is there. It saves the copy under the name the list
gives, so a registered file stored under another name in the clone still arrives under the name
the later steps use.

**What you should see:** one line per file copied, such as
`'<path in the clone>/baseline_wing_only_refined.vsp3' -> 'baseline_wing_only_refined.vsp3'`,
naming all six files of the list. If the clone holds identical copies of a file in more than one
folder, that file is copied once per copy, so you see more than six lines. That is harmless:
every copy has the registered checksum, and the last one simply replaces an identical file.
Step 1d confirms what you end up with.

The command was checked on a folder of 27 `.vsp3` files in subfolders. It held
`baseline_corrected.vsp3`, a second identical copy of `baseline_wing_only_refined.vsp3`, one
registered file under a different name, and a file named `mcv2_i002_c01.vsp3` with other content.
It printed seven lines, left out the file with the wrong content, and gave the six registered
files their registered names.

**If it does not:**

1. `find: '.../design-study': No such file or directory`, or no output at all: the folder is
   not where you cloned the repository, or it holds none of the six files with the registered
   content. Check that `$HOME/design-study` is where you cloned it.
2. Only F, H and W appear: that is expected from the clone, which does not hold B, C or M.
   Copy those three from the project drive. If F, H or W is missing too, the repository no
   longer holds it with the registered content; take it from the project drive as well.

**1d. Check what you have:**

```bash
cd $GEOM/vsp3
sha256sum -c --ignore-missing vsp3.sha256
```

**What you should see:** one line ending in `OK` for each file you have, for example
`baseline_wing_only_refined.vsp3: OK`. **Count them:** six `OK` lines if you want all six
geometries. `--ignore-missing` skips, without a word, every name in the list that is not in the
folder, so a missing line is the only sign that a file is missing or under another name.

**If it does not:**

1. A line says `FAILED`: that file is not the one the project built from. The design study
   repository may hold another version under the same name. Do not continue with it; ask the
   project which file to use.
2. Fewer `OK` lines than you expected: a file you need is missing, or it is in the folder under
   a different name (this can happen with files you downloaded one by one). List the checksums
   of what you have:

   ```bash
   sha256sum *.vsp3
   ```

   If a file's checksum is in `vsp3.sha256`, rename the file to the name listed beside that
   checksum (`mv <your file> <name in the list>`) and run the check again. If no file has the
   checksum, the file is missing: get it as in 1b or 1c.
3. It prints only `sha256sum: vsp3.sha256: no file was verified`: none of your files has one of
   the names in the list, so nothing was checked. Do item 2.

### Step 2. Keep only the wing (strip everything else by name)

An OpenVSP model can hold many components besides the wing: slats, flaps, tails, slat tracks.
Only the wing may go into the meshing surface, and you select it **by name**.

List the wing components in the file:

```bash
cd $GEOM/B
python3 $REPO/scripts/vsp3_sections.py ../vsp3/baseline_wing_only_refined.vsp3 --list
```

**What you should see:** one line, `cruise_wing`. All six project files contain exactly one
wing component, called `cruise_wing`.

**If it does not:** if it lists several names, the file is a multi-component model. The full
wind-tunnel model of this wing (not in this repository) lists twelve: `trapezoidal_reference`, `cruise_wing`, `main_wing`, `slat_inboard`,
`vane_inboard`, `flap_inboard`, `slat_outboard`, `vane_outboard`, `flap_outboard`,
`main_and_flap`, `htail`, `vtail`. If you give such a file to the next step without naming the
wing, it stops with (shortened here):

```text
ValueError: <file>: need exactly one wing geom, found ['trapezoidal_reference', 'cruise_wing', ...]
```

That is the safe behaviour. Always pass `--geom cruise_wing` (or whatever the wing is called in
your file) in Step 4, so the choice is written down in the command rather than left to chance.

If you only have an STL exported by OpenVSP, the non-wing parts are separate named solids
inside the STL. Section 7 shows how to strip them.

### Step 3. Also read the placement values (you need them in Step 5)

The wing's dihedral and incidence are **not** stored on the wing component. They are stored on
its parent, a component called `wing_comps`. Print the parent's values:

```bash
grep -oE '<(X_Rel_Location|Y_Rel_Location|Z_Rel_Location|X_Rel_Rotation|Y_Rel_Rotation|Z_Rel_Rotation) Value="[^"]*"' \
    ../vsp3/baseline_wing_only_refined.vsp3 | head -6
```

**What you should see**, identical for all six project files (the parent's block comes first
in each file, which is why the first six matches are the ones you want):

```text
<X_Rel_Location Value="3.859999999999999876e+00"
<X_Rel_Rotation Value="5.000000000000000000e+00"
<Y_Rel_Location Value="0.000000000000000000e+00"
<Y_Rel_Rotation Value="1.872999999999999998e+00"
<Z_Rel_Location Value="-4.000000000000000222e-01"
<Z_Rel_Rotation Value="0.000000000000000000e+00"
```

These are 5.000 deg dihedral (rotation about x), 1.873 deg incidence (rotation about y), no
rotation about z, and a translation of (3.860, 0, -0.400) **feet**. The placement script in
Step 5 has exactly these values built in.

**If it does not:** if your file shows different numbers, stop. The built-in translation cannot
be changed from the command line, and the project has to be asked.

### Step 4. Loft the wing and convert feet to metres

**The units trap.** OpenVSP models and OpenVSP STL exports for this wing are in **feet**.
OpenFOAM assumes metres. A surface left in feet is 3.28 times too large in every direction, and
nothing in OpenFOAM will warn you. The conversion is feet times 0.3048 = metres.

`scripts/vsp3_to_stl.py` builds the wing surface directly from the section shapes stored in the
`.vsp3` (OpenVSP itself is not needed) and **applies the 0.3048 scaling by default** (its
`--scale` option defaults to `0.3048`; leave it alone).

```bash
cd $GEOM/B
python3 $REPO/scripts/vsp3_to_stl.py ../vsp3/baseline_wing_only_refined.vsp3 --geom cruise_wing \
    --out baseline_oml.stl --nchord 200 --nspan 240 --blunt-te 0.0063 \
    --cap delaunay --cap-edge 0.0030 > step1.json
```

What the options mean:

1. `--geom cruise_wing`: the wing component, by name (Step 2).
2. `--out baseline_oml.stl`: the output. A summary is also written beside it as
   `baseline_oml.json`, and printed (here into `step1.json`).
3. `--nchord 200`: points per surface around each section (cosine spacing, dense at leading
   and trailing edge).
4. `--nspan 240`: spanwise stations. Every station defined in the file (19 for the baseline) is
   kept; the rest are placed in between.
5. `--blunt-te 0.0063`: restores the blunt trailing edge. OpenVSP closes the trailing edge to
   a sharp point; this undoes that. 0.0063 is the setting of the correction, chosen to match the
   0.63 % chord blunt trailing edge measured on the project's reference aerofoil section. It is
   not the number you will measure on the finished 3D surface: the registry records a measured
   gap of about 0.68 % of the local chord on surfaces built this way (`te_gap_pct_c: 0.6769`),
   and `check_oml_gates.py` (Section 10c) reports 0.78 % at the root and 0.85 % at the tip by
   its own definition. All of these are expected.
6. `--cap delaunay`: the tip cap. **Always pass it.** The default, `fan`, is kept only so old
   files can be regenerated, and it is a measured defect as a meshing surface (one tip vertex
   shared by 398 triangles, tip triangle aspect ratios up to 35,212).
7. `--cap-edge 0.0030`: target tip-cap edge length, 3.0 mm. Variant A only; leave it out for
   H, F and W (Section 3).

It takes under 20 seconds and prints nothing while it works (the summary goes into
`step1.json`). Check the summary:

```bash
grep -E '"(half_span_m|half_area_m2|full_area_m2|root_chord_m|tip_chord_m|scale_applied|triangles)"' step1.json
```

**What you should see** for the baseline:

```text
    "triangles": 191400
    "scale_applied": 0.3048,
  "half_span_m": 1.8287989056393246,
  "half_area_m2": 0.6204632184259208,
  "full_area_m2": 1.2409264368518416,
  "root_chord_m": 0.6400810263333012,
  "tip_chord_m": 0.15240002665074773,
```

How to read it, with the frame of each number stated:

1. `half_span_m` 1.8288 m is the span of one half-wing measured along the wing before any
   dihedral, in metres. It is half of the reference span `bref` 3.6576 m.
2. `full_area_m2` 1.24093 m² is the planform area of the whole wing (both halves). It matches
   the project's primary reference area `Sref` 1.24092 m² (called the DSO basis in the
   project's files; DSO is the design-optimisation study the candidate wings came from).
   `half_area_m2` 0.620463 m² matches the half-model reference area 0.620462 m² that the force
   coefficients use.
3. Root chord 0.64008 m and tip chord 0.15240 m.

**If it does not:** if `half_span_m` reads about 6.0, **the file is in feet.** That happens if
`--scale 1.0` was passed. Note that the summary then still labels the value `half_span_m` even though the number
is in feet: the label does not prove the unit, the value does. Rerun without `--scale`.

This surface is an in-house loft: it matches the `.vsp3` section shapes exactly at the stations
defined in the file and interpolates linearly between them. It is not OpenVSP's own tessellation.
The summary says so in its `limitation` field.

### Step 5. Apply the placement transform (dihedral and incidence)

The loft from Step 4 is built in the wing's own frame: root leading edge at the origin, no
dihedral, no incidence. The dihedral and incidence live on the parent component (Step 3), which
the loft does not read. Step 5 applies them.

**Why this step exists.** The project's first 3D surface skipped it. That surface had exactly
the right span (3.657598 m) and area (1.240926 m²), so every planform check passed, and it was
still the wrong wing: no dihedral and no incidence. The registry records it as superseded for
that reason. A span or area check cannot catch this; Section 10c shows the check that does.

```bash
python3 $REPO/scripts/apply_placement_transform.py --loft baseline_oml.stl \
    --out baseline_oml_placed.stl --order xy --json-out placement.json
```

1. `--order xy` applies the incidence first and the dihedral second (R = Rx(5.000) Ry(1.873)).
   The project chose this order by measurement against an OpenVSP export of the same model.
   **Always pass `--order xy`.** The default, `auto`, needs an OpenVSP export given with
   `--reference`, and without one it stops with `--order auto requires --reference`.
2. The translation, (3.860, 0, -0.400) ft = (1.176528, 0, -0.12192) m, is built in.

**What you should see:** a summary (also saved in `placement.json`) that contains these lines:

```text
    "rotation_centre": "wing root LE = loft origin (verified, |root_LE| = 0.00e+00 m)",
    "order_applied": "xy",
    "translation_m": [
      1.176528,
```

**If it does not:**

1. `--order auto requires --reference`: `--order xy` is missing. Add it (item 1 above).
2. `RuntimeError: loft root LE is at [1.176528, 0.0, -0.12192], not the origin`: the input was
   already placed. **Do not run this step twice.** The script checks that its input still has
   the root leading edge at the origin, and stops on an already placed file, which is what you
   want. Use the Step 4 output as `--loft`.

### Step 6. Close the root (extend it through the symmetry plane)

After Step 5 the surface is still open at the root, and the tilted root outline no longer lies
in `y = 0`: it runs from 3.692 mm below the plane to 4.679 mm above it. This is exactly the slot
from Section 1.

`scripts/close_wing_surface.py` drops a skirt from every root edge straight down to
`y = -0.020 m` and caps it there with a flat plate. The cap is **not** at `y = 0`: a cap lying
on the symmetry plane would make the mesher choose between wing and symmetry at the same place.
With the wing extending 20 mm through the plane, the mesh boundary at `y = 0` cuts through solid
wing. It does not move any existing triangle.

Where the added surface ends up: the cap, and the part of the skirt below `y = 0`, lie outside
the fluid domain (the recipes' domain starts at `y = 0`) and are cut away by the mesher. Where
the root outline sits above the plane (up to 4.679 mm), the upper part of the skirt is what fills
the slot. That part is inside the domain and becomes a small piece of wetted wing surface: for
the baseline, 2.2e-3 m² of the 6.6e-2 m² the step adds, about 0.17 % of the half-wing's wetted
area of 1.2794 m².

```bash
python3 $REPO/scripts/close_wing_surface.py baseline_oml_placed.stl --write
```

`--write` replaces the file in place, and only after the script's own checks pass. Without
`--write` it only reports.

**What you should see:**

```text
  baseline_oml_placed.stl                398 rim edges skirted, 396 cap triangles, 192592 facets total
```

For H, F and W the last number is `195266`.

**If it does not:**

1. `already closed`: the input was already closed and nothing was done. The usual reason is
   that you ran Step 6 on this file before (it changes the file in place); that is fine, go on
   to Step 7. For OpenVSP STL exports, see Section 7.
2. `REFUSED`: read the reason and go back to Step 4.

### Step 7. Weld the root seam

The skirt from Step 6 was placed against the root outline, not joined into it: each root point
now exists twice, at positions that differ by less than 1 micrometre. The case-level Python
checker in this repository rounds coordinates to the nearest micrometre before matching points,
and calls this surface closed. **OpenFOAM matches points exactly**, so to OpenFOAM the surface
is still open along the root (Step 8 shows the evidence). `scripts/weld_surface.py` merges the
duplicated points.

```bash
python3 $REPO/scripts/weld_surface.py baseline_oml_placed.stl --in-place
```

For H, F and W add `--tol 1e-6` (Section 3).

The script runs `surfaceCheck` before and after welding, writes the welded surface to a
temporary file first, and replaces your file (`--in-place`) **only if** `surfaceCheck` then
calls it closed, one zone, one orientation. It takes 10 to 15 seconds.

**What you should see** for the baseline, in your Python terminal from Section 4 (OpenFOAM-12
loaded):

```text
=== baseline_oml_placed.stl ===
  surfaceCheck: OpenFOAM build 12-86e126a7bc4d, already loaded (WM_PROJECT_DIR=/opt/openfoam12)
  before: {'closed': False, 'zones': 2, 'one_orientation': False, 'illegal': False}
  merged 399 vertices (96695 -> 96296), max displacement 9989.2 nm
  triangles 192592 -> 192588 (4 degenerate dropped)
  after:  {'closed': True, 'zones': 1, 'one_orientation': True, 'illegal': False}
  PROMOTED over baseline_oml_placed.stl
```

The `surfaceCheck:` line says which OpenFOAM ran `surfaceCheck` (item 1 below). In a terminal where
OpenFOAM-12 is not loaded, it ends `sourced /opt/openfoam12/etc/bashrc (default)` instead. Both
were run, and they gave the same output otherwise and the same checksum. The characters after
`12-` identify the exact build and may differ on your machine.

`'illegal': False` means `surfaceCheck` reported no illegal triangles. For H, F and W with
`--tol 1e-6` the middle lines read `merged 397 vertices (98032 -> 97635), max displacement
778.3 nm` and `triangles 195266 -> 195266 (0 degenerate dropped)`.

The `9989.2 nm` and the `4 degenerate dropped` are **not** the root seam. The seam is 397 of the
399 merges, all under 1 micrometre. The default tolerance (`1e-5` m, 10 micrometres) also merges
two pairs of genuinely separate points about 10 micrometres apart near the tip trailing edge,
and that removes 4 triangles. The registered B, C and M surfaces include this merge, which is why
their route uses the default. With `--tol 1e-6` only the seam is welded (`merged 397 vertices`,
`max displacement 778.3 nm`, `0 degenerate dropped`); the registered H, F and W surfaces were
built that way.

**If it does not:** `FAILED: surfaceCheck could not be run ...` (exit status 2) is item 2 below.
`NOT PROMOTED` (exit status 1) is at the end of item 2.

Things to know about this step:

1. **Which OpenFOAM it uses.** The script starts `surfaceCheck` itself, and only ever under
   OpenFOAM-12. It takes, in this order: an OpenFOAM-12 already loaded in your terminal; else
   the `etc/bashrc` named by the environment variable `ARGUS_OF_BASHRC`; else
   `/opt/openfoam12/etc/bashrc`. It then reads the version back from `surfaceCheck`'s own
   `Build` line and refuses anything that is not 12. It never falls back to another OpenFOAM,
   even when one is installed beside OpenFOAM-12. So with OpenFOAM-12 installed in
   `/opt/openfoam12` as in [01_setup.md](./01_setup.md), this step works in any terminal: it
   was run with no OpenFOAM loaded, with OpenFOAM-12 loaded, and with OpenFOAM-org 7 loaded, and
   all three used OpenFOAM-12 and gave the same checksum. Step 8 still needs OpenFOAM-12 loaded.
   If your OpenFOAM-12 is somewhere else, either load it first or tell the script where it is:

   ```bash
   export ARGUS_OF_BASHRC=<your OpenFOAM-12 folder>/etc/bashrc
   ```

   The `surfaceCheck:` line then ends `(ARGUS_OF_BASHRC)`. An OpenFOAM-12 loaded in the
   terminal still comes first: the line then ends `already loaded (...); ARGUS_OF_BASHRC ignored`.

2. **If `surfaceCheck` cannot be run**, the script stops before it changes anything, with exit
   status 2. Your file is left as it was, and no temporary file is left behind. It prints, for
   example:

   ```text
   FAILED: surfaceCheck could not be run, so nothing was checked and nothing was promoted.
     reason: OpenFOAM-12 is not loaded in this shell and its bashrc was not found at <path> (ARGUS_OF_BASHRC)
     if OpenFOAM-12 is missing: load it (source /opt/openfoam12/etc/bashrc) or set ARGUS_OF_BASHRC to its etc/bashrc, and rerun.
   ```

   The `reason:` line says what went wrong. Each of these has been reproduced:

   1. `its bashrc was not found at <path> (default)`: OpenFOAM-12 is not loaded and is not at
      `/opt/openfoam12`. Load it, or set `ARGUS_OF_BASHRC` (item 1).
   2. `its bashrc was not found at <path> (ARGUS_OF_BASHRC)`: `ARGUS_OF_BASHRC` names a file that
      does not exist. A wrong value is an error: the script does not then try
      `/opt/openfoam12`. Correct it, or remove it with `unset ARGUS_OF_BASHRC`.
   3. `surfaceCheck ran under OpenFOAM build 7-..., not OpenFOAM-12`: `ARGUS_OF_BASHRC` points at
      another OpenFOAM version. Point it at OpenFOAM-12.
   4. `surfaceCheck on <file> did not complete ...`, with `invalid fileName` further along the
      line: OpenFOAM does not accept the file name, for example because it contains a space.
      Rename the file; the names in this guide are safe.

   `NOT PROMOTED. surfaceCheck still refuses it; the temporary is at <file>` (exit status 1) is
   a different case: `surfaceCheck` did run, and the welded surface still fails it. The
   `after:` line shows how. For example, on a surface that was never closed (Step 6 skipped) it
   read `after:  {'closed': False, 'zones': 1, 'one_orientation': True, 'illegal': False}`. Your
   file is left as it was. The temporary file (in `/tmp`, or in `$TMPDIR` if you set it) can be
   deleted. Go back to Step 6.
3. The welded file is created readable by you only (permissions `-rw-------`). That is harmless.
4. The script leaves `surfaceCheck` debris in the current folder: `badFaces`, `problemFaces`,
   `baseline_oml_placed_0.obj`, `baseline_oml_placed_1.obj`, and next to the STL
   `zone_baseline_oml_placed.vtk`. You can delete them. This is why you work in `$GEOM`, not in
   the repository.

### Step 8. Run `surfaceCheck`: the input gate

This is the step that was missing when the 3D campaign was lost (Section 1). Run it on the
finished surface. `-checkSelfIntersection` adds a self-intersection test; it took about two
seconds here.

```bash
surfaceCheck -checkSelfIntersection baseline_oml_placed.stl > surfaceCheck.log 2>&1
grep -E "Build|Triangles|Vertices|Bounding Box|illegal|Surface is|unconnected parts|Number of zones|orientation|nearby" surfaceCheck.log
```

**What you should see.** For the baseline the `grep` prints exactly these ten lines:

```text
Build  : 12-86e126a7bc4d
Triangles    : 192588
Vertices     : 96296
Bounding Box : (1.17653 -0.02 -0.175402) (2.33626 1.82541 0.0136577)
Surface has no illegal triangles.
Found 0 nearby points.
Surface is closed. All edges connected to two faces.
Number of unconnected parts : 1
Number of zones (connected area with consistent normal) : 1
Surface is not self-intersecting
```

The first line must start with `Build  : 12-`. Anything else means the check ran under a
different OpenFOAM version: go back to Section 4 and load OpenFOAM-12. (The characters after
`12-` identify the exact build and may differ on your machine.)

For reference, this is the relevant part of the full log, `surfaceCheck.log`, for the same
surface (the banner, the per-bin quality table and the edge-length lines are left out):

```text
Reading surface from "baseline_oml_placed.stl" ...

Statistics:
Triangles    : 192588
Vertices     : 96296
Bounding Box : (1.17653 -0.02 -0.175402) (2.33626 1.82541 0.0136577)

Region	Size
------	----
wing_closed	192588


Surface has no illegal triangles.

Triangle quality (equilateral=1, collapsed=0):

    min 6.23599e-13 for triangle 192364
    max 0.812855 for triangle 590

Dumping bad quality faces to "badFaces"
Paste this into the input for surfaceSubset

Checking for points less than 1e-6 of bounding box ((1.15973 1.84541 0.18906) metre) apart.
Found 0 nearby points.

Surface is closed. All edges connected to two faces.

Number of unconnected parts : 1

Number of zones (connected area with consistent normal) : 1

Checking self-intersection.
Surface is not self-intersecting
```

Line by line, what must be true:

1. `Surface has no illegal triangles.`
2. `Surface is closed. All edges connected to two faces.`
3. `Number of unconnected parts : 1` (one piece).
4. `Number of zones (connected area with consistent normal) : 1` (all normals agree), and no
   line saying `More than one normal orientation.`
5. `Found 0 nearby points.`
6. `Surface is not self-intersecting`.
7. **The bounding box is in metres and in the placed frame** (half model, x streamwise, y
   spanwise, z up). For every one of the six geometries it is close to
   `(1.17653 -0.02 -0.175402) (2.33626 1.82541 0.0137)` (W, whose tip is twisted, reads
   `(2.33617 1.82542 0.0147544)` for the maximum):
   1. minimum x 1.17653 m is the root leading edge after the placement translation;
   2. minimum y **-0.02** m is the root extension from Step 6 (if it reads about `-0.0037`, the
      root was never closed; if about `0`, the surface was never placed);
   3. maximum y **1.82541** m is the tip. It is a little less than the 1.8288 m half-span of
      Step 4 because y is measured straight across the flow, and the 5 degree dihedral tilts
      the span upward. It is not exactly 1.8288 times cos 5 degrees (1.8218 m) because the
      incidence rotation, made about the root leading edge, has also lowered the swept-back
      tip, and tilting the span moves points that sit low a little outward in y. A value near
      **6** means the file is in **feet** (Section 6).

The quality minimum of `6.23599e-13` is expected. The faces listed in `badFaces` (five of them
for the baseline) all lie in the root cap at `y = -0.020 m`, which is outside the fluid domain.

**If it does not.** These are real outputs from the same wing at earlier steps.

A surface that was never closed (Step 6 skipped). Note that it is still one piece with one
normal orientation, which is why it can look healthy at a glance:

```text
Surface is not closed since not all edges connected to two faces:
    connected to one face : 398
    connected to >2 faces : 0
...
Number of unconnected parts : 1

Number of zones (connected area with consistent normal) : 1
```

A surface that was closed but not welded (Step 7 skipped). The case-level checker and the
cluster mesh job's built-in gate both call this surface closed (Section 10); `surfaceCheck`
does not:

```text
Found 395 nearby points.

Surface is not closed since not all edges connected to two faces:
    connected to one face : 796
    connected to >2 faces : 0
...
Number of unconnected parts : 2
...
Number of zones (connected area with consistent normal) : 2
More than one normal orientation.

Checking self-intersection.
Surface is self-intersecting at 216 locations.
```

The self-intersection here is at the unwelded root seam (all 216 points written to
`selfInterPoints.obj` lie within 5 mm of `y = 0`) and goes away after Step 7. The surface that
was never closed reports `Surface is not self-intersecting`.

**Anything other than the good pattern means: do not mesh this surface.** Go back to the step
that produces the property that failed.

`surfaceCheck` also writes diagnostic files. Into the folder you run it from: `badFaces` when
some triangles have poor quality (the good baseline has five, so you will see it), `problemFaces`
when edges are open, `<name>_0.obj`, `<name>_1.obj` and so on (one per piece) when the surface
is in more than one piece, and `selfInterPoints.obj` with `-checkSelfIntersection`. Next to the
STL: `zone_<name>.vtk`, also when the surface is in more than one piece. They are diagnostics
only and can be deleted.

Optional: to look at the surface, open it in ParaView ([01_setup.md](./01_setup.md), section A7;
it needs a display):

```bash
paraview baseline_oml_placed.stl   # NOT VERIFIED: the program exists, but no display was available to start it
```

### Step 9. Compare the checksum with the registry

For the surface you built in Steps 4 to 8:

```bash
sha256sum baseline_oml_placed.stl
echo "20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641  baseline_oml_placed.stl" | sha256sum -c -
```

**What you should see:**

```text
20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641  baseline_oml_placed.stl
baseline_oml_placed.stl: OK
```

**If it does not:** the second command prints

```text
baseline_oml_placed.stl: FAILED
sha256sum: WARNING: 1 computed checksum did NOT match
```

Your surface is not the registered one. Work through these in order:

1. Run Step 8 on it. A step was probably skipped: for example, a closed but unwelded surface
   (Step 7 skipped) gives exactly this `FAILED`.
2. Compare the commands you ran with the ones in Section 11, flag by flag, for the letter you
   are building.
3. Check your Python, numpy and scipy versions (Section 4, item 2).
4. Then see the last paragraph of Section 11.

If you have finished surfaces under the registry's file names (Section 3 says how that can
happen), check them all at once. Go to the folder that holds them (shown here as
`<folder with the STL files>`) and save the finished-surface list from Section 3 there:

```bash
cd <folder with the STL files>
cat > surfaces.sha256 <<'EOF'
20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006c641  baseline_oml_placed_laddercap.stl
790c3ded3568fa877bf8daf1a7296d7f66e7e865e42f4a1e9b90a14501c45d3b  cte_i002_c04_oml_placed_laddercap.stl
dae22d0c4bf3cba849e39a3b59d84d1412ffed3430aa4e10c0ebdf92556941dc  mcv2_i002_c01_oml_placed_laddercap.stl
18329aa33c0c1d971b85b40bffe588cd45eb25cba83683715816e5be29186fc3  chc_g02_c06_oml_placed_laddercap.stl
82bf2e1e570187ab5b3dfa987f122142ab01415995ef5c47fd27838901d0037b  cfft_b02_c01_oml_placed_laddercap.stl
200fa64ca9d6eecc54bbecaa62c384edeab95c2560ce5f514c4a9ffc830de692  cffw_b01_c01_oml_placed_laddercap.stl
EOF
sha256sum -c --ignore-missing surfaces.sha256
```

**What you should see:** one line ending in `OK` per file you have, for example
`baseline_oml_placed_laddercap.stl: OK`. As in Step 1d, `no file was verified` means none of your
files has one of these names; use the name-independent search below instead.

Whatever a file is called, you can find out which registry entry it belongs to by searching the
registry and the run cards for its checksum:

```bash
cd $REPO
grep -n "$(sha256sum $GEOM/B/baseline_oml_placed.stl | cut -d' ' -f1)" registry/candidates.yaml data/run_cards/*.json | cut -c1-110
```

**What you should see** (the lines are long, so they are cut at 110 characters):

```text
registry/candidates.yaml:70:  geometry_3d_laddercap: {file: geometry/derived/wing3d/baseline_oml_placed_ladder
data/run_cards/SWB_r2_trim.json:13:    "sha256": "20bd6acdee770b8241930d3f81dedac4dedd33f78f57dc0abc09dca8b006
```

**If it does not:** no output means the file is not a registered geometry. The registry's rule
is that no case runs on an unregistered geometry, so do not use it for results until the project
has registered it.

Copies kept in the delivered mesh cases on the cluster are not all the registered files. For F,
H and W, the surface in each mesh case is the closed but unwelded version, the input of Step 7
(checked on HPC12 by reading the files). It has the same 195,266 triangles but 98,032 points
against the registered file's 97,635. Its points differ from the registered file only along the
root seam, by at most 7.8e-07 m (the weld's `max displacement 778.3 nm` in Step 7), or 5.0e-07 m
in any one coordinate, which is the figure the registry quotes. Exact point matching finds 796
open edges there, which is what `surfaceCheck` reports for such a surface in Step 8. The `STATUS`
notes in `registry/candidates.yaml` describe these copies as the same surface rewritten at 7
significant figures. They are not: they carry 10 significant figures, like the registered file.
Never mesh such a copy. Use the original file, or a copy whose checksum you have checked as
above, and run Step 8 on it. At the production surface level, the delivered F, H and W meshes
leaked into the wing interior through this seam part-way through meshing; the consequences are
in the [README](../README.md) (Known gaps, item 8).

---

## 6. Is my surface in feet? How to tell, and how to fix it

The quickest test is the bounding box from `surfaceCheck` (Step 8). Real values for this wing:

| Surface | Units | Bounding box maximum y |
|---|---|---|
| loft from Step 4, before placement | metres | 1.8288 |
| the same loft built with `--scale 1.0` | feet | 6 |
| finished meshing surface (placed, with root extension) | metres | 1.82541 |
| OpenVSP's own STL export, +y half | feet | 5.98885 |
| the same export after scaling by 0.3048 | metres | 1.8254 |

The ratio between feet and metres is 1/0.3048 = 3.28. For this wing the maximum y should be
about 1.83 (metres). If it is about 6, the file is in feet.

**Fix for a file you built yourself:** rerun Step 4 without `--scale`.

**Fix for a file you were given in feet** (for example an OpenVSP STL export): scale it with
OpenFOAM-12's `surfaceTransformPoints`.

```bash
surfaceTransformPoints "scale=(0.3048 0.3048 0.3048)" wing_in_feet.stl wing_in_metres.stl
surfaceCheck wing_in_metres.stl | grep "Bounding Box"
```

Expected: the log ends with `Writing surf to "wing_in_metres.stl" ...` and `End`, and the new
bounding box is the old one times 0.3048. On the feet version of the Step 4 loft this gave
`(0 0 -0.0471814) (1.1601 1.8288 0.0479245)`, identical to the loft built in metres.

Scale **once**. Scaling a file that is already in metres makes it 3.28 times too small; the
maximum y would then read about 0.56 (1.8254 times 0.3048).

---

## 7. If you only have an OpenVSP STL export

The production route is the `.vsp3` route above. An STL exported by OpenVSP is **not** a
substitute for it, for four reasons found on the project's own export:

1. Its trailing edge is sharp; the production surfaces have the 0.63 % chord blunt trailing
   edge (Step 4). The registry records the two as not interchangeable.
2. The project's 6,240-triangle export is far too coarse to mesh (the registry calls it not
   meshable).
3. Each half-wing is already capped at the root, but the cap lies on `y = 0` (tilted by the
   dihedral, reaching only 3.655 mm below the plane). That is the "capped on the plane" failure
   from Step 6. `scripts/close_wing_surface.py` only closes open roots, so on this file it prints
   `already closed` and changes nothing.
4. It contains non-wing solids.

It is still useful as a **cross-check of the placement**: after the steps below, the export's
bounding box agreed with the placed loft from Step 5 to within 0.72 mm in every direction, while
the unplaced loft was 1.18 m away in x.

To strip non-wing solids from an ASCII STL (one whose first five characters are `solid`; only
ASCII files were tested):

1. List the solids:

   ```bash
   head -c 5 cruise_wing.stl; echo
   grep -a "^solid" cruise_wing.stl
   ```

   On the project's export this printed `solid`, then:

   ```text
   solid 1_cruise_wing_S_Surf0
   solid 2_cruise_wing_S_Surf1
   solid 3_slat_out_track_S_Surf0
   solid 4_slat_out_track_S_Surf1
   ```

   OpenVSP names each solid `<number>_<component>_S_Surf<k>`. The `slat_out_track` pieces are
   not wing. For `cruise_wing`, `Surf0` and `Surf1` are the two mirror halves. For
   `slat_out_track` they are two separate tracks, both on the positive-y side, so do not assume
   `Surf<k>` means a side. Select by component name first, then pick the wing half by its
   bounding box (item 3 of this list).

2. Split the file into one STL per solid:

   ```bash
   surfaceSplitByPatch cruise_wing.stl
   ```

   Expected: one line per solid, such as
   `Writing patch 1_cruise_wing_S_Surf0 to file "cruise_wing_1_cruise_wing_S_Surf0.stl"`.

3. Keep the wing half on the positive-y side. Check each wing piece's bounding box:

   ```bash
   surfaceCheck cruise_wing_1_cruise_wing_S_Surf0.stl | grep -E "Bounding Box|Surface is"
   ```

   The piece to keep has y running from about 0 to about +6 (feet). On the project's export
   that was `Surf0`: `(3.86 -0.0119906 -0.573113) (7.66486 5.98885 0.04382)`. The two slat-track
   pieces were open 16-triangle sheets: the whole export failed `surfaceCheck` with
   `connected to one face : 32`, and all 32 of those edges were on the slat tracks.

4. Scale it to metres (Section 6) and run `surfaceCheck` again.

What to do next: get the `.vsp3` from the design study repository (Step 1) and use the route
in Section 5.

`scripts/export_vsp_stl.py` is the project's tool for making such exports from a `.vsp3`. It
needs OpenVSP's own Python module (`openvsp`; the registry records OpenVSP 3.51.2), which is
not part of this repository, and the launcher its own instructions use, `scripts/bin/vsppython`,
is not shipped either. Its options were read from the script, but no export has been run with
it, so it is **NOT VERIFIED** here. It is not part of the meshing route: its own header still
calls it the production route, but the registry records it as a verification tool only,
because its exports have a sharp trailing edge.

---

## 8. Putting the surface into a case

You do this in [03](./03_build_a_case.md), not here. Step 3 of that guide creates the case
folder and sets `$CASE`, and its Step 5 copies the surface in and checks it. `$CASE` is not set
in this guide, so do not run case commands from here. The rules that matter:

1. The mesh recipes read the surface as `constant/triSurface/wing.stl`. The file name must be
   `wing.stl`, whatever your file was called.
2. After copying, its checksum must still be the one from Step 9; copying does not change it.
3. The copy inside the case is the file that will be meshed, so that is the copy to run
   `surfaceCheck` on (03 does this).

**Do not create `constant/geometry/`.** OpenFOAM-12 looks for surfaces in `constant/geometry/`
if that folder exists and only falls back to `constant/triSurface/` if it does not. This was
tested: with an empty `constant/geometry/` present, `surfaceFeatures` stopped with
`Cannot read ".../constant/geometry/wing.stl"`. The same message, naming
`constant/geometry/`, also appears when neither folder exists, that is, when the surface was
never copied into the case (also tested). Worse, the project's surface gates
(`scripts/assert_surface_closed.py` and the mesh job's built-in gate in `hpc/mesh_hpc12.pbs`)
only look in `constant/triSurface/`. With a surface in both folders, OpenFOAM would mesh the
one in `constant/geometry/` while the gates checked a different file.

When you reach feature extraction in [03](./03_build_a_case.md), the output of `surfaceFeatures`
is one more input check. With the CR recipe's `surfaceFeaturesDict` and the baseline surface it
reported `open edges : 0` and `points : 1273`. The 1273 is the number of points on the
feature-edge mesh it writes (`wing.eMesh`); the separate `feature points` line reads `4`, so do
not confuse the two. The recipe notes record 1273 `wing.eMesh` points for B and 2067 for W.
W's extra points come from `open edges : 796`, the unwelded root seam of the surface its mesh
case held (Step 9): an `open edges` count above 0 means the surface is not the registered,
welded one. Stop and replace it.

Do not rotate the surface to set the angle of attack. The surface stays as built; the project
sets the angle of attack by rotating the inflow direction (see
[06](./06_trim_to_target_cl.md)). The same surface is used for every flow condition: no geometry
in `data/rans_forces.json` carries two different `geometry_sha256` values across its cases
(three early-cruise baseline records carry none).

---

## 9. Scripts that look relevant but are not for this

1. `scripts/regenerate_oml.py` regenerates the project's older, superseded surfaces (fan tip
   cap by default, root extension added, but never welded) from folders that are not in this
   repository. It does not produce the meshing surfaces.
2. `scripts/export_vsp_stl.py`: see Section 7.
3. The scripts print internal record references in brackets in some messages and help texts
   (for example after a gate name). The records they refer to are not part of this repository;
   the messages themselves are complete without them.

---

## 10. The project's own surface checks (extra, not a replacement for `surfaceCheck`)

These scripts are in `scripts/`. They are useful, but **none of them replaces Step 8**.
`assert_surface_closed.py` and the cluster mesh job's built-in gate (in `hpc/mesh_hpc12.pbs`)
round coordinates to the nearest micrometre before matching points, and both passed the
closed-but-unwelded surface that `surfaceCheck` rejects (Step 8). Only `surfaceCheck` and gate
G6 of `check_surface_closed.py` see that seam.

### 10a. `assert_surface_closed.py`: the case-level gate

It checks every `constant/triSurface/*.stl` in a case folder: closed, consistent orientation,
normals pointing outward, and with `--require-overhang` also that the surface reaches at least
10 mm below `y = 0`. Always pass `--require-overhang`; the root extension is part of a correct
surface.

First prove the checker can fail (four synthetic surfaces, one that must pass and three that
must fail for different reasons):

```bash
python3 $REPO/scripts/assert_surface_closed.py --self-test
```

Expected: four lines starting `[ok ]`. The summary line then says `5/5 fixtures behaved as
specified`; in this repository only four fixtures run, because the fifth needs a project file
that is not shipped. Four `[ok ]` lines is the correct result.

Then run it on your surface. The script wants a case folder, so make a small test folder that
holds only the surface in the place a case keeps it (you will run the same check on the real
case in [03](./03_build_a_case.md), Step 6):

```bash
mkdir -p $GEOM/B/gatecase/constant/triSurface
cp $GEOM/B/baseline_oml_placed.stl $GEOM/B/gatecase/constant/triSurface/wing.stl
python3 $REPO/scripts/assert_surface_closed.py $GEOM/B/gatecase --require-overhang
echo "exit status $?"
```

Expected: `<path>/gatecase: all surfaces closed, outward, and crossing y=0`, then
`exit status 0`. On the surface from Step 5 (never closed) it printed `REFUSED` and
`constant/triSurface/wing.stl: OPEN: 398 edges with one neighbour (vertices unified at 1e-6 m)`.
On the OpenVSP export from Section 7 it printed
`DOES NOT CROSS THE SYMMETRY PLANE: y_min -3.655 mm, needs < -10.0 mm`. On a case folder with no
surface it refuses rather than passing.

Do not use its `--audit-all` option here: it looks for a `cases/` folder that this repository
does not have, checks zero surfaces, and still reports success.

### 10b. `check_surface_closed.py`: the file-level gate, gate by gate

```bash
python3 $REPO/scripts/check_surface_closed.py --glob "$GEOM/B/baseline_oml_placed.stl"
```

It reports gates G1 to G6 for each file. **In this repository its overall verdict is always
`FAIL`** and it exits with status 1, for two reasons that are not about your surface: gate G5
compares against a manifest file that is not shipped, and it also checks for two retired
project files that are not shipped (`RETIREMENT BROKEN` lines, then `PARTITION FAILURE`).
So read the per-file line instead:

1. Good surface: the only failed gate listed is `G5(NO PRE-FIX MANIFEST ENTRY, ...)`.
2. Closed but unwelded surface: `G6(796 open edges exact vs 0 at 6 dp  <-- CLOSED TO THIS GATE,
   OPEN TO snappyHexMesh)` appears as well. G6 is the one Python check that sees the unwelded
   seam.

Do not use `--snapshot` (it writes into a project folder that does not exist here and stops
with `FileNotFoundError`). `--negative-control` only works on files written by Steps 4 to 6; on a
welded file from Step 7 it wrongly reports `MISSED` because it cannot find a triangle to
remove, so do not read that as a fault in your surface.

### 10c. `check_oml_gates.py`: is the placement applied?

This measures dihedral, incidence, chords, sweep and tip twist from the triangles themselves.

```bash
python3 $REPO/scripts/check_oml_gates.py $GEOM/B/baseline_oml_placed.stl
```

Expected for the finished baseline (the first lines of the output; exit status 0):

```text
TE reference: midpoint
GATE                       MEASURED       TARGET        DELTA      TOL
1 DIHEDRAL PRESENT         0.130494     0.139000    -0.008506    0.010  PASS  (m)
2 INCIDENCE PRESENT        1.852651     1.873000    -0.020349    0.250  PASS  (deg)
3 ROOT CHORD               0.639125     0.640081    -0.000956    0.002  PASS  (m)
3 TIP CHORD                0.152151     0.152400    -0.000249    0.002  PASS  (m)
4 LE SWEEP                28.844652    28.861990    -0.017338    0.250  PASS  (deg)
4 TIP TWIST               -4.122956    -4.114950    -0.008006    0.250  PASS  (deg)

OVERALL: PASS
```

After these lines it prints further values under `reported, not gated`, which are for
information only. Gate 1 measures how far the tip quarter-chord point rises above the root
quarter-chord point, in metres (target 0.139 m), not the dihedral angle in degrees.

On the unplaced loft from Step 4 it fails exactly the two gates a planform check cannot see,
while the chords, sweep and twist still pass (exit status 1):

```text
TE reference: midpoint
GATE                       MEASURED       TARGET        DELTA      TOL
1 DIHEDRAL PRESENT        -0.000037     0.139000    -0.139037    0.010  FAIL  (m)
2 INCIDENCE PRESENT       -0.020060     1.873000    -1.893060    0.250  FAIL  (deg)
```

The targets are the baseline's. For W, the twist concept, gate `4 TIP TWIST` fails by design:
it measured -5.127872 deg against the baseline's -4.114950 deg, a difference of -1.013 deg,
which is W's -1.00 deg tip twist. C, M, H and F pass all six gates.

---

## 11. The known-answer test, and all six build commands

To confirm that your tools produce the project's surfaces exactly, rebuild one and compare its
checksum with the registry. All six have been rebuilt from their `.vsp3` files with the
scripts in this repository and the commands below, and **all six reproduced the registered
sha256 checksums exactly** (Python 3.11.8, numpy 1.26.4, scipy 1.11.4).
They were built in a terminal with no OpenFOAM loaded, so Step 7 found OpenFOAM-12 by itself
(Step 7, item 1).

Each block is Steps 4 to 7. It needs `$REPO` and `$GEOM` set as in Section 4 and the `.vsp3`
in `$GEOM/vsp3/`, and it creates and enters `$GEOM/<letter>` itself. The output is named
`<ID>_oml_placed.stl`, not the registry's `..._laddercap.stl` name; the names differ, and only
the checksum says whether the content is the same.

**Variant A: B, C and M.** For B use `L=B; ID=baseline_wing_only_refined` as shown; for C use
`L=C; ID=cte_i002_c04`; for M use `L=M; ID=mcv2_i002_c01`.

```bash
L=B; ID=baseline_wing_only_refined
mkdir -p $GEOM/$L && cd $GEOM/$L
python3 $REPO/scripts/vsp3_to_stl.py ../vsp3/$ID.vsp3 --geom cruise_wing --out ${ID}_oml.stl \
    --nchord 200 --nspan 240 --blunt-te 0.0063 --cap delaunay --cap-edge 0.0030 > step1.json
python3 $REPO/scripts/apply_placement_transform.py --loft ${ID}_oml.stl --out ${ID}_oml_placed.stl --order xy
python3 $REPO/scripts/close_wing_surface.py ${ID}_oml_placed.stl --write
python3 $REPO/scripts/weld_surface.py ${ID}_oml_placed.stl --in-place
sha256sum ${ID}_oml_placed.stl
```

**Variant B: H, F and W.** Use `L=H; ID=chc_g02_c06` for H, `L=F; ID=cfft_b02_c01` for F,
`L=W; ID=cffw_b01_c01` for W.

```bash
L=H; ID=chc_g02_c06
mkdir -p $GEOM/$L && cd $GEOM/$L
python3 $REPO/scripts/vsp3_to_stl.py ../vsp3/$ID.vsp3 --geom cruise_wing --out ${ID}_oml.stl \
    --nchord 200 --nspan 240 --blunt-te 0.0063 --cap delaunay > step1.json
python3 $REPO/scripts/apply_placement_transform.py --loft ${ID}_oml.stl --out ${ID}_oml_placed.stl --order xy
python3 $REPO/scripts/close_wing_surface.py ${ID}_oml_placed.stl --write
python3 $REPO/scripts/weld_surface.py ${ID}_oml_placed.stl --tol 1e-6 --in-place
sha256sum ${ID}_oml_placed.stl
```

Compare the printed checksum with the list in Section 3. Then run Step 8 on the result. On all
six, `surfaceCheck` reported a closed surface, one piece, one zone and no illegal triangles, and
`assert_surface_closed.py --require-overhang` passed.

Three notes on the variants:

1. The route text stored in `registry/candidates.yaml` for H, F and W leaves out two things:
   `--write` on `close_wing_surface.py`, and the weld step. Followed literally, it leaves the
   root open, because without `--write` Step 6 only reports. With `--write` added but still no
   weld, it gives the earlier, unwelded versions of those files, not the registered ones. The
   commands above are the ones that reproduce the registered checksums.
2. With the default weld tolerance instead of `--tol 1e-6`, H, F and W come out with different
   checksums (the weld then merges more points and drops 4 to 16 triangles). Use `--tol 1e-6`
   for them.
3. For a geometry that is not in the table, there is no registered answer and the project has not
   recorded which variant to use. Use variant A, which is the baseline's, unless the project tells
   you otherwise, and record which you used. Know what that choice includes: on B, C and M,
   variant A's default weld tolerance also merged two pairs of points about 10 micrometres apart
   near the tip trailing edge and dropped 4 triangles (Step 7), while variant B welds only the
   root seam. Check the weld output of your own geometry for the same thing. Build every
   geometry you intend to compare by the same route.

If your checksum differs, first run Step 8. A surface that passes Step 8 but has a different
checksum is still not a registered geometry: note your Python, numpy and scipy versions and ask
the project before using it for results.

---

## 12. Troubleshooting

| What you see | Cause | What to do |
|---|---|---|
| Step 1c prints nothing, or some of the six names never appear | wrong folder, or the design study repository does not hold those files with the registered content | check the folder you cloned into; otherwise ask the project (Step 1c) |
| Step 1c prints more than six lines | the clone holds identical copies of a file in more than one folder | harmless; check with Step 1d |
| `FAILED` from `sha256sum -c` in Step 1d | not the file the project built from (possibly another version under the same name) | do not use it; ask the project which file to use |
| fewer than six `OK` lines in Step 1d | a file is missing, or is in the folder under another name | Step 1d, item 2 |
| `ValueError: ... need exactly one wing geom, found [...]` | the `.vsp3` holds more than one wing component | run Step 2's `--list` and pass `--geom <wing name>` |
| `half_span_m` about 6.0 in Step 4 | built in feet (`--scale 1.0` was passed) | rerun Step 4 without `--scale` |
| bounding box maximum y about 6 | the surface is in feet | Section 6 |
| bounding box maximum y about 0.56 | scaled twice | start again from the unscaled file |
| `--order auto requires --reference` | `--order xy` missing in Step 5 | add `--order xy` |
| `RuntimeError: loft root LE is at [1.176528, 0.0, -0.12192], not the origin` | Step 5 run on an already placed file | use the Step 4 output as `--loft` |
| `already closed` in Step 6 | the input was already closed | if you already ran Step 6 on this file, this is expected: go on to Step 7. Otherwise a different file was passed (the Step 5 output is always open); check the name. For OpenVSP exports see Section 7 |
| `FAILED: surfaceCheck could not be run, so nothing was checked and nothing was promoted.` in Step 7 (exit status 2) | OpenFOAM-12 not found where the script looked, `ARGUS_OF_BASHRC` wrong or pointing at another OpenFOAM version, or a file name OpenFOAM refuses (for example one with a space); the `reason:` line says which | Step 7, items 1 and 2 |
| `NOT PROMOTED` in Step 7 with real results in `after:` (exit status 1) | the surface is still not closed after welding | if Step 6 was skipped, go back to it; otherwise rerun from Step 4 and check the flags against Section 11 |
| `surfaceCheck`: `connected to one face : 398` | root never closed | Step 6 |
| `surfaceCheck`: `connected to one face : 796`, 2 zones, `More than one normal orientation.` | closed but not welded | Step 7 |
| `surfaceCheck`: `Surface is self-intersecting at 216 locations.` together with `connected to one face : 796` | closed but not welded | Step 7 (it goes away after welding) |
| `<file>: FAILED` from `sha256sum -c -` in Step 9 | the surface is not the registered one (for example Step 7 skipped, or a flag differs) | Step 9, "If it does not" |
| bounding box minimum y about -0.0037 | not closed | Step 6 |
| bounding box minimum x about 0 | not placed | Step 5 |
| `Cannot read ".../constant/geometry/wing.stl"` from `surfaceFeatures` | a `constant/geometry/` folder exists, or the surface was never copied into `constant/triSurface/` | remove `constant/geometry/` if it exists; make sure the surface is at `constant/triSurface/wing.stl` (Section 8) |
| `no file was verified` from `sha256sum -c` | none of your files has a name in the list | Step 1d, item 2 (for `.vsp3` files); in Step 9, use the name-independent search |
| checksum not found in the registry | not a registered geometry | Section 11, last paragraph |
| `foamVersion` does not print `OpenFOAM-12`, `foamRun` not found, or the `surfaceCheck` log's `Build` line does not start with `12-` | your shell has loaded a different OpenFOAM version | Section 4 |
| `surfaceCheck` not found | no OpenFOAM loaded | Section 4 |
