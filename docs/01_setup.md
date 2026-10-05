# 01. Setting up your environment

> **At a glance**
> 1. **What this guide gets done:** it installs and checks OpenFOAM-org 12, Open MPI, ParaView, Python and LaTeX on your local Linux machine, and checks OpenFOAM-org 12 and Python on the HPC12 cluster.
> 2. **Before you start:** a clone of this repository (`argus-production`), a local Linux machine where you can use `sudo` (the install commands are for Ubuntu), and an HPC12 account. This is the first guide; the overview is [../README.md](../README.md).
> 3. **At the end:** every line of the checklist at the bottom of this page matches.
> 4. **How long:** the local smoke tests take seconds. The OpenFOAM-12 build on HPC12 runs for hours, as a batch job.

This guide prepares two machines:

1. **Your local Linux machine.** You prepare and check cases here, look at results and build the report.
2. **The HPC12 cluster.** The full-size meshes are built and solved here.

**OpenFOAM-org 12 is the only supported OpenFOAM version for every guide in this repository.** Every OpenFOAM command in these guides (`foamRun`, `blockMesh`, `surfaceFeatures`, `snappyHexMesh`, `checkMesh`, `surfaceCheck`, `foamDictionary`, `decomposePar`, ...) is an OpenFOAM-org 12 command. Some recipe files carry older version labels in their header banners (`Version:  7`, `v1912` or `v2006`). Those banners sit inside a comment and select nothing.

> **NO GEOMETRY IS SHIPPED.** This repository has no wing surface: no STL file and no `.vsp3` file. The wing geometry comes as `.vsp3` files: F, H and W from the design study's public GitHub repository, https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow (folder `results/geometries/low_speed/`); B, C and M, which that repository does not hold, from the ARGUS project drive, `\\tudelft.net\staff-umbrella\ARGUS Morphing wing CFD\argus_hpc12_backup\geometry_vsp3\`, which holds all six with a `SHA256SUMS` file. [02_geometry_from_openvsp.md](./02_geometry_from_openvsp.md) explains how to check those files and build every meshing STL from them yourself. You do not need geometry for anything in this guide.

Back to the overview: [../README.md](../README.md).

## Two conventions used in every guide

`$REPO` is your clone of this repository. The repository is kept on GitHub as `argus-production`, a private repository. Cloning it needs a GitHub account that the project has given access to; the project also gives you the clone address and tells you how to sign in. Clone it once:

```bash
git clone <argus-production clone URL> $HOME/argus-production   # NOT VERIFIED: the address is not yet available
```

The first line it prints starts with `Cloning into` and names the folder. Keep the full path of the repository free of spaces and commas, here and on HPC12: on HPC12 the paths of its scripts go into the job's `qsub -v` list, which cannot carry either.

Then set `$REPO` once in each new terminal:

```bash
export REPO=$HOME/argus-production   # or wherever you cloned it
```

On HPC12, [04_send_to_hpc.md](./04_send_to_hpc.md), Step 4, copies the repository into your cluster home. There too, set `$REPO` to that copy, as a full path starting with `/` (as `$HOME/...` does), in each new login shell.

`$CASE` is the case directory you are building or running. [03_build_a_case.md](./03_build_a_case.md) introduces it.

Unless a guide says otherwise, paths are relative to the repository root.

## What you will have at the end

| Item | Where | Version the project used |
|---|---|---|
| OpenFOAM-org 12 | local, `/opt/openfoam12` | Ubuntu package `openfoam12`, version `20250206` |
| ParaView (for looking at results) | local, `/opt/paraviewopenfoam510` | 5.10.1, installed together with OpenFOAM |
| Open MPI (for parallel runs) | local, `/usr/bin/mpirun` | 4.0.3 |
| Python 3 with six packages | local | Python 3.11.8 |
| pdflatex and latexmk (for the report) | local | TeX Live 2024 |
| OpenFOAM-org 12 | HPC12, `$HOME/OpenFOAM/OpenFOAM-12` | built from source; reports `Build: 12` |
| Python 3 (the system one; nothing to install) | HPC12, `/usr/bin/python3` | 3.6.8, with numpy 1.12.1 and no jsonschema |

The project's local machine ran Ubuntu 20.04.6 LTS under WSL2 on Windows, with 12 physical cores (24 hardware threads, which is what `nproc` counts) and 23 GB of memory visible to Linux. That machine cannot build or solve the production meshes. For example, the wake-refined (r2) baseline (B) mesh at Condition CR has 99,111,506 cells (see `data/run_cards/SWB_r2_trim.json`, field `mesh.cells`). Those meshes are built and solved on HPC12.

---

## Part A. Your local machine

### A1. Install OpenFOAM-org 12

The OpenFOAM Foundation publishes OpenFOAM-org 12 as an Ubuntu package. The project's machine was set up this way. Its apt configuration has three parts:

1. The Foundation's package repository, as the line `deb http://dl.openfoam.org/ubuntu focal main` in `/etc/apt/sources.list`. `focal` is the codename of Ubuntu 20.04.
2. The OpenFOAM Foundation's signing key, in apt's keyring.
3. The package itself, installed with `apt install openfoam12`. That command also installed ParaView (`paraviewopenfoam510`) automatically.

**Step 1. Add the repository and its signing key.** Follow the OpenFOAM Foundation's download page for your Ubuntu release: https://openfoam.org/download (this address appears in OpenFOAM-12's own `bin/paraFoam` script; it could not be opened while this guide was written). **NOT VERIFIED:** this step was done on the project's machine before this guide was written, and the exact commands were not recorded, so they are not repeated here.

**Step 2. Install the package.**

```bash
sudo apt update
sudo apt install openfoam12
```

**Step 3. Confirm that it installed.**

```bash
dpkg -l openfoam12 | tail -1
```

Expected output (the project's machine):

```
ii  openfoam12     20250206     amd64        OpenFOAM is the leading free, open source software for
```

The `ii` at the start means "installed". If you see `dpkg-query: no packages found matching openfoam12` instead, the install did not happen. Check that apt knows about the Foundation's repository:

```bash
grep -h "dl.openfoam.org" /etc/apt/sources.list /etc/apt/sources.list.d/*.list
```

On the project's machine this prints `deb http://dl.openfoam.org/ubuntu focal main` (plus a commented-out `# deb-src` line). If it prints nothing, Step 1 did not take effect.

If you are not on Ubuntu, the same download page lists other routes. None of them was tested for this guide.

### A2. Load OpenFOAM-12 into your terminal

OpenFOAM only works in a terminal that has loaded ("sourced") its environment file:

```bash
source /opt/openfoam12/etc/bashrc
```

If it works, it usually prints nothing. If it prints `x86_64-conda-linux-gnu-cc: error: unrecognized command-line option '--showme:link'`, an Anaconda installation is ahead of the system tools on your `PATH`: OpenFOAM still loaded, but read A5 before you run anything in parallel.

You need to do this in every new terminal. To make it happen automatically, first check whether your `~/.bashrc` already loads a different OpenFOAM:

```bash
grep -n -i openfoam ~/.bashrc
```

If a line loads a different version (for example `source /opt/openfoam7/etc/bashrc`), open `~/.bashrc` in a text editor and put a `#` at the start of that line. Then add the OpenFOAM-12 line at the end of the file:

```bash
echo 'source /opt/openfoam12/etc/bashrc' >> ~/.bashrc
```

If you use Anaconda, its set-up block sits higher up in `~/.bashrc`, so every new terminal will now print the `--showme:link` message above. That is the MPI trap in A5.

Close the terminal and open a new one. A new terminal does not remember `$REPO`, so set it again (see "Two conventions" above). Then continue with A3.

### A3. Check the version

```bash
foamVersion
echo $WM_PROJECT_VERSION
which foamRun
foamRun -help | grep -E "^(Using|Build):"
```

Expected output (the project's machine):

```
OpenFOAM-12
12
/opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/foamRun
Using: OpenFOAM-12 (see https://openfoam.org)
Build: 12-86e126a7bc4d
```

Things that can look wrong but are not:

1. `foamVersion` is a shell function that the OpenFOAM environment file defines, not a program. So `which foamVersion` prints nothing even on a correct setup. Use `foamVersion` itself, or `echo $WM_PROJECT_VERSION`.
2. The `Build:` line identifies the exact source build. `12-86e126a7bc4d` is what package version `20250206` reports. If your package version differs, your `Build:` line may differ too; the `Using:` line must say `OpenFOAM-12`. (The build on HPC12 prints `Build: 12`; see B4.)
3. Do not use `which simpleFoam` as a version test. OpenFOAM-12 still ships a program called `simpleFoam`: it prints a notice that it has been replaced and then runs `foamRun -solver incompressibleFluid` itself. So `simpleFoam` exists under both versions and says nothing about which one you loaded.

Now check the utilities the guides use:

```bash
which blockMesh surfaceFeatures surfaceCheck snappyHexMesh checkMesh foamDictionary decomposePar reconstructPar
```

Expected: eight lines, each starting with `/opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/`. A missing line means that utility is not on your `PATH`; go back to A2.

### A4. Trap: a terminal that loaded a different OpenFOAM

**Symptom.** You type an OpenFOAM-12 command and get:

```
foamRun: command not found
```

In a script, the same message is preceded by the script's name and line number, for example `./myjob.sh: line 2: foamRun: command not found`. Also, `foamVersion` prints another version (for example `OpenFOAM-7`), or `echo $WM_PROJECT_VERSION` prints something other than `12`. This happened on the project's machine, where `~/.bashrc` loaded an older version. The same `command not found` message appears in a terminal that has loaded no OpenFOAM at all.

**Fix.** Run `source /opt/openfoam12/etc/bashrc` in that terminal. It replaces the other version's settings. This was checked on a machine with an older version installed beside OpenFOAM-12: after loading OpenFOAM-12 on top of the older one, none of the older version's directories was left on `PATH` or `LD_LIBRARY_PATH`. Then fix `~/.bashrc` as described in A2 so that new terminals start correctly.

### A5. Trap: Anaconda's MPI hiding the system Open MPI

The OpenFOAM-12 package runs in parallel through the system's Open MPI (`/usr/bin/mpirun`, version 4.0.3). Anaconda ships its own `mpirun` (MPICH). If Anaconda's directory comes first on your `PATH`, parallel OpenFOAM runs break. This was reproduced on the project's machine.

**Symptom 1**, printed while you source OpenFOAM:

```
x86_64-conda-linux-gnu-cc: error: unrecognized command-line option '--showme:link'
```

**Check which `mpirun` you have:**

```bash
mpirun --version | head -1
```

Good: `mpirun (Open MPI) 4.0.3`. Bad: `HYDRA build details:`, which is MPICH.

**Symptom 2**, in a parallel run. Every rank prints:

```
--> FOAM FATAL ERROR:
bool IPstream::init(int& argc, char**& argv) : attempt to run parallel on 1 processor
```

and the run then either stops within seconds or hangs. Both were seen on the project's machine: one test was still hanging after 120 seconds, and another returned after about 1 second. If it hangs, press Ctrl+C. Either way, look for the message in the log with the `grep` shown in A6.

**Fix.** Put the system directory first, then load OpenFOAM again:

```bash
export PATH=/usr/bin:$PATH
source /opt/openfoam12/etc/bashrc
mpirun --version | head -1
```

The last line should now print `mpirun (Open MPI) 4.0.3`, and the `--showme:link` message no longer appears. This also means that `python3` in that terminal is now the system Python, which lacks the packages from A8.

The problem only affects parallel runs. With Anaconda first on `PATH` and OpenFOAM-12 loaded, serial utilities such as `surfaceCheck` still worked in the test, and so did the Python scripts that call them. A practical arrangement is two terminals, which the later guides refer to by these names:

1. **Your OpenFOAM terminal**, for parallel OpenFOAM runs: a terminal with the fix above.
2. **Your Python terminal**, for the Python scripts: a terminal with your Python environment (A8) first on `PATH` and OpenFOAM-12 also loaded. Ignore the `--showme:link` message there, and do not start parallel runs from it. Scripts that need numpy, such as the geometry scripts and `scripts/install_postpro.py`, run here.

### A6. Smoke test: run an OpenFOAM tutorial, serial and parallel

This proves that the solver runs and that parallel runs work before you touch a real case. It takes seconds. It works in a folder outside the repository.

Serial run:

```bash
source /opt/openfoam12/etc/bashrc
mkdir -p $HOME/of12-smoke && cd $HOME/of12-smoke
cp -r $FOAM_TUTORIALS/incompressibleFluid/cavity cavity-serial
cd cavity-serial
blockMesh > log.blockMesh 2>&1
foamRun > log.foamRun 2>&1
tail -3 log.foamRun
```

Expected: the word `End` (between two blank lines). If you see `FOAM FATAL ERROR` instead, read the lines under it in `log.foamRun`. A `command not found` means you are in the situation of A4.

Parallel run on 2 cores, in a fresh copy:

```bash
cd $HOME/of12-smoke
cp -r $FOAM_TUTORIALS/incompressibleFluid/cavity cavity-parallel
cd cavity-parallel
cp $FOAM_ETC/caseDicts/preProcessing/decomposeParDict system/
foamDictionary -writePrecision 12 -entry numberOfSubdomains -set 2 system/decomposeParDict
foamDictionary -writePrecision 12 -entry method -set scotch system/decomposeParDict
blockMesh > log.blockMesh 2>&1
decomposePar > log.decomposePar 2>&1
mpirun -np 2 foamRun -parallel > log.foamRun 2>&1
grep "^nProcs" log.foamRun
tail -3 log.foamRun
```

Expected: the two `foamDictionary` lines print `New entry numberOfSubdomains 2;` and `New entry method          scotch;`, the `grep` prints `nProcs : 2`, and the last lines of the log are `End` and `Finalising parallel run`. On the project's machine the `mpirun` line took 5 to 7 seconds.

**`foamDictionary` reads every file name relative to the folder you are in, even one that starts with `/`.** So run it from inside the case, with a relative name such as `system/decomposeParDict`, as above. Given a full path, it stops with `FOAM FATAL IO ERROR` and `file "<the folder you are in>/<the full path>" does not exist`. This was reproduced here with this OpenFOAM-12 package, and the same failure has been seen on HPC12.

If the `mpirun` line has not returned after a minute, it has probably hung. Because its output goes into `log.foamRun`, a hang shows nothing on the screen. Press Ctrl+C. If it returned but the first `grep` printed no `nProcs : 2`, something also went wrong. In either case run:

```bash
grep "attempt to run parallel" log.foamRun
```

If that prints anything, you have the MPI trap: go to A5.

You can delete `$HOME/of12-smoke` afterwards.

### A7. ParaView, for looking at meshes and results

ParaView 5.10.1 was installed with the `openfoam12` package (package `paraviewopenfoam510`). Loading OpenFOAM-12 puts it on your `PATH`:

```bash
pvbatch --version
paraview --version
```

Expected: `paraview version 5.10.1` from each. If `pvbatch` is not found, ParaView was not installed: the `openfoam12` package only *recommends* `paraviewopenfoam510`, so an apt set up to skip recommended packages leaves it out. Install it the same way as OpenFOAM:

```bash
sudo apt install paraviewopenfoam510
```

`pvbatch` works without a screen. `paraview` is the graphical program and needs a display. On the project's machine (Windows with WSL2) Windows provided the display. Without a display, `paraview` fails with `qt.qpa.xcb: could not connect to display`.

To open a case, go to the case directory and run `paraFoam`. Or create a small file that ParaView can open, and open it from ParaView's File menu:

```bash
paraFoam -touch
```

This prints `Created '<case-name>.OpenFOAM'`. `.OpenFOAM` files are ignored by git in this repository.

Two notes:

1. `scripts/trefftz.py` is written to run under `pvbatch`, not `python3`. Its own header says so.
2. The project did not use ParaView on HPC12. When OpenFOAM was built there, its ParaView reader did not build, because the cluster lacks the VTK headers it needs. Look at results on your local machine.

### A8. Python

The scripts in `scripts/`, `report/scripts/` and `report/fig/` (117 files) need, on your own machine, **Python 3.11**: the project ran them with Python 3.11.8. The package metadata set the lowest version at 3.10 (duckdb 1.5.5 requires 3.10 or newer; numpy, scipy and matplotlib at the pinned versions require 3.9 or newer). No Python other than 3.11 was tried with the pinned packages: NOT VERIFIED. Ubuntu 20.04's own `python3` (3.8.10) is not enough, for three further reasons:

1. It cannot parse `scripts/render_delta_maps_of12.py`, which uses syntax added in Python 3.9.
2. It has no numpy, scipy or matplotlib, and no pip.
3. Its jsonschema (3.2.0) is too old for `scripts/validate_run_card.py`, which needs `jsonschema.Draft202012Validator`.

**Python on HPC12 is not this Python.** The HPC12 login node's `python3` is 3.6.8, with numpy 1.12.1 and no jsonschema (checked at an HPC12 login shell). That is enough for `scripts/assert_case_postpro.py`, the post-processing check you run there before a solve, which needs only the standard library: it gave `8 of 8` there on each of the six delivered Condition CR solve cases. If a case reaches HPC12 without its post-processing block, `scripts/install_postpro.py` can be run there as a fallback (its `--dry-run` ran there with numpy 1.12.1, but its writing step has not been run on HPC12: NOT VERIFIED). Everything else in this section, including `scripts/validate_run_card.py`, runs on your own machine.

Every `import` in those 117 files was read, including imports inside functions and imports of other scripts in the repository. These are the third-party packages they need:

| Package to install | Import name | Scripts that need it (of 117) | Version tested |
|---|---|---|---|
| `numpy` | `numpy` | 77, including the geometry chain `vsp3_to_stl.py`, `apply_placement_transform.py`, `close_wing_surface.py` and `weld_surface.py`, and `install_postpro.py` | 1.26.4 |
| `matplotlib` | `matplotlib` | 46: all 29 `deck_*.py` and `make_*.py` in `report/fig/`, `report/scripts/build_span_efficiency.py`, and 16 plotting scripts in `scripts/` | 3.8.0 |
| `scipy` | `scipy` | 16, including `vsp3_to_stl.py` and `weld_surface.py` | 1.11.4 |
| `PyYAML` | `yaml` | 4: `validate_run_card.py`, `build_2d_case.py`, `regenerate_oml.py`, `build_results_db.py` | 6.0.1 |
| `jsonschema` | `jsonschema` | 1: `validate_run_card.py` | 4.19.2 |
| `duckdb` | `duckdb` | 1: `build_results_db.py` | 1.5.5 |

How the 117 files split up:

1. **83 files** need at least one of the six packages in the table.
2. **32 files** need only the standard library. They include `scripts/assert_case_postpro.py`, `scripts/verify_recipe.py`, `scripts/derive_run_card_fields.py`, `scripts/harvest_rans_forces.py` and most of `report/scripts/`.
3. **2 files** need, besides the standard library, only a package that is **not** installed with pip:
   1. `scripts/export_vsp_stl.py` imports `openvsp`. The project's `.vsp3`-to-STL route, `scripts/vsp3_to_stl.py`, reads the `.vsp3` file itself and does not need OpenVSP.
   2. `scripts/trefftz.py` imports `paraview`. It runs under `pvbatch` (A7).

**Seven files cannot run to the end from this repository, whatever you install.** They import a module that is not in this repository, or look for `scripts/` in the wrong folder, and stop with `ModuleNotFoundError` when they reach that import. This is not a fault in your Python set-up.

1. At start-up: `scripts/check_synthetic_gates.py` (missing `make_synthetic_wing`), `scripts/collect_2d_grid.py` and `scripts/collect_trim_deltas.py` (missing `collect_results`), and `report/scripts/build_span_efficiency.py` (looks for `case_derived_quantities` outside the repository).
2. Part-way through a run, if they get that far: `scripts/surface_integrate.py` (missing `m6_pressure_excess`; its `--self-test` option does not reach that import and passes), `report/fig/make_spanwise.py` and `report/fig/make_hinge.py` (both look for `scripts/` outside the repository).

None of these guides runs the four files in `scripts/`. [07_postprocessing_and_report.md](./07_postprocessing_and_report.md) covers the three report files.

**Install.** You need a Python 3.11 that is not Ubuntu 20.04's own. Ubuntu 20.04's standard package sources offer only `python3.8` and `python3.9` (checked with `apt-cache pkgnames`), so it has to come from elsewhere. This guide does not cover installing one. The project used an Anaconda installation with Python 3.11.8; if you use Anaconda too, read A5 first. If you already have `conda`, one way to get a separate environment with the tested Python version is:

```bash
conda create -n argus python=3.11
conda activate argus
```

Then check that `python3` is the new one: `python3 --version` must print `Python 3.11.` and a number, not Ubuntu's `Python 3.8.10`. `python3 -m pip` installs into whichever `python3` comes first on your `PATH`. Install the six packages at the versions the project tested:

```bash
python3 -m pip install numpy==1.26.4 scipy==1.11.4 matplotlib==3.8.0 PyYAML==6.0.1 jsonschema==4.19.2 duckdb==1.5.5
```

The versions are pinned because the bit-for-bit surface reproduction in [02_geometry_from_openvsp.md](./02_geometry_from_openvsp.md) was only tested with numpy 1.26.4 and scipy 1.11.4.

**NOT VERIFIED:** none of these three install lines was run in a fresh environment for this guide, because they download packages. The `conda create` options were checked with `conda create --help`. The package names and versions were read from the installed package metadata of the environment the project used, and an offline dry run of the `pip` line there (with `--dry-run --no-index` added) accepted it and found every package already at the pinned version.

**Check:**

```bash
python3 --version
python3 -c "import numpy, scipy, matplotlib, yaml, jsonschema, duckdb; print('python packages ok')"
python3 -c "import jsonschema; jsonschema.Draft202012Validator; print('jsonschema ok')"
```

Expected: `Python 3.11.` and a number (the project's was `Python 3.11.8`), then `python packages ok`, then `jsonschema ok`. What failure looks like, from Ubuntu 20.04's own `python3`:

```
ModuleNotFoundError: No module named 'numpy'
AttributeError: module 'jsonschema' has no attribute 'Draft202012Validator'
```

A quick test that a real script starts:

```bash
cd $REPO
python3 scripts/vsp3_to_stl.py --help
```

Expected first line: `usage: vsp3_to_stl.py [-h] [--geom GEOM] --out OUT [--nchord NCHORD]`. If you get `python3: can't open file '<your home>/scripts/vsp3_to_stl.py'`, `$REPO` is not set in this terminal, so `cd $REPO` went to your home folder (see "Two conventions").

A test that the run-card check works end to end. It reads the run-card schema, `scripts/run_card.schema.json`, and needs PyYAML and jsonschema 4:

```bash
cd $REPO
python3 scripts/validate_run_card.py --self-test
python3 scripts/validate_run_card.py data/run_cards/*.json
```

Expected:

```
self-test PASS (example accepted; 4 corrupted variants rejected)
PASS data/run_cards/SWB_r2_trim.json
PASS data/run_cards/SWC_r2_trim.json
PASS data/run_cards/SWF_r2_trim.json
PASS data/run_cards/SWH_r2_trim.json
PASS data/run_cards/SWM_r2_trim.json
PASS data/run_cards/SWW_r2_trim.json
```

If it ends with `AttributeError: module 'jsonschema' has no attribute 'Draft202012Validator'`, the `python3` in this terminal is not the one you installed the packages into.

**`scripts/weld_surface.py` runs `surfaceCheck` itself, and only under OpenFOAM-12.** It uses the OpenFOAM-12 already loaded in your terminal. If none is loaded, it uses the `etc/bashrc` named by the variable `ARGUS_OF_BASHRC`, and otherwise `/opt/openfoam12/etc/bashrc`. It never uses another OpenFOAM version. Under the first `=== <file> ===` line it prints, one line says which one it used (later surfaces do not repeat it), for example:

```
  surfaceCheck: OpenFOAM build 12-86e126a7bc4d, sourced /opt/openfoam12/etc/bashrc (default)
```

In a terminal where OpenFOAM-12 is loaded, the end of that line reads `already loaded (WM_PROJECT_DIR=/opt/openfoam12)` instead. If it cannot run `surfaceCheck` (for example, `ARGUS_OF_BASHRC` names a file that does not exist), it prints `FAILED: surfaceCheck could not be run, so nothing was checked and nothing was promoted.` followed by a `reason:` line, and exits with status 2. [02_geometry_from_openvsp.md](./02_geometry_from_openvsp.md), Step 7, shows it in use. The HPC12 job scripts read a variable with the same name, `ARGUS_OF_BASHRC`, but with a different default (B2).

**Do not run the seven older 2D-study scripts** (`solver_recipe.py`, `recipe_trial.py`, `overnight_recipe_study.py`, `wallresolved_transonic_routes.py`, `trim2d_transonic.py`, `build_rae2822_transonic.py`, `build_rae2822_wolf.py`): they are OpenFOAM-org 7 tools that no guide uses.

### A9. LaTeX, for the report

The report needs `pdflatex` and `latexmk`, plus the 13 LaTeX packages that `report/src/argus_rans_validation_report.tex` loads. The project's machine had TeX Live 2024, installed with the upstream TeX Live installer rather than with apt (`dpkg -l` lists no TeX Live packages). The installer's log, `/usr/local/texlive/2024/install-tl.log`, records how it was run:

```bash
./install-tl --no-interaction
```

The log records that this installed the full scheme into `/usr/local/texlive/2024`. The full scheme is large: on the project's machine `du -sh /usr/local/texlive/2024` prints `8.8G` for that folder, so check your free disk space first. The log ends by asking you to add `/usr/local/texlive/2024/bin/x86_64-linux` to your `PATH`, and the project did so with this line in `~/.profile`:

```bash
export PATH=/usr/local/texlive/2024/bin/x86_64-linux:$PATH
```

**NOT VERIFIED:** how the installer was downloaded and unpacked was not recorded, and the installation was not repeated for this guide. The TeX Live web site named in the installer's log, https://tug.org/texlive/, describes how to get it. If you install a later release, the year in these paths may differ: use the directory that the last lines of your own installer output name. Any TeX installation that passes the checks below has the tools and every package the report loads.

```bash
which pdflatex latexmk
pdflatex --version | head -1
latexmk --version | grep Version
kpsewhich geometry.sty booktabs.sty graphicx.sty siunitx.sty amsmath.sty caption.sty subcaption.sty longtable.sty hyperref.sty microtype.sty enumitem.sty titlesec.sty tikz.sty | wc -l
```

Expected (the project's machine): two paths (there, under `/usr/local/texlive/2024/bin/x86_64-linux/`), then `pdfTeX 3.141592653-2.6-1.40.26 (TeX Live 2024)`, then a line ending `7 Apr. 2024. Version 4.85`, then `13`. The last number counts the packages found. If it is below 13, run the `kpsewhich` command without `| wc -l`: the package whose path is not printed is the missing one.

This guide only checks the tools. How to build the report is in [07_postprocessing_and_report.md](./07_postprocessing_and_report.md).

### A10. Check that the recipes arrived intact

Each folder in `recipes/` carries a `MANIFEST.sha256` file with a checksum for each recipe file (exactly which files is stated below). Checking it shows whether any listed file was changed or damaged since the recipes were published. This needs only `sha256sum`, not OpenFOAM:

```bash
cd $REPO/recipes/L11_wallResolved_CR && sha256sum -c MANIFEST.sha256
```

Observed output:

```
0.orig/U: OK
0.orig/k: OK
0.orig/nut: OK
0.orig/omega: OK
0.orig/p: OK
RECIPE_FILES.json: OK
constant/momentumTransport: OK
constant/physicalProperties: OK
system/blockMeshDict: OK
system/controlDict: OK
system/decomposeParDict: OK
system/fvSchemes: OK
system/fvSolution: OK
system/layerSettings: OK
system/meshLevels: OK
system/meshQualityDict: OK
system/planeSample: OK
system/snappyHexMeshDict: OK
system/surfaceFeaturesDict: OK
```

All four recipe folders at once, printing only failures plus one summary line per folder:

```bash
cd $REPO/recipes
for d in L11_wallResolved_CR L11_wallResolved_WT M6_wallModelled_CRUISE_EARLY M6_wallModelled_CRUISE_LATE; do
  (cd $d && sha256sum -c --quiet MANIFEST.sha256 && echo "$d: all files OK")
done
```

Observed output:

```
L11_wallResolved_CR: all files OK
L11_wallResolved_WT: all files OK
M6_wallModelled_CRUISE_EARLY: all files OK
M6_wallModelled_CRUISE_LATE: all files OK
```

What a failure looks like. This was produced on a copy of the recipe with one line added to `system/controlDict`:

```
system/controlDict: FAILED
sha256sum: WARNING: 1 computed checksum did NOT match
```

If you see that in `recipes/`, the recipe is not the published one. Get a clean copy of the repository before building anything from it.

What the check covers: the two `L11_` manifests list all 19 files in their folders. The two `M6_` manifests list 22 of the 23 files in their folders; `RECIPE_FILES.json` is not in them, so a change to that one file would not be detected. No manifest lists itself.

Run this check on the recipe folders in the repository only. A case copied from a recipe changes by design in two ways, and both make `system/controlDict` report `FAILED`:

1. You install the post-processing block, which adds an include line to `system/controlDict` and writes a new file, `system/argusPostPro`, that is not in the manifest at all. [04_send_to_hpc.md](./04_send_to_hpc.md), Step 5e, shows when and how.
2. You set the angle of attack ([06_trim_to_target_cl.md](./06_trim_to_target_cl.md), section 6.4.2), which also changes `0.orig/U`.

Before either change, a copied case reports every line `OK`. `hpc/mesh_hpc12.pbs` does not compare `system/controlDict` or `0.orig/U` with the manifest; it checks them for a consistent angle instead. See [05_submit_and_run.md](./05_submit_and_run.md), section 5, "What the mesh job checks before it meshes".

**A second check, for the two `L11_` recipes.** `scripts/verify_recipe.py` re-checks the manifest, compares 14 mesh settings in the recipe's dictionaries with the values pinned in `hpc/v5_recipe.json`, and checks that the flow direction and the lift and drag directions give one angle. It needs only the Python standard library. Run it from the repository root; from any other folder it stops with `FileNotFoundError: ... 'hpc/v5_recipe.json'`.

```bash
cd $REPO
python3 scripts/verify_recipe.py recipes/L11_wallResolved_CR
```

Observed output (speed `|U|` in m/s, angle in degrees):

```
  manifest: 19 files verified
  no suffixed variants: the folder is unambiguous
  pinned scalars compared: 14
  alpha triple agrees: |U| 34.0000 m/s, alpha 1.3000 deg, axis lag -0.000000 deg, orthogonal, unit
RECIPE OK
```

Three things to know about it:

1. The two cruise recipes (`M6_...`) end with `RECIPE REFUSED: 4 problems`. That is by design: the pinned layer settings are those of the `L11_` recipes. The manifest check above is the check for the cruise recipes.
2. `--self-test` does not work in this repository's layout: it prints `self-test SKIPPED: recipe folder absent` and exits with status 1.
3. `L11_wallResolved_WT` passes both checks here, but `hpc/mesh_hpc12.pbs` refuses even a fresh copy of it: its lift and drag directions are written to 6 significant figures, which gives an axis lag of -3.66e-06 degrees against the mesh job's tolerance of 1e-6 degrees (`dragDir implies alpha 2.102797, Uinf implies 2.102800`). It is not one of the three delivered conditions. The Condition CR and the two cruise recipes pass that check.

---

## Part B. The HPC12 cluster

> The HPC12 route is three files: `hpc/mesh_hpc12.pbs` (the mesh job), `hpc/stage_solve_case.sh` (turns a finished mesh case into a solve case at the login prompt; it needs only bash, coreutils, sed and grep, not OpenFOAM or Python) and `hpc/solve_hpc12.pbs` (the solve job). `hpc/README_HPC.md`, `hpc/submit.sh` and `hpc/push_to_cluster.sh` are from an earlier DelftBlue/SLURM setup and are not the HPC12 route.

### B1. What you need on HPC12

1. **An HPC12 account.** You log in as `<your-netid>@hpc12.tudelft.net`; setting up SSH access and keys for it is your own set-up. [04_send_to_hpc.md](./04_send_to_hpc.md) checks the login and covers moving files.
2. **Access to the `fpt-*` queues.** The shipped scripts submit to `fpt-medium` (solve) and `fpt-large` (mesh). These queues accept only members of the group `lr-fpt`. The `asm-*` queues require `lr-asm`, which the project's account did not have (notes inside `hpc/solve_hpc12.pbs`). If you are not in `lr-fpt`, ask the cluster administrators. You can read each queue's limits at the login prompt:

   ```bash
   qstat -Qf fpt-medium fpt-large fpt-small | grep -E "Queue:|max_user|resources_max|resources_default|max_running|max_queuable|acl_group"
   ```

   Observed output (HPC12, September 2026):

   ```
   Queue: fpt-medium
       max_user_queuable = 3
       resources_max.nodes = 4
       resources_max.walltime = 72:00:00
       resources_default.walltime = 72:00:00
       acl_group_enable = True
       acl_groups = lr-fpt,delftblue-groupadmins
   Queue: fpt-large
       max_user_queuable = 1
       resources_max.walltime = 72:00:00
       resources_default.walltime = 72:00:00
       acl_group_enable = True
       acl_groups = lr-fpt,delftblue-groupadmins
   Queue: fpt-small
       max_user_queuable = 10
       resources_max.nodes = 1
       resources_max.walltime = 72:00:00
       resources_default.walltime = 72:00:00
       acl_group_enable = True
       acl_groups = lr-fpt,delftblue-groupadmins
   ```

   In words: you can have at most 3 jobs (queued or running) in `fpt-medium` and 1 in `fpt-large`, and no job may ask for more than 72 hours of walltime. [05_submit_and_run.md](./05_submit_and_run.md), section 3, explains the queues and node types.
3. **OpenFOAM-org 12 built in `$HOME/OpenFOAM/OpenFOAM-12`.** Both job scripts look there by default. See B3 and B4.
4. **A copy of this repository in your cluster home.** [04_send_to_hpc.md](./04_send_to_hpc.md) makes it. Set `$REPO` to it in each login shell, as a full path (see "Two conventions"): the job scripts and `hpc/stage_solve_case.sh` are run from there, and `hpc/stage_solve_case.sh` refuses to work when the path of your copy contains a space or a comma.
5. **`python3` on HPC12.** You install nothing: the system `python3` is 3.6.8 on the login node (with numpy 1.12.1 and no jsonschema), and 3.6 on the compute nodes. Both job scripts run small Python checks that use only the standard library. The project's solve jobs ran them on HPC12 nodes: the convergence verdicts quoted in `data/run_cards/*.json` came from them.
6. **`ARGUS_GATE` in the submission line of every solve leg.** Before each solve leg, `hpc/solve_hpc12.pbs` runs the post-processing check, `python3 <checker> --before <case>`, if it finds the checker file. The check asks whether the case will write everything the post-processing needs ([07_postprocessing_and_report.md](./07_postprocessing_and_report.md)). If it fails, the leg stops with `REFUSED: post-processing gate refused` before it runs anything on the case. It leaves only an empty `log/leg<N>` folder, which moves the next leg's number on by one ([05_submit_and_run.md](./05_submit_and_run.md), section 10a, point 5). The job looks for the checker at the path in the variable `ARGUS_GATE`. Without `ARGUS_GATE`, it looks at `$HOME/argus_hpc12/assert_case_postpro.py`, which does not exist on a new account; it then prints `NOTE: <path> absent, post-processing gate NOT run` and runs the leg unchecked. So:
   1. Put `ARGUS_GATE=$REPO/scripts/assert_case_postpro.py` (the checker shipped in your copy of this repository) in the `-v` list of every solve leg you submit. You do not need to create anything in `$HOME/argus_hpc12/`. A real PBS job has received `ARGUS_GATE` through `-v` on HPC12 and ran the check (2026-09-24, `8 of 8 requirements satisfied`). The gate lines of the job script, run unchanged on staged cases, ran the check when given the checker's path and printed the `NOTE` line when not. Once a leg has started, its live output, `$HOME/<job number>.hpc12.hpc.OU`, shows which happened; the output file `<job name>.o<job number>` appears in the folder you ran `qsub` from only when the job ends ([05_submit_and_run.md](./05_submit_and_run.md), sections 7b and 11b).
   2. Submit every leg yourself. On HPC12 a job cannot submit its own next leg: its `qsub` fails with `could not connect to trqauthd`. See [05_submit_and_run.md](./05_submit_and_run.md), section 10a.
   3. A case passes the check only once its post-processing block, `system/argusPostPro`, is installed. The normal place to install it is your own machine, before you send the case to HPC12 ([04_send_to_hpc.md](./04_send_to_hpc.md), Step 5e). The installer, `scripts/install_postpro.py`, needs numpy, so you run it in your Python terminal (A5).

### B2. What the shipped job scripts load

`hpc/solve_hpc12.pbs` (lines 130 to 145) and `hpc/mesh_hpc12.pbs` (lines 36 to 50) set up their environment with these lines (condensed here: the real scripts also switch `set -u` off around `source`, see below). This is what the scripts do; you do not type it:

```bash
module load devtoolset/11     2>/dev/null || source /opt/rh/devtoolset-11/enable 2>/dev/null
module load mpi/openmpi-4.1.2 2>/dev/null
OF_BASHRC="${ARGUS_OF_BASHRC:-$HOME/OpenFOAM/OpenFOAM-12/etc/bashrc}"
source "$OF_BASHRC"
```

After that, they check the result rather than trusting it:

1. `hpc/solve_hpc12.pbs` stops with `REFUSED: OpenFOAM-12 not found at ...` if the environment file is missing, and with `REFUSED: foamRun not on PATH after sourcing ...` if loading it did not work. When it works, the job log shows `foamRun: <path>` and `WM_OPTIONS=<value>`.
2. `hpc/mesh_hpc12.pbs` stops with `REFUSED: <utility> not on PATH after sourcing ...` unless `snappyHexMesh`, `blockMesh`, `decomposePar`, `reconstructPar`, `surfaceFeatures` and `checkMesh` are all found. When they are, the job log shows `OpenFOAM: <path>  WM_PROJECT_VERSION=12`.

Loading OpenFOAM-12 on HPC12 prints `dirname: missing operand` and `Try 'dirname --help' for more information.`, so every mesh and solve job log contains that pair just before its `foamRun:` or `OpenFOAM:` line (HPC12 job logs). It is harmless (B4).

If your OpenFOAM-12 is somewhere else, set `ARGUS_OF_BASHRC` to its `etc/bashrc` when you submit. [05_submit_and_run.md](./05_submit_and_run.md) shows the submission line. Because you submit every solve leg by hand on HPC12 ([05_submit_and_run.md](./05_submit_and_run.md), section 10a), put `ARGUS_OF_BASHRC` in the `-v` list of every leg.

**If you write your own job script:** OpenFOAM's `etc/bashrc` reads a variable that is not set under bash, so a script running with `set -u` dies when it sources the file:

```
/opt/openfoam12/etc/bashrc: line 46: ZSH_NAME: unbound variable
```

That message was reproduced on the local package. The project's build job hit the same failure on HPC12. The shipped scripts turn `set -u` off just around the `source` line (`set +u`, then `source`, then `set -u`).

### B3. Building OpenFOAM-12 on HPC12

The project built OpenFOAM-org 12 from source in its home directory on HPC12. **The build script is not in this repository.** What the project's records establish about the result:

1. The source is in `$HOME/OpenFOAM/OpenFOAM-12`, with the third-party packages in `$HOME/OpenFOAM/ThirdParty-12` beside it.
2. The programs are in `$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin` (148 of them). `WM_OPTIONS` is therefore `linux64GccDPInt32Opt`.
3. The third-party build provides the scotch decomposition library. The project's mesh jobs used it: `data/run_cards/SWB_r2_trim.json` records, under `mesh.mesh_generation`, the method `scotch` at 192 ranks, read from the mesh case's `system/decomposeParDict`. The solves that followed used the method `hierarchical` at 128 ranks (fields `decomposition.method` and `decomposition.n_subdomains` in the same run card).

The outline below describes what the project's build job did. **NOT VERIFIED:** the copy of the build job available when this guide was written is known to differ from the copy that ran on the cluster, and this outline was not re-run. Treat it as a starting point, and run it as a batch job, not on the shared login node, because the build keeps 20 cores busy for hours.

```bash
#!/bin/bash
# NOT VERIFIED: outline of the project's OpenFOAM-12 build job on HPC12
#PBS -N argus-of12-build
#PBS -q fpt-small
#PBS -l nodes=1:ppn=20,walltime=12:00:00
#PBS -j oe
module load devtoolset/11 || source /opt/rh/devtoolset-11/enable
module load mpi/openmpi-4.1.2
mkdir -p $HOME/OpenFOAM && cd $HOME/OpenFOAM
wget -O - http://dl.openfoam.org/source/12 | tar xz
wget -O - http://dl.openfoam.org/third-party/12 | tar xz
mv OpenFOAM-12-version-12 OpenFOAM-12
mv ThirdParty-12-version-12 ThirdParty-12
export LIBRARY_PATH=$(mpicc -show | tr ' ' '\n' | sed -n 's/^-L//p' | head -1):$LIBRARY_PATH
source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc
cd $HOME/OpenFOAM/OpenFOAM-12
./Allwmake -j 20 > $HOME/OpenFOAM/Allwmake.log 2>&1
```

To run it, save the outline in your cluster home as a file, for example `build_of12.pbs`, and submit it from there (NOT VERIFIED, like the outline):

```bash
cd $HOME
qsub build_of12.pbs
```

`qsub` prints the job identifier (a number followed by the server name). PBS writes the job's own messages to `argus-of12-build.o<job id>` in the folder you ran `qsub` from; [05_submit_and_run.md](./05_submit_and_run.md), section 11, shows how to follow a job. The build's own messages go to `$HOME/OpenFOAM/Allwmake.log` (the project's job wrote a log with a date in its name in the same folder). This is "the build log" in B4.

Lessons from the project's build:

1. **The compiler.** The cluster's default compiler is too old for OpenFOAM-12 (gcc 4.8.5, according to the build notes). `module load devtoolset/11` provides a newer one: gcc 11.2.1 at the HPC12 login prompt. B4 shows the check.
2. **`ld: cannot find -lmpi`.** This happened while building the third-party packages. The MPI module sets `LD_LIBRARY_PATH` (used when a program runs) but not `LIBRARY_PATH` (used when a program is linked). The `export LIBRARY_PATH=...` line above fixes it.
3. **A non-zero exit code is not the verdict.** `./Allwmake` finished with exit code 2 while every program had been built. The only part that failed was ParaView's OpenFOAM reader, which needs VTK headers the cluster does not have. Judge the build by the check in B4, not by the exit code.
4. **No source patch.** A patch to `src/OpenFOAM/db/dictionary/entry/entryIO.C` circulates for building OpenFOAM-11 on HPC12. Do not apply it to OpenFOAM-12: the code it changes is not in OpenFOAM-12's version of that file.
5. **`set -u`.** See the end of B2.

### B4. Check the cluster installation

Log in to HPC12 (`ssh <your-netid>@hpc12.tudelft.net`; [04_send_to_hpc.md](./04_send_to_hpc.md), Step 2, checks the login). Before anything else, OpenFOAM is not loaded: `command -v foamDictionary` prints nothing. Type these three lines at the prompt, one at a time. They are the lines the shipped job scripts run (B2), without the fallback after `||` and without hiding error messages:

```bash
module load devtoolset/11
module load mpi/openmpi-4.1.2
source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc
```

The two `module` lines print nothing. The `source` line prints two lines:

```
dirname: missing operand
Try 'dirname --help' for more information.
```

**These two lines are harmless.** OpenFOAM-12 is loaded after them, and the same pair appears in every mesh and solve job log of the project's HPC12 runs (B2).

Type the three lines directly, not through a pipe. A command such as `source ... | grep -v dirname` runs in a separate shell, whose settings are thrown away when it ends, so OpenFOAM is not loaded afterwards. If the `module` lines were piped as well, the `source` line also prints `.../OpenFOAM-12/etc/config.sh/mpi: line 46: mpicc: command not found`.

Now check the result, in the same shell:

```bash
echo $WM_PROJECT_VERSION
gcc --version | head -1
command -v foamRun potentialFoam snappyHexMesh blockMesh decomposePar reconstructPar surfaceFeatures checkMesh foamDictionary
echo $WM_OPTIONS
foamDictionary -help | grep -E "^(Using|Build):"
```

Observed output (HPC12, September 2026; the paths show your real home folder where this shows `$HOME`):

```
12
gcc (GCC) 11.2.1 20220127 (Red Hat 11.2.1-9)
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/foamRun
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/potentialFoam
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/snappyHexMesh
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/blockMesh
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/decomposePar
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/reconstructPar
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/surfaceFeatures
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/checkMesh
$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/foamDictionary
linux64GccDPInt32Opt
Using: OpenFOAM-12 (see https://openfoam.org)
Build: 12
```

How to read it:

1. `12`: OpenFOAM-12 is loaded. An empty line means the `source` line did not take effect in this shell.
2. The `gcc` line must show version 9 or newer; the project's build job refused anything older. If it shows 4.8.5, the cluster's default compiler, the `devtoolset/11` module did not load. The job scripts then try `source /opt/rh/devtoolset-11/enable` instead, but hide its error messages, so this is the place to see them.
3. Nine paths: the eight OpenFOAM programs the two job scripts call, and `foamDictionary`, which the guides use at the prompt. Fewer than nine means the build is incomplete: look in `$HOME/OpenFOAM/Allwmake.log` for the first error, and see lesson 3 in B3 before concluding the build failed. A missing `potentialFoam` does not stop a solve job: the job prints `potentialFoam returned non-zero (continuing; it is an initialiser)` and solves without that starting step, so check for it here.
4. `linux64GccDPInt32Opt`.
5. `Using: OpenFOAM-12`. The `Build:` line here is just `12`, where your local package prints `12-86e126a7bc4d` (A3). Both are OpenFOAM-12.

Next, check the system Python (B1, item 5). You install nothing, and the output is the same with or without the three lines above:

```bash
command -v python3
python3 --version
python3 -c "import numpy; print(numpy.__version__)"
```

Observed output (HPC12, September 2026):

```
/usr/bin/python3
Python 3.6.8
1.12.1
```

The post-processing check, `scripts/assert_case_postpro.py`, needs only the standard library, and this `python3` runs it (A8). The numpy line matters only for the fallback of installing the post-processing block on HPC12 (A8). If it prints `ModuleNotFoundError: No module named 'numpy'`, install the block on your own machine instead ([04_send_to_hpc.md](./04_send_to_hpc.md), Step 5e).

You do not need OpenFOAM loaded at the login prompt to submit a job, because the job scripts load it themselves. You need it there for `foamDictionary`, which later guides use on HPC12 (remember the relative-path rule in A6).

Last, check that the build reads dictionaries that use `#include` and `#includeFunc`. The recipes use both (`#include` in `system/snappyHexMeshDict`, `#includeFunc` in the cruise recipes' `system/controlDict`). The older OpenFOAM-11, built on this cluster, crashed on `#includeFunc` (even in `foamDictionary`) until its source was patched. Lesson 4 in B3 explains why OpenFOAM-12 needs no patch; this check shows that your build handles both directives:

```bash
mkdir -p $HOME/of12-dict-test && cd $HOME/of12-dict-test
printf 'FoamFile { version 2.0; format ascii; class dictionary; object inc; }\nalpha 1.234;\n' > inc
printf 'FoamFile { version 2.0; format ascii; class dictionary; object t1; }\n#include "inc"\nbeta 5.678;\n' > t1
printf 'FoamFile { version 2.0; format ascii; class dictionary; object t2; }\nfunctions\n{\n    #includeFunc residuals\n}\n' > t2
foamDictionary -expand -entry alpha -value t1
foamDictionary -expand -entry functions/residuals/type -value t2
```

Expected:

```
1.234
residuals
```

The first line shows that `#include "inc"` was read, the second that `#includeFunc residuals` was expanded. A broken directive gives `FOAM FATAL IO ERROR`, or a crash, instead. The `printf` and `foamDictionary` lines were run on HPC12, after the three `module` and `source` lines (in a temporary folder rather than `$HOME/of12-dict-test`), and printed exactly these two lines. They were also run verbatim on the local OpenFOAM-12 package, with two failure checks: a missing include file and an unknown function name each stopped with `FOAM FATAL IO ERROR`. Keep `-expand`: without it, `foamDictionary` does not expand these directives, so both lines stop with `FOAM FATAL IO ERROR` (`Cannot find entry ...`) even on a working build. You can delete `$HOME/of12-dict-test` afterwards.

---

## Verify your setup: checklist

Every command below was run for this guide: items 1 to 13 on the project's local machine, item 14 at an HPC12 login shell. The expected output is the pass condition. Where a value in brackets is marked "the project's", it is what the project's machine printed, and yours may differ (A3). Run the local OpenFOAM items in your OpenFOAM terminal, where OpenFOAM-12 is loaded (A2) and `mpirun` is Open MPI (A5). Run the Python items in your Python terminal (A5), with your Python 3.11 (A8), from `$REPO`.

| # | Command | Expected output |
|---|---|---|
| 1 | `dpkg -l openfoam12 \| tail -1` | a line starting `ii  openfoam12` (the project's version: `20250206`) |
| 2 | `foamVersion` | `OpenFOAM-12` |
| 3 | `which foamRun` | `/opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/foamRun` |
| 4 | `foamRun -help \| grep -E "^(Using\|Build):"` | `Using: OpenFOAM-12 (see https://openfoam.org)`, then a `Build:` line (the project's: `Build: 12-86e126a7bc4d`) |
| 5 | `mpirun --version \| head -1` | a line starting `mpirun (Open MPI)`, not `HYDRA build details:` (the project's: `mpirun (Open MPI) 4.0.3`) |
| 6 | the parallel smoke test in A6, then `grep "^nProcs" log.foamRun` | `nProcs : 2` |
| 7 | `pvbatch --version` | `paraview version 5.10.1` |
| 8 | `python3 -c "import numpy, scipy, matplotlib, yaml, jsonschema, duckdb; print('python packages ok')"` | `python packages ok` |
| 9 | `python3 -c "import jsonschema; jsonschema.Draft202012Validator; print('jsonschema ok')"` | `jsonschema ok` |
| 10 | `python3 scripts/validate_run_card.py --self-test` | `self-test PASS (example accepted; 4 corrupted variants rejected)` |
| 11 | the `kpsewhich ... \| wc -l` line in A9 | `13` |
| 12 | the four-folder loop in A10 | four lines ending `: all files OK` |
| 13 | `python3 scripts/verify_recipe.py recipes/L11_wallResolved_CR` | last line `RECIPE OK` |
| 14 | on HPC12, the three lines at the top of B4, then the three check blocks in B4 | `12`; a `gcc` version of 9 or newer (11.2.1 on HPC12); nine paths under `$HOME/OpenFOAM/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/`; `linux64GccDPInt32Opt`; `Using: OpenFOAM-12 (see https://openfoam.org)` and `Build: 12`; then `/usr/bin/python3`, `Python 3.6.8` and `1.12.1`; then `1.234` and `residuals` |

When every line matches, go on to [02_geometry_from_openvsp.md](./02_geometry_from_openvsp.md).
