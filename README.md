# ARGUS morphing wing: RANS validation of the VLM design study

OpenFOAM-org 12 RANS workflows for validating the OpenVSP/VSPAERO outer-wing morphing design
study of the ARGUS project at TU Delft: case recipes, mesh and solve jobs for the HPC12 cluster,
trim to target lift, post-processing, run cards, and the validation report. It is the viscous
counterpart of the
[ARGUS Aerodynamic Morphing-Wing Workflow](https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow)
(Liming Zheng, TU Delft), whose geometries and vortex-lattice results it checks, and the
maintained starting point for further CFD work on the ARGUS wing.

This repository was produced within the ARGUS research project by **Tyler Buchanan** (TU Delft).
Original software and documentation are available under the MIT License
([LICENSE](LICENSE)); third-party software, NASA-derived material, design-study data and
branding retain their own terms. See [NOTICE.md](NOTICE.md) for provenance and license scope,
and [CITATION.cff](CITATION.cff) to cite it.

This repository provides computational tools, methods, and workflows developed by Delft
University of Technology (TU Delft) within Work Package 1 (Modelling and Simulation),
specifically Task 1.3 (CFD Simulation-Based Design of Aerodynamic Forms for Demonstrator), of the
ARGUS project (Adaptive and Responsive Morphing Surfaces for Green Aviation Using Shape Memory
Alloys). Led by the University of Twente, the project brings together TU Delft, the Fraunhofer
Innovation Platform for Advanced Manufacturing at the University of Twente (FIP-AM@UT),
Nederlands Lucht- en Ruimtevaartlaboratorium (NLR), and the Materials innovation institute
(M2i). ARGUS aims to advance wing morphing for more sustainable aviation by combining advanced
materials, computational modelling, and intelligent control to improve aerodynamic efficiency,
reduce drag, and lower aviation's environmental impact. The repository supports the public
dissemination, reuse, and further development of these research outputs. This research (project
number LITO&O24-03725850) is part of the Netherlands National Growth Fund program Luchtvaart in
Transitie, which is partly financed by the Dutch Research Council NWO
([www.nwo.nl](https://www.nwo.nl)). It is also part of the Partnership Program of the Materials
innovation institute M2i ([www.m2i.nl](https://www.m2i.nl)) (with project number N24014).

This repository is the curated project handover, not a dump of every solver run. It contains
the case recipes, the Python and shell tooling, the run cards, audited summary data, the report
and its sources, the group presentation, and seven guides that take a new researcher from
installation to a trimmed, post-processed case. Meshes, solution fields and complete case
folders are in the project archive on the ARGUS project drive (`argus_hpc12_backup/`), and the
`.vsp3` geometries come from the design study (docs/02, Step 1).

ARGUS asks how the outer wing of a transport aircraft should adapt during cruise. A low-order
design study (a design-space optimisation, "DSO", run with the VSPAERO vortex-lattice method,
"VLM") searched that design space and selected morphing candidates on induced drag alone. This
repository is the viscous check of that selection.

Six wings were meshed with snappyHexMesh and solved with steady RANS: one rigid baseline and five
morphing candidates of the NASA EET AR-12 cruise wing (the wing tested at low speed in NASA
TP-1580), each as a half model. They were solved at three operating points, and every geometry
was trimmed to each operating point's own target lift coefficient, so geometries are compared at
the same C_L, not at the same angle of attack. The trim tolerance is 1 count of C_L
(|dC_L| < 1e-4), and all 18 trims meet it; the widest is `LCW_trim` at 0.895 counts.
`data/rans_forces.json` records every residual under `_trim_quality`, with the drag bias it
implies (0.036 counts at most).

The report is `report/argus_rans_validation_report.pdf` (36 pages, 9 October 2026). Start there if you only want the results, but
read [Which mesh each result is on](#which-mesh-each-result-is-on-read-this-before-using-any-number)
before quoting any number from it.

## Release of 9 October 2026: what changed, and why

This release replaces every early- and late-cruise result. If you read the group deck of
23 September or the report of 15 September, the cruise conclusions there are withdrawn.

1. **What was wrong.** The first cruise meshes of F, H and W were built from an unwelded copy of
   their surfaces (Known gaps, item 8). Under the 40,000,000-cell cap the leaked wing interior
   used up the cell budget and starved those meshes near the wing and in the wake. They read
   +50.1 to +56.5 counts of C_D against B.
2. **What was done.** All twelve cruise trims (six geometries, two cruise points) were meshed
   again from the registered surfaces with the shipped cruise recipes and the 130,000,000-cell
   cap (73.64 to 73.75 million cells), and re-trimmed to |dC_L| < 1e-4. All six geometries were
   rebuilt, not only F, H and W, so every within-condition difference is taken inside one mesh
   generation.
3. **How it was checked.** B, C and M had sound meshes before, so a correct rebuild must
   reproduce them, and it does: they move by 0.10 to 0.34 counts. F, H and W move by 51.0 to
   55.3 counts. Their earlier penalties were the mesh, not the shapes.
4. **What it changes** (C_D counts, DSO basis, M 0.78; early cruise at C_L 0.529297087450, late
   cruise at C_L 0.544114800210; each difference against that condition's own B trim):

   | Geometry | Early cruise | Late cruise |
   |---|---|---|
   | B (baseline) | 177.839 | 187.959 |
   | W minus B | -0.769 | -0.781 |
   | F minus B | +1.796 | +1.348 |
   | H minus B | +2.244 | +1.927 |
   | C minus B | +12.871 | +12.214 |
   | M minus B | +14.554 | +13.728 |

   W is the only candidate that saves drag at cruise, and the cheapest geometry at all three
   operating points. The cruise order W, F, H, C, M differs from the Condition CR order only in
   the C and H pair, 0.18 counts apart at Condition CR. C and M carry a strong aft shock at
   eta 0.80 and separate behind it; F and H a moderate shock, attached; W none.
5. **What is new in the report.** Late cruise now has its spanwise loading, section cuts and
   surface maps, and its root bending; a table of every trim leg; and the cruise hinge moments
   and root bending on the rebuilt meshes. Condition CR results are unchanged.
6. **Where the record is.** `data/rans_forces.json` keeps each replaced cruise value under
   `superseded_entry` and `superseded_trim_entry`; `report/src/sections/gen/tab_generation_null.tex`
   is the earlier-against-rebuilt table behind item 3. The report itself states only the final
   numbers.

## Supported toolchain

**OpenFOAM-org 12** (the OpenFOAM Foundation release; its solver application is `foamRun`) on
your own Linux machine, and on the **HPC12** cluster (login `<your-netid>@hpc12.tudelft.net`)
with **PBS** (`qsub`, `qstat`). The guides in `docs/` cover nothing else. SSH access to HPC12,
and your SSH keys, are your own set-up.

Install it first, on both machines: [docs/01_setup.md](docs/01_setup.md). Then load the
version 12 environment. On your own machine (packaged install under `/opt/openfoam12`):

```bash
source /opt/openfoam12/etc/bashrc
```

On HPC12 you need OpenFOAM-org 12 built from source in `$HOME/OpenFOAM/OpenFOAM-12`. The job
that builds it is **not in this repository**; [docs/01_setup.md](docs/01_setup.md) (section B3)
gives an outline of the project's build job, marked NOT VERIFIED there because it was not
re-run. Once OpenFOAM-12 exists there, load the same compiler and MPI
modules and the same OpenFOAM location that `hpc/mesh_hpc12.pbs` and `hpc/solve_hpc12.pbs` use.
Type the three lines at the HPC12 login prompt, one after the other, in the same shell. Do not
run them through a pipe: a pipe runs them in a separate shell, and what they set is lost.

```bash
module load devtoolset/11
module load mpi/openmpi-4.1.2
source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc
```

On HPC12 the `source` line prints two lines, `dirname: missing operand` and
`Try 'dirname --help' for more information.`. They are harmless; the same pair appears in
every job's output file.

Then check it, on either machine:

```bash
echo "WM_PROJECT_VERSION=$WM_PROJECT_VERSION"
command -v foamRun
```

**What you should see:** `WM_PROJECT_VERSION=12`, then a path ending in `/bin/foamRun` (under
`/opt/openfoam12/` on your own machine, under `$HOME/OpenFOAM/OpenFOAM-12/` on HPC12). On
HPC12 this was checked on the login node, after the three lines above in one shell; before
them, the shell does not find `foamDictionary`. You do not need OpenFOAM loaded to submit a
job, because the job scripts load it themselves. If you built OpenFOAM-12 somewhere else on
HPC12, add `ARGUS_OF_BASHRC=<path to its etc/bashrc>` to the `-v` list of every job.

**If it does not:** if the first line shows another version or is empty, or `command -v foamRun`
prints nothing, your shell has loaded a different OpenFOAM version (or none). Source the
version 12 file again in that shell (it replaces the other version) and repeat the check. If
the file you source does not exist, OpenFOAM-org 12 is not installed at that location: go back
to [docs/01_setup.md](docs/01_setup.md).

Several of the Python scripts, the post-processing installer among them, need numpy;
[docs/01_setup.md](docs/01_setup.md) says how to set up Python.

Every guide refers to your clone of this repository as `$REPO`:

```bash
export REPO=$HOME/argus-production   # or wherever you cloned it
```

On HPC12, `$REPO` is the copy of this repository that
[docs/04_send_to_hpc.md](docs/04_send_to_hpc.md), Step 4, puts there, not a clone. Two rules
apply to it, because every solve leg's `qsub -v` list carries paths under it:

1. its full path must contain no space and no comma;
2. `$REPO` must be an absolute path, set in the shell you submit from. The job looks for the
   files those paths name only after it has moved into the case folder. A path starting with
   `$HOME` is absolute.

## What is in this repository

| Path | What it holds |
|---|---|
| `README.md` | This file. |
| `LICENSE`, `NOTICE.md`, `CITATION.cff` | The MIT License for the original software and documentation, the provenance and scope of that license, and how to cite this repository. |
| `docs/` | The seven guides, in reading order (see the quickstart below). |
| `recipes/` | The OpenFOAM-org 12 case setups, one folder per condition (`0.orig/`, `constant/`, `system/`, plus `MANIFEST.sha256` and `RECIPE_FILES.json`). The two cruise folders also hold `CONDITION.json` (the cruise flight state and what the recipe derives from it, including the Reynolds-number matching) and `variants/fvSchemes.nolimiter` (the same schemes without the gradient limiter; `recipes/README.md` says when to use it), both covered by their manifests. No geometry and no mesh. `recipes/README.md` explains them, including the section "Decomposition and run length, aligned with the r2 cases". `L11_wallResolved_WT` is a tunnel-matched condition that is not one of the three delivered conditions and was not given the wake refinement. |
| `hpc/` | The HPC12 (PBS) route is three files: `mesh_hpc12.pbs` (the mesh job), `stage_solve_case.sh` (turns a finished mesh case into a solve case) and `solve_hpc12.pbs` (the solve job, one leg per submission). `templates/argusPostPro` is the post-processing block that `scripts/install_postpro.py` writes into a case. `v5_recipe.json` holds the mesh settings that `scripts/verify_recipe.py` checks the two `L11` recipes against. `README_HPC.md`, `submit.sh` and `push_to_cluster.sh` are from an earlier DelftBlue/SLURM setup, depend on files not shipped here, and are not the HPC12 route. |
| `scripts/` | 71 Python scripts (geometry conversion, input checks, the post-processing installer and checker, post-processing, and run-card derivation and validation), plus `run_card.schema.json`, the schema the run cards are validated against. |
| `registry/` | `candidates.yaml`: binds each candidate id to its geometry files and their sha256 checksums. |
| `data/` | Small result files (mostly JSON): `rans_forces.json` (forces for every converged case), `span_efficiency.json` (Trefftz-plane span efficiency at Condition CR, built from `span_eff_r2trim_raw.json`), `vlm.json` (the extracted VLM results), `dso_vlm_claims.json` (the design study's own VLM statements, each with its operating point and frame), `lineage_measured.json` (the measured difference between B's aerofoil lineage and the corrected one, behind the report's lineage figure), case provenance, `spanwise/` (M's dense spanwise loading at Condition CR), `phaseA_run_card_exemption.json` (a working record about an earlier campaign; none of its cases is among the delivered results and no guide uses it), and `run_cards/` with one run card per trim, all eighteen. A run card is the JSON record of one run: geometry checksum, mesh statistics, solver settings, cost and results. |
| `report/` | `argus_rans_validation_report.pdf`; its LaTeX source in `report/src/`; figures, the scripts that draw them and `pipeline_diagram.tex` in `report/fig/`; table and data builders in `report/scripts/`. |
| `presentations/` | The group slide deck of 9 October 2026 with the delivered results, `argus_group_2026-10-09`, as `.pptx` and `.pdf` (30 slides). It replaces the deck of 23 September. |
| `template_case/airfoil2d/` | The OpenFOAM-org 7 numerics that `scripts/build_2d_case.py` copies for the 2D section studies. Not part of the guides' OpenFOAM-12 route. |

All eighteen run cards pass `scripts/validate_run_card.py`
([docs/07_postprocessing_and_report.md](docs/07_postprocessing_and_report.md), section 7.5). It
needs Python with jsonschema 4 or later. The HPC12 login node's `python3` has no jsonschema, so
validate on your own machine.

Two terms used throughout:

1. **Span efficiency e** comes from the induced drag measured on planes of constant x behind
   the wing (Trefftz planes). For a planar wake it cannot exceed 1.
2. **y+** is the distance of the first cell centre from the wall, in viscous units. y+ about 1
   puts the first cell inside the viscous sublayer, so the boundary layer is resolved rather than
   modelled.

## The three operating points

| Condition | Case prefix | Recipe folder | Mach | U (m/s) | Target C_L | Turbulence model | Wall treatment |
|---|---|---|---|---|---|---|---|
| CR (low speed) | `SW` | `L11_wallResolved_CR` | 0.10 | 34.0 | 0.428277635108 | kOmegaSST | resolved, y+ about 1 |
| Early cruise | `CMP` | `M6_wallModelled_CRUISE_EARLY` | 0.78 | 233.934028 | 0.529297087450 | Spalart-Allmaras | modelled (wall functions) |
| Late cruise | `LC` | `M6_wallModelled_CRUISE_LATE` | 0.78 | 230.501789 | 0.544114800210 | Spalart-Allmaras | modelled (wall functions) |

1. **"Condition CR" is the M 0.10 low-speed point** at which the design study selected its
   candidates. CR does not stand for cruise. "Early cruise" and "late cruise" are the two M 0.78
   states.
2. U is the free-stream speed each recipe's `forceCoeffs` carries (`magUInf`).
3. The cruise states are defined at aircraft scale but meshed at model scale. The recipes lower
   the viscosity so that the Reynolds number on the model chord `lRef` 0.393957 m equals the
   aircraft Reynolds number, which `CONDITION.json` computes on an aircraft chord of 3.196 m.
   Those two chords are not the same fraction of the wing: scaled down by the ratio of the spans
   (b_ref 3.6576 m against the aircraft's 38.36 m), 3.196 m corresponds to about 0.305 m, not to
   0.394 m. `recipes/M6_wallModelled_CRUISE_EARLY/CONDITION.json` and its `_LATE` twin record
   this Reynolds-matching choice, and its alternatives, under `scale_resolution`, as not formally
   approved; `recipes/README.md` lists it under "What is NOT settled".
4. Reynolds number as run, on the model reference chord `lRef` 0.393957 m: CR 917,434 (the same
   flow is 0.79e6 on the TP-1580 chord 0.3404 m); early cruise 21,176,454; late cruise
   16,107,442.
5. Condition CR is solved incompressible (`foamRun` solver `incompressibleFluid`); both cruise
   states are solved compressible (`fluid`).
6. The Condition CR wall conditions in `recipes/L11_wallResolved_CR/0.orig/` have
   `WallFunction` in their type names (`nutLowReWallFunction`, `kLowReWallFunction`,
   `omegaWallFunction`). They are not wall functions in the modelled sense:
   `nutLowReWallFunction` sets the turbulent viscosity to zero at the wall, which is the
   resolved-wall condition. The cruise recipes use `nutUSpaldingWallFunction`, which is a wall
   function.
7. The Condition CR target C_L is defined on the DSO basis. The RANS cases apply all three
   targets on the DSO basis; for the cruise targets that is an open question (see item 7 of the
   next section).

**The three operating points are never differenced against each other.** Each carries its own
baseline trim. A candidate's drag change is always taken against baseline B at the same
condition, never against a number from another condition. This is the easiest mistake to make
with these results.

## Reference quantities and units

Reference quantities are **dual**. State which basis a coefficient is on every time you quote it.

1. **DSO basis (primary)**: S_ref 1.24092 m², C_ref 0.39396 m. Every RANS coefficient in
   `data/` and the report is on this basis.
2. **TP-1580 basis (reported alongside)**: S_ref 1.1148 m², c_ref 0.3404 m.
3. b_ref 3.6576 m is common to both.
4. **Half model**: every case solves half a wing. `forceCoeffs` uses A_ref 0.620462 m² (half of
   the DSO S_ref) and `lRef` 0.393957 m, so a half-model coefficient equals the full-wing
   coefficient on the DSO basis.
5. The two reference areas differ by a factor of 1.1131, i.e. 11.31 %. The report puts that at
   about 28 times a morphing effect of order 0.4 %. Never compare coefficients across the two
   bases without converting, at any condition.
6. The VLM results in `data/vlm.json` keep the reference quantities read from each VSPAERO input
   deck: the Condition CR decks are on the DSO basis, the early- and late-cruise decks on the
   TP-1580 basis. Each run's `frame_conversion` block gives the factor to the other basis.
7. **Open frame question for the cruise states.** The RANS cruise cases trim to the cruise target
   C_L on the DSO basis (`recipes/M6_*/CONDITION.json`), while `data/vlm.json` states that the
   design study's cruise targets are on the TP-1580 basis. `data/rans_forces.json` records this
   disagreement as unresolved (`_open_cross_file_frame_conflict`). Until it is settled, do not put
   a cruise RANS coefficient beside a cruise VLM coefficient.
8. Drag is quoted in counts: 1 count = 0.0001 in C_D.

**Units trap.** OpenVSP files (the `.vsp3` design files and any STL exported from them) are in
**feet**. Scale by 0.3048 to metres before meshing, and strip non-wing solids (for example slat
tracks) by solid name. The meshing surfaces that
[docs/02_geometry_from_openvsp.md](docs/02_geometry_from_openvsp.md) builds are in metres; a raw
OpenVSP export is not.

## Which mesh each result is on (read this before using any number)

There are three mesh generations:

1. **published**: the earlier mesh. Its cell cap (`maxGlobalCells` 40,000,000) was reached on
  every mesh. On W, H and F the cap was used up early, by the refinement of a wing interior that
  leaked in through an unwelded root seam (Known gaps, item 8), before snappyHexMesh applied its
  other refinement regions (the wake box included), so their wake core carried 62.5 mm cells
  where the B mesh already had 7.8125 mm.
2. **wake-refined**: the cap raised to 130,000,000 and two refinement boxes added over the aft
  part of the wing and the near wake (x from 1.80 m to 3.40 m), giving 7.8125 mm cells in the
  wake core. `recipes/README.md` records the three changes and why.
3. **welded**: the wake-refined recipe built from the registered (welded) surfaces, used for every
  early- and late-cruise trim. Same recipe and cap as generation 2; only the surface differs for
  F, H and W.

The rules:

1. **Condition CR**: the trimmed result of each of the six geometries (case names ending
   `_trim`; run cards `data/run_cards/SW?_r2_trim.json`, where `r2` means wake-refined) is on
   the wake-refined mesh. These **supersede** the published CR results.
2. At Condition CR, the angle-of-attack sweep cases that bracket the trim angle (named like
   `SWB_a1p75`) were run on the published mesh and were not re-run. The report's lift-curve
   slopes (Table 17's dC_L/dalpha and implied-dC_L columns, the "Slope" sentence of section 5.3,
   and the slope in `report/fig/cl_alpha_condition_CR.pdf`) are fitted over these published
   sweep cases **plus** the wake-refined trim, so they mix two generations (rule 5).
3. **Early and late cruise**: every trim, and every off-trim leg that reached it, is on the
   welded mesh, all six geometries. These **supersede** the published cruise results. The
   published cruise alpha-bracket cases (named like `CMPB_a1p20`) remain in `data/` labelled
   `published` and are used by no cruise result, fit or figure.
4. `recipes/` ships the **wake-refined recipe for all three conditions**. The welded cruise meshes
   were built from it: in the B, F and W mesh cases at both cruise points, 13 of the 14 recipe
   files are byte-identical to `recipes/M6_wallModelled_CRUISE_EARLY` and `_LATE`, and
   `controlDict` differs only by the staging steps (`purgeWrite` 0, the post-processing include,
   and the case's own force axes). At Condition CR the recipe is set up as the delivered cases were: the mesh stage
   splits the mesh into 192 pieces with method `scotch`, as for all six wake-refined meshes; the
   first solve leg ends at `endTime` 4000 and saves every 500 iterations (`writeInterval`); and
   `hpc/stage_solve_case.sh` gives the solve case the delivered solves' 128-rank split
   (`hierarchical`, 8 x 2 x 8). No full-size case has been rebuilt from the recipe and compared
   with the delivered numbers, so do not use the delivered CR numbers as a regression test either.
5. **Never difference, fit or rank across generations**, and never put numbers from two of them
   in one table without showing each one's generation.
6. The published Condition CR span efficiencies **W 1.2861, H 1.2974, F 1.3028** are an artefact
   of the unresolved wake on the starved mesh (e above 1 is physically impossible for a planar
   wake). **Do not quote them.**
7. **The surface figures of `report/argus_rans_validation_report.pdf` still show M on its published mesh at
   Condition CR.** M's panels come from its published-mesh case in the pressure and
   skin-friction section figures (`report/fig/fig_cp_sections_CR.pdf`,
   `fig_cf_sections_upper.pdf`, `fig_cf_sections_lower.pdf`) and in the surface delta maps
   (`delta_map_CR_cp.png`, `delta_map_CR_cf_s.png`). B, C, F, H and W are wake-refined there.
   M's forces and span efficiency in `data/`, its spanwise loading in `data/spanwise/`, and its
   row of the report's trim-angle table (Table 17: 1.257423 deg, the value
   `report/scripts/make_rank_table.py` prints) are wake-refined.
8. **The report's hinge-moment and root-bending results are wake-refined at Condition CR and
   welded at early and late cruise.** This covers the tables and figures of the report
   subsections "Hinge moment about the morph line" and "Root bending, against the low-order
   method and against its screen" (`report/fig/fig_hinge_integrated.pdf`,
   `fig_hinge_sections.pdf`, `fig_rootbending.pdf`, and the values behind them in
   `report/fig/hinge_moment_values.csv` and `rootbending_values.csv`). At Condition CR all six
   geometries were integrated at iteration 6000 of the wake-refined trims of rule 1, at the
   wake-refined trim angles (B 1.747609, C 1.268311, F 1.475456, H 1.552310, M 1.257423,
   W 1.520545 deg); at cruise at each welded trim's final iteration. `report/fig/derived_provenance.json`
   records which source served each of the 18 rows: 6 wake-refined, 12 welded. So a table that
   sets Condition CR loads beside cruise loads crosses a mesh generation as well as an operating
   point (rule 5).
9. **Within a cruise condition, all six geometries are on one generation**, so F, H and W can
   now be compared with B, C and M there. The published-mesh cruise values (+50.1 to +56.5 counts
   for F, H and W against B) are withdrawn; see the release notes above.

Results are in `report/argus_rans_validation_report.pdf` and `data/`; this README deliberately carries no results
table. Where the mesh generation is recorded:

1. `data/run_cards/SW?_r2_trim.json`: the six Condition CR trims, all wake-refined; and
   `data/run_cards/CMP?_welded_trim.json` and `LC?_welded_trim.json`: the twelve cruise trims, all
   welded. Each card's `mesh.mesh_generation` block records how the mesh was built (192 ranks,
   method `scotch`, machine); it is **not** a generation label. The solve split is in
   `decomposition` and `hpc` (128 ranks, `hierarchical`).
2. `data/rans_forces.json`: every row under `cases` carries `mesh_generation`, `published`,
   `wake_refined` or `welded`; the six `SW?_trim` rows are the `wake_refined` ones, and the twelve
   cruise `_trim` rows and the fifteen `<case>_welded_leg<n>` rows are `welded`, each naming its
   cluster case in `source_case`. Three provenance
   entries in each of those six still cite the published case archive (guide 07, section 5.5);
   take results from `trims`. The summary rows under
   `trims` mostly do not (M's Condition CR row does); each names its case in `trim_case`, and the
   label is on that case. In the Condition CR block of `trims`, the fields
   `superseded_trim_entry` (all six geometries) and `superseded_published_row` (M) keep the
   earlier published-mesh values for the record. Do not quote them.
3. `data/span_efficiency.json`: its `frame` block and all six geometries carry
   `mesh_generation`, spelt `wake-refined`.
4. `data/rans_forces_current.json` and `data/case_provenance_current.json` carry
   `mesh_generation` on every row. Every trim row in both is the current case: the force values
   equal `data/rans_forces.json` row for row, and the eighteen trim rows of the provenance file
   were re-read from the current cases (its `current_cases_2026_10_05` block names the sources).
   Replaced values are kept under `superseded_entry`.

Never infer a generation from a missing label. Every trim, at all three operating points, has
a run card.

**Quote window means.** The quotable Condition CR drag of a geometry is the window mean on the
`REPORT THESE` line quoted in the `notes` field of its run card. The window is the last 200
force samples, the one the solve job's convergence gate tests
([docs/06_trim_to_target_cl.md](docs/06_trim_to_target_cl.md), section 6.6.1), and the window
mean is the average over it. The card's `results.CD_counts` is not the window mean
([docs/07_postprocessing_and_report.md](docs/07_postprocessing_and_report.md), section 9.2).
For B, from the repository root:

```bash
grep -o 'REPORT THESE: Cd [0-9.]*' data/run_cards/SWB_r2_trim.json
```

**What you should see:** one line, `REPORT THESE: Cd 0.0188031` (C_D on the DSO basis, 188.031
counts), while the same card's `results.CD_counts` reads 188.030265.

**If it does not:** if it prints `No such file or directory` instead, you are not in the
repository root: run `cd "$REPO"` and repeat the command.

## Geometry letters

Derived from the `candidate` and `geometry.sha256` fields of each run card in
`data/run_cards/`, and checked against the sha256 checksums in `registry/candidates.yaml` and the
geometry table in the report.

| Letter | Registry key | Concept | Design-study delivery | Meshing surface sha256 (first 8) | Condition CR run card |
|---|---|---|---|---|---|
| B | `baseline` | rigid baseline | none (reference wing) | `20bd6acd` | `data/run_cards/SWB_r2_trim.json` |
| C | `cte_i002_c04` | continuous trailing-edge camber | 29 July 2026 | `790c3ded` | `data/run_cards/SWC_r2_trim.json` |
| F | `cfft_b02_c01` | continuous trailing-edge camber | 19 August 2026 | `82bf2e1e` | `data/run_cards/SWF_r2_trim.json` |
| H | `chc_g02_c06` | conventional hinged | 19 August 2026 | `18329aa3` | `data/run_cards/SWH_r2_trim.json` |
| M | `mcv2_i002_c01` | continuous trailing-edge camber | 29 July 2026 | `dae22d0c` | `data/run_cards/SWM_r2_trim.json` |
| W | `cffw_b01_c01` | distributed twist | 19 August 2026 | `200fa64c` | `data/run_cards/SWW_r2_trim.json` |

1. A case name is the condition prefix, the letter, then a suffix: `SWB_trim` is B at Condition
   CR, `CMPW_trim` is W at early cruise, `LCH_trim` is H at late cruise.
2. The August delivery is the design study's reselection after it moved from near-field to
   far-field induced drag. The registry describes M as a RANS geometry-validation case rather than
   a feasible optimum; it was kept for traceability.
3. **The registry's status notes contradict this table for three candidates.** The `STATUS`
   text under `chc_g02_c06`, `cfft_b02_c01` and `cffw_b01_c01` in `registry/candidates.yaml`
   names the wrong case (it says `SWW_trim`, `SWH_trim` and `SWF_trim` respectively), and the
   `cte_i002_c04` entry's `description` still says C was not run. The registry's sha256
   values are correct. Where the registry's prose disagrees, trust each run card's `candidate`
   and `geometry.sha256` fields.
4. B is the as-delivered wing (older lineage: 9 distinct aerofoils across 19 defining stations).
   All five candidates are built on a corrected lineage. The report's section "Two lineages of
   the same wing" measures the geometric difference (`data/lineage_measured.json`): at most
   0.407 percentage points of t/c and 0.209 % of chord in surface ordinate, and zero at and
   inboard of eta 0.450. It does **not**
   estimate the difference's effect on drag. That effect is common to all five candidates, so it
   drops out of a candidate-to-candidate difference (to first order: each wing responds from its
   own flow), but it is present, unquantified, in every candidate-minus-B drag change.
5. The registry also holds an entry called `baseline_corrected`: it is **not** B, and no case
   was meshed on it.
6. The `.vsp3` design file each surface is built from, and its checksum, are listed in
   [docs/02_geometry_from_openvsp.md](docs/02_geometry_from_openvsp.md), section 3.
7. For F, H and W the checksum in the table is that of the registered, welded surface, which
   their trim solve cases (names ending `_trim`) and run cards carry, and which their cruise
   meshes were built from. Their Condition CR mesh cases on HPC12 hold an earlier copy of the
   same surface, closed but not welded
   ([docs/02_geometry_from_openvsp.md](docs/02_geometry_from_openvsp.md), Step 7 not applied).
   See item 8 of [Known gaps in this release](#known-gaps-in-this-release).

## Quickstart

Read the guides in this order.

1. [docs/01_setup.md](docs/01_setup.md): install OpenFOAM-org 12 on your machine and on HPC12, set up Python, and check the environment.
2. [docs/02_geometry_from_openvsp.md](docs/02_geometry_from_openvsp.md): check the `.vsp3` design files, build each wing surface in metres, and pass `surfaceCheck`.
3. [docs/03_build_a_case.md](docs/03_build_a_case.md): set up baseline B at Condition CR on your own machine, from the surface you built for B, and go through each meshing step before any cluster step.
4. [docs/04_send_to_hpc.md](docs/04_send_to_hpc.md): build the production case on your machine, then copy it to HPC12, where the full-size meshes are built.
5. [docs/05_submit_and_run.md](docs/05_submit_and_run.md): submit the mesh job, make the solve case, submit the solve legs with `qsub`, watch them and restart them.
6. [docs/06_trim_to_target_cl.md](docs/06_trim_to_target_cl.md): trim each geometry to its target C_L.
7. [docs/07_postprocessing_and_report.md](docs/07_postprocessing_and_report.md): forces, run cards, figures and the report.

**The HPC12 route in brief** (each step is explained in
[docs/05_submit_and_run.md](docs/05_submit_and_run.md)):

1. `hpc/mesh_hpc12.pbs` meshes the case at 192 ranks. Its header asks for queue `fpt-large`,
   6 nodes x 32 cores (`typej`) and 12 hours.
2. `hpc/stage_solve_case.sh` copies the finished mesh case into a new solve case at 128 ranks
   (`hierarchical`, 8 x 2 x 8, as the delivered Condition CR solves), checks the copy file by
   file and prints ready-to-run `qsub` lines. It needs only bash and standard command-line
   tools: no OpenFOAM and no Python.
3. `hpc/solve_hpc12.pbs` runs one solve leg per job. You submit every leg yourself, each with
   `ARGUS_GATE` in its `-v` list, so that the leg checks the post-processing block before it
   starts.

**Install the post-processing block before a case leaves your machine.** A case built from
`recipes/` alone writes the force coefficients (`forceCoeffs1`, plus a few extras in the cruise
recipes), but no Trefftz planes, no wing-surface data and no drag split.
`scripts/install_postpro.py` writes `system/argusPostPro` from `hpc/templates/argusPostPro` and
adds one include line to `system/controlDict`. The solver then also writes the pressure and
viscous parts of the forces, the wing-surface fields (p, cp, y+ and the wall-shear-stress
vector), four Trefftz planes behind the wing, and the residuals.
`scripts/assert_case_postpro.py --before` checks a case: a bare recipe gives
`1 of 8 requirements satisfied`, an installed case `8 of 8`. The install needs numpy and the
wing surface as an ASCII STL, `constant/triSurface/wing.stl`. When and where to run both:
[docs/04_send_to_hpc.md](docs/04_send_to_hpc.md) and
[docs/05_submit_and_run.md](docs/05_submit_and_run.md); what the block writes:
[docs/07_postprocessing_and_report.md](docs/07_postprocessing_and_report.md).

### Check that your copy of a recipe is intact

Before building anything, check the recipe folder in your clone:

```bash
cd "$REPO/recipes/L11_wallResolved_CR"
sha256sum -c --quiet MANIFEST.sha256 && echo "recipe intact"
```

**What you should see:** `recipe intact` and nothing else. The same check passes on all four
recipe folders.

**If it does not:** any line containing `FAILED` names a file that differs from the shipped
recipe (`FAILED open or read` means the file is missing), and `recipe intact` is not printed. If
you cloned the repository with git, restore the named file from inside the same recipe folder
(this example restores `system/fvSchemes`; use the path printed before `FAILED`):

```bash
git checkout -- system/fvSchemes
```

Otherwise download the repository again. Then repeat the check.

Run this check on the recipe folders in your clone, not on a case. In a case copied from a
recipe, `system/controlDict` reports `FAILED` once the post-processing block is installed, and
`0.orig/U` too once you set the angle of attack. That is expected, and the mesh job allows for
it ([docs/04_send_to_hpc.md](docs/04_send_to_hpc.md)).

For the two `L11` recipes there is a second check, which compares the recipe's dictionaries with
the settings pinned in `hpc/v5_recipe.json`. It needs only the Python standard library. Run it
from the repository root:

```bash
cd "$REPO"
python3 scripts/verify_recipe.py recipes/L11_wallResolved_CR
```

**What you should see:** the last line is `RECIPE OK`, after lines that include
`manifest: 19 files verified` and `pinned scalars compared: 14`.

**If it does not:** each problem is printed on a `FAIL:` line, and the last line reads
`RECIPE REFUSED: N problems`. The two cruise recipes are always refused by this check, by
design, because its pinned layer settings are the `L11` ones; use the manifest check for them.

## What is not in this repository

> **NO GEOMETRY IS SHIPPED.** There is no `.vsp3` and no STL file anywhere in this repository.
> You download the `.vsp3` design files and build the wing surfaces yourself.

1. **Geometry.** The six wings come as OpenVSP `.vsp3` files: F, H and W from the design study's public GitHub repository, https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow (folder `results/geometries/low_speed/`); B, C and M, which that repository does not hold, from the ARGUS project drive, `\\tudelft.net\staff-umbrella\ARGUS Morphing wing CFD\argus_hpc12_backup\geometry_vsp3\`, which holds all six with a `SHA256SUMS` file.
   Check each downloaded `.vsp3` against the checksum table in
   [docs/02_geometry_from_openvsp.md](docs/02_geometry_from_openvsp.md), section 3. Then build
   each meshing surface with the same guide: its section 11 gives the commands that reproduce
   all six registered surfaces bit for bit, and section 3 lists their full sha256 checksums
   (the first 8 characters are in the geometry-letter table above). Compare by checksum, not by
   file name: a surface you build is named differently from the registry's file. The surfaces
   are half models in metres (52.6 to 53.3 MB each), with the root extended 20 mm through the
   symmetry plane and capped, so a file's own area is not the wing's wetted area. For B the
   registry records the file's area as 1.3497 m² against a wetted area of 1.2794 m², 5.49 %
   more. Take wetted area from the solver's wing patch, never from the STL. Do not follow the
   `route` text in `registry/candidates.yaml`: for H, F and W it leaves out `--write` on
   `close_wing_surface.py` and the `weld_surface.py` step, and followed literally it leaves the
   root open (02, section 11, note 1).
2. **Meshes.** None are shipped. The full-size meshes (99,111,506 to 99,208,470 cells for the
   six Condition CR trims, field `mesh.cells` of each run card) do not fit on a workstation. They
   are built on HPC12 from `recipes/` with `hpc/mesh_hpc12.pbs`
   ([docs/04_send_to_hpc.md](docs/04_send_to_hpc.md),
   [docs/05_submit_and_run.md](docs/05_submit_and_run.md)).
3. **Raw fields**: time directories, `processor*` directories and `postProcessing/`. None are
   shipped (`.gitignore` excludes them). The numbers harvested from them are in `data/`.
4. **Solver and mesher logs.** None are shipped. What was read from them is recorded in the run
   cards and in `data/`.
5. **Part of the HPC12 campaign tooling.** Shipped: `hpc/mesh_hpc12.pbs`, `hpc/solve_hpc12.pbs`
   and `hpc/stage_solve_case.sh`, which does the step from a finished mesh case to a solve case.
   Not in this repository: the job that builds OpenFOAM-org 12 in your HPC12 home directory; the
   scripts that reconstructed and archived finished cases; the campaign scheduler and watcher;
   and the builder of the trim cases.
6. **The source repository's working records.** Some shipped files (`registry/candidates.yaml`,
   `recipes/README.md`, script comments) cite internal record identifiers and source-repository
   paths (`hpc/recipes/`, `geometry/derived/`, `cases/validation/`, `dso_reference/`) that do not
   exist here. In `recipes/README.md`, read `hpc/recipes/` as `recipes/`. Some script output
   does the same: the explanation lines that `scripts/assert_case_postpro.py` prints under each
   result name scripts that are not shipped.
7. **The design-study deliveries** (VSPAERO inputs and outputs). The VLM results extracted from
   them are in `data/vlm.json`. `report/scripts/extract_vlm.py`, which extracts them, needs
   those deliveries in a `dso_reference/` folder, which is not shipped.

## Known gaps in this release

1. **The report source does not compile in place, and two of its inputs are older than the
   shipped PDF.** `report/src/argus_rans_validation_report.tex` looks for its figures in `fig/` next to itself, but
   they are in `report/fig/`. Build it as
   [docs/07_postprocessing_and_report.md](docs/07_postprocessing_and_report.md) says, following
   sections 8.2 to 8.4 in full: copy into a separate build directory, regenerate the figures,
   rewrite the span-efficiency paragraph, then build. `latexmk` then finishes with 34 pages and
   no undefined references. Stopping after section 8.2 also gives 34 pages and no undefined
   references, but brings back M's superseded Condition CR result from two shipped inputs:
   1. `report/fig/fig_trimmed_performance.pdf`: C_D 191.25 counts, a drag change of +0.240
      counts and a trim angle of 1.252 deg for M;
   2. `report/fig/cl_alpha_condition_CR.pdf`: a trim angle of 1.2520 deg for M.

   The shipped `report/argus_rans_validation_report.pdf` has M's current values in these places: 188.26 counts,
   +0.232 counts and 1.257 deg in the first, 1.2574 deg in the second.
   Never quote those two files as shipped. Or use the shipped `report/argus_rans_validation_report.pdf`.
2. **`recipes/L11_wallResolved_WT` is refused by `hpc/mesh_hpc12.pbs`**, even as a fresh copy.
   Its `liftDir` and `dragDir` are written to 6 significant figures, which gives an axis lag of
   -3.66e-06 deg against the job's tolerance of 1e-6 deg (`dragDir implies alpha 2.102797,
   Uinf implies 2.102800`). The CR, early-cruise and late-cruise recipes pass. WT is not one of
   the three delivered conditions. To mesh it anyway, set its angle first as
   [docs/05_submit_and_run.md](docs/05_submit_and_run.md), section 5 (Gate A), describes.
3. **`scripts/assert_case_postpro.py --after` reports false negatives on OpenFOAM-12 output.** On
   a correctly run case it ends `---- 5 of 8 requirements satisfied ----`, with
   `wing_surface_p`, `wing_surface_tau_VECTOR` and `wing_surface_yplus` reported `MISSING`. It
   looks for the field names only at the start of each wing-surface file, and OpenFOAM-12 writes
   them after the points and faces. Check those fields directly, as
   [docs/07_postprocessing_and_report.md](docs/07_postprocessing_and_report.md) shows. The usage
   lines at the top of the script (not its `--help` output) also show a `--time` option that
   does not exist (`error: unrecognized arguments`). Its `compressible_cp_frame` requirement
   cannot fail on an OpenFOAM-12 case, so it cannot catch a wrong pressure-coefficient frame on
   a cruise case.
4. **You submit every solve leg by hand.** `hpc/solve_hpc12.pbs` tries to submit its own next
   leg with `qsub` from inside the job. On HPC12 that fails (`could not connect to trqauthd`),
   and it would not pass `ARGUS_GATE` on anyway. Submit each leg yourself with the full `-q`,
   `-l` and `-v` ([docs/05_submit_and_run.md](docs/05_submit_and_run.md), section 10a). For a
   cruise case, the job's own raise of `endTime` by 1500 (10000 to 11500) is not a multiple of
   the cruise `writeInterval` 2000, so the job refuses it: give `ARGUS_ENDTIME` a multiple of
   2000 instead.
5. **The post-processing gate runs only if you point the job at it.** Without `ARGUS_GATE` in
   the `-v` list, `hpc/solve_hpc12.pbs` looks for `$HOME/argus_hpc12/assert_case_postpro.py`,
   which does not exist on a new account, prints
   `NOTE: <path> absent, post-processing gate NOT run` and runs the leg unchecked. Add
   `ARGUS_GATE=$REPO/scripts/assert_case_postpro.py` to every leg's `-v` list, as the lines that
   `hpc/stage_solve_case.sh` prints do. A real PBS job has received `ARGUS_GATE` through `-v`
   on HPC12 and ran the check (2026-09-24, `8 of 8 requirements satisfied`).
6. **Some comments and pinned values still describe an older set-up.** The header comments of
   `hpc/solve_hpc12.pbs` describe 192 solve ranks (hierarchical 8 x 3 x 8) and an
   incompressible `endTime` of 4500, and the `solve` and `ranks` entries of
   `hpc/v5_recipe.json` give 192 solve ranks, `endTime` 4500 and walltime 24:00:00. What runs
   is set by the case: the job reads the rank count and `endTime` from it, and refuses an
   allocation with fewer slots than the case has ranks. `scripts/splice_postpro.py` marks a case
   quotable only once it has reached iteration 4500, so it reports a Condition CR case that
   converged at `endTime` 4000 as not quotable. The verdict that counts is the case's
   `.converged` marker ([docs/07_postprocessing_and_report.md](docs/07_postprocessing_and_report.md)).
7. **Not yet run on HPC12** (NOT VERIFIED there):
   1. a comparison showing the post-processing block does not change a production mesh. A
      192-rank mesh job with the block installed has run on HPC12, but there was no mesh of the
      same surface without it to compare. At tutorial resolution on a local machine, the block
      did not change the mesh;
   2. the writing step of `scripts/install_postpro.py` at the HPC12 login prompt. Its
      derivations ran there with `--dry-run`;
   3. the manual fallback for a solve case without the block: copying `system/argusPostPro`
      into it and adding the include line to `system/controlDict` by hand
      ([docs/05_submit_and_run.md](docs/05_submit_and_run.md), section 7a).
8. **The Condition CR F, H and W meshes were built from the unwelded surface.** In the F, H and
   W mesh cases of the published generation (all three conditions) and of the wake-refined
   generation (Condition CR), `constant/triSurface/wing.stl` is the closed but unwelded version
   of the registered surface; the welded cruise meshes are not affected (sub-item 2). It has the same 195,266 triangles, and its points differ from the registered file
   only along the root seam, by at most 7.8e-07 m (metre frame of `wing.stl`). Exact point
   matching finds 796 open edges there. The mesh job's surface gate rounds points to the
   nearest micrometre and accepts such a surface
   ([docs/02_geometry_from_openvsp.md](docs/02_geometry_from_openvsp.md), section 10).
   `data/rans_forces.json` and the run cards record the registered, welded checksum for the
   Condition CR cases, and the trim solve cases (names ending `_trim`) carry the welded file. The
   B, C and M mesh cases hold the registered files in every generation.

   **What it changed, measured from the HPC12 meshing logs.** The open edges are the root rim,
   where the 20 mm root extension meets the wing. The two copies of each rim point lie 0.07 to
   0.78 micrometres apart, and about half of the rim lies inside the mesh domain (y from 0 to
   4.68 mm). At the production surface level that gap lets snappyHexMesh's first cell-removal
   step keep the wing interior: it did so in every F, H and W mesh build (Condition CR on both
   generations, early cruise, late cruise) and in none of the B, C and M builds. A later
   removal step deletes the interior again, so **no delivered mesh is flooded**: the wing face
   counts match B, C and M to 0.01 % and the wing patch area matches the surface to 0.004 %.
   But with the 40,000,000-cell cap of the published mesh, refining that leaked interior used
   up the cell budget (43.5 to 44.3 million cells after the first shell refinement step,
   against 29.6 to 33.8 million for B, C and M) and starved the near field and the wake of F, H
   and W. That, not feature refinement, is why those meshes were starved.
   1. **Condition CR:** the wake-refined rebuild removed the effect. F, H and W minus B moved
      from +5.733, +6.580 and +5.734 counts (published) to -0.815, -0.009 and -0.868 counts
      (wake-refined). The seam still leaves local marks at the root of the wake-refined F, H
      and W meshes (a few points up to 0.58 mm below the symmetry plane, and higher y+ within
      0.5 mm of the root), with a measured effect of a few hundredths of a count. The delivered
      Condition CR results stand.
   2. **Early and late cruise:** rebuilt from the registered surfaces for all six geometries (the
      welded generation). In those mesh cases `constant/triSurface/wing.stl` was hashed and matches
      the registered file for every geometry. See the release notes above.
   3. **Future meshes:** build F, H and W only from the registered, welded surfaces, and require
      `surfaceCheck` to report `Surface is closed` before meshing
      ([docs/02_geometry_from_openvsp.md](docs/02_geometry_from_openvsp.md), Step 8). The mesh
      job's own surface gate cannot see this defect.
9. **Two small loose ends in the report's hinge-moment section.** The bin-count sweep sentence
   (report section 5.6, "SWB_trim gives M_h = 0.9375, ... 0.8630 N m") does not state which
   mesh generation or iteration it was run on, and its 1600-bin value (0.8630 N m) matches
   neither Table 18's B total (0.8668 N m) nor its pressure term (0.8633 N m). The q ratio in
   section 5.7 is printed as 1.24431, while `report/fig/rootbending_values.csv` implies
   1.24425. Neither changes a table, a percentage or a verdict.

## Citation and contact

Cite this repository with [CITATION.cff](CITATION.cff). The RANS validation workflow was
developed by Tyler Buchanan at TU Delft for the ARGUS project. The geometries and vortex-lattice
results it validates come from the ARGUS Aerodynamic Morphing-Wing Workflow by Liming Zheng,
https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow; cite that work with its own
`CITATION.cff` when you use the design study. Project-specific questions should go through the
current ARGUS work-package lead.
