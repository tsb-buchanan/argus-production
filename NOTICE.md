# Project provenance and third-party notice

## Project provenance

The ARGUS RANS validation workflow was developed by Tyler Buchanan at Delft University of
Technology as part of the ARGUS research project. The repository records the OpenFOAM case
recipes, HPC12 mesh and solve jobs, trim procedure, post-processing, run cards, audited summary
data and validation report produced for that project.

It validates the design study of the ARGUS Aerodynamic Morphing-Wing Workflow by Liming Zheng
(Delft University of Technology), https://github.com/Liming-Zheng/ARGUS_Aerodynamic_Workflow.
Please cite this repository using CITATION.cff, and cite that workflow with its own
CITATION.cff when you use its geometries or vortex-lattice results. References to ARGUS, Delft
University of Technology, or project collaborators describe provenance only and do not imply
endorsement of derived work.

## License scope

Original software and documentation in this repository are released under the MIT
License in LICENSE.

The license does not relicense third-party software, source data, geometry, or branding. In
particular:

- OpenFOAM, OpenVSP and VSPAERO are external software and retain their own licenses. None of
  them is redistributed here; the files under `recipes/` and `template_case/` are input
  dictionaries for OpenFOAM.
- The design-study geometries (`.vsp3`) are not included. The vortex-lattice results in
  `data/vlm.json` and `data/dso_vlm_claims.json` are extracted from the design study's
  deliveries and keep that provenance.
- NASA-derived geometry and reference material (the EET AR-12 wing, NASA TP-1580) retain their
  original provenance and applicable terms.
- The presentation in `presentations/` uses the TU Delft slide template; its branding is not
  covered by the license.
- Names, logos, and trademarks of ARGUS, TU Delft, NASA, OpenFOAM, OpenVSP, and project
  partners remain the property of their respective owners.

Users are responsible for checking the applicable terms before redistributing third-party or
derived assets.
