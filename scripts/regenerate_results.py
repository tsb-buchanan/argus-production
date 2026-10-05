#!/usr/bin/env python3
"""regenerate_results.py: rebuild every results/ artifact from the case data.

Single entry point for project record-keeping rule 3 ("all plots regenerated
by script"). Each step is a collector that reads case output and writes CSV/JSON
plus its figures; nothing here computes physics of its own.

Steps that need an external tool (XFOIL) are skipped with a notice unless the
tool is available, so this runs end-to-end on a clean checkout.

Usage:
  python3 scripts/regenerate_results.py [--xfoil PATH] [--skip-db]
"""

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
XFOIL_DEFAULT = Path.home() / "tools/xfoil-6.99/xfoil"

STEPS = [
    # FIRST, and deliberately so: the derived STLs are git-ignored (*.stl), so a
    # clean checkout has the recipe and the checksum but not the binary. This
    # step regenerates both production OMLs and FAILS LOUDLY on any checksum
    # mismatch against registry/candidates.yaml. The ruling (D050) is that the
    # binary is a cache and the recipe plus checksum is the artifact, which only
    # holds as long as this gate passes. If it ever fails, the geometry is not
    # reproducible and the promotion is void until it is.
    ("3D OML regeneration gate (D050)",
     [sys.executable, "scripts/regenerate_oml.py", "--check"]),
    ("TMR Layer-1 forces + Richardson + figures",
     [sys.executable, "validation_data/naca0012_tmr/collect_ogrid.py"]),
    ("TMR Layer-1 figures",
     [sys.executable, "validation_data/naca0012_tmr/plot_ogrid.py"]),
    ("TMR Layer-1 verification memo",
     [sys.executable, "validation_data/naca0012_tmr/make_layer1_memo.py"]),
    ("EET 2D grid-convergence family",
     [sys.executable, "scripts/collect_2d_grid.py",
      "--glob", "cases/baseline/CR2d_grid_*", "--out", "results/eet2d_cr"]),
    ("Task-3 fixed-cl morphed deltas (D019)",
     [sys.executable, "scripts/collect_trim_deltas.py"]),
    ("Task-4 fixed-alpha sweep run cards",
     [sys.executable, "scripts/collect_alpha_sweep.py"]),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--xfoil", type=Path, default=XFOIL_DEFAULT)
    ap.add_argument("--skip-db", action="store_true")
    ap.add_argument("--keep-going", action="store_true",
                    help="report failures instead of stopping at the first one")
    args = ap.parse_args()

    steps = list(STEPS)
    if args.xfoil.exists():
        steps.append(("Task-4 XFOIL code-to-code cross-check",
                      [sys.executable, "scripts/xfoil_crosscheck.py",
                       "--xfoil", str(args.xfoil)]))
    else:
        print(f"NOTE: xfoil not found at {args.xfoil}; skipping the cross-check step "
              f"(build recipe in docs/environment.md). Existing artifacts kept.")
    steps.append(("Results index + RESULTS.md",
                  [sys.executable, "scripts/consolidate_results.py"]))
    steps.append(("Phase 0 report",
                  [sys.executable, "scripts/make_phase0_report.py"]))
    if not args.skip_db:
        steps.append(("DuckDB derived query layer (D016)",
                      [sys.executable, "scripts/build_results_db.py"]))

    failures = []
    for i, (label, cmd) in enumerate(steps, 1):
        print(f"\n=== [{i}/{len(steps)}] {label}")
        p = subprocess.run(cmd, cwd=REPO)
        if p.returncode != 0:
            failures.append(label)
            if not args.keep_going:
                sys.exit(f"FAILED at step {i}: {label}")
    print()
    if failures:
        sys.exit(f"{len(failures)} step(s) failed: {'; '.join(failures)}")
    print(f"all {len(steps)} steps OK")


if __name__ == "__main__":
    main()
