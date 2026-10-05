#!/usr/bin/env python3
"""validate_run_card.py: validate ARGUS run cards against scripts/run_card.schema.json.

Checks each card against the JSON schema and, when the registry is available,
cross-checks that the candidate is registered and the geometry sha256 appears
in that candidate's registry entry (or baseline's, for baseline geometry).

Usage:
    python3 scripts/validate_run_card.py cases/<...>/run_card.json [...]
    python3 scripts/validate_run_card.py --self-test
"""

import argparse
import json
import sys
from pathlib import Path

import jsonschema
import yaml

REPO = Path(__file__).resolve().parents[1]
SCHEMA = REPO / "scripts/run_card.schema.json"
REGISTRY = REPO / "registry/candidates.yaml"


def collect_sha256(node, out):
    if isinstance(node, dict):
        for k, v in node.items():
            if k == "sha256" and isinstance(v, str):
                out.add(v)
            else:
                collect_sha256(v, out)
    elif isinstance(node, list):
        for v in node:
            collect_sha256(v, out)


def registry_errors(card, registry):
    errs = []
    cand = card.get("candidate")
    entry = registry.get(cand) or registry.get("benchmarks", {}).get(cand)
    if entry is None:
        return [f"candidate '{cand}' not in registry/candidates.yaml "
                f"(top level or benchmarks section)"]
    shas = set()
    collect_sha256(entry, shas)
    collect_sha256(registry.get("baseline", {}), shas)
    sha = card.get("geometry", {}).get("sha256")
    if sha and sha not in shas:
        errs.append(f"geometry sha256 {sha[:12]}... not registered for "
                    f"'{cand}' (or baseline)")
    return errs


def cross_field_errors(card):
    """Constraints JSON Schema cannot express, because they relate two fields.

    THE RANK-COUNT DUPLICATION. `decomposition` and `hpc` both carry the solve rank count
    and decomposition method: `decomposition` because it is present for laptop AND cluster
    runs, `hpc` because the cluster generator writes its cards from the job log. A NUMBER
    RECORDED IN TWO PLACES WITHOUT A CHECK IS A NUMBER THAT WILL EVENTUALLY DISAGREE, and a
    disagreement here is exactly the unstated-frame failure both blocks exist to prevent.
    `decomposition` is canonical; when `hpc` is non-null the two must agree."""
    errs = []
    dec, hpc = card.get("decomposition"), card.get("hpc")
    if isinstance(dec, dict) and isinstance(hpc, dict):
        if hpc.get("solve_ranks") is not None and dec.get("n_subdomains") is not None \
           and hpc["solve_ranks"] != dec["n_subdomains"]:
            errs.append(f"hpc/solve_ranks {hpc['solve_ranks']} disagrees with "
                        f"decomposition/n_subdomains {dec['n_subdomains']}; "
                        f"decomposition is canonical")
        a, b = hpc.get("solve_decomposition_method"), dec.get("method")
        if a and b and a != b:
            errs.append(f"hpc/solve_decomposition_method {a!r} disagrees with "
                        f"decomposition/method {b!r}; decomposition is canonical")
    # A PAIR MEMBER MUST MATCH ITS PARTNER'S RANK COUNT (D079). Rank count does not change a
    # converged steady answer, but an unstated difference between two halves of a delta is
    # the D068 failure, so it is recorded and checked rather than assumed.
    if isinstance(dec, dict) and dec.get("pair_partner"):
        if dec.get("pair_partner_n_subdomains") != dec.get("n_subdomains"):
            errs.append(f"pair member rank mismatch: this run {dec.get('n_subdomains')} vs "
                        f"partner {dec.get('pair_partner_n_subdomains')} "
                        f"({dec['pair_partner']})")
    # A SERIAL MESH MUST SAY SO IN method, not merely leave ranks null (schema requires it
    # in prose; enforced here).
    mg = card.get("mesh", {}).get("mesh_generation")
    if isinstance(mg, dict) and mg.get("ranks") is None:
        if "serial" not in str(mg.get("method", "")).lower():
            errs.append("mesh/mesh_generation: ranks is null but method does not say "
                        "'serial'; a null rank count is only meaningful if the method "
                        "states the mesh was built serially")
    if isinstance(mg, dict) and mg.get("inferred") and not mg.get("inferred_from"):
        errs.append("mesh/mesh_generation: inferred is true but inferred_from is absent; "
                    "an inferred value without its evidence is not auditable")
    return errs


def validate(card, validator, registry):
    errs = [f"{'/'.join(str(p) for p in e.absolute_path) or '(root)'}: {e.message}"
            for e in sorted(validator.iter_errors(card), key=lambda e: list(e.absolute_path))]
    errs += cross_field_errors(card)
    if not errs and registry is not None:
        errs += registry_errors(card, registry)
    return errs


def example_card():
    """Schema-conformant EXAMPLE (not a run record) for --self-test."""
    return {
        "schema_version": 1,
        "candidate": "baseline",
        "case": {"path": "cases/baseline/CR_M0p10_CL0p428_2d_eta0p80",
                 "dimensionality": "2D", "condition": "CR", "eta": 0.80},
        "geometry": {"file": "geometry/derived/eet_section/eet_cruise_faired.dat",
                     "sha256": "0" * 64},
        "provenance": {"git_commit": "0123abc", "git_dirty": False,
                       "openfoam_version": "openfoam-org-7", "date": "2026-07-21"},
        "mesh": {"generator": "blockMesh+snappyHexMesh+extrudeMesh", "cells": 120000,
                 "max_skewness": 2.1, "max_non_orthogonality": 60.2,
                 "yplus": {"min": 0.2, "mean": 0.8, "max": 1.6},
                 "layers": {"n_surface_layers": 28, "first_layer_thickness_m": 8e-6,
                            "expansion_ratio": 1.2, "coverage_fraction": 0.98},
                 "checkMesh_pass": True,
                 # mesh_generation, wall_treatment and hpc became MANDATORY in commit
                 # d1abb50, whose message asserted "Every card now names nut, k and
                 # omega". No generator was changed and THIS FIXTURE WAS NOT UPDATED
                 # EITHER, so the validator's own self-test has failed from that commit
                 # onward with nothing surfacing it: the one check that would have
                 # caught the schema/card divergence was itself broken by the same
                 # commit. A fixture that cannot satisfy its own schema is not a
                 # fixture (GEO-089). Repaired 2026-08-21, ARG-170.
                 # ranks is part of the MESH RECIPE, not the solve: snappyHexMesh
                 # does not produce the same mesh at different rank counts.
                 "mesh_generation": {"ranks": 12, "method": "scotch", "inferred": False,
                                     "inferred_from": "authored by the case builder, "
                                                      "which writes decomposeParDict "
                                                      "and this block from one constant",
                                     "machine": "laptop-WSL2"}},
        "model": {"solver": "simpleFoam", "turbulence_model": "kOmegaSST",
                  "wall_treatment": "low-Re resolved"},
        "bc": {"airfoil": {"U": "noSlip"}, "inlet": {"U": "fixedValue 34.0 at alpha"}},
        "numerics": {"schemes_summary": "steadyState, linearUpwind, SIMPLEC",
                     "convergence_criteria": {"p": 1e-5, "U": 1e-6, "k|omega": 1e-6}},
        "convergence": {"iterations": 2400,
                        "final_residuals": {"p": 9e-6, "U": 8e-7, "k": 6e-7, "omega": 4e-7},
                        "converged": True},
        "conditions": {"U_mag_m_s": 34.0, "alpha_deg": 1.234, "nu_m2_s": 1.46e-5,
                       "Re_ref": 7.9e5, "mach_nominal": 0.10, "rho_ref_kg_m3": 1.225,
                       "reference": {"Aref_m2": 0.3404, "lRef_m": 0.3404,
                                     "CofR_m": [0.0851, 0.0, 0.0],
                                     "convention": "per unit span: Aref = chord x 1 m span"}},
        "results": {"CL": 0.4283, "CD": 0.0091, "CD_pressure": 0.0034,
                    "CD_viscous": 0.0057, "Cm": -0.083},
        "trim": {"target_CL": 0.428277635108, "tolerance": 1e-4,
                 "alpha_trim_deg": 1.234, "converged": True,
                 "history": [{"alpha_deg": 0.5, "CL": 0.35}, {"alpha_deg": 1.5, "CL": 0.455},
                             {"alpha_deg": 1.234, "CL": 0.42831}]},
        "cost": {"wall_clock_s": 1800, "n_cores": 8, "core_hours": 4.0,
                 "machine": "laptop-WSL2"},
        "wall_treatment": {"nut": "nutLowReWallFunction", "k": "fixedValue",
                           "omega": "omegaWallFunction", "yplus_max": 1.6,
                           "yplus_convention": "cell-centre"},
        "hpc": None,
        "notes": "EXAMPLE CARD (validator self-test), not a run record.",
    }


def self_test(validator):
    good = example_card()
    errs = validate(good, validator, None)
    if errs:
        print("self-test FAIL: example card should validate:")
        print("\n".join(f"  {e}" for e in errs))
        return 1
    bad_cases = []
    b = json.loads(json.dumps(good)); del b["geometry"]["sha256"]
    bad_cases.append(("missing geometry.sha256", b))
    b = json.loads(json.dumps(good)); b["model"]["wall_treatment"] = "wall functions"
    bad_cases.append(("wall functions without decision_ref", b))
    b = json.loads(json.dumps(good)); b["results"].pop("CD_viscous")
    bad_cases.append(("missing CD split", b))
    b = json.loads(json.dumps(good)); b["typo_field"] = 1
    bad_cases.append(("unknown top-level field", b))
    for name, card in bad_cases:
        if not validate(card, validator, None):
            print(f"self-test FAIL: '{name}' should be rejected")
            return 1
    print("self-test PASS (example accepted; 4 corrupted variants rejected)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cards", nargs="*", type=Path)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--no-registry-check", action="store_true")
    args = ap.parse_args()

    schema = json.loads(SCHEMA.read_text())
    validator = jsonschema.Draft202012Validator(schema)

    if args.self_test:
        sys.exit(self_test(validator))
    if not args.cards:
        ap.error("no run cards given (or use --self-test)")

    registry = None
    if not args.no_registry_check and REGISTRY.exists():
        registry = yaml.safe_load(REGISTRY.read_text())

    n_bad = 0
    for path in args.cards:
        errs = validate(json.loads(path.read_text()), validator, registry)
        if errs:
            n_bad += 1
            print(f"FAIL {path}")
            print("\n".join(f"  {e}" for e in errs))
        else:
            print(f"PASS {path}")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
