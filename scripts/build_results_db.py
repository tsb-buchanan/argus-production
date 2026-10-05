#!/usr/bin/env python3
"""build_results_db.py: (re)build results/argus.duckdb from the files in place.

D016: DuckDB is a DERIVED QUERY LAYER, not a store of record. The files
(run cards, results CSVs, registry/candidates.yaml) remain the source of
truth; delete the .duckdb at any time and rerun this script.

Loads:
1. run_cards_raw:       read_json_auto over cases/**/run_card.json
2. forces_tmr:          results/tmr_naca0012/ogrid_forces.csv
3. registry_geometries: flattened export of registry/candidates.yaml
                        (candidate, kind, file, sha256, status)
then executes scripts/views.sql (run_cards, registry, forces, runs_full,
convergence).

Example queries (duckdb results/argus.duckdb):
1. All runs per candidate:
     SELECT candidate, case_name, alpha_deg, cells, cl, cd, converged
     FROM run_cards ORDER BY candidate, alpha_deg, cells;
2. Grid-convergence table per family (deltas between successive levels):
     SELECT * FROM convergence WHERE candidate = 'tmr_naca0012';
3. Cost summary:
     SELECT candidate, count(*) AS runs, round(sum(core_hours), 2) AS core_h,
            round(sum(wall_clock_s)/3600, 2) AS wall_h
     FROM run_cards GROUP BY candidate;
4. Registry gate audit (runs whose geometry is not registered):
     SELECT case_name, geometry_file FROM runs_full
     WHERE registry_candidate IS NULL;
"""

import sys
from pathlib import Path

import duckdb
import yaml

REPO = Path(__file__).resolve().parents[1]
DB = REPO / "results/argus.duckdb"


def flatten_registry():
    reg = yaml.safe_load((REPO / "registry/candidates.yaml").read_text())
    rows = []
    def walk(cand, kind, node):
        if isinstance(node, dict):
            if "sha256" in node and node.get("sha256"):
                rows.append({"candidate": cand, "kind": kind,
                             "file": node.get("file"), "sha256": node["sha256"],
                             "status": node.get("status")})
            else:
                for k, v in node.items():
                    walk(cand, k if kind == "" else f"{kind}.{k}", v)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(cand, f"{kind}[{i}]", v)
    for cand, entry in reg.items():
        walk(cand, "", entry)
    return rows


def main():
    DB.parent.mkdir(exist_ok=True)
    DB.unlink(missing_ok=True)
    import os
    os.chdir(REPO)
    con = duckdb.connect(str(DB))
    con.execute(f"""
        CREATE TABLE run_cards_raw AS
        SELECT * FROM read_json_auto('{REPO}/cases/**/run_card.json',
                                     filename=true, union_by_name=true)
    """)
    con.execute(f"""
        CREATE TABLE forces_tmr AS
        SELECT * FROM read_csv_auto('{REPO}/results/tmr_naca0012/ogrid_forces.csv')
    """)
    rows = flatten_registry()
    con.execute("""CREATE TABLE registry_geometries
                   (candidate VARCHAR, kind VARCHAR, file VARCHAR,
                    sha256 VARCHAR, status VARCHAR)""")
    con.executemany("INSERT INTO registry_geometries VALUES (?, ?, ?, ?, ?)",
                    [(r["candidate"], r["kind"], r["file"], r["sha256"], r["status"])
                     for r in rows])
    con.execute((REPO / "scripts/views.sql").read_text())
    n_cards = con.execute("SELECT count(*) FROM run_cards").fetchone()[0]
    n_geo = con.execute("SELECT count(*) FROM registry_geometries").fetchone()[0]
    print(f"built {DB.relative_to(REPO)}: {n_cards} run cards, {n_geo} registered "
          f"geometries, views: run_cards, registry, forces, runs_full, convergence")
    con.close()


if __name__ == "__main__":
    main()
