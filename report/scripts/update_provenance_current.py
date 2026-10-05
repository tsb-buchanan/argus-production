#!/usr/bin/env python3
"""Point data/case_provenance_current.json's eighteen trim rows at the current cases.

WHY. Every row was harvested from the published case archives (argus_hpc12_backup/cases/), so
the six Condition CR trims and the twelve cruise trims described cases that no result in the
report now comes from. This re-reads each trim from the archive of the case that IS the
result: argus_hpc12_backup/r2/cases/SW<L>_trim_r2 at Condition CR and
argus_hpc12_backup/cruise_welded/cases/<prefix><L>_welded at cruise, through the harvester's
own harvest_one, so the per-file reading rules are the same ones the rest of the file used.

WHAT THE CURRENT ARCHIVES DO NOT CARRY, AND WHERE IT COMES FROM INSTEAD. Neither archive has a
CASE_PROVENANCE.json (only the published case builder wrote one). The geometry identity is
therefore taken from each trim's run card, whose geometry block was hashed from the case
(results/run_cards_r2, results/run_cards_welded), and the recipe from the case's own
RECIPE_FILES.json and MANIFEST.sha256. Each row says which source served each field.

Every row also gets mesh_generation from data/rans_forces.json, so the alpha-bracket rows, which
stay as harvested, are labelled too. The previous row is kept under superseded_entry; the file
is backed up first and written through a temporary name. Report first, exit after.
"""
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from report_paths import data_root  # noqa: E402
import harvest_case_provenance as H  # noqa: E402

D = data_root() / "data"
REPO = HERE.parents[3]
U = Path("/mnt/u/argus_hpc12_backup")
CARDS = {"SW": (REPO / "results/run_cards_r2", "%s%s_r2_trim.json"),
         "CMP": (REPO / "results/run_cards_welded", "%s%s_welded_trim.json"),
         "LC": (REPO / "results/run_cards_welded", "%s%s_welded_trim.json")}
# Cruise rows are read from a copy of the four small files each HPC12 solve case holds
# (constant/momentumTransport, system/controlDict, RECIPE_FILES.json, MANIFEST.sha256), fetched
# into the directory given as argv[1]: three of the twelve were not yet on U: when this ran, and
# one source for all twelve beats two. The row records the HPC12 path, which is what was read.
FETCHED = Path(sys.argv[1]) if len(sys.argv) > 1 else None
ARCH = {"SW": (U / "r2/cases", "%s%s_trim_r2", "wake_refined"),
        "CMP": (FETCHED, "%s%s_welded", "welded"),
        "LC": (FETCHED, "%s%s_welded", "welded")}
CLUSTER = "hpc12:/home/scratch/$USER/argus_wfix/solve/%s%s_welded"
# The r2 archives do not carry RECIPE_FILES.json; that absence is recorded, not refused.
ALLOWED_ABSENT = {"SW": {"recipe_files"}, "CMP": set(), "LC": set()}


def main():
    if not U.is_dir():
        sys.exit("REFUSED: %s not mounted" % U)
    if FETCHED is None or not FETCHED.is_dir():
        sys.exit("usage: update_provenance_current.py <dir of fetched cruise case files>")
    p = D / "case_provenance_current.json"
    doc = json.loads(p.read_text())
    cases = doc["cases"]
    gen = {n: e.get("mesh_generation") for n, e in
           json.loads((D / "rans_forces.json").read_text())["cases"].items()}
    done, problems = [], []
    for pre in ("SW", "CMP", "LC"):
        for L in "BCFHMW":
            name = "%s%s_trim" % (pre, L)
            adir = ARCH[pre][0] / (ARCH[pre][1] % (pre, L))
            if not adir.is_dir():
                problems.append("%s: archive %s missing" % (name, adir)); continue
            rec, notes = H.harvest_one(str(adir), name)
            if pre != "SW":
                # harvest_one records the path it READ; for the fetched copy that is a local
                # temporary directory, so every such path is restated as the HPC12 file it is a
                # copy of. Nothing local and temporary may reach the record.
                for k, v in list(rec.items()):
                    if isinstance(v, str) and v.startswith(str(adir)):
                        rec[k] = (CLUSTER % (pre, L)) + v[len(str(adir)):]
            bad = {k: v for k, v in notes.items() if k != "provenance" and v != "OK"
                   and not (v == "ABSENT" and k in ALLOWED_ABSENT[pre])}
            if bad:
                problems.append("%s: %s" % (name, bad)); continue
            card = json.loads((CARDS[pre][0] / (CARDS[pre][1] % (pre, L))).read_text())
            rec["prov_geometry_sha256"] = card["geometry"]["sha256"]
            rec["prov_geometry_source"] = card["geometry"]["file"]
            rec["prov_alpha_deg_in_case"] = card["conditions"]["alpha_deg"]
            rec["prov_speed_m_s"] = card["conditions"]["U_mag_m_s"]
            man = adir / "MANIFEST.sha256"
            if man.is_file():
                rec["prov_recipe_sha"] = hashlib.sha256(man.read_bytes()).hexdigest()
            rec["prov_source"] = {
                "case_files": (str(adir).replace("/mnt/u/", "U:/") if pre == "SW"
                               else CLUSTER % (pre, L)),
                "RECIPE_FILES.json": notes.get("recipe_files"),
                "geometry_alpha_speed": "run card %s" % (CARDS[pre][1] % (pre, L)),
                "turbulence_controlDict_recipe": "the case archive, by harvest_case_provenance.harvest_one",
                "CASE_PROVENANCE.json": "not carried by this archive"}
            rec["mesh_generation"] = ARCH[pre][2]
            if gen.get(name) != ARCH[pre][2]:
                problems.append("%s: rans_forces.json says %s, archive is %s"
                                % (name, gen.get(name), ARCH[pre][2])); continue
            old = cases.get(name)
            if old is not None:
                rec["superseded_entry"] = old.get("superseded_entry", old)
            cases[name] = rec
            done.append(name)
    for n, e in cases.items():
        if n not in done and gen.get(n):
            e["mesh_generation"] = gen[n]
    unlabelled = [n for n, e in cases.items() if not e.get("mesh_generation")]
    print("  trims re-pointed: %d of 18" % len(done))
    for x in problems:
        print("  PROBLEM %s" % x)
    print("  rows %d, unlabelled %d %s" % (len(cases), len(unlabelled), unlabelled[:5]))
    if problems or len(done) != 18:
        print("REFUSED: nothing written")
        return 1
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    bak = p.with_name(p.name + ".bak_" + now.replace("-", "").replace(":", ""))
    shutil.copy2(p, bak)
    doc["current_cases_2026_10_05"] = dict(
        what="the eighteen trim rows re-read from the archives of the cases the results come from",
        sources={"condition_CR": "U:/argus_hpc12_backup/r2/cases/SW<L>_trim_r2",
                 "cruise": "hpc12:/home/scratch/$USER/argus_wfix/solve/<prefix><L>_welded"},
        geometry_from="the trims' run cards", trims=sorted(done), backup=bak.name,
        generator="scripts/update_provenance_current.py")
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(doc, indent=1))
    json.loads(tmp.read_text())
    os.replace(tmp, p)
    print("  written; backup %s" % bak.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
