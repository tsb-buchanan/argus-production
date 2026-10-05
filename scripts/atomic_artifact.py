#!/usr/bin/env python3
"""Write an artifact ONLY if it validates. Nothing invalid reaches a final filename.

GEO-081, a project decision. THE DEFECT IS THE SURVIVAL MECHANISM, NOT THE MISSING
FIELD. A collector that writes to the final path and validates afterwards presents
a failure as a CRASH while leaving the invalid artifact on disk under its real
name. Anything downstream reads the file, not the exit status, and cannot tell it
was rejected. That is how 22 of 22 run cards stayed invalid from d1abb50 onward.

THE FIX IS ORDERING, NOT DILIGENCE:
    VALIDATE BEFORE WRITE, OR WRITE TO A TEMPORARY PATH AND PROMOTE ONLY ON PASS.

write_if_valid() does the second, which also covers validators that need a real
file on disk. The temporary lives in the DESTINATION DIRECTORY so the promoting
os.replace is atomic (same filesystem), and it is removed on failure.

SECOND RULE, same commit: THE RUN-CARD SCHEMA IS PARTITIONED INTO DERIVED AND
RECORDED-AT-RUN-TIME FIELDS, and the partition lives in the schema so it covers
fields nobody has added yet. preserve_run_time_fields() reads the class off the
schema, so patching one more field by hand is never the answer.
    DERIVED             regenerate every time from the case: solver, turbulence
                        model, wall treatment, mesh stats, residual targets, forces
    RECORDED-AT-RUN-TIME immutable on regeneration: git_commit, date, machine,
                        wall-clock, core-hours, OpenFOAM build
The `date` field already had this protection with a comment explaining it, and the
comment did not reach the field beside it. Patching git_commit alone would have
left machine, wall_clock, core_hours and the build waiting to repeat it.
"""
import json
import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCHEMA = REPO / "scripts" / "run_card.schema.json"

# Schema keyword marking the class of a field. Custom `x-` keywords are ignored by
# JSON Schema validators, so this is inert to validation and readable by us.
CLASS_KEY = "x-argus-field-class"
DERIVED = "derived"
RUN_TIME = "recorded-at-run-time"


class ArtifactRejected(RuntimeError):
    """The candidate artifact failed validation, so it was never written."""


def write_if_valid(path, content, validator, encoding="utf-8"):
    """Write `content` to `path` only if `validator(tmp_path)` passes.

    validator: callable taking the temporary file's Path. It must RAISE or return
    a falsy-but-explanatory value on failure. Anything it returns that is a
    non-empty string or list is treated as an error report.

    On failure the temporary is removed and ArtifactRejected is raised, so the
    previous contents of `path` survive untouched. A rejected artifact never
    exists under its real name, not even briefly.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix="." + path.name + ".",
                               suffix=".candidate")
    tmp = Path(tmp)
    try:
        with os.fdopen(fd, "w", encoding=encoding) as fh:
            fh.write(content)
        problems = validator(tmp)
        if problems:
            raise ArtifactRejected(
                "%s rejected by its validator and NOT written:\n  %s"
                % (path, problems if isinstance(problems, str) else
                   "\n  ".join(str(p) for p in problems)))
        os.replace(str(tmp), str(path))          # atomic within one filesystem
        tmp = None
    finally:
        if tmp is not None and tmp.exists():
            tmp.unlink()
    return path


def run_card_validator(tmp_path):
    """Validator for run cards: shells out to the project's own validator."""
    import subprocess
    r = subprocess.run([sys.executable, str(REPO / "scripts" / "validate_run_card.py"),
                        str(tmp_path)], capture_output=True, text=True, cwd=str(REPO))
    return "" if r.returncode == 0 else (r.stdout + r.stderr).strip()


def _walk_schema(node, path=()):
    """Yield (dotted_path, subschema) for every named property in the schema."""
    if not isinstance(node, dict):
        return
    for name, sub in (node.get("properties") or {}).items():
        yield ".".join(path + (name,)), sub
        yield from _walk_schema(sub, path + (name,))


def run_time_fields(schema=None):
    """Dotted paths of every field the SCHEMA marks as recorded-at-run-time.

    Read from the schema, never hard-coded here: a field added later inherits the
    protection without anyone remembering this module exists.
    """
    schema = schema if schema is not None else json.loads(SCHEMA.read_text())
    return [p for p, sub in _walk_schema(schema)
            if isinstance(sub, dict) and sub.get(CLASS_KEY) == RUN_TIME]


def _get(d, dotted):
    for k in dotted.split("."):
        if not isinstance(d, dict) or k not in d:
            return None, False
        d = d[k]
    return d, True


def _set(d, dotted, value):
    keys = dotted.split(".")
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


def preserve_run_time_fields(card, prior, schema=None):
    """Copy every RECORDED-AT-RUN-TIME field from `prior` into `card`.

    Regeneration may recompute anything DERIVED, and must not touch anything
    recorded when the run happened. Returns the list of preserved dotted paths so
    the caller can report rather than assume.
    """
    if not prior:
        return []
    kept = []
    for dotted in run_time_fields(schema):
        val, present = _get(prior, dotted)
        if present:
            _set(card, dotted, val)
            kept.append(dotted)
    return kept


def self_test():
    import shutil
    tmpdir = Path(tempfile.mkdtemp(prefix="argus_atomic_"))
    ok = True
    try:
        target = tmpdir / "artifact.json"
        target.write_text('{"original": true}\n')

        # 1. A rejected artifact must NOT reach the final name, and must not
        #    disturb what is already there.
        try:
            write_if_valid(target, '{"bad": true}\n', lambda p: "schema says no")
            print("  rejected artifact written   : *** FAIL, it was written ***")
            ok = False
        except ArtifactRejected:
            survived = json.loads(target.read_text()).get("original") is True
            leftovers = [p.name for p in tmpdir.iterdir() if p.name != "artifact.json"]
            print("  rejected artifact NOT written: %s" % ("PASS" if survived else "FAIL"))
            print("  prior contents survived      : %s" % ("PASS" if survived else "FAIL"))
            print("  no candidate left behind     : %s"
                  % ("PASS" if not leftovers else "FAIL %s" % leftovers))
            ok &= survived and not leftovers

        # 2. A valid artifact is promoted.
        write_if_valid(target, '{"original": false, "new": true}\n', lambda p: "")
        promoted = json.loads(target.read_text()).get("new") is True
        print("  valid artifact promoted      : %s" % ("PASS" if promoted else "FAIL"))
        ok &= promoted

        # 3. The schema partition is READ FROM THE SCHEMA, and a field added to the
        #    schema later is protected without touching this module.
        fields = run_time_fields()
        want = {"provenance.git_commit", "provenance.date", "cost.machine",
                "cost.wall_clock_s", "cost.core_hours"}
        got = set(fields)
        print("  run-time fields from schema  : %s (%d found)"
              % ("PASS" if want <= got else "FAIL missing %s" % sorted(want - got),
                 len(got)))
        ok &= want <= got

        card = {"provenance": {"git_commit": "NEW", "openfoam_version": "openfoam-org-7"},
                "cost": {"machine": "NEW", "core_hours": 9.9}}
        prior = {"provenance": {"git_commit": "OLD", "openfoam_version": "openfoam-org-7"},
                 "cost": {"machine": "OLD", "core_hours": 1.1}}
        kept = preserve_run_time_fields(card, prior)
        held = (card["provenance"]["git_commit"] == "OLD"
                and card["cost"]["machine"] == "OLD"
                and card["cost"]["core_hours"] == 1.1)
        # and a DERIVED field must NOT be preserved
        derived_free = card["provenance"]["openfoam_version"] == "openfoam-org-7"
        print("  run-time fields immutable    : %s (%d preserved)"
              % ("PASS" if held else "FAIL", len(kept)))
        print("  derived fields still derived : %s" % ("PASS" if derived_free else "FAIL"))
        ok &= held and derived_free
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    print("\nSELF TEST: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(self_test())
