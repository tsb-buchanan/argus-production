#!/usr/bin/env python3
"""Build Q_PINF_TABLE.csv for a sections directory by READING that directory's own logs.

WHY THIS EXISTS. The published data/sections_nflag/Q_PINF_TABLE.csv carries q, p_inf, the
fit residuals, alpha, dragDir and the shear sign for every case. Those are all values that
section_slices_from_surface.py PRINTS and does not write, so the table was assembled by
hand once. A hand-assembled table beside a regenerated set of CSVs is the half-updated
artefact this project keeps paying for, so the table is now derived from the logs the
extraction itself wrote (GEO-080: the claim about a file's contents comes from reading the
file, never from the hand that caused it).

VALIDATED AGAINST THE PUBLISHED TABLE. Run with --check it re-derives every row of
data/sections_nflag/Q_PINF_TABLE.csv from that directory's logs and exits non-zero on any
field that disagrees. That is a known-answer null (D058 item 2), not a self-comparison.

COVERAGE IS ASSERTED, NOT ASSUMED (GEO-087). Every *.log in the directory is enumerated
into exactly one of PARSED or UNPARSEABLE, and the two are asserted to sum to the file
count. A log the parser cannot read is a failure, never a silent omission.
"""
import argparse
import csv
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent

# The condition label is a property of the case name, and it is the one thing in this table
# that is not in the log. CR is the M 0.10 wall-resolved set, CMP the M 0.78 wall-modelled
# one; a case matching neither prefix is a failure rather than a default.
PREFIX_CONDITION = {"SW": "CR (M 0.10)", "CMP": "early cruise (M 0.78)", "LC": "late cruise (M 0.78)"}

FIELDS = ["case", "condition", "surface_time", "n_faces", "q_from_file", "p_inf_from_file",
          "fit_resid_max_over_q", "fit_resid_p99_over_q", "alpha_deg", "dragDir",
          "shear_sign_multiplier"]

# MULTILINE, because these anchors are matched against the WHOLE log text. Without it
# `^` binds only to the start of the string and every faces line reads as absent, which
# is exactly how the first run of this parser scored 12 of 12 logs UNPARSEABLE.
RE_READ = re.compile(r"^reading (.+)$", re.M)
RE_FACES = re.compile(r"^\s+(\d+) faces, fields:", re.M)
RE_Q = re.compile(r"DERIVED FROM THIS FILE: q = (\S+), p_inf = (\S+) "
                  r"\(fit residual: max (\S+) of q, p99 (\S+)\)")
RE_DRAG = re.compile(r"freestream dragDir \[(.+?)\] \(alpha (\S+) deg\)")
RE_SHEAR = re.compile(r"multiplier (-?\d+) \(attached is positive\)")


def condition_for(case):
    for p, c in sorted(PREFIX_CONDITION.items(), key=lambda kv: -len(kv[0])):
        if case.startswith(p):
            return c
    raise SystemExit("FATAL: %s matches no known condition prefix %s"
                     % (case, sorted(PREFIX_CONDITION)))


def fmt_g(x, sig=6):
    """Match the published table's %g-style rendering of q and p_inf."""
    return "%g" % float("%.*g" % (sig, float(x)))


def parse_log(path, surface_time):
    """Return one table row, or raise KeyError naming the field that was absent.

    Nothing here is defaulted. A log missing any of the five lines is UNPARSEABLE and the
    caller counts it as such; a row assembled from four of five would be a partial record
    wearing a complete record's shape.
    """
    txt = path.read_text(errors="replace")
    m_f = RE_FACES.search(txt)
    m_q = RE_Q.search(txt)
    m_d = RE_DRAG.search(txt)
    m_s = RE_SHEAR.search(txt)
    for name, m in (("n_faces", m_f), ("q", m_q), ("dragDir", m_d), ("shear sign", m_s)):
        if m is None:
            raise KeyError(name)
    case = path.name[:-len(".log")]
    # The surface time is in the VTK path when the extraction read it out of a case
    # tree, and absent when it read a hand-staged copy. Derive it where it exists and
    # fall back to what the caller supplied, rather than trusting the caller over the
    # file.
    m_r = RE_READ.search(txt)
    if m_r:
        parts = m_r.group(1).split("/")
        if len(parts) >= 2 and parts[-2].isdigit():
            surface_time = parts[-2]
    dx, dy, dz = [float(v) for v in m_d.group(1).split()]
    alpha = math.degrees(math.atan2(dz, dx))
    if abs(alpha - float(m_d.group(2))) > 5e-6:
        raise SystemExit("FATAL: %s alpha %.9f disagrees with the logged %s"
                         % (path, alpha, m_d.group(2)))
    return {
        "case": case,
        "condition": condition_for(case),
        "surface_time": surface_time,
        "n_faces": m_f.group(1),
        "q_from_file": fmt_g(m_q.group(1)),
        "p_inf_from_file": fmt_g(m_q.group(2)),
        "fit_resid_max_over_q": "%.2e" % float(m_q.group(3)),
        "fit_resid_p99_over_q": "%.2e" % float(m_q.group(4)),
        "alpha_deg": "%.6f" % alpha,
        "dragDir": "%.9f %.9f %.9f" % (dx, dy, dz),
        "shear_sign_multiplier": m_s.group(1),
    }


def build(sect_dir, times):
    logs = sorted(sect_dir.glob("*.log"))
    parsed, unparseable = [], []
    for lg in logs:
        case = lg.name[:-len(".log")]
        try:
            parsed.append(parse_log(lg, times.get(case, "")))
        except KeyError as e:
            unparseable.append((lg.name, str(e)))
    # GEO-087: enumerated must equal adjudicated, asserted rather than observed.
    assert len(parsed) + len(unparseable) == len(logs), "bucket sum does not close"
    return parsed, unparseable, len(logs)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sections", required=True, type=Path,
                    help="sections directory holding <case>.log files")
    ap.add_argument("--time", action="append", default=[], metavar="CASE=TIME",
                    help="surface time for a case; the log does not carry it")
    ap.add_argument("--check", action="store_true",
                    help="re-derive and diff against the directory's existing table")
    a = ap.parse_args()

    times = dict(t.split("=", 1) for t in a.time)
    rows, bad, n = build(a.sections, times)
    # REPORT FIRST, EXIT AFTER (GEO-089 item 3).
    for name, field in bad:
        print("UNPARSEABLE %s: no %s line" % (name, field), file=sys.stderr)
    print("logs: %d total = %d PARSED + %d UNPARSEABLE" % (n, len(rows), len(bad)),
          file=sys.stderr)
    if bad:
        return 1

    out = a.sections / "Q_PINF_TABLE.csv"
    if a.check:
        if not out.exists():
            print("FATAL: --check needs an existing %s" % out, file=sys.stderr)
            return 1
        have = {r["case"]: r for r in csv.DictReader(out.open())}
        diffs = 0
        for r in rows:
            old = have.get(r["case"])
            if old is None:
                print("MISSING from published table: %s" % r["case"], file=sys.stderr)
                diffs += 1
                continue
            for k in FIELDS:
                if k == "surface_time" and not r[k]:
                    continue  # not supplied on the command line; nothing to compare
                if k == "dragDir":
                    # COMPARE THE VECTOR, NOT ITS RENDERING. The published table was
                    # transcribed from the log text, which prints numpy's "0." and drops a
                    # trailing zero; this writes a fixed 9-decimal field. A string compare
                    # would report eleven disagreements where no component differs at all,
                    # which is a formatting difference wearing a data difference's clothes.
                    if any(abs(float(u) - float(v)) > 5e-10 for u, v in
                           zip(old[k].split(), r[k].split())):
                        print("DIFF %s.%s: published %r derived %r"
                              % (r["case"], k, old[k], r[k]), file=sys.stderr)
                        diffs += 1
                    continue
                if str(old[k]).strip() != str(r[k]).strip():
                    print("DIFF %s.%s: published %r derived %r"
                          % (r["case"], k, old[k], r[k]), file=sys.stderr)
                    diffs += 1
        print("check: %d rows, %d field disagreements" % (len(rows), diffs), file=sys.stderr)
        return 1 if diffs else 0

    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote %s (%d rows)" % (out, len(rows)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
