#!/usr/bin/env python3
"""Harvest every converged 3D RANS result on HPC12 into one JSON blob.

Runs ON HPC12 (python 3.6 compatible). Emits JSON on stdout.

Design notes, per project standing rules:

1. GEO-080 (nothing vouches for a step it did not observe). Every field is
   DERIVED by reading the case: alpha comes from the forceCoeffs header's
   dragDir (what the solver actually used), Aref/magUInf/lRef from the same
   header, endTime from system/controlDict. Nothing is authored here.

2. D058 (a null needs an independently known answer). The window mean for a
   LEGACY marker has to be recomputed. The SAME code path is run on every
   NEW-format marker, whose window mean is independently known because the
   gate wrote it into the marker. Any disagreement is reported as a residual,
   so the legacy recomputation is validated rather than assumed.

3. GEO-087 / GEO-089 (a gate must partition its input set). Every case
   directory lands in exactly one named bucket and the buckets are asserted
   to sum to the enumerated total.

4. Two marker formats exist and a reader that handles only one silently drops
   converged cases (this already happened with CMPC_a1p20). Both are parsed
   and the format is recorded per case.
"""

import json
import math
import os
import re
import sys

# Case root on HPC12, with $USER expanded at run time. Set ARGUS_SOLVE to read another root.
SOLVE = os.environ.get("ARGUS_SOLVE") or os.path.expandvars("/home/scratch/$USER/argus/solve")
ARCHIVE = os.path.expanduser("~/argus_archive")
WINDOW = 200

FAMILY = {"B": "baseline_EET_AR12", "C": "continuous_TE_camber",
          "F": "F", "H": "H", "M": "M", "W": "W"}
CONDITION = {"SW": "condition_CR", "CMP": "early_cruise", "LC": "late_cruise"}


def case_prefix(case):
    """Split a case name into (condition prefix, family letter)."""
    for pre in ("CMP", "SW", "LC"):
        if case.startswith(pre):
            return pre, case[len(pre)]
    return None, None


# ---------------------------------------------------------------- markers

RE_NEW_GATE = re.compile(
    r"CONVERGENCE GATE over\s+(\d+)\s+samples:\s*MEAN drift\s+"
    r"Cd\s+([-\d.eE+]+)\s*ct.*?"
    r"Cl\s+([-\d.eE+]+)\s*ct.*?"
    r"Cl SPAN\s+([-\d.eE+]+)\s*ct.*?->\s*(\S+)", re.I)
RE_NEW_REPORT = re.compile(
    r"REPORT THESE:\s*Cd\s+([-\d.eE+]+)\s+Cl\s+([-\d.eE+]+)", re.I)
RE_NEW_FINAL = re.compile(r"finalised\s+(\S+)\s+on\s+(\S+?),\s*leg\s*(\d+)", re.I)
RE_LEGACY = re.compile(
    r"case\s+(\S+):\s*reached\s+([\d.]+)\s+of\s+endTime\s+([\d.]+);\s*"
    r"drift\s+Cd\s+([-\d.eE+]+)\s*ct.*?"
    r"Cl\s+([-\d.eE+]+)\s*ct.*?"
    r"Cl span\s+([-\d.eE+]+)\s*ct.*?->\s*(\S+)", re.I)


def parse_marker(path):
    """Return a dict of what the marker itself states. Never guesses."""
    with open(path) as fh:
        txt = fh.read()
    out = {"marker_path": path, "marker_raw": txt.strip()}

    m = RE_NEW_GATE.search(txt)
    r = RE_NEW_REPORT.search(txt)
    if m and r:
        out["marker_format"] = "window_means"
        out["gate_samples"] = int(m.group(1))
        out["drift_Cd_ct"] = float(m.group(2))
        out["drift_Cl_ct"] = float(m.group(3))
        out["Cl_span_ct"] = float(m.group(4))
        out["gate_verdict"] = m.group(5)
        out["marker_Cd"] = float(r.group(1))
        out["marker_Cl"] = float(r.group(2))
        f = RE_NEW_FINAL.search(txt)
        if f:
            out["finalised_utc"] = f.group(1)
            out["finalised_node"] = f.group(2)
            out["leg"] = int(f.group(3))
        out["marker_iterations"] = None
        out["marker_endTime"] = None
        return out

    m = RE_LEGACY.search(txt)
    if m:
        out["marker_format"] = "legacy_single_line"
        out["gate_samples"] = None
        out["marker_iterations"] = float(m.group(2))
        out["marker_endTime"] = float(m.group(3))
        out["drift_Cd_ct"] = float(m.group(4))
        out["drift_Cl_ct"] = float(m.group(5))
        out["Cl_span_ct"] = float(m.group(6))
        out["gate_verdict"] = m.group(7)
        out["marker_Cd"] = None
        out["marker_Cl"] = None
        return out

    out["marker_format"] = "UNPARSEABLE"
    return out


# ------------------------------------------------------------ forceCoeffs

RE_VEC = re.compile(r"\(([-\d.eE+\s]+)\)")


def parse_fc_header(path):
    hdr = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if ":" not in line:
                continue
            key, val = line[1:].split(":", 1)
            key = key.strip()
            val = val.strip()
            v = RE_VEC.search(val)
            if v:
                hdr[key] = [float(x) for x in v.group(1).split()]
            else:
                try:
                    hdr[key] = float(val)
                except ValueError:
                    hdr[key] = val
    return hdr


def read_fc_series(case_dir):
    """Merge every forceCoeffs1 time directory into one Time->row series.

    Restart handling: directories are applied in ASCENDING start time, so a
    later restart OVERWRITES the overlapping tail of an earlier run. SWM_a1p20
    needs this (its 0/ dir runs to 3430 while the 3000/ restart covers
    3000-4000); taking the union without overwriting would blend two different
    solution paths over 3000-3430.
    """
    root = os.path.join(case_dir, "postProcessing", "forceCoeffs1")
    if not os.path.isdir(root):
        return None, None, [], None
    dirs = []
    for name in os.listdir(root):
        f = os.path.join(root, name, "forceCoeffs.dat")
        if os.path.isfile(f):
            try:
                dirs.append((float(name), f))
            except ValueError:
                pass
    if not dirs:
        return None, None, [], None
    dirs.sort()

    series = {}
    headers = []
    for _, f in dirs:
        headers.append((f, parse_fc_header(f)))
        with open(f) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                if len(p) < 4:
                    continue
                try:
                    t = float(p[0])
                    series[t] = (float(p[1]), float(p[2]), float(p[3]))
                except ValueError:
                    continue
    # header of record = the last restart leg, i.e. what the solver last used
    hdr_path, hdr = headers[-1]
    # every leg must agree on the frame; a disagreement is reported, not hidden
    mismatches = []
    for f, h in headers[:-1]:
        for k in ("dragDir", "Aref", "magUInf", "lRef"):
            if h.get(k) != hdr.get(k):
                mismatches.append({"file": f, "field": k,
                                   "value": h.get(k), "header_of_record": hdr.get(k)})
    return hdr, hdr_path, sorted(series.items()), mismatches


def window_mean(series, n=WINDOW):
    """Mean of the last n samples. The single recomputation route, used for
    legacy markers and run on new-format markers as its own known-answer null."""
    tail = series[-n:]
    cd = sum(v[1] for _, v in tail) / float(len(tail))
    cl = sum(v[2] for _, v in tail) / float(len(tail))
    return cd, cl, len(tail), tail[0][0], tail[-1][0]


def read_endtime(case_dir):
    cd = os.path.join(case_dir, "system", "controlDict")
    if not os.path.isfile(cd):
        return None, None
    with open(cd) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("endTime"):
                tok = s.rstrip(";").split()
                if len(tok) >= 2:
                    try:
                        return float(tok[1]), cd
                    except ValueError:
                        return None, cd
    return None, cd


def case_status(case_dir):
    """Classify a case with no .converged marker. Derived from the tree."""
    has_proc = any(n.startswith("processor") for n in os.listdir(case_dir))
    logdir = os.path.join(case_dir, "log")
    nlog = len(os.listdir(logdir)) if os.path.isdir(logdir) else 0
    has_zero = os.path.isdir(os.path.join(case_dir, "0"))
    if has_proc and nlog > 0:
        return "running", {"processor_dirs": True, "log_files": nlog, "zero_dir": has_zero}
    return "staged", {"processor_dirs": has_proc, "log_files": nlog, "zero_dir": has_zero}


# ----------------------------------------------------------------- main

def main():
    solve = sorted(d for d in os.listdir(SOLVE)
                   if os.path.isdir(os.path.join(SOLVE, d)) and not d.startswith("_"))
    arch = sorted(d for d in os.listdir(ARCHIVE)
                  if os.path.isdir(os.path.join(ARCHIVE, d)) and not d.startswith("_"))
    enumerated = sorted(set(solve) | set(arch))

    cases = {}
    buckets = {"converged": [], "running": [], "staged": [], "unparseable_marker": []}
    validations = []

    for case in enumerated:
        sdir = os.path.join(SOLVE, case)
        adir = os.path.join(ARCHIVE, case)
        smark = os.path.join(sdir, ".converged")
        amark = os.path.join(adir, ".converged")

        if os.path.isfile(smark):
            marker, case_dir, src = smark, sdir, "scratch_solve"
        elif os.path.isfile(amark):
            marker, case_dir, src = amark, adir, "argus_archive"
        else:
            base = sdir if os.path.isdir(sdir) else adir
            st, ev = case_status(base)
            buckets[st].append(case)
            pre, fam = case_prefix(case)
            cases[case] = {
                "case": case, "family": FAMILY.get(fam, fam),
                "family_letter": fam,
                "condition": CONDITION.get(pre),
                "status": st, "evidence": ev, "case_dir": base,
                "converged": False,
                "alpha_deg": "NOT FOUND", "Cd_counts": "NOT FOUND",
                "Cl": "NOT FOUND", "L_over_D": "NOT FOUND",
            }
            continue

        m = parse_marker(marker)
        if m["marker_format"] == "UNPARSEABLE":
            buckets["unparseable_marker"].append(case)
        else:
            buckets["converged"].append(case)

        hdr, hdr_path, series, mism = read_fc_series(case_dir)
        # fall back to the other copy if this one was purged of postProcessing
        alt = adir if case_dir == sdir else sdir
        fc_dir_used = case_dir
        if not series and os.path.isdir(alt):
            hdr, hdr_path, series, mism = read_fc_series(alt)
            fc_dir_used = alt

        pre, fam = case_prefix(case)
        rec = {
            "case": case,
            "family": FAMILY.get(fam, fam),
            "family_letter": fam,
            "condition": CONDITION.get(pre),
            "is_trim_case": case.endswith("_trim"),
            "status": "converged",
            "converged": m.get("gate_verdict", "").upper() == "CONVERGED",
            "gate_verdict": m.get("gate_verdict"),
            "marker_format": m["marker_format"],
            "marker_source": src,
            "drift_Cd_ct": m.get("drift_Cd_ct"),
            "drift_Cl_ct": m.get("drift_Cl_ct"),
            "Cl_span_ct": m.get("Cl_span_ct"),
            "gate_samples": m.get("gate_samples"),
            "sources": {"marker": marker},
            "marker_text": m["marker_raw"],
        }
        for k in ("finalised_utc", "finalised_node", "leg"):
            if k in m:
                rec[k] = m[k]

        if not series:
            rec.update({
                "alpha_deg": "NOT FOUND", "Cd_counts": "NOT FOUND",
                "Cl": "NOT FOUND", "L_over_D": "NOT FOUND",
                "iterations": "NOT FOUND", "endTime": "NOT FOUND",
                "note": "no forceCoeffs1 output found in either scratch or archive",
            })
            cases[case] = rec
            continue

        drag = hdr.get("dragDir")
        alpha = math.degrees(math.atan2(drag[2], drag[0])) if drag else None
        cd_w, cl_w, nw, t0, t1 = window_mean(series)

        if m["marker_format"] == "window_means":
            cd, cl = m["marker_Cd"], m["marker_Cl"]
            cd_prov = marker + " (gate window mean, REPORT THESE)"
            cl_prov = cd_prov
            # known-answer null: the same recomputation route, checked against
            # the value the gate itself wrote.
            validations.append({
                "case": case,
                "marker_Cd": cd, "recomputed_Cd": cd_w,
                "resid_Cd_ct": (cd_w - cd) * 1e4,
                "marker_Cl": cl, "recomputed_Cl": cl_w,
                "resid_Cl_ct": (cl_w - cl) * 1e4,
                "n_samples_used": nw,
            })
        else:
            cd, cl = cd_w, cl_w
            cd_prov = (os.path.join(fc_dir_used, "postProcessing/forceCoeffs1/*/forceCoeffs.dat")
                       + " (mean of last %d samples, t %g to %g; legacy marker carries no window mean)"
                       % (nw, t0, t1))
            cl_prov = cd_prov

        endtime, cdict = read_endtime(case_dir)
        if endtime is None and fc_dir_used != case_dir:
            endtime, cdict = read_endtime(fc_dir_used)

        rec.update({
            "alpha_deg": alpha,
            "dragDir": drag,
            "Cd": cd,
            "Cd_counts": cd * 1e4,
            "Cl": cl,
            "L_over_D": cl / cd if cd else None,
            "Aref_m2": hdr.get("Aref"),
            "lRef_m": hdr.get("lRef"),
            "magUInf_m_s": hdr.get("magUInf"),
            "window_samples_used": nw,
            "window_t_first": t0,
            "window_t_last": t1,
            "iterations": series[-1][0],
            "n_force_samples_total": len(series),
            "endTime": endtime,
            "endTime_marker_stated": m.get("marker_endTime"),
            "iterations_marker_stated": m.get("marker_iterations"),
            "restart_legs": len(set(x for x in os.listdir(
                os.path.join(fc_dir_used, "postProcessing", "forceCoeffs1")))),
            "header_mismatches_across_legs": mism,
        })
        rec["sources"].update({
            "Cd": cd_prov, "Cl": cl_prov,
            "alpha_deg": hdr_path + " (# dragBir/dragDir header vector)",
            "Aref_m2": hdr_path, "lRef_m": hdr_path, "magUInf_m_s": hdr_path,
            "iterations": hdr_path,
            "endTime": cdict if cdict else "NOT FOUND",
            "drift_Cd_ct": marker, "drift_Cl_ct": marker, "Cl_span_ct": marker,
        })
        rec["sources"]["alpha_deg"] = hdr_path + " (dragDir header vector)"
        cases[case] = rec

    total = len(enumerated)
    part_sum = sum(len(v) for v in buckets.values())
    out = {
        "enumerated_case_dirs": total,
        "partition": {k: sorted(v) for k, v in buckets.items()},
        "partition_counts": dict((k, len(v)) for k, v in buckets.items()),
        "partition_sums_to_total": part_sum == total,
        "partition_sum": part_sum,
        "cases": cases,
        "window_mean_route_validation": validations,
        "roots": {"scratch_solve": SOLVE, "argus_archive": ARCHIVE},
    }
    json.dump(out, sys.stdout, indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
