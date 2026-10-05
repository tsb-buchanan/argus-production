#!/usr/bin/env python3
"""Assemble data/span_efficiency.json for the report from the validated span_efficiency() output.

WHAT THIS IS. Trefftz-plane induced drag at four downstream stations and the span
efficiency derived from it, for every geometry at Condition CR. The numbers come from
span_efficiency() in scripts/case_derived_quantities.py, UNCHANGED, reading the solver's
own function-object plane VTKs. This file only SELECTS which evaluation each geometry
publishes and records the frame it was taken in.

-------------------------------------------------------------------------------
THE REFUSAL THIS FILE EXISTS FOR (D068)
-------------------------------------------------------------------------------
C_Di and e are both normalised quantities and both shift with the grid. Putting six of
them in one figure asserts they share a frame. They share a frame only if they were taken
on ONE MESH GENERATION at ONE OPERATING POINT, so that is ASSERTED here and the assertion
is the reason the script exists: it REFUSES to write the data file if the selected
evaluations do not agree on both. A figure that mixes generations would look entirely
normal, every curve plausible, and the whole comparison would be meaningless.

The same refusal covers the trim state. An evaluation taken at a geometry's pre-nudge
alpha carries a uniform freestream-w residual over the whole plane, because
span_efficiency removes w = U*sin(alpha) using the alpha it is handed. Publishing one of
those beside four trimmed evaluations is the same class of error one level down, so the
selected entry must be a TRIMMED leg for every geometry.

-------------------------------------------------------------------------------
EVERY STATION IS REPORTED, BECAUSE THERE IS NO PLATEAU
-------------------------------------------------------------------------------
Crossflow kinetic energy is conserved in an inviscid wake, so a Trefftz integral should
not decay downstream at all; what decay remains is the grid dissipating the trailing
vortex system. Quoting the last station alone would hide that, so all four stations and
the spread across them go into the data file, and the report plots the decay rather than
a point value. The spread IS the uncertainty on C_Di, not scatter to be averaged away.

e > 1 IS PHYSICALLY IMPOSSIBLE and is therefore a self-certifying gate rather than a
tolerance. Any geometry whose published entry exceeds 1 is REFUSED here, not flagged.

FRAME, travelling with the numbers (D065):
  C_Di on the HALF wing, Aref 0.620462 m2, matching each case's own forceCoeffs.
  e on AR_FULL from case_derived_quantities, DSO basis (bref 3.6576, Sref 1.24092).
  Condition CR, U_inf 34.0 m/s, every geometry at its own trimmed alpha.
"""
import json
import sys
from pathlib import Path

# THE DATA ROOT IS FOUND, NOT ASSUMED. This was parent.parent, which is the report
# root only when data/, fig/ and scripts/ are siblings. The delivery repo lifts data/
# to its own root, so the assumption held in one layout and failed silently in the
# other. data_root() walks up for the directory whose data/ carries rans_forces.json
# and RAISES if there is none, rather than returning a plausible wrong tree.
# In the working repo and in an assembled build directory it returns exactly what
# parent.parent returned, so this changes no behaviour there.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_paths import data_root  # noqa: E402
HERE = data_root()
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import case_derived_quantities as cdq  # noqa: E402

AREF_HALF = 0.620462
U_CR = 34.0
CL_TARGET = 0.428277635108
LETTER = {"SWB": "B", "SWC": "C", "SWF": "F", "SWH": "H", "SWM": "M", "SWW": "W"}

# Which raw evaluation each geometry publishes, as <GEOM>_<tag>_t<time> keys of the driver
# output. The tag IS the mesh generation and the trim state, so it is the thing asserted
# below rather than a label carried for convenience:
#   r2trim  wake-refined mesh, trimmed leg     <- the only publishable combination
#   r2      wake-refined mesh, published alpha <- untrimmed, refused
#   pub     original campaign mesh             <- different generation, refused
GENERATION = {"r2trim": "wake-refined", "r2": "wake-refined", "pub": "original"}
TRIMMED = {"r2trim": True, "r2": False, "pub": True}


def halt(msg):
    raise SystemExit("REFUSED: " + msg)


def latest_for(raw, geom, tag):
    """The highest-time evaluation of one geometry under one tag, or None."""
    keys = [k for k in raw if k.startswith(geom + "_" + tag + "_t")]
    if not keys:
        return None
    return max(keys, key=lambda k: float(k.rsplit("_t", 1)[1]))


def main():
    ap = __import__("argparse").ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", required=True, type=Path,
                    help="span_eff_r2.json written by the span_efficiency driver")
    ap.add_argument("--tag", default="r2trim",
                    help="which evaluation to publish for every geometry")
    ap.add_argument("--out", type=Path, default=HERE / "data/span_efficiency.json")
    # A SUBSET IS PERMITTED ONLY IF IT IS STATED (GEO-087). The refusal exists to stop a
    # figure captioned as the whole field being drawn over part of it, not to stop a
    # deliberate partial figure. So the bypass is available, it REQUIRES A REASON, and the
    # reason and the excluded set are written into the data file so the caption cannot
    # quietly omit them. Same shape as --allow requiring --reason elsewhere in this repo:
    # A BYPASS THAT DOES NOT RECORD ITSELF IS A BYPASS WITHOUT CONDITIONS.
    ap.add_argument("--allow-subset", action="store_true",
                    help="emit over the geometries that have a trimmed evaluation, "
                         "naming the others; requires --reason")
    ap.add_argument("--reason", default="",
                    help="why a subset is acceptable here; recorded in the data file")
    a = ap.parse_args()

    raw = json.loads(a.raw.read_text())
    if a.tag not in GENERATION:
        halt("unknown tag %r; known %s" % (a.tag, sorted(GENERATION)))
    if not TRIMMED[a.tag]:
        halt("tag %r is an UNTRIMMED leg. span_efficiency removes w = U*sin(alpha) using "
             "the alpha it is handed, so a pre-nudge evaluation leaves a uniform residual "
             "over the whole plane and is not publishable beside trimmed ones." % a.tag)

    # REPORT FIRST, EXIT AFTER (GEO-089 item 3): every geometry is named into a bucket
    # before anything is refused, so the message says WHICH are missing, not how many.
    present, missing = {}, []
    for geom in sorted(LETTER):
        k = latest_for(raw, geom, a.tag)
        if k is None:
            missing.append(geom)
        else:
            present[geom] = k
    print("geometries: %d of %d present under tag %r" % (len(present), len(LETTER), a.tag))
    for geom in sorted(present):
        print("  %s <- %s" % (geom, present[geom]))
    for geom in missing:
        print("  %s MISSING: no %s evaluation in %s" % (geom, a.tag, a.raw.name))
    assert len(present) + len(missing) == len(LETTER), "bucket sum does not close"
    if missing and not a.allow_subset:
        halt("%d of %d geometries have no %s evaluation: %s. A figure over %d of them "
             "captioned as the field is a verdict over an unstated subset (GEO-087). "
             "Pass --allow-subset with --reason to emit over the %d deliberately."
             % (len(missing), len(LETTER), a.tag, ", ".join(missing), len(present),
                len(present)))
    if missing and a.allow_subset and not a.reason.strip():
        halt("--allow-subset requires --reason. An unrecorded bypass is indistinguishable "
             "from the check not having run.")

    geoms, bad_e, bad_cl = {}, [], []
    for geom, key in sorted(present.items()):
        r = raw[key]
        planes = sorted(r["planes"], key=lambda p: p["x_m"])
        e_last = planes[-1]["span_efficiency"]
        if e_last is None or e_last > 1.0:
            bad_e.append((geom, e_last))
        cl = float(r["CL"])
        if abs(cl - CL_TARGET) * 1e4 > 1.0:
            bad_cl.append((geom, (cl - CL_TARGET) * 1e4))
        first, last = planes[0], planes[-1]
        geoms[LETTER[geom]] = {
            "case": geom + "_trim",
            "mesh_generation": GENERATION[a.tag],
            "trimmed": TRIMMED[a.tag],
            "time": r["time_req"],
            "alpha_deg": r["alpha_deg"],
            "CL": cl,
            "CL_minus_target_counts": (cl - CL_TARGET) * 1e4,
            "Cd_counts": float(r["Cd"]) * 1e4,
            "stations": [{"plane": p["plane"], "x_m": p["x_m"],
                          "CDi_counts": p["CDi_counts"],
                          "span_efficiency": p["span_efficiency"],
                          "n_faces": p["n_faces"]} for p in planes],
            "CDi_counts_first": first["CDi_counts"],
            "CDi_counts_last": last["CDi_counts"],
            "decay_first_to_last_percent":
                100.0 * (first["CDi_counts"] - last["CDi_counts"]) / first["CDi_counts"],
            "e_last": e_last,
        }

    for geom, e in bad_e:
        print("  %s e at the last station = %s" % (geom, e))
    if bad_e:
        halt("span efficiency exceeds 1 for %s. That is physically impossible and reports "
             "an unresolved wake, so it is refused rather than plotted."
             % ", ".join(g for g, _ in bad_e))
    for geom, d in bad_cl:
        print("  %s CL is %+.3f counts off target" % (geom, d))
    if bad_cl:
        halt("%s are not trimmed to within the 1.0-count gate; C_Di compared across "
             "different lift is not one comparison." % ", ".join(g for g, _ in bad_cl))

    # THE ASSERTION THIS FILE EXISTS FOR.
    gens = sorted({g["mesh_generation"] for g in geoms.values()})
    if len(gens) != 1:
        halt("the selected evaluations span %d mesh generations, %s. C_Di and e both shift "
             "with the grid, so one figure over two generations is not one comparison "
             "(D068)." % (len(gens), gens))

    xs = [tuple(round(s["x_m"], 3) for s in g["stations"]) for g in geoms.values()]
    if len(set(xs)) != 1:
        halt("the geometries were sampled at different stations: %s. A decay rate is only "
             "comparable across a common set of x." % sorted(set(xs)))

    out = {
        "what": "Trefftz-plane induced drag and span efficiency at Condition CR",
        "generated_by": "docs/report/all_geometry_2026-09-15/scripts/build_span_efficiency.py",
        "source": str(a.raw),
        "produced_by": "span_efficiency() in scripts/case_derived_quantities.py, unchanged",
        "frame": {
            "Aref_m2": AREF_HALF, "Aref_note": "half wing, matching each case's forceCoeffs",
            "AR_full": cdq.AR_FULL, "AR_note": "DSO basis, bref 3.6576 m, Sref 1.24092 m2",
            "U_inf_m_s": U_CR, "condition": "CR (M 0.10)", "CL_target": CL_TARGET,
            "mesh_generation": gens[0],
            "stations_x_m": list(xs[0]),
        },
        "gate": {
            "rule": "e = CL^2/(pi*AR*CDi) cannot exceed 1",
            "max_e_last_station": max(g["e_last"] for g in geoms.values()),
        },
        "geometries": geoms,
    }
    if missing:
        # CARRIED IN THE DATA so the figure and its caption cannot omit it.
        out["subset"] = {
            "emitted": sorted(LETTER[g] for g in present),
            "excluded": sorted(LETTER[g] for g in missing),
            "excluded_cases": sorted(missing),
            "reason": a.reason.strip(),
            "note": "This file covers a STATED SUBSET. Any caption drawn from it must "
                    "name the excluded geometries and why.",
        }
    # ---- THE EXCLUDED SET, ON ITS OWN EVALUATION, IN ITS OWN BLOCK ------------------
    # M is on the span-efficiency slide, by request, even on the published evaluation.
    # It CANNOT go in `geometries`: the assertion above refuses a mixed-generation set,
    # and rightly, because C_Di and e both shift with the grid. So the older evaluation
    # gets a SEPARATE, LABELLED block. The partition is preserved rather than widened,
    # which is the same move as an exclusion pinning its expected values instead of
    # quietly being dropped from the comparison.
    #
    # IT IS GENERATED HERE RATHER THAN HAND-ADDED because this script rewrites the whole
    # file: anything written in by hand would be silently dropped on the next rebuild, and
    # nobody would notice until a slide went out with a missing bar.
    pub = {}
    for geom in sorted(missing):
        key = "%s_pub_t4000" % geom
        r = raw.get(key)
        if not r or not r.get("planes"):
            continue
        planes = r["planes"]
        first, last = planes[0], planes[-1]
        pub[LETTER[geom]] = {
            "case": geom,
            "mesh_generation": "published",
            "trimmed": True,
            "time": 4000,
            "alpha_deg": r.get("alpha_deg"),
            "CL": r.get("CL"),
            "CL_minus_target_counts": (float(r["CL"]) - CL_TARGET) * 1e4,
            "Cd_counts": float(r["Cd"]) * 1e4,
            "stations": [{"plane": q["plane"], "x_m": q["x_m"],
                          "CDi_counts": q["CDi_counts"],
                          "span_efficiency": q["span_efficiency"],
                          "n_faces": q["n_faces"]} for q in planes],
            "CDi_counts_first": first["CDi_counts"],
            "CDi_counts_last": last["CDi_counts"],
            "decay_first_to_last_percent":
                100.0 * (first["CDi_counts"] - last["CDi_counts"]) / first["CDi_counts"],
            "e_last": last["span_efficiency"],
        }
    if pub:
        out["published_evaluation"] = {
            "what": "the EXCLUDED geometries on their own earlier evaluation, for "
                    "completeness only",
            "mesh_generation": "published",
            "MUST_NOT": "be plotted on the same axis as `geometries` without a visual "
                        "distinction and a label naming the evaluation. C_Di and e both "
                        "shift with the grid, so the two sets are not one comparison "
                        "(D068). Draw it hatched, or in a separate panel.",
            "produced_by": "the SAME span_efficiency() call and the SAME four planes as "
                           "`geometries`; only the input evaluation differs",
            "geometries": pub,
        }
        print("\npublished_evaluation block: %s" % ", ".join(sorted(pub)))
        for L, g in sorted(pub.items()):
            print("  %s  CDi %7.3f -> %7.3f ct  (decay %5.2f%%)  e_last %.4f"
                  % (L, g["CDi_counts_first"], g["CDi_counts_last"],
                     g["decay_first_to_last_percent"], g["e_last"]))

    a.out.write_text(json.dumps(out, indent=2) + "\n")
    print("\nwrote %s" % a.out)
    print("  mesh generation: %s   stations: %s" % (gens[0], list(xs[0])))
    for L in sorted(geoms, key=lambda L: geoms[L]["e_last"]):
        g = geoms[L]
        print("  %s  CDi %7.3f -> %7.3f ct  (decay %5.2f%%)   e_last %.4f"
              % (L, g["CDi_counts_first"], g["CDi_counts_last"],
                 g["decay_first_to_last_percent"], g["e_last"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
