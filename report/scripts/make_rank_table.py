#!/usr/bin/env python3
"""Emit the selection-test table (tab:qv:rank) and its prose values, read from the data.

WHY GENERATE RATHER THAN TYPE. This table and the three paragraphs around it are the
report's answer to "does the low-order search order the candidates correctly", and the
answer is a FUNCTION OF THE RANS DRAG INCREMENTS. Those increments changed: F, W and H
moved from +5.733, +5.734 and +6.580 counts to negative values, which reorders the RANS
column and therefore every concordance count, every Kendall tau and every sentence that
quotes them.

THE CONCLUSION REVERSED, WHICH IS WHY THIS MUST NOT BE HAND-EDITED. On the previous
increments the near-field metric reproduced the RANS order exactly (tau +1.00, 0 pairs
swapped) and the far-field metric did not (tau 0.00, 5 swapped). On the corrected
increments it is the other way round. A table that size, updated by hand against a
reversed conclusion, is where a stale cell survives and quietly contradicts the paragraph
above it.

RANKS AND STATISTICS ARE RECOMPUTED, NEVER CARRIED. Kendall tau, Spearman rho, the
concordant and discordant counts and the swapped-pair list all come from
fig/fig_vlmerr_values.json, which fig/make_vlmerr.py writes from the force data. Run
make_vlmerr.py first; this reads its output.

TIES ARE REPORTED, NOT ASSUMED. The old caption said "RANS ranks 3 and 4, F and W, are
separated by 0.001 ct and are a tie". On the corrected numbers that separation is 0.094 ct
and they are not tied. The tie threshold is stated here and applied, rather than described
in a caption that nothing checks.
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
LETTERS = ["C", "M", "F", "W", "H"]
TIE_CT = 0.01  # increments closer than this are reported as a tie, and the pair is named

# Early-cruise increments, which tab:perf:vlm carries beside the Condition CR ones. They are
# NOT recomputed here: early cruise is a different operating point and nothing in this
# session re-ran it, so these are the published values, carried unchanged and labelled.
EARLY_CRUISE = {"C": 13.105, "M": 14.646, "F": 55.619, "W": 50.119, "H": 54.246}



def _fmt(v, nd):
    return "--" if v is None else ("%%.%df" % nd) % v


def _vlm_span_efficiency():
    """VLM near- and far-field span efficiency per geometry at Condition CR.

    Keyed by the report letter through the same geometry mapping the increments use, so a
    deck that maps to no meshed geometry contributes nothing rather than being guessed at.
    """
    v = json.loads((HERE / "data/vlm.json").read_text())
    name_to_letter = {"cte_i002_c04": "C", "mcv2_i002_c01": "M", "cfft_b02_c01": "F",
                      "cffw_b01_c01": "W", "chc_g02_c06": "H",
                      "baseline_wing_only_refined": "B"}
    E, Ew = {}, {}
    for r in v.get("runs", []):
        L = name_to_letter.get(r.get("geometry_name"))
        pol = r.get("polar_on_DSO_basis") or {}
        if L is None or not pol:
            continue
        if abs(float(pol.get("Mach", 0)) - 0.10) > 1e-6:
            continue          # Condition CR only; never mix operating points (D076)
        E[L], Ew[L] = pol.get("E"), pol.get("Ew")
    return E, Ew


def _rans_span_efficiency():
    """RANS span efficiency at the last Trefftz station, for every geometry that has one.

    BOTH BLOCKS ARE READ. A geometry the data file carries under `published_evaluation`
    rather than under `geometries` still HAS a measured value, and printing "--" beside the
    other five says the quantity was never obtained rather than that it came from the
    earlier evaluation. Project decision, 2026-09-23: the table carries the result; the mesh
    generation is recorded in data/span_efficiency.json and fig/fig_spaneff_values.csv.

    THE PRECONDITION IS ASSERTED, NOT ASSUMED (GEO-092). A published entry is taken only if
    it is a TRIMMED leg; an untrimmed one carries a uniform freestream-w residual over the
    plane and is not the same quantity. Anything failing that is left out and named, which
    is what "--" is actually for.
    """
    f = HERE / "data/span_efficiency.json"
    if not f.exists():
        return {}
    d = json.loads(f.read_text())
    out = {L: e["e_last"] for L, e in d["geometries"].items()}
    for L, e in d.get("published_evaluation", {}).get("geometries", {}).items():
        if L in out:
            continue
        if not e.get("trimmed"):
            print("%% %s published span efficiency skipped: not a trimmed leg" % L,
                  file=sys.stderr)
            continue
        out[L] = e["e_last"]
    return out


def ranks(values):
    """1 = smallest. Returns {id: rank}."""
    order = sorted(values, key=lambda k: values[k])
    return {k: i + 1 for i, k in enumerate(order)}


def main():
    src = HERE / "fig/fig_vlmerr_values.json"
    if not src.exists():
        sys.exit("HALT: no %s. Run fig/make_vlmerr.py first." % src)
    d = json.loads(src.read_text())
    e = d["errors_vs_mapped_B"]
    rk = d["ranking"]

    missing = [L for L in LETTERS if L not in e]
    if missing:
        sys.exit("HALT: %s absent from %s" % (", ".join(missing), src.name))

    rans = {L: e[L]["rans_dCd_total_ct"] for L in LETTERS}
    near = {L: e[L]["vlm_dCDi_near_ct"] for L in LETTERS}
    far = {L: e[L]["vlm_dCDiw_far_ct"] for L in LETTERS}
    r_rans, r_near, r_far = ranks(rans), ranks(near), ranks(far)

    print("%% table body for tab:qv:rank, generated by scripts/make_rank_table.py")
    for L in sorted(LETTERS, key=lambda k: r_rans[k]):
        print("    %s & %d & %d & %d \\\\" % (L, r_rans[L], r_near[L], r_far[L]))
    print("    \\midrule")
    nf, ff = rk["near_field"], rk["far_field"]
    print("    \\multicolumn{2}{l}{Kendall $\\tau$, all five}  & $%+.2f$ & $%+.2f$ \\\\"
          % (nf["kendall"], ff["kendall"]))
    print("    \\multicolumn{2}{l}{Spearman $\\rho$, all five} & $%+.2f$ & $%+.2f$ \\\\"
          % (nf["spearman"], ff["spearman"]))
    print("    \\multicolumn{2}{l}{Pairs swapped, of 10}      & %d & %d \\\\"
          % (len(nf["swapped_pairs"]), len(ff["swapped_pairs"])))

    # ---- tab:qv:increments, section 5 -------------------------------------------------
    # THE SAME NUMBERS APPEAR IN THREE TABLES ACROSS TWO SECTIONS. A stale-value sweep
    # found two of them still carrying the superseded increments after the ranking table
    # had been corrected, which is the half-updated-document failure in its usual form:
    # the table everyone looks at gets fixed and its two siblings do not. All three are
    # emitted here from one computation.
    rng = max(rans.values()) - min(rans.values())

    # "% OF RANGE" IS GONE, by project convention: nobody reads a percentage of a
    # 1.108-count spread as anything meaningful, and it invited confusion with span
    # efficiency, which is a different quantity entirely. The table now carries what the
    # reader actually wants to compare: each method's PREDICTION, the DELTA against the
    # measurement, and the SPAN EFFICIENCIES beside them.
    #
    # E and Ew are the VLM's own near- and far-field span efficiencies from
    # polar_on_DSO_basis. e_x4 is the RANS value at the last Trefftz station, and is
    # absent for any geometry whose wake-plane data is not yet available; absent is
    # printed as "--" rather than omitted, so the gap is visible in the table.
    vlm_E, vlm_Ew = _vlm_span_efficiency()
    rans_e = _rans_span_efficiency()
    print("\n%% table body for tab:qv:increments (section 5)")
    for L in sorted(LETTERS, key=lambda k: r_rans[k]):
        en, ef = near[L] - rans[L], far[L] - rans[L]
        print("    %s & $%+.3f$ & $%+.3f$ & $%+.3f$ & $%+.3f$ & $%+.3f$ & %s & %s & %s \\\\"
              % (L, rans[L], near[L], en, far[L], ef,
                 _fmt(vlm_E.get(L), 4), _fmt(vlm_Ew.get(L), 4), _fmt(rans_e.get(L), 4)))

    # ---- tab:perf:vlm, section 6 ------------------------------------------------------
    ec = EARLY_CRUISE
    r_ec = ranks(ec) if all(L in ec for L in LETTERS) else None
    print("\n%% table body for tab:perf:vlm (section 6)")
    for L in sorted(LETTERS, key=lambda k: r_rans[k]):
        tail = (" & $%+.3f$ & %d" % (ec[L], r_ec[L])) if r_ec else " & -- & --"
        print("%s & $%+.3f$ & %d & $%+.3f$ & %d & $%+.3f$ & %d%s \\\\"
              % (L, near[L], r_near[L], far[L], r_far[L], rans[L], r_rans[L], tail))

    # ---- tab:lo:selection, the ONE selection table of the restructured report ----------
    # The 22 September restructure collapsed tab:qv:rank, tab:qv:increments and
    # tab:perf:vlm into a single table, because three copies of one computation is how two
    # of them came to hold superseded numbers. The three blocks above are kept so that
    # nothing reading them breaks; THIS block is what sections/s6_low_order.tex carries.
    # Column order: RANS increment and rank; near-field increment, rank, error; far-field
    # increment, rank, error; VLM E, VLM Ew, RANS e. Twelve columns, one row per candidate.
    print("\n%% table body for tab:lo:selection (section 6 of the restructured report)")
    for L in sorted(LETTERS, key=lambda k: r_rans[k]):
        en, ef = near[L] - rans[L], far[L] - rans[L]
        print("    %s & $%+.3f$ & %d & $%+.3f$ & %d & $%+.3f$ & $%+.3f$ & %d & $%+.3f$ & %s & %s & %s \\\\"
              % (L, rans[L], r_rans[L], near[L], r_near[L], en, far[L], r_far[L], ef,
                 _fmt(vlm_E.get(L), 4), _fmt(vlm_Ew.get(L), 4), _fmt(rans_e.get(L), 4)))
    print("    \\midrule")
    print("    \\multicolumn{3}{l}{Kendall $\\tau$, all five}  & \\multicolumn{3}{c}{$%+.2f$} "
          "& \\multicolumn{3}{c}{$%+.2f$} & & & \\\\" % (nf["kendall"], ff["kendall"]))
    print("    \\multicolumn{3}{l}{Spearman $\\rho$, all five} & \\multicolumn{3}{c}{$%+.2f$} "
          "& \\multicolumn{3}{c}{$%+.2f$} & & & \\\\" % (nf["spearman"], ff["spearman"]))
    print("    \\multicolumn{3}{l}{Pairs swapped, of 10}      & \\multicolumn{3}{c}{%d} "
          "& \\multicolumn{3}{c}{%d} & & & \\\\"
          % (len(nf["swapped_pairs"]), len(ff["swapped_pairs"])))

    # ---- tab:qv:lift, the trim-angle table --------------------------------------------
    # ADDED 22 September BECAUSE THE TABLE WAS FOUND STALE BY THE RESTRUCTURE: its B, F, W
    # and H rows still carried the trim angles from before those four cases were re-trimmed
    # (B 1.742369 against 1.747609 in this file), and so did the paragraph above it. The
    # table had no generator, which is how it stayed wrong while fig_vlmerr_correlation.pdf,
    # drawn from the same JSON, was right. This block emits it; pasting is a manual decision.
    lift = d.get("lift", {})
    slope = d.get("rans_lift_slope_fit", {})
    alt = lift.get("B_alternative_baseline_corrected")
    print("\n%% table body for tab:qv:lift (section 6.2), generated by scripts/make_rank_table.py")

    def _lift_row(label, L, rec):
        sl = slope.get(L, {}).get("dCl_dalpha_per_deg")
        print("    %s & %.6f & %.5f & $%+.4f$ & $%+.2f$ & %s \\\\"
              % (label, rec["alpha_rans_deg"], rec["alpha_vlm_deg"], rec["d_alpha_deg"],
                 rec["implied_dCL_pct_of_target"], "%.6f" % sl if sl is not None else "--"))
    for L in ["B", "C", "M", "F", "W", "H"]:
        if L in lift:
            _lift_row(L, L, lift[L])
    if alt:
        print("    \\midrule")
        _lift_row("B$^{*}$", "B", alt)
    cand = {L: abs(lift[L]["d_alpha_deg"]) for L in LETTERS if L in lift}
    if cand:
        worst = max(cand, key=cand.get)
        print("%% trim-angle prose: candidates agree to within %.3f deg (worst %s); implied "
              "dCL at most %.2f%% of target; B %+.3f deg, B* %s"
              % (cand[worst], worst,
                 max(abs(lift[L]["implied_dCL_pct_of_target"]) for L in cand),
                 lift["B"]["d_alpha_deg"],
                 ("%+.3f deg" % alt["d_alpha_deg"]) if alt else "--"), file=sys.stderr)

    print("\n%% RANS increment range at Condition CR: %.3f ct (was 6.729 before the "
          "re-trim)" % rng)

    # ---- prose values, so the paragraphs quote a computation and not a memory ----------
    def order_str(r):
        return " < ".join(sorted(r, key=lambda k: r[k]))

    print("\n%% ---- values for the prose ----", file=sys.stderr)
    print("RANS order by dCd      : %s" % order_str(r_rans), file=sys.stderr)
    print("VLM near-field order   : %s" % order_str(r_near), file=sys.stderr)
    print("VLM far-field order    : %s" % order_str(r_far), file=sys.stderr)
    print("near-field  kendall %+.2f  spearman %+.2f  concordant %d  discordant %d  "
          "swapped %d of 10: %s"
          % (nf["kendall"], nf["spearman"], nf["concordant"], nf["discordant"],
             len(nf["swapped_pairs"]),
             ", ".join("-".join(p) for p in nf["swapped_pairs"]) or "none"), file=sys.stderr)
    print("far-field   kendall %+.2f  spearman %+.2f  concordant %d  discordant %d  "
          "swapped %d of 10: %s"
          % (ff["kendall"], ff["spearman"], ff["concordant"], ff["discordant"],
             len(ff["swapped_pairs"]),
             ", ".join("-".join(p) for p in ff["swapped_pairs"]) or "none"), file=sys.stderr)

    # WHICH METRIC WINS, stated as a comparison rather than left to the reader, because
    # this is the sentence that reversed.
    better = "near-field" if nf["kendall"] > ff["kendall"] else "far-field"
    print("\nTHE METRIC THAT TRACKS THE RANS ORDER IS THE %s ONE "
          "(kendall %+.2f against %+.2f)."
          % (better.upper(), max(nf["kendall"], ff["kendall"]),
             min(nf["kendall"], ff["kendall"])), file=sys.stderr)

    # Ties among the RANS increments, named rather than asserted in a caption.
    ties = []
    ordered = sorted(LETTERS, key=lambda k: rans[k])
    for a, b in zip(ordered, ordered[1:]):
        gap = abs(rans[a] - rans[b])
        if gap < TIE_CT:
            ties.append((a, b, gap))
    if ties:
        for a, b, g in ties:
            print("TIE: RANS separates %s and %s by %.3f ct (< %.2f)" % (a, b, g, TIE_CT),
                  file=sys.stderr)
    else:
        adj = min(((a, b, abs(rans[a] - rans[b])) for a, b in zip(ordered, ordered[1:])),
                  key=lambda t: t[2])
        print("NO TIES among the RANS increments; the closest adjacent pair is %s-%s at "
              "%.3f ct. The old caption's \"F and W are separated by 0.001 ct and are a "
              "tie\" no longer holds." % (adj[0], adj[1], adj[2]), file=sys.stderr)

    # The provisional set, named, because a row that is not a measurement must not be
    # quoted as one (the same rule the deck's trim gate applies).
    forces = json.loads((HERE / "data/rans_forces.json").read_text())
    g = forces["trims"]["condition_CR"]["geometries"]
    prov = [L for L in LETTERS if not g[L].get("trim_within_tolerance", False)]
    if prov:
        print("\nPROVISIONAL: %s %s not trimmed to within the gate, so every count above "
              "that involves %s is subject to change."
              % (", ".join(prov), "is" if len(prov) == 1 else "are",
                 "it" if len(prov) == 1 else "them"), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
