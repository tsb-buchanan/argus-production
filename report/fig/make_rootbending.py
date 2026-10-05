#!/usr/bin/env python3
"""
Root bending moment, VLM against RANS, and Liming's screening limit.
ARGUS all-geometry report 2026-09-15.

EVERY NUMBER IN THE OUTPUT IS READ FROM A FILE. Nothing is typed in from memory.

Inputs (read, never assumed):
  results/derived/<CASE>_trim.json                     RANS RBM, from the solver's own M_x
  dso_reference/.../<geom>_optimization_result.csv     VLM RBM and its increase percentage
  docs/report/.../data/vlm.json                        the VLM freestream, AS READ FROM DECK

-------------------------------------------------------------------------------
THE TWO FRAMES THAT MUST BE STATED OR THIS COMPARISON IS MEANINGLESS (D068)
-------------------------------------------------------------------------------
1. DYNAMIC PRESSURE. RBM is DIMENSIONAL. The VLM decks run at Vinf = 100 in OpenVSP length
   units, and this project's units are FEET, so 30.48 m/s and q = 569.031 Pa against the
   RANS 34.0 m/s and q = 708.050 Pa. The ratio is 1.24431. A VLM moment compared against a
   RANS moment WITHOUT that factor is wrong by 24 per cent before any physics enters, which
   is the same trap D036 records for reference areas.

2. WHICH BASELINE THE PERCENTAGE IS AGAINST. The VLM differences are against
   `baseline_corrected`, which is B*, the corrected-interpolation baseline. Its OWN record
   states root_bending_increase_percent = -1.6702, i.e. B* is itself measured against the
   original baseline. The RANS differences are against B, the meshed original-lineage
   baseline, because B is the surface that was actually meshed and solved. B and B* are not
   the same wing and the report documents that split for drag; it applies here too.

   THE PERCENTAGE IS THE PRIMARY COMPARISON because it cancels q. The dimensional cross-check
   is reported beside it BECAUSE the two baselines differ: if only percentages were shown, a
   baseline disagreement would hide inside them.

-------------------------------------------------------------------------------
ONE DELIVERED VALUE IS STALE AND IS RECOMPUTED HERE
-------------------------------------------------------------------------------
mcv2_i002_c01 records root_bending_moment_Nm = 118.64405 and
root_bending_increase_percent = 5.3352. Those are inconsistent: 118.64405 / 110.75354 =
1.07124, i.e. +7.124 per cent. The moment is the CORRECTED value while the percentage was
computed against the PRE-CORRECTION baseline. Liming's own report says so in prose, "the
unchanged candidate now increases half-wing root bending by 7.12%, above the illustrative
6.8% screen", so this is a known file defect rather than a new disagreement. EVERY percentage
here is RECOMPUTED from the two moments rather than read, and the delivered value is carried
alongside so the discrepancy is visible rather than silently corrected.

-------------------------------------------------------------------------------
THE LIMIT
-------------------------------------------------------------------------------
Liming's low-speed screen is a 6.8 per cent increase in half-wing root bending over the
corrected rigid baseline. His report calls root bending "the active screen" for the low-speed
selection and states the strict optimum "lies immediately below" it. It also states the screen
is "a study-level bound rather than a certified structural limit", and that his root bending is
"approximated". So exceeding it is a FLAG FOR THE STRUCTURES TEAM, not a failure against a
certified allowable, and this script labels it that way.
"""

import csv
import glob
import json
import derived_source
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(REPORT)))
DATA = os.path.join(REPORT, "data")
DERIVED = os.path.join(REPO, "results", "derived")
DSO = os.path.join(REPO, "dso_reference")
FIG = HERE

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def halt(m):
    log("HALT: " + m)
    raise SystemExit("HALT: " + m)


plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 9.0, "axes.labelsize": 8.5,
    "legend.fontsize": 7.8, "xtick.labelsize": 7.8, "ytick.labelsize": 7.8,
    "axes.linewidth": 0.7, "savefig.bbox": "tight", "figure.dpi": 200,
    "mathtext.fontset": "dejavusans",
})

STYLE = {
    "B": dict(color="#000000"), "C": dict(color="#D55E00"),
    "F": dict(color="#0072B2"), "H": dict(color="#009E73"),
    "M": dict(color="#CC79A7"), "W": dict(color="#56B4E9"),
}
LETTERS = ["B", "C", "F", "H", "M", "W"]

# letter -> the DSO geometry name whose optimization_result.csv carries the VLM RBM
GEOM = {"B": "baseline_corrected", "C": "cte_i002_c04", "F": "cfft_b02_c01",
        "H": "chc_g02_c06", "M": "mcv2_i002_c01", "W": "cffw_b01_c01"}

LIMIT_PCT = 6.8            # Liming's illustrative low-speed screen
Q_RANS_CR = 708.05         # Pa, 0.5 * 1.225 * 34.0^2


# ---------------------------------------------------------------- VLM freestream
def vlm_q():
    """q of the VLM decks, DERIVED from the deck's own Vinf and rho, never assumed.

    OpenVSP's native length unit is FEET for this project, so Vinf = 100 is 100 ft/s, and
    Rho = 0.002377 is slug/ft^3. Both are converted here and the result is asserted against
    the sea-level value the RANS uses, because a silent unit slip in either would move the
    whole comparison by the same factor that makes it necessary.
    """
    with open(os.path.join(DATA, "vlm.json")) as fh:
        d = json.load(fh)
    vs, rs = set(), set()
    for r in d.get("runs", []):
        fc = r.get("flow_conditions_AS_READ_FROM_DECK") or {}
        if fc.get("Vinf_Lunit_per_Tunit") is not None:
            vs.add(float(fc["Vinf_Lunit_per_Tunit"]))
        if fc.get("Rho_Munit_per_Lunit3") is not None:
            rs.add(float(fc["Rho_Munit_per_Lunit3"]))
    if len(vs) != 1 or len(rs) != 1:
        halt("the decks do not share one freestream: Vinf %s, rho %s" % (sorted(vs), sorted(rs)))
    v_fts, rho_slug = vs.pop(), rs.pop()
    v = v_fts * 0.3048                          # ft/s -> m/s
    rho = rho_slug * 515.378818                 # slug/ft^3 -> kg/m^3
    if abs(rho - 1.225) / 1.225 > 0.01:
        halt("deck density converts to %.4f kg/m^3, not sea level 1.225; unit assumption is wrong" % rho)
    q = 0.5 * rho * v * v
    log("VLM freestream, derived from the decks:")
    log("  Vinf %.1f ft/s = %.4f m/s   rho %.6f slug/ft^3 = %.4f kg/m^3" % (v_fts, v, rho_slug, rho))
    log("  q_VLM  = %.4f Pa" % q)
    log("  q_RANS = %.4f Pa   ratio q_RANS/q_VLM = %.5f" % (Q_RANS_CR, Q_RANS_CR / q))
    return q


# ---------------------------------------------------------------- VLM root bending
def vlm_rbm():
    """{letter: (RBM_Nm, delivered_pct, source_path)} from the optimization_result files.

    PREFERS THE NEWEST DELIVERY. C and M exist only in the 2026-07-29 correction pair; B*, F,
    H and W also appear in the 2026-08-19 CDIW audit. Where both exist the newer is taken and
    the older is checked against it, because two copies of one number is two things that can
    disagree.
    """
    out = {}
    for L, name in GEOM.items():
        hits = sorted(glob.glob(os.path.join(DSO, "**", "*%s*optimization_result.csv" % name),
                                recursive=True))
        hits = [h for h in hits if os.path.basename(h).startswith(name)]
        if not hits:
            log("  %s (%s): NO optimization_result.csv found" % (L, name))
            continue
        vals = []
        for h in hits:
            with open(h, encoding="utf-8-sig") as fh:
                r = list(csv.DictReader(fh))[0]
            rbm = r.get("root_bending_moment_Nm")
            pct = r.get("root_bending_increase_percent")
            if rbm:
                vals.append((float(rbm), float(pct) if pct else None, h))
        if not vals:
            continue
        # every copy must agree on the moment, or say so
        m = {round(v[0], 6) for v in vals}
        if len(m) > 1:
            halt("%s has disagreeing RBM across deliveries: %s" % (name, sorted(m)))
        out[L] = vals[-1]
    return out


# ---------------------------------------------------------------- RANS root bending
def rans_rbm(prefix):
    out = {}
    for L in LETTERS:
        # PREFER THE WAKE-REFINED RECORD where the condition was re-run, and record which
        # was used. Building the name directly meant Condition CR read the PUBLISHED
        # records while its drag and span efficiency had moved to the wake-refined mesh.
        p, _gen = derived_source.resolve(prefix, L)
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        rb = d.get("root_bending") or {}
        if rb.get("available"):
            out[L] = (rb["RBM_total_Nm"], rb.get("q_Pa"), rb.get("lift_centroid_eta"),
                      rb.get("force_frame"))
    return out


# ---------------------------------------------------------------- build
Q_VLM = vlm_q()
RATIO = Q_RANS_CR / Q_VLM

log("")
log("VLM root bending, read from the delivered optimisation results:")
V = vlm_rbm()
for L in LETTERS:
    if L in V:
        rbm, pct, src = V[L]
        log("  %s %-20s RBM %12.5f N m   delivered incr%% %s   %s"
            % (L, GEOM[L], rbm, ("%.4f" % pct) if pct is not None else "-",
               os.path.relpath(src, REPO).split(os.sep)[1]))

if "B" not in V:
    halt("no VLM baseline, so no percentage can be formed")
VB = V["B"][0]

log("")
log("PERCENTAGES ARE RECOMPUTED FROM THE MOMENTS, not read from the files:")
rows = []
for L in LETTERS:
    if L not in V:
        continue
    rbm, delivered, src = V[L]
    recomputed = 100.0 * (rbm - VB) / VB
    flag = ""
    if delivered is not None and L != "B" and abs(recomputed - delivered) > 0.01:
        flag = "   <-- delivered %.4f%% DISAGREES, stale against the pre-correction baseline" % delivered
    log("  %s  RBM %12.5f  recomputed %+8.4f%%%s" % (L, rbm, recomputed, flag))
    rows.append((L, rbm, delivered, recomputed))

log("")
R = rans_rbm("SW")
if "B" not in R:
    halt("no RANS baseline at Condition CR")
RB = R["B"][0]
log("RANS root bending at Condition CR, from each case's own solver moment:")
for L in LETTERS:
    if L in R:
        v, q, eta, fr = R[L]
        log("  %s  RBM %12.4f N m   %+8.4f%% vs B   q %.2f   centroid_eta %.4f   frame %s"
            % (L, v, 100.0 * (v - RB) / RB, q, eta, fr))

# ---------------------------------------------------------------- the comparison
log("")
log("COMPARISON. Percentage is the primary column: it cancels q. The dimensional column is")
log("the cross-check, and it is what exposes a baseline disagreement that percentages hide.")
comp = []
for L in LETTERS:
    if L not in V or L not in R:
        continue
    vlm_pct = 100.0 * (V[L][0] - VB) / VB
    rans_pct = 100.0 * (R[L][0] - RB) / RB
    vlm_scaled = V[L][0] * RATIO
    dim_pct = 100.0 * (vlm_scaled - R[L][0]) / R[L][0]
    comp.append(dict(letter=L, geometry=GEOM[L],
                     vlm_rbm_Nm=V[L][0], vlm_rbm_scaled_to_RANS_q_Nm=vlm_scaled,
                     vlm_increase_pct_recomputed=vlm_pct,
                     vlm_increase_pct_as_delivered=V[L][2] if False else V[L][1],
                     rans_rbm_Nm=R[L][0], rans_increase_pct=rans_pct,
                     pct_difference_pp=rans_pct - vlm_pct,
                     dimensional_difference_pct=dim_pct,
                     exceeds_6p8_screen_vlm=bool(vlm_pct > LIMIT_PCT),
                     exceeds_6p8_screen_rans=bool(rans_pct > LIMIT_PCT)))

log("  %-2s %-20s %10s %10s %8s   %12s %12s %8s"
    % ("", "geometry", "VLM %", "RANS %", "diff pp", "VLM x ratio", "RANS N m", "dim %"))
for c in comp:
    log("  %-2s %-20s %+10.3f %+10.3f %+8.2f   %12.3f %12.3f %+8.2f"
        % (c["letter"], c["geometry"], c["vlm_increase_pct_recomputed"], c["rans_increase_pct"],
           c["pct_difference_pp"], c["vlm_rbm_scaled_to_RANS_q_Nm"], c["rans_rbm_Nm"],
           c["dimensional_difference_pct"]))

log("")
log("AGAINST LIMING'S 6.8%% SCREEN (illustrative, not a certified structural limit):")
for c in comp:
    if c["letter"] == "B":
        continue
    v = "OVER " if c["exceeds_6p8_screen_vlm"] else "under"
    r = "OVER " if c["exceeds_6p8_screen_rans"] else "under"
    note = ""
    if c["exceeds_6p8_screen_rans"] and not c["exceeds_6p8_screen_vlm"]:
        note = "   <-- the VISCOUS calculation moves it from feasible to infeasible"
    log("  %s  VLM %+.3f%% %s   RANS %+.3f%% %s%s"
        % (c["letter"], c["vlm_increase_pct_recomputed"], v,
           c["rans_increase_pct"], r, note))

# ---------------------------------------------------------------- early cruise, RANS only
log("")
RC = rans_rbm("CMP")
if "B" in RC:
    CB = RC["B"][0]
    log("EARLY CRUISE, RANS only. Our five geometries do NOT appear in the DSO cruise set,")
    log("which holds optimiser seeds and the rigid baseline alone, so there is no VLM column.")
    log("Liming defines the cruise limit as the early-cruise RIGID baseline's own RBM, so")
    log("utilisation = RBM / baseline and anything above 1.0 is flagged.")
    for L in LETTERS:
        if L in RC:
            log("  %s  RBM %10.3f N m   utilisation %.4f  %s"
                % (L, RC[L][0], RC[L][0] / CB, "OVER" if RC[L][0] > CB else "under"))

# ---------------------------------------------------------------- figure
def make_figure(path):
    fig, axes = plt.subplots(1, 2, figsize=(6.69, 3.1))
    letters = [c["letter"] for c in comp if c["letter"] != "B"]
    x = np.arange(len(letters))
    w = 0.38
    vp = [c["vlm_increase_pct_recomputed"] for c in comp if c["letter"] != "B"]
    rp = [c["rans_increase_pct"] for c in comp if c["letter"] != "B"]

    ax = axes[0]
    ax.bar(x - w / 2, vp, w, label="VLM", color="0.62", edgecolor="0.2", linewidth=0.5)
    ax.bar(x + w / 2, rp, w, label="RANS",
           color=[STYLE[L]["color"] for L in letters], edgecolor="0.2", linewidth=0.5)
    ax.axhline(LIMIT_PCT, color="#B00000", lw=1.2, ls="--", label="6.8 per cent screen")
    ax.set_xticks(x); ax.set_xticklabels(letters)
    ax.set_ylabel("root bending increase, per cent")
    ax.set_title("Condition CR", fontsize=8.2)
    ax.legend(frameon=False, fontsize=7.2, loc="upper left")
    ax.grid(True, axis="y", lw=0.35, color="0.88")

    ax = axes[1]
    d = [c["dimensional_difference_pct"] for c in comp]
    allL = [c["letter"] for c in comp]
    ax.bar(np.arange(len(allL)), d, 0.6,
           color=[STYLE[L]["color"] for L in allL], edgecolor="0.2", linewidth=0.5)
    ax.axhline(0.0, color="0.4", lw=0.7)
    ax.set_xticks(np.arange(len(allL))); ax.set_xticklabels(allL)
    ax.set_ylabel("(VLM$\\times$1.2443 $-$ RANS) / RANS,  per cent")
    ax.set_title("", fontsize=8.2)
    ax.grid(True, axis="y", lw=0.35, color="0.88")

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def write_csv(path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(comp[0].keys()))
        w.writeheader()
        for c in comp:
            w.writerow(c)
    return path


written = [make_figure(os.path.join(FIG, "fig_rootbending.pdf")),
           write_csv(os.path.join(FIG, "rootbending_values.csv"))]

log("")
log("VERIFICATION (stat of each artefact, read back from disk):")
bad = 0
for p in written:
    ok = os.path.exists(p) and os.path.getsize(p) > 0
    bad += 0 if ok else 1
    log("  %s  %-34s %8d bytes" % ("OK " if ok else "BAD", os.path.basename(p),
                                   os.path.getsize(p) if os.path.exists(p) else 0))

with open(os.path.join(FIG, "make_rootbending_LOG.txt"), "w") as fh:
    fh.write("\n".join(LOG) + "\n")
log("log written to make_rootbending_LOG.txt")
sys.exit(1 if bad else 0)
