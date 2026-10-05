#!/usr/bin/env python3
"""THE reference-quantity guard. Import this before dividing by an aspect ratio.

    from reference_quantities import aspect_ratio, assert_ar_consistent

WHY THIS EXISTS (D062)
----------------------
This programme carries THREE aspect ratios for one wing, because it carries
three reference areas, and two of them differ by 11.3%:

    AR 10.7807   Bref^2 / 13.3572 ft^2   DSO / WETTED planform. This is the
                 basis our RANS coefficients and every DSO comparison use.
    AR 12.0000   Bref^2 / 12.0 ft^2      TP-1580's stated reference area. This
                 is the basis the tunnel data's own coefficients use.
    AR 12.0032   Bref^2 / 11.996793      the trapezoidal planform, which is
                 where the wing's "AR 12" DESIGNATION comes from.

The 12.0 versus 12.0032 distinction is 0.027% and harmless. THE DANGEROUS PAIR
IS 10.7807 VERSUS 12.0, an 11.3% error, and it enters silently: CDi = CL^2 /
(pi AR e) is dimensionally fine with any AR, so taking CL from one normalisation
and AR from the other inflates the induced term with no error and no warning.
Two workstreams walked into it from opposite sides in a single night.

THE RULE: AR MUST BE Bref^2 / Sref COMPUTED FROM THE SAME Sref THE COEFFICIENTS
IN HAND ARE NORMALISED ON. Never hard-code it. Read Sref from the case, the
.polar header or the run card, and derive AR from that.

Prior art, and both do this correctly inline: scripts/analyse_noise_floor.py
(aspect_ratio(), reads Aref from the case and refuses a mismatch above 0.5%) and
scripts/analyse_A2.py (AR from the cases' own Aref). This module exists so the
guard is importable rather than duplicated, and so a third consumer cannot be
written without it.
"""
import math

BREF_FT = 12.0
BREF_M = 3.6576

SREF_DSO_FT2 = 13.3572          # wetted planform; our RANS and all DSO work
SREF_TP1580_FT2 = 12.0          # TP-1580's stated reference area
SREF_TRAPEZOID_FT2 = 11.996793  # source of the "AR 12" designation

AR_DSO = BREF_FT ** 2 / SREF_DSO_FT2            # 10.7807
AR_TP1580 = BREF_FT ** 2 / SREF_TP1580_FT2      # 12.0000
AR_TRAPEZOID = BREF_FT ** 2 / SREF_TRAPEZOID_FT2  # 12.0032


def aspect_ratio(sref, bref=None, half_model=False):
    """AR = bref^2 / sref, in whatever consistent units are passed.

    half_model=True doubles sref first, for a half-model Aref. Units of bref
    default to feet; pass bref explicitly when working in metres.
    """
    if bref is None:
        bref = BREF_FT if sref > 5.0 else BREF_M
    s = 2.0 * sref if half_model else sref
    if s <= 0:
        raise ValueError("sref must be positive, got %r" % sref)
    return bref ** 2 / s


def assert_ar_consistent(ar, sref, bref=None, half_model=False, tol=0.005, label=""):
    """Raise unless `ar` is bref^2/sref for THIS consumer's own sref.

    tol is a RELATIVE tolerance and defaults to 0.5%, which admits the harmless
    12.0/12.0032 distinction and rejects the 11.3% DSO/TP-1580 confusion.
    """
    expect = aspect_ratio(sref, bref, half_model)
    rel = abs(ar - expect) / expect
    if rel > tol:
        raise ValueError(
            "%sAR %.4f is not Bref^2/Sref for this consumer's own Sref %.6f, which "
            "gives %.4f (%.2f%% apart). The programme carries AR 10.7807 on the DSO "
            "wetted area and AR 12.0 on TP-1580's reference area, 11.3%% apart; using "
            "CL from one with AR from the other silently inflates CDi = CL^2/(pi AR e)."
            % (label + ": " if label else "", ar, sref, expect, 100 * rel))
    return expect


def basis_of(ar, tol=0.005):
    """Name the basis an AR belongs to, or None. For diagnostics and messages."""
    for name, val in (("DSO / wetted (Sref 13.3572 ft2)", AR_DSO),
                      ("TP-1580 reference (Sref 12.0 ft2)", AR_TP1580),
                      ("trapezoidal designation (Sref 11.996793 ft2)", AR_TRAPEZOID)):
        if abs(ar - val) / val <= tol:
            return name
    return None




# ---------------------------------------------------------------------------
# CHORD REGISTER (D068). ONE CHORD PER PURPOSE.
# ---------------------------------------------------------------------------
# This wing has FIVE defensible chords and they span 0.3048 to 0.3940 m, a 29%
# range. Every one of them is a correct number; what makes them dangerous is
# that none of them announces its purpose. TP-1580 itself prints three in a
# single sentence (report page 4, read visually at 600 DPI):
#
#   "The aspect-ratio-12 wing had a mean geometric chord of 33.02 cm (13.00 in.)
#    and the aspect-ratio-10 wing, 34.90 cm (13.74 in.); however, an approximate
#    average value of 34.04 cm (13.40 in.) was used as the reference mean
#    geometric chord for both ... This average value also corresponded to the
#    value of the local wing chord at the wing trailing-edge break station."
#
# Note the last clause: 34.04 is NOT a two-wing mean. (33.02+34.90)/2 = 33.96.
# It is the local chord at the trailing-edge break station, and TP-1580 says so.
#
# CHORDS[purpose] -> (metres, source, note). Callers MUST name a purpose.
CHORDS = {
    "re_vs_tp1580": (
        0.3404, "TP-1580 p.4, 34.04 cm (13.40 in.)",
        "TP-1580's OWN reference mean geometric chord, used for BOTH wings and equal "
        "to the local chord at the trailing-edge break station. THE ONLY chord for "
        "Reynolds matching against TP-1580, because it is the chord its tabulated Re "
        "is built on."),
    "similarity": (
        0.3302, "TP-1580 p.4, 33.02 cm (13.00 in.)",
        "the AR-12 wing's OWN mean geometric chord. Use for geometric-similarity "
        "arguments ONLY. 0.625 x 33.02 = 20.6375 cm reproduces TM-83111's printed "
        "20.64 cm; against the 34.04 reference chord the ratio comes out 0.6063 and "
        "looks like a similarity failure when it is a frame collision (D067)."),
    "dso_internal": (
        0.3393, "S/b on Sref 13.3572 ft2 = 1.1131 ft",
        "mean geometric chord on the DSO WETTED area. DSO-internal use only; never "
        "against TP-1580, whose coefficients are on the 12.0 ft2 tunnel area."),
    "moment_ref": (
        0.393957, "DSO Cref (MAC), 1.29251 ft",
        "the MEAN AERODYNAMIC chord, not a mean geometric chord. This is the chord "
        "forceCoeffs divides Cm by. It differs from dso_internal by 16% and the two "
        "are routinely confused (D040 item 4)."),
    "trapezoidal_sb": (
        0.3048, "S/b on the 12.0 ft2 trapezoidal area = 1.0 ft",
        "NO SANCTIONED USE. Recorded so that a value of exactly 1.0 ft or 0.3048 m "
        "appearing as a chord is recognised as this quantity and rejected, rather "
        "than mistaken for a feet-to-metres conversion factor."),
    # THE 0.3048 QUADRUPLE COLLISION (GEO-084). The Phase 0.5 synthetic wings make
    # 0.3048 m legitimate for the first time, and the temptation is to relax the
    # refusal above. DO NOT. A guard loosened once to admit a legitimate case admits
    # the illegitimate one next. Two NEW NAMED purposes instead, so the refusal on
    # trapezoidal_sb stays intact for every other caller.
    "synthetic_rectangular_chord": (
        0.3048, "syn_rectangular constant chord, S/b with S = b^2/12",
        "the SYNTHETIC rectangular wing's actual chord. Phase 0.5 method validation "
        "ONLY; never against an ARGUS or TP-1580 quantity."),
    "synthetic_reference_mgc": (
        0.3048, "mean geometric chord of BOTH synthetic wings, S/b",
        "mean geometric chord of both synthetic wings, equal by construction because "
        "they share S and b. Phase 0.5 method validation ONLY."),
}
UNSANCTIONED = {"trapezoidal_sb"}

# FOUR DISTINCT MEANINGS SHARE THE NUMERAL 0.3048, which is worse than the 1.1131
# collision because one of the four is a UNIT CONVERSION (GEO-084):
#   1. syn_rectangular's chord                        (synthetic_rectangular_chord)
#   2. the mean geometric chord of BOTH synthetic wings (synthetic_reference_mgc)
#   3. the feet-to-metres conversion factor            <- the live hazard
#   4. S/b on the 12.0 ft2 trapezoidal area            (trapezoidal_sb, REFUSED)
# A REFACTOR THAT DE-DUPLICATES THIS NUMBER IS A DEFECT. Any code or reviewer that
# sees 0.3048 and infers a unit conversion is wrong on three of the four.
COLLIDING_NUMERAL_0P3048 = 0.3048


def chord_for(purpose):
    """Return the chord in metres for a NAMED purpose. Raises otherwise.

    There is no default and there will not be one: a caller that does not know
    which chord it needs cannot be given one safely.
    """
    if purpose not in CHORDS:
        raise ValueError(
            "chord_for() requires a NAMED purpose. Got %r. Valid: %s. This wing has "
            "five defensible chords spanning 0.3048 to 0.3940 m (29%%); picking one "
            "without stating why is how D067's 0.6063-versus-0.6251 collision happened."
            % (purpose, sorted(CHORDS)))
    if purpose in UNSANCTIONED:
        raise ValueError(
            "chord %r has NO SANCTIONED USE in this programme: %s"
            % (purpose, CHORDS[purpose][2]))
    return CHORDS[purpose][0]


def reynolds_on_chord(re, from_purpose, to_purpose):
    """Convert a Reynolds number from one chord basis to another. Re scales with chord."""
    return re * chord_for(to_purpose) / chord_for(from_purpose)


# TP-1580's own stated (Mach, Re) couple, report p.9, read as an AUTHORITATIVE
# PAIR and never re-derived from Mach with our nu (D042 item 3c):
#   "Reynolds numbers of 0.97 to 1.63 x 10^6 based on the reference mean geometric
#    chord of 34.04 cm, with corresponding Mach numbers of 0.12 to 0.20"
TP1580_RE_COUPLE = ((0.12, 0.97e6), (0.20, 1.63e6))


def assert_re_basis_vs_tp1580(re, mach, chord_purpose, tol=0.05, label=""):
    """Raise unless `re`, once converted to TP-1580's chord basis, matches its
    tabulated Re at this Mach.

    Re is proportional to Mach in an atmospheric tunnel, so the couple is
    interpolated linearly in Mach. tol is relative and defaults to 5%, wide
    enough for the read precision and far tighter than the 15.7% error of using
    the DSO MAC without converting.
    """
    (m0, r0), (m1, r1) = TP1580_RE_COUPLE
    expect = r0 + (re_i := (mach - m0) / (m1 - m0)) * (r1 - r0)
    got = reynolds_on_chord(re, chord_purpose, "re_vs_tp1580")
    rel = abs(got - expect) / expect
    if rel > tol:
        raise ValueError(
            "%sRe %.4g on chord basis %r converts to %.4g on TP-1580's 0.3404 m chord, "
            "against its tabulated %.4g at Mach %.3f (%.1f%% apart). Using the DSO MAC "
            "0.393957 m without converting is a 15.7%% error and looks like a physical "
            "disagreement."
            % (label + ": " if label else "", re, chord_purpose, got, expect, mach, 100 * rel))
    return expect




# ---------------------------------------------------------------------------
# AREA-BASIS LOCK (D068 item 3). ASSERTION, NOT A COMMENT.
# ---------------------------------------------------------------------------
# TP-1580 and TM-83111 coefficients are on the 12.0 ft^2 TUNNEL area. DSO and our
# RANS are on the 13.3572 ft^2 WETTED area. The factor is 1.1131, i.e. 11.31%,
# against a morphing effect of order 0.4%. A missed conversion is 28x the signal.
#
# THE D036 TRAP LIVES EXACTLY HERE, and it is worth restating because the numeral
# is shared: the area ratio 13.3572/12.0 = 1.113100 is DIMENSIONLESS, and the DSO
# mean geometric chord S/b = 13.3572/12.0 = 1.113100 ft is a LENGTH. Same digits,
# because bref = 12.0 ft and Sref_TP1580 = 12.0 ft^2 share the numeral 12.
# Anything landing on 1.1131 must say which one it means.
AREA_RATIO_DSO_OVER_TP1580 = SREF_DSO_FT2 / SREF_TP1580_FT2   # 1.113100, DIMENSIONLESS
# DELIBERATELY A SECOND, SEPARATE DEFINITION. Do NOT merge this with the line
# above as duplication: this one is a LENGTH IN FEET, that one is a RATIO. They
# share the numeral 1.113100 only because bref = 12.0 ft and Sref_TP1580 = 12.0
# ft^2 share the numeral 12. regression_test_area_factor() asserts they remain
# separately defined and fails if either is aliased to the other.
DSO_MEAN_GEOM_CHORD_FT = SREF_DSO_FT2 / BREF_FT               # 1.113100 FEET, a LENGTH
VALID_BASES = ("dso_wetted", "tp1580_tunnel")


def assert_same_basis(basis_a, basis_b, label=""):
    """Raise unless two quantities declare the SAME reference-area basis.

    Callers must declare a basis for each side. There is no inference: a
    coefficient does not carry its basis in its value, which is the entire
    reason this exists.
    """
    for b in (basis_a, basis_b):
        if b not in VALID_BASES:
            raise ValueError("%sunknown reference-area basis %r; valid: %s"
                             % (label + ": " if label else "", b, list(VALID_BASES)))
    if basis_a != basis_b:
        raise ValueError(
            "%srefusing to compare a %r quantity with a %r one. The reference areas "
            "differ by a factor of %.6f (%.2f%%) against a morphing effect of order "
            "0.4%%, so an unconverted comparison is about 28x the signal. Convert "
            "explicitly with scripts/tp1580_normalisation.py."
            % (label + ": " if label else "", basis_a, basis_b,
               AREA_RATIO_DSO_OVER_TP1580, 100 * (AREA_RATIO_DSO_OVER_TP1580 - 1)))
    return basis_a


def regression_test_area_factor():
    """Fails if the 1.1131 conversion is removed or neutered anywhere.

    Three independent checks, so weakening any one of them still trips the test:
      1. the factor is 1.113100 from the two registered areas;
      2. it is NOT 1.0, which is what removing the conversion would make it;
      3. tp1580_normalisation actually APPLIES it, checked round trip;
      4. THE TWO 1.113100 CONSTANTS REMAIN SEPARATELY DEFINED. This is the check
         the numeral coincidence actually invites: a refactor sees 1.113100
         defined twice, calls it duplication, and merges them. Both call sites
         keep working, and checks 1 to 3 ALL KEEP PASSING, because the value is
         right either way. What is lost is that one is a dimensionless ratio and
         the other is a length in feet.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import tp1580_normalisation as tn

    ok = True
    if abs(AREA_RATIO_DSO_OVER_TP1580 - 1.113100) > 1e-6:
        print("  FAIL area factor is %.6f, expected 1.113100" % AREA_RATIO_DSO_OVER_TP1580); ok = False
    if abs(AREA_RATIO_DSO_OVER_TP1580 - 1.0) < 1e-9:
        print("  FAIL area factor collapsed to 1.0, the conversion has been removed"); ok = False
    cl_tp = 0.4
    cl_dso = tn.cl_tp1580_to_dso(cl_tp)
    if abs(cl_dso - cl_tp) < 1e-9:
        print("  FAIL cl_tp1580_to_dso is the identity; the factor is not applied"); ok = False
    if abs(tn.cl_dso_to_tp1580(cl_dso) - cl_tp) > 1e-9:
        print("  FAIL round trip does not return the original CL"); ok = False

    # 4. cross-purpose aliasing. The two constants must be SEPARATE objects with
    #    separate derivations, not one name pointing at the other.
    import ast, inspect
    src = inspect.getsource(sys.modules[__name__])
    tree = ast.parse(src)
    defs = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            t = node.targets[0]
            if isinstance(t, ast.Name) and t.id in ("AREA_RATIO_DSO_OVER_TP1580",
                                                    "DSO_MEAN_GEOM_CHORD_FT"):
                defs[t.id] = ast.dump(node.value)
    missing = {"AREA_RATIO_DSO_OVER_TP1580", "DSO_MEAN_GEOM_CHORD_FT"} - set(defs)
    if missing:
        print("  FAIL constant(s) no longer separately defined: %s" % sorted(missing)); ok = False
    elif len(set(defs.values())) < 2:
        print("  FAIL the two 1.113100 constants have been merged or aliased; one is a "
              "RATIO and the other is a LENGTH IN FEET"); ok = False
    else:
        for name in ("AREA_RATIO_DSO_OVER_TP1580", "DSO_MEAN_GEOM_CHORD_FT"):
            if "Name" in defs[name] and "BinOp" not in defs[name]:
                print("  FAIL %s is aliased to another name rather than derived" % name); ok = False
        print("  the two 1.113100 constants are separately DERIVED (ratio vs length) : OK")
    print("  area factor %.6f, cl 0.4 tp1580 -> %.6f dso, round trip exact : %s"
          % (AREA_RATIO_DSO_OVER_TP1580, cl_dso, "PASS" if ok else "FAIL"))
    return ok




# ---------------------------------------------------------------------------
# OPERATING-POINT LOCK (D077). THE LAST FRAME INSTANCE THAT LIVED IN PROSE.
# ---------------------------------------------------------------------------
# Instance 6 of the frame class was the only one with no code path guarding it,
# and prose is what passes every automated check. This is the guard.
#
# TWO DIFFERENT M 0.1 STATES EXIST and they are not interchangeable:
#   (0.1, 0.529297)         the DSO report's four-concept matched study; the
#                           report's own words, "a low-speed method-sensitivity
#                           case, not an A320 weight-supporting flight state"
#   (0.1, 0.428277635108)   the ARGUS validation point; every number in this repo
#
# Same Mach, same concept names, different operating point. A percentage from one
# may not be quoted against a count from the other.
OPERATING_POINTS = {
    "argus_validation": (0.10, 0.428277635108),
    "dso_method_sensitivity": (0.10, 0.529297),
}


def operating_point(mach, cl):
    """Normalise an (mach, cl) pair and name it if it is a known state."""
    for name, (m, c) in OPERATING_POINTS.items():
        if abs(mach - m) < 1e-9 and abs(cl - c) < 1e-9:
            return name, (mach, cl)
    return None, (mach, cl)


def assert_same_operating_point(a, b, tol_mach=1e-6, tol_cl=1e-6, label=""):
    """Raise unless two quantities carry the SAME (mach, cl). No inference.

    a and b are (mach, cl) tuples. There is no default operating point and there
    will not be one: a quantity that does not know its own state cannot be
    compared safely, exactly as chord_for() refuses an unqualified chord.
    """
    for x in (a, b):
        if not (isinstance(x, (tuple, list)) and len(x) == 2):
            raise ValueError("%soperating point must be a (mach, cl) pair, got %r"
                             % (label + ": " if label else "", x))
    dm, dc = abs(a[0] - b[0]), abs(a[1] - b[1])
    if dm > tol_mach or dc > tol_cl:
        na, _ = operating_point(*a)
        nb, _ = operating_point(*b)
        raise ValueError(
            "%srefusing to compare a quantity at (M %.4f, CL %.9f)%s with one at "
            "(M %.4f, CL %.9f)%s. Same Mach with a different CL is still a "
            "different operating point: the DSO report's four-concept study sits "
            "at CL 0.529297 and every number in this repo sits at 0.428277635108. "
            "State BOTH CLs inline or do not compare them."
            % (label + ": " if label else "", a[0], a[1], " [%s]" % na if na else "",
               b[0], b[1], " [%s]" % nb if nb else ""))
    return a


def self_test():
    ok = True
    print("AR_DSO        %.4f  <- our RANS and every DSO comparison" % AR_DSO)
    print("AR_TP1580     %.4f  <- the tunnel data's own coefficients" % AR_TP1580)
    print("AR_TRAPEZOID  %.4f  <- the 'AR 12' designation only" % AR_TRAPEZOID)
    print("dangerous pair %.2f%% apart; harmless pair %.3f%% apart"
          % (100 * (AR_TP1580 / AR_DSO - 1), 100 * (AR_TRAPEZOID / AR_TP1580 - 1)))
    print()
    # the harmless distinction must PASS
    try:
        assert_ar_consistent(AR_TRAPEZOID, SREF_TP1580_FT2, label="harmless 12.0032 vs 12.0")
        print("harmless 12.0032 against Sref 12.0 : PASS (accepted, 0.027%)")
    except ValueError as e:
        ok = False
        print("harmless pair REJECTED, tolerance too tight:", e)
    # the dangerous confusion must RAISE
    try:
        assert_ar_consistent(AR_TP1580, SREF_DSO_FT2, label="AR 12 with DSO Sref")
        print("dangerous AR 12.0 against Sref 13.3572 : NOT CAUGHT")
        ok = False
    except ValueError:
        print("dangerous AR 12.0 against Sref 13.3572 : CAUGHT")
    # half-model form, which is how OpenFOAM cases store Aref
    half = SREF_DSO_FT2 / 2.0
    got = aspect_ratio(half, half_model=True)
    print("half-model Aref %.4f ft2 -> AR %.4f %s"
          % (half, got, "OK" if abs(got - AR_DSO) < 1e-9 else "FAIL"))
    ok &= abs(got - AR_DSO) < 1e-9
    # --- chord register (D068)
    print("CHORD REGISTER")
    for k, (m, src, _) in CHORDS.items():
        print("  %-16s %.6f m   %s%s" % (k, m, src, "   [NO SANCTIONED USE]" if k in UNSANCTIONED else ""))
    try:
        chord_for("cref"); print("  unnamed purpose 'cref' : NOT CAUGHT"); ok = False
    except ValueError:
        print("  unnamed purpose 'cref' : CAUGHT")
    try:
        chord_for("trapezoidal_sb"); print("  unsanctioned chord     : NOT CAUGHT"); ok = False
    except ValueError:
        print("  unsanctioned chord     : CAUGHT")

    # --- Reynolds basis gate, two runs
    print()
    print("REYNOLDS BASIS GATE (TP-1580 couple 0.97e6 @ M0.12, 1.63e6 @ M0.20 on 0.3404 m)")
    for mach, re_tp in ((0.12, 0.97e6), (0.20, 1.63e6)):
        re_mac = re_tp * CHORDS["moment_ref"][0] / CHORDS["re_vs_tp1580"][0]
        try:
            assert_re_basis_vs_tp1580(re_mac, mach, "moment_ref", label="M %.2f" % mach)
            print("  M %.2f  Re on MAC %.4g -> converts and MATCHES" % (mach, re_mac))
        except ValueError as e:
            print("  M %.2f  unexpected failure: %s" % (mach, e)); ok = False
        try:
            assert_re_basis_vs_tp1580(re_mac, mach, "re_vs_tp1580")
            print("  M %.2f  same number claimed already on 0.3404 m : NOT CAUGHT" % mach); ok = False
        except ValueError:
            print("  M %.2f  same number claimed already on 0.3404 m : CAUGHT (15.7%% error)" % mach)

    print()
    print("AREA-BASIS LOCK (D068)")
    ok &= regression_test_area_factor()
    try:
        assert_same_basis("dso_wetted", "tp1580_tunnel", label="mixed comparison")
        print("  mixed-basis comparison : NOT CAUGHT"); ok = False
    except ValueError:
        print("  mixed-basis comparison : CAUGHT")
    try:
        assert_same_basis("dso_wetted", "wetted")
        print("  undeclared basis       : NOT CAUGHT"); ok = False
    except ValueError:
        print("  undeclared basis       : CAUGHT")

    print()
    print("OPERATING-POINT LOCK (D077)")
    argus=OPERATING_POINTS["argus_validation"]; dso=OPERATING_POINTS["dso_method_sensitivity"]
    try:
        assert_same_operating_point(argus, argus); print("  same state, same state : accepted")
    except ValueError: print("  same state REJECTED"); ok=False
    try:
        assert_same_operating_point(argus, dso, label="exec summary vs our pair")
        print("  M 0.1 CL 0.4283 vs M 0.1 CL 0.5293 : NOT CAUGHT"); ok=False
    except ValueError:
        print("  M 0.1 CL 0.4283 vs M 0.1 CL 0.5293 : CAUGHT (same Mach, different CL)")
    try:
        assert_same_operating_point((0.1, 0.4283), 0.4283); print("  malformed pair : NOT CAUGHT"); ok=False
    except ValueError: print("  malformed pair                     : CAUGHT")
    print("  named states:", ", ".join("%s (M %.2f, CL %.9f)"%(k,v[0],v[1]) for k,v in OPERATING_POINTS.items()))

    print()
    print("SELF TEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(self_test())
