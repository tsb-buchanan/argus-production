#!/usr/bin/env python3
"""hinge_moment.py: the ARGUS hinge moment, and a REFUSAL for every other route to it.

RULING (project decision, 2026-07-29), and it is a refusal rather than a note because the wrong
route is the one a careful person would reach for first:

    OUR HINGE MOMENT COMES FROM RANS PRESSURE INTEGRATION ABOUT x_h/c = 0.62
    AND NOTHING ELSE. No cmy, no transfer formula, no exception.

WHY THE OBVIOUS ROUTE IS FORBIDDEN. The delivered spanwise-loads files carry a column
named sectional_cmy_native. It is exactly the quantity someone would reach for, it is
per-strip, and it is already in SI-labelled files. It is ALSO flagged
"native-unverified" in the author's own DATA_DICTIONARY.csv, with the instruction
"check VSPAERO documentation before dimensional use". Independently, D053 could not
close the delivered Cm from cmy under ANY simple area weighting. Two independent
reasons, one from the author and one from our own arithmetic, and neither has been
resolved. A quantity that cannot be closed against a known total must not be
integrated into a new one.

THE WITHDRAWN 1050 N m BOUND IS NOT A CONSTRAINT HERE. The hinge moment is an OUTPUT
reported with its full normalising frame inline (q, reference length, and model versus
aircraft scale), unconstrained until Twente supplies a traceable limit.
"""

X_H_OVER_C = 0.62
ETA_RANGE = (0.60, 1.00)

_FORBIDDEN = ("cmy", "sectional_cmy", "sectional_cmy_native", "Cmy", "CMY",
              "cmyo", "cmyi", "cmy*c/cref")


class HingeMomentRouteError(RuntimeError):
    """Raised when a forbidden route to the hinge moment is attempted."""


def refuse_cmy_route(source: str, column: str = None):
    """Call this at the top of ANY hinge-moment code path that touches delivered data.

    It raises. That is the whole point: the note version of this rule was already in
    the record and would not have stopped anyone."""
    raise HingeMomentRouteError(
        f"FORBIDDEN ROUTE TO THE HINGE MOMENT: source={source!r} column={column!r}.\n"
        f"sectional_cmy_native is flagged native-unverified by the author's own data "
        f"dictionary, and D053 could not close the delivered Cm from cmy under any "
        f"simple area weighting.\n"
        f"The ONLY sanctioned route is RANS pressure integration about "
        f"x_h/c = {X_H_OVER_C} over eta {ETA_RANGE[0]} to {ETA_RANGE[1]}. "
        f"No cmy, no transfer formula, no exception (project decision, 2026-07-29).")


def guard_columns(columns):
    """Raise if a caller hands us a table containing a forbidden moment column while
    computing a hinge moment. Use on every delivered-data read in this path."""
    hit = [c for c in columns if str(c).strip() in _FORBIDDEN]
    if hit:
        refuse_cmy_route(source="delivered spanwise loads", column=hit[0])
    return columns


def sanctioned_route_description():
    return (f"RANS surface-pressure integration of the local aerodynamic moment about "
            f"the chordwise line x_h/c = {X_H_OVER_C}, over {ETA_RANGE[0]} <= eta <= "
            f"{ETA_RANGE[1]}, reported as an OUTPUT with q, reference length and "
            f"model-versus-aircraft scale stated inline.")


if __name__ == "__main__":
    print(sanctioned_route_description())
    for bad in ("sectional_cmy_native", "cmy"):
        try:
            guard_columns(["eta", "chord_m", bad])
        except HingeMomentRouteError as e:
            print(f"\nREFUSAL FIRED on {bad!r}:\n{e}")
