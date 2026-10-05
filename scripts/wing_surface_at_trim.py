#!/usr/bin/env python3
"""wing_surface_at_trim.py: one wing surface AT the trim condition, for experiment comparison.

    python3 scripts/wing_surface_at_trim.py

WHY THIS EXISTS. The surface an experimentalist compares against pressure taps has to be at
the SAME LIFT COEFFICIENT as the experiment, and no computed case is: every case runs at a
fixed angle and lands where it lands. The surfaces previously delivered were worse than
that. They were `KR_wake` at C_L 0.455054 and `MR_mcv2` at C_L 0.492932, which are +268 and
+647 counts from the target, on the superseded OpenFOAM-7 chain, at 188k cells, and
carrying skin friction as a MAGNITUDE so no direction and no separation line can be taken
from them.

WHAT THIS PRODUCES. The three converged OpenFOAM-12 angles interpolated field by field to
the trim C_L, giving a single surface at C_L 0.428277635108 exactly, 5,276,769 faces, with
wall shear as a VECTOR.

THE PRECONDITION, ASSERTED RATHER THAN ASSUMED. Interpolating field by field is only valid
if the three surfaces are the SAME MESH in the same order. Point count, polygon count and
the polygon connectivity block are compared byte for byte across the three cases, and the
run HALTS if they differ. A silent mismatch here would blend fields from different cells
and produce a plausible, wrong surface.

INTERPOLATION IS AGAINST C_L, NOT ALPHA, matching every other quantity in this delivery:
the target is defined on C_L, so interpolating in C_L lands on it exactly.

OUTPUT IS BINARY LEGACY VTK, which is about half the size of the ASCII the solver writes
and is read by ParaView, VisIt and Tecplot without conversion. Legacy VTK binary is
BIG-ENDIAN by specification regardless of host, which is the one thing implementations get
wrong; numpy's byteswap is applied explicitly rather than relied on.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
TARGET_CL = 0.428277635108
CASES = ["SWB_a1p60", "SWB_a1p75", "SWB_a1p90"]


def surface_path(case):
    hits = sorted((REPO / "results/postpro_latest" / case /
                   "surfaces/conv/argusWingSurface").glob("*/wingSurface.vtk"))
    if not hits:
        sys.exit("no converged surface for %s" % case)
    return hits[-1]


def parse(path):
    """Return (header_bytes_through_polygons, npoints, ncells, {field: array})."""
    raw = path.read_bytes()
    txt = raw.decode("ascii", errors="replace")

    mpts = re.search(r"POINTS (\d+) float\n", txt)
    mpol = re.search(r"POLYGONS (\d+) (\d+)\n", txt)
    mcd = re.search(r"CELL_DATA (\d+)\n", txt)
    if not (mpts and mpol and mcd):
        sys.exit("%s: not the expected sampleSurface layout" % path)
    npts, ncells = int(mpts.group(1)), int(mcd.group(1))

    # geometry block: from the start through the end of POLYGONS, verbatim
    geom = raw[:mcd.start()]

    fields = {}
    pos = mcd.end()
    mfield = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*) (\d+) (\d+) float\n", re.M)
    for m in mfield.finditer(txt, pos):
        name, ncomp, n = m.group(1), int(m.group(2)), int(m.group(3))
        if n != ncells:
            continue
        nxt = mfield.search(txt, m.end())
        block = txt[m.end(): nxt.start() if nxt else len(txt)]
        a = np.fromstring(block, sep=" ", dtype=np.float64)
        want = ncells * ncomp
        if a.size < want:
            sys.exit("%s: field %s short, %d of %d" % (path, name, a.size, want))
        fields[name] = a[:want].reshape(ncells, ncomp) if ncomp > 1 else a[:want]
    return geom, npts, ncells, fields


def main():
    at = json.loads((REPO / "results/derived/at_target_cl.json").read_text())
    cl = {p["case"]: p["cl"] for p in at["points"]}

    geom0 = npts0 = ncells0 = None
    data, xs = {}, []
    for c in CASES:
        p = surface_path(c)
        print("  reading %-12s %s" % (c, p.relative_to(REPO)))
        geom, npts, ncells, f = parse(p)
        if geom0 is None:
            geom0, npts0, ncells0 = geom, npts, ncells
        else:
            # THE PRECONDITION. Same mesh, same order, byte for byte.
            if npts != npts0 or ncells != ncells0 or geom != geom0:
                sys.exit("HALTED: %s is not the same mesh as %s. Field-by-field "
                         "interpolation would blend different cells." % (c, CASES[0]))
        data[c] = f
        xs.append(cl[c])
        print("    %d points, %d cells, fields: %s" % (npts, ncells, ", ".join(f)))

    print("  mesh identity asserted across all three: %d points, %d cells"
          % (npts0, ncells0))

    names = [n for n in data[CASES[0]] if all(n in data[c] for c in CASES)]
    xs = np.asarray(xs)
    order = np.argsort(xs)
    out = {}
    print("  interpolating to the trim C_L %.12f (bracket %.6f to %.6f)"
          % (TARGET_CL, xs.min(), xs.max()))
    for n in names:
        stack = np.stack([data[c][n] for c in CASES])[order]
        A = np.vstack([xs[order], np.ones(len(xs))]).T
        coef, *_ = np.linalg.lstsq(A, stack.reshape(len(xs), -1), rcond=None)
        out[n] = (coef[0] * TARGET_CL + coef[1]).reshape(stack.shape[1:])
        print("    %-18s %s" % (n, "x".join(str(s) for s in out[n].shape)))

    dst = REPO / "results/vtk_export/of12_SWB_wing_surface_at_trim.vtk"
    dst.parent.mkdir(parents=True, exist_ok=True)
    hdr = (b"# vtk DataFile Version 3.0\n"
           b"ARGUS baseline wing AT TRIM: C_L 0.428277635108, alpha 1.7295 deg, "
           b"Condition CR M 0.10, OpenFOAM-12. Interpolated against C_L across "
           b"SWB_a1p60/75/90 on an identical mesh. wallShearStress is a VECTOR and is "
           b"KINEMATIC (divide-free of rho); p is kinematic. cp is dimensionless.\n"
           b"BINARY\n")
    body = geom0.split(b"\n", 3)[3]        # drop the source's own 3 header lines
    with open(dst, "wb") as fh:
        fh.write(hdr)
        fh.write(body)
        fh.write(b"CELL_DATA %d\n" % ncells0)
        fh.write(b"FIELD attributes %d\n" % len(names))
        for n in names:
            a = np.asarray(out[n], dtype=">f4")
            ncomp = 1 if a.ndim == 1 else a.shape[1]
            fh.write(b"%s %d %d float\n" % (n.encode(), ncomp, ncells0))
            fh.write(a.tobytes())
            fh.write(b"\n")
    print("  wrote %s  (%.1f MB)" % (dst.relative_to(REPO), dst.stat().st_size / 1048576))

    (REPO / "results/derived/wing_surface_at_trim.json").write_text(json.dumps(dict(
        target_cl=TARGET_CL, cases=CASES, case_cl=cl, n_points=npts0, n_cells=ncells0,
        fields=names, output=str(dst.relative_to(REPO)),
        note=("Interpolated against C_L on an identical mesh, asserted byte for byte. "
              "Supersedes the OF-7 KR_wake / MR_mcv2 surfaces, which were at C_L 0.455054 "
              "and 0.492932, +268 and +647 counts from the target, and carried skin "
              "friction as a magnitude.")), indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
