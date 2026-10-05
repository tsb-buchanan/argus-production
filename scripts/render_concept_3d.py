#!/usr/bin/env python3
"""render_concept_3d.py: the wing in 3D, top and bottom, coloured by what each concept moves.

Review, 2026-08-22, of the section plots: they do not show the wing; 3D top and bottom
views of the surface were wanted.

Fair. Sections answer "how much and where along the chord"; they do not show the wing.
This writes a VTK per concept carrying the DISPLACEMENT FROM THE BASELINE as a point
field, then renders top and bottom with one shared camera and one shared colour scale.

WHY NEAREST-POINT AND NOT VERTEX INDEX. mcv2 and cte share the baseline's tessellation
(192,588 triangles) but the three newly ingested concepts do not (195,266), because the
delaunay tip cap adds facets. Matching by index would silently compare unrelated
vertices on exactly the three surfaces that matter most. A KD-tree nearest-point query
makes no assumption about the mesh and works for all five.

SIGNED, NOT ABSOLUTE. The sign is the physics: a trailing-edge droop moves the surface
DOWN, and |d| would hide that while making an upward twist at the leading edge look the
same as a downward one. The sign is taken along the baseline's own outward normal.

ONE COLOUR SCALE ACROSS EVERY PANEL. Per-panel autoscaling is the classic way to make
four different magnitudes look identical; the range is fixed from the largest concept so
the panels are comparable by eye, which is the only reason to put them side by side.
"""
import argparse, gzip, re, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
# PANELS CARRY THE CAMPAIGN TAG LETTER, not just the concept description, because the
# results tables name geometries as SWC / CMPC / LCC and a reader cannot otherwise connect
# this figure to them. The letter is the suffix shared by all three conditions (SW =
# Condition CR, CMP = early cruise, LC = late cruise).
#
# THE MAPPING IS TAKEN FROM EACH CASE'S OWN CASE_PROVENANCE.json, i.e. from the STL the
# mesh was actually built on, and cross-checked against the registry sha256. It is NOT
# taken from any prose document, because two of them disagree:
# docs/report/argus_geometry_report.tex says "SWC & cfft_b02_c01", which is wrong;
# docs/CAMPAIGN_V5_PLAN.md says SWC = cte, SWF = cfft, which the provenance confirms.
# Witness: SWF's wing.stl sha256 begins 82bf2e1e570187ab, the registry hash for
# cfft_b02_c01_oml_placed_laddercap.stl (registry/candidates.yaml:190).
#
# ORDERED BY TAG LETTER, so the panels run in the same order as the geometry tables in the
# report. The previous order was the ingestion order, which matches nothing a reader has.
# The shared colour range is computed over all five regardless, so order does not move it.
CONCEPTS = [("cte_i002_c04",  "C: TE camber (retained)"),
            ("cfft_b02_c01",  "F: TE camber A/c 0.0139"),
            ("chc_g02_c06",   "H: hinged"),
            ("mcv2_i002_c01", "M: TE camber A/c 0.0337"),
            ("cffw_b01_c01",  "W: TWIST")]


def tris(name):
    txt = gzip.open(REPO / ("geometry/stl/%s.stl.gz" % name), "rt", errors="replace").read()
    return np.array(re.findall(r"vertex\s+(-?\d+\.?\d*(?:[eE][-+]?\d+)?)\s+"
                               r"(-?\d+\.?\d*(?:[eE][-+]?\d+)?)\s+"
                               r"(-?\d+\.?\d*(?:[eE][-+]?\d+)?)", txt),
                    dtype=float).reshape(-1, 3, 3)


def unique_points(T):
    P = T.reshape(-1, 3)
    key = np.round(P, 7)
    _, idx, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    return P[idx], inv.reshape(-1, 3)


def normals_at_points(T, pts, faces):
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.where(ln > 0, n / np.maximum(ln, 1e-30), 0.0)
    acc = np.zeros_like(pts)
    for k in range(3):
        np.add.at(acc, faces[:, k], n)
    ln = np.linalg.norm(acc, axis=1, keepdims=True)
    return np.where(ln > 0, acc / np.maximum(ln, 1e-30), 0.0)


def write_vtk(path, pts, faces, field, fname):
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\nARGUS concept displacement\nASCII\n"
                "DATASET POLYDATA\nPOINTS %d float\n" % len(pts))
        for p in pts:
            f.write("%.6f %.6f %.6f\n" % tuple(p))
        f.write("POLYGONS %d %d\n" % (len(faces), 4 * len(faces)))
        for t in faces:
            f.write("3 %d %d %d\n" % tuple(t))
        f.write("POINT_DATA %d\nSCALARS %s float 1\nLOOKUP_TABLE default\n"
                % (len(pts), fname))
        for v in field:
            f.write("%.6f\n" % v)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="results/figures/deck")
    ap.add_argument("--size", nargs=2, type=int, default=[1500, 700])
    # FOR PRINT, THE IN-IMAGE CAPTION IS WORSE THAN NONE. Placed in a report at page width
    # it renders as four lines of unreadable grey type directly above the LaTeX caption
    # that says the same thing. On a slide, where there is no caption, it is the only place
    # the frame can live, so it stays the default and is suppressed explicitly.
    ap.add_argument("--no-caption", action="store_true",
                    help="omit the in-image caption; use when a document supplies one")
    a = ap.parse_args()
    from scipy.spatial import cKDTree

    Tb = tris("baseline")
    Pb, Fb = unique_points(Tb)
    Nb = normals_at_points(Tb, Pb, Fb)
    tree = cKDTree(Pb)

    work = Path(tempfile.mkdtemp(prefix="argus_c3d_"))
    made, vmax = [], 0.0
    for name, lab in CONCEPTS:
        T = tris(name)
        P, F = unique_points(T)
        d, i = tree.query(P, k=1)
        signed = np.einsum("ij,ij->i", P - Pb[i], Nb[i]) * 1e3      # mm along the normal
        vmax = max(vmax, float(np.percentile(np.abs(signed), 99.5)))
        out = work / ("%s.vtk" % name)
        write_vtk(out, P, F, signed, "displacement_mm")
        made.append((name, lab, out))
        print("  %-16s %7d pts  displacement %+.2f .. %+.2f mm"
              % (name, len(P), signed.min(), signed.max()))
    print("  shared colour range: +/- %.2f mm (99.5th percentile)" % vmax)

    # ---- RENDER WITH MATPLOTLIB, NOT PARAVIEW ---------------------------------
    # pvbatch aborts here with "xcb_xlib_threads_sequence_lost" even under
    # --force-offscreen-rendering and a sequential VTK SMP backend. It has cost this
    # project time twice. A PLANFORM view is also the better picture for this question:
    # the wing is thin, so a perspective 3D view hides the aft region behind foreshorten-
    # ing, while looking straight down shows exactly where on the wing the surface moves.
    # Upper and lower are separated by the sign of the facet normal, which is exact.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    outdir = REPO / a.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, len(made), figsize=(3.2 * len(made), 7.4),
                             constrained_layout=True)
    im = None
    for j, (name, lab, _) in enumerate(made):
        T = tris(name)
        P, F = unique_points(T)
        d, i2 = tree.query(P, k=1)
        signed = np.einsum("ij,ij->i", P - Pb[i2], Nb[i2]) * 1e3
        fn = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
        upper_face = fn[:, 2] >= 0
        for row, sel, title in ((0, upper_face, "upper"), (1, ~upper_face, "lower")):
            ax = axes[row, j]
            faces = F[sel]
            # SPAN HORIZONTAL, CHORD VERTICAL. Plotting (x, y) directly puts span
            # up the page, which is not how a planform is read and which my first
            # caption then described wrongly as horizontal. Swapping the axes here
            # makes the picture and the words agree.
            tri = Triangulation(P[:, 1], P[:, 0], faces)
            im = ax.tripcolor(tri, signed, cmap="coolwarm", vmin=-vmax, vmax=vmax,
                              shading="gouraud", rasterized=True)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            if row == 0:
                ax.set_title(lab, fontsize=10)
            if j == 0:
                ax.set_ylabel(title.upper(), fontsize=10)
    cb = fig.colorbar(im, ax=axes, shrink=0.55, pad=0.01)
    cb.set_label("displacement from baseline, along the outward normal [mm]", fontsize=9)
    fig.suptitle("Each concept on the wing: where the surface actually moves, "
                 "seen from above and below", fontsize=13)
    # WRAPPED EXPLICITLY. As one run-on string this overran the right page edge mid-word,
    # so the last clause was unreadable in the report. Matplotlib does not wrap figure text
    # by default and constrained_layout does not police it, so the breaks are authored here.
    if not a.no_caption:
     fig.text(0.005, 0.004,
             "Panel letters are the campaign geometry suffix: C is SWC / CMPC / LCC, and so on, "
             "where SW = Condition CR, CMP = early cruise, LC = late cruise.\n"
             "Planform view: SPAN runs left to right (root at left), CHORD up the page, flow from "
             "top to bottom. Signed along the BASELINE's outward normal, so the same\n"
             "colour means the same physical direction on both surfaces. One scale across every "
             "panel (+/- %.2f mm, 99.5th percentile) because per-panel autoscaling\n"
             "makes different magnitudes look identical. Matched by nearest point, not vertex "
             "index: the three newly ingested concepts carry 195,266 triangles\n"
             "against the baseline's 192,588." % vmax, fontsize=7, va="bottom")
    out = outdir / "concept_displacement_3d.png"
    fig.savefig(out, dpi=150)
    # relative_to RAISES on any --outdir outside the repo, which is the normal case when
    # rendering to a scratch directory before promoting. It fired AFTER savefig, so the
    # figure was written and the script still exited non-zero: a success reported as a
    # crash, which is the inverse of the vouched-success defect and just as misleading to
    # anything reading the exit status.
    try:
        shown = out.relative_to(REPO)
    except ValueError:
        shown = out
    print("  wrote %s" % shown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
