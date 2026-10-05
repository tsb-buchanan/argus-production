#!/usr/bin/env python3
"""Where each concept moves the surface: signed displacement from the baseline along the
baseline's outward normal, upper and lower, one shared scale. Same computation as
scripts/render_concept_3d.py (nearest point through a KD-tree, 99.5th-percentile range);
drawn at slide size with large labels and the colour bar along the bottom."""
import gzip
import os
import re
import numpy as np
import deck_common as dc
from matplotlib.tri import Triangulation
from scipy.spatial import cKDTree

dc.style()
REPO = os.path.normpath(os.path.join(dc.R, "..", "..", ".."))
CONCEPTS = [("C", "cte_i002_c04"), ("F", "cfft_b02_c01"), ("H", "chc_g02_c06"),
            ("M", "mcv2_i002_c01"), ("W", "cffw_b01_c01")]
NUM = r"(-?\d+\.?\d*(?:[eE][-+]?\d+)?)"


def tris(name):
    txt = gzip.open(os.path.join(REPO, "geometry/stl/%s.stl.gz" % name), "rt", errors="replace").read()
    return np.array(re.findall(r"vertex\s+" + NUM + r"\s+" + NUM + r"\s+" + NUM, txt), dtype=float).reshape(-1, 3, 3)


def unique_points(T):
    P = T.reshape(-1, 3)
    _, idx, inv = np.unique(np.round(P, 7), axis=0, return_index=True, return_inverse=True)
    return P[idx], inv.reshape(-1, 3)


Tb = tris("baseline"); Pb, Fb = unique_points(Tb)
n = np.cross(Tb[:, 1] - Tb[:, 0], Tb[:, 2] - Tb[:, 0]); n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)
acc = np.zeros_like(Pb)
for k in range(3):
    np.add.at(acc, Fb[:, k], n)
Nb = acc / np.maximum(np.linalg.norm(acc, axis=1, keepdims=True), 1e-30)
tree = cKDTree(Pb)
data = []
for L, name in CONCEPTS:
    T = tris(name); P, F = unique_points(T)
    _, i = tree.query(P, k=1)
    signed = np.einsum("ij,ij->i", P - Pb[i], Nb[i]) * 1e3
    fn = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    data.append((L, P, F, signed, fn[:, 2] >= 0))
vmax = max(float(np.percentile(np.abs(d[3]), 99.5)) for d in data)
fig, axes = dc.plt.subplots(2, 5, figsize=(dc.WIDE_W, 4.5))
fig.set_layout_engine("none")
fig.subplots_adjust(left=0.035, right=0.995, top=0.995, bottom=0.2, wspace=0.04, hspace=0.0)
for col, (L, P, F, signed, up) in enumerate(data):
    for row, sel in enumerate((up, ~up)):
        ax = axes[row, col]
        im = ax.tripcolor(Triangulation(P[:, 1], P[:, 0], F[sel]), signed, cmap="coolwarm",
                          vmin=-vmax, vmax=vmax, shading="gouraud", rasterized=True)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        if row == 0:
            ax.text(0.03, 0.95, L, transform=ax.transAxes, fontsize=20, fontweight="bold",
                    color=dc.COLOURS[L], va="top")
        if col == 0:
            ax.set_ylabel(("upper", "lower")[row], fontsize=17)
cax = fig.add_axes([0.3, 0.125, 0.4, 0.035])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("displacement from B along the outward normal  (mm)", fontsize=16)
cb.ax.tick_params(labelsize=14)
dc.save(fig, "deck_concept", dc.WIDE_W, 4.5)
