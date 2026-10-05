#!/usr/bin/env python3
"""Shared loaders, colours and sizing for the group-deck figures (deck_*.py).

EVERY FIGURE IS BUILT AT THE SIZE IT IS PLACED AT, in inches, and saved at 200 dpi WITHOUT
bbox 'tight', so the PNG's pixel size is exactly (w*200, h*200) and every font lands on the
slide at the point size it was drawn at (brief rule 3: axis labels 16 pt or larger). save()
asserts the pixel size after writing.

TEXT ON A FIGURE IS RESTRICTED to axis labels, tick labels, legend and numeric data labels
(brief rule 4). No titles, no notes. Panel identity, where a figure has to carry two panels,
is expressed through the axis label.

COLOURS are the report's Okabe-Ito map from make_cp.py, one per geometry, on every slide.
"""
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.dirname(HERE)
DATA = os.path.join(R, "data")
OUT = os.path.join(HERE, "deck")
DPI = 200

COLOURS = {"B": "#000000", "C": "#D55E00", "F": "#0072B2",
           "H": "#009E73", "M": "#CC79A7", "W": "#E69F00"}
ORDER = ["B", "C", "F", "H", "M", "W"]
CAND = ["C", "F", "H", "M", "W"]
# TU Delft theme colours, read off the template's theme1.xml
TUD_CYAN, TUD_DARK, TUD_GREY = "#00A6D6", "#017188", "#A7A7A7"

# Placed sizes, inches. RESULT is the right-hand figure box of a result slide (rule 2: at
# least half the 13.333 in slide width, full content height).
RESULT_W, RESULT_H = 6.7, 5.15
HALF_W, HALF_H = 6.0, 4.5   # two figures side by side under two bullets
NARROW_W = 3.25         # two rank panels side by side inside the result box
WIDE_W, WIDE_H = 12.1, 4.55

CONDITIONS = {"condition_CR": "Condition CR", "early_cruise": "early cruise",
              "late_cruise": "late cruise"}


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 15, "axes.labelsize": 17, "xtick.labelsize": 15,
        "ytick.labelsize": 15, "legend.fontsize": 14, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25,
        "grid.linewidth": 0.6, "legend.frameon": False, "figure.constrained_layout.use": True,
    })


def jload(rel):
    with open(os.path.join(R, rel)) as fh:
        return json.load(fh)


def cload(rel):
    with open(os.path.join(R, rel)) as fh:
        return list(csv.DictReader(fh))


def forces():
    return jload("data/rans_forces.json")


def trims(cond):
    """letter -> trim entry for one condition, straight from rans_forces.json."""
    return forces()["trims"][cond]["geometries"]


def vlmerr():
    return jload("fig/fig_vlmerr_values.json")


def spaneff():
    return jload("data/span_efficiency.json")


def vlm():
    return jload("data/vlm.json")


def save(fig, name, w, h):
    """Write fig/deck/<name>.png at exactly (w, h) inches x 200 dpi and assert it."""
    os.makedirs(OUT, exist_ok=True)
    fig.set_size_inches(w, h)
    path = os.path.join(OUT, name + ".png")
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    with open(path, "rb") as fh:
        head = fh.read(24)
    pw, ph = int.from_bytes(head[16:20], "big"), int.from_bytes(head[20:24], "big")
    want = (round(w * DPI), round(h * DPI))
    if (pw, ph) != want:
        raise SystemExit("HALT: %s is %dx%d px, expected %dx%d" % (path, pw, ph, *want))
    print("wrote %s  %dx%d px  (%.2f x %.2f in)" % (os.path.relpath(path, R), pw, ph, w, h))
    return path


def bar_labels(ax, bars, fmt="%+.3f", size=14, pad=3):
    """Numeric data labels on bars, above positive bars and below negative ones."""
    for b in bars:
        v = b.get_height()
        ax.annotate(fmt % v, (b.get_x() + b.get_width() / 2, v),
                    xytext=(0, pad if v >= 0 else -pad), textcoords="offset points",
                    ha="center", va="bottom" if v >= 0 else "top", fontsize=size)


def read_section_csv(path):
    """x/c, cp, cf_s, is_upper from a section CSV, columns by name (see make_cp.py)."""
    with open(path) as fh:
        rows = [r for r in csv.reader(fh) if r and not r[0].startswith("#")]
    col = {h: i for i, h in enumerate(rows[0])}
    body = rows[1:]
    x = np.array([float(r[col["x_over_c"]]) for r in body])
    cp = np.array([float(r[col["cp"]]) for r in body])
    cfs = np.array([float(r[col["cf_s"]]) for r in body])
    up = np.array([r[col["surface"]] == "upper" for r in body])
    return x, cp, cfs, up


def binned_median(x, y, n=200, min_faces=5):
    """Median of y in n equal x/c bins, NaN where fewer than min_faces (as make_cp/make_cf)."""
    edges = np.linspace(0.0, 1.0, n + 1)
    idx = np.digitize(x, edges) - 1
    xc = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(n, np.nan)
    for k in range(n):
        m = idx == k
        if m.sum() >= min_faces:
            med[k] = np.median(y[m])
    return xc, med


def section_file(case, eta_int):
    """Resolve the section CSV for a case through the report's own rule (section_source)."""
    import sys
    sys.path.insert(0, HERE)
    from section_source import resolve
    d, _, gen = resolve(case)
    return os.path.join(d, "%s_eta%03d.csv" % (case, eta_int)), gen
