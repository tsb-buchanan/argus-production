"""Shared loader for the binned planform maps (results/wing_maps_of12) and the planform
mesh they are painted on, copied from scripts/render_delta_maps_of12.py so the deck draws
the deltas on the same support the report does."""
import json
import os
import numpy as np
import deck_common as dc

MAPS = os.path.normpath(os.path.join(dc.R, "..", "..", "..", "results", "wing_maps_of12"))


def load_map(case):
    npz = os.path.join(MAPS, "maps_%s.npz" % case)
    js = os.path.join(MAPS, "maps_%s.json" % case)
    if not (os.path.exists(npz) and os.path.exists(js)):
        raise SystemExit("HALT: no binned map for %s at %s" % (case, npz))
    return dict(np.load(npz)), json.load(open(js))


def planform_mesh(M):
    """Bin-centre coordinates regularised to a monotone (eta, x/c) mesh, as the report does."""
    X, Y = np.array(M["x_m"], float), np.array(M["y_m"], float)
    n_eta, n_xc = X.shape
    xcc = (np.arange(n_xc) + 0.5) / n_xc
    sl = np.full(n_eta, np.nan); ic = np.full(n_eta, np.nan); yv = np.full(n_eta, np.nan)
    for i in range(n_eta):
        g = np.isfinite(X[i])
        if g.sum() >= 2:
            sl[i], ic[i] = np.polyfit(xcc[g], X[i][g], 1)
        gy = np.isfinite(Y[i])
        if gy.any():
            yv[i] = np.median(Y[i][gy])
    idx = np.arange(n_eta)
    for a in (sl, ic, yv):
        f = np.isfinite(a)
        if not f.any():
            raise SystemExit("HALT: planform mesh has no finite row")
        a[~f] = np.interp(idx[~f], idx[f], a[f])
    Xf = ic[:, None] + sl[:, None] * xcc[None, :]
    Yf = np.repeat(yv[:, None], n_xc, axis=1)
    return Xf, Yf
