#!/usr/bin/env python3
"""Title slide: the baseline wing seen from above, coloured by c_p, nothing else. Drawn from
the binned planform map results/wing_maps_of12/maps_SWB_trim.npz (the same bins the report's
surface maps use, with their physical bin coordinates), painted on the planform mesh the way
render_delta_maps_of12.py does. Sized to the template's title-page picture placeholder."""
import numpy as np
import deck_common as dc
from deck_maps_common import planform_mesh, load_map

dc.style()
M, prov = load_map("SWB_trim")
X, Y = planform_mesh(M)
fin = M["cp_upper"][np.isfinite(M["cp_upper"])]
lo, hi = np.percentile(fin, [0.5, 99.5])
fig = dc.plt.figure()
fig.set_layout_engine("none")
ax = fig.add_axes([0.02, 0.16, 0.96, 0.83])
pm = ax.pcolormesh(X, Y, np.ma.masked_invalid(M["cp_upper"]), cmap="coolwarm_r", vmin=lo, vmax=hi,
                   shading="gouraud")
ax.set_aspect("equal"); ax.axis("off")
cax = fig.add_axes([0.25, 0.105, 0.5, 0.03])
cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
cb.set_label(r"$c_p$, upper surface, Condition CR", fontsize=18)
cb.ax.tick_params(labelsize=15)
dc.save(fig, "deck_title_cp", 8.0, 7.5)
