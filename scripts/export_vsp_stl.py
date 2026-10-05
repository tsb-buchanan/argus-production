#!/usr/bin/env python3
"""Export a high-tessellation STL from a .vsp3 using OpenVSP's own tessellator.

Runs INSIDE the OpenVSP 3.51.2 container:

    OPENVSP_SCRATCH=<dir> scripts/bin/vsppython scripts/export_vsp_stl.py \
        --vsp3 geometry/source/baseline_wing_only_refined.vsp3 \
        --outdir <dir> --name baseline

This is the production geometry route. It supersedes scripts/vsp3_to_stl.py,
which walked the wing geom's own parms and therefore never saw the placement
transform carried by the PARENT Blank geom 'wing_comps' (X_Rel_Rotation = 5 deg
dihedral, Y_Rel_Rotation = 1.873 deg incidence). See docs/decisions.md D048.

UNITS: OpenVSP writes STL in model units, which for this model are FEET. Both
artifacts are emitted and checksummed: the raw export in feet (provenance) and
a scaled copy in metres (x 0.3048) for meshing. Nothing is left implicit.
"""
import argparse
import hashlib
import json
import os
import struct
import sys

import openvsp as vsp

FT_TO_M = 0.3048


def drain_errors():
    """Return and clear queued VSP errors.

    ErrorMgr.PopErrorAndPrint wants a C FILE*, which the Python binding cannot
    supply from sys.stdout, so errors are popped individually instead.
    """
    errs = vsp.ErrorMgrSingleton.getInstance()
    out = []
    while errs.GetNumTotalErrors() > 0:
        out.append(errs.PopLastError().GetErrorString())
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_stl(path):
    """Return (triangles, is_binary). triangles is a list of 3x3 float tuples."""
    with open(path, "rb") as fh:
        head = fh.read(5)
        fh.seek(0)
        if head[:5] == b"solid":
            tris, cur = [], []
            for line in fh:
                s = line.strip().split()
                if s and s[0] == b"vertex":
                    cur.append(tuple(float(v) for v in s[1:4]))
                    if len(cur) == 3:
                        tris.append(cur)
                        cur = []
            return tris, False
        fh.seek(80)
        n = struct.unpack("<I", fh.read(4))[0]
        tris = []
        for _ in range(n):
            d = struct.unpack("<12fH", fh.read(50))
            tris.append([d[3:6], d[6:9], d[9:12]])
        return tris, True


def write_stl_ascii(path, tris, name="wing"):
    with open(path, "w") as fh:
        fh.write("solid %s\n" % name)
        for t in tris:
            fh.write("  facet normal 0 0 0\n    outer loop\n")
            for v in t:
                fh.write("      vertex %.9e %.9e %.9e\n" % tuple(v))
            fh.write("    endloop\n  endfacet\n")
        fh.write("endsolid %s\n" % name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vsp3", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--tess-w", type=int, default=401,
                    help="points around the airfoil; 401 gives ~200 per surface")
    ap.add_argument("--tess-u-total", type=int, default=240,
                    help="target spanwise stations summed over all sections")
    args = ap.parse_args()

    vsp.ClearVSPModel()
    vsp.ReadVSPFile(args.vsp3)

    geoms = vsp.FindGeoms()
    wings = [g for g in geoms if vsp.GetGeomTypeName(g) == "Wing"]
    if len(wings) != 1:
        sys.exit("expected exactly 1 Wing geom, found %d" % len(wings))
    wing = wings[0]

    # Record the full parent transform chain. This is the data the in-house
    # loft dropped, so it is provenance, not decoration.
    drain_errors()
    chain = []
    node = wing
    # GetGeomParent returns the string "NONE" (not "") at the root, so testing
    # truthiness alone walks one step too far and raises spurious errors.
    while node and node != "NONE":
        entry = {"id": node, "name": vsp.GetGeomName(node),
                 "type": vsp.GetGeomTypeName(node), "parms": {}}
        for pname in ("X_Rel_Rotation", "Y_Rel_Rotation", "Z_Rel_Rotation",
                      "X_Rel_Location", "Y_Rel_Location", "Z_Rel_Location"):
            pid = vsp.GetParm(node, pname, "XForm")
            if pid:
                entry["parms"][pname] = vsp.GetParmVal(pid)
        chain.append(entry)
        node = vsp.GetGeomParent(node)

    # Tessellation. Chordwise via Tess_W on the geom; spanwise via SectTess_U
    # per section, distributed in proportion to section span so that panel
    # density is uniform rather than per-panel-constant.
    pid_w = vsp.GetParm(wing, "Tess_W", "Shape")
    vsp.SetParmVal(pid_w, args.tess_w)
    vsp.Update()
    tess_w_set = vsp.GetParmVal(pid_w)

    xsurf = vsp.GetXSecSurf(wing, 0)
    nsec = vsp.GetNumXSec(xsurf)
    spans = []
    for i in range(1, nsec):
        pid = vsp.FindParm(wing, "Span", "XSec_%d" % i)
        spans.append(vsp.GetParmVal(pid) if pid else 0.0)
    total = sum(spans)

    sect_tess = []
    for i in range(1, nsec):
        want = max(2, int(round(args.tess_u_total * spans[i - 1] / total)))
        pid = vsp.FindParm(wing, "SectTess_U", "XSec_%d" % i)
        if pid:
            vsp.SetParmVal(pid, want)
        sect_tess.append(want)
    vsp.Update()

    os.makedirs(args.outdir, exist_ok=True)
    ft_path = os.path.join(args.outdir, "%s_vsp_ft.stl" % args.name)
    m_path = os.path.join(args.outdir, "%s_vsp_m.stl" % args.name)

    drain_errors()
    vsp.ExportFile(ft_path, vsp.SET_ALL, vsp.EXPORT_STL)
    export_errors = drain_errors()
    if export_errors:
        print("VSP reported errors during export:", file=sys.stderr)
        for e in export_errors:
            print("   ", e, file=sys.stderr)

    tris, was_binary = read_stl(ft_path)
    tris_m = [[tuple(c * FT_TO_M for c in v) for v in t] for t in tris]
    write_stl_ascii(m_path, tris_m, name=args.name)

    meta = {
        "source_vsp3": args.vsp3,
        "source_vsp3_sha256": sha256(args.vsp3),
        "exporter": "OpenVSP %s ExportFile(EXPORT_STL) via python API" % vsp.GetVSPVersion(),
        # OPENVSP_IMAGE is not forwarded into the container, so a default here
        # would record a GUESSED image name as provenance. Record it only when
        # it is actually known, and say so when it is not.
        "container_image": os.environ.get("OPENVSP_IMAGE", "UNRECORDED (not passed into container)"),
        "geom": vsp.GetGeomName(wing),
        "transform_chain_root_last": chain,
        "tessellation": {
            "Tess_W": tess_w_set,
            "n_chordwise_per_surface_approx": int((tess_w_set - 1) / 2),
            "SectTess_U_per_section": sect_tess,
            "sum_SectTess_U": sum(sect_tess),
            "n_sections": nsec - 1,
        },
        "units": {
            "native": "feet",
            "raw_export": "feet (unscaled, as written by OpenVSP)",
            "scale_applied_to_metre_copy": FT_TO_M,
        },
        "triangles": len(tris),
        "export_errors": export_errors,
        "files": {
            "feet": {"path": ft_path, "sha256": sha256(ft_path),
                     "format": "binary" if was_binary else "ascii"},
            "metres": {"path": m_path, "sha256": sha256(m_path), "format": "ascii"},
        },
    }
    with open(os.path.join(args.outdir, "%s_vsp_export.json" % args.name), "w") as fh:
        json.dump(meta, fh, indent=2)

    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
