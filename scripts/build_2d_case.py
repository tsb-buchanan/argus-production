#!/usr/bin/env python3
"""build_2d_case.py: registry-driven 2D O-grid case builder (ARGUS sections).

Builds cases/<candidate>/<COND>2d_<tag>_<level>_a<alpha>/ from
template_case/airfoil2d numerics dicts (verbatim) + a D015 O-grid on a
REGISTERED geometry. Inherits the three Layer-1 lessons (D012 verdict):
1. 500c farfield for lifting cases (uncorrected farfield error ~ cl^2).
2. Fully-turbulent enforcement: turbulent-branch k init (internal Tu ~2%),
   farfield BC keeps the condition's freestream spec.
3. Wall spacing y+ ~0.5 at the fine level, refining with the family
   (t1 8.8e-6 / 6.2e-6 / 4.4e-6 at Condition CR; y+ ~1.0/0.7/0.5) - D018.1;
   500c farfield default D018.2; turbulent-branch k init D018.3.

Conditions per docs/ARGUS_reference_data.md s9 (D006): CR: U 34.0, WT: U 40.8;
nu 1.46e-5; turbulence per template comments (CR k 1.73e-3, omega 22).

Usage:
  python3 scripts/build_2d_case.py --geometry baseline:geometry_2d_faired \
      --condition CR --level fine --alpha 2.0 [--np 8] [--tag grid]
"""

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
CHORD = 0.3404
LEVELS = {"coarse": (320, 1.2, 8.8e-6), "medium": (452, 1.137, 6.2e-6),
          "fine": (640, 1.0955, 4.4e-6)}
CONDITIONS = {
    "CR": {"U": 34.0, "nu": 1.46e-5, "k": 1.73e-3, "omega": 22.0,
           "mach": 0.10, "Re": 34.0 * CHORD / 1.46e-5},
    "WT": {"U": 40.8, "nu": 1.46e-5, "k": 2.50e-3, "omega": 27.0,
           "mach": 0.12, "Re": 40.8 * CHORD / 1.46e-5},
    # ---- NASA CR-2215 benchmark conditions (D012 Layer 2). ----
    # CHORD IS 0.6096 m HERE, NOT 0.3404, AND IT IS NOT COSMETIC. Re and M cannot both
    # be matched on a different chord: holding Re 2.83e6 on the 0.3404 m chord forces
    # U = 121 m/s and M = 0.36, which breaks the incompressible assumption the solver
    # rests on. So these conditions carry their own chord and --chord must be passed.
    # k and omega are scaled from CR: k as U^2, omega as U/c (dimensionally consistent).
    "C2215A06": {"U": 67.78, "nu": 1.46e-5, "k": 6.875e-3, "omega": 24.48,
                 "mach": 0.199, "Re": 2.83e6, "chord": 0.6096,
                 "note": "M 0.201 Re 2.83e6, transition FIXED 0.10c both surfaces (#80 "
                         "grit). The like-for-like case: our SST is fully turbulent, so a "
                         "tripped experiment is the correct comparison. Wake-rake drag. "
                         "CAVEAT: M 0.2 is borderline incompressible, PG scaling advised."},
    "C2215A05": {"U": 26.34, "nu": 1.46e-5, "k": 1.038e-3, "omega": 9.51,
                 "mach": 0.0775, "Re": 1.10e6, "chord": 0.6096,
                 "note": "M 0.080 Re 1.10e6, FREE transition. Closest condition to ARGUS "
                         "CR (M 0.10, Re 7.9e5). CAVEAT: free transition against a fully "
                         "turbulent SST solve -- expect the solve to over-read cd, the "
                         "same laminar-run gap seen against XFOIL ncrit 9."},
}

FIELDS = {
"U": """FoamFile {{ version 2.0; format ascii; class volVectorField; object U; }}
dimensions [0 1 -1 0 0 0 0];
internalField uniform ({ux} 0 {uz});
boundaryField
{{
    farfield {{ type freestream; freestreamValue uniform ({ux} 0 {uz}); }}
    front    {{ type empty; }}
    back     {{ type empty; }}
    airfoil  {{ type noSlip; }}
}}
""",
"p": """FoamFile {{ version 2.0; format ascii; class volScalarField; object p; }}
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{{
    farfield {{ type freestreamPressure; freestreamValue uniform 0; }}
    front    {{ type empty; }}
    back     {{ type empty; }}
    airfoil  {{ type zeroGradient; }}
}}
""",
"k": """/* Freestream k per condition (D006/template); internal init 1.0 =
   turbulent-branch selection (D012 Layer-1 lesson 2). */
FoamFile {{ version 2.0; format ascii; class volScalarField; object k; }}
dimensions [0 2 -2 0 0 0 0];
internalField uniform 1.0;
boundaryField
{{
    farfield {{ type inletOutlet; inletValue uniform {k}; value uniform {k}; }}
    front    {{ type empty; }}
    back     {{ type empty; }}
    airfoil  {{ type fixedValue; value uniform 1e-12; }}
}}
""",
"omega": """FoamFile {{ version 2.0; format ascii; class volScalarField; object omega; }}
dimensions [0 0 -1 0 0 0 0];
internalField uniform {omega};
boundaryField
{{
    farfield {{ type inletOutlet; inletValue uniform {omega}; value uniform {omega}; }}
    front    {{ type empty; }}
    back     {{ type empty; }}
    airfoil  {{ type omegaWallFunction; value uniform {omega}; }}
}}
""",
"nut": """FoamFile {{ version 2.0; format ascii; class volScalarField; object nut; }}
dimensions [0 2 -1 0 0 0 0];
internalField uniform 1e-8;
boundaryField
{{
    farfield {{ type calculated; value uniform 1e-8; }}
    front    {{ type empty; }}
    back     {{ type empty; }}
    airfoil  {{ type nutLowReWallFunction; value uniform 0; }}
}}
""",
}

ALLRUN = """#!/bin/sh
cd "${0%/*}" || exit
set -e
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
case "$(mpirun --version 2>/dev/null | head -1)" in
    *"Open MPI"*) ;;
    *) PATH=/usr/bin:$PATH; export PATH ;;
esac
runApplication renumberMesh -overwrite
runApplication checkMesh
[ -d 0 ] || cp -r 0.orig 0
[ -f warmstart_meta.json ] || runApplication potentialFoam
NP=$(grep -oE 'numberOfSubdomains [0-9]+' system/decomposeParDict | awk '{print $2}')
if [ "$NP" = 1 ]; then
    runApplication simpleFoam
else
    runApplication decomposePar
    runParallel simpleFoam
    runApplication reconstructPar -latestTime
fi
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--geometry", required=True,
                    help="registry spec <candidate>:<kind>, e.g. baseline:geometry_2d_faired")
    ap.add_argument("--condition", required=True, choices=list(CONDITIONS))
    ap.add_argument("--level", required=True, choices=list(LEVELS))
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--np", type=int, default=4)
    ap.add_argument("--rfar", type=float, default=500.0)
    ap.add_argument("--tag", default="grid")
    ap.add_argument("--endtime", type=int, default=15000)
    args = ap.parse_args()

    parts = args.geometry.split(":")
    cand, kind = parts[0], parts[1]
    reg = yaml.safe_load((REPO / "registry/candidates.yaml").read_text())
    entry = reg[cand][kind]
    if len(parts) == 3:
        entry = entry[int(parts[2])]   # list entries, e.g. geometry_2d_morphed:1
    dat = REPO / entry["file"]
    sha = hashlib.sha256(dat.read_bytes()).hexdigest()
    if sha != entry["sha256"]:
        sys.exit(f"registry gate: sha256 mismatch for {entry['file']}")

    cond = CONDITIONS[args.condition]
    # PER-CONDITION CHORD. A benchmark section is not our section; using the
    # project chord for it would silently change its Reynolds number.
    chord = cond.get("chord", CHORD)
    n, er, t1 = LEVELS[args.level]
    atag = f"a{args.alpha:g}".replace(".", "p").replace("-", "m")
    case = REPO / f"cases/{cand}/{args.condition}2d_{args.tag}_{args.level}_{atag}"
    if case.exists():
        shutil.rmtree(case)
    (case / "system").mkdir(parents=True)
    (case / "constant").mkdir()
    (case / "0.orig").mkdir()

    tpl = REPO / "template_case/airfoil2d"
    for f in ("controlDict", "fvSchemes", "fvSolution"):
        shutil.copy(tpl / f"system/{f}", case / f"system/{f}")
    shutil.copy(tpl / "constant/turbulenceProperties", case / "constant/turbulenceProperties")
    (case / "constant/transportProperties").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object transportProperties; }\n"
        "transportModel Newtonian;\n"
        f"nu [0 2 -1 0 0 0 0] {cond['nu']};   // {args.condition} (D006)\n")

    a = math.radians(args.alpha)
    ca, sa = math.cos(a), math.sin(a)
    subs = {"ux": f"{cond['U']*ca:.8f}", "uz": f"{cond['U']*sa:.8f}",
            "k": cond["k"], "omega": cond["omega"]}
    for name, txt in FIELDS.items():
        (case / f"0.orig/{name}").write_text(txt.format(**subs))

    cd = case / "system/controlDict"
    t = cd.read_text()
    t = re.sub(r"liftDir \(0 0 1\); dragDir \(1 0 0\);",
               f"liftDir ({-sa:.12f} 0 {ca:.12f}); dragDir ({ca:.12f} 0 {sa:.12f});", t, count=1)
    t = t.replace("endTime 4000", f"endTime {args.endtime}")
    t = t.replace("magUInf 34.0", f"magUInf {cond['U']}")
    # RESCALE THE FORCE REFERENCE WITH THE CHORD. The template hard-codes lRef/Aref 0.3404,
    # the EET section chord, and its own comment says "rescale with chord". Nothing did.
    # CR-2215 is a 0.6096 m chord, so every coefficient it produced was too large by
    # 0.6096/0.3404 = 1.7908: cl at alpha 0 read 0.489 against a measured 0.28, and the
    # lift slope read 0.2034 /deg, which is nearly double thin-aerofoil theory and
    # therefore impossible. The chord was ALREADY KNOWN here (line above sets `chord` from
    # the condition and passes it to the mesher); it simply never reached the force
    # reference. Ninth instance of the frame class: a correct number under an unstated frame.
    t = re.sub(r"lRef 0\.3404; Aref 0\.3404;",
               f"lRef {chord}; Aref {chord};", t, count=1)
    t = re.sub(r"CofR \(0\.0851 0 0\);",
               f"CofR ({0.25*chord:.6g} 0 0);", t, count=1)
    cd.write_text(t)

    fs = case / "system/fvSolution"
    t = fs.read_text()
    t = t.replace('residualControl { p 1e-5; U 1e-6; "(k|omega)" 1e-6; }',
                  'residualControl { p 1e-7; U 1e-8; "(k|omega)" 1e-8; }  // grid-study grade')
    t = t.replace('    "(U|k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }',
                  '''    "(U|k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
    Phi { solver PCG; preconditioner DIC; tolerance 1e-6; relTol 0.01; }   // potentialFoam init''')
    t = t.replace("cache { grad(U); }",
                  "potentialFlow { nNonOrthogonalCorrectors 10; PhiRefCell 0; PhiRefValue 0; }\ncache { grad(U); }")
    fs.write_text(t)

    (case / "system/decomposeParDict").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }\n"
        f"numberOfSubdomains {args.np};\nmethod scotch;\n")

    subprocess.run([sys.executable, str(REPO / "scripts/ogrid2d.py"),
                    "--dat", str(dat), "--chord", str(chord), "--n", str(n),
                    "--er", str(er), "--t1", str(t1),
                    "--rfar-chords", str(args.rfar), "--case", str(case)], check=True)

    (case / "Allrun").write_text(ALLRUN)
    (case / "Allrun").chmod(0o755)
    meta = {
        "candidate": cand, "geometry": entry["file"], "geometry_sha256": sha,
        "chord_m": chord, "alpha_deg": args.alpha, "condition": args.condition,
        "U_m_s": cond["U"], "nu_m2_s": cond["nu"], "Re": cond["Re"],
        "mach_nominal": cond["mach"], "level": args.level,
        "mesh": f"D015 O-grid, level {args.level} (n {n}, ER {er}, t1 {t1}), rfar ~{args.rfar:g}c",
        "lessons": "D012 Layer-1: 500c farfield, turbulent-branch init, y+~0.5 fine",
    }
    (case / "build_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"built {case.relative_to(REPO)}")


if __name__ == "__main__":
    main()
