#!/usr/bin/env python3
"""trefftz.py: induced drag from a Trefftz plane, the way a VLM defines it.

Run with pvbatch:
    pvbatch scripts/trefftz.py --case <caseDir> --alpha <deg> [--x 1.5 2.5 4.0]

WHY THIS EXISTS. Our surface integration gives C_Dp, which contains induced drag AND form
drag together. The VLM reports C_Di, induced ALONE. Comparing them is not like-for-like,
and the comparison currently produces a SIGN VIOLATION:

    VLM  C_Di (induced only)          98.79 ct
    RANS C_Dp (induced + form)        88.29 ct     <- should be the LARGER of the two

A Trefftz plane isolates induced drag using the same definition the VLM uses, so the two
become directly comparable and the violation can be explained rather than noted.

THE FORMULA, and its frame.
    D_i = (rho/2) * INTEGRAL over the plane of (v'^2 + w'^2) dS
    C_Di = D_i / (0.5 rho U^2 A_ref)
         = INTEGRAL (v'^2 + w'^2) dS / (U^2 A_ref)

    v', w' are CROSSFLOW PERTURBATIONS, so the freestream must be removed first:
        v' = U_y - 0
        w' = U_z - U*sin(alpha)
    Forgetting the w' offset adds a uniform U^2 sin^2(alpha) over the whole plane and
    swamps the answer -- at alpha 2.1 that is 0.0013*U^2 per unit area against a signal
    concentrated in the vortex core.

    A_ref is the HALF-wing 0.620462 m2, matching the case's forceCoeffs, because the
    plane covers the half model only.

SAMPLED AT SEVERAL DOWNSTREAM STATIONS ON PURPOSE. Too close and the wake has not rolled
up; too far and numerical diffusion has spread it and the integral decays. Reporting one
station would hide that; reporting the trend shows whether a plateau exists.
"""
import argparse
import os
import sys

try:
    from paraview.simple import *          # noqa: F401,F403
    from paraview import servermanager as sm
except ImportError:
    sys.exit("run with pvbatch, not python3")

import math

paraview.simple._DisableFirstRenderCameraReset()
AREF_HALF = 0.620462


def mpi_rank_size():
    """(rank, size) of this pvbatch process.

    THIS DISTINGUISHES TWO SITUATIONS THAT LOOK IDENTICAL IN A LOG. Under `srun -n 48`,
    an MPI-ENABLED pvbatch is ONE parallel job with size 48; a NON-MPI build is 48
    INDEPENDENT SERIAL ParaViews, each reading the whole 113M-cell case and each printing
    its own full result. The second is slower than plain serial and its interleaved
    output looks like a working parallel run. Read it, do not assume it.
    """
    try:
        from paraview import servermanager
        c = servermanager.vtkProcessModule.GetProcessModule().GetGlobalController()
        return c.GetLocalProcessId(), c.GetNumberOfProcesses()
    except Exception:
        return 0, 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--uinf", type=float, default=40.8)
    ap.add_argument("--x", type=float, nargs="+", default=None,
                    help="absolute x of the planes; default = TE + [0.5,1,2,4,8] chords")
    ap.add_argument("--reconstructed", action="store_true")
    ap.add_argument("--te", type=float, default=None,
                    help="trailing-edge x [m]. Given explicitly, the wing patch is never "
                    "read: one fewer full pass over a 113M-cell case, and one fewer way "
                    "to fail silently.")
    ap.add_argument("--chord", type=float, default=None,
                    help="reference chord [m] for expressing plane offsets")
    a = ap.parse_args()
    if (a.te is None) != (a.chord is None):
        sys.exit("--te and --chord are given together or not at all")

    rank, nproc = mpi_rank_size()
    ntasks = int(os.environ.get("SLURM_NTASKS", "0"))
    if rank == 0:
        print("pvbatch: MPI rank %d of %d (SLURM_NTASKS=%s)" % (rank, nproc, ntasks or "-"))
    # A NON-MPI BUILD UNDER srun IS A SILENT PERFORMANCE AND CORRECTNESS TRAP: every rank
    # would integrate the WHOLE domain and print a full, plausible, duplicated answer.
    # Refuse it rather than let 48 copies of the same number reach the log.
    if ntasks > 1 and nproc == 1:
        sys.exit("FATAL: launched under srun with SLURM_NTASKS=%d but this pvbatch reports "
                 "a single MPI process, i.e. it is NOT MPI-enabled. %d independent serial "
                 "ParaViews would each read the entire case and print duplicate results. "
                 "Run with `pvbatch` directly (serial) or use an MPI-enabled build."
                 % (ntasks, ntasks))
    if rank != 0:
        # only rank 0 reports; the integral itself is global
        sys.stdout = open(os.devnull, "w")

    stub = os.path.join(a.case, "case.foam")
    if not os.path.exists(stub):
        open(stub, "w").close()
    src = OpenFOAMReader(FileName=stub)
    src.CaseType = "Reconstructed Case" if a.reconstructed else "Decomposed Case"
    src.UpdatePipelineInformation()
    times = list(src.TimestepValues or [])
    if not times:
        sys.exit("FATAL: no time directories (CaseType=%s)" % src.CaseType)
    t = times[-1]
    src.MeshRegions = ["internalMesh"]
    src.CellArrays = ["U"]
    src.UpdatePipeline(t)
    b = src.GetDataInformation().GetBounds()
    print("case %s   latest time %s" % (a.case, t))
    print("  domain x %.2f..%.2f  y %.2f..%.2f  z %.2f..%.2f" % b)

    if a.te is not None:
        te, chord = a.te, a.chord
        print("  TE %.4f, ref chord %.4f (given, wing patch not read)" % (te, chord))
    else:
        # wing trailing edge, for placing the planes in chords
        wing = OpenFOAMReader(FileName=stub)
        wing.CaseType = src.CaseType
        wing.UpdatePipelineInformation()
        avail = list(wing.GetProperty("MeshRegions"))
        if "patch/wing" not in avail:
            sys.exit("FATAL: no 'patch/wing' region. Available: %s\n"
                     "  (a DECOMPOSED case carries the wing only in processor*/; the "
                     "top-level constant/polyMesh/boundary is the blockMesh background "
                     "and lists no wing at all)" % avail[:12])
        wing.MeshRegions = ["patch/wing"]
        wing.UpdatePipeline(t)
        wb = wing.GetDataInformation().GetBounds()
        # AN EMPTY VTK DATASET REPORTS BOUNDS OF +/-1e308, NOT AN ERROR. Without this the
        # run continues on garbage: every plane lands outside the domain, every one is
        # "skipped", and the script prints its footer and EXITS 0 having computed nothing.
        # That is exactly what job 10579302 did for 2h55m. The same guard already existed
        # in pv_postprocess.py after the identical failure and was not carried here.
        if not all(abs(v) < 1e30 for v in wb) or wb[1] <= wb[0]:
            sys.exit("FATAL: wing patch read returned an EMPTY dataset (bounds %r). "
                     "Nothing downstream of this is meaningful; pass --te/--chord to "
                     "place the planes without reading the patch." % (wb,))
        te, chord = wb[1], wb[1] - wb[0]
        print("  wing x %.4f..%.4f (TE %.4f, ref chord %.4f)" % (wb[0], wb[1], te, chord))

    xs = a.x if a.x else [te + f * chord for f in (0.5, 1.0, 2.0, 4.0, 8.0)]
    wfs = a.uinf * math.sin(math.radians(a.alpha))          # freestream w
    print("  freestream w removed: %.5f m/s" % wfs)
    print()
    print("  %10s %8s %14s %12s" % ("x", "chords", "CDi [counts]", "area m2"))
    got = []
    for x in xs:
        if not (b[0] < x < b[1]):
            print("  %10.4f   outside the domain, skipped" % x)
            continue
        sl = Slice(Input=src)
        sl.SliceType = "Plane"
        sl.SliceType.Origin = [x, 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5])]
        sl.SliceType.Normal = [1, 0, 0]
        # crossflow kinetic energy per unit area, freestream removed
        c = Calculator(Input=sl)
        c.AttributeType = "Point Data"
        c.ResultArrayName = "q2"
        c.Function = "U_Y*U_Y + (U_Z-%.9f)*(U_Z-%.9f)" % (wfs, wfs)
        integ = IntegrateVariables(Input=c)
        # PORTABILITY GUARD, EQUIVALENCE PROVEN NOT ASSUMED. HPC12's only ParaView is
        # /usr/bin/pvbatch 4.4.0, which has no DivideCellDataByVolume; the property was
        # added in 5.x and setting it there REQUESTS SUMMATION. 4.4 sums unconditionally:
        # measured with a known answer, a 2x3 Plane carrying 77 cells integrates to
        # Area = 6.0 exactly, not 6/77. So omitting it where absent leaves this integral
        # bit-for-bit identical. Without the guard the script aborts with AttributeError
        # before any plane is evaluated, which is how this was found.
        try:
            integ.DivideCellDataByVolume = 0
        except AttributeError:
            pass
        integ.UpdatePipeline(t)
        d = sm.Fetch(integ)
        pd = d.GetPointData().GetArray("q2")
        cd = d.GetCellData().GetArray("Area")
        if pd is None or cd is None:
            print("  %10.4f   integration returned nothing" % x)
            Delete(integ); Delete(c); Delete(sl); continue
        q2 = pd.GetValue(0)                 # already area-integrated by IntegrateVariables
        area = cd.GetValue(0)
        cdi = q2 / (a.uinf ** 2 * AREF_HALF)
        print("  %10.4f %8.2f %14.3f %12.4f"
              % (x, (x - te) / chord, 1e4 * cdi, area))
        got.append((x, 1e4 * cdi))
        Delete(integ); Delete(c); Delete(sl)
    print()
    print("  frame: C_Di on A_ref %.6f m2 (HALF wing, matches forceCoeffs), U %.1f m/s."
          % (AREF_HALF, a.uinf))
    print("  Compare against VSPAERO C_Di, which is induced drag ALONE on the same basis.")
    # A RUN THAT INTEGRATED NOTHING MUST NOT EXIT 0. The footer above prints whether or
    # not any plane produced a number, so the exit status is the only thing that can
    # distinguish "no induced drag found" from "no plane was ever evaluated".
    if not got:
        sys.exit("FATAL: every plane was skipped; no Trefftz integral was computed.")
    print("\n  %d of %d planes integrated; CDi spread %.3f to %.3f counts"
          % (len(got), len(xs), min(v for _, v in got), max(v for _, v in got)))
    return 0


if __name__ == "__main__":
    rc = main()
    # HARD EXIT, DELIBERATELY. Under `srun -n 48`, pvbatch 5.9.1 printed every result for
    # the first case and then NEVER RETURNED: the interpreter hung in MPI finalize during
    # shutdown, so the surrounding loop never reached the second case and the job burned
    # its 3h wall with the answer already sitting in the log. os._exit skips interpreter
    # cleanup and the MPI teardown with it. stdout is flushed FIRST, because os._exit does
    # not flush and the whole point is that the results survive.
    sys.stdout.flush()
    sys.stderr.flush()
    # THIS DID NOT CURE THE HANG AND THE HONEST RECORD MATTERS MORE THAN THE TIDY ONE.
    # os._exit runs on EVERY rank here (the rank-0 guard above covers printing only), and
    # job 10580838 STILL sat RUNNING for 35 minutes after printing its complete result at
    # 3.5 minutes, until cancelled. So the hang is below Python -- in pvbatch/MPI teardown
    # or in srun's step reaping -- and is UNRESOLVED.
    #
    # THE MITIGATION THAT ACTUALLY WORKS IS OPERATIONAL, NOT CODE: run ONE CASE PER JOB and
    # give it a short wall. The numbers are flushed and complete long before the hang, so a
    # job killed at its limit has already delivered; a two-case job loses the second case
    # entirely, which is what happened to 10579302 (2h55m, nothing).
    os._exit(rc or 0)
