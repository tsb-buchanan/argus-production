# airfoil2d template (OpenFOAM-org 7)

2D mirror of wing3d (D009/D011 conventions: org-7 syntax, simpleFoam + kOmegaSST,
y+~1 low-Re with the same absolute 28-layer stack, second-order turbulence schemes,
SIMPLEC). Plane: x streamwise, z vertical, y spanwise (single cell, empty patches),
so the trim script's inlet-vector/force-direction rotation is identical to 3D.

## Mesh chain (Allmesh)
snappyHexMesh has no 2D mode, so the mesh is built in a spanwise slab and flattened:
1. blockMesh: slab x [-17, 25.5] x z [-17, 17] m (50c up / 75c down / 50c v for
   c = 0.3404), ONE base cell (0.425 m cube) across y.
2. surfaceFeatures (org-7 utility) on constant/triSurface/airfoil.stl.
3. snappyHexMesh (parallel) castellate + snap + ABSOLUTE layers in the slab.
4. reconstructParMesh, then extrudeMesh: the front patch (y=0) is re-extruded to a
   single 1 m cell; side faces inherit inlet/outlet/farfield/airfoil from edge
   adjacency; front/back become the 2D planes.
5. changeDictionary sets front/back to empty; checkMesh.

Slab cost note: octree refinement is isotropic, so the slab holds 2^level cells
across y in refined zones; peak slab size ~5-8M cells at surface level (8 9),
inside the laptop ceiling but minutes, not seconds. The extruded 2D mesh is
~100-150k cells + 28x~800 layer cells.

## Per-case parameters (set by case builder / trim script)
1. Geometry: dat2stl.py --dat <registered .dat> --chord <local chord m> writes
   constant/triSurface/airfoil.stl. Registered sha256 only (registry gate).
2. Alpha: inlet vector (Ux = U cos a, Uz = U sin a) in 0.orig/U + rotated
   liftDir/dragDir in controlDict forceCoeffs. Geometry stays fixed.
3. Condition (D006): CR defaults set (U 34.0, k 1.73e-3, omega 22); WT values
   (U 40.8, k 2.50e-3, omega 27) in the file comments.
4. forceCoeffs is PER UNIT SPAN: extrusion thickness 1 m, Aref = chord x 1 m,
   lRef = chord, CofR = quarter chord. Rescale all three when the chord changes.

## GRID-STUDY VARIABLES (not frozen, D012 Layer 2 gates the freeze)
Surface level (8 9), box levels, nSurfaceLayers 28, firstLayerThickness 8e-6,
base cell 0.425 m. Deviations from the frozen recipe get a decision entry.

## Checklist before first run
STL regenerated from a REGISTERED .dat at the correct chord; locationInMesh
outside the airfoil; per-condition U/k/omega/magUInf consistent across 0.orig
and controlDict; run card written and validated (scripts/validate_run_card.py).
