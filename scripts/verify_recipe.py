#!/usr/bin/env python3
"""Verify a recipe folder: manifest intact, no suffixed variants, dicts agree with
the pinned scalars.

THE FAILURE THIS CATCHES. hpc/v5_recipe.json pins the recipe as scalars a human can
read. hpc/recipes/<name>/ holds the dicts OpenFOAM actually reads. Nothing keeps
those two in step, and a pinned value that has quietly stopped describing the dict
beside it is worse than no pin: it is a confident, checkable-looking claim that is
wrong. That is not hypothetical here either. v5_recipe.json first pinned
surfaceLevel 8, read from meshLevels.full, while the mesh that ran was
surfaceLevel 11 from the unsuffixed meshLevels.

So this asserts, in both directions, that the number in the JSON is the number in
the dict. Run it before any mesh job and in the pre-commit sweep.

CHECK 4 IS THE ONE THAT NEEDS ITS OWN NEGATIVE CONTROL (GEO-089 fixture rule): a
comparison that silently finds no keys to compare would report PASS over an empty
set. It therefore asserts a minimum key count and --self-test mutates a real dict
to prove the comparison can fail.
"""
import hashlib, json, pathlib, re, shutil, sys, tempfile

BANNED = (".wallResolved", ".wallModelled", ".full", ".sandbox")

# pinned-JSON path -> (dict file, key as written in the dict)
PINNED = {
    ("layers", "nSurfaceLayers"):      ("system/layerSettings", "nSurfaceLayers"),
    ("layers", "expansionRatio"):      ("system/layerSettings", "expansionRatio"),
    ("layers", "firstLayerThickness"): ("system/layerSettings", "firstLayerThickness"),
    ("layers", "minThickness"):        ("system/layerSettings", "minThickness"),
    ("layers", "featureAngle"):        ("system/layerSettings", "featureAngle"),
    ("levels", "surfaceLevelMin"):     ("system/meshLevels", "surfaceLevelMin"),
    ("levels", "surfaceLevelMax"):     ("system/meshLevels", "surfaceLevelMax"),
    ("levels", "levelOuter"):          ("system/meshLevels", "levelOuter"),
    ("levels", "levelMid"):            ("system/meshLevels", "levelMid"),
    ("levels", "levelNear"):           ("system/meshLevels", "levelNear"),
    ("levels", "levelWing"):           ("system/meshLevels", "levelWing"),
    ("levels", "levelWake"):           ("system/meshLevels", "levelWake"),
    ("levels", "levelCone"):           ("system/meshLevels", "levelCone"),
    ("levels", "levelMorph"):          ("system/meshLevels", "levelMorph"),
}
MIN_COMPARED = 12


def scalar(path, key):
    """Read `key <value>;` from an OpenFOAM dict. Returns None if absent.
    Comments are stripped first: several of these dicts carry `// 3.906 mm` after
    the value and a naive match picks up the comment."""
    txt = re.sub(r"//.*", "", pathlib.Path(path).read_text(errors="ignore"))
    m = re.search(rf"\b{re.escape(key)}\s+([-\d.eE+]+)\s*;", txt)
    if not m:
        m = re.search(rf"\b{re.escape(key)}\s+(\S+?)\s*;", txt)
    return m.group(1) if m else None


def check(folder, pinned_json):
    folder = pathlib.Path(folder)
    pin = json.loads(pathlib.Path(pinned_json).read_text())
    fails, notes = [], []

    # 1. manifest
    man = folder / "MANIFEST.sha256"
    if not man.is_file():
        fails.append("MANIFEST.sha256 missing")
    else:
        listed = {}
        for line in man.read_text().splitlines():
            if line.strip():
                h, _, rel = line.partition("  ")
                listed[rel] = h
        for rel, h in listed.items():
            p = folder / rel
            if not p.is_file():
                fails.append(f"manifest lists {rel} but it is absent")
            elif hashlib.sha256(p.read_bytes()).hexdigest() != h:
                fails.append(f"{rel} does not match its manifest hash")
        on_disk = {str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file()}
        extra = on_disk - set(listed) - {"MANIFEST.sha256", "RECIPE_FILES.json"}
        if extra:
            fails.append(f"files present but unlisted: {sorted(extra)}")
        notes.append(f"manifest: {len(listed)} files verified")

    # 2. no suffixed variants. The structural guarantee.
    stray = [str(p.relative_to(folder)) for p in folder.rglob("*")
             if p.is_file() and any(str(p).endswith(b) for b in BANNED)]
    if stray:
        fails.append(f"suffixed variants present: {stray}")
    else:
        notes.append("no suffixed variants: the folder is unambiguous")

    # 3/4. pinned scalars vs dicts
    compared = 0
    for (sec, jkey), (rel, dkey) in PINNED.items():
        p = folder / rel
        if not p.is_file():
            fails.append(f"{rel} absent, cannot check {sec}.{jkey}")
            continue
        if sec not in pin or jkey not in pin[sec]:
            fails.append(f"{pinned_json} has no {sec}.{jkey} to check against")
            continue
        got, want = scalar(p, dkey), pin[sec][jkey]
        if got is None:
            fails.append(f"{rel} does not define {dkey}")
            continue
        try:
            same = abs(float(got) - float(want)) <= 1e-12 * max(1.0, abs(float(want)))
        except ValueError:
            same = str(got) == str(want)
        compared += 1
        if not same:
            fails.append(f"{sec}.{jkey}: pinned {want}, {rel} says {got}")
    if compared < MIN_COMPARED:
        fails.append(f"only {compared} scalars compared, expected at least "
                     f"{MIN_COMPARED}. A comparison over an empty set is not a pass.")
    notes.append(f"pinned scalars compared: {compared}")

    # 5. THE ALPHA TRIPLE. alpha is written into THREE places that must agree:
    #    0.orig/U's Uinf vector, and forceCoeffs1's liftDir and dragDir. They are
    #    regenerated per case, so this is the check most likely to catch a real
    #    campaign error rather than a stale template.
    #    D019 IS THE PRECEDENT: forceCoeffs axes lagged the U rotation by up to
    #    0.0083 deg and biased cd by cl * lag. A wrong-by-a-degree dragDir is a
    #    LARGE drag error that no residual, force history or mesh check can see,
    #    because every one of them is perfectly happy with it.
    import math
    u = folder / "0.orig/U"
    cd_ = folder / "system/controlDict"
    if u.is_file() and cd_.is_file():
        ut = u.read_text(errors="ignore")
        ct = re.sub(r"//.*", "", cd_.read_text(errors="ignore"))
        mu = re.search(r"Uinf\s*\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", ut)
        md = re.search(r"dragDir\s*\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", ct)
        ml = re.search(r"liftDir\s*\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", ct)
        mm = re.search(r"magUInf\s+([-\d.eE+]+)\s*;", ct)
        if not (mu and md and ml and mm):
            fails.append("could not read the alpha triple (Uinf, dragDir, liftDir, magUInf); "
                         "that is the ABSENCE of a check, not a pass")
        else:
            ux, _, uz = (float(g) for g in mu.groups())
            dx, _, dz = (float(g) for g in md.groups())
            lx, _, lz = (float(g) for g in ml.groups())
            mag = math.hypot(ux, uz)
            a_u = math.degrees(math.atan2(uz, ux))
            a_d = math.degrees(math.atan2(dz, dx))
            a_l = math.degrees(math.atan2(-lx, lz))
            # (a) magUInf must equal |Uinf|, else every coefficient is scaled wrong
            if abs(mag - float(mm.group(1))) > 1e-3:
                fails.append(f"magUInf {mm.group(1)} != |Uinf| {mag:.6f}. Every coefficient "
                             f"would be scaled by the square of that ratio.")
            # (b) the axes must follow U to well inside the D019 lag
            for lab, a in (("dragDir", a_d), ("liftDir", a_l)):
                if abs(a - a_u) > 1e-3:
                    fails.append(f"{lab} implies alpha {a:.6f} deg but Uinf implies "
                                 f"{a_u:.6f} deg, a lag of {a-a_u:+.6f} deg (D019)")
            # (c) lift must be perpendicular to drag, and both unit
            if abs(lx * dx + lz * dz) > 1e-6:
                fails.append(f"liftDir.dragDir = {lx*dx+lz*dz:.2e}, not orthogonal")
            for lab, n in (("dragDir", math.hypot(dx, dz)), ("liftDir", math.hypot(lx, lz))):
                if abs(n - 1.0) > 1e-5:
                    fails.append(f"|{lab}| = {n:.9f}, not unit; forceCoeffs does not "
                                 f"normalise it, so the coefficient is scaled by 1/{n:.9f}")
            # ONLY claim agreement if nothing above disagreed. The first version
            # appended this note unconditionally and printed "alpha triple agrees"
            # directly beside a FAIL line saying the axes lag by 0.8 deg. An output
            # that simultaneously asserts and disclaims the same result is the exact
            # shape record-keeping rule 5 forbids in the decision log, and it has no
            # more business in a gate.
            if not any("dragDir" in x or "liftDir" in x or "magUInf" in x for x in fails):
                notes.append(f"alpha triple agrees: |U| {mag:.4f} m/s, alpha {a_u:.4f} deg, "
                             f"axis lag {a_d-a_u:+.6f} deg, orthogonal, unit")
            else:
                notes.append(f"alpha triple CHECKED AND REJECTED: |U| {mag:.4f} m/s, "
                             f"alpha from Uinf {a_u:.4f} deg (see FAIL lines)")
    else:
        fails.append("0.orig/U or system/controlDict absent; the alpha triple was not checked")

    return fails, notes


def self_test():
    """NEGATIVE CONTROL. A verifier that has only ever passed is doing no work (D058)."""
    src = pathlib.Path("hpc/recipes/L11_wallResolved_WT")
    pin = "hpc/v5_recipe.json"
    if not src.is_dir():
        print("self-test SKIPPED: recipe folder absent"); return 1
    f, _ = check(src, pin)
    ok0 = not f
    print(f"  1. unmodified recipe verifies: {'PASS' if ok0 else 'FAIL ' + str(f)}")
    results = [ok0]
    with tempfile.TemporaryDirectory() as td:
        # (a) a changed VALUE must be caught by the pinned-scalar comparison
        d = pathlib.Path(td) / "mutval"; shutil.copytree(src, d)
        t = (d / "system/layerSettings").read_text().replace("1.450", "1.180")
        (d / "system/layerSettings").write_text(t)
        f2, _ = check(d, pin)
        hit = any("expansionRatio" in x for x in f2)
        print(f"  2. expansionRatio 1.45->1.18 is caught: {'PASS' if hit else 'FAIL'}")
        results.append(hit)
        # (b) ... and independently by the manifest, which is the point of having both
        hit_m = any("manifest hash" in x for x in f2)
        print(f"  3. the same edit trips the manifest independently: {'PASS' if hit_m else 'FAIL'}")
        results.append(hit_m)
        # (c) a reintroduced suffixed variant must be caught
        d3 = pathlib.Path(td) / "stray"; shutil.copytree(src, d3)
        shutil.copy2(d3 / "system/layerSettings", d3 / "system/layerSettings.wallResolved")
        f3, _ = check(d3, pin)
        hit3 = any("suffixed variants" in x for x in f3)
        print(f"  4. a reintroduced .wallResolved is caught: {'PASS' if hit3 else 'FAIL'}")
        results.append(hit3)
        # (d) a deleted dict must not silently shrink the compared set into a pass
        d4 = pathlib.Path(td) / "gone"; shutil.copytree(src, d4)
        (d4 / "system/meshLevels").unlink()
        f4, _ = check(d4, pin)
        hit4 = any("only" in x and "compared" in x for x in f4) or any("absent" in x for x in f4)
        print(f"  5. a deleted meshLevels fails rather than narrowing: {'PASS' if hit4 else 'FAIL'}")
        results.append(hit4)
        # (e) THE D019 CASE: a dragDir left at the previous alpha while U moved on.
        #     This is the single most likely real error in an 18-case campaign, and
        #     nothing else in the pipeline can see it.
        d5 = pathlib.Path(td) / "lag"; shutil.copytree(src, d5)
        t5 = (d5 / "system/controlDict").read_text().replace(
            "dragDir         ( 0.999327 0 0.0366925 )",
            "dragDir         ( 0.999742 0 0.0227 )")
        (d5 / "system/controlDict").write_text(t5)
        f5, _ = check(d5, pin)
        hit5 = any("dragDir implies alpha" in x for x in f5)
        print(f"  6. a dragDir left at a stale alpha is caught: {'PASS' if hit5 else 'FAIL'}")
        results.append(hit5)
        # (f) magUInf disagreeing with |Uinf| scales every coefficient
        d6 = pathlib.Path(td) / "mag"; shutil.copytree(src, d6)
        t6 = (d6 / "system/controlDict").read_text().replace("magUInf         40.8;",
                                                             "magUInf         34.0;")
        (d6 / "system/controlDict").write_text(t6)
        f6, _ = check(d6, pin)
        hit6 = any("magUInf" in x for x in f6)
        print(f"  7. magUInf 40.8 -> 34.0 (the WT/CR mixup) is caught: {'PASS' if hit6 else 'FAIL'}")
        results.append(hit6)
    print(f"self-test {sum(results)}/{len(results)}")
    return 0 if all(results) else 1


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(self_test())
    folder = sys.argv[1] if len(sys.argv) > 1 else "hpc/recipes/L11_wallResolved_WT"
    pinned = sys.argv[2] if len(sys.argv) > 2 else "hpc/v5_recipe.json"
    fails, notes = check(folder, pinned)
    for n in notes:
        print(f"  {n}")
    # REPORT FIRST, EXIT AFTER (GEO-087 item 5).
    for f in fails:
        print(f"  FAIL: {f}")
    print("RECIPE OK" if not fails else f"RECIPE REFUSED: {len(fails)} problems")
    sys.exit(0 if not fails else 1)
