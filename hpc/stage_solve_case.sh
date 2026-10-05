#!/bin/bash
# stage_solve_case.sh -- turn a finished MESH case into a SOLVE case for hpc/solve_hpc12.pbs.
#
# usage: hpc/stage_solve_case.sh <mesh-case-dir> <solve-case-dir> [--ranks N]
#
#   <mesh-case-dir>   a case meshed by hpc/mesh_hpc12.pbs from a recipe folder: the
#                     reconstructed constant/polyMesh, 0.orig/, system/, the other
#                     constant/ dictionaries, constant/triSurface/, MANIFEST.sha256,
#                     RECIPE_FILES.json, and CONDITION.json for a cruise recipe.
#   <solve-case-dir>  must NOT exist yet. It is created; nothing is ever overwritten.
#   --ranks N         the solve rank count. Default 128. ARGUS_RANKS=N in the
#                     environment also sets it; the flag wins over the environment.
#
# THE RANK RULE. The six wake-refined Condition CR cases (the r2 cases B, C, F, H, M, W)
# were meshed at 192 ranks and SOLVED AT 128 RANKS, method hierarchical. 128 is the
# default so that a case staged here is split the way the delivered numbers were:
#     128             hierarchical (8 2 8), delta 0.001, order xyz: the r2 solves
#     192             hierarchical (8 3 8), delta 0.001, order xyz: the published solves'
#                     case definitions
#     any other N     scotch, with a printed note that this departs from both
# This is the SOLVE decomposition only. The mesh case keeps the recipe's own
# decomposeParDict, the mesh-stage one, because parallel snappyHexMesh builds a
# different mesh on a different decomposition; re-decomposing a finished mesh does not.
# The 128 block is the one read back from the six r2 solve cases on HPC12: every one of
# them carries 128 subdomains, method hierarchical, n (8 2 8), delta 0.001, order xyz.
# WHY 128: fpt-medium takes at most 4 nodes per job and 4 x typej is 128 cores, with up
# to three jobs per user; fpt-large takes ONE job per user. A different rank count is a
# different partition, so the default is the one the delivered results came from.
#
# WHAT IT DOES. Nothing is created until every input gate has passed.
#   1. Gates the mesh case as an INPUT: checkMesh passed every flooded mesh this project
#      built, because a gate on the output cannot see a wrong input. It also gates every
#      path the printed qsub line carries in its -v list (the solve case, solve_hpc12.pbs
#      and the post-processing checker), resolved to the absolute paths it prints.
#   2. Copies into <solve-case-dir>.incomplete_<UTC stamp>, a sibling directory, so that a
#      copy which dies half way never sits under the real name.
#   3. Writes system/decomposeParDict for the rank rule above.
#   4. Verifies FROM THE FILES: every mesh file by size, and by sha256 up to 64 MiB; every
#      other copied file by sha256; the staged file set equal to the source file set; no
#      symlinks; the decomposition read back with solve_hpc12.pbs's own sed expression;
#      and the conditions solve_hpc12.pbs checks on a cold start.
#   5. Writes SOLVE_SETUP.json from what it read back, renames the directory to
#      <solve-case-dir>, and prints the next steps: the post-processing block (and the
#      exact fix when the case lacks it), the angle, the iterations, a cd into the solve
#      case, and the qsub lines, every -v list carrying ARGUS_GATE.
#
# COPIES, NEVER SYMLINKS. A symlinked constant/polyMesh lets a purge inside one case
# follow the link and destroy a mesh that other cases share. Symlinks inside the mesh
# case are followed and their targets copied.
#
# NEVER DELETES. A failure after the gates leaves the .incomplete_ directory in place and
# names it; remove it yourself once you have read why it failed.
#
# NOT COPIED, deliberately: 0/ (solve_hpc12.pbs regenerates 0 from 0.orig on a cold
# start, and a leftover 0/ would run whatever it holds), processor* (the mesh job's own
# decomposition; the first solve leg decomposes afresh), log/, and anything else at the
# top of the mesh case. Every top-level entry is printed as copied or excluded.
#
# MANIFEST.sha256 is copied as the record of the recipe the mesh was built from. Its
# system/decomposeParDict entry describes the MESH decomposition and is expected not to
# match the file written here; SOLVE_SETUP.json records whether it does.
#
# Exit status: 0 staged; 2 refused (usage or an input gate, nothing created);
# 1 failed after staging began (the partial directory is named and left in place).
#
# Needs bash, coreutils, sed and grep only. It runs on the HPC12 login node and needs
# neither OpenFOAM nor python.

set -uo pipefail

HASH_MAX=67108864                  # mesh files up to 64 MiB are compared by sha256
REQMESH="boundary faces neighbour owner points"
USAGE="usage: $0 <mesh-case-dir> <solve-case-dir> [--ranks N]"

say()    { echo "  $*"; }
note()   { echo "  NOTE: $*"; }
refuse() { echo "REFUSED: $*" >&2; echo "Nothing was created." >&2; exit 2; }
STAGE=""
die() {
    echo "FAILED: $*" >&2
    if [ -n "$STAGE" ] && [ -e "$STAGE" ]; then
        echo "The partial case is left at $STAGE; nothing was deleted." >&2
        echo "Remove it yourself once you have read why it failed." >&2
    fi
    exit 1
}
sha() { sha256sum -- "$1" 2>/dev/null | cut -d' ' -f1; }
jstr() {
    local s
    s=$(printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g; s/[[:cntrl:]]/?/g')
    printf '"%s"' "$s"
}

# ---------------------------------------------------------------- arguments
NPOS=0; MESH=""; DST=""; FLAG_SET=0; FLAG_VAL=""
add_pos() {
    NPOS=$((NPOS + 1))
    if [ "$NPOS" -eq 1 ]; then MESH="$1"; elif [ "$NPOS" -eq 2 ]; then DST="$1"; fi
}
while [ $# -gt 0 ]; do
    case "$1" in
        --ranks)
            [ $# -ge 2 ] || refuse "--ranks needs a value. $USAGE"
            [ "$FLAG_SET" -eq 0 ] || refuse "--ranks given more than once"
            FLAG_SET=1; FLAG_VAL="$2"; shift 2 ;;
        --ranks=*)
            [ "$FLAG_SET" -eq 0 ] || refuse "--ranks given more than once"
            FLAG_SET=1; FLAG_VAL="${1#--ranks=}"; shift ;;
        -h|--help)
            echo "$USAGE"; exit 0 ;;
        --)
            shift; while [ $# -gt 0 ]; do add_pos "$1"; shift; done ;;
        -*)
            refuse "unknown option '$1'. $USAGE" ;;
        *)
            add_pos "$1"; shift ;;
    esac
done
[ "$NPOS" -eq 2 ] || refuse "expected a mesh case and a solve case, got $NPOS path(s). $USAGE"

# ---------------------------------------------------------------- rank count
# THE FLAG WINS, THEN THE ENVIRONMENT, THEN 128. The source is printed with the value,
# because a rank count nobody can trace is how a case ends up on a different partition.
if [ "$FLAG_SET" -eq 1 ]; then
    RAW="$FLAG_VAL"; RSRC="--ranks"; RWHY="from --ranks"
    [ -n "${ARGUS_RANKS:-}" ] && RWHY="from --ranks (ARGUS_RANKS='${ARGUS_RANKS}' ignored: the flag wins)"
elif [ -n "${ARGUS_RANKS:-}" ]; then
    RAW="$ARGUS_RANKS"; RSRC="ARGUS_RANKS"; RWHY="from the environment variable ARGUS_RANKS"
else
    RAW="128"; RSRC="default"; RWHY="the default, as the r2 solves; neither --ranks nor ARGUS_RANKS was given"
fi
case "$RAW" in
    ''|*[!0-9]*) refuse "rank count '$RAW' ($RSRC) is not a positive integer" ;;
esac
[ "${#RAW}" -le 6 ] || refuse "rank count '$RAW' ($RSRC) is not a plausible rank count"
NR=$((10#$RAW))
[ "$NR" -ge 1 ] || refuse "rank count '$RAW' ($RSRC) is below 1"

case "$NR" in
    128) METHOD=hierarchical; NX=8; NY=2; NZ=8
         RULE="128 ranks: hierarchical (8 2 8), the r2 solves' configuration" ;;
    192) METHOD=hierarchical; NX=8; NY=3; NZ=8
         RULE="192 ranks: hierarchical (8 3 8), the published solves' case definitions" ;;
    *)   METHOD=scotch; NX=""; NY=""; NZ=""
         RULE="$NR ranks: scotch, which departs from the delivered configuration" ;;
esac

echo "=== stage_solve_case.sh $(date -u '+%Y-%m-%dT%H:%M:%SZ') ==="
say "ranks: $NR ($RWHY)"
say "rule:  $RULE"
if [ "$METHOD" = scotch ]; then
    note "$NR RANKS DEPARTS FROM THE DELIVERED CONFIGURATION. The r2 solves ran at 128"
    note "(hierarchical 8 2 8) and the published solve definitions carry 192 (hierarchical"
    note "8 3 8). A result on a different decomposition is a result on a different"
    note "partition. scotch splits into any count but needs an OpenFOAM build that carries"
    note "the real scotch library."
fi

# ---------------------------------------------------------------- paths
while [ "${DST%/}" != "$DST" ] && [ "$DST" != "/" ]; do DST="${DST%/}"; done
case "$MESH$DST" in
    *[[:cntrl:]]*) refuse "a path contains a control character" ;;
esac

# INPUT GATE 1: the mesh case exists.
[ -d "$MESH" ] || refuse "mesh case '$MESH' does not exist or is not a directory"

# INPUT GATE 2: the solve case does NOT exist. Never overwrite, never delete.
if [ -e "$DST" ] || [ -L "$DST" ]; then
    refuse "solve case '$DST' already exists. This script never overwrites and never
         deletes: choose a new name, or move the old case aside yourself."
fi
MESH_ABS=$(cd "$MESH" && pwd -P) || refuse "cannot enter mesh case '$MESH'"
DST_ABS=$(realpath -m -- "$DST") || refuse "cannot resolve solve case path '$DST'"
case "$DST_ABS/" in
    "$MESH_ABS"/*) refuse "solve case '$DST' is the mesh case or lies inside it" ;;
esac

# INPUT GATE 2a: every path the printed qsub line puts in its -v list. That list is
# comma-separated and the printed line does not quote it, so a comma or a space in any
# of these paths breaks the submission. THE RESOLVED PATHS ARE CHECKED, because they are
# what gets printed: a relative solve-case path is clean as typed and can still resolve
# under a folder whose name holds a space or a comma.
SCRIPT_DIR=$(cd "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P) \
    || refuse "cannot resolve the folder this script is in"
REPO_DIR=$(dirname -- "$SCRIPT_DIR")
SOLVE_PBS="$SCRIPT_DIR/solve_hpc12.pbs"
GATE_PY="$REPO_DIR/scripts/assert_case_postpro.py"
vlist_ok() {
    case "$1" in *[[:space:],]*|*[[:cntrl:]]*) return 1 ;; esac
    return 0
}
vlist_ok "$DST_ABS" || refuse "the solve case's full path '$DST_ABS' contains a
         space or a comma. The qsub line carries it as ARGUS_CASE=... in its
         comma-separated -v list, which cannot hold either. Choose a solve case whose
         full path has neither."
for vp in "ARGUS_SCRIPT|$SOLVE_PBS" "ARGUS_GATE|$GATE_PY"; do
    vn=${vp%%|*}; vf=${vp#*|}
    [ -f "$vf" ] || continue
    vlist_ok "$vf" || refuse "'$vf' contains a space or a comma. The qsub line carries
         it as $vn=... in its comma-separated -v list, which cannot hold either.
         Clone the repository to a folder whose full path has neither, and run this
         script from there."
done

# ---------------------------------------------------------------- file walk
# Lists regular files under a directory, symlinks followed, one relative path per line.
# A dangling link, a special file or a control character in a name is collected as BAD,
# so the gates refuse it before anything is copied instead of cp failing half way.
WALK_OUT=""; WALK_BAD=""; WALK_LINKS=0
walk() {
    local d="$1" rel="$2" depth="$3" e b r
    if [ "$depth" -gt 24 ]; then
        WALK_BAD="$WALK_BAD
    nested deeper than 24 levels (a symlink loop?): $rel"
        return
    fi
    for e in "$d"/* "$d"/.[!.]* "$d"/..?*; do
        [ -e "$e" ] || [ -L "$e" ] || continue
        b=${e##*/}; r="$rel$b"
        case "$r" in
            *[[:cntrl:]]*) WALK_BAD="$WALK_BAD
    control character in a file name: $(printf '%q' "$r")"; continue ;;
        esac
        if [ -L "$e" ]; then
            WALK_LINKS=$((WALK_LINKS + 1))
            if [ ! -e "$e" ]; then
                WALK_BAD="$WALK_BAD
    dangling symlink: $r"; continue
            fi
        fi
        if [ -d "$e" ]; then
            if [ ! -r "$e" ] || [ ! -x "$e" ]; then
                WALK_BAD="$WALK_BAD
    unreadable directory: $r"; continue
            fi
            walk "$e" "$r/" $((depth + 1))
        elif [ -f "$e" ]; then
            if [ ! -r "$e" ]; then
                WALK_BAD="$WALK_BAD
    unreadable file: $r"; continue
            fi
            WALK_OUT="$WALK_OUT$r
"
        else
            WALK_BAD="$WALK_BAD
    neither a regular file nor a directory: $r"
        fi
    done
}
nlines() { [ -z "$1" ] && { echo 0; return; }; printf '%s' "$1" | grep -c ''; }

# ---------------------------------------------------------------- input gates on the mesh
PM="$MESH/constant/polyMesh"
NPROC=0
for p in "$MESH"/processor*; do [ -d "$p" ] && NPROC=$((NPROC + 1)); done

MISS=""; EMPTY=""; GZ=""
for f in $REQMESH; do
    if [ ! -f "$PM/$f" ]; then
        MISS="$MISS $f"
        [ -f "$PM/$f.gz" ] && GZ="$GZ $f.gz"
    elif [ ! -s "$PM/$f" ]; then
        EMPTY="$EMPTY $f"
    fi
done

# INPUT GATE 3: processor* WITHOUT a reconstructed mesh. After `snappyHexMesh -parallel
# -overwrite` the fine mesh exists only inside processor*, so a case in that state has
# nothing to copy. mesh_hpc12.pbs collects it with `reconstructPar -constant`.
if [ "$NPROC" -gt 0 ] && { [ ! -d "$PM" ] || [ -n "$MISS$EMPTY" ]; }; then
    refuse "mesh case '$MESH' has $NPROC processor* directories but no reconstructed
         constant/polyMesh (missing:${MISS:- none}; zero bytes:${EMPTY:- none}). The mesh
         lives only in processor*. Collect it first (mesh_hpc12.pbs runs
         'reconstructPar -constant'), then stage."
fi
# INPUT GATE 4: the reconstructed mesh is present, complete and non-empty.
[ -d "$PM" ] || refuse "mesh case '$MESH' has no constant/polyMesh"
[ -z "$GZ" ] || refuse "constant/polyMesh holds compressed files ($GZ ). The recipes write
         uncompressed meshes (writeCompression off) and this script does not handle
         compressed ones."
[ -z "$MISS" ] || refuse "constant/polyMesh in '$MESH' is missing:$MISS"
[ -z "$EMPTY" ] || refuse "constant/polyMesh in '$MESH' has ZERO-BYTE:$EMPTY"
# INPUT GATE 5: the mesh carries a wall patch. The background blockMesh has none, so a
# polyMesh without one is the box that blockMesh writes first, not the snappy mesh.
# The same test mesh_hpc12.pbs applies after reconstruction.
grep -qE '^[[:space:]]+wing$|type[[:space:]]+wall' "$PM/boundary" \
    || refuse "constant/polyMesh/boundary has no wall patch. That is the background
         blockMesh, not the snappy mesh; the mesh job did not finish."

# INPUT GATE 6: what solve_hpc12.pbs reads at a cold start is present.
[ -d "$MESH/0.orig" ] || refuse "mesh case '$MESH' has no 0.orig/"
REQ=$(ls "$MESH/0.orig" 2>/dev/null | grep -v uniform)
[ -n "$REQ" ] || refuse "0.orig in '$MESH' is empty; solve_hpc12.pbs derives the solved
         field set from it and refuses an empty one"
[ -d "$MESH/system" ] || refuse "mesh case '$MESH' has no system/"
CD="$MESH/system/controlDict"
[ -f "$CD" ] || refuse "mesh case '$MESH' has no system/controlDict"
grep -q "^startFrom *latestTime" "$CD" \
    || refuse "system/controlDict is not 'startFrom latestTime'; solve_hpc12.pbs refuses
         that, because a restart would silently begin at 0"
CASE_END=$(sed -n 's/^endTime *\([0-9]\+\) *;.*/\1/p' "$CD" | head -1)
[ -n "$CASE_END" ] || refuse "cannot read an integer endTime from system/controlDict, the
         way solve_hpc12.pbs reads it"
[ -d "$MESH/constant" ] || refuse "mesh case '$MESH' has no constant/"

# INPUT GATE 7: the recipe record travels with the case.
for f in MANIFEST.sha256 RECIPE_FILES.json; do
    [ -f "$MESH/$f" ] || refuse "mesh case '$MESH' has no $f; a case meshed by
         mesh_hpc12.pbs from a recipe folder always carries it"
done

# INPUT GATE 8: every file to be copied is readable as a regular file.
walk "$MESH/0.orig" "0.orig/" 1
walk "$MESH/system" "system/" 1
for e in "$MESH"/constant/* "$MESH"/constant/.[!.]* "$MESH"/constant/..?*; do
    [ -e "$e" ] || [ -L "$e" ] || continue
    b=${e##*/}
    [ "$b" = polyMesh ] && continue
    if [ -L "$e" ] && [ ! -e "$e" ]; then WALK_BAD="$WALK_BAD
    dangling symlink: constant/$b"; continue; fi
    [ -L "$e" ] && WALK_LINKS=$((WALK_LINKS + 1))
    if [ -d "$e" ]; then
        if [ -r "$e" ] && [ -x "$e" ]; then walk "$e" "constant/$b/" 2
        else WALK_BAD="$WALK_BAD
    unreadable directory: constant/$b"
        fi
    elif [ -f "$e" ] && [ -r "$e" ]; then WALK_OUT="${WALK_OUT}constant/$b
"
    elif [ -f "$e" ]; then WALK_BAD="$WALK_BAD
    unreadable file: constant/$b"
    else WALK_BAD="$WALK_BAD
    neither a regular file nor a directory: constant/$b"
    fi
done
for f in MANIFEST.sha256 RECIPE_FILES.json CONDITION.json; do
    if [ -L "$MESH/$f" ] && [ ! -e "$MESH/$f" ]; then WALK_BAD="$WALK_BAD
    dangling symlink: $f"; continue; fi
    [ -e "$MESH/$f" ] || continue
    [ -L "$MESH/$f" ] && WALK_LINKS=$((WALK_LINKS + 1))
    if [ -f "$MESH/$f" ] && [ -r "$MESH/$f" ]; then WALK_OUT="$WALK_OUT$f
"
    else WALK_BAD="$WALK_BAD
    not a readable regular file: $f"
    fi
done
DICTLIST="$WALK_OUT"
WALK_OUT=""
[ -L "$PM" ] && WALK_LINKS=$((WALK_LINKS + 1))
walk "$PM" "constant/polyMesh/" 2
MESHLIST="$WALK_OUT"
SRC_LINKS=$WALK_LINKS
[ -z "$WALK_BAD" ] || refuse "the mesh case holds entries that cannot be copied as files:$WALK_BAD"
N_DICT=$(nlines "$DICTLIST")
N_MESH=$(nlines "$MESHLIST")

# bytes to copy, from the files themselves
NEED=0
set -f
IFS_SAVE="$IFS"; IFS='
'
for rel in $DICTLIST $MESHLIST; do
    s=$(stat -L -c%s -- "$MESH/$rel" 2>/dev/null) || { IFS="$IFS_SAVE"; set +f; refuse "cannot stat $rel"; }
    NEED=$((NEED + s))
done
IFS="$IFS_SAVE"; set +f

# INPUT GATE 9: somewhere to put it. The nearest existing ancestor must be a writable
# directory with room for the copy, so the copy does not die half way on a full disk.
# A filesystem quota is not visible to df, so this is necessary, not sufficient.
ANC=$(dirname -- "$DST_ABS")
while [ ! -e "$ANC" ]; do ANC=$(dirname -- "$ANC"); done
[ -d "$ANC" ] || refuse "'$ANC' is not a directory, so '$DST' cannot be created under it"
[ -w "$ANC" ] || refuse "no write permission in '$ANC'"
DFL=$(df -Pk -- "$ANC" 2>/dev/null | sed -n '2p')
set -f
# shellcheck disable=SC2086
set -- $DFL
set +f
AVAIL_K="${4:-}"
NEED_K=$(( (NEED + 1023) / 1024 ))
case "$AVAIL_K" in
    ''|*[!0-9]*) note "could not read the free space under '$ANC'; not checked" ;;
    *) [ "$AVAIL_K" -ge "$NEED_K" ] \
           || refuse "the copy needs $NEED_K KiB and '$ANC' has $AVAIL_K KiB free" ;;
esac

# ---------------------------------------------------------------- the plan, stated
NC_SRC=$(head -c 4096 "$PM/owner" 2>/dev/null | tr -d '\000' | grep -ao 'nCells:[0-9]*' | head -1)
say "mesh case:  $MESH"
say "solve case: $DST"
say "mesh:       ${NC_SRC:-nCells:unknown}, $N_MESH files in constant/polyMesh, $NEED bytes to copy in all"
[ "$SRC_LINKS" -gt 0 ] && note "$SRC_LINKS symlink(s) in the mesh case; their TARGETS are copied, never the links"

# EVERY TOP-LEVEL ENTRY LANDS IN EXACTLY ONE BUCKET, copied or excluded with a reason,
# because an assertion over a set must partition the set.
TOP_N=0; TOP_COPIED=0; TOP_EXCL=0; EXCL_JSON=""; EXCL_PROC=0
say "top level of the mesh case:"
for e in "$MESH"/* "$MESH"/.[!.]* "$MESH"/..?*; do
    [ -e "$e" ] || [ -L "$e" ] || continue
    b=${e##*/}
    TOP_N=$((TOP_N + 1))
    case "$b" in
        0.orig|system|constant|MANIFEST.sha256|RECIPE_FILES.json|CONDITION.json)
            TOP_COPIED=$((TOP_COPIED + 1))
            if [ "$b" = constant ]; then
                say "  copied        constant/ (polyMesh verified file by file; the rest by sha256)"
            else
                say "  copied        $b"
            fi
            continue ;;
        processor*)
            TOP_EXCL=$((TOP_EXCL + 1)); EXCL_PROC=$((EXCL_PROC + 1)); continue ;;
        0)
            why="the solve job regenerates 0 from 0.orig on a cold start" ;;
        log)
            why="the mesh job's logs stay with the mesh case" ;;
        *)
            why="not part of a solve case; left in the mesh case" ;;
    esac
    TOP_EXCL=$((TOP_EXCL + 1))
    say "  NOT copied    $(printf '%q' "$b"): $why"
    EXCL_JSON="$EXCL_JSON${EXCL_JSON:+, }$(jstr "$b")"
done
if [ "$EXCL_PROC" -gt 0 ]; then
    say "  NOT copied    processor* ($EXCL_PROC directories): the mesh job's decomposition;"
    say "                the first solve leg runs decomposePar at $NR ranks"
    EXCL_JSON="$EXCL_JSON${EXCL_JSON:+, }$(jstr "processor* ($EXCL_PROC)")"
fi
# The case statement's last branch catches every other name, so the partition holds by
# construction; the counts are printed, not asserted, since that assertion could not fail.
say "  $TOP_N entries = $TOP_COPIED copied + $TOP_EXCL excluded"
if [ -d "$MESH/0" ]; then
    note "the mesh case holds a leftover 0/. It is NOT copied: solve_hpc12.pbs creates 0/"
    note "from 0.orig/ on a cold start, and only when 0/ is absent, so a copied 0/ would run"
    note "whatever it holds instead of 0.orig."
fi

# ---------------------------------------------------------------- build
STAMP=$(date -u '+%Y%m%dT%H%M%SZ')
STAGE_NAME="${DST_ABS}.incomplete_${STAMP}"
if [ -e "$STAGE_NAME" ] || [ -L "$STAGE_NAME" ]; then
    refuse "'$STAGE_NAME' already exists"
fi
# Past this point a refusal can no longer promise that nothing was created: mkdir -p
# may have made some of the missing parents before failing.
mkdir -p -- "$(dirname -- "$DST_ABS")" \
    || { echo "REFUSED: cannot create the parent of '$DST'" >&2; exit 2; }
mkdir -- "$STAGE_NAME" \
    || { echo "REFUSED: cannot create '$STAGE_NAME'" >&2; exit 2; }
STAGE="$STAGE_NAME"
say "staging in $STAGE"

cpl() { cp -RL --preserve=mode,timestamps -- "$1" "$2" || die "copy failed: $1 -> $2"; }
cpl "$MESH/0.orig" "$STAGE/0.orig"
cpl "$MESH/system" "$STAGE/system"
mkdir -- "$STAGE/constant" || die "cannot create $STAGE/constant"
for e in "$MESH"/constant/* "$MESH"/constant/.[!.]* "$MESH"/constant/..?*; do
    [ -e "$e" ] || continue
    b=${e##*/}
    [ "$b" = polyMesh ] && continue
    cpl "$e" "$STAGE/constant/$b"
done
for f in MANIFEST.sha256 RECIPE_FILES.json CONDITION.json; do
    [ -e "$MESH/$f" ] && cpl "$MESH/$f" "$STAGE/$f"
done
say "copying constant/polyMesh ($N_MESH files; the slow part)"
S0=$SECONDS
mkdir -- "$STAGE/constant/polyMesh" || die "cannot create $STAGE/constant/polyMesh"
cpl "$PM/." "$STAGE/constant/polyMesh/"
say "mesh copied in $((SECONDS - S0)) s"

# ---------------------------------------------------------------- decomposeParDict
# DETERMINISTIC: the text depends on the rank count alone (no date, no user, no path),
# so every case staged at the same count carries a byte-identical file and one sha256.
# The comments never carry the subdomain-count keyword followed by a space (older
# tooling rewrote every such line, comments included) and never a comment delimiter.
DP="$STAGE/system/decomposeParDict"
{
    cat <<'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// Solve-case decomposition, written by hpc/stage_solve_case.sh. The mesh case's own
// decomposition, the one the mesh was built with, is recorded in SOLVE_SETUP.json.
// MANIFEST.sha256 describes the recipe as meshed, so its entry for this file is
// expected not to match; SOLVE_SETUP.json records whether it does.
//
EOF
    case "$NR" in
        128) cat <<'EOF'
// RULE APPLIED: 128 ranks gives hierarchical (8 2 8). This is the configuration of the
// six wake-refined Condition CR solves (the r2 cases B, C, F, H, M, W). Read back on
// HPC12, every one of their solve cases carries exactly this: 128 subdomains, method
// hierarchical, n (8 2 8), delta 0.001, order xyz. 8 by 2 by 8 is 128 exactly.
// WHY 128: fpt-medium takes at most 4 nodes of 32 cores per job and up to three jobs
// per user; fpt-large takes one job per user.
EOF
            ;;
        192) cat <<'EOF'
// RULE APPLIED: 192 ranks gives hierarchical (8 3 8), the case definitions of the
// published 192-rank solves. It is NOT the r2 configuration, which solved at 128.
// 8 by 3 by 8 is 192 exactly.
EOF
            ;;
        *) cat <<'EOF'
// RULE APPLIED: a rank count other than 128 or 192 gives method scotch, which splits
// into any count. THIS DEPARTS FROM THE DELIVERED CONFIGURATION: the r2 solves ran at
// 128 (hierarchical 8 2 8) and the published solve definitions carry 192 (hierarchical
// 8 3 8). A result on a different decomposition is a result on a different partition.
// scotch needs an OpenFOAM build that carries the real scotch library.
EOF
            ;;
    esac
    echo ""
    echo "numberOfSubdomains $NR;"
    printf 'method          %s;\n' "$METHOD"
    if [ "$METHOD" = hierarchical ]; then
        echo ""
        echo "hierarchicalCoeffs"
        echo "{"
        echo "    n               ($NX $NY $NZ);"
        echo "    delta           0.001;"
        echo "    order           xyz;"
        echo "}"
    fi
} > "$DP" || die "cannot write $DP"

# ---------------------------------------------------------------- verify, from the files
say "verifying"
# (a) the decomposition, read back with solve_hpc12.pbs's own expression (its line 109)
RB_N=$(sed -n 's/^ *numberOfSubdomains *\([0-9]\+\) *;.*/\1/p' "$DP" 2>/dev/null | head -1)
RB_M=$(sed -n 's/^ *method *\([A-Za-z]\+\) *;.*/\1/p' "$DP" | head -1)
RB_NV=$(sed -n 's/^ *n *( *\([0-9]\+\) \+\([0-9]\+\) \+\([0-9]\+\) *) *;.*/\1 \2 \3/p' "$DP" | head -1)
[ "$RB_N" = "$NR" ] || die "decomposeParDict reads back as '$RB_N' subdomains, not $NR"
[ "$RB_M" = "$METHOD" ] || die "decomposeParDict reads back method '$RB_M', not $METHOD"
KW=$(grep -c 'numberOfSubdomains ' "$DP")
[ "$KW" -eq 1 ] || die "the subdomain-count keyword plus a space appears $KW times, not once"
DELIM=$(grep -c '/\*\|\*/' "$DP")
[ "$DELIM" -eq 0 ] || die "decomposeParDict carries a comment delimiter"
if [ "$METHOD" = hierarchical ]; then
    [ -n "$RB_NV" ] || die "cannot read hierarchicalCoeffs n back"
    set -f
    # shellcheck disable=SC2086
    set -- $RB_NV
    set +f
    [ "$1 $2 $3" = "$NX $NY $NZ" ] || die "n reads back as ($RB_NV), not ($NX $NY $NZ)"
    [ $(($1 * $2 * $3)) -eq "$NR" ] || die "n ($RB_NV) multiplies to $(($1 * $2 * $3)), not $NR"
    N_JSON="[$1, $2, $3]"
    say "  decomposeParDict: $RB_N subdomains, $RB_M ($RB_NV), product asserted"
else
    [ -z "$RB_NV" ] || die "a scotch dictionary should carry no hierarchicalCoeffs n"
    N_JSON="null"
    say "  decomposeParDict: $RB_N subdomains, $RB_M"
fi
DP_SHA=$(sha "$DP")

# (b) every file: the staged set equals the source set, dictionaries by sha256, mesh
#     files by size and by sha256 up to 64 MiB. The buckets must sum to the count.
set -f
IFS_SAVE="$IFS"; IFS='
'
# same_sha: true only when BOTH digests were computed and agree. Two failed reads would
# otherwise yield two empty strings, which compare equal and vouch for a copy nobody read.
same_sha() {
    local h1 h2
    h1=$(sha "$1"); h2=$(sha "$2")
    [ ${#h1} -eq 64 ] && [ "$h1" = "$h2" ]
}
D_OK=0; D_BAD=""; D_REPL=0
for rel in $DICTLIST; do
    if [ "$rel" = "system/decomposeParDict" ]; then D_REPL=1; continue; fi
    if same_sha "$MESH/$rel" "$STAGE/$rel"; then D_OK=$((D_OK + 1)); else D_BAD="$D_BAD $rel"; fi
done
M_HASH=0; M_SIZE=0; M_BAD=""
for rel in $MESHLIST; do
    s1=$(stat -L -c%s -- "$MESH/$rel" 2>/dev/null); s2=$(stat -c%s -- "$STAGE/$rel" 2>/dev/null)
    if [ -z "$s2" ] || [ "$s1" != "$s2" ]; then
        M_BAD="$M_BAD $rel(src=$s1,dst=${s2:-absent})"; continue
    fi
    if [ "$s1" -le "$HASH_MAX" ]; then
        if same_sha "$MESH/$rel" "$STAGE/$rel"; then M_HASH=$((M_HASH + 1)); else M_BAD="$M_BAD $rel(sha256)"; fi
    else
        M_SIZE=$((M_SIZE + 1))
    fi
done
IFS="$IFS_SAVE"; set +f
[ -z "$D_BAD" ] || die "copied dictionary files differ from the mesh case:$D_BAD"
[ -z "$M_BAD" ] || die "mesh files did not copy identically:$M_BAD"
# The counts below come from grep -c over the list, the loops from IFS word-splitting of
# the same list: two parsings, so the sum is a real cross-check of what was adjudicated.
[ $((D_OK + D_REPL)) -eq "$N_DICT" ] || die "dictionary partition broken: $N_DICT = $D_OK + $D_REPL"
[ $((M_HASH + M_SIZE)) -eq "$N_MESH" ] || die "mesh partition broken: $N_MESH = $M_HASH + $M_SIZE"
say "  mesh: $N_MESH files = $M_HASH identical by sha256 + $M_SIZE identical by size (over 64 MiB)"
say "  other files: $N_DICT = $D_OK identical by sha256 + $D_REPL decomposeParDict replaced"

WALK_OUT=""; WALK_BAD=""; WALK_LINKS=0
walk "$STAGE" "" 1
GOT=$(printf '%s' "$WALK_OUT" | sort)
WANT=$(printf '%s%ssystem/decomposeParDict\n' "$DICTLIST" "$MESHLIST" | sort -u)
[ -z "$WALK_BAD" ] || die "the staged case holds unexpected entries:$WALK_BAD"
[ "$WALK_LINKS" -eq 0 ] || die "the staged case holds $WALK_LINKS symlink(s); it must hold copies only"
[ "$GOT" = "$WANT" ] || die "the staged file set is not the source file set plus decomposeParDict"
say "  staged file set equals the source set; no symlinks"

# (c) what solve_hpc12.pbs checks at a cold start, read from the staged case
SREQ=$(ls "$STAGE/0.orig" 2>/dev/null | grep -v uniform)
[ -n "$SREQ" ] || die "staged 0.orig is empty"
for f in $SREQ; do
    head -c 64 "$STAGE/0.orig/$f" > /dev/null 2>&1 || die "staged 0.orig/$f is unreadable"
done
for x in 0 processor0 log .converged; do
    [ -e "$STAGE/$x" ] && die "staged case holds '$x', which would not start cold"
done
grep -q "^startFrom *latestTime" "$STAGE/system/controlDict" || die "staged controlDict is not startFrom latestTime"
S_END=$(sed -n 's/^endTime *\([0-9]\+\) *;.*/\1/p' "$STAGE/system/controlDict" | head -1)
S_WI=$(sed -n 's/^writeInterval *\([0-9]\+\) *;.*/\1/p' "$STAGE/system/controlDict" | head -1)
[ -n "$S_END" ] || die "cannot read endTime from the staged controlDict"
say "  cold start: 0.orig holds $(echo "$SREQ" | wc -w) fields ($(echo "$SREQ" | tr '\n' ' '| sed 's/ $//')); no 0/, processor*, log/ or .converged"

NC=$(head -c 4096 "$STAGE/constant/polyMesh/owner" 2>/dev/null | tr -d '\000' | grep -ao 'nCells:[0-9]*' | head -1 | sed 's/nCells://')
case "$NC" in ''|*[!0-9]*) NC_JSON="null"; note "nCells not found in the owner header" ;; *) NC_JSON="$NC" ;; esac
[ "$NC_JSON" != null ] && say "  nCells $NC (owner header)"

# ---------------------------------------------------------------- SOLVE_SETUP.json
MDP="$MESH/system/decomposeParDict"
if [ -f "$MDP" ]; then
    M_N=$(sed -n 's/^ *numberOfSubdomains *\([0-9]\+\) *;.*/\1/p' "$MDP" | head -1)
    M_M=$(sed -n 's/^ *method *\([A-Za-z]\+\) *;.*/\1/p' "$MDP" | head -1)
    MESH_DP_JSON="{\"subdomains\": ${M_N:-null}, \"method\": $( [ -n "$M_M" ] && jstr "$M_M" || echo null ), \"sha256\": $(jstr "$(sha "$MDP")")}"
else
    MESH_DP_JSON="null"
fi
MAN_DP=$(grep '  system/decomposeParDict$' "$STAGE/MANIFEST.sha256" | head -1 | cut -d' ' -f1)
if [ -n "$MAN_DP" ]; then
    [ "$MAN_DP" = "$DP_SHA" ] && MAN_MATCH=true || MAN_MATCH=false
    MAN_JSON="{\"sha256\": $(jstr "$MAN_DP"), \"matches_written_file\": $MAN_MATCH}"
else
    MAN_JSON="null"
fi
SELF_SHA=$(sha "${BASH_SOURCE[0]}")
NOW=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

SJ="$STAGE/SOLVE_SETUP.json"
cat > "$SJ.part" <<EOF || die "cannot write $SJ.part"
{
  "written_by": "hpc/stage_solve_case.sh",
  "script_sha256": $(jstr "${SELF_SHA:-unknown}"),
  "date_utc": "$NOW",
  "mesh_case": $(jstr "$MESH"),
  "mesh_case_absolute": $(jstr "$MESH_ABS"),
  "solve_case": $(jstr "$DST"),
  "ranks": $RB_N,
  "ranks_source": $(jstr "$RSRC"),
  "rule": $(jstr "$RULE"),
  "method": $(jstr "$RB_M"),
  "n": $N_JSON,
  "decomposeParDict_sha256": $(jstr "$DP_SHA"),
  "nCells": $NC_JSON,
  "mesh_case_decomposeParDict": $MESH_DP_JSON,
  "manifest_entry_decomposeParDict": $MAN_JSON,
  "files_verified": {"mesh_sha256": $M_HASH, "mesh_size_only": $M_SIZE, "other_sha256": $D_OK},
  "not_copied": [$EXCL_JSON]
}
EOF
mv -- "$SJ.part" "$SJ" || die "cannot promote $SJ"
grep -q "^  \"ranks\": $NR,\$" "$SJ" || die "SOLVE_SETUP.json does not read back ranks $NR"

# ---------------------------------------------------------------- promote
if [ -e "$DST_ABS" ] || [ -L "$DST_ABS" ]; then
    die "'$DST' appeared while staging; not overwriting it"
fi
mv -T -- "$STAGE" "$DST_ABS" || die "cannot rename $STAGE to $DST_ABS"
{ [ -f "$DST_ABS/SOLVE_SETUP.json" ] && [ ! -e "$STAGE" ]; } || die "the rename did not land"
STAGE=""
say "STAGED: $DST_ABS ($(du -sh -- "$DST_ABS" 2>/dev/null | cut -f1))"

# ---------------------------------------------------------------- next steps
# SCRIPT_DIR, SOLVE_PBS and GATE_PY were resolved and gated with the inputs (gate 2a). A
# file missing from the clone is printed as a placeholder to fill in, never as a path
# that does not exist.
if [ ! -f "$SOLVE_PBS" ]; then
    note "$SOLVE_PBS not found; use the absolute path of hpc/solve_hpc12.pbs in your clone"
    SOLVE_PBS="<absolute path of hpc/solve_hpc12.pbs>"
fi
if [ ! -f "$GATE_PY" ]; then
    note "$GATE_PY not found; use the absolute path of scripts/assert_case_postpro.py in your clone"
    GATE_PY="<absolute path of scripts/assert_case_postpro.py>"
fi
INSTALL_PY="$REPO_DIR/scripts/install_postpro.py"
[ -f "$INSTALL_PY" ] || INSTALL_PY="<absolute path of scripts/install_postpro.py>"
JOB="argus-$(basename -- "$DST_ABS")"
VLIST="ARGUS_CASE=$DST_ABS,ARGUS_SCRIPT=$SOLVE_PBS,ARGUS_GATE=$GATE_PY"

# THE BLOCK COUNTS AS PRESENT ONLY WHEN BOTH HALVES ARE: a live include line in
# controlDict and the file it names. assert_case_postpro.py --before pulls an included
# file in only if it exists, so an include without the file fails the gate as surely as
# no include at all, and a mention inside a comment is not an include.
PP_INC=0; PP_FILE=0
grep -Eq '^[[:space:]]*#include[[:space:]]+"argusPostPro"' "$DST_ABS/system/controlDict" && PP_INC=1
[ -f "$DST_ABS/system/argusPostPro" ] && PP_FILE=1

echo
echo "NEXT STEPS"
echo "1. Post-processing block."
if [ "$PP_INC" -eq 1 ] && [ "$PP_FILE" -eq 1 ]; then
    echo "   system/controlDict already includes argusPostPro."
    echo "   Check it here before you submit. The check needs only the Python standard"
    echo "   library, and it must end \"---- 8 of 8 requirements satisfied ----\":"
    echo "       python3 $GATE_PY --before $DST_ABS"
else
    if [ "$PP_INC" -eq 1 ]; then
        echo "   WARNING: system/controlDict includes argusPostPro, but the file"
        echo "   system/argusPostPro is missing."
    else
        echo "   WARNING: system/controlDict does NOT include argusPostPro: the mesh case"
        echo "   was built without the block."
    fi
    echo "   Without the block this case does not write the pressure and viscous forces,"
    echo "   the wing-surface fields or the Trefftz planes, and with ARGUS_GATE in the -v"
    echo "   list (step 4) the solve job refuses every leg. Install it now, before the"
    echo "   first leg, here at the HPC12 login prompt:"
    echo "       python3 $INSTALL_PY $DST_ABS"
    echo "       python3 $GATE_PY --before $DST_ABS"
    echo "   The check must end \"---- 8 of 8 requirements satisfied ----\"."
    echo "   install_postpro.py needs numpy, which the HPC12 login node has (1.12.1). Its"
    echo "   writing step has not been run on HPC12: NOT VERIFIED."
    echo "   It reads system/controlDict, constant/polyMesh/boundary (and processor0's,"
    echo "   if the first names no wall patch), 0.orig/U, constant/momentumTransport"
    echo "   and the ASCII STL in constant/triSurface; for a compressible (cruise) case"
    echo "   also 0.orig/p, 0.orig/T and constant/physicalProperties (or"
    echo "   thermophysicalProperties); and the template hpc/templates/argusPostPro from"
    echo "   the repository. Once a leg has run, it reads 0/U in place of 0.orig/U, and"
    echo "   tries a sampled wing surface under postProcessing/ before the STL. It writes"
    echo "   system/argusPostPro and adds the include line to system/controlDict, where"
    echo "   it also removes any wingSurface, wallShearStress, yPlus or forces1 function"
    echo "   object, because the block replaces them."
    echo "   If the install fails here, run it on your own machine's copy of the case"
    echo "   and copy only that copy's system/argusPostPro into"
    echo "   $DST_ABS/system/."
    if [ "$PP_INC" -eq 1 ]; then
        echo "   The include line is already in system/controlDict."
    else
        echo "   Add this line directly under the { that follows functions in"
        echo "   $DST_ABS/system/controlDict:"
        echo "       #include \"argusPostPro\""
    fi
    echo "   Then repeat the check (NOT VERIFIED on HPC12)."
fi
echo "2. Angle of attack. The case carries the angle the mesh case was given:"
grep -E '^Uinf' "$DST_ABS/0.orig/U" 2>/dev/null | sed 's/^/       /'
grep -E '^ *(liftDir|dragDir)' "$DST_ABS/system/controlDict" 2>/dev/null | sed 's/^ */       /'
echo "   To run another angle, set it as in docs/06_trim_to_target_cl.md, section 6.4.2"
echo "   (a cold start), and check it with section 6.4.4."
echo "3. Iterations. system/controlDict: endTime $S_END, writeInterval ${S_WI:-unread}. The r2"
echo "   solves ran their cold leg to endTime 4000 with writeInterval 500. This script"
echo "   changes neither."
if [ -n "$S_WI" ] && [ "$S_WI" -gt 0 ] && [ $((S_END % S_WI)) -ne 0 ]; then
    echo "   WARNING: endTime $S_END is not a multiple of writeInterval $S_WI. solve_hpc12.pbs"
    echo "   refuses that unless ARGUS_ENDTIME in the -v list is a multiple of $S_WI."
fi
echo "4. Submit ($NR ranks) from the HPC12 login node, from inside the solve case: PBS"
echo "   copies the job's output file, $JOB.o<job id>, into the folder qsub runs from when"
echo "   the job ends. While it runs, HPC12 keeps that text in \$HOME/<job id>.hpc12.hpc.OU."
echo "       cd $DST_ABS"
if [ "$NR" -eq 128 ]; then
    echo "       qsub -N $JOB -q fpt-medium -l nodes=4:ppn=32:typej,walltime=72:00:00 \\"
    echo "            -v $VLIST \\"
    echo "            $SOLVE_PBS"
    echo "   Faster, if you have no other fpt-large job (fpt-large holds one job per user),"
    echo "   in place of that qsub line:"
    echo "       qsub -N $JOB -q fpt-large -l nodes=1:ppn=128:typen,walltime=72:00:00 \\"
    echo "            -v $VLIST \\"
    echo "            $SOLVE_PBS"
    echo "   The r2 legs ran at about 18.2 to 19.9 s per iteration on 4 typej hosts and"
    echo "   about 10.1 s per iteration on one typen host (solver ClockTime)."
else
    echo "   No -l line for $NR ranks was used in the delivered runs. The allocation must"
    echo "   provide at least $NR slots: solve_hpc12.pbs refuses fewer, and idles any surplus."
    echo "       qsub -N $JOB -q <queue> -l <resources giving at least $NR slots>,walltime=72:00:00 \\"
    echo "            -v $VLIST \\"
    echo "            $SOLVE_PBS"
fi
echo "   ARGUS_GATE makes solve_hpc12.pbs run the step 1 check at the start of every leg,"
echo "   before the write probe, decomposePar and foamRun, and stop a leg that fails it"
echo "   with \"REFUSED: post-processing gate refused\". Without ARGUS_GATE the job looks"
echo "   for \$HOME/argus_hpc12/assert_case_postpro.py, which does not exist on a new"
echo "   account, and runs the leg unchecked. Once the job has started:"
echo "       grep -E \"requirements satisfied|gate NOT run|gate refused\" \$HOME/<job id>.hpc12.hpc.OU"
echo "   should print \"---- 8 of 8 requirements satisfied ----\". After the job ends, grep"
echo "   $JOB.o<job id> in the solve case instead."
echo "   Add ARGUS_OF_BASHRC=<path> to the -v list if OpenFOAM-12 is not in"
echo "   \$HOME/OpenFOAM/OpenFOAM-12."
echo "5. Later legs. Submit every leg yourself, from the same folder, with the same -v"
echo "   list, ARGUS_GATE included. A job's own qsub for its next leg fails on HPC12"
echo "   (could not connect to trqauthd), and would not pass ARGUS_GATE on anyway."
echo "   Add ARGUS_ENDTIME where docs/05_submit_and_run.md, section 10a, says so. It"
echo "   must be a multiple of writeInterval (${S_WI:-unread}), or solve_hpc12.pbs refuses the leg."
exit 0
