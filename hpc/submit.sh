#!/usr/bin/env bash
# hpc/submit.sh -- assemble a job script at SUBMISSION time and submit it, with
# every run-time-read file hash-pinned into the case.
#
# Usage:
#   hpc/submit.sh <template> <case-dir> [extra sbatch args...]
#
#   <template>  either a .slurm (submitted as-is) or a .header (assembled with
#               the body named in its `# ARGUS-BODY:` line)
#
# =============================================================================
# WHAT THIS EXISTS TO CLOSE: THE GAP BETWEEN SNAPSHOT AND RESOLUTION
# =============================================================================
# `sbatch` COPIES the submitted script into SLURM's spool at submission, which is
# what makes a queued job reproducible: editing the template afterwards cannot
# reach it. VERIFIED, not assumed, with `scontrol write batch_script <id> -`.
#
# ANYTHING THE SCRIPT READS AT RUN TIME ESCAPES THAT GUARANTEE. It is read from
# the filesystem when the job STARTS, which on a deep queue is days later. Three
# consequences, in increasing order of nastiness:
#   1. A `source` of a sibling file does not even resolve, because the spool copy
#      has no siblings. Leg B died in 3 seconds on exactly this.
#   2. If it did resolve, an edit between submission and start would silently
#      change what a queued job executes.
#   3. THE INSTRUMENT CANNOT SEE IT. `scontrol write batch_script` returns the
#      source line, not the body. So the one check you would use to confirm what
#      a queued job will do returns a POINTER, and the pointer resolves at a
#      different time than the check ran.
#
# Two fixes, and they are different fixes for different things:
#   ASSEMBLY   the body is concatenated at submission, so it is inside the
#              snapshot and `scontrol write batch_script` returns the real bytes.
#   PINNING    what genuinely cannot be assembled -- env.sh, which every template
#              sources, and the python the jobs invoke -- is HASHED at submission
#              into the case, and re-checked by the job at start. An edit then
#              fails LOUDLY instead of substituting silently.
#
# Pinning is strictly weaker than assembly and is used only where assembly is not
# possible. It detects drift; it does not prevent it.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"

die() { echo "submit: FATAL: $*" >&2; exit 1; }
note() { echo "submit: $*"; }

# --dry-run assembles and validates but does not submit, so the assembly itself is
# testable on a machine with no SLURM. An assembler that can only be exercised by
# submitting a job is an assembler nothing checks.
DRY=0
if [ "${1:-}" = "--dry-run" ]; then DRY=1; shift; fi
[ $# -ge 2 ] || die "usage: hpc/submit.sh [--dry-run] <template> <case-dir> [sbatch args...]"
TEMPLATE="$1"; shift
CASE="$1"; shift

[ -f "$TEMPLATE" ] || die "no such template: ${TEMPLATE}"
[ -d "$CASE" ] || die "no such case directory: ${CASE}"

# ------------------------------------------------------------- 1. assemble ----
ASSEMBLED="$(mktemp "${TMPDIR:-/tmp}/argus_submit.XXXXXX.slurm")"
BODY=""

# env.sh IS INLINED, NOT PINNED. It is sourced by every template and read at JOB
# START, so it has been outside the sbatch snapshot for EVERY JOB THIS PROJECT HAS
# EVER RUN, including job 1. The sourced body enlarged that exposure and made it
# visible; it did not create it.
#
# Leg B failed LOUDLY only by luck of path construction: ${BASH_SOURCE[0]} pointed
# at SLURM's spool, which has no siblings, so nothing resolved. An ABSOLUTE path,
# or one relative to $SLURM_SUBMIT_DIR, resolves perfectly well from a compute
# node -- and that is the SILENT variant. env.sh is written exactly that way.
#
# Inlining beats pinning wherever it is possible: pinning DETECTS drift, inlining
# leaves nothing to drift. One mechanism instead of two, and the single definition
# survives because the inlined text is assembled from the one file.
#
# INSERTED AT THE DIRECTIVE BOUNDARY, which is the first line that is neither
# blank nor a comment. SLURM reads #SBATCH only until the first executable line,
# so inserting there keeps every directive above it effective. Nothing is
# filtered or rewritten: this is an INSERTION, and the template's own bytes are
# all present in order.
insert_env() {
    local src="$1" out="$2"
    awk -v envf="${REPO}/hpc/env.sh" '
        BEGIN { done = 0 }
        {
            if (!done && $0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/) {
                print "# ---- BEGIN INLINED hpc/env.sh (assembled by hpc/submit.sh) ----"
                print "ARGUS_ENV_INLINED=1"
                while ((getline line < envf) > 0) print line
                close(envf)
                print "# ---- END INLINED hpc/env.sh ----"
                done = 1
            }
            print
        }
        END {
            if (!done) {
                print "# ---- BEGIN INLINED hpc/env.sh (assembled by hpc/submit.sh) ----"
                print "ARGUS_ENV_INLINED=1"
                while ((getline line < envf) > 0) print line
                close(envf)
                print "# ---- END INLINED hpc/env.sh ----"
            }
        }
    ' "$src" > "$out"
}

if [ "${TEMPLATE##*.}" = "header" ]; then
    BODY="$(sed -n 's/^# ARGUS-BODY:[[:space:]]*//p' "$TEMPLATE" | head -1)"
    [ -n "$BODY" ] || die "${TEMPLATE} is a .header but declares no '# ARGUS-BODY:' line.
The pairing lives in the header on purpose: a header that does not say what body
it needs cannot be assembled correctly by anything, including a careful reader."
    [ -f "${REPO}/${BODY}" ] || die "declared body ${BODY} does not exist under ${REPO}"
    insert_env "$TEMPLATE" "${ASSEMBLED}.hdr"
    cat "${ASSEMBLED}.hdr" "${REPO}/${BODY}" > "$ASSEMBLED"
    rm -f "${ASSEMBLED}.hdr"
    note "assembled ${TEMPLATE} + inlined env.sh + ${BODY}"
else
    insert_env "$TEMPLATE" "$ASSEMBLED"
    note "assembled ${TEMPLATE} + inlined env.sh"
fi

# DERIVED: the assembled file must actually contain env.sh, not merely have been
# through a function named insert_env. The assembling step does not get to vouch
# for its own result.
grep -q "END INLINED hpc/env.sh" "$ASSEMBLED" || die "env.sh was NOT inlined into the assembled script"
grep -q "^argus_assert_runtime_pins" "$ASSEMBLED" || die "the inlined env.sh does not define argus_assert_runtime_pins; the wrong file was inlined"

ASM_SHA="$(sha256sum "$ASSEMBLED" | cut -d' ' -f1)"
note "assembled sha256 ${ASM_SHA}"

# --------------------------------------------------------------- 2. pin -------
# WHAT A QUEUED JOB READS AT RUN TIME, ENUMERATED RATHER THAN REMEMBERED. Every
# entry here is outside the sbatch snapshot. env.sh has ALWAYS been outside it,
# for every job this project has ever submitted including job 1 -- the sourced
# body made the exposure larger, it did not create it.
# env.sh IS NO LONGER HERE. It is inlined above, so there is nothing left to
# drift. Pinning is the fallback for what cannot be inlined: a python file cannot
# sensibly be pasted into a bash script, so it keeps a pin, and the pin is checked
# BEFORE it is invoked.
RUNTIME_READS=(
    "hpc/measure_mesh_surface.py"
)
PINFILE="${CASE}/RUNTIME_PINS.sha256"
: > "${PINFILE}.tmp"
for f in "${RUNTIME_READS[@]}"; do
    if [ -f "${REPO}/${f}" ]; then
        printf '%s  %s\n' "$(sha256sum "${REPO}/${f}" | cut -d' ' -f1)" "$f" >> "${PINFILE}.tmp"
    else
        printf 'MISSING  %s\n' "$f" >> "${PINFILE}.tmp"
    fi
done
mv "${PINFILE}.tmp" "$PINFILE"
note "pinned $(wc -l < "$PINFILE") run-time-read files into ${PINFILE}"

# The assembled script's own hash is recorded beside the pins so a run card can
# state what actually ran. DERIVED from the file, not from the intent.
#
# THE ENVIRONMENT IS RECORDED TOO, AND IT IS THE HALF THAT WAS MISSING (2026-08-21).
# The script hash pins the INVARIANT part. Every CONSEQUENTIAL setting is in the
# ENVIRONMENT: ARGUS_ENDTIME decides how far the run goes, ARGUS_WALL decides the wall
# treatment (off-limits to get wrong), ARGUS_SKIP_MESH and ARGUS_RESTART decide whether
# the mesh and the existing times survive. A record of the script alone therefore says
# almost nothing about what the job will do.
#
# MEASURED COST OF THE OMISSION: jobs 10680519/20 (L11 pair) sat queued for 14 h each
# against cases whose controlDict read endTime 300 at latestTime 300, i.e. a guaranteed
# zero-iteration run. Whether ARGUS_ENDTIME had been pinned COULD NOT BE ESTABLISHED from
# anything on disk or in SLURM, because the only record was the script hash and the
# variable is read at line 1143 of the assembled body from an environment nobody kept.
# The diagnosis was impossible, not merely difficult, and that is a record defect rather
# than a debugging one.
_ARGUS_ENV_JSON="$(env | sed -n 's/^\(ARGUS_[A-Z_]*\)=\(.*\)$/    "\1": "\2",/p' | sed '$ s/,$//')"
cat > "${CASE}/SUBMITTED_SCRIPT.json" <<EOF
{
  "template": "${TEMPLATE}",
  "body": "${BODY:-none}",
  "assembled_sha256": "${ASM_SHA}",
  "note": "The assembled file is what sbatch snapshotted. Retrieve it with scontrol write batch_script <jobid> - and it must hash to assembled_sha256.",
  "sbatch_args": "$(printf '%s ' "$@" | sed 's/"/\\"/g; s/ *$//')",
  "argus_env": {
${_ARGUS_ENV_JSON}
  },
  "env_note": "ARGUS_* as seen by submit.sh. sbatch propagates these via --export; anything passed ONLY in --export=ALL,VAR=val on the command line appears in sbatch_args above, not here. Read BOTH.",
  "submitted_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF

# ------------------------------------------------------------- 3. submit ------
# The assembled script must PARSE. Catching a syntax error here costs nothing;
# catching it on a compute node costs a queue slot and a node-hour.
bash -n "$ASSEMBLED" || die "the assembled script does not parse"

if [ "$DRY" -eq 1 ]; then
    note "DRY RUN: assembled and validated, not submitted"
    note "  assembled sha256 ${ASM_SHA}"
    cp "$ASSEMBLED" "${ARGUS_DRYRUN_KEEP:-/dev/null}" 2>/dev/null || true
    rm -f "$ASSEMBLED"
    exit 0
fi

JOBID="$(sbatch --parsable "$@" "$ASSEMBLED")"
note "submitted job ${JOBID}"

# DERIVED, NOT ASSERTED BY THE SUBMITTING STEP. Read the snapshot back out of
# SLURM and confirm it hashes to what was submitted. If this ever disagrees, the
# assembly story is wrong and everything downstream of it is unsupported.
SPOOL_SHA="$(scontrol write batch_script "$JOBID" - 2>/dev/null | sha256sum | cut -d' ' -f1 || echo unavailable)"
if [ "$SPOOL_SHA" = "$ASM_SHA" ]; then
    note "VERIFIED: SLURM's snapshot hashes to the assembled file"
elif [ "$SPOOL_SHA" = "unavailable" ]; then
    echo "submit: WARNING: could not read the snapshot back; assembly is UNVERIFIED" >&2
else
    echo "submit: WARNING: snapshot hash ${SPOOL_SHA} != assembled ${ASM_SHA}." >&2
    echo "submit: SLURM may normalise trailing newlines. Compare before trusting." >&2
fi

rm -f "$ASSEMBLED"
echo "$JOBID"
