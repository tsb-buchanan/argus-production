#!/usr/bin/env bash
# hpc/push_to_cluster.sh -- LAPTOP SIDE. Put hpc/ on DelftBlue's /home and create
# the /scratch working directories. This is the prerequisite for job 0.
#
# THE USER RUNS THIS. It needs your credentials, and automation running it has none: no
# password is stored, read, or accepted by this script, and none should ever be
# typed into a chat or ticket. Authentication is whatever ssh already does
# for you, i.e. your SSH key or an interactive password prompt.
#
# UNTESTED against DelftBlue, like everything else in hpc/.
#
# Usage:
#   hpc/push_to_cluster.sh [--remote <sshhost>] [--dry-run] [--bastion <netid@host>]
#
#   --remote    ssh destination. Default 'delftblue', which assumes the
#               ~/.ssh/config block in README_HPC.md section 4. Without that
#               block, pass e.g. --remote <your-netid>@login.delftblue.tudelft.nl
#   --bastion   route via a jump host, for off-campus without EduVPN:
#               --bastion <netid>@linux-bastion.tudelft.nl
#   --dry-run   show what would transfer, transfer nothing
#   --check     report whether the cluster copy is current, then stop.
#               Exits 0 if in sync, 1 if stale or absent. Run this before
#               any `sbatch --test-only`: a header fixed after the last push
#               is invisible from the cluster side, and the test then
#               validates the OLD file while looking entirely normal.
#   --interactive  allow ssh to PROMPT you for a password instead of requiring a
#               key. Use only if you have not set up a key yet. The password is
#               typed at ssh's own prompt: this script never reads, stores or
#               accepts one, and none should ever be pasted into a chat, a
#               command line or a file.
#
# =============================================================================
# WHAT GETS TRANSFERRED, AND WHY IT IS ONLY hpc/
# =============================================================================
# 192 KB. Every cluster-side script sources ${REPO}/hpc/env.sh and nothing else
# outside hpc/. The analysis scripts (scripts/analyse_coverage_ladder.py and the
# rest) run on the LAPTOP after pull_results.sh, so they never need to be there.
#
# NOT the whole repo: the working tree minus cases/ and .git/ is 700 MB, almost
# all of it untracked geometry and results. /home is a 30.00 GiB hard quota with
# 25.02 GiB already used, i.e. under 5 GiB free, so 700 MB of files nothing on
# the cluster reads would be a seventh of the remaining space.
#
# NOT via git clone, which would otherwise be the clean route: as of this
# writing NOTHING IS PUSHED (every commit including hpc/ is local-only) and
# github.com/<owner>/argus-validation returns 404 unauthenticated, so a
# clone on the cluster would need a token as well as getting nothing useful.
# If that changes, `git clone` on the cluster is better than this script and
# makes the COMMIT stamping below unnecessary.
#
# NOT the case bundles. Those are 61 MB each, they go to /scratch, not /home,
# and package_case.sh has only been dry-run so far: none exists yet. See
# README_HPC.md sections 5 and 6.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"
# shellcheck source=hpc/env.sh
source "${HERE}/env.sh"

die() { echo "push_to_cluster: FATAL: $*" >&2; exit 1; }
note() { echo "push_to_cluster: $*"; }

# THE SILENT-DEATH TRAP. `set -e` exits with no message at all, and a script
# that vanishes leaves the operator with nothing to act on. This happened here:
# a failed command substitution inside an assignment killed the script after the
# last note() and before the next one, so the visible output simply stopped.
#
# Every unhandled non-zero status now names the line it died on. The explicit
# die() and FATAL paths below are still the right way to fail; this only ensures
# that the ones nobody anticipated are LOUD rather than invisible.
_argus_err_trap() {
    local rc=$1 line=$2
    echo >&2
    echo "push_to_cluster: ABORTED at line ${line} with status ${rc}." >&2
    echo "push_to_cluster: No diagnosis was printed before this, which means the" >&2
    echo "push_to_cluster: failure was NOT anticipated. That is a bug in this" >&2
    echo "push_to_cluster: script, not necessarily in your setup. Re-run with:" >&2
    echo "push_to_cluster:     bash -x $0 $*" >&2
    echo "push_to_cluster: and report the last few lines." >&2
    exit "$rc"
}
trap '_argus_err_trap $? $LINENO' ERR

REMOTE="delftblue"
BASTION=""
DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --remote)  [ $# -ge 2 ] || die "--remote needs a host"; REMOTE="$2"; shift 2 ;;
        --bastion) [ $# -ge 2 ] || die "--bastion needs netid@host"; BASTION="$2"; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        --check)   CHECK_ONLY=1; shift ;;
        --interactive) ARGUS_INTERACTIVE=1; export ARGUS_INTERACTIVE; shift ;;
        *) die "unknown option: $1" ;;
    esac
done
CHECK_ONLY="${CHECK_ONLY:-0}"

SSH_OPTS=()
RSYNC_SSH=()
if [ -n "$BASTION" ]; then
    SSH_OPTS=(-J "$BASTION")
    RSYNC_SSH=(-e "ssh -J ${BASTION}")
    note "routing via bastion ${BASTION}"
fi

echo "=========================================================================="
echo "push_to_cluster: ${REPO}/hpc  ->  ${REMOTE}:argus-validation/hpc"
[ "$DRY" -eq 1 ] && echo "                 DRY RUN"
echo "=========================================================================="

command -v rsync >/dev/null 2>&1 || die "rsync not found on this machine"

# ---------------------------------------------------- 1. stamp the commit ----
# git_commit is a MANDATORY run-card field and there will be no git repository
# on the cluster, so it is captured here, at the moment of transfer, and read
# back by argus_provenance in env.sh.
#
# A DIRTY TREE IS RECORDED, NOT REFUSED. Refusing would be wrong: it is entirely
# normal to push a work-in-progress hpc/ to test something. But a run card that
# says "commit abc1234" while the transferred files did not match abc1234 is a
# false provenance record, so the dirty flag travels with the hash.
GIT_COMMIT="$(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
if [ -n "$(git -C "$REPO" status --porcelain 2>/dev/null)" ]; then
    GIT_DIRTY=true
else
    GIT_DIRTY=false
fi
HPC_DIRTY=false
if [ -n "$(git -C "$REPO" status --porcelain -- hpc 2>/dev/null)" ]; then
    HPC_DIRTY=true
fi

note "commit ${GIT_COMMIT}, repo dirty ${GIT_DIRTY}, hpc/ dirty ${HPC_DIRTY}"
if [ "$HPC_DIRTY" = "true" ]; then
    cat >&2 <<'WARN'
push_to_cluster: WARNING: hpc/ has uncommitted changes. They WILL be transferred
push_to_cluster: and the stamped commit will NOT describe them. The stamp records
push_to_cluster: this as dirty so the run card cannot claim clean provenance, but
push_to_cluster: committing first is better than relying on that flag.
WARN
fi

if [ "$DRY" -eq 0 ]; then
    cat > "${HERE}/COMMIT" <<EOF
git_commit          ${GIT_COMMIT}
git_dirty           ${GIT_DIRTY}
git_dirty_hpc_only  ${HPC_DIRTY}
stamped_utc         $(date -u +%Y-%m-%dT%H:%M:%SZ)
stamped_from        $(hostname)
EOF
    note "stamped ${HERE}/COMMIT"
else
    note "dry run: would stamp ${HERE}/COMMIT with the values above"
fi

# ------------------------------------------------------ 2. reachability ------
# STAGED DIAGNOSIS, because "cannot reach X" is not a diagnosis.
#
# The first version of this check ran one ssh, threw its stderr away with
# >/dev/null 2>&1, and printed a generic three-cause message that LED WITH
# EduVPN. On its first real use the cause was a missing ~/.ssh/config alias, the
# host was reachable the whole time, and ssh's own message ("Could not resolve
# hostname delftblue") had been discarded. That is this repo's own no-silent-
# error-swallowing rule broken by the error handler itself: the tool that exists
# to explain the failure deleted the explanation.
#
# The three causes are DISTINGUISHABLE and are now distinguished:
#   name does not resolve          -> wrong alias or hostname, NOT the network
#   resolves but TCP 22 refused    -> VPN or network
#   TCP fine but authentication    -> key or password
note "diagnosing the path to ${REMOTE}"

# The bare hostname, with any user@ stripped, for the DNS and TCP probes.
PROBE_HOST="${REMOTE##*@}"

# An ssh config alias resolves inside ssh, not in DNS, so ask ssh what it will
# actually connect to before concluding anything about the name.
#
# THIS BLOCK IS DELIBERATELY VERBOSE, because the previous version of it was
# `RESOLVED_HOST="$(ssh -G "$REMOTE" 2>/dev/null | awk ...)"` and that had TWO
# defects, the second worse than the first:
#   1. 2>/dev/null discarded ssh's own explanation. Exactly the fault this
#      section was rewritten to remove, reintroduced one line below the comment
#      describing it.
#   2. Under `set -euo pipefail` a command substitution that exits non-zero
#      makes the ASSIGNMENT non-zero, which set -e turns into an immediate exit.
#      With stderr already discarded, the script died printing NOTHING AT ALL.
#      A wrong diagnosis is bad; a silent exit is worse, because there is
#      nothing to be wrong about.
# So: keep both streams and treat a non-zero status as a finding rather than as
# a reason to vanish.
#
# THE `|| SSH_G_RC=$?` FORM IS LOAD-BEARING, and `set +e` around it is NOT a
# substitute. An ERR trap fires INDEPENDENTLY OF errexit: `set +e` suppresses the
# exit but not the trap, so the first attempt at this fix still died, just with a
# better message. Bash exempts the left-hand side of a `||` list from the ERR
# trap, which is what actually makes a tolerated failure tolerated here.
SSH_G_RC=0
SSH_G_OUT="$(ssh -G "$REMOTE" 2>&1)" || SSH_G_RC=$?

if [ "$SSH_G_RC" -ne 0 ]; then
    cat >&2 <<EOF
push_to_cluster: FATAL: ssh cannot even parse its own configuration for '${REMOTE}'.
ssh said:
EOF
    printf '    %s\n' "$SSH_G_OUT" >&2
    cat >&2 <<EOF

  THIS IS YOUR ~/.ssh/config, not the network and not your key. ssh refuses to
  read the file, so nothing that uses ssh will work until it is valid.

  Look at it:
      cat ~/.ssh/config

  The block for this host must be exactly six lines, with FOUR SPACES of indent
  on the last five and no stray characters:

      Host delftblue
          HostName login.delftblue.tudelft.nl
          User <your-netid>
          IdentityFile ~/.ssh/id_ed25519
          ServerAliveInterval 60
          ServerAliveCountMax 10

  Editing it by hand is more reliable than pasting a printf:
      nano ~/.ssh/config
      chmod 600 ~/.ssh/config

  Or skip the alias entirely, which needs no config file at all:
      $0 --remote <your-netid>@login.delftblue.tudelft.nl
EOF
    exit 1
fi

RESOLVED_HOST="$(printf '%s\n' "$SSH_G_OUT" | awk '/^hostname /{print $2; exit}' || true)"
[ -n "$RESOLVED_HOST" ] || RESOLVED_HOST="$PROBE_HOST"
note "  ssh will connect to: ${RESOLVED_HOST}"

# A resolved host identical to the alias means ssh found no Host block for it,
# so it is about to try connecting to a literal machine called 'delftblue'.
if [ "$RESOLVED_HOST" = "$REMOTE" ] && [ "$REMOTE" = "delftblue" ]; then
    cat >&2 <<EOF
push_to_cluster: FATAL: ssh has no 'Host delftblue' block, so it would try to
push_to_cluster: connect to a machine literally named 'delftblue'.
  Add the block shown in README_HPC.md section 4, or bypass it:
      $0 --remote <your-netid>@login.delftblue.tudelft.nl
EOF
    exit 1
fi

if [ -z "$BASTION" ]; then
    if ! getent hosts "$RESOLVED_HOST" >/dev/null 2>&1; then
        cat >&2 <<EOF
push_to_cluster: FATAL: '${RESOLVED_HOST}' does not resolve.

  THIS IS NOT A NETWORK OR VPN PROBLEM. The name itself is unknown.

  If you passed no --remote, the default is the alias 'delftblue', which only
  works once you have the ~/.ssh/config block from README_HPC.md section 4.
  You do not have it yet. Two ways forward:

  (a) RECOMMENDED, do it once and never type a password again:
        ssh-keygen -t ed25519 -C "argus-delftblue"
        ssh-copy-id <netid>@login.delftblue.tudelft.nl
      then add to ~/.ssh/config:
        Host delftblue
            HostName login.delftblue.tudelft.nl
            User <netid>
            IdentityFile ~/.ssh/id_ed25519
            ServerAliveInterval 60
            ServerAliveCountMax 10
      then: chmod 600 ~/.ssh/config

  (b) skip the alias and name the host outright:
        $0 --remote <netid>@login.delftblue.tudelft.nl
EOF
        exit 1
    fi
    note "  DNS OK"

    if ! timeout 10 bash -c "cat < /dev/null > /dev/tcp/${RESOLVED_HOST}/22" 2>/dev/null; then
        cat >&2 <<EOF
push_to_cluster: FATAL: ${RESOLVED_HOST} resolves but port 22 is not reachable.

  THIS IS THE NETWORK. The DHPC docs: "A direct SSH to DelftBlue from outside
  of the university network is impossible!"

  (a) connect EduVPN (Institute Access), then re-run; or
  (b) route via the bastion:
        $0 --bastion <netid>@linux-bastion.tudelft.nl
      (students: student-linux.tudelft.nl)
EOF
        exit 1
    fi
    note "  TCP 22 OK"
fi

# Now authentication, and ssh's OWN message is shown rather than discarded.
# BatchMode=yes deliberately: it makes an ssh that WOULD have prompted fail
# instead, so a missing key is reported as a missing key rather than hanging on
# a prompt inside a script. If it fails this way the fallback below runs
# interactively, where a password prompt is legitimate.
SSH_ERR="$(ssh "${SSH_OPTS[@]}" -o ConnectTimeout=15 -o BatchMode=yes \
             "$REMOTE" 'echo ok' 2>&1 >/dev/null)" && AUTH_OK=1 || AUTH_OK=0

if [ "$AUTH_OK" -eq 0 ]; then
    echo "push_to_cluster: key-based login did not succeed. ssh said:" >&2
    printf '    %s\n' "$SSH_ERR" >&2
    case "$SSH_ERR" in
        *"Permission denied"*|*"publickey"*)
            cat >&2 <<EOF

  The host is reachable and it is refusing your credentials, not your route.
  It offers publickey and password. You have no usable key for it yet.

  STRONGLY RECOMMENDED, one interactive password entry, then never again:
      ssh-keygen -t ed25519 -C "argus-delftblue"
      ssh-copy-id <netid>@login.delftblue.tudelft.nl
  Re-run this script afterwards.

  To proceed on password auth instead, re-run with --interactive and ssh will
  prompt you. Type it at the prompt. NEVER paste a password into a chat, a
  script, a command line or a file: command lines are visible to every process
  on the machine and land in your shell history.
EOF
            ;;
    esac
    if [ "${ARGUS_INTERACTIVE:-0}" != "1" ]; then
        exit 1
    fi
    note "--interactive given: retrying, ssh will prompt on your terminal"
    ssh "${SSH_OPTS[@]}" -o ConnectTimeout=30 "$REMOTE" 'echo ok' >/dev/null \
        || die "interactive ssh also failed"
fi
note "ssh OK"

# --------------------------------------------------- 2b. sync-state report ---
# WHY THIS EXISTS. On 2026-07-29 a push ran at 17:47:28 and `sbatch --test-only`
# ran 29 seconds later against a header that had been fixed AFTER that push. The
# output was byte-identical to the pre-fix run, including the warning the fix had
# removed, and NOTHING ON THE CLUSTER SIDE COULD HAVE REVEALED IT. The operator
# reasonably assumed the file was current.
#
# The cluster carries hpc/COMMIT, so the drift is computable rather than
# guessable. Reported on every run, and `--check` reports it and stops.
# `|| RC=$?` rather than `|| true`, and not because of style. A missing COMMIT
# on the cluster is an EXPECTED state (nothing pushed yet) and is reported as
# ABSENT below; a failed ssh is NOT, and `|| true` would render the two
# indistinguishable. The || form is also what exempts this from the ERR trap,
# which `set +e` would not do.
RC_RC=0
REMOTE_COMMIT="$(ssh "${SSH_OPTS[@]}" "$REMOTE" \
    'sed -n "s/^git_commit *//p" $HOME/argus-validation/hpc/COMMIT 2>/dev/null' 2>&1)" || RC_RC=$?
if [ "$RC_RC" -ne 0 ] && [ "$RC_RC" -ne 1 ]; then
    echo "push_to_cluster: WARNING: could not read the cluster's COMMIT stamp (ssh exit ${RC_RC})." >&2
    printf '    %s\n' "$REMOTE_COMMIT" >&2
    echo "push_to_cluster: sync state is UNKNOWN, not 'in sync'. Treat the cluster as stale." >&2
    REMOTE_COMMIT=""
fi
REMOTE_COMMIT="$(printf '%s' "$REMOTE_COMMIT" | tr -d '[:space:]')"

echo
echo "--- sync state ---"
if [ -z "$REMOTE_COMMIT" ]; then
    echo "    cluster:  NO hpc/COMMIT -- nothing pushed yet, or pushed before stamping existed"
    SYNC="ABSENT"
elif [ "$REMOTE_COMMIT" = "$GIT_COMMIT" ]; then
    echo "    cluster:  ${REMOTE_COMMIT}"
    echo "    laptop:   ${GIT_COMMIT}"
    echo "    IN SYNC"
    SYNC="OK"
else
    echo "    cluster:  ${REMOTE_COMMIT}"
    echo "    laptop:   ${GIT_COMMIT}"
    echo "    STALE. Commits on the laptop that the cluster does not have:"
    git -C "$REPO" log --oneline "${REMOTE_COMMIT}..HEAD" -- hpc 2>/dev/null \
        | sed 's/^/        /' || echo "        (cannot list; the cluster commit is not in this repo)"
    echo "    Files that differ:"
    git -C "$REPO" diff --stat "${REMOTE_COMMIT}..HEAD" -- hpc 2>/dev/null \
        | sed 's/^/        /' || true
    SYNC="STALE"
fi

if [ "$CHECK_ONLY" -eq 1 ]; then
    echo
    case "$SYNC" in
        OK)     note "--check: cluster is current. Nothing to do."; exit 0 ;;
        STALE)  note "--check: cluster is STALE. Re-run without --check to update it."; exit 1 ;;
        ABSENT) note "--check: nothing on the cluster. Re-run without --check to push."; exit 1 ;;
    esac
fi

# ---------------------------------------------------------- 3. transfer ------
# rsync -av is the DHPC-documented form. -z added because this is all text and
# compresses well; --delete so a file removed on the laptop does not linger on
# the cluster and get sourced by a job.
RSYNC_OPTS=(-avz --delete --progress)
[ "$DRY" -eq 1 ] && RSYNC_OPTS+=(--dry-run)

note "creating the remote layout"
if [ "$DRY" -eq 0 ]; then
    ssh "${SSH_OPTS[@]}" "$REMOTE" \
        'mkdir -p $HOME/argus-validation && mkdir -p /scratch/$USER/argus/bundles /scratch/$USER/argus/cases'
else
    note "dry run: would mkdir \$HOME/argus-validation and /scratch/\$USER/argus/{bundles,cases}"
fi

note "rsync hpc/"
rsync "${RSYNC_OPTS[@]}" "${RSYNC_SSH[@]}" \
    "${REPO}/hpc/" "${REMOTE}:argus-validation/hpc/"

if [ "$DRY" -eq 1 ]; then
    echo
    note "DRY RUN COMPLETE. Nothing transferred, nothing stamped."
    exit 0
fi

# --------------------------------------------------------- 4. verify ---------
# WHAT MAKES THIS FAIL: hpc/env.sh missing or unreadable on the cluster, the
# COMMIT stamp not having arrived, or the scripts landing without +x. Each has
# happened to somebody; none is visible from the rsync summary alone.
echo
echo "--- verifying what landed ---"
ssh "${SSH_OPTS[@]}" "$REMOTE" 'bash -s' <<'REMOTE_CHECK'
set -euo pipefail
H="$HOME/argus-validation/hpc"
fail=0
for f in env.sh slurm/smoke.slurm slurm/solve.slurm slurm/mesh.slurm \
         package_case.sh unpack_case.sh pull_results.sh selftest.sh COMMIT; do
    if [ -e "$H/$f" ]; then
        printf '    ok      %s\n' "$f"
    else
        printf '    MISSING %s\n' "$f"; fail=1
    fi
done
echo "    --- stamped provenance ---"
sed 's/^/        /' "$H/COMMIT" 2>/dev/null || { echo "        NO COMMIT FILE"; fail=1; }
echo "    --- storage ---"
printf '    /home used: %s\n' "$(du -sh "$HOME" 2>/dev/null | cut -f1)"
printf '    scratch dirs: %s\n' "$(ls -d /scratch/$USER/argus/* 2>/dev/null | tr '\n' ' ')"
if ! bash -n "$H/slurm/smoke.slurm"; then
    echo "    smoke.slurm FAILED bash -n on the cluster"; fail=1
fi
[ "$fail" -eq 0 ] || { echo "    VERIFICATION FAILED"; exit 1; }
echo "    all present, smoke.slurm parses"
REMOTE_CHECK

echo
echo "=========================================================================="
echo "push_to_cluster: DONE"
echo "  commit stamped   ${GIT_COMMIT} (dirty: ${GIT_DIRTY})"
echo
echo "  THIS SCRIPT TRANSFERS SCRIPTS ONLY, never case bundles: those are 61 MB"
echo "  each and go to /scratch, not /home. Build with package_case.sh, then"
echo "    rsync -avz --progress ~/argus-bundles/<name>.tar.gz \\"
echo "          ${REMOTE}:${ARGUS_BUNDLE_STAGE}/"
echo
echo "  NEXT STEP depends on where you are in hpc/jobs/job_index.md."
if ssh "${SSH_OPTS[@]}" "$REMOTE" \
     'ls /scratch/$USER/argus/smoke/*/cavity/log/icoFoam.log' >/dev/null 2>&1; then
    echo "  Job 0 has run on this cluster (a smoke work directory exists)."
    echo "  If it passed, the next thing is job 1, the migration gate:"
    echo "    ssh ${REMOTE}"
    echo "    bash \$HOME/argus-validation/hpc/unpack_case.sh \\"
    echo "         ${ARGUS_BUNDLE_STAGE}/A2_base_a.tar.gz"
    echo "  It must print MESH ASSERTION PASSED before anything is submitted."
else
    echo "  No job-0 smoke directory found on the cluster. Start there:"
    echo "    ssh ${REMOTE}"
    echo "    sbatch \$HOME/argus-validation/hpc/slurm/smoke.slurm"
fi
echo "=========================================================================="
