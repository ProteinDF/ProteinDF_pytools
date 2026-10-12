#!/usr/bin/env bash
# Reviewer-side verification of a working branch (run by Claude, from any worktree).
# Runs check.sh in the branch's worktree and saves the full output for the review record.
#
# Usage: devtool/review.sh [--verbose] <branch>
#   --verbose  print the whole output (default: header, then the summary only)
#
# The full output is always saved under <git-common-dir>/review-logs/.
set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
verbose=0
while (($#)); do
    case "$1" in
        --verbose) verbose=1; shift ;;
        -*) sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
        *) break ;;
    esac
done
(($# == 1)) || { sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
branch=$1

wt=$("$here/flow.sh" path "$branch") || exit 1
logdir="$(git rev-parse --path-format=absolute --git-common-dir)/review-logs"
mkdir -p "$logdir"
log="$logdir/$(date +%Y%m%d-%H%M%S)_${branch//\//-}.log"

run() {
    echo "branch:   $branch"
    echo "worktree: $wt"
    echo
    echo "--- commits (develop..$branch)"
    git -C "$wt" log --oneline "develop..$branch"
    echo
    echo "--- diffstat"
    git -C "$wt" diff --stat "develop...$branch"
    echo
    echo "--- worktree status (uncommitted / untracked)"
    git -C "$wt" status --short
    echo
    echo "--- check.sh"
    # use the reviewer's (develop) check.sh so an edited copy on the branch cannot weaken it
    (cd "$wt" && "$here/check.sh")
    check_status=$?
    echo
    echo "=== review: check.sh exit $check_status"
    ((check_status == 0))
}

if ((verbose)); then
    run 2>&1 | tee "$log"
    status=${PIPESTATUS[0]}
else
    run >"$log" 2>&1
    status=$?
    sed -n '1,/^--- worktree status/p' "$log" | sed '$d'
    sed -n '/^=== summary/,$p' "$log"
    sed -n '/^--- worktree status/,/^--- check.sh/p' "$log" | sed '1d;$d' | grep . | sed 's/^/uncommitted: /'
fi
echo "log: $log"
exit "$status"
