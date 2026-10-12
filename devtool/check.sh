#!/usr/bin/env bash
# Definition-of-done checks for a working branch. Run it in the branch's worktree.
# Both the implementing agent and the reviewer run this and paste its output verbatim.
#
# Usage: devtool/check.sh [--base <ref>] [--no-tests] [--verbose]
#   --base      ref to diff against (default: develop)
#   --no-tests  only run the diff checks
#   --verbose   print the whole test output (default: the last lines; the full output
#               is in the temporary log whose path is printed)
# Tests need PDF_HOME (ProteinDF source tree) and a Python with proteindf_bridge.
# Python: PDFTOOLS_PYTHON, else <git-common-dir>/check-venv/bin/python, else python3
# (it must be able to import proteindf_bridge). Tests are the doctests and unit tests of run_test.sh.
set -uo pipefail

base=develop
do_tests=1
verbose=0
while (($#)); do
    case "$1" in
        --base) base=$2; shift 2 ;;
        --no-tests) do_tests=0; shift ;;
        --verbose) verbose=1; shift ;;
        *) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
    esac
done
common_dir=$(git rev-parse --path-format=absolute --git-common-dir) || exit 1
if [[ -n "${PDFTOOLS_PYTHON:-}" ]]; then
    py=$PDFTOOLS_PYTHON
elif [[ -x "$common_dir/check-venv/bin/python" ]]; then
    py=$common_dir/check-venv/bin/python
else
    py=python3
fi

top=$(git rev-parse --show-toplevel) || exit 1
cd "$top"
merge_base=$(git merge-base "$base" HEAD) || exit 1

declare -a results=()
failed=0
record() { # <name> <PASS|FAIL|SKIP> [note]
    results+=("$(printf '%-14s %-4s %s' "$1" "$2" "${3:-}")")
    [[ "$2" == FAIL ]] && failed=1
    return 0
}

echo "branch: $(git branch --show-current)  HEAD: $(git rev-parse --short HEAD)  base: $base ($(git rev-parse --short "$merge_base"))"
echo "commits:"
git log --oneline "$merge_base..HEAD" | sed 's/^/  /'
mapfile -t changed < <(git diff --name-only --diff-filter=d "$merge_base" HEAD)
mapfile -t changed_py < <(printf '%s\n' "${changed[@]}" | grep -E '\.py$')
echo "changed files: ${#changed[@]} (Python: ${#changed_py[@]})"
[[ -n "$(git status --porcelain --untracked-files=no)" ]] \
    && echo "WARNING: uncommitted changes exist; they are tested but not part of the diff checks"

# 1. whitespace errors in the branch diff
echo; echo "--- whitespace (git diff --check)"
if git diff --check "$merge_base" HEAD; then
    record whitespace PASS
else
    record whitespace FAIL
fi

# 2. syntax of the changed Python files
echo; echo "--- syntax (py_compile)"
if ((${#changed_py[@]} == 0)); then
    record syntax SKIP "no Python changes"
elif "$py" -I -c 'import sys, py_compile
bad = 0
for f in sys.argv[1:]:
    try:
        py_compile.compile(f, doraise=True, cfile="/dev/null")
    except Exception as e:
        print(e); bad = 1
sys.exit(bad)' "${changed_py[@]}"; then
    record syntax PASS
else
    record syntax FAIL
fi

# 3. black on changed lines only (the existing code is not fully formatted)
echo; echo "--- black (changed files, diff only)"
if ((${#changed_py[@]} == 0)); then
    record black SKIP "no Python changes"
elif ! "$py" -m black --version >/dev/null 2>&1; then
    record black SKIP "black not installed"
else
    out=$("$py" -m black --check --diff -q -l 100 "${changed_py[@]}" 2>&1)
    if [[ -z "$out" ]]; then
        record black PASS
    else
        echo "$out" | head -60
        record black FAIL "formatting differs (existing files may already differ; check only the changed lines)"
    fi
fi

# 4. tests: the doctests and unit tests listed in run_test.sh, run against this worktree's
#    sources (PYTHONPATH). Needs PDF_HOME (the ProteinDF source tree; data/basis2 is read).
#    doctest.testmod() exits 0 even on failure, so doctests are judged by their output.
if ((do_tests)); then
    echo; echo "--- tests (doctest + unittest from run_test.sh)"
    log=$(mktemp -t pdftools-check-tests.XXXXXX)
    if [[ -z "${PDF_HOME:-}" ]]; then
        record tests SKIP "PDF_HOME is not set"
    elif ! PYTHONPATH=$top "$py" -I -c 'import sys; sys.path.insert(0, sys.argv[1]); import proteindf_bridge' "$top" 2>/dev/null; then
        record tests SKIP "$py cannot import proteindf_bridge"
    else
        mapfile -t doctests < <(sed -n '/^DOCTESTS=/,/^  "$/p' run_test.sh | grep -oE '^\s+[a-z0-9_]+' | tr -d ' ')
        mapfile -t unittests < <(sed -n '/^UNITTESTS=/,/^  "$/p' run_test.sh | grep -oE '^\s+[a-z0-9_]+' | tr -d ' ')
        n_fail=0; n_run=0
        for m in "${doctests[@]}"; do
            n_run=$((n_run + 1))
            echo ">>>> doctest: $m" >>"$log"
            out=$(PYTHONPATH=$top "$py" -m "proteindf_tools.$m" 2>&1); st=$?
            printf '%s\n' "$out" >>"$log"
            if ((st != 0)) || grep -qE 'Failed example|Traceback' <<<"$out"; then
                echo "FAILED: doctest $m"; n_fail=$((n_fail + 1))
            fi
        done
        for m in "${unittests[@]}"; do
            n_run=$((n_run + 1))
            echo ">>>> unittest: $m" >>"$log"
            if ! PYTHONPATH=$top "$py" -m unittest "tests.test_$m" >>"$log" 2>&1; then
                echo "FAILED: unittest $m"; n_fail=$((n_fail + 1))
            fi
        done
        ((verbose)) && cat "$log"
        if ((n_fail == 0)); then
            record tests PASS "$n_run module(s); log: $log"
        else
            record tests FAIL "$n_fail of $n_run module(s) failed; log: $log"
        fi
    fi
fi

echo; echo "=== summary"
printf '%s\n' "${results[@]}"
exit $failed
