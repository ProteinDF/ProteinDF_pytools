#!/usr/bin/env bash
# Create the virtualenv that devtool/check.sh uses: <git-common-dir>/check-venv.
#
# Usage: devtool/setup-check-venv.sh [--recreate]
#   PDFTOOLS_BRIDGE_DIR  source directory of ProteinDF_bridge (default: pip package "proteindf_bridge");
#                        a wheel file also works. Needs proteindf_bridge >= 2026.10.6.
#   PDFTOOLS_BASE_PYTHON Python 3.12 interpreter to use (default: python3.12 in PATH)
#
# Notes:
#   - Do not reuse a virtualenv that has a sitecustomize.py: it can patch proteindf_tools at
#     startup and hide bugs.
#   - proteindf_bridge >= 2026.10.6 works with numpy 2 (older releases do not).
set -euo pipefail

recreate=0
[[ "${1:-}" == --recreate ]] && recreate=1
common_dir=$(git rev-parse --path-format=absolute --git-common-dir)
venv=$common_dir/check-venv
py=${PDFTOOLS_BASE_PYTHON:-python3.12}
"$py" -c "import venv" 2>/dev/null || { echo "error: $py is not a working Python; set PDFTOOLS_BASE_PYTHON" >&2; exit 1; }

if [[ -e "$venv" ]]; then
    ((recreate)) || { echo "error: $venv exists (use --recreate)" >&2; exit 1; }
    rm -rf "$venv"
fi
"$py" -m venv "$venv"
"$venv/bin/pip" install -q numpy h5py matplotlib "${PDFTOOLS_BRIDGE_DIR:-proteindf_bridge}"
"$venv/bin/python" -c "import proteindf_bridge, h5py, numpy, matplotlib"
echo "created $venv"
find "$venv" -name sitecustomize.py | grep . && echo "warning: sitecustomize.py found" >&2 || true
