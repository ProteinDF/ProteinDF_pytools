#!/usr/bin/env bash
# Create the virtualenv that devtool/check.sh uses: <git-common-dir>/check-venv.
#
# Usage: devtool/setup-check-venv.sh [--recreate]
#   PDFTOOLS_BRIDGE_DIR  source directory of ProteinDF_bridge (default: pip package "proteindf_bridge")
#   PDFTOOLS_BASE_PYTHON Python 3.12 interpreter to use (default: python3.12 in PATH)
#
# Notes:
#   - Do not reuse a virtualenv that has a sitecustomize.py: it can patch proteindf_tools at
#     startup and hide bugs.
#   - numpy must be < 2: proteindf_bridge uses ndarray.tostring / numpy.fromstring.
set -euo pipefail

recreate=0
[[ "${1:-}" == --recreate ]] && recreate=1
common_dir=$(git rev-parse --path-format=absolute --git-common-dir)
venv=$common_dir/check-venv
py=${PDFTOOLS_BASE_PYTHON:-python3.12}
command -v "$py" >/dev/null || { echo "error: $py not found; set PDFTOOLS_BASE_PYTHON" >&2; exit 1; }

if [[ -e "$venv" ]]; then
    ((recreate)) || { echo "error: $venv exists (use --recreate)" >&2; exit 1; }
    rm -rf "$venv"
fi
"$py" -m venv "$venv"
"$venv/bin/pip" install -q "numpy<2" h5py matplotlib "${PDFTOOLS_BRIDGE_DIR:-proteindf_bridge}"
"$venv/bin/python" -c "import proteindf_bridge, h5py, numpy, matplotlib"
echo "created $venv"
find "$venv" -name sitecustomize.py | grep . && echo "warning: sitecustomize.py found" >&2 || true
