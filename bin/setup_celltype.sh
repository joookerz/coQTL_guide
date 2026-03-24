#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <celltype_name> <work_root>"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GUIDE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CELLTYPE="$1"
TARGET_DIR="$2"

mkdir -p "$(dirname "${TARGET_DIR}")"
WORK_ROOT="$(cd "$(dirname "${TARGET_DIR}")" && pwd)/$(basename "${TARGET_DIR}")"

mkdir -p "${WORK_ROOT}/cluster_identification/method2_sc_hurdle"
mkdir -p "${WORK_ROOT}/pcQTL"

cp "${GUIDE_ROOT}/templates/config.template.R" "${WORK_ROOT}/config.R"
cp "${GUIDE_ROOT}/scripts/method2/"*.R "${WORK_ROOT}/cluster_identification/method2_sc_hurdle/"
cp "${GUIDE_ROOT}/scripts/pcqtl/"*.R "${WORK_ROOT}/pcQTL/"

sed -i "s|CELL_TYPE <- \"example_celltype\"|CELL_TYPE <- \"${CELLTYPE}\"|" "${WORK_ROOT}/config.R"
sed -i "s|WORK_ROOT <- \"/abs/path/to/coQTL_runs/example_celltype\"|WORK_ROOT <- \"${WORK_ROOT}\"|" "${WORK_ROOT}/config.R"

chmod +x "${GUIDE_ROOT}/bin/"*.sh

echo "Created work directory: ${WORK_ROOT}"
echo "Edit ${WORK_ROOT}/config.R before running the pipeline."
