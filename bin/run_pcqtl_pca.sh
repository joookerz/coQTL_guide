#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <config.R>"
  exit 1
fi

CONFIG_PATH="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
PCQTL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../scripts/pcqtl" && pwd)"

export COQTL_CONFIG="${CONFIG_PATH}"

Rscript "${PCQTL_DIR}/step2_cluster_pca.R"
Rscript "${PCQTL_DIR}/step3_prepare_saige_inputs.R"
