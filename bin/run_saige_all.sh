#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <config.R>"
  exit 1
fi

CONFIG_PATH="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
PCQTL_DIR="$(Rscript -e "cfg_env <- new.env(); sys.source('${CONFIG_PATH}', envir = cfg_env); cat(file.path(cfg_env\$WORK_ROOT, 'pcQTL'))")"
PC_MAP="${PCQTL_DIR}/step3_saige/cluster_pc_map.tsv"

if [[ ! -f "${PC_MAP}" ]]; then
  echo "Missing ${PC_MAP}. Run run_pcqtl_pca.sh first."
  exit 1
fi

mapfile -t CLUSTERS < <(awk -F'\t' 'NR>1 {print $1}' "${PC_MAP}" | awk '!seen[$0]++')
for CLUSTER_ID in "${CLUSTERS[@]}"; do
  "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_saige_cluster.sh" "${CONFIG_PATH}" "${CLUSTER_ID}"
done
