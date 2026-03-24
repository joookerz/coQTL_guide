#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <config.R>"
  exit 1
fi

CONFIG_PATH="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
METHOD2_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../scripts/method2" && pwd)"

export COQTL_CONFIG="${CONFIG_PATH}"

Rscript "${METHOD2_DIR}/step1_filter_sparse_genes.R"
Rscript "${METHOD2_DIR}/calculate_chunks.R"

CHUNK_PLAN="${METHOD2_DIR}/chunk_plan.tsv"
while IFS=$'\t' read -r chromosome chunk_id total_chunks gene_start gene_end n_genes_chunk; do
  [[ "${chromosome}" == "chromosome" ]] && continue
  Rscript "${METHOD2_DIR}/step2_calculate_sc_associations_chunked.R" \
    --chr "${chromosome}" \
    --gene_start "${gene_start}" \
    --gene_end "${gene_end}" \
    --chunk_id "${chunk_id}"
done < "${CHUNK_PLAN}"

Rscript "${METHOD2_DIR}/step2b_merge_chunks.R" --chr 0

for chr in $(seq 1 22); do
  Rscript "${METHOD2_DIR}/step3_identify_clusters_single_chr.R" "${chr}"
done

Rscript "${METHOD2_DIR}/step4_merge_clusters.R"
