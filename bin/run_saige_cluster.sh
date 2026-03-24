#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <config.R> <cluster_id>"
  exit 1
fi

CONFIG_PATH="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
CLUSTER_ID="$2"

cfg() {
  Rscript -e "cfg_env <- new.env(); sys.source('${CONFIG_PATH}', envir = cfg_env); cat(${1})"
}

PCQTL_DIR="$(cfg "file.path(cfg_env\$WORK_ROOT, 'pcQTL')")"
PC_MAP="${PCQTL_DIR}/step3_saige/cluster_pc_map.tsv"
if [[ ! -f "${PC_MAP}" ]]; then
  echo "Missing ${PC_MAP}. Run run_pcqtl_pca.sh first."
  exit 1
fi

CHR="$(awk -F'\t' -v cid="${CLUSTER_ID}" 'NR>1 && $1==cid {print $4; exit}' "${PC_MAP}")"
if [[ -z "${CHR}" ]]; then
  echo "Cluster not found in ${PC_MAP}: ${CLUSTER_ID}"
  exit 1
fi

PHENO_FILE="${PCQTL_DIR}/step2_pca/${CLUSTER_ID}/pheno_with_pcs.tsv"
REGION_FILE="${PCQTL_DIR}/step3_saige/regions/${CLUSTER_ID}_region.txt"
SAIGE_DIR="${PCQTL_DIR}/step3_saige"
mkdir -p "${SAIGE_DIR}/step1/${CLUSTER_ID}" "${SAIGE_DIR}/step2/${CLUSTER_ID}" "${SAIGE_DIR}/step3/${CLUSTER_ID}"

SAIGE_PREFIX="$(cfg "cfg_env\$SAIGE_CMD_PREFIX")"
STEP1_CMD="$(cfg "cfg_env\$SAIGE_STEP1_CMD")"
STEP2_CMD="$(cfg "cfg_env\$SAIGE_STEP2_CMD")"
STEP3_CMD="$(cfg "cfg_env\$SAIGE_STEP3_CMD")"
VR_PLINK="$(cfg "cfg_env\$VR_PLINK_PREFIX")"
COVARS="$(cfg "paste(cfg_env\$SAIGE_COVAR_COLS, collapse = ',')")"
BED_FILE="$(cfg "gsub('\\{chr\\}', '${CHR}', cfg_env\$GENO_BED_PATTERN)")"
BIM_FILE="$(cfg "gsub('\\{chr\\}', '${CHR}', cfg_env\$GENO_BIM_PATTERN)")"
FAM_FILE="$(cfg "gsub('\\{chr\\}', '${CHR}', cfg_env\$GENO_FAM_PATTERN)")"
MIN_MAF="$(cfg "cfg_env\$MIN_MAF")"
MARKERS_PER_CHUNK="$(cfg "cfg_env\$MARKERS_PER_CHUNK")"

mapfile -t PCS < <(awk -F'\t' -v cid="${CLUSTER_ID}" 'NR>1 && $1==cid {print $2}' "${PC_MAP}")

run_cmd() {
  if [[ -n "${SAIGE_PREFIX}" ]]; then
    local quoted_args=""
    printf -v quoted_args ' %q' "$@"
    bash -lc "${SAIGE_PREFIX}${quoted_args}"
  else
    "$@"
  fi
}

for PC in "${PCS[@]}"; do
  STEP1_PFX="${SAIGE_DIR}/step1/${CLUSTER_ID}/${PC}"
  STEP2_OUT="${SAIGE_DIR}/step2/${CLUSTER_ID}/${PC}"
  STEP3_OUT="${SAIGE_DIR}/step3/${CLUSTER_ID}/${PC}_genePval"

  run_cmd "${STEP1_CMD}" \
    --useSparseGRMtoFitNULL=FALSE \
    --useGRMtoFitNULL=FALSE \
    --phenoFile="${PHENO_FILE}" \
    --phenoCol="${PC}" \
    --covarColList="${COVARS}" \
    --sampleCovarColList="${COVARS}" \
    --sampleIDColinphenoFile=individual \
    --traitType=quantitative \
    --outputPrefix="${STEP1_PFX}" \
    --invNormalize=TRUE \
    --skipVarianceRatioEstimation=FALSE \
    --isRemoveZerosinPheno=FALSE \
    --isCovariateOffset=FALSE \
    --isCovariateTransform=TRUE \
    --skipModelFitting=FALSE \
    --tol=0.00001 \
    --plinkFile="${VR_PLINK}" \
    --IsOverwriteVarianceRatioFile=TRUE

  run_cmd "${STEP2_CMD}" \
    --bedFile="${BED_FILE}" \
    --bimFile="${BIM_FILE}" \
    --famFile="${FAM_FILE}" \
    --SAIGEOutputFile="${STEP2_OUT}" \
    --chrom="${CHR}" \
    --minMAF="${MIN_MAF}" \
    --LOCO=FALSE \
    --GMMATmodelFile="${STEP1_PFX}.rda" \
    --varianceRatioFile="${STEP1_PFX}.varianceRatio.txt" \
    --rangestoIncludeFile="${REGION_FILE}" \
    --markers_per_chunk="${MARKERS_PER_CHUNK}"

  run_cmd "${STEP3_CMD}" \
    --assocFile="${STEP2_OUT}" \
    --geneName="${PC}" \
    --genePval_outputFile="${STEP3_OUT}"
done
