# Example config for one cell type.
# Copy and modify this file for your own dataset.

CELL_TYPE <- "B_in"

# Input files
COUNT_FILE <- "/data/project/input/B_in_readcounts.tsv.gz"
GENE_INFO_FILE <- "/data/project/reference/gene_info_with_location.tsv"

# Output root for this cell type
WORK_ROOT <- "/data/project/coqtl_runs/B_in"

# Optional column names in the count matrix
CELL_ID_COL_CANDIDATES <- c("CellID", "barcode")
INDIVIDUAL_COL_CANDIDATES <- c("individual", "IndividualID")
CELLTYPE_COL_CANDIDATES <- c("CellType")

# Method 2
NONZERO_CUTOFF <- 0.01
P_THRESHOLD <- 0.05
CHUNK_SIZE_GENES <- 50
SECS_PER_PAIR <- 0.5
MAX_WINDOW_SIZE <- 50
MIN_WINDOW_SIZE <- 2
CLUSTER_THRESHOLD <- 0.70

# PCA and phenotype generation
COVAR_COLS <- c("age", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pf1", "pf2")
VARIANCE_THRESHOLD <- 0.95
CIS_WINDOW <- 500000

# SAIGE-QTL
SAIGE_COVAR_COLS <- COVAR_COLS

# If SAIGE-QTL is installed directly and scripts are in PATH:
SAIGE_CMD_PREFIX <- ""
SAIGE_STEP1_CMD <- "step1_fitNULLGLMM_qtl.R"
SAIGE_STEP2_CMD <- "step2_tests_qtl.R"
SAIGE_STEP3_CMD <- "step3_gene_pvalue_qtl.R"

# Example if using Singularity/Apptainer instead:
# SAIGE_CMD_PREFIX <- "singularity exec --bind /data:/data --cleanenv /path/to/saigeqtl.sif"

VR_PLINK_PREFIX <- "/data/project/genotype/pruned/variance_ratio_markers"
GENO_BED_PATTERN <- "/data/project/genotype/raw/full_genome_chr{chr}.bed"
GENO_BIM_PATTERN <- "/data/project/genotype/raw/full_genome_chr{chr}.bim"
GENO_FAM_PATTERN <- "/data/project/genotype/raw/full_genome_chr{chr}.fam"
MIN_MAF <- 0.05
MARKERS_PER_CHUNK <- 10000
