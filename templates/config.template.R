# Copy this file to <your_celltype_workdir>/config.R and edit the values.

CELL_TYPE <- "example_celltype"

# Input files
COUNT_FILE <- "/abs/path/to/example_celltype_readcounts.tsv.gz"
GENE_INFO_FILE <- "/abs/path/to/gene_info_with_location.tsv"

# Output root for this cell type
WORK_ROOT <- "/abs/path/to/coQTL_runs/example_celltype"

# Optional: column names in the count matrix
CELL_ID_COL_CANDIDATES <- c("CellID", "barcode")
INDIVIDUAL_COL_CANDIDATES <- c("individual", "IndividualID")
CELLTYPE_COL_CANDIDATES <- c("CellType")

# Method 2: cluster identification
NONZERO_CUTOFF <- 0.01
P_THRESHOLD <- 0.05
CHUNK_SIZE_GENES <- 50
SECS_PER_PAIR <- 0.5
MAX_WINDOW_SIZE <- 50
MIN_WINDOW_SIZE <- 2
CLUSTER_THRESHOLD <- 0.70

# pcQTL: PCA on each cluster
COVAR_COLS <- c("age", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pf1", "pf2")
VARIANCE_THRESHOLD <- 0.95
CIS_WINDOW <- 500000

# SAIGE-QTL
SAIGE_COVAR_COLS <- COVAR_COLS
SAIGE_STEP1_CMD <- "step1_fitNULLGLMM_qtl.R"
SAIGE_STEP2_CMD <- "step2_tests_qtl.R"
SAIGE_STEP3_CMD <- "step3_gene_pvalue_qtl.R"
SAIGE_CMD_PREFIX <- ""
VR_PLINK_PREFIX <- "/abs/path/to/variance_ratio_plink_prefix"
GENO_BED_PATTERN <- "/abs/path/to/full_genome_chr{chr}.bed"
GENO_BIM_PATTERN <- "/abs/path/to/full_genome_chr{chr}.bim"
GENO_FAM_PATTERN <- "/abs/path/to/full_genome_chr{chr}.fam"
MIN_MAF <- 0.05
MARKERS_PER_CHUNK <- 10000
