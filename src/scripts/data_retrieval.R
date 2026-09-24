suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(edgeR)
  library(sva)
  library(ggplot2)
  library(EnsDb.Hsapiens.v86)
})

# ==============================================================================
# User Parameters
# ==============================================================================
TCGA_PROJECT_NAME <- "COAD"            # TCGA Colon Adenocarcinoma
GTEX_TISSUE_NAME  <- "COLON"           # GTEx Colon
SRA_PROJECT_NAME  <- "SRP029880"
MIN_RIN_SCORE     <- 6.0               # GTEx RIN threshold
EXCLUDE_HARDY_4   <- TRUE              # GTEx Exclude Hardy 4 donors

data_dir <- "data/COAD-COLON-SRP"

cat("=== Starting recount3 TCGA, GTEx & SRP029880 Preprocessing Pipeline ===\n\n")

# ==============================================================================
# [Step 1] Fetch & Extract TCGA Counts
# ==============================================================================
cat("[Step 1] Fetching TCGA Data...\n")
human_projects <- available_projects(organism = "human")

tcga_info <- subset(human_projects, toupper(project) == toupper(TCGA_PROJECT_NAME) & file_source == "tcga")
if (nrow(tcga_info) == 0) stop("TCGA project not found!")

rse_tcga <- create_rse(tcga_info)
counts_tcga_full <- compute_read_counts(rse_tcga)

# Filter TCGA Primary Tumors using tcga.cgc_sample_sample_type
tcga_cd <- colData(rse_tcga)
tcga_keep <- tcga_cd$tcga.cgc_sample_sample_type %in% c("Primary Tumor", "Solid Tissue Normal")

counts_tcga <- counts_tcga_full[, tcga_keep]
cat(sprintf("  -> Retained %d TCGA samples (Total reads: %.2e)\n\n", 
            ncol(counts_tcga), sum(counts_tcga)))

# ==============================================================================
# [Step 2] Fetch & Extract GTEx Counts
# ==============================================================================
cat("[Step 2] Fetching GTEx Data...\n")
gtex_info <- subset(human_projects, toupper(project) == toupper(GTEX_TISSUE_NAME) & file_source == "gtex")
if (nrow(gtex_info) == 0) stop("GTEx tissue project not found!")

rse_gtex <- create_rse(gtex_info)
counts_gtex_full <- compute_read_counts(rse_gtex)

# Filter GTEx Quality Controls
gtex_cd    <- colData(rse_gtex)
gtex_rin   <- as.numeric(gtex_cd$gtex.smrin)
gtex_hardy <- as.numeric(gtex_cd$gtex.dthhrdy)

gtex_keep  <- (!is.na(gtex_rin) & gtex_rin >= MIN_RIN_SCORE) & 
  (is.na(gtex_hardy) | gtex_hardy != 4)

counts_gtex <- counts_gtex_full[, gtex_keep]
cat(sprintf("  -> Retained %d GTEx Control samples (Total reads: %.2e)\n\n", 
            ncol(counts_gtex), sum(counts_gtex)))

# ==============================================================================
# [Step 2b] Fetch & Extract SRP029880 (Balanced Cohort)
# ==============================================================================
cat("[Step 2b] Fetching SRP029880 Data...\n")
srp_info <- subset(human_projects, project == SRA_PROJECT_NAME & project_home == "data_sources/sra")
if (nrow(srp_info) == 0) stop("SRP029880 project not found!")

rse_srp <- create_rse(srp_info)
counts_srp_full <- compute_read_counts(rse_srp)

# Parse metadata to separate Primary Tumors from Normal Epithelium, excluding Liver Metastases
srp_attrs <- colData(rse_srp)$sra.sample_attributes
srp_is_normal <- grepl("normal", srp_attrs, ignore.case = TRUE)
srp_is_tumor  <- grepl("primary", srp_attrs, ignore.case = TRUE) & !grepl("metastasis", srp_attrs, ignore.case = TRUE)

srp_keep <- srp_is_normal | srp_is_tumor
counts_srp <- counts_srp_full[, srp_keep]

srp_conditions <- ifelse(srp_is_tumor[srp_keep], "Tumor", "Control")

cat(sprintf("  -> Retained %d SRP029880 samples (%d Tumor, %d Control)\n\n", 
            ncol(counts_srp), sum(srp_conditions == "Tumor"), sum(srp_conditions == "Control")))

# ==============================================================================
# [Step 3] Align Genes and initial filtering
# ==============================================================================
cat("[Step 3] Aligning Genes and Running TMM Normalization...\n")

# Strip Ensembl version
rownames(counts_tcga) <- sub("\\..*", "", rownames(counts_tcga))
rownames(counts_gtex) <- sub("\\..*", "", rownames(counts_gtex))
rownames(counts_srp)  <- sub("\\..*", "", rownames(counts_srp))

# Remove duplicate gene IDs caused by stripping versions
counts_tcga <- counts_tcga[!duplicated(rownames(counts_tcga)), ]
counts_gtex <- counts_gtex[!duplicated(rownames(counts_gtex)), ]
counts_srp  <- counts_srp[!duplicated(rownames(counts_srp)), ]

# Find common genes across all three datasets
common_genes <- intersect(intersect(rownames(counts_tcga), rownames(counts_gtex)), rownames(counts_srp))
cat(sprintf("  -> Aligned %d common genes across TCGA, GTEx, and SRP029880\n", length(common_genes)))

# Merge into one dataframe
counts_merged <- cbind(counts_tcga[common_genes, ], 
                       counts_gtex[common_genes, ], 
                       counts_srp[common_genes, ])

# Query the local database directly (no internet connection needed)
gene_annotations <- genes(EnsDb.Hsapiens.v86, return.type = "data.frame")

# Filter for protein coding genes
cat("  -> Fetching protein-coding gene annotations via local EnsDb...\n")
pc_genes_info <- subset(gene_annotations, gene_biotype == "protein_coding")
keep_ids <- intersect(rownames(counts_merged), pc_genes_info$gene_id)
counts_merged <- counts_merged[keep_ids, ]
pc_genes_info <- pc_genes_info[match(keep_ids, pc_genes_info$gene_id), ]

# Extract symbols and handle missing ones (replace blank/NA with Ensembl ID)
gene_symbols <- pc_genes_info$gene_name
missing_sym <- is.na(gene_symbols) | gene_symbols == ""
gene_symbols[missing_sym] <- pc_genes_info$gene_id[missing_sym]

# Aggregate counts for duplicate symbols (sums the reads of duplicate mappings)
counts_merged <- rowsum(counts_merged, group = gene_symbols)

cat(sprintf("  -> Retained %d unique protein-coding gene symbols\n", nrow(counts_merged)))

# ==============================================================================
# [Step 4] Filtering samples
# ==============================================================================
# Define Sample Group Annotations
tcga_sample_types <- tcga_cd$tcga.cgc_sample_sample_type[tcga_keep]
tcga_conditions   <- ifelse(tcga_sample_types == "Primary Tumor", "Tumor", "Control")

# Extract raw TCGA stage data
tcga_stage_raw <- tcga_cd$tcga.gdc_cases.diagnoses.tumor_stage[tcga_keep]
if(is.null(tcga_stage_raw)) {
  tcga_stage_raw <- tcga_cd$tcga.cgc_case_pathologic_stage[tcga_keep]
}

# Clean and convert Roman numeral stages to numeric 1-4
s_clean <- gsub("stage\\s*", "", tolower(as.character(tcga_stage_raw)))
tcga_numeric_stage <- rep(NA, length(s_clean))
tcga_numeric_stage[grepl("^iv", s_clean)] <- 4
tcga_numeric_stage[grepl("^iii", s_clean)] <- 3
tcga_numeric_stage[grepl("^ii", s_clean) & !grepl("^iii", s_clean)] <- 2
tcga_numeric_stage[grepl("^i", s_clean) & !grepl("^ii", s_clean) & !grepl("^iv", s_clean)] <- 1

sample_info <- data.frame(
  sample_id    = c(colnames(counts_tcga), colnames(counts_gtex), colnames(counts_srp)),
  dataset      = c(rep("TCGA", ncol(counts_tcga)), rep("GTEx", ncol(counts_gtex)), rep("SRP029880", ncol(counts_srp))),
  condition    = c(tcga_conditions, rep("Control", ncol(counts_gtex)), srp_conditions),
  cancer_stage = c(tcga_numeric_stage, rep(NA, ncol(counts_gtex)), rep(NA, ncol(counts_srp))),
  stringsAsFactors = FALSE
)
sample_info$cancer_stage[sample_info$condition == "Control"] <- 0
rownames(sample_info) <- sample_info$sample_id

# Keep samples that are NOT (Condition == Tumor AND Stage == NA), strictly enforcing this on TCGA only
keep_mask <- !(sample_info$condition == "Tumor" & is.na(sample_info$cancer_stage) & sample_info$dataset == "TCGA")

sample_info   <- sample_info[keep_mask, ]
counts_merged <- counts_merged[, keep_mask]

cat(sprintf("  -> Dropped %d un-staged TCGA tumor samples\n", sum(!keep_mask)))

# ==============================================================================
# [Step 5] Filtering Genes
# ==============================================================================
# Adaptive low-expression filtering via edgeR
design <- model.matrix(~ dataset + condition, data = sample_info)
dge <- DGEList(counts = counts_merged, samples = sample_info)
keep_genes <- filterByExpr(dge, design = design)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge, method = "TMM")
tmm_log2cpm <- cpm(dge, log = TRUE, prior.count = 1)

cat(sprintf("  -> Retained %d expressed genes after TMM normalization\n\n", nrow(tmm_log2cpm)))

# ------------------------------------------------------------------------------
# [Step 6] edgeR DEG Analysis (Uses raw counts + TMM factors + GLM)
# ------------------------------------------------------------------------------
dge_deg <- estimateDisp(dge, design = design)
fit_deg <- glmQLFit(dge_deg, design = design)
tr <- glmTreat(fit_deg, coef = "conditionTumor", lfc = 1)
deg_results <- topTags(tr, n = Inf)$table
deg_df <- data.frame(
  Gene_ID = rownames(deg_results),
  deg_results,
  check.names = FALSE
)

# ==============================================================================
# [Step 7] ComBat Batch Effect Correction
# ==============================================================================
# PCA verification plot (Pre-ComBat)
pca_res <- prcomp(t(tmm_log2cpm), scale. = TRUE)
pca_df  <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Dataset = sample_info$dataset,
  Condition = sample_info$condition
)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Dataset)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(title = "PCA of Pre-ComBat Expression (TMM log2CPM)",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)"))

cat("[Step 4] Performing ComBat Batch Correction...\n")

mod <- model.matrix(~ condition, data = sample_info)
combat_expr <- ComBat(dat = tmm_log2cpm, batch = sample_info$dataset, mod = mod)

# PCA verification plot (Post-ComBat)
pca_res_combat <- prcomp(t(combat_expr), scale. = TRUE)
pca_df_combat  <- data.frame(
  PC1 = pca_res_combat$x[, 1],
  PC2 = pca_res_combat$x[, 2],
  Dataset = sample_info$dataset,
  Condition = sample_info$condition
)
pca_plot_combat <- ggplot(pca_df_combat, aes(x = PC1, y = PC2, color = Condition, shape = Dataset)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(title = "PCA of ComBat-Corrected Expression (TMM log2CPM)",
       x = paste0("PC1 (", round(summary(pca_res_combat)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res_combat)$importance[2,2]*100, 1), "%)"))

# ==============================================================================
# [Step 8] Export csv files
# ==============================================================================
cat("[Step 5] Exporting CSV Files...\n")

# Ensure output directory exists before saving
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

expr_df <- data.frame(Gene_ID = rownames(tmm_log2cpm), tmm_log2cpm, check.names = FALSE)
combat_expr_df <- data.frame(Gene_ID = rownames(combat_expr), combat_expr, check.names = FALSE)

write.csv(combat_expr_df, paste(data_dir, "combat_tmm_log2cpm_expression.csv", sep="/"), row.names = FALSE)
write.csv(expr_df, paste(data_dir, "tmm_log2cpm_expression.csv", sep="/"), row.names = FALSE)
write.csv(sample_info, paste(data_dir, "sample_info.csv", sep="/"), row.names = FALSE)
write.csv(deg_df, paste(data_dir, "deg_results.csv", sep="/"), row.names = FALSE)

cat("\n=== Processing Complete! ===\n")
