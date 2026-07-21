# ==============================================================================
# 0. DEPENDENCY MANAGEMENT & LIBRARIES
# ==============================================================================
required_cran <- c("tidyverse", "matrixStats")
required_bioc <- c("limma", "sva", "biomaRt")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

for (pkg in required_cran) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}

for (pkg in required_bioc) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(pkg, update = FALSE)
        library(pkg, character.only = TRUE)
    }
}

# ==============================================================================
# 1. DATA LOADING & DEDUPLICATION
# ==============================================================================
gse_id <- "GSE42568"
cat("Loading microarray expression dataset", gse_id, "...\n")

df_microarray <- read.csv(file = file.path("data", gse_id, "raw", "counts.csv"), check.names = FALSE)
metadata <- read.csv(file = file.path("data", gse_id, "raw", "samples.csv"), check.names = TRUE)

# Determine the feature/probe column name
if ("probe_id" %in% colnames(df_microarray)) {
    id_col <- "probe_id"
} else if ("ensembl_gene_id" %in% colnames(df_microarray)) {
    id_col <- "ensembl_gene_id"
} else {
    id_col <- colnames(df_microarray)[1] # Default to first column
}

# Clean missing/empty IDs
df_microarray <- df_microarray %>%
    filter(!is.na(.data[[id_col]]) & .data[[id_col]] != "")

# Aggregate duplicate probe IDs by taking the mean across numeric sample columns
if (any(duplicated(df_microarray[[id_col]]))) {
    cat("Found duplicate probe/feature IDs. Averaging signal values across duplicates...\n")
    df_microarray <- df_microarray %>%
        group_by(.data[[id_col]]) %>%
        summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
        ungroup()
}

# Now safe to set row names
df_microarray <- df_microarray %>% column_to_rownames(id_col)

# Align metadata with expression matrix columns
metadata <- metadata[match(colnames(df_microarray), metadata$sample_id), ]
df_matrix <- as.matrix(df_microarray)

# Ensure data is log2 transformed
if (max(df_matrix, na.rm = TRUE) > 100) {
    cat("Data appears untransformed. Applying log2 transformation...\n")
    df_matrix <- log2(df_matrix + 1)
}

# ==============================================================================
# 2. FILTER SINGLE-SAMPLE BATCHES & LOW VARIANCE PROBES
# ==============================================================================
batch_counts <- table(metadata$batch)
valid_batches <- names(batch_counts[batch_counts > 1])

samples_to_drop <- sum(batch_counts[batch_counts == 1])
cat(paste("Dropping", samples_to_drop, "samples belonging to single-sample batches.\n"))

metadata <- metadata %>% filter(batch %in% valid_batches)
df_matrix <- df_matrix[, metadata$sample_id]

# Microarray Filtering: Drop non-informative low-variance probes (bottom 20%)
row_variances <- rowVars(df_matrix)
variance_threshold <- quantile(row_variances, probs = 0.20, na.rm = TRUE)
keep_probes <- row_variances > variance_threshold

df_matrix_filtered <- df_matrix[keep_probes, ]
cat(paste("Features before filtering:", nrow(df_matrix), "\n"))
cat(paste("Features after variance filtering:", nrow(df_matrix_filtered), "\n"))

# Save intermediate filtered matrices
dir.create(file.path("data", gse_id, "processed"), recursive = TRUE, showWarnings = FALSE)
write.csv(df_matrix_filtered, file = file.path("data", gse_id, "processed", "matrix_filtered.csv"), quote = FALSE)
write.csv(metadata, file = file.path("data", gse_id, "processed", "samples_filtered.csv"), row.names = FALSE)

# ==============================================================================
# 3. BATCH CORRECTION FOR CONTINUOUS DATA (limma::removeBatchEffect)
# ==============================================================================
cat("Running batch correction using limma::removeBatchEffect...\n")

batches <- droplevels(as.factor(metadata$batch))
groups <- droplevels(as.factor(metadata$characteristics_ch1))

# Design matrix to preserve biological group variance while removing batch effect
design_matrix <- model.matrix(~groups)

# Apply removeBatchEffect (designed specifically for continuous log-intensity values)
df_matrix_corrected <- removeBatchEffect(
    x = df_matrix_filtered,
    batch = batches,
    design = design_matrix
)

write.csv(df_matrix_corrected, file = file.path("data", gse_id, "processed", "matrix_corrected.csv"), quote = FALSE)

# ==============================================================================
# 4. ANNOTATE PROBES & DEDUPLICATE TO ENSEMBL IDs
# ==============================================================================
cat("Annotating probe IDs to Ensembl IDs...\n")

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

feature_ids <- rownames(df_matrix_corrected)

# Query Ensembl mapping (adjust filter attribute if using a different array)
annotations <- getBM(
    attributes = c("affy_hg_u133_plus_2", "ensembl_gene_id", "hgnc_symbol"),
    filters    = "affy_hg_u133_plus_2",
    values     = feature_ids,
    mart       = mart
)

# Filter out unmapped probes or empty Ensembl IDs
annotations <- annotations %>%
    filter(!is.na(ensembl_gene_id) & ensembl_gene_id != "")

# Convert corrected matrix to long format / dataframe with probe IDs
df_corrected_df <- as.data.frame(df_matrix_corrected) %>%
    rownames_to_column("affy_hg_u133_plus_2")

# Merge annotations with corrected expression values
df_mapped <- annotations %>%
    inner_join(df_corrected_df, by = "affy_hg_u133_plus_2")

# ------------------------------------------------------------------------------
# DEDUPLICATION: Collapse multiple probes per Ensembl ID
# Select the probe with the highest mean expression across all samples
# ------------------------------------------------------------------------------
cat("Collapsing duplicate probes to unique Ensembl IDs (selecting max mean signal)...\n")

sample_cols <- metadata$sample_id

df_collapsed <- df_mapped %>%
    # Calculate mean expression per probe across samples
    mutate(mean_expr = rowMeans(across(all_of(sample_cols)), na.rm = TRUE)) %>%
    # Group by Ensembl ID and keep the probe with highest mean expression
    group_by(ensembl_gene_id) %>%
    slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(-mean_expr, -affy_hg_u133_plus_2)

# Set Ensembl ID as the clean matrix row index / rownames
df_final_matrix <- df_collapsed %>%
    dplyr::select(-hgnc_symbol) %>%
    column_to_rownames("ensembl_gene_id")

# ==============================================================================
# 5. SAVE FINAL OUTPUTS
# ==============================================================================
cat("Saving final outputs...\n")

saveRDS(
    list(
        "expression_matrix" = df_final_matrix,
        "metadata" = metadata
    ),
    file = file.path("data", gse_id, "processed", "matrix_final.rds")
)

write.csv(
    df_final_matrix,
    file = file.path("data", gse_id, "processed", "matrix_final.csv"),
    row.names = TRUE,
    quote = FALSE
)

cat("Pipeline completed successfully! Clean matrix indexed by Ensembl IDs.\n")
