# ==============================================================================
# 0. DEPENDENCY MANAGEMENT & LIBRARIES
# ==============================================================================
required_cran <- c("tidyverse")
required_bioc <- c("edgeR", "biomaRt", "sva", "GenomicRanges", "IRanges")

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
# 1. DATA LOADING & ENSEMBL VERSION CLEANUP
# ==============================================================================
cat("Loading datasets...\n")
df_rna <- read.csv("data/GSE42568/raw/counts.csv", check.names = FALSE)
metadata <- read.csv("data/GSE42568/raw/samples.csv", check.names = TRUE)

# Strip Ensembl ID version numbers (e.g., ENSG00000000003.14 -> ENSG00000000003)
# df_rna$ensembl_gene_id <- sub("\\.\\d+$", "", df_rna$ensembl_gene_id)

# Deduplicate if stripping version numbers resulted in non-unique Ensembl IDs
if (any(duplicated(df_rna$ensembl_gene_id))) {
  cat("Consolidating duplicate Ensembl IDs post-version stripping...\n")
  df_rna <- df_rna %>%
    group_by(ensembl_gene_id) %>%
    summarise(across(everything(), sum)) %>%
    ungroup()
}

df_rna <- df_rna %>% column_to_rownames("ensembl_gene_id")

# Align metadata with expression matrix columns
metadata <- metadata[match(colnames(df_rna), metadata$sample_id), ]
df_rna_matrix <- as.matrix(df_rna)

# ==============================================================================
# 2. FILTER SINGLE-SAMPLE BATCHES & LOW-EXPRESSION GENES
# ==============================================================================
batch_counts <- table(metadata$batch)
valid_batches <- names(batch_counts[batch_counts > 1])

samples_to_drop <- sum(batch_counts[batch_counts == 1])
cat(paste("Dropping", samples_to_drop, "samples belong to single-sample sequencing batches.\n"))

metadata <- metadata %>% filter(batch %in% valid_batches)
df_rna_matrix <- df_rna_matrix[, metadata$sample_id]

batches <- droplevels(as.factor(metadata$batch))
groups <- droplevels(as.factor(metadata$characteristics_ch1))

cat(paste("Transcripts before filterByExpr:", nrow(df_rna_matrix), "\n"))
dge_temp <- DGEList(counts = df_rna_matrix, group = groups)
keep_transcripts <- filterByExpr(dge_temp)
df_rna_matrix_filtered <- df_rna_matrix[keep_transcripts, ]
cat(paste("Transcripts after filterByExpr:", nrow(df_rna_matrix_filtered), "\n"))

# Save intermediate filtered matrices
write.csv(df_rna_matrix_filtered, "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/rna_matrix_filtered.csv", quote = FALSE)
write.csv(metadata, "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/samples_filtered.csv", row.names = FALSE)

# ==============================================================================
# 3. BATCH CORRECTION (ComBat_seq)
# ==============================================================================
cat("Running ComBat_seq batch correction...\n")
chunk_size <- 1000
num_transcripts <- nrow(df_rna_matrix_filtered)
num_chunks <- ceiling(num_transcripts / chunk_size)

corrected_chunks <- list()

for (i in 1:num_chunks) {
  start_row <- ((i - 1) * chunk_size) + 1
  end_row <- min(i * chunk_size, num_transcripts)

  cat(paste("Processing Chunk", i, "of", num_chunks, "(Rows", start_row, "to", end_row, ")...\n"))

  current_chunk <- df_rna_matrix_filtered[start_row:end_row, , drop = FALSE]
  chunk_corrected <- ComBat_seq(counts = current_chunk, batch = batches, group = groups)
  corrected_chunks[[i]] <- chunk_corrected
  gc()
}

df_rna_corrected <- do.call(rbind, corrected_chunks)
write.csv(df_rna_corrected, "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/rna_matrix_corrected.csv", quote = FALSE)

# ==============================================================================
# 4. EXON-BASED NON-OVERLAPPING GENE LENGTH CALCULATION
# ==============================================================================
cat("Fetching exon boundaries from Ensembl to calculate non-overlapping exonic lengths...\n")

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Retrieve individual exon start and end positions
exon_coords <- getBM(
  attributes = c("ensembl_gene_id", "exon_chrom_start", "exon_chrom_end"),
  filters    = "ensembl_gene_id",
  values     = rownames(df_rna_corrected),
  mart       = mart
)

# Use GenomicRanges / IRanges to merge overlapping exons per gene
cat("Reducing overlapping exon coordinates per gene...\n")
gr <- GRanges(
  seqnames = "chr1", # Dummy seqname to construct GRanges
  ranges   = IRanges(start = exon_coords$exon_chrom_start, end = exon_coords$exon_chrom_end),
  gene_id  = exon_coords$ensembl_gene_id
)

# Split by gene_id, reduce overlaps, and calculate total exonic length
gr_split <- split(gr, gr$gene_id)
gr_reduced <- reduce(gr_split)
exonic_lengths <- sum(width(gr_reduced))

# Construct a clean mapping dataframe
gene_length_df <- data.frame(
  gene_id = names(exonic_lengths),
  gene_length = as.vector(exonic_lengths),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 5. GeTMM NORMALIZATION (FPK -> TMM -> CPM)
# ==============================================================================
cat("Calculating GeTMM normalized expression...\n")

# Merge corrected count matrix with true exonic lengths
df_rna_corrected_df <- as.data.frame(df_rna_corrected)
df_rna_corrected_df$gene_id <- rownames(df_rna_corrected_df)

df_rna_combined <- df_rna_corrected_df %>%
  inner_join(gene_length_df, by = "gene_id")

# Set clean rownames and extract sample columns
rownames(df_rna_combined) <- df_rna_combined$gene_id
sample_cols <- metadata$sample_id

# 1. Calculate length in kilobases
length_kb <- df_rna_combined$gene_length / 1000

# 2. Calculate FPK (Fragments Per Kilobase) using sweep for matrix-vector division
fpk_matrix <- sweep(as.matrix(df_rna_combined[, sample_cols]), 1, length_kb, FUN = "/")

# 3. Apply TMM normalization across samples & calculate CPM
dge_fpk <- DGEList(counts = fpk_matrix)
dge_fpk <- calcNormFactors(dge_fpk, method = "TMM")
rna_getmm <- cpm(dge_fpk)

# Format final dataframe
rna_getmm_df <- as.data.frame(rna_getmm) %>%
  rownames_to_column("gene_id")

# ==============================================================================
# 6. SAVE OUTPUTS
# ==============================================================================
cat("Saving final outputs...\n")
saveRDS(
  list("rna_getmm" = rna_getmm_df, "metadata" = metadata),
  file = "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/gtex_gene_getmm.rds"
)

write.csv(
  rna_getmm_df,
  file = "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/data_gtex_gene_getmm.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("Pipeline completed successfully!\n")
