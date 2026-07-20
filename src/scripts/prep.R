if (!require("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse")
}
library(tidyverse)

if (!require("edgeR", quietly = TRUE)) {
    BiocManager::install("edgeR")
}
library(edgeR)

if (!require("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
}
library(biomaRt)

if (!require("sva", quietly = TRUE)) {
    BiocManager::install("sva")
}
library(sva)

df_rna <- read.csv("data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/raw_expression.csv", check.names = FALSE)
metadata <- read.csv("data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/samples.csv")

rownames(df_rna) <- df_rna$gene_id
df_rna <- df_rna %>% dplyr::select(-gene_id)

metadata <- metadata[match(colnames(df_rna), metadata$SAMPID), ]

df_rna_matrix <- as.matrix(df_rna)
batch_counts <- table(metadata$SMGEBTCH)
valid_batches <- names(batch_counts[batch_counts > 1])

samples_to_drop <- sum(batch_counts[batch_counts == 1])
print(paste("Dropping", samples_to_drop, "samples because they are the only sample in their sequencing batch."))

metadata <- metadata %>%
    filter(SMGEBTCH %in% valid_batches)

df_rna_matrix <- df_rna_matrix[, metadata$SAMPID]

batches <- droplevels(as.factor(metadata$SMGEBTCH))
groups <- droplevels(as.factor(metadata$SMTS))

print(paste("Transcripts before filtering:", nrow(df_rna_matrix)))

dge_temp <- DGEList(counts = df_rna_matrix, group = groups)

keep_transcripts <- filterByExpr(dge_temp)

df_rna_matrix_filtered <- df_rna_matrix[keep_transcripts, ]

print(paste("Transcripts after filtering:", nrow(df_rna_matrix_filtered)))

write.csv(df_rna_matrix_filtered, file = "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/rna_matrix_filtered.csv", row.names = TRUE, quote = FALSE)
write.csv(metadata, file = "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/samples_filtered.csv", row.names = FALSE, quote = TRUE)

print("Running ComBat_seq")
chunk_size <- 1000
num_transcripts <- nrow(df_rna_matrix_filtered)
num_chunks <- ceiling(num_transcripts / chunk_size)

corrected_chunks <- list()

for (i in 1:num_chunks) {
    start_row <- ((i - 1) * chunk_size) + 1
    end_row <- min(i * chunk_size, num_transcripts)

    print(paste("Processing Chunk", i, "of", num_chunks, "(Rows", start_row, "to", end_row, ")..."))

    current_chunk <- df_rna_matrix_filtered[start_row:end_row, , drop = FALSE]
    chunk_corrected <- ComBat_seq(counts = current_chunk, batch = batches, group = groups)
    corrected_chunks[[i]] <- chunk_corrected
    gc()
}
df_rna_corrected <- do.call(rbind, corrected_chunks)
write.csv(df_rna_corrected, file = "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/rna_matrix_corrected.csv", row.names = TRUE, quote = FALSE)
print("Chunked ComBat_seq completed")

mart <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    mirror = "useast"
)
mapping <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "start_position", "end_position"),
    filters = "ensembl_gene_id",
    values = rownames(df_rna_corrected),
    mart = mart
)

mapping$gene_length <- mapping$end_position - mapping$start_position + 1

df_rna_corrected_df <- as.data.frame(df_rna_corrected)
df_rna_corrected_df$gene_id <- rownames(df_rna_corrected_df)
df_rna_combined <- df_rna_corrected_df %>%
    inner_join(mapping, by = c("gene_id" = "ensembl_gene_id"))
rownames(df_rna_combined) <- df_rna_combined$gene_id

df_rna_combined$length_kb <- df_rna_combined$gene_length / 1000

sample_cols <- metadata$SAMPID
fpk <- df_rna_combined[, sample_cols] / df_rna_combined$length_kb

rna_getmm <- DGEList(counts = fpk)
rna_getmm <- calcNormFactors(rna_getmm, method = "TMM")
rna_getmm <- cpm(rna_getmm)

rna_getmm <- as.data.frame(rna_getmm)
rna_getmm$gene_id <- rownames(rna_getmm)

rna_getmm <- rna_getmm %>% dplyr::select(gene_id, dplyr::everything())

saveRDS(list("rna_getmm" = rna_getmm, "metadata" = metadata), file = "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/gtex_gene_getmm.rds")
write.csv(rna_getmm, file = "data/for_SWAMP/expression_datasets/GTEx_agg_gene_rna/data_gtex_gene_getmm.csv", row.names = FALSE, quote = FALSE)

print("GeTMM Pipeline completed")
