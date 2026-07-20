library(GEOquery)
library(Biobase)
library(biomaRt)

# Define your GEO Accession Number (e.g., "GSE42568" or an RNA-seq GSE)
gse_id <- "GSE42568"

cat("Fetching dataset from GEO...\n")

# 1. Fetch GEO ExpressionSet & Sample Metadata
gse_data <- getGEO(gse_id, GSEMatrix = TRUE)

# Select the primary ExpressionSet object (usually the first platform)
eset <- gse_data[[1]]

# Extract Sample Metadata (PhenoData)
sample_metadata <- pData(eset)

# Extract Expression / Count Matrix
counts_matrix <- exprs(eset)

cat("Successfully retrieved matrix with", nrow(counts_matrix), "rows and", ncol(counts_matrix), "samples.\n")

# 2. Annotate Feature IDs to Ensembl IDs & Gene Symbols
# Connect to Ensembl BioMart
cat("Connecting to Ensembl BioMart for gene annotation...\n")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract probe/feature identifiers from the expression matrix
feature_ids <- rownames(counts_matrix)

# Map Affymetrix/Probe IDs or Entrez IDs to Ensembl IDs & Gene Symbols
# Note: Adjust 'filters' attribute if your input IDs are from a different array or standard Entrez IDs
annotation_results <- getBM(
  attributes = c("affy_hg_u133_plus_2", "ensembl_gene_id", "hgnc_symbol", "description"),
  filters    = "affy_hg_u133_plus_2",
  values     = feature_ids,
  mart       = mart
)

# 3. Merge Annotation Table with Expression Counts
# Convert expression matrix to a data frame with Probe IDs as a column
counts_df <- data.frame(
  affy_hg_u133_plus_2 = rownames(counts_matrix),
  counts_matrix,
  check.names = FALSE
)

# Merge matrix with gene annotations
annotated_counts <- merge(
  annotation_results,
  counts_df,
  by = "affy_hg_u133_plus_2",
  all.y = TRUE
)

annotated_counts <- annotated_counts[, !(names(annotated_counts) %in% c(
  "ensembl_gene_id",
  "hgnc_symbol",
  "description"
))]

# 4. Save Everything to CSV Files
cat("Saving data to CSV...\n")

# Matches 1 to 2 digits for day/month and 2 digits for year at the end of the string
sample_metadata$batch <- sub(".*_(\\d{1,2}_\\d{1,2}_\\d{2})$", "\\1", sample_metadata$title)

# Save data
dir.create(file.path("data", gse_id, "raw"), recursive = TRUE, showWarnings = FALSE)
write.csv(sample_metadata, file = file.path("data", gse_id, "raw", "samples.csv"), row.names = FALSE)
write.csv(annotated_counts, file = file.path("data", gse_id, "raw", "counts.csv"), row.names = FALSE)

cat("Done! Output saved as CSV files in your current working directory.\n")
