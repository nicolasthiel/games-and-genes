library(DESeq2)

path <- "GSE141910"

# read in data
data_counts <- read.csv(file.path("data", path, "raw", "counts.csv"), row.names = 1)
data_sample <- read.csv(file.path("data", path, "raw", "sample.csv"), row.names = 1)

# only keep NF and DCM samples
data_sample <- data_sample[data_sample$etiology %in% c("DCM", "NF"), ]
data_counts <- data_counts[, rownames(data_sample)]

# create DESeq object
dds <- DESeqDataSetFromMatrix(
    countData = data_counts,
    colData = data_sample,
    design = ~etiology
)

smallestGroupSize <- 166
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

# estimate size factors
dds <- estimateSizeFactors(dds)

# normalize
norm <- counts(dds, normalize = TRUE)
data_norm <- as.data.frame(norm)

# log2 transform
data_counts_log <- log2(data_norm + 1)
data_norm_log <- log2(data_norm + 1)

# save normalised data
write.csv(data_norm, file = file.path("data", path, "processed", "counts_norm.csv"))
write.csv(data_counts_log, file = file.path("data", path, "processed", "counts_log.csv"))
write.csv(data_sample, file = file.path("data", path, "processed", "sample.csv"))


# Perform differential expression analysis
dds <- DESeq(dds)

# Extract results for DCM vs NF
res <- results(dds)

# Order results by adjusted p-value
res <- res[order(res$padj, na.last = NA), ]

# Save results to a CSV file
write.csv(as.data.frame(res), file = file.path("out", path, "deg.csv"))
