source("functions.R")

path <- "GSE141910"

data_norm <- read.csv(file.path("data", path, "processed", "counts_norm.csv"), row.names = 1)
data_sample <- read.csv(file.path("data", path, "processed", "sample.csv"), row.names = 1)

A_SR <- data_norm[, data_sample$etiology == "NF"] # expression levels of reference samples
A_SD <- data_norm[, data_sample$etiology == "DCM"] # expression levels of disease samples
SR <- colnames(A_SR) # reference samples
SD <- colnames(A_SD) # disease samples
num_SR <- length(SR) # number of reference samples
num_SD <- length(SD) # number of disease samples
genes <- rownames(data_norm) # genes
n <- length(genes) # number of genes

# create boolean expression matrix showing which genes are differentially expressed in which disease samples
B <- boolean_expression_matrix(A_SD, A_SR, 0, 100)

# determine the support of B
# each element in sp_B corresponds to a disease sample and contains the set of genes that are differentially expressed in that sample
sp_B <- support_of(B)

# find all unique coalitions from the support sets
coalitions <- find_coalitions(sp_B)
num_coalitions <- length(coalitions)

# compute lambda values for each coalition
lambda_bar <- numeric(num_coalitions)
for (j in 1:num_coalitions) {
    S <- coalitions[[j]]
    for (k in 1:length(SD)) {
        T <- sp_B[[k]]
        if (identical(T, S)) {
            lambda_bar[j] <- lambda_bar[j] + 1
        }
    }
}
lambda <- lambda_bar / num_SD

# compute Shapley values for each gene
shapley_value <- numeric(n)
names(shapley_value) <- genes
for (gene in genes) {
    cat(sprintf("\rIteration %d/%d", which(genes == gene), n))
    flush.console()
    for (i in 1:num_coalitions) {
        S <- coalitions[[i]]
        if (gene %in% S) {
            shapley_value[gene] <- shapley_value[gene] + lambda[i] / length(S)
        }
    }
}
shapley_value <- sort(shapley_value, decreasing = TRUE)

# save the Shapley values
dir.create(file.path("out", path), recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(shapley = shapley_value), file = file.path("out", path, "shapley.csv"), row.names = TRUE)
