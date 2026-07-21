# Load necessary library
library(dplyr)
library(ggplot2)

path <- "GSE141910"

# Load the data
shapley <- read.csv(file.path("out", path, "shapley.csv"), row.names = 1)
deg <- read.csv(file.path("out", path, "deg.csv"), row.names = 1)

# Merge the data on rownames
merged_data <- merge(shapley, deg, by = "row.names", all = TRUE)
rownames(merged_data) <- merged_data$Row.names
merged_data$Row.names <- NULL
