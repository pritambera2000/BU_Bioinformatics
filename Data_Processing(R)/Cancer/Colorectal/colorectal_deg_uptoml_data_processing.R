library(affy)
library(affyio)
library(GEOquery)
library(dplyr)
library(tidyverse)
library(limma)
# Reading the cel files
raw_data_colorectal <- ReadAffy(celfile.path = "cel_files_colorectal/")
# Performing normalization
normalized_data_colorectal <- rma(raw_data_colorectal)
# Get expression estimates
normalized_expr_colorectal <- as.data.frame(exprs(normalized_data_colorectal))
# Map probe IDs to gene symbols
gse_colorectal <- getGEO("GSE44076", GSEMatrix = TRUE)
# Fetch feature data to get ID - gene symbol mapping
eset_colorectal <- gse_colorectal[[1]]
# Extract metadata
metadata_colorectal <- pData(eset_colorectal)
# subset the metadata
metadata_colorectal <- metadata_colorectal[,c(2,10)]
# Extract the GSM IDs from the col names
gsm_ids <- sapply(colnames(normalized_expr_colorectal), function(x) strsplit(x, "_")[[1]][1])
# Replace column names
colnames(normalized_expr_colorectal) <- gsm_ids
# Convert the metadata to a data frame
metadata_df <- as.data.frame(metadata_colorectal)
# Remove the "sample type:" prefix from the characteristics_ch1 column
metadata_df$characteristics_ch1 <- sapply(strsplit(metadata_df$characteristics_ch1, ": "), `[`, 2)
# Remove the geo_accession column from the metadata
metadata_df <- subset(metadata_df, select = -geo_accession)

# making sure the row names in colData matches to column names in counts_data
all(colnames(normalized_expr_colorectal) %in% rownames(metadata_df)) ## TRUE
# are they in the same order?
all(colnames(normalized_expr_colorectal) == rownames(metadata_df)) ## TRUE

## DEG
# Create design matrix directly from metadata dataframe
design <- model.matrix(~ 0 + factor(metadata_df$characteristics_ch1))
# Set column names of design matrix
colnames(design) <- levels(factor(metadata_df$characteristics_ch1))
# Create contrast matrix
contrasts <- makeContrasts(
  Tumor_vs_Normal = Tumor - Normal,
  Mucosa_vs_Normal = Mucosa - Normal,
  Mucosa_vs_Tumor = Mucosa - Tumor,
  levels = colnames(design)
)

# Fit the linear model
fit <- lmFit(normalized_expr_colorectal, design)

# Apply empirical Bayes smoothing and create DGELRT objects
tumor_vs_normal_dgelrt <- eBayes(contrasts.fit(fit, contrasts[, "Tumor_vs_Normal"]))
mucosa_vs_normal_dgelrt <- eBayes(contrasts.fit(fit, contrasts[, "Mucosa_vs_Normal"]))
mucosa_vs_tumor_dgelrt <- eBayes(contrasts.fit(fit, contrasts[, "Mucosa_vs_Tumor"]))

# Get top tables
tumor_vs_normal_top <- topTable(tumor_vs_normal_dgelrt, n = Inf)
mucosa_vs_normal_top <- topTable(mucosa_vs_normal_dgelrt, n = Inf)
mucosa_vs_tumor_top <- topTable(mucosa_vs_normal_dgelrt, n = Inf)

# Filter DEGs by significance (example: using adjusted p-value)
sig_cutoff <- 0.05
tumor_vs_normal_sig <- tumor_vs_normal_top[tumor_vs_normal_top$adj.P.Val < sig_cutoff, ]
mucosa_vs_normal_sig <- mucosa_vs_normal_top[mucosa_vs_normal_top$adj.P.Val < sig_cutoff, ]
mucosa_vs_tumor_sig <- mucosa_vs_tumor_top[mucosa_vs_tumor_top$adj.P.Val < sig_cutoff, ]

# Select top significant genes based on adjusted p-value
num_features <- 100
top_tumor_vs_normal_genes <- head(rownames(tumor_vs_normal_sig), num_features)
top_mucosa_vs_normal_genes <- head(rownames(mucosa_vs_normal_sig), num_features)
top_mucosa_vs_tumor_genes <- head(rownames(mucosa_vs_tumor_sig), num_features)

# Combine significant genes from all comparisons
significant_genes <- unique(c(top_tumor_vs_normal_genes, top_mucosa_vs_normal_genes, top_mucosa_vs_tumor_genes))

# Subset normalized expression data to include only significant genes
significant_expr <- normalized_expr_colorectal[significant_genes, ]
# Transpose the significant expression data
significant_expr <- t(significant_expr)

# Define the directory where you want to save the files
output_directory <- "colorectal_deg"

# Save significant genes to a file
write.csv(significant_expr, file.path(output_directory, "significant_genes.csv"), row.names = TRUE)

## Process the DEG Data for ML ##

## Re annotate the metadata_df with the sample_type column
metadata_df_annot <- as.data.frame(metadata_colorectal)
# Remove the "sample type:" prefix from the characteristics_ch1 column
metadata_df_annot$characteristics_ch1 <- sapply(strsplit(metadata_df_annot$characteristics_ch1, ": "), `[`, 2)
# Match the GSM IDs with the metadata
matched_metadata <- metadata_df_annot[match(gsm_ids, metadata_df_annot$geo_accession), ]

# Convert transposed_normalized_expr_colorectal to a data frame
significant_expr_df <- as.data.frame(significant_expr)

# Add the sample types to the transposed expression data
significant_expr_df$sample_type <- matched_metadata$characteristics_ch1

# Save significant genes to a file
write.csv(significant_expr_df, file.path(output_directory, "colorectal_ml_data_deg.csv"), row.names = FALSE)
