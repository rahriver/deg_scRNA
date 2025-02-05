# # Load necessary libraries
# library(Seurat)
# library(Signac)
# 
# # Define the file path
# file_path <- "/home/core/hyp/hyp.rds"
# 
# 
# if (!requireNamespace("Signac", quietly = TRUE)) {
#     install.packages("remotes")
#     remotes::install_github("timoast/signac")
# }
# 
# 
# # Load the Seurat object
# hyp <- readRDS(file_path)
# 
# # Print the structure of the object
# print(hyp)
# 
# # View metadata
# head(hyp@meta.data)
# 
# Assays(hyp)
# 
# head(GetAssayData(hyp, slot = "counts"))
# 
# Idents(hyp)
# 
# summary(hyp@meta.data)
#================================================================#
# CANT RUN TOO MUCH RAM USAGE
#================================================================#
# library(Seurat)
# library(DESeq2)
# library(tidyverse)
# 
# install.packages(c("Matrix", "data.table", "pheatmap"))
# 
# 
# file_path <- "/home/core/hyp/hyp.rds"
# hyp <- readRDS(file_path)
# 
# # extract metadata
# meta_data <- hyp@meta.data
# 
# # extract raw RNA counts
# rna_counts <- GetAssayData(hyp, assay = "RNA", slot = "counts")
# 
# # group metadata by species
# species_labels <- meta_data$Species  # Modify this if another column contains species names
# 
# # pseudo-bulk counts (sum counts per species)
# pseudo_bulk_counts <- aggregate(Matrix::t(rna_counts), by = list(Species = species_labels), FUN = sum)
# rownames(pseudo_bulk_counts) <- pseudo_bulk_counts$Species
# pseudo_bulk_counts <- pseudo_bulk_counts[, -1]  # Remove species column
# 
# pseudo-bulk count matrix
# print(dim(pseudo_bulk_counts))
# 
# 
# library(DESeq2)
# 
# meta_df <- data.frame(Species = rownames(pseudo_bulk_counts))
# rownames(meta_df) <- meta_df$Species
# 
# dds <- DESeqDataSetFromMatrix(
#     countData = t(pseudo_bulk_counts),  # Transpose so genes are rows
#     colData = meta_df,
#     design = ~ Species  # Model species-specific effects
# )
# 
# dds <- DESeq(dds)
# 
# # Mbuna vs Non-Mbuna
# dds$Mbuna <- ifelse(dds$Species %in% c("MbunaSpecies1", "MbunaSpecies2"), "Mbuna", "NonMbuna")
# 
# # Re-run DESeq2 with new contrast
# dds <- DESeq(dds)
# 
# res <- results(dds, contrast = c("Mbuna", "Mbuna", "NonMbuna"))
# 
# res <- res[order(res$padj), ]
# # Top ones
# head(res)
# 
# library(ggplot2)
# 
# res_df <- as.data.frame(res)
# 
# # Volcano plot
# ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
#     geom_point(aes(color = padj < 0.05), alpha = 0.6) +
#     theme_minimal() +
#     ggtitle("Differential Expression: Mbuna vs Non-Mbuna") +
#     xlab("Log2 Fold Change") + ylab("-log10 Adjusted P-value") +
#     scale_color_manual(values = c("grey", "red")) +
#     theme(legend.position = "none")
#================================================================#
#================================================================#
library(Seurat)
library(Matrix)
library(data.table)
library(DESeq2)
library(limma)
library(ggplot2)
library(pheatmap)
library(future.apply)

# just for my laptop spec
plan(multicore, workers = 4)

file_path <- "/home/core/hyp/hyp.rds"
hyp <- readRDS(file_path)

# raw RNA counts as a sparse matrix (genes x cells)
rna_counts <- GetAssayData(hyp, assay = "RNA", slot = "counts")
meta_data <- hyp@meta.data

# species labels (one label per cell)
species_labels <- factor(meta_data$Species)

################################################################################
# Pseudo-bulk by species
################################################################################
# Summary workflow:
# we're aggregating our single-cell data into pseudobulk by summing the RNA counts across all cells within each species.
# we're also doing sparse matrix: rows are genes and columns are cells.
# adding a new column 'Species' to the matrix
# we group all the cells belonging to the same species and sums up their RNA counts for each gene.
# for each species we have a single pseudobulk profile that counts from all cells of that species are summed.

# transposing to get a matrix of dimension (cells x genes)
rna_counts_t <- as.data.table(as.matrix(t(rna_counts)))

# adding a 'Species' column
rna_counts_t[, Species := species_labels]

# using sum counts within each species
pseudo_bulk_counts <- rna_counts_t[, lapply(.SD, sum), by = Species]

# converting to a matrix. its now (species x genes)
rownames(pseudo_bulk_counts) <- pseudo_bulk_counts$Species
pseudo_bulk_counts <- pseudo_bulk_counts[, -1, with = FALSE]
pseudo_bulk_counts <- as(as.matrix(pseudo_bulk_counts), "dgCMatrix")
print(dim(pseudo_bulk_counts))

# so DESeq2 needs a matrix with genes in rows and samples in columns.
# then I transpose to get (genes x species).
pseudo_bulk_counts <- t(pseudo_bulk_counts)

# metadata coldata
meta_df <- data.frame(Species = colnames(pseudo_bulk_counts))
rownames(meta_df) <- meta_df$Species
stopifnot(all(colnames(pseudo_bulk_counts) == rownames(meta_df)))

dds <- DESeqDataSetFromMatrix(
    countData = pseudo_bulk_counts,
    colData   = meta_df,
    design    = ~ Species
)

dds <- DESeq(dds)

################################################################################
# Mbuna vs Non-Mbuna
################################################################################

# Define Mbuna vs NonMbuna
dds$Mbuna <- ifelse(dds$Species %in% c("MbunaSpecies1", "MbunaSpecies2"), 
                    "Mbuna", 
                    "NonMbuna")

# again DESeq2 but with Mbuna factor
dds <- DESeq(dds)

# DGE
res <- results(dds, contrast = c("Mbuna", "Mbuna", "NonMbuna"))

# filter
res <- res[order(res$padj), ]

head(res)

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < 0.05), alpha = 0.6) +
    theme_minimal() +
    ggtitle("Mbuna vs Non-Mbuna") +
    xlab("Log2 Fold Change") + ylab("-log10 Adjusted P-value") +
    scale_color_manual(values = c("grey", "red")) +
    theme(legend.position = "none")