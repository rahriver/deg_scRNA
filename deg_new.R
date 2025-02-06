library(Seurat)
library(dplyr)
library(ggplot2)

# Ensure the active assay is RNA
DefaultAssay(seurat_obj) <- "RNA"

# Check unique values in Mbuna_nonMbuna
unique(seurat_obj@meta.data$Mbuna_nonMbuna)

# Set identities based on Mbuna vs. Non-Mbuna
Idents(seurat_obj) <- "Mbuna_nonMbuna"

# Perform DGE: Compare Mbuna vs. Non-Mbuna
dge_results <- FindMarkers(
    seurat_obj,
    ident.1 = "mbuna",
    ident.2 = "non-mbuna",
    test.use = "MAST",  # Change to "MAST" if needed
    min.pct = 0.1,
    logfc.threshold = 0.25
)

# Calculate the proportion of zero counts per gene
zero_fraction <- rowSums(GetAssayData(seurat_obj, slot = "counts") == 0) / ncol(seurat_obj)

# Plot histogram of zero fractions
hist(zero_fraction, breaks = 50, col = "blue", main = "Proportion of Zero-Expression Genes", xlab = "Fraction of Zeros")

# Find overlapping genes between Wilcoxon and MAST results
wilcox_deg <- rownames(FindMarkers(seurat_obj, ident.1 = "mbuna", ident.2 = "non-mbuna", test.use = "wilcox"))
mast_deg <- rownames(FindMarkers(seurat_obj, ident.1 = "mbuna", ident.2 = "non-mbuna", test.use = "MAST"))

# Find genes detected by only one method
wilcox_only <- setdiff(wilcox_deg, mast_deg)
mast_only <- setdiff(mast_deg, wilcox_deg)
both_methods <- intersect(wilcox_deg, mast_deg)

# Print summary
cat("Wilcoxon Only DEGs:", length(wilcox_only), "\n")
cat("MAST Only DEGs:", length(mast_only), "\n")
cat("Shared DEGs:", length(both_methods), "\n")



# View top results
head(dge_results)

# Save results
write.csv(dge_results, "DGE_results_mbuna_vs_nonmbuna_MAST.csv", row.names = TRUE)

# Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(dge_results,
                lab = rownames(dge_results),
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Differential Expression: Mbuna vs. Non-Mbuna",
                pCutoff = 0.05,
                FCcutoff = 0.5)

# Violin Plot for top differentially expressed gene
top_gene <- rownames(dge_results)[1]  # Get the top gene
VlnPlot(seurat_obj, features = top_gene, group.by = "Mbuna_nonMbuna")
