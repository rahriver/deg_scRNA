# # ################################################################################
# # # 1) Load libraries
# # ################################################################################
# # library(Seurat)
# # library(glmmSeq)  # Make sure you've installed glmmSeq and dependencies
# # 
# # ################################################################################
# # # 2) Load Seurat object
# # ################################################################################
# # file_path <- "/home/core/hyp/hyp.rds"
# # hyp <- readRDS(file_path)
# # 
# # ################################################################################
# # # 3) Extract raw (single-cell) count data for glmmSeq
# # ################################################################################
# # # Seurat typically stores (features x cells):
# # rna_counts <- GetAssayData(hyp, assay = "RNA", slot = "counts")
# # # Make it a regular matrix (or it can stay sparse, but glmmSeq often expects dense):
# # countdata <- as.matrix(rna_counts)
# # 
# # # Confirm dimensions:
# # #  - rows = genes
# # #  - cols = cells
# # dim(countdata)
# # 
# # ################################################################################
# # # 4) Prepare cell-level metadata
# # ################################################################################
# # # The Seurat meta.data has rows = cells, so we match that order to countdata columns:
# # metadata <- hyp@meta.data
# # 
# # # Reorder metadata rows so they align exactly to colnames(countdata):
# # # (Often this is already the same, but good to be explicit.)
# # metadata <- metadata[colnames(countdata), ]
# # 
# # # Example: define a new factor "Mbuna" vs "NonMbuna"
# # # (If you already have a Species column, you can recode it, or keep it as is.)
# # metadata$Mbuna <- ifelse(metadata$Species %in% c("MbunaSpecies1", "MbunaSpecies2"), 
# #                          "Mbuna", 
# #                          "NonMbuna")
# # 
# # metadata$cellID <- colnames(countdata)  # Each cell gets a unique ID
# # ################################################################################
# # # 5) (Optional) Compute library-size offsets
# # ################################################################################
# # # glmmSeq can automatically handle offsets via 'offset = log(...)'.
# # # This is helpful if library sizes differ greatly among cells:
# # library_sizes <- colSums(countdata)
# # offset_vec    <- log(library_sizes)
# # 
# # ################################################################################
# # # 6) Fit glmmSeq model
# # ################################################################################
# # # In glmmSeq, we specify:
# # #   - formula   : e.g. count ~ Mbuna (fixed effect)
# # #   - countdata : matrix with genes in rows, cells in columns
# # #   - metadata  : one row per column of 'countdata'
# # #   - id        : the grouping factor for random intercepts
# # #                 If no repeated measures or no grouping, set id=1
# # #   - offset    : optional, log of library size (if desired)
# # #   - method    : "LRT", "wald", etc. for significance test
# # #   - parallel  : whether to run in parallel
# # #   - cores     : number of CPU cores to use
# # #
# # # Below, we treat "Mbuna" as our main effect of interest, 
# # # and assume no random grouping factor (id=1).
# # #
# # # NOTE: For a large single-cell dataset, glmmSeq can be quite slow. 
# # #       Adjust 'min.count', 'min.cells', or subset genes to speed up.
# # 
# # res_glmm <- glmmSeq(
# #     modelFormula   = count ~ Mbuna + (1 | cellID),  # Adding a random effect
# #     countdata      = countdata,       # genes x cells
# #     metadata       = metadata,        # matching order of columns(cells)
# #     method         = "lme4",          # fitting engine ("lme4" or "glmmTMB")
# #     test.statistic = "LRT",           # "LRT" or "Wald"
# #     offset         = offset_vec,      
# #     min.count      = 1,              
# #     min.cells      = 1,              
# #     remove.0.cells = TRUE,           
# #     parallel       = TRUE,            
# #     cores          = 4                
# # )
# # 
# # colnames(metadata)
# # variables <- as.character(variables)
# # 
# # 
# # 
# # ################################################################################
# # # 7) Inspect DE results
# # ################################################################################
# # # The main results table is typically found in the 'results' slot:
# # deg_results <- res_glmm$results
# # 
# # # 'deg_results' contains estimates, p-values, etc. sorted by p-value
# # head(deg_results)
# # 
# # # For example, you might check:
# # #   deg_results$logFC   (log fold-change)
# # #   deg_results$Pvalue  (raw p-value)
# # #   deg_results$FDR     (adjusted p-value)
# # 
# # ################################################################################
# # # Done!
# # ################################################################################
# 
# ################################################################################
# # 1) Load libraries
# ################################################################################
# 
# library(Seurat)
# library(glmmSeq)
# 
# ################################################################################
# # 2) Load Seurat object
# ################################################################################
# 
# file_path <- "/home/core/hyp/hyp.rds"
# hyp <- readRDS(file_path)
# 
# ################################################################################
# # 3) Extract raw (single-cell) count data
# ################################################################################
# 
# # Seurat stores counts as (features x cells):
# rna_counts <- GetAssayData(hyp, assay = "RNA", slot = "counts")
# 
# # Convert to a normal matrix (if sparse):
# countdata <- as.matrix(rna_counts)
# 
# ################################################################################
# # 4) Prepare metadata (must match countdata)
# ################################################################################
# 
# # Get metadata from Seurat
# metadata <- hyp@meta.data
# 
# # Ensure rownames(metadata) match colnames(countdata)
# metadata <- metadata[colnames(countdata), ]
# 
# # Define Mbuna vs NonMbuna groups
# metadata$Mbuna <- ifelse(metadata$Species %in% c("MbunaSpecies1", "MbunaSpecies2"), 
#                          "Mbuna", 
#                          "NonMbuna")
# 
# # Add a unique 'cellID' as a random effect (required by glmmSeq)
# metadata$cellID <- colnames(countdata)  # Each cell gets a unique ID
# 
# ################################################################################
# # 5) Define Variables and Subset Metadata
# ################################################################################
# 
# # Specify the columns we need
# variables <- c("Species", "Mbuna", "cellID")
# 
# # (Optional) Check that all required columns are present
# missing_vars <- setdiff(variables, colnames(metadata))
# if (length(missing_vars) > 0) {
#     warning("⚠️ The following columns are missing from metadata: ", 
#             paste(missing_vars, collapse = ", "))
#     variables <- intersect(variables, colnames(metadata))
# }
# 
# # Convert metadata to a pure data frame, then subset
# metadata <- as.data.frame(metadata)
# metadata_subset <- metadata[, variables, drop = FALSE]
# 
# ################################################################################
# # 6) Library-size offsets (optional but improves model accuracy)
# ################################################################################
# 
# library_sizes <- colSums(countdata)
# offset_vec    <- log(library_sizes)
# 
# ################################################################################
# # 7) Run glmmSeq with the Corrected Model Formula
# ################################################################################
# 
# # Note: We now use a model formula without the response variable.
# res_glmm <- glmmSeq(
#     modelFormula   = ~ Mbuna + (1 | cellID),  # Only covariates & random effect
#     countdata      = countdata,               # genes x cells
#     metadata       = metadata_subset,         # matching order of columns (cells)
#     method         = "lme4",                  # fitting engine ("lme4" or "glmmTMB")
#     test.statistic = "LRT",                   # "LRT" or "Wald"
#     offset         = offset_vec,      
#     min.count      = 1,              
#     min.cells      = 1,              
#     remove.0.cells = TRUE,           
#     parallel       = TRUE,            
#     cores          = 4                
# )
# 
# table(metadata$Species)
# table(metadata$cellID)  # Less likely, but sometimes cellID can collapse.
# 
# table(metadata$Mbuna)
# 
# # 1) Confirm which species you have in your metadata
# unique(metadata$Species)
# 
# # 2) Confirm how many cells per Mbuna level
# table(metadata$Mbuna)
# 
# all(rownames(metadata) == colnames(countdata))  # Should return TRUE
# 

# deg_results <- res_glmm$results
# head(deg_results)

library(Seurat)
library(glmmSeq)

file_path <- "/home/core/hyp/hyp.rds"
hyp <- readRDS(file_path)

################################################################################
# 3) Extract raw (single-cell) count data
################################################################################

rna_counts <- GetAssayData(hyp, assay = "RNA", slot = "counts")
countdata  <- as.matrix(rna_counts)

################################################################################
# 4) Prepare metadata
################################################################################

metadata <- hyp@meta.data
metadata <- metadata[colnames(countdata), ]

# "LF" and "MC" are Mbuna species maybe?
mbuna_species <- c("LF", "MC")
metadata$Mbuna <- ifelse(
    metadata$Species %in% mbuna_species,
    "Mbuna",
    "NonMbuna"
)
metadata$Mbuna <- factor(metadata$Mbuna, levels = c("Mbuna", "NonMbuna"))

# Random effect: one cell ID per cell
metadata$cellID <- colnames(countdata)

################################################################################
# Filtering out low-count genes
################################################################################

# Keeping genes with >=10 total counts
keep_genes <- rowSums(countdata) >= 10
countdata_filtered <- countdata[keep_genes, ]

################################################################################
# 6) Subseting metadata columns
################################################################################

variables <- c("Species", "Mbuna", "cellID")
missing_vars <- setdiff(variables, colnames(metadata))
if (length(missing_vars) > 0) {
    warning("Missing from metadata: ", paste(missing_vars, collapse = ", "))
    variables <- intersect(variables, colnames(metadata))
}

metadata <- as.data.frame(metadata)
metadata_subset <- metadata[, variables, drop = FALSE]

################################################################################
# 7) Library-size offsets
################################################################################

library_sizes <- colSums(countdata_filtered)
offset_vec    <- log(library_sizes)

################################################################################
# 8) Run glmmSeq
################################################################################

res_glmm <- glmmSeq(
    modelFormula   = ~ Mbuna + (1 | cellID),
    countdata      = countdata_filtered, 
    metadata       = metadata_subset,
    method         = "glmmTMB", # either lme4 or glmmTMB
    test.statistic = "LRT",
    offset         = offset_vec,      
    min.count      = 5,       # Raise the threshold from 1 to 5
    min.cells      = 5,       # Raise the threshold from 1 to 5
    remove.0.cells = TRUE,
    parallel       = TRUE,
    cores          = 4
)

deg_results <- res_glmm$results
head(deg_results)

