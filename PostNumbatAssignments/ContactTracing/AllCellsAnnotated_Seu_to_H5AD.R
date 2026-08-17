required_packages <- c("Seurat", "SeuratDisk")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    sprintf(
      "Missing required package(s): %s",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

# Using the:
# 1. the full integrated Seurat object
# 2. the Numbat-derived non-malignant Seurat object
#
# Inputs are generated in:
# - PreProcessing/Integration.R.ipynb
# - PostNumbatAssignments/NonTumorObjAnalysis-Copy1.ipynb
#
# The broad non-malignant labels in BI_broad_assignments_v2 were assigned after
# manual review of each cluster together with the DEGs/marker genes enriched in
# those clusters.

# Required packages: Seurat, SeuratDisk
# Update these paths if running from a different working directory.
all_cells_rds <- "data/sarcoma_all/data_sarcoma_all_merged_obj.rds"
non_malignant_rds <- "data/sarcoma_all/data_numbat_sarcoma_all_merged_obj_non_malignant.rds"
output_h5seurat <- "PostNumbatAssignments/RefinedCellTyping/all_cells_merged_bi_assignments.h5Seurat"
output_h5ad <- "PostNumbatAssignments/RefinedCellTyping/all_cells_merged_bi_assignments.h5ad"

bi_broad_assignments_v2_levels <- function() {
  c(
    "Myeloid",
    "T-cell",
    "Fibroblast",
    "B-cell/Plasma cell",
    "Endothelial cell",
    "Myofibroblast",
    "Malignant"
  )
}

assign_bi_broad_v2 <- function(
  seu,
  cluster_field = "integrated_snn_res.0.2",
  assignment_field = "BI_broad_assignments_v2"
) {
  if (assignment_field %in% colnames(seu@meta.data)) {
    return(seu)
  }

  if (!cluster_field %in% colnames(seu@meta.data)) {
    stop(
      sprintf(
        "Missing '%s' in non-malignant metadata; cannot reconstruct '%s'.",
        cluster_field,
        assignment_field
      )
    )
  }

  # Cluster-to-label mapping from the refined non-malignant cell typing step.
  broad_assignments_v2 <- c(
    "Myeloid",
    "T-cell",
    "Fibroblast",
    "B-cell/Plasma cell",
    "T-cell",
    "Myeloid",
    "Endothelial cell",
    "T-cell",
    "T-cell",
    "B-cell/Plasma cell",
    "Myeloid",
    "B-cell/Plasma cell",
    "Myofibroblast",
    "Myeloid"
  )

  Idents(object = seu) <- cluster_field
  cluster_levels <- levels(Idents(seu))

  if (length(cluster_levels) != length(broad_assignments_v2)) {
    stop(
      sprintf(
        "Expected %d cluster levels in '%s', found %d.",
        length(broad_assignments_v2),
        cluster_field,
        length(cluster_levels)
      )
    )
  }

  seu <- RenameIdents(
    object = seu,
    new.names = setNames(broad_assignments_v2, cluster_levels)
  )
  seu[[assignment_field]] <- Idents(object = seu)
  Idents(object = seu) <- cluster_field
  seu
}

merge_bi_assignments_into_all_cells <- function(
  all_cells,
  non_malignant,
  assignment_field = "BI_broad_assignments_v2"
) {
  # Ensure the refined labels exist on the non-malignant object.
  non_malignant <- assign_bi_broad_v2(
    seu = non_malignant,
    assignment_field = assignment_field
  )

  # Merge the non-malignant labels onto the full object by barcode.
  non_malignant_meta <- non_malignant@meta.data[, assignment_field, drop = FALSE]
  non_malignant_meta$barcode <- rownames(non_malignant_meta)

  all_meta <- all_cells@meta.data
  all_meta$barcode <- rownames(all_meta)

  merged_df <- merge(
    x = all_meta,
    y = non_malignant_meta,
    all.x = TRUE,
    by = "barcode",
    sort = FALSE
  )

  rownames(merged_df) <- merged_df$barcode
  merged_df <- merged_df[rownames(all_cells@meta.data), , drop = FALSE]

  if (!identical(rownames(merged_df), rownames(all_cells@meta.data))) {
    stop("Merged metadata row order no longer matches the full all-cells Seurat object.")
  }

  # Cells not present in the non-malignant object are labeled Malignant.
  merged_df[[assignment_field]] <- ifelse(
    is.na(merged_df[[assignment_field]]),
    "Malignant",
    as.character(merged_df[[assignment_field]])
  )

  all_cells@meta.data[[assignment_field]] <- factor(
    merged_df[[assignment_field]],
    levels = bi_broad_assignments_v2_levels()
  )

  all_cells
}

reconstruct_all_cells_merged_bi_assignments <- function(
  all_cells_rds,
  non_malignant_rds,
  output_h5seurat,
  output_h5ad = sub("\\.h5Seurat$", ".h5ad", output_h5seurat),
  assignment_field = "BI_broad_assignments_v2"
) {
  dir.create(dirname(output_h5seurat), recursive = TRUE, showWarnings = FALSE)

  # Load the full integrated object and the Numbat-derived non-malignant object.
  all_cells <- readRDS(all_cells_rds)
  non_malignant <- readRDS(non_malignant_rds)

  # Transfer refined broad labels onto the full object.
  all_cells <- merge_bi_assignments_into_all_cells(
    all_cells = all_cells,
    non_malignant = non_malignant,
    assignment_field = assignment_field
  )

  # Export to h5Seurat, then convert to h5ad.
  SaveH5Seurat(all_cells, filename = output_h5seurat, overwrite = TRUE)
  Convert(output_h5seurat, dest = "h5ad", overwrite = TRUE)

  default_h5ad <- sub("\\.h5Seurat$", ".h5ad", output_h5seurat)
  if (!identical(normalizePath(default_h5ad, mustWork = FALSE), normalizePath(output_h5ad, mustWork = FALSE))) {
    if (file.exists(output_h5ad)) {
      file.remove(output_h5ad)
    }
    ok <- file.rename(default_h5ad, output_h5ad)
    if (!ok) {
      stop(sprintf("Failed to move '%s' to '%s'.", default_h5ad, output_h5ad))
    }
  }

  invisible(all_cells)
}

# Run the reconstruction when the script is sourced or executed.
reconstruct_all_cells_merged_bi_assignments(
  all_cells_rds = all_cells_rds,
  non_malignant_rds = non_malignant_rds,
  output_h5seurat = output_h5seurat,
  output_h5ad = output_h5ad
)
