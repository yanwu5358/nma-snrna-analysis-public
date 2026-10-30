# Marker-defined motor neurons using the original SCT corrected counts.
identify_motor_neurons <- function(object, celltype_column = "celltype_res05") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Install SeuratObject")
  md <- object[[]]
  if (!celltype_column %in% names(md)) stop("Missing cell-type metadata")
  genes <- c("SLC5A7", "ISL1", "ISL2", "MNX1", "STMN2", "VSX2", "SIM1")
  if (!"SCT" %in% SeuratObject::Assays(object)) stop("Missing SCT assay")
  missing <- setdiff(genes, rownames(object[["SCT"]]))
  if (length(missing)) stop("Missing markers: ", paste(missing, collapse = ", "))
  x <- SeuratObject::GetAssayData(object, assay = "SCT", layer = "counts")[genes, , drop = FALSE]
  keep <- !is.na(md[[celltype_column]]) & md[[celltype_column]] == "Neuron" &
    x["SLC5A7", ] > 0 &
    (x["ISL1", ] > 0 | x["ISL2", ] > 0 | x["MNX1", ] > 0 | x["STMN2", ] > 5) &
    x["VSX2", ] < 1 & x["SIM1", ] < 1
  setNames(as.logical(keep), colnames(object))
}
