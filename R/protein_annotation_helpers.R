.validate_protein_annotation <- function(lfqdata, protein_annotation) {
  protein_id <- protein_annotation$pID
  row_annot <- protein_annotation$row_annot
  quant_data <- lfqdata$data_long()

  stopifnot(length(protein_id) >= 1L)
  stopifnot(anyDuplicated(row_annot[, protein_annotation$pID, drop = FALSE]) == 0L)
  # A site-level annotation legitimately carries the deeper hierarchy keys: that
  # is what makes it one row per site rather than per protein.

  annotated_quant <- dplyr::right_join(
    row_annot,
    quant_data,
    by = protein_id,
    multiple = "all"
  )
  stopifnot(nrow(annotated_quant) == nrow(quant_data))
  invisible(protein_annotation)
}
