.validate_protein_annotation <- function(lfqdata, protein_annotation) {
  protein_id <- protein_annotation$pID
  row_annot <- protein_annotation$row_annot
  quant_data <- lfqdata$data_long()

  stopifnot(length(protein_id) == 1L)
  stopifnot(anyDuplicated(row_annot[[protein_id]]) == 0L)
  stopifnot(
    length(intersect(
      setdiff(lfqdata$relevant_hierarchy_keys(), protein_id),
      colnames(row_annot)
    )) ==
      0L
  )

  annotated_quant <- dplyr::right_join(
    row_annot,
    quant_data,
    by = protein_id,
    multiple = "all"
  )
  stopifnot(nrow(annotated_quant) == nrow(quant_data))
  invisible(protein_annotation)
}
