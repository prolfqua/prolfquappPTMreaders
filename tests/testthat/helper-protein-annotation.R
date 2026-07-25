expect_protein_annotation_contract <- function(result) {
  protein_annotation <- result$protein_annotation
  protein_id <- protein_annotation$pID
  row_annot <- protein_annotation$row_annot
  quant_data <- result$lfqdata$data_long()

  protein_annotation$annotate_contaminants()
  annotated_quant <- dplyr::right_join(
    protein_annotation$row_annot,
    quant_data,
    by = protein_id,
    multiple = "all"
  )

  expect_identical(anyDuplicated(row_annot[[protein_id]]), 0L)
  expect_length(
    intersect(
      setdiff(result$lfqdata$relevant_hierarchy_keys(), protein_id),
      colnames(row_annot)
    ),
    0L
  )
  expect_equal(nrow(annotated_quant), nrow(quant_data))
  expect_identical(sum(is.na(annotated_quant$description)), 0L)
  expect_identical(sum(is.na(annotated_quant$gene_name)), 0L)
  expect_identical(sum(is.na(annotated_quant$CON)), 0L)
}
