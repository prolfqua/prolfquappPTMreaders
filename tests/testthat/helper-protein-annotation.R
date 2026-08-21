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

  # The annotation is one row per key, where the key is protein_Id for a
  # protein-level reader and protein_Id + site for a site-level one. Either way
  # the join below must not multiply the quant rows.
  expect_identical(anyDuplicated(row_annot[, protein_id, drop = FALSE]), 0L)
  expect_true(all(protein_id %in% result$lfqdata$relevant_hierarchy_keys()))
  expect_identical(protein_id[[1]], result$lfqdata$relevant_hierarchy_keys()[[1]])
  expect_equal(nrow(annotated_quant), nrow(quant_data))
  expect_identical(sum(is.na(annotated_quant$description)), 0L)
  expect_identical(sum(is.na(annotated_quant$gene_name)), 0L)
  expect_identical(sum(is.na(annotated_quant$CON)), 0L)
}
