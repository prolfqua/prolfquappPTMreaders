fixture <- function() {
  path <- system.file("extdata", "MSSTATS_site", package = "prolfquappPTMreaders")
  files <- get_MSSTATS_site_files(path)
  annot <- dataset_template_MSSTATS_site(files)
  annot$Group <- rep(c("A", "B"), length.out = nrow(annot))
  annot$Control <- ifelse(annot$Group == "A", "C", "T")
  list(files = files, annotation = prolfquapp::read_annotation(annot))
}

test_that("the reader is registered and callable the way prolfquapp calls it", {
  expect_true("MSSTATS_site" %in% names(prolfqua_preprocess_functions))
  # prolfquapp builds the call from these five by name, and adds nr_peptides
  # only for readers that declare it (.forward_nr_peptides). This reader does
  # not: a site-level table has one quantity per site, so a minimum-peptides
  # filter has nothing to act on.
  expect_identical(
    names(formals(preprocess_MSSTATS_site)),
    c("quant_data", "fasta_file", "annotation", "pattern_contaminants", "pattern_decoys")
  )
})

test_that("get_MSSTATS_site_files finds the table and the fasta", {
  files <- fixture()$files
  expect_match(basename(files$data), "\\.csv$")
  expect_match(basename(files$fasta), "\\.fasta$")
})

test_that("PeptideSequence defaults to the index", {
  sites <- read_MSSTATS_site(fixture()$files$data)
  expect_identical(sites$PeptideSequence, sites$Index)
})

test_that("a table without the required columns fails, naming them", {
  incomplete <- tempfile(fileext = ".csv")
  readr::write_csv(data.frame(ProteinName = "PA", Run = "run1"), incomplete)
  expect_error(read_MSSTATS_site(incomplete), "Index")
  expect_error(read_MSSTATS_site(incomplete), "Intensity")
})

test_that("preprocessing yields site-level lfqdata and site annotation", {
  f <- fixture()
  result <- preprocess_MSSTATS_site(f$files$data, f$files$fasta, f$annotation)

  expect_named(result, c("lfqdata", "protein_annotation"))
  expect_equal(result$lfqdata$relevant_hierarchy_keys(), c("protein_Id", "site"))
  # four sites in four runs
  expect_equal(nrow(result$lfqdata$data_long()), 16)

  row_annot <- result$protein_annotation$row_annot
  expect_true(all(c("modAA", "posInProtein", "SequenceWindow") %in% names(row_annot)))
  expect_equal(nrow(row_annot), 4)
})

test_that("the site annotation is read out of the index and cut from the fasta", {
  f <- fixture()
  result <- preprocess_MSSTATS_site(f$files$data, f$files$fasta, f$annotation)
  row_annot <- result$protein_annotation$row_annot
  row_annot <- row_annot[order(row_annot$Index), ]

  expect_equal(row_annot$Index, c("PA_S10", "PA_S31", "PB_T10", "PB_Y31"))
  expect_equal(row_annot$modAA, c("S", "S", "T", "Y"))
  expect_equal(row_annot$posInProtein, c(10L, 31L, 10L, 31L))
  # The window is the fasta's own residues, so its centre is the modified one.
  expect_true(all(nchar(row_annot$SequenceWindow) == 15))
  expect_equal(substr(row_annot$SequenceWindow, 8, 8), row_annot$modAA)
})

test_that("an annotation naming runs the table does not have fails loudly", {
  f <- fixture()
  annot <- f$annotation
  annot$annot[[annot$atable$file_name]] <- paste0("absent_", seq_len(nrow(annot$annot)))
  expect_error(
    preprocess_MSSTATS_site(f$files$data, f$files$fasta, annot),
    "no run in the site table matches"
  )
})
