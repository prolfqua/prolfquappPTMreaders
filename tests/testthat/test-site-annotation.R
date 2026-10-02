test_that("get_sequence_windows centres the window on the site", {
  d <- data.frame(
    sequence = "MKFLVLLFNILCLFPVLAADNH",
    posInProtein = c(5L, 12L)
  )

  res <- get_sequence_windows(d, pos_in_protein = "posInProtein")

  expect_equal(res$sequence_window, c("XXXMKFLVLLFNILC", "VLLFNILCLFPVLAA"))
})


test_that("get_sequence_windows leaves the window missing when the sequence is", {
  # A site whose protein is not in the FASTA reaches this with sequence NA.
  # paste0() used to render that as the literal "NA", giving a window of X's
  # around two letters that reads like a real one.
  d <- data.frame(
    sequence = c("MKFLVLLFNILCLFPVLAADNH", NA_character_, NA_character_),
    posInProtein = c(5L, 5L, 108L)
  )

  res <- get_sequence_windows(d, pos_in_protein = "posInProtein")

  expect_equal(res$sequence_window[1], "XXXMKFLVLLFNILC")
  expect_true(all(is.na(res$sequence_window[2:3])))
})


test_that("get_sequence_windows leaves the window missing when the position is", {
  d <- data.frame(sequence = "MKFLVLLFNILCLFPVLAADNH", posInProtein = NA_integer_)

  res <- get_sequence_windows(d, pos_in_protein = "posInProtein")

  expect_true(is.na(res$sequence_window))
})


test_that("parse_site_index reads the residue and its position", {
  res <- parse_site_index(c("A0A1B0GTU1_S108", "P02545_T3"))

  expect_equal(res$modAA, c("S", "T"))
  expect_equal(res$posInProtein, c(108L, 3L))
})


test_that("parse_site_index returns NA for an index naming no site", {
  # The multi-site file writes the unmodified form of a peptide this way.
  res <- parse_site_index(c("A2AKB9_354_369_1_0", "nonsense"))

  expect_true(all(is.na(res$modAA)))
  expect_true(all(is.na(res$posInProtein)))
})


test_that("the Spectronaut site reader hands over residue, position and window", {
  path <- system.file("extdata", "BGS_site", package = "prolfquappPTMreaders")
  files <- get_BGS_site_files(path)
  annot <- dataset_template_BGS_site(files)
  annot$Group <- "A"
  annot$Control <- "C"
  annotation <- prolfquapp::read_annotation(annot)

  result <- suppressMessages(preprocess_BGS_site(files$data, files$fasta, annotation))
  row_annot <- result$protein_annotation$row_annot
  sites <- tidyr::separate(
    row_annot,
    "site",
    into = c("p", "key", "aa", "location", "multiplicity"),
    sep = "~",
    remove = FALSE
  )

  expect_true(all(c("modAA", "posInProtein", "SequenceWindow") %in% names(row_annot)))
  expect_equal(row_annot$modAA, sites$aa)
  expect_equal(row_annot$posInProtein, as.integer(sites$location))
  expect_true(all(nchar(row_annot$SequenceWindow) == 15L))
  expect_equal(substr(row_annot$SequenceWindow, 8, 8), row_annot$modAA)
})
