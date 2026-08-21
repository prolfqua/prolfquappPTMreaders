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
