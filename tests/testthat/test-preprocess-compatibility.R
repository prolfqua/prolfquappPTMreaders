test_that("FragPipe multisite preprocessing accepts AnalysisConfiguration", {
  path <- system.file(
    "extdata",
    "FP_multisite",
    package = "prolfquappPTMreaders"
  )
  files <- get_FP_multi_site_files(path)
  annot <- dataset_template_FP_multi_site(files)
  annot$Group <- ifelse(grepl("^WT", annot$Name), "WT", "KO")
  annot$Subject <- sub(".*_", "", annot$Name)
  annot$CONTROL <- ifelse(annot$Group == "WT", "C", "T")
  annot$isotopeLabel <- "light"
  annotation <- prolfquapp::read_annotation(annot)

  expect_s3_class(annotation$atable, "AnalysisConfiguration")

  result <- preprocess_FP_multi_site(
    files$data,
    files$fasta,
    annotation,
    sitetype = "multisite"
  )

  expect_s3_class(result$lfqdata, "LFQData")
  expect_gt(nrow(result$lfqdata$data_long()), 0)
})

test_that("FragPipe combined STY preprocessing accepts AnalysisConfiguration", {
  path <- system.file(
    "extdata",
    "FP_combined_STY",
    package = "prolfquappPTMreaders"
  )
  files <- get_FP_combined_STY_files(path)
  annot <- dataset_template_FP_combined_STY(files)
  annot$Group <- ifelse(grepl("37C", annot$Name), "37C", "42C")
  annot$Control <- ifelse(annot$Group == "37C", "C", "T")
  annot$isotopeLabel <- "light"
  annotation <- prolfquapp::read_annotation(annot)

  expect_s3_class(annotation$atable, "AnalysisConfiguration")

  expect_snapshot(
    {
      invisible(capture.output(
        result <- preprocess_FP_combined_STY(
          files$data,
          files$fasta,
          annotation
        )
      ))
    },
    cran = TRUE
  )

  expect_s3_class(result$lfqdata, "LFQData")
  expect_gt(nrow(result$lfqdata$data_long()), 0)
})

test_that("Spectronaut site preprocessing accepts AnalysisConfiguration", {
  path <- system.file(
    "extdata",
    "BGS_site",
    package = "prolfquappPTMreaders"
  )
  files <- get_BGS_site_files(path)
  annot <- dataset_template_BGS_site(files)
  annot$Group <- "A"
  annot$Control <- "C"
  annot$isotopeLabel <- "light"
  annotation <- prolfquapp::read_annotation(annot)

  expect_s3_class(annotation$atable, "AnalysisConfiguration")

  result <- preprocess_BGS_site(
    files$data,
    files$fasta,
    annotation
  )

  expect_s3_class(result$lfqdata, "LFQData")
  expect_gt(nrow(result$lfqdata$data_long()), 0)
})
