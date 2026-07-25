# FragPipe combined STY preprocessing accepts AnalysisConfiguration

    Code
      invisible(capture.output(result <- preprocess_FP_combined_STY(files$data, files$
        fasta, annotation)))
    Condition
      Warning in `dplyr::left_join()`:
      Detected an unexpected many-to-many relationship between `x` and `y`.
      i Row 1 of `x` matches multiple rows in `y`.
      i Row 1 of `y` matches multiple rows in `x`.
      i If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.
    Message
      column sampleName already exists, using :sampleName
      completing cases
      completing cases done
      setup done
      completing cases
    Condition
      Warning in `dplyr::inner_join()`:
      Detected an unexpected many-to-many relationship between `x` and `y`.
      i Row 1 of `x` matches multiple rows in `y`.
      i Row 1 of `y` matches multiple rows in `x`.
      i If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.
