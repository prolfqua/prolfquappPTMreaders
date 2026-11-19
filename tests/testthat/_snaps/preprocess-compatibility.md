# FragPipe combined STY preprocessing accepts AnalysisConfiguration

    Code
      invisible(capture.output(result <- preprocess_FP_combined_STY(files$data, files$
        fasta, annotation)))
    Message
      remove_small_intensities: 924 values set to NA (e.g. A0A2R8Y4L2~A0A2R8Y4L2_S309~NQGGYGGSsSSSSYGSGR, A0A2R8Y4L2~A0A2R8Y4L2_S310~NQGGYGGSSsSSSYGSGR, A0A2R8Y4L2~A0A2R8Y4L2_S311~NQGGYGGSSSsSSYGSGR, A0A2R8Y4L2~A0A2R8Y4L2_S316~NQGGYGGSSSSSSYGsGR, A0A2R8Y4L2~A0A2R8Y4L2_Y314~NQGGYGGSSSSSSyGSGR, ...): intensity below 1
      See $provenance().

