#' preprocess functions implemented in this package
#' @export
prolfqua_preprocess_functions <- list(
  FP_multisite = list(
    extra_args = "list(sitetype = 'multisite')",
    preprocess = "prolfquappPTMreaders::preprocess_FP_multi_site",
    get_files = "prolfquappPTMreaders::get_FP_multi_site_files",
    dataset = "prolfquappPTMreaders::dataset_template_FP_multi_site"
  ),
  FP_singlesite = list(
    extra_args = "list(sitetype = 'singlesite')",
    preprocess = "prolfquappPTMreaders::preprocess_FP_multi_site",
    get_files = "prolfquappPTMreaders::get_FP_single_site_files",
    dataset = "prolfquappPTMreaders::dataset_template_FP_multi_site"
  ),
  FP_combined_STY = list(
    extra_args = "list(annotation_join_by = 'SampleName')",
    preprocess = "prolfquappPTMreaders::preprocess_FP_combined_STY",
    get_files = "prolfquappPTMreaders::get_FP_combined_STY_files",
    dataset = "prolfquappPTMreaders::dataset_template_FP_combined_STY"
  ),
  BGS_site = list(
    extra_args = "list()",
    preprocess = "prolfquappPTMreaders::preprocess_BGS_site",
    get_files = "prolfquappPTMreaders::get_BGS_site_files",
    dataset = "prolfquappPTMreaders::dataset_template_BGS_site"
  )

)


