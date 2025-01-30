# Function implementing reading routines for PTM data

This package extends the functionality of the $prolfquapp$ package to support PTM analysis.


After installing this R package you can pass the funcitons to the prolfquapp package using this section in the `config.yaml` file to parse PTM centric output formats generated by FragPipe.

When starting the analysis from `abundance_multi-site_None.tsv` files.

```
ext_reader:
  extra_args: list(pattern_contaminants = "^zz|^CON|Cont_", pattern_decoys = "^REV_|^rev_")
  preprocess: prolfquappPTMreaders::preprocess_FP_multi_site
  get_files: prolfquappPTMreaders::get_FP_multi_site_files
```

When starting the analysis from `^combined_site_STY_.+\\.tsv`

```
ext_reader:
  extra_args: list(pattern_contaminants = "^zz|^CON|Cont_", pattern_decoys = "^REV_|^rev_",  annotation_join_by = c("raw.file", "Name")")
  preprocess: prolfquappPTMreaders::preprocess_FP_combined_STY
  get_files: prolfquappPTMreaders::get_FP_combined_STY_files
```


# Installation instructions

```
remotes::install_github("prolfqua/prolfquappPTMreaders", dependencies = TRUE)
```
