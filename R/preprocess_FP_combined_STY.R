


#' get report.tsv and fasta file location in folder
#' @return list with paths to data and fasta
#' @export
get_FP_combined_STY_files <- function(path){
  psm_file <- dir(path = path, pattern = "^combined_site_STY_.+\\.tsv", recursive = TRUE, full.names = TRUE)
  fasta.files <- grep("*.fasta$|*.fas$", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  fp.manifest <- grep("*.fp-manifest",  dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = psm_file, fasta = fasta.files, fp.manifest = fp.manifest))
}


#' reads combined_site_STY file and converts to long format.
#' @return tidy data table
#' @export
#'
read_combined_STY_file <- function(file){
  xd <- readr::read_tsv(file, show_col_types = FALSE)
  colnames(xd) <- gsub("Localization Probability", "Localization_Probability", colnames(xd))
  colnames(xd) <- gsub("MaxLFQ Intensity", "MaxLFQ_Intensity", colnames(xd))
  xd <- xd |> dplyr::rename(BLP = "Best Localization_Probability")
  tidy_data <- xd |>
    tidyr::pivot_longer(
      cols = tidyselect::contains("Localization_Probability") | tidyselect::contains("Intensity") | tidyselect::contains("MaxLFQ_Intensity"),
      names_to = c("SampleName", ".value"),
      names_sep = " "
    )
  tidy_data <- tidy_data |> dplyr::rename(ProteinID = !!rlang::sym("Protein ID"))
  return(tidy_data)
}


#' create dataset tempalte for FP combined STY file
#' @return data.frame
#' @export
#'
dataset_template_FP_combined_STY_v2 <- function(files){
  res_data <- prolfquappPTMreaders::read_combined_STY_file(files$data)

  manifest <- readr::read_tsv(files$fp.manifest, col_names = FALSE)
  colnames(manifest) <- c("raw.file", "Name", "Experiment", "Data_type")
  manifest$Experiment <- NULL
  manifest$Data_type <- NULL

  res_data <- dplyr::inner_join(manifest, res_data, by = c(Name = "SampleName"))

  datasetannot <- res_data |>
    dplyr::select(tidyselect::all_of(c("raw.file", "Name"))) |>
    dplyr::distinct()
  datasetannot$Group <- ""
  datasetannot$Subject <- ""
  datasetannot$Control <- ""
  return(datasetannot)
}

#' create dataset tempalte for FP combined STY file
#' @return data.frame
#' @export
#'

dataset_template_FP_combined_STY <- function(files){
  res_data <- prolfquappPTMreaders::read_combined_STY_file(files$data)

  # manifest <- readr::read_tsv(files$fp.manifest, col_names = FALSE)
  # colnames(manifest) <- c("raw.file", "Name", "Experiment", "Data_type")
  # manifest$Experiment <- NULL
  # manifest$Data_type <- NULL

  # res_data <- dplyr::inner_join(manifest, res_data, by = c(Name = "SampleName"))

  datasetannot <- res_data |>
    dplyr::select(file = "SampleName", Name = "SampleName") |>
    dplyr::distinct()
  datasetannot$Group <- ""
  datasetannot$Subject <- ""
  datasetannot$Control <- ""
  return(datasetannot)
}

#' preprocess FP multisite, filter by purity_threshold and PeptideProphetProb
#' @return list with lfqdata and protein annotation
#' @export
#' @param annotation_join_by column in annotation file
preprocess_FP_combined_STY <- function(
    quant_data,
    fasta_file,
    annotation,
    pattern_contaminants = "^zz|^CON|Cont_",
    pattern_decoys = "^REV_|^rev_",
    annotation_join_by = "SampleName"
){

  #annotation_join_by <- match.arg(annotation_join_by)
  # pattern_contaminants = "^zz|^CON"

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(normalize_path(annot[[atable$fileName]])))
    ))

  multiSite_long <- prolfquappPTMreaders::read_combined_STY_file(quant_data)
  # join with anno again this should work now with Name # if not all samples are used in the dataset they would be removed here (to be tested)
  by = annotation_join_by
  names(by) = annotation$atable$fileName
  multiSite_long <- dplyr::inner_join(x = annotation$annot, y = multiSite_long, by = by)

  # add missing required parameters (qvalue)
  multiSite_long$qValue <- 0
  multiSite_long$nr_children  <- 1


  # Setup configuration manually for peptide analysis (phospho)
  atable$ident_Score = "Localization_Probability"
  atable$ident_qValue = "qValue"
  atable$nr_children = "nr_children"
  atable$hierarchy[["protein_Id"]] <- c("ProteinID")
  atable$hierarchy[["site"]] <- c("Index", "Peptide")
  atable$set_response("Intensity")
  atable$hierarchyDepth <- 2

  # Preprocess data - aggregate proteins.
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(multiSite_long, config)

  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities(threshold = 1)


  # Create fasta annotation
  # Create Site Annotation
  site_annot <- multiSite_long |>
    dplyr::select(c("Index", "Protein", "ProteinID", "Peptide", "BLP")) |>
    dplyr::distinct()

  phosSite <- site_annot |> dplyr::rowwise() |> dplyr::mutate(siteinfo = gsub(ProteinID, "", Index))
  phosSite <- phosSite |> dplyr::rowwise() |> dplyr::mutate(PhosSites = gsub("^_", "", siteinfo))
  phosSite <- phosSite |>
    tidyr::extract("PhosSites", into = c("modAA", "posInProtein"), regex = "([A-Z])(\\d+)", convert = TRUE, remove = FALSE)
  nrPep_exp <- multiSite_long |>
    dplyr::select(ProteinID, Peptide) |>
    dplyr::distinct() |>
    dplyr::group_by(ProteinID) |>
    dplyr::summarize(nrPeptides = dplyr::n()) |> dplyr::ungroup()

  fasta_annot <- prolfquapp::get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)

  fasta_annot <- dplyr::left_join(nrPep_exp, fasta_annot, by = c(ProteinID = "proteinname"), multiple = "all")
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)
  fasta_annot2 <- dplyr::inner_join(fasta_annot, phosSite, by = "ProteinID")

  # Make names to match lfqdata
  fasta_annot2 <- fasta_annot2 |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!rlang::sym("ProteinID"))
  fasta_annot2 <- fasta_annot2 |> dplyr::mutate(!!lfqdata$config$table$hierarchy_keys_depth()[2] := paste(!!rlang::sym("Index"),!!rlang::sym("Peptide"), sep = "~"))
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata ,
    fasta_annot2,
    description = "description",
    cleaned_ids = "protein_Id",
    full_id = "protein_Id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )

  # Verify lfqdata and protein_annotation match on hierarchy keys
  stopifnot(nrow(dplyr::inner_join(prot_annot$row_annot, lfqdata$data, by = lfqdata$config$table$hierarchy_keys_depth())) > 0)
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}

