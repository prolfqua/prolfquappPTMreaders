#' Modified Residue and Position from a FragPipe Site Index
#'
#' A `combined_site_*` index names the protein, the modified residue and its
#' position in that protein, e.g. `A0A1B0GTU1_S108`. The site file carries no
#' other positional column -- unlike the TMT `abundance_single-site` file, which
#' ships `Start` and `SequenceWindow` -- so everything downstream that needs a
#' position or a residue is derived from here.
#'
#' @param index character vector of `Index` values.
#' @return data.frame with `modAA` and `posInProtein`; both `NA` for an index
#'   that does not end in a residue and a position.
#' @export
#' @examples
#' parse_site_index(c("A0A1B0GTU1_S108", "P02545_T3", "nonsense"))
parse_site_index <- function(index) {
  m <- regmatches(index, regexpr("_[A-Za-z][0-9]+$", index))
  hit <- regexpr("_[A-Za-z][0-9]+$", index) > 0
  modAA <- rep(NA_character_, length(index))
  pos <- rep(NA_integer_, length(index))
  modAA[hit] <- substr(m, 2, 2)
  pos[hit] <- as.integer(substr(m, 3, nchar(m)))
  data.frame(modAA = modAA, posInProtein = pos, stringsAsFactors = FALSE)
}

#' get report.tsv and fasta file location in folder
#' @param path directory path to search for files
#' @return list with paths to data and fasta
#' @export
get_FP_combined_STY_files <- function(path) {
  psm_file <- dir(
    path = path,
    pattern = "^combined_site_STY_.+\\.tsv",
    recursive = TRUE,
    full.names = TRUE
  )
  fasta.files <- grep(
    "*.fasta$|*.fas$",
    dir(path = path, recursive = TRUE, full.names = TRUE),
    value = TRUE
  )
  fp.manifest <- grep(
    "*.fp-manifest",
    dir(path = path, recursive = TRUE, full.names = TRUE),
    value = TRUE
  )
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
#' @param file path to combined_site_STY file
#' @return tidy data table
#' @export
#'
read_combined_STY_file <- function(file) {
  xd <- readr::read_tsv(file, show_col_types = FALSE)
  colnames(xd) <- gsub(
    "Localization Probability",
    "Localization_Probability",
    colnames(xd)
  )
  colnames(xd) <- gsub("MaxLFQ Intensity", "MaxLFQ_Intensity", colnames(xd))
  xd <- xd |> dplyr::rename(BLP = "Best Localization_Probability")
  tidy_data <- xd |>
    tidyr::pivot_longer(
      cols = tidyselect::contains("Localization_Probability") |
        tidyselect::contains("Intensity") |
        tidyselect::contains("MaxLFQ_Intensity"),
      names_to = c("SampleName", ".value"),
      names_sep = " "
    )
  tidy_data <- tidy_data |>
    dplyr::rename(ProteinID = !!rlang::sym("Protein ID"))
  tidy_data <- dplyr::bind_cols(
    tidy_data,
    prolfquappPTMreaders::parse_site_index(tidy_data$Index)
  )
  return(tidy_data)
}


#' create dataset template for FP combined STY file (v2 with manifest)
#' @param files list with data, fasta and fp.manifest paths
#' @return data.frame with annotation template
#' @export
#'
dataset_template_FP_combined_STY_v2 <- function(files) {
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

#' create dataset template for FP combined STY file
#' @param files list with data and fasta paths
#' @return data.frame with annotation template
#' @export
#'
dataset_template_FP_combined_STY <- function(files) {
  res_data <- prolfquappPTMreaders::read_combined_STY_file(files$data)

  datasetannot <- res_data |>
    dplyr::select(file = "SampleName", Name = "SampleName") |>
    dplyr::distinct()
  datasetannot$Group <- ""
  datasetannot$Subject <- ""
  datasetannot$Control <- ""
  return(datasetannot)
}

#' preprocess FP combined STY file
#' @param quant_data path to combined_site_STY file
#' @param fasta_file path to FASTA file
#' @param annotation annotation object from prolfquapp::read_annotation()
#' @param pattern_contaminants regex pattern to identify contaminant proteins
#' @param pattern_decoys regex pattern to identify decoy proteins
#' @param annotation_join_by column in annotation file to join on
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#' # Use package test data
#' path <- system.file("extdata", "FP_combined_STY", package = "prolfquappPTMreaders")
#' files <- get_FP_combined_STY_files(path)
#' annot_template <- dataset_template_FP_combined_STY(files)
#' annot_template$Group <- ifelse(grepl("37C", annot_template$Name), "37C", "42C")
#' annot_template$Control <- ifelse(annot_template$Group == "37C", "C", "T")
#' annotation <- prolfquapp::read_annotation(annot_template)
#' result <- preprocess_FP_combined_STY(files$data, files$fasta, annotation)
#' stopifnot(nrow(result$lfqdata$data) > 0)
#' stopifnot(nrow(result$protein_annotation$row_annot) > 0)
#' head(result$protein_annotation$row_annot$SequenceWindow)

preprocess_FP_combined_STY <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  annotation_join_by = "SampleName"
) {
  annot <- annotation$annot
  config <- annotation$atable$clone(deep = TRUE)
  annot <- annot |>
    dplyr::mutate(
      raw.file = gsub(
        "^x|.d.zip$|.raw$",
        "",
        (basename(prolfquapp::normalize_path(annot[[config$file_name]])))
      )
    )

  multiSite_long <- prolfquappPTMreaders::read_combined_STY_file(quant_data)
  # join with anno again; this should work now with Name. If not all samples are used in the
  # dataset they would be removed here (to be tested)
  by <- annotation_join_by
  names(by) <- annotation$atable$file_name
  multiSite_long <- dplyr::inner_join(
    x = annotation$annot,
    y = multiSite_long,
    by = by
  )

  fasta_annot_early <- prolfquapp::get_annot_from_fasta(
    fasta_file,
    pattern_decoys = pattern_decoys
  )
  multiSite_long$Protein <- multiSite_long$ProteinID

  # add missing required parameters (qvalue)
  multiSite_long$qValue <- 0
  multiSite_long$nr_children <- 1

  # Setup configuration manually for peptide analysis (phospho)
  config$ident_score <- "Localization_Probability"
  config$ident_q_value <- "qValue"
  config$nr_children <- "nr_children"
  config$hierarchy[["protein_Id"]] <- c("Protein")
  config$hierarchy[["site"]] <- c("Index", "Peptide")
  config$set_response("Intensity")
  config$hierarchy_depth <- 2

  # Preprocess data - aggregate proteins.
  adata <- prolfqua::setup_analysis(multiSite_long, config)

  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities(threshold = 1)

  nrPep_exp <- multiSite_long |>
    dplyr::select(Protein, Peptide) |>
    dplyr::distinct() |>
    dplyr::group_by(Protein) |>
    dplyr::summarize(nrPeptides = dplyr::n()) |>
    dplyr::ungroup()

  fasta_annot <- dplyr::left_join(
    nrPep_exp,
    fasta_annot_early,
    by = c(Protein = "proteinname"),
    multiple = "all"
  )
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)

  protein_id <- lfqdata$relevant_hierarchy_keys()[1]
  fasta_annot <- fasta_annot |>
    dplyr::rename(!!protein_id := !!rlang::sym("Protein"))

  fasta_annot <- site_row_annotation(
    lfqdata = lfqdata,
    config = config,
    long = multiSite_long,
    protein_annot = fasta_annot,
    fasta_file = fasta_file,
    pattern_decoys = pattern_decoys
  )

  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    fasta_annot,
    description = "description",
    cleaned_ids = "protein_Id",
    full_id = "fasta.id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )

  .validate_protein_annotation(lfqdata, prot_annot)
  return(list(lfqdata = lfqdata, protein_annotation = prot_annot))
}
