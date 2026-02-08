
#' get report.tsv and fasta file location in folder
#' @param path directory path to search for files
#' @return list with paths to data and fasta
#' @export
#' @examples
#' identical(formals(get_FP_multi_site_files), formals(prolfquapp::get_dummy_files))
#'
get_FP_single_site_files <- function(path){
  psm_file <- dir(path = path, pattern = "abundance_single-site_None.tsv", recursive = TRUE, full.names = TRUE)
  fasta.files <- grep("*.fasta$|*.fas", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = psm_file, fasta = fasta.files))
}


#' get report.tsv and fasta file location in folder
#' @param path directory path to search for files
#' @return list with paths to data and fasta
#' @export
#' @examples
#' identical(formals(get_FP_multi_site_files), formals(prolfquapp::get_dummy_files))
#'
get_FP_multi_site_files <- function(path){
  psm_file <- dir(path = path, pattern = "abundance_multi-site_None.tsv", recursive = TRUE, full.names = TRUE)
  fasta.files <- grep("*.fasta$|*.fas", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }
  return(list(data = psm_file, fasta = fasta.files))
}



#' read multisite file
#' @param quant_data path to multisite quantification file
#' @return tidy long-format data frame
#' @export
read_FP_multisite_to_long <- function(quant_data) {
  xx <- readr::read_tsv(quant_data)
  quant_idx_start <- grep(pattern = "ReferenceIntensity", x = colnames(xx))  + 1
  multiSite_long <- xx |>
    tidyr::pivot_longer(cols = all_of(quant_idx_start:ncol(xx)), values_to = "abundance", names_to = "channel")

  return(multiSite_long)
}


#' create dataset template from FP multi_site
#' @param files list with data and fasta paths
#' @return data.frame with annotation template
#' @export
dataset_template_FP_multi_site <- function(files){
  xx <- read_FP_multisite_to_long(files$data)
  channels <- unique(xx$channel)
  dataset <- data.frame(raw.file = channels, Name = channels ,Group = NA, Subject = NA, CONTROL = NA)
  return(dataset)
}


#' preprocess FP multisite, filter by purity_threshold and PeptideProphetProb
#' @param quant_data path to multisite quantification file
#' @param fasta_file path to FASTA file
#' @param annotation annotation object from prolfquapp::read_annotation()
#' @param sitetype type of site analysis ("singlesite" or "multisite")
#' @param pattern_contaminants regex pattern to identify contaminant proteins
#' @param pattern_decoys regex pattern to identify decoy proteins
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#' identical(names(formals(preprocess_FP_multi_site)), names(formals(prolfquapp::preprocess_dummy)))
#'
preprocess_FP_multi_site <- function(
    quant_data,
    fasta_file,
    annotation,
    sitetype = c("singlesite", "multisite"),
    pattern_contaminants = "^zz|^CON|Cont_",
    pattern_decoys = "^REV_|^rev_"){
  sitetype <- match.arg(sitetype)
  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))

  multiSite_long <- read_FP_multisite_to_long(quant_data)
  bw <- "channel"

  names(bw) <- annotation$atable$fileName

  # join with anno again this should work now with Name # if not all samples are used in the dataset they would be removed here (to be tested)

  multiSite_long <- dplyr::inner_join(x = annotation$annot, y = multiSite_long, by = bw)
  # add missing required parameters (qvalue)
  multiSite_long$qValue <- 1 - multiSite_long$MaxPepProb
  multiSite_long$nr_children  <- 1

  # Map short ProteinID to full FASTA ID (sp|ACC|NAME) for consistency with preprocess_FP_PSM
  fasta_annot_early <- prolfquapp::get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
  id_map <- fasta_annot_early |> dplyr::select(proteinname, fasta.id) |> dplyr::distinct()
  multiSite_long <- dplyr::left_join(multiSite_long, id_map, by = c("ProteinID" = "proteinname"))
  # Use full FASTA ID where available, fall back to ProteinID
  multiSite_long$Protein <- ifelse(is.na(multiSite_long$fasta.id), multiSite_long$ProteinID, multiSite_long$fasta.id)

  # Setup configuration manually for peptide analysis (phospho)
  atable$ident_Score = "MaxPepProb"
  atable$ident_qValue = "qValue"
  atable$nr_children = "nr_children"
  atable$hierarchy[["protein_Id"]] <- c("Protein")
  atable$hierarchy[["site"]] <- c("Index", "Peptide")
  atable$set_response("abundance")
  atable$hierarchyDepth <- 2

  # Preprocess data - aggregate proteins.
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(multiSite_long, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities(threshold = 1)


  # Create fasta annotation
  # Create Site Annotation
  site_annot <- multiSite_long |>
    dplyr::select(c("Index", "ProteinID", "Peptide", "SequenceWindow", "Start", "End", "MaxPepProb", "ReferenceIntensity")) |>
    dplyr::distinct()
  phosSite <- site_annot |> dplyr::rowwise() |> dplyr::mutate(siteinfo = gsub(ProteinID, "", Index))

  if(sitetype == "singlesite"){
    phosSite <- phosSite |>
      tidyr::separate_wider_delim(siteinfo, names = c(NA, "PhosSites"), delim = "_",
                                  too_few = "align_start")
    phosSite <- phosSite |>
      tidyr::extract("PhosSites", into = c("modAA", "posInProtein"), regex = "([A-Z])(\\d+)", convert = TRUE, remove = FALSE)
  } else if(sitetype == "multisite") {
    phosSite <- phosSite |>
      tidyr::separate_wider_delim(siteinfo, names = c(NA, "startModSite", "endModSite", "NumPhos", "LocalizedNumPhos", "PhosSites"), delim = "_",
                                  too_few = "align_start")
    split_codes <- function(x) {
      if (is.na(x)) return(NA)
      return(gsub("([A-Z]\\d+)(?=[A-Z]\\d+)", "\\1;", x, perl = TRUE))
    }
    phosSite$PhosSites <- sapply(phosSite$PhosSites, split_codes)

  }



  nrPep_exp <- multiSite_long |>
    dplyr::select(Protein, Peptide) |>
    dplyr::distinct() |>
    dplyr::group_by(Protein) |>
    dplyr::summarize(nrPeptides = dplyr::n()) |> dplyr::ungroup()

  fasta_annot <- dplyr::left_join(nrPep_exp, fasta_annot_early, by = c(Protein = "fasta.id"), multiple = "all")
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)
  fasta_annot2 <- dplyr::inner_join(fasta_annot, phosSite, by = c("proteinname" = "ProteinID"))

  # Make names to match lfqdata
  fasta_annot2 <- fasta_annot2 |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!rlang::sym("Protein"))
  fasta_annot2 <- fasta_annot2 |> dplyr::mutate(!!lfqdata$config$table$hierarchy_keys_depth()[2] := paste(!!rlang::sym("Index"),!!rlang::sym("Peptide"), sep = "~"))
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata ,
    fasta_annot2,
    description = "description",
    cleaned_ids = "proteinname",
    full_id = "protein_Id",
    exp_nr_children = "nrPeptides",
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )

  # Verify lfqdata and protein_annotation match on hierarchy keys
  stopifnot(nrow(dplyr::inner_join(prot_annot$row_annot, lfqdata$data, by = lfqdata$config$table$hierarchy_keys_depth())) > 0)
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}
