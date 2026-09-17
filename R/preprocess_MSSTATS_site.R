#' Read a Site-Level Table in MSstats Long Format
#'
#' The MSstats long format is what the PTM statistics literature publishes in:
#' MSstatsPTM's converters emit it, its simulations are distributed in it, and
#' every tool in that family reads it. This reader is the way into our pipeline
#' for data that arrives that way rather than as a search engine's own output.
#'
#' One row is one site in one run. Required columns:
#'
#' \describe{
#'   \item{ProteinName}{the parent protein.}
#'   \item{Index}{the site, as `<protein>_<residue><position>`, for example
#'     `P07947_K235` -- the convention [parse_site_index()] reads and the one
#'     MSstatsPTM uses for its PTM protein names.}
#'   \item{Run}{the run, joined against the annotation's file name column.}
#'   \item{Intensity}{on the linear scale, not log2.}
#' }
#'
#' `PeptideSequence` is optional and defaults to `Index`: a site-level table has
#' one quantity per site, so there is no peptide to distinguish, and the site
#' hierarchy then keys on the index alone.
#'
#' As with every reader here, the sequence window is cut from the FASTA rather
#' than taken from the input, so `modAA`, `posInProtein` and `SequenceWindow`
#' follow the same convention as the FragPipe and Spectronaut readers. The FASTA
#' therefore has to contain the proteins named in `ProteinName`, with the
#' residue named in `Index` at the position named in `Index`.
#'
#' @param path directory path to search for files
#' @return list with paths to data and fasta
#' @export
#' @examples
#' path <- system.file("extdata", "MSSTATS_site", package = "prolfquappPTMreaders")
#' files <- get_MSSTATS_site_files(path)
#' basename(unlist(files))
get_MSSTATS_site_files <- function(path) {
  data_file <- grep(
    "msstats.*site.*\\.csv$|.*sites.*\\.csv$",
    dir(path = path, recursive = TRUE, full.names = TRUE),
    value = TRUE,
    ignore.case = TRUE
  )
  fasta.files <- grep(
    "\\.fasta$|\\.fas$",
    dir(path = path, recursive = TRUE, full.names = TRUE),
    value = TRUE
  )
  return(list(data = data_file, fasta = fasta.files))
}

#' Read the site table
#'
#' @param file path to the table; comma or tab separated.
#' @return data frame with `PeptideSequence` filled in where absent.
#' @export
#' @examples
#' path <- system.file("extdata", "MSSTATS_site", package = "prolfquappPTMreaders")
#' sites <- read_MSSTATS_site(get_MSSTATS_site_files(path)$data)
#' head(sites)
read_MSSTATS_site <- function(file) {
  x <- readr::read_delim(
    file,
    delim = if (grepl("\\.tsv$|\\.txt$", file[[1]])) "\t" else ",",
    show_col_types = FALSE,
    progress = FALSE
  )
  required <- c("ProteinName", "Index", "Run", "Intensity")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop(
      "the site table is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"PeptideSequence" %in% colnames(x)) {
    # A site-level table has one quantity per site; the hierarchy still wants a
    # second key, and the index is the honest value for it.
    x$PeptideSequence <- x$Index
  }
  x
}

#' Annotation template for a site table
#'
#' @param files output of [get_MSSTATS_site_files()].
#' @return data.frame with annotation template
#' @export
#' @examples
#' path <- system.file("extdata", "MSSTATS_site", package = "prolfquappPTMreaders")
#' dataset_template_MSSTATS_site(get_MSSTATS_site_files(path))
dataset_template_MSSTATS_site <- function(files) {
  sites <- read_MSSTATS_site(files$data)
  data.frame(
    Relative.Path = unique(sites$Run),
    Name = unique(sites$Run),
    Group = "",
    Subject = "",
    Control = "",
    stringsAsFactors = FALSE
  )
}

#' Preprocess a site-level table in MSstats long format
#'
#' @param quant_data path to the site table.
#' @param fasta_file path to the FASTA the windows are cut from.
#' @param annotation output of [prolfquapp::read_annotation()].
#' @param pattern_contaminants contaminant pattern.
#' @param pattern_decoys decoy pattern.
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#' path <- system.file("extdata", "MSSTATS_site", package = "prolfquappPTMreaders")
#' files <- get_MSSTATS_site_files(path)
#' annot <- dataset_template_MSSTATS_site(files)
#' annot$Group <- rep(c("A", "B"), length.out = nrow(annot))
#' annot$Control <- ifelse(annot$Group == "A", "C", "T")
#' annotation <- prolfquapp::read_annotation(annot)
#' result <- preprocess_MSSTATS_site(
#'   quant_data = files$data,
#'   fasta_file = files$fasta,
#'   annotation = annotation
#' )
#' stopifnot(nrow(result$lfqdata$data) > 0)
#' head(result$protein_annotation$row_annot$SequenceWindow)
preprocess_MSSTATS_site <- function(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_"
) {
  config <- annotation$atable$clone(deep = TRUE)

  site_long <- read_MSSTATS_site(quant_data)

  by <- "Run"
  names(by) <- config$file_name
  site_long <- dplyr::inner_join(x = annotation$annot, y = site_long, by = by)
  if (nrow(site_long) == 0) {
    stop(
      "no run in the site table matches the annotation's ",
      config$file_name,
      " column",
      call. = FALSE
    )
  }

  # The columns site_row_annotation() and the site hierarchy key on.
  site_long$Protein <- site_long$ProteinName
  site_long$Peptide <- site_long$PeptideSequence
  site_long$qValue <- 0
  site_long$nr_children <- 1
  site_long$isotopeLabel <- "light"

  config$isotope_label <- "isotopeLabel"
  config$ident_q_value <- "qValue"
  config$nr_children <- "nr_children"
  config$hierarchy[["protein_Id"]] <- c("Protein")
  config$hierarchy[["site"]] <- c("Index", "Peptide")
  config$set_response("Intensity")
  config$hierarchy_depth <- 2

  adata <- prolfqua::setup_analysis(site_long, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities(threshold = 1)

  nrPep_exp <- site_long |>
    dplyr::select("Protein", "Peptide") |>
    dplyr::distinct() |>
    dplyr::group_by(.data$Protein) |>
    dplyr::summarize(nrPeptides = dplyr::n()) |>
    dplyr::ungroup()

  fasta_annot <- prolfquapp::get_annot_from_fasta(
    fasta_file,
    pattern_decoys = pattern_decoys
  )
  fasta_annot <- dplyr::left_join(
    nrPep_exp,
    fasta_annot,
    by = c(Protein = "proteinname"),
    multiple = "all"
  )
  fasta_annot <- fasta_annot |> dplyr::rename(description = "fasta.header")

  protein_id <- lfqdata$relevant_hierarchy_keys()[1]
  fasta_annot <- fasta_annot |>
    dplyr::rename(!!protein_id := !!rlang::sym("Protein"))

  fasta_annot <- site_row_annotation(
    lfqdata = lfqdata,
    config = config,
    long = site_long,
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
