

#' get Spectronaut report.tsv and fasta file location in folder
#' @param path directory path to search for files
#' @return list with paths to data and fasta
#' @export
#' @examples
#' identical(formals(get_BGS_site_files), formals(prolfquapp::get_dummy_files))
#' # Use package test data
#' path <- system.file("extdata", "BGS_site", package = "prolfquappPTMreaders")
#' files <- get_BGS_site_files(path)
#' files$data
#' files$fasta
#'
get_BGS_site_files <- function(path){
  # Find Spectronaut PTM report files (typically named *Report_WithProteinRollup.tsv)
  psm_file <- dir(path = path, pattern = "Report.*\\.tsv$", recursive = TRUE, full.names = TRUE)
  fasta.files <- grep("*.fasta$|*.fas$", dir(path = path, recursive = TRUE,full.names = TRUE), value = TRUE)
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_warn("No fasta file found!")
  }
  return(list(data = psm_file, fasta = fasta.files))
}


#' read Spectronaut PTM report file
#' @export
#' @param quant_data path to Spectronaut PTM report file
#' @param min_site_sample_loc minimum site probability threshold for individual observations (default: 0.00)
#' @param min_site_loc minimum site probability threshold - keep only sites where at least one sample exceeds this (default: 0.99)
#' @return data table with Spectronaut column names (dots replaced with underscores),
#'   filtered to keep only single phosphorylation sites (PTM_Multiplicity == 1) and high confidence observations
#' @examples
#' # Use package test data
#' path <- system.file("extdata", "BGS_site", package = "prolfquappPTMreaders")
#' files <- get_BGS_site_files(path)
#' data <- read_BGS_site(files$data, min_site_loc = 0.75)
#' nrow(data)
#' range(data$PTM_SiteProbability)
#' unique(data$PTM_ModificationTitle)
#'
read_BGS_site <- function(quant_data, min_site_sample_loc=0.3, min_site_loc = 0.95) {
  xx <- readr::read_tsv(quant_data, show_col_types = FALSE)

  # Replace dots with underscores in column names for R compatibility
  colnames(xx) <- gsub("\\.", "_", colnames(xx))

  # Filter for phosphorylation sites only with multiplicity = 1 (single phosphorylation)
  xx <- xx |>
    dplyr::filter(grepl("Phospho", PTM_ModificationTitle)) |>
    dplyr::filter(PTM_Multiplicity == 1)

  # Two-step filtering for site probability:
  # 1. Remove very low confidence observations
  # 2. Keep sites where at least one sample has very high confidence
  xx <- xx |>
    dplyr::filter(PTM_SiteProbability >= min_site_sample_loc) |>
    dplyr::group_by(PTM_CollapseKey) |>
    dplyr::filter(max(PTM_SiteProbability) > min_site_loc) |>
    dplyr::ungroup()

  return(xx)
}


#' create dataset template from BGS site file
#' @export
#' @param files list with data and fasta paths
#' @return data.frame with annotation template
#' @examples
#' # Use package test data
#' path <- system.file("extdata", "BGS_site", package = "prolfquappPTMreaders")
#' files <- get_BGS_site_files(path)
#' annot_template <- dataset_template_BGS_site(files)
#' head(annot_template)
#'
dataset_template_BGS_site <- function(files){
  xx <- read_BGS_site(files$data)
  channels <- unique(xx$R_FileName)
  dataset <- data.frame(raw.file = channels, Name = channels, Group = NA, Control = NA)
  return(dataset)
}


#' preprocess Spectronaut PTM site data
#' @param quant_data path to Spectronaut PTM report file
#' @param fasta_file path to FASTA file
#' @param annotation annotation object from prolfquapp::read_annotation()
#' @param pattern_contaminants regex pattern to identify contaminant proteins
#' @param pattern_decoys regex pattern to identify decoy proteins
#' @return list with lfqdata and protein annotation
#' @export
#' @examples
#' identical(names(formals(preprocess_BGS_site)), names(formals(prolfquapp::preprocess_dummy)))
#' # Use package test data (single sample subset)
#' path <- system.file("extdata", "BGS_site", package = "prolfquappPTMreaders")
#' files <- get_BGS_site_files(path)
#'
#' # Create annotation
#' annot <- dataset_template_BGS_site(files)
#' # Fill in Group and Control columns (test data has single sample)
#' annot$Group <- "A"
#' annot$Control <- "C"
#' annotation <- prolfquapp::read_annotation(annot)
#'
#' # Preprocess
#' result <- preprocess_BGS_site(
#'   quant_data = files$data,
#'   fasta_file = files$fasta,
#'   annotation = annotation
#' )
#' result$lfqdata
#' result$protein_annotation
#' }
preprocess_BGS_site <- function(
    quant_data,
    fasta_file,
    annotation,
    pattern_contaminants = "^zz|^CON|Cont_",
    pattern_decoys = "^REV_|^rev_"){

  annot <- annotation$annot
  atable <- annotation$atable
  annot <- annot |> dplyr::mutate(
    raw.file = gsub("^x|.d.zip$|.raw$","",
                    (basename(annot[[atable$fileName]]))
    ))

  site_long <- read_BGS_site(quant_data)

  # Rename only what's needed for prolfqua hierarchy and joining

  bw <- "R_FileName"
  names(bw) <- annotation$atable$fileName

  # join with annotation
  site_long <- dplyr::inner_join(x = annotation$annot, y = site_long, by = bw)

  # add missing required parameters (qvalue)
  # Use 1 - PTM_SiteProbability as qValue (lower is better)
  site_long$qValue <- 1 - site_long$PTM_SiteProbability
  site_long$nr_children  <- 1


  # Setup configuration for site-level analysis (phospho)
  # Note: BGS data is already aggregated at site level by Spectronaut
  # PTM_Group contains which peptides were used but there's only 1 quantity per site-sample
  atable$ident_Score = "PTM_SiteProbability"
  atable$ident_qValue = "qValue"
  atable$nr_children = "nr_children"
  atable$hierarchy[["protein_Id"]] <- c("PTM_ProteinId")
  atable$hierarchy[["site"]] <- c("PTM_ProteinId","PTM_CollapseKey", "PTM_SiteAA", "PTM_SiteLocation", "PTM_Multiplicity")

  atable$set_response("PTM_Quantity")
  atable$hierarchyDepth <- 2

  # Preprocess data - setup analysis
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(site_long, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities(threshold = 1)


  # Create Site Annotation - one row per site with max PTM_SiteProbability
  site_annot <- site_long |>
    dplyr::group_by(PTM_ProteinId, PTM_CollapseKey, PTM_FlankingRegion,
                    PTM_SiteAA, PTM_SiteLocation,
                    PTM_ModificationTitle, PTM_Multiplicity) |>
    dplyr::summarize(PTM_SiteProbability = max(PTM_SiteProbability, na.rm = TRUE), .groups = "drop") |>
    dplyr::distinct()

  # Rename columns to match expected naming convention and create PhosSites


  # Count peptides per protein
  nrPep_exp <- site_long |>
    dplyr::select(PTM_ProteinId, PTM_Group) |>
    dplyr::distinct() |>
    dplyr::group_by(PTM_ProteinId) |>
    dplyr::summarize(nrPeptides = dplyr::n()) |> dplyr::ungroup()

  fasta_annot <- prolfquapp::get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
  fasta_annot <- dplyr::left_join(nrPep_exp, fasta_annot, by = c(PTM_ProteinId = "proteinname"), multiple = "all")
  fasta_annot <- fasta_annot |> dplyr::rename(description = fasta.header)
  fasta_annot2 <- dplyr::inner_join(fasta_annot, site_annot, by = "PTM_ProteinId")

  # Make names to match lfqdata - must unite same columns as setup_analysis does
  fasta_annot2 <- fasta_annot2 |> dplyr::rename(!!lfqdata$config$table$hierarchy_keys_depth()[1] := !!rlang::sym("PTM_ProteinId"))
  fasta_annot2 <- fasta_annot2 |> tidyr::unite(
    !!lfqdata$config$table$hierarchy_keys_depth()[2],
    c("protein_Id", "PTM_CollapseKey", "PTM_SiteAA", "PTM_SiteLocation", "PTM_Multiplicity"),
    sep = "~", remove = FALSE
  )

  fasta_annot2 <- fasta_annot2 |>
    dplyr::rename(
      modAA = PTM_SiteAA,
      posInProtein = PTM_SiteLocation,
      SequenceWindow = PTM_FlankingRegion
    ) |>
    dplyr::mutate(PhosSites = paste0(modAA, posInProtein))


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

  #verify
  stopifnot(nrow(dplyr::inner_join(prot_annot$row_annot, lfqdata$data, by = lfqdata$config$table$hierarchy_keys_depth())) > 0)
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}
