#' Per-Site Row Annotation for a PTM Site Reader
#'
#' Both FragPipe site readers describe the same thing -- a modified residue in a
#' protein -- so they hand over the same columns, and nothing downstream has to
#' know which quantification the data came from. One row per protein and site,
#' carrying the modified residue, its position, and the sequence window the motif
#' analyses read.
#'
#' The window is always cut from the FASTA, even when the input file ships one,
#' because the shipped ones do not agree: the TMT single-site file carries a
#' 15-mer centred on the residue, the multi-site file carries a peptide with its
#' flanking regions. Cutting it here gives every reader the same convention --
#' `flank_size` residues either side of the site -- which is what
#' `prophosqua::validate_sequence_window()` and the motif analyses expect.
#'
#' The position likewise comes from the index rather than FragPipe's `Start`
#' column, which disagrees with the index for some proteins (`O75531_S4` is
#' reported at `Start = 3`).
#'
#' @param lfqdata The [prolfqua::LFQData] built by the reader.
#' @param config The analysis table annotation the reader set up.
#' @param long The reader's long table, before `setup_analysis`.
#' @param protein_annot Protein-level annotation, keyed on the protein id.
#' @param fasta_file Path to the FASTA, used when the window has to be computed.
#' @param pattern_decoys Decoy pattern, so a window is never cut out of a
#'   reversed sequence.
#' @return `protein_annot` extended to one row per protein and site, with
#'   `modAA`, `posInProtein` and `SequenceWindow`.
#' @keywords internal
site_row_annotation <- function(
  lfqdata,
  config,
  long,
  protein_annot,
  fasta_file,
  pattern_decoys = NULL
) {
  protein_id <- lfqdata$relevant_hierarchy_keys()[[1]]

  site <- long |>
    dplyr::select(dplyr::any_of(c("Protein", "Index", "Peptide"))) |>
    dplyr::distinct()

  parsed <- prolfquappPTMreaders::parse_site_index(site$Index)
  site$modAA <- parsed$modAA
  site$posInProtein <- parsed$posInProtein

  sequences <- prolfquapp::get_annot_from_fasta(
    fasta_file,
    pattern_decoys = pattern_decoys,
    include_seq = TRUE
  )
  if (!is.null(pattern_decoys) && nzchar(pattern_decoys)) {
    # A decoy shares its cleaned id with the protein it was reversed from, so
    # an unfiltered lookup silently cuts windows out of reversed sequences.
    sequences <- sequences |>
      dplyr::filter(!grepl(pattern_decoys, .data$fasta.id))
  }
  sequences <- sequences |>
    dplyr::select(Protein = "proteinname", "sequence") |>
    dplyr::distinct(.data$Protein, .keep_all = TRUE)

  site <- site |>
    dplyr::left_join(sequences, by = "Protein") |>
    prolfquappPTMreaders::get_sequence_windows(
      sequence = "sequence",
      pos_in_protein = "posInProtein"
    ) |>
    dplyr::rename(SequenceWindow = "sequence_window") |>
    dplyr::select(-dplyr::all_of("sequence"))

  site <- site |> dplyr::rename(!!protein_id := !!rlang::sym("Protein"))

  # The site key is the hierarchy's own concatenation of Index and Peptide;
  # separate_hierarchy splits it back rather than assuming the separator.
  site_keys <- lfqdata$data_long() |>
    dplyr::select(dplyr::all_of(lfqdata$relevant_hierarchy_keys())) |>
    dplyr::distinct()
  site_keys <- prolfqua::separate_hierarchy(site_keys, config)

  # Join on the protein as well: one Index and Peptide pair can belong to more
  # than one protein, and joining without it multiplies the annotation.
  site_keys <- dplyr::left_join(
    site_keys,
    site,
    by = c(protein_id, "Index", "Peptide")
  )

  # Many-to-many by design: a protein carries many sites, and the FASTA
  # annotation still holds a forward and a decoy row per cleaned id. The decoy
  # rows are dropped by ProteinAnnotation's decoy-aware resolution, which owns
  # that decision.
  dplyr::left_join(
    site_keys,
    protein_annot,
    by = protein_id,
    relationship = "many-to-many"
  )
}
