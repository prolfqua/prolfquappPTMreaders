# Global variables to avoid R CMD check NOTEs for NSE in dplyr/tidyr
# These are column names used in dplyr/tidyr operations
utils::globalVariables(c(

  # Common columns
  "Index", "Peptide", "ProteinID", "fasta.header", "sequence_window", "siteinfo",
  "modAA", "posInProtein",


  # BGS-specific columns
  "PTM_ProteinId", "PTM_CollapseKey", "PTM_FlankingRegion", "PTM_SiteAA",
 "PTM_SiteLocation", "PTM_ModificationTitle", "PTM_Multiplicity",
  "PTM_SiteProbability", "PTM_Group",

  # tidyselect functions
  "all_of",

  # rlang operators
  ":="
))
