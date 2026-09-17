# prolfquappPTMreaders 0.3.0

- New `MSSTATS_site` reader: a site-level table in MSstats long format
  (`ProteinName`, `Index`, `Run`, `Intensity`) now enters the pipeline
  directly. That is the format the PTM statistics literature publishes in --
  MSstatsPTM's converters emit it and its simulations are distributed in it --
  so data from that world no longer has to be turned into a search engine's
  output format first. As with the other site readers, the sequence window is
  cut from the FASTA, so the modified residue, its position and the window
  follow the same convention as the FragPipe and Spectronaut readers.
- Both FragPipe site readers now hand over the same per-site annotation --
  modified residue, position in the protein, and sequence window -- so a
  downstream analysis no longer has to tell the TMT and LFQ quantifications
  apart. The window is always cut from the FASTA, using one flank convention
  rather than the differing ones the input files ship.
- Sequence windows are reported as missing for a site whose protein is absent
  from the FASTA, instead of a window of placeholder residues that read like a
  real one.
- New `get_sequence_windows()` (moved here from `prophosqua`, which is where the
  FASTA is read) and `parse_site_index()` for reading the residue and position
  out of a FragPipe site index.
- Package installation now declares its R 4.1 minimum and no longer installs
  the unused `prophosqua` dependency.
- PTM preprocessors now keep FASTA-derived `ProteinAnnotation` objects unique
  by protein ID, allowing protein metadata to propagate to every quantified
  site without treating site metadata as protein annotation.
- PTM preprocessors now work with the `AnalysisConfiguration` objects returned
  by current `prolfquapp` releases.
- Began tracking user-visible changes in `NEWS.md`. For changes before this version, see the git history.
