# prolfquappPTMreaders 0.3.0

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
