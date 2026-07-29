# prolfquappPTMreaders 0.3.0

- Package installation now declares its R 4.1 minimum and no longer installs
  the unused `prophosqua` dependency.
- PTM preprocessors now keep FASTA-derived `ProteinAnnotation` objects unique
  by protein ID, allowing protein metadata to propagate to every quantified
  site without treating site metadata as protein annotation.
- PTM preprocessors now work with the `AnalysisConfiguration` objects returned
  by current `prolfquapp` releases.
- Began tracking user-visible changes in `NEWS.md`. For changes before this version, see the git history.
