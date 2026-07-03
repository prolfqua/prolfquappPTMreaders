# TODO: expose the stripped (unmodified) peptide so the nr_peptides filter works for PTM/site data

> Add the stripped-sequence peptide to the PTM/site readers so prolfquapp's
> minimum-peptides-per-parent-protein filter can be applied to site-level data.

## Requirements

- **Why:** prolfquapp is adding a minimum-peptides-per-parent-protein filter —
  drop a protein (and its child peptides/sites) when its parent protein is
  supported by fewer than **N distinct stripped (unmodified) peptides**. See
  `prolfquapp/TODO/TODO_expose_nr_peptides.md`. It runs generically on the
  `LFQData`'s `protein_Id` → stripped-peptide hierarchy, counting distinct
  stripped peptides per protein.
- For standard readers the `LFQData` child **is** the stripped peptide, so the
  count is a direct group-by. **PTM/site readers are the exception**: their
  hierarchy is `protein_Id` → **site**, not `protein_Id` → stripped peptide, so
  there is no stripped-sequence child to count — counting children would count
  *sites*, giving the wrong cut.
- **Goal:** PTM/site `LFQData` carries the stripped peptide so "distinct stripped
  peptides per parent protein" is computable, and the same filter can drop sites
  whose parent protein has `< N` stripped peptides.
- Behaviour unchanged at the default `N = 1`.
- **Near-term (before the stripped-peptide element exists):** prolfquapp will
  forward `nr_peptides` into readers via `formals()` inspection and *warn* for
  readers that don't declare it (see `prolfquapp/TODO/TODO_expose_nr_peptides.md`).
  So each PTM/site reader should **accept an `nr_peptides = 1` parameter now**
  (accepted and ignored) — this keeps the forwarded `do.call` safe and silences
  the unsupported-reader warning until the real filter is implemented here.

## Current state (verified in code)

- `R/preprocess_FP_combined_STY.R:155-158`:
  `config$hierarchy[["protein_Id"]] <- c("Protein")`,
  `config$hierarchy[["site"]] <- c("Index", "Peptide")`, `hierarchy_depth <- 2`;
  and `nr_children <- 1` is hardcoded (`:148`).
- The FP multisite/singlesite (`R/preprocess_FP_multisite.R`) and BGS site
  (`R/preprocess_BGS_site.R`) readers use the same site-keyed hierarchy.
- A `Peptide` column exists inside the site key, but it is (likely) the
  **modified** peptide, and there is no per-protein distinct-stripped-peptide
  count.

## Design (sketch)

- Add the **stripped (unmodified) peptide sequence** as a hierarchy element (or,
  at minimum, a carried row column) in each PTM/site reader, alongside the
  existing `site` key.
- prolfquapp's generic filter should take the **stripped-peptide column name as a
  parameter** (default = the standard peptide hierarchy key). PTM readers then
  just need to expose that column; standard readers keep using the peptide key.
- Decide whether the stripped sequence is a true hierarchy key (affects
  aggregation / site rollup) or only a filter-support column.

## Implementation plan (later)

- [ ] For each PTM/site reader (`preprocess_FP_multi_site`,
      `preprocess_FP_combined_STY`, `preprocess_BGS_site`), determine the source
      of the (modified) peptide and derive the stripped sequence.
- [ ] Add the stripped-sequence column/hierarchy element to the `config` /
      long table before `setup_analysis()`.
- [ ] Coordinate the column name with prolfquapp's generic filter parameter
      (`TODO_expose_nr_peptides.md`).
- [ ] Tests: a PTM fixture where a parent protein has a single stripped peptide
      (possibly several sites) is dropped at `N = 2` and kept at `N = 1`.

## Open questions

- Is the `Peptide` in the site key the modified or the stripped sequence? If
  modified, what is the stripping rule (strip modifications / charge)?
- Stripped sequence as a real hierarchy key vs. a filter-only column — which,
  given site aggregation?
- Does `prophosqua` have its own PTM path that needs the same treatment?
