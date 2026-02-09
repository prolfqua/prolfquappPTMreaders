# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

prolfquappPTMreaders is an R package that extends prolfquapp to support Post-Translational Modification (PTM) analysis. It provides reader functions and preprocessing pipelines for PTM output formats from FragPipe and Spectronaut (BGS).

## Development Commands

```r
devtools::document()   # Build documentation (roxygen2)
devtools::install()    # Install locally
devtools::check()      # R CMD check
```

A `Makefile` is provided with convenient targets. **Run `make check-fast` (or `make check`) periodically** to catch errors, warnings, and NOTEs early:

```bash
make check-fast   # R CMD check without vignettes (quick)
make check        # Full R CMD check including vignettes
make test         # Run testthat tests only
make lint         # Run lintr
```

## Supported Formats

The package handles four PTM data sources (see `R/prolfqua_preprocess_functions.R` for config):

| Format | File Pattern | Reader | Source |
|--------|--------------|--------|--------|
| FP_multisite | `abundance_multi-site_None.tsv` | `preprocess_FP_multi_site.R` | FragPipe |
| FP_singlesite | `abundance_single-site_None.tsv` | `preprocess_FP_multi_site.R` | FragPipe |
| FP_combined_STY | `combined_site_STY_*.tsv` | `preprocess_FP_combined_STY.R` | FragPipe |
| BGS_site | `Report*.tsv` | `preprocess_BGS_site.R` | Spectronaut |

## Architecture

Each format follows the same pattern with four exported functions:

1. **`get_*_files(path)`** - Locate data files and FASTA in a directory
2. **`read_*(...)`** - Parse raw output to tidy long format
3. **`dataset_template_*(files)`** - Generate annotation template with Group/Subject/Control columns
4. **`preprocess_*(quant_data, fasta_file, annotation, ...)`** - Full preprocessing pipeline returning `list(lfqdata, protein_annotation)`

### Data Flow

```
get_*_files() → dataset_template_*() → [user fills annotation] → preprocess_*() → lfqdata + protein_annotation
```

### Hierarchy Configuration

All PTM data uses 2-level hierarchy:
- Level 1: Protein ID
- Level 2: Site (composite key: `Index~Peptide` or equivalent)

## Function Signature Compatibility

**IMPORTANT**: Functions must maintain signature compatibility with prolfquapp:

- `preprocess_*()` must match `prolfquapp::preprocess_dummy()` signature
- `get_*_files()` must match `prolfquapp::get_dummy_files()` signature (single `path` parameter, returns list with `data`, `fasta`, optionally `fp.manifest`)

Verify compatibility in examples:
```r
identical(names(formals(preprocess_BGS_site)), names(formals(prolfquapp::preprocess_dummy)))
identical(formals(get_FP_multi_site_files), formals(prolfquapp::get_dummy_files))
```

## Integration with prolfquapp

Reference in `config.yaml`:

```yaml
# FragPipe multi-site
ext_reader:
  extra_args: list(sitetype = 'multisite')
  preprocess: prolfquappPTMreaders::preprocess_FP_multi_site
  get_files: prolfquappPTMreaders::get_FP_multi_site_files

# FragPipe combined STY
ext_reader:
  extra_args: list(annotation_join_by = 'SampleName')
  preprocess: prolfquappPTMreaders::preprocess_FP_combined_STY
  get_files: prolfquappPTMreaders::get_FP_combined_STY_files

# Spectronaut BGS
ext_reader:
  extra_args: list()
  preprocess: prolfquappPTMreaders::preprocess_BGS_site
  get_files: prolfquappPTMreaders::get_BGS_site_files
```

## Key Dependencies

- **prolfquapp/prolfqua**: Core proteomics framework
- **prophosqua**: Phosphoproteomics utilities (e.g., `get_sequence_windows()`)
