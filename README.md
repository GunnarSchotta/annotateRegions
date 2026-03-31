# annotate.regions

`annotate_regions.R` is an R script for annotating genomic intervals from a BED file against gene features and RepeatMasker repeat annotations.

It supports both mouse and human genomes and produces tabular summaries plus publication-friendly plots.

## What The Script Does

Given a BED file of genomic regions, the script:

- reads UCSC-style BED coordinates (`chr1`, `chr2`, etc.; 0-based start)
- annotates each region relative to gene structure
- assigns nearby or overlapping gene information
- summarizes transposable element and repeat overlap
- compares the observed regions against matched randomized genomic background regions
- writes result tables and PNG plots

The repeat analysis excludes several non-TE repeat classes such as simple repeats, low-complexity regions, satellite annotations, and small RNA-related categories so the repeat summaries stay focused on TE-derived annotations.

## Supported Genomes

- `mm10`
- `mm39`
- `hg19`
- `hg38`

## Input Requirements

The input must be a BED file with at least 3 tab-separated columns:

1. chromosome
2. start
3. end

Important assumptions:

- BED coordinates are expected to be standard UCSC/BED format with 0-based starts.
- Chromosome names should use UCSC-style naming such as `chr1`, `chrX`, or `chrY`.
- The script expects standard chromosomes for the selected genome assembly.

Additional BED columns are preserved and carried through to the annotated output table.

## Usage

```bash
Rscript annotate_regions.R <input.bed> <output_prefix> [top_n_repeat_names] \
  [--genome mm10|mm39|hg19|hg38] \
  [--repeat-summary-mode dominant|all] \
  [--background-iterations N]
```

### Positional Arguments

- `<input.bed>`: input BED file
- `<output_prefix>`: prefix used for all generated output files
- `[top_n_repeat_names]`: optional positive integer controlling how many repeat names are shown in the top-repeat plot; default is `10`

### Optional Flags

- `--genome`: genome assembly to use; default is `mm10`
- `--repeat-summary-mode`: how repeat class/name plots are summarized
- `--background-iterations`: number of randomized background samplings; default is `100`

`--repeat-summary-mode` options:

- `dominant`: one dominant repeat assignment per region is used for class/name plots
- `all`: all repeat overlaps are used for class/name plots

## Example

```bash
Rscript annotate_regions.R peaks.bed peaks_mm39 10 \
  --genome mm39 \
  --repeat-summary-mode dominant \
  --background-iterations 100
```

## Output Files

Using an output prefix such as `peaks_mm39`, the script generates:

- `peaks_mm39.annotated.tsv`: main region-level annotation table
- `peaks_mm39.repeat_detail.tsv`: detailed repeat overlap information
- `peaks_mm39.background_summary.tsv`: observed vs. background summary statistics
- `peaks_mm39.feature_annotation_percent.png`: feature annotation percentages
- `peaks_mm39.repeat_association_percent.png`: repeat-associated vs. non-repeat-associated fractions
- `peaks_mm39.repeat_class_frequency.png`: repeat class frequency plot
- `peaks_mm39.repeat_name_topN_frequency.png`: top repeat-name frequency plot

## Dependencies

The script uses CRAN and Bioconductor packages, including:

- `ggplot2`
- `GenomicRanges`
- `IRanges`
- `S4Vectors`
- `GenomeInfoDb`
- `GenomicFeatures`
- `AnnotationDbi`
- `AnnotationHub`
- `UCSCRepeatMasker`
- assembly-specific `TxDb` and `org.*.eg.db` packages

By default, the script attempts to install missing packages automatically.

## Notes

- Promoters are defined in the script using a window of 2000 bp upstream and 200 bp downstream.
- The script validates that `end >= start` for all BED intervals.
- `top_n_repeat_names` and `--background-iterations` must both be positive integers.

## Repository Contents

- `annotate_regions.R`: main annotation script
