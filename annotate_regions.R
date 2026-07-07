#!/usr/bin/env Rscript

# Automatic annotation of mouse or human genomic regions from a BED file.
# Input BED is assumed to be UCSC/BED 0-based with chromosome names like chr1.
#
# Supported genome assemblies:
#   mm10, mm39, hg19, hg38
#
# Outputs:
#   <prefix>.annotated.tsv
#   <prefix>.repeat_detail.tsv
#   <prefix>.background_summary.tsv
#   <prefix>.feature_annotation_percent.png / .pdf
#   <prefix>.repeat_class_frequency.png / .pdf
#   <prefix>.repeat_name_topN_frequency.png / .pdf
#
# Usage:
#   Rscript annotate_regions_background_repeatmode.R regions.bed output_prefix [top_n_repeat_names]
#       [--genome mm10|mm39|hg19|hg38]
#       [--repeat-summary-mode dominant|all]
#       [--background-iterations N]
# Example:
#   Rscript annotate_regions_background_repeatmode.R peaks.bed peaks_mm39 10 \
#       --genome mm39 --repeat-summary-mode dominant --background-iterations 100

VERSION <- "1.1.0"

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

parse_args <- function(args) {
  opts <- list(
    genome = "mm10",
    repeat_summary_mode = "dominant",
    background_iterations = 100L
  )
  positional <- character()

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--version", "-v")) {
      cat(paste0("annotate_regions.R version ", VERSION, "\n"))
      quit(save = "no", status = 0)
    } else if (arg == "--genome") {
      if (i == length(args)) stop("Missing value after --genome")
      opts$genome <- args[[i + 1L]]
      i <- i + 2L
    } else if (grepl("^--genome=", arg)) {
      opts$genome <- sub("^--genome=", "", arg)
      i <- i + 1L
    } else if (arg == "--repeat-summary-mode") {
      if (i == length(args)) stop("Missing value after --repeat-summary-mode")
      opts$repeat_summary_mode <- args[[i + 1L]]
      i <- i + 2L
    } else if (grepl("^--repeat-summary-mode=", arg)) {
      opts$repeat_summary_mode <- sub("^--repeat-summary-mode=", "", arg)
      i <- i + 1L
    } else if (arg == "--background-iterations") {
      if (i == length(args)) stop("Missing value after --background-iterations")
      opts$background_iterations <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (grepl("^--background-iterations=", arg)) {
      opts$background_iterations <- as.integer(sub("^--background-iterations=", "", arg))
      i <- i + 1L
    } else if (startsWith(arg, "--")) {
      stop("Unknown option: ", arg)
    } else {
      positional <- c(positional, arg)
      i <- i + 1L
    }
  }

  list(positional = positional, options = opts)
}

args_raw <- commandArgs(trailingOnly = TRUE)
parsed <- parse_args(args_raw)
args <- parsed$positional
opts <- parsed$options

if (length(args) < 2) {
  stop(
    paste0(
      "Usage: Rscript annotate_regions_background_repeatmode.R <input.bed> <output_prefix> [top_n_repeat_names] ",
      "[--genome mm10|mm39|hg19|hg38] [--repeat-summary-mode dominant|all] [--background-iterations N]
",
      "Example: Rscript annotate_regions_background_repeatmode.R peaks.bed peaks_mm39 10 ",
      "--genome mm39 --repeat-summary-mode dominant --background-iterations 100"
    )
  )
}

bed_file <- args[[1]]
out_prefix <- args[[2]]
top_n_repeat_names <- if (length(args) >= 3) as.integer(args[[3]]) else 10L
if (is.na(top_n_repeat_names) || top_n_repeat_names < 1L) {
  stop("top_n_repeat_names must be a positive integer.")
}

GENOME_CONFIG <- list(
  mm10 = list(
    species = "Mus musculus",
    genome = "mm10",
    txdb_pkg = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    txdb_obj = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    orgdb_pkg = "org.Mm.eg.db",
    orgdb_obj = "org.Mm.eg.db",
    standard_chroms = c(paste0("chr", c(1:19, "X", "Y")))
  ),
  mm39 = list(
    species = "Mus musculus",
    genome = "mm39",
    txdb_pkg = "TxDb.Mmusculus.UCSC.mm39.knownGene",
    txdb_obj = "TxDb.Mmusculus.UCSC.mm39.knownGene",
    orgdb_pkg = "org.Mm.eg.db",
    orgdb_obj = "org.Mm.eg.db",
    standard_chroms = c(paste0("chr", c(1:19, "X", "Y")))
  ),
  hg19 = list(
    species = "Homo sapiens",
    genome = "hg19",
    txdb_pkg = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    txdb_obj = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    orgdb_pkg = "org.Hs.eg.db",
    orgdb_obj = "org.Hs.eg.db",
    standard_chroms = c(paste0("chr", c(1:22, "X", "Y")))
  ),
  hg38 = list(
    species = "Homo sapiens",
    genome = "hg38",
    txdb_pkg = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    txdb_obj = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    orgdb_pkg = "org.Hs.eg.db",
    orgdb_obj = "org.Hs.eg.db",
    standard_chroms = c(paste0("chr", c(1:22, "X", "Y")))
  )
)

selected_genome <- tolower(opts$genome)
if (!(selected_genome %in% names(GENOME_CONFIG))) {
  stop(
    "--genome must be one of: ",
    paste(names(GENOME_CONFIG), collapse = ", "),
    ". Received: ", opts$genome
  )
}
genome_cfg <- GENOME_CONFIG[[selected_genome]]

# Focus repeat analyses on TE-derived annotations by excluding non-TE repeat classes
# before any overlap calculations. This affects repeat_any, dominant repeat
# assignment, detailed overlap tables, and all observed/background summaries.
EXCLUDED_REPEAT_CLASSES <- c(
  "Simple_repeat",
  "Low_complexity",
  "Satellite",
  "rRNA",
  "scRNA",
  "snRNA",
  "srpRNA",
  "tRNA",
  "LTR?",
  "DNA?",
  "RNA",
  "RC",
  "SINE?",
  "LINE?",
  "RC?"
)

repeat_summary_mode <- tolower(opts$repeat_summary_mode)
if (!(repeat_summary_mode %in% c("dominant", "all"))) {
  stop("--repeat-summary-mode must be either 'dominant' or 'all'.")
}
background_iterations <- as.integer(opts$background_iterations)
if (is.na(background_iterations) || background_iterations < 1L) {
  stop("--background-iterations must be a positive integer.")
}

PROMOTER_UPSTREAM <- 2000L
PROMOTER_DOWNSTREAM <- 200L
AUTO_INSTALL <- TRUE
STANDARD_CHROMS <- genome_cfg$standard_chroms

cran_pkgs <- c("ggplot2")
bioc_pkgs <- c(
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "GenomeInfoDb",
  "GenomicFeatures",
  "AnnotationDbi",
  genome_cfg$txdb_pkg,
  genome_cfg$orgdb_pkg,
  "AnnotationHub",
  "UCSCRepeatMasker"
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) == 0) return(invisible(TRUE))
  if (!AUTO_INSTALL) {
    stop(
      "Missing packages: ", paste(missing_pkgs, collapse = ", "),
      ". Set AUTO_INSTALL <- TRUE or install them manually."
    )
  }
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(GenomeInfoDb)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(AnnotationHub)
  library(ggplot2)
})

suppressPackageStartupMessages({
  library(genome_cfg$txdb_pkg, character.only = TRUE)
  library(genome_cfg$orgdb_pkg, character.only = TRUE)
})

message("Reading BED file: ", bed_file)
if (!file.exists(bed_file)) stop("Input BED file does not exist: ", bed_file)

bed_df <- utils::read.delim(
  bed_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(bed_df) < 3) stop("BED file must contain at least 3 columns: chrom, start, end")

orig_colnames <- paste0("bed_col", seq_len(ncol(bed_df)))
colnames(bed_df) <- orig_colnames
colnames(bed_df)[1:3] <- c("chrom", "start0", "end")

if (any(is.na(bed_df$start0)) || any(is.na(bed_df$end))) {
  stop("BED start/end columns contain NA values.")
}
if (any(bed_df$end < bed_df$start0)) {
  stop("Found BED rows where end < start.")
}

regions <- GRanges(
  seqnames = bed_df$chrom,
  ranges = IRanges(start = as.integer(bed_df$start0) + 1L, end = as.integer(bed_df$end))
)

region_ids <- paste0(bed_df$chrom, ":", bed_df$start0, "-", bed_df$end)
region_ids <- make.unique(region_ids)
names(regions) <- region_ids
mcols(regions)$region_id <- region_ids
mcols(regions)$input_order <- seq_along(regions)

# Load selected gene models.
txdb <- getExportedValue(genome_cfg$txdb_pkg, genome_cfg$txdb_obj)
orgdb <- getExportedValue(genome_cfg$orgdb_pkg, genome_cfg$orgdb_obj)

message("Preparing transcript, promoter, gene-body, and TSS annotations...")

tx_by_gene <- transcriptsBy(txdb, by = "gene", use.names = FALSE)
valid_gene_names <- names(tx_by_gene)
valid_gene_names <- valid_gene_names[!is.na(valid_gene_names) & nzchar(valid_gene_names)]
tx_by_gene <- tx_by_gene[valid_gene_names]

transcripts_gr <- unlist(tx_by_gene, use.names = FALSE)
mcols(transcripts_gr)$gene_id <- rep(names(tx_by_gene), elementNROWS(tx_by_gene))

if (is.null(names(transcripts_gr)) || all(!nzchar(names(transcripts_gr)))) {
  tx_names <- mcols(transcripts_gr)$tx_name
  if (!is.null(tx_names)) {
    names(transcripts_gr) <- as.character(tx_names)
  }
}

promoters_gr <- promoters(
  transcripts_gr,
  upstream = PROMOTER_UPSTREAM,
  downstream = PROMOTER_DOWNSTREAM
)
mcols(promoters_gr)$gene_id <- mcols(transcripts_gr)$gene_id

gene_ranges_list <- range(tx_by_gene)
gene_spans <- unlist(gene_ranges_list, use.names = FALSE)
gene_id_per_span <- rep(names(tx_by_gene), elementNROWS(gene_ranges_list))
mcols(gene_spans)$gene_id <- gene_id_per_span

tss_gr <- promoters(transcripts_gr, upstream = 0L, downstream = 1L)
mcols(tss_gr)$gene_id <- mcols(transcripts_gr)$gene_id

all_gene_ids <- unique(as.character(c(
  mcols(transcripts_gr)$gene_id,
  mcols(gene_spans)$gene_id,
  mcols(tss_gr)$gene_id
)))
all_gene_ids <- all_gene_ids[!is.na(all_gene_ids)]

symbol_map <- AnnotationDbi::mapIds(
  orgdb,
  keys = all_gene_ids,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

annotate_primary <- function(gr, gene_spans, promoters_gr) {
  primary_annotation <- rep("intergenic", length(gr))

  gene_hits <- findOverlaps(gr, gene_spans, ignore.strand = TRUE)
  if (length(gene_hits) > 0L) {
    primary_annotation[unique(queryHits(gene_hits))] <- "gene_body"
  }

  promoter_hits <- findOverlaps(gr, promoters_gr, ignore.strand = TRUE)
  if (length(promoter_hits) > 0L) {
    primary_annotation[unique(queryHits(promoter_hits))] <- "promoter"
  }

  primary_annotation
}

primary_annotation <- annotate_primary(regions, gene_spans, promoters_gr)

message("Computing nearest TSS and gene mapping...")
nearest_hits <- distanceToNearest(regions, tss_gr, ignore.strand = TRUE)
closest_tss_dist <- rep(NA_integer_, length(regions))
closest_gene_id <- rep(NA_character_, length(regions))
closest_gene_symbol <- rep(NA_character_, length(regions))
closest_tx_name <- rep(NA_character_, length(regions))

if (length(nearest_hits) > 0L) {
  qh <- queryHits(nearest_hits)
  sh <- subjectHits(nearest_hits)
  closest_tss_dist[qh] <- mcols(nearest_hits)$distance
  closest_gene_id[qh] <- as.character(mcols(tss_gr)$gene_id[sh])
  closest_tx_name[qh] <- as.character(mcols(tss_gr)$tx_name[sh])
  mapped_symbols <- unname(symbol_map[closest_gene_id[qh]])
  mapped_symbols[is.na(mapped_symbols) | mapped_symbols == ""] <- closest_gene_id[qh][is.na(mapped_symbols) | mapped_symbols == ""]
  closest_gene_symbol[qh] <- mapped_symbols
}

message("Retrieving ", genome_cfg$genome, " RepeatMasker annotations from AnnotationHub...")
get_genome_repeats <- function(genome_cfg) {
  ah <- AnnotationHub()
  q <- query(ah, "RepeatMasker")
  meta <- as.data.frame(S4Vectors::mcols(q))
  keep <- rep(TRUE, nrow(meta))

  if ("species" %in% colnames(meta)) {
    keep <- keep & meta$species == genome_cfg$species
  }
  if ("genome" %in% colnames(meta)) {
    keep <- keep & meta$genome == genome_cfg$genome
  }
  if ("dataprovider" %in% colnames(meta)) {
    keep <- keep & meta$dataprovider == "UCSC"
  }
  if ("rdataclass" %in% colnames(meta)) {
    keep <- keep & meta$rdataclass == "GRanges"
  }

  q2 <- q[keep]
  if (length(q2) == 0L) {
    stop(
      "Could not find a UCSC ", genome_cfg$species, " ", genome_cfg$genome,
      " RepeatMasker GRanges resource in AnnotationHub."
    )
  }
  q2[[1]]
}

repeats_gr <- get_genome_repeats(genome_cfg)
repeat_cols <- colnames(S4Vectors::mcols(repeats_gr))
find_first_col <- function(candidates, available) {
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0L) return(NA_character_)
  hit[[1]]
}
rep_name_col <- find_first_col(c("repName", "name"), repeat_cols)
rep_class_col <- find_first_col(c("repClass", "class"), repeat_cols)
rep_family_col <- find_first_col(c("repFamily", "family"), repeat_cols)
rep_id_col <- find_first_col(c("id", "repID", "repeat_id"), repeat_cols)

if (is.na(rep_name_col) || is.na(rep_class_col)) {
  stop(
    "RepeatMasker metadata columns were not found as expected. Available columns: ",
    paste(repeat_cols, collapse = ", ")
  )
}

rep_name <- as.character(S4Vectors::mcols(repeats_gr)[[rep_name_col]])
rep_class <- as.character(S4Vectors::mcols(repeats_gr)[[rep_class_col]])
rep_family <- if (!is.na(rep_family_col)) as.character(S4Vectors::mcols(repeats_gr)[[rep_family_col]]) else rep("", length(repeats_gr))
rep_link_id <- if (!is.na(rep_id_col)) as.character(S4Vectors::mcols(repeats_gr)[[rep_id_col]]) else as.character(seq_along(repeats_gr))
rep_name[is.na(rep_name)] <- ""
rep_class[is.na(rep_class)] <- ""
rep_family[is.na(rep_family)] <- ""
rep_link_id[is.na(rep_link_id) | rep_link_id == ""] <- as.character(which(is.na(rep_link_id) | rep_link_id == ""))

keep_repeat_idx <- !(rep_class %in% EXCLUDED_REPEAT_CLASSES)
message(
  "Filtering RepeatMasker annotations before overlap calculation. Kept ",
  sum(keep_repeat_idx), " of ", length(keep_repeat_idx),
  " entries after excluding classes: ", paste(EXCLUDED_REPEAT_CLASSES, collapse = ", ")
)
repeats_gr <- repeats_gr[keep_repeat_idx]
rep_name <- rep_name[keep_repeat_idx]
rep_class <- rep_class[keep_repeat_idx]
rep_family <- rep_family[keep_repeat_idx]
rep_link_id <- rep_link_id[keep_repeat_idx]

compute_repeat_details <- function(gr, repeats_gr, rep_link_id, rep_name, rep_class, rep_family) {
  repeat_hits <- findOverlaps(gr, repeats_gr, ignore.strand = TRUE)

  repeat_any <- rep(FALSE, length(gr))
  dominant_repeat_name <- rep(NA_character_, length(gr))
  dominant_repeat_class <- rep(NA_character_, length(gr))
  dominant_repeat_family <- rep(NA_character_, length(gr))
  dominant_repeat_id <- rep(NA_character_, length(gr))
  dominant_repeat_overlap_bp <- rep(0L, length(gr))
  dominant_repeat_overlap_fraction <- rep(0, length(gr))

  repeat_detail_df <- data.frame(
    region_id = character(),
    region_index = integer(),
    repeat_link_id = character(),
    repeat_name = character(),
    repeat_class = character(),
    repeat_family = character(),
    overlap_bp = integer(),
    overlap_fraction = numeric(),
    stringsAsFactors = FALSE
  )

  if (length(repeat_hits) > 0L) {
    qh <- queryHits(repeat_hits)
    sh <- subjectHits(repeat_hits)
    ov_width <- width(pintersect(gr[qh], repeats_gr[sh]))

    detail_raw <- data.frame(
      region_id = names(gr)[qh],
      region_index = qh,
      repeat_link_id = rep_link_id[sh],
      repeat_name = rep_name[sh],
      repeat_class = rep_class[sh],
      repeat_family = rep_family[sh],
      overlap_bp = as.integer(ov_width),
      stringsAsFactors = FALSE
    )

    repeat_detail_df <- stats::aggregate(
      overlap_bp ~ region_id + region_index + repeat_link_id + repeat_name + repeat_class + repeat_family,
      data = detail_raw,
      FUN = sum
    )
    repeat_detail_df$overlap_fraction <- repeat_detail_df$overlap_bp / width(gr[repeat_detail_df$region_index])

    split_idx <- split(seq_len(nrow(repeat_detail_df)), repeat_detail_df$region_index)
    for (idx_chr in names(split_idx)) {
      idx <- split_idx[[idx_chr]]
      sub <- repeat_detail_df[idx, , drop = FALSE]
      sub <- sub[order(-sub$overlap_bp, -sub$overlap_fraction, sub$repeat_name, sub$repeat_class), , drop = FALSE]
      best <- sub[1, , drop = FALSE]
      i <- as.integer(idx_chr)
      repeat_any[i] <- TRUE
      dominant_repeat_id[i] <- best$repeat_link_id
      dominant_repeat_name[i] <- best$repeat_name
      dominant_repeat_class[i] <- best$repeat_class
      dominant_repeat_family[i] <- best$repeat_family
      dominant_repeat_overlap_bp[i] <- best$overlap_bp
      dominant_repeat_overlap_fraction[i] <- best$overlap_fraction
    }
  }

  list(
    repeat_any = repeat_any,
    dominant_repeat_name = dominant_repeat_name,
    dominant_repeat_class = dominant_repeat_class,
    dominant_repeat_family = dominant_repeat_family,
    dominant_repeat_id = dominant_repeat_id,
    dominant_repeat_overlap_bp = dominant_repeat_overlap_bp,
    dominant_repeat_overlap_fraction = dominant_repeat_overlap_fraction,
    repeat_detail_df = repeat_detail_df
  )
}

message("Computing repeat overlaps...")
repeat_res <- compute_repeat_details(regions, repeats_gr, rep_link_id, rep_name, rep_class, rep_family)
repeat_any <- repeat_res$repeat_any
dominant_repeat_name <- repeat_res$dominant_repeat_name
dominant_repeat_class <- repeat_res$dominant_repeat_class
dominant_repeat_family <- repeat_res$dominant_repeat_family
dominant_repeat_id <- repeat_res$dominant_repeat_id
dominant_repeat_overlap_bp <- repeat_res$dominant_repeat_overlap_bp
dominant_repeat_overlap_fraction <- repeat_res$dominant_repeat_overlap_fraction
repeat_detail_df <- repeat_res$repeat_detail_df

out_df <- bed_df
out_df$region_id <- names(regions)
out_df$width_bp <- width(regions)
out_df$primary_annotation <- primary_annotation
out_df$repeat_any <- ifelse(repeat_any, "yes", "no")
out_df$repeat_name <- dominant_repeat_name
out_df$repeat_class <- dominant_repeat_class
out_df$repeat_family <- dominant_repeat_family
out_df$repeat_link_id <- dominant_repeat_id
out_df$repeat_overlap_bp <- dominant_repeat_overlap_bp
out_df$repeat_overlap_fraction <- round(dominant_repeat_overlap_fraction, 4)
out_df$closest_gene <- closest_gene_symbol
out_df$closest_gene_entrezid <- closest_gene_id
out_df$closest_transcript <- closest_tx_name
out_df$distance_to_closest_tss_bp <- closest_tss_dist
out_df$repeat_summary_mode_used_for_class_name_plots <- repeat_summary_mode
out_df$genome <- selected_genome
out_df$excluded_repeat_classes <- paste(EXCLUDED_REPEAT_CLASSES, collapse = ",")

get_seqinfo_lengths <- function(txdb, repeats_gr, regions, allowed_seqlevels = NULL) {
  lens <- GenomeInfoDb::seqlengths(txdb)
  if (is.null(lens) || all(is.na(lens))) {
    lens <- GenomeInfoDb::seqlengths(repeats_gr)
  }
  if (is.null(lens) || all(is.na(lens))) {
    lens <- GenomeInfoDb::seqlengths(regions)
  }
  lens <- lens[!is.na(lens)]
  if (!is.null(allowed_seqlevels)) {
    lens <- lens[names(lens) %in% allowed_seqlevels]
  }
  lens
}

sample_background_regions <- function(observed_regions, seq_lengths, iter_index) {
  chroms <- as.character(seqnames(observed_regions))
  widths_bp <- width(observed_regions)
  sampled <- observed_regions

  for (chr in unique(chroms)) {
    idx <- which(chroms == chr)
    if (!(chr %in% names(seq_lengths))) {
      stop("Chromosome ", chr, " is not present in available sequence lengths for background sampling.")
    }
    chr_len <- as.integer(seq_lengths[[chr]])
    region_widths <- widths_bp[idx]
    if (any(region_widths > chr_len)) {
      stop("At least one region is longer than chromosome ", chr, " and cannot be sampled for background.")
    }
    max_starts <- chr_len - region_widths + 1L
    starts <- vapply(max_starts, function(ms) sample.int(ms, 1L), integer(1))
    ends <- starts + region_widths - 1L
    ranges(sampled)[idx] <- IRanges(start = starts, end = ends)
  }

  bg_ids <- paste0("bg", iter_index, "_", chroms, ":", start(sampled) - 1L, "-", end(sampled))
  names(sampled) <- bg_ids
  mcols(sampled)$region_id <- bg_ids
  mcols(sampled)$input_order <- seq_along(sampled)
  sampled
}

feature_levels <- c("promoter", "gene_body", "intergenic")
repeat_any_levels <- c("yes", "no")

count_levels <- function(vals, levels) {
  vals <- as.character(vals)
  fac <- factor(vals, levels = levels)
  as.integer(table(fac))
}

OTHER_UNKNOWN_CLASSES <- c("Other", "Unknown")
OTHER_UNKNOWN_LABEL <- "Other/Unknown"

collapse_other_unknown <- function(class_vals) {
  class_vals[class_vals %in% OTHER_UNKNOWN_CLASSES] <- OTHER_UNKNOWN_LABEL
  class_vals
}

extract_repeat_plot_data <- function(out_df, repeat_detail_df, mode, top_n_repeat_names) {
  repeat_regions_df <- out_df[out_df$repeat_any == "yes", , drop = FALSE]

  if (mode == "dominant") {
    class_vals <- repeat_regions_df$repeat_class
    class_vals <- class_vals[!is.na(class_vals) & class_vals != ""]
    name_vals <- repeat_regions_df$repeat_name
    name_vals <- name_vals[!is.na(name_vals) & name_vals != ""]
  } else {
    class_vals <- repeat_detail_df$repeat_class
    class_vals <- class_vals[!is.na(class_vals) & class_vals != ""]
    name_vals <- repeat_detail_df$repeat_name
    name_vals <- name_vals[!is.na(name_vals) & name_vals != ""]
  }
  class_vals <- collapse_other_unknown(class_vals)

  class_counts <- sort(table(class_vals), decreasing = TRUE)
  name_counts <- sort(table(name_vals), decreasing = TRUE)
  name_total_all <- sum(name_counts)  # total across ALL detected names, before trimming to top N
  if (length(name_counts) > top_n_repeat_names) {
    name_counts <- name_counts[seq_len(top_n_repeat_names)]
  }

  list(
    class_counts = class_counts,
    name_counts = name_counts,
    name_total_all = name_total_all
  )
}

message("Sampling genomic background with ", background_iterations, " matched randomization(s)...")
seq_lengths <- get_seqinfo_lengths(txdb, repeats_gr, regions, allowed_seqlevels = STANDARD_CHROMS)
input_chroms <- unique(as.character(seqnames(regions)))
missing_chroms <- setdiff(input_chroms, names(seq_lengths))
if (length(missing_chroms) > 0L) {
  stop(
    "Some input chromosomes are unavailable for background sampling for the selected genome: ",
    paste(missing_chroms, collapse = ", "),
    ". The script currently supports standard mm10 chromosomes when estimating background."
  )
}

feature_bg_mat <- matrix(0L, nrow = background_iterations, ncol = length(feature_levels), dimnames = list(NULL, feature_levels))
repeat_any_bg_mat <- matrix(0L, nrow = background_iterations, ncol = length(repeat_any_levels), dimnames = list(NULL, repeat_any_levels))
repeat_class_bg_list <- vector("list", background_iterations)
repeat_name_bg_list <- vector("list", background_iterations)
repeat_name_bg_total_list <- vector("list", background_iterations)

for (iter in seq_len(background_iterations)) {
  bg_regions <- sample_background_regions(regions, seq_lengths, iter)
  bg_primary <- annotate_primary(bg_regions, gene_spans, promoters_gr)
  bg_repeat_res <- compute_repeat_details(bg_regions, repeats_gr, rep_link_id, rep_name, rep_class, rep_family)
  bg_out_df <- data.frame(
    primary_annotation = bg_primary,
    repeat_any = ifelse(bg_repeat_res$repeat_any, "yes", "no"),
    repeat_class = bg_repeat_res$dominant_repeat_class,
    repeat_name = bg_repeat_res$dominant_repeat_name,
    stringsAsFactors = FALSE
  )

  feature_bg_mat[iter, ] <- count_levels(bg_out_df$primary_annotation, feature_levels)
  repeat_any_bg_mat[iter, ] <- count_levels(bg_out_df$repeat_any, repeat_any_levels)

  bg_repeat_plot_data <- extract_repeat_plot_data(
    bg_out_df,
    bg_repeat_res$repeat_detail_df,
    mode = repeat_summary_mode,
    top_n_repeat_names = top_n_repeat_names
  )
  repeat_class_bg_list[[iter]] <- bg_repeat_plot_data$class_counts
  repeat_name_bg_list[[iter]] <- bg_repeat_plot_data$name_counts
  repeat_name_bg_total_list[[iter]] <- bg_repeat_plot_data$name_total_all
}

feature_observed_counts <- count_levels(out_df$primary_annotation, feature_levels)
repeat_any_observed_counts <- count_levels(out_df$repeat_any, repeat_any_levels)
feature_background_mean <- colMeans(feature_bg_mat)
repeat_any_background_mean <- colMeans(repeat_any_bg_mat)

feature_plot_df <- data.frame(
  category = rep(feature_levels, 2L),
  source = rep(c("Observed", "Background"), each = length(feature_levels)),
  n = c(feature_observed_counts, feature_background_mean),
  total = length(regions),
  stringsAsFactors = FALSE
)
feature_plot_df$percent <- feature_plot_df$n / feature_plot_df$total

repeat_assoc_plot_df <- data.frame(
  category = rep(repeat_any_levels, 2L),
  source = rep(c("Observed", "Background"), each = length(repeat_any_levels)),
  n = c(repeat_any_observed_counts, repeat_any_background_mean),
  total = length(regions),
  stringsAsFactors = FALSE
)
repeat_assoc_plot_df$percent <- repeat_assoc_plot_df$n / repeat_assoc_plot_df$total

observed_repeat_plot_data <- extract_repeat_plot_data(out_df, repeat_detail_df, repeat_summary_mode, top_n_repeat_names)
observed_class_counts <- observed_repeat_plot_data$class_counts
observed_name_counts <- observed_repeat_plot_data$name_counts

mean_named_table <- function(count_list, categories) {
  out <- setNames(numeric(length(categories)), categories)
  if (length(categories) == 0L) return(out)
  for (x in count_list) {
    if (length(x) == 0L) next
    nm <- intersect(names(x), categories)
    out[nm] <- out[nm] + as.numeric(x[nm])
  }
  out / length(count_list)
}

class_categories <- union(names(observed_class_counts), unique(unlist(lapply(repeat_class_bg_list, names))))
class_categories <- class_categories[nzchar(class_categories)]
class_categories <- class_categories[order(-mean_named_table(repeat_class_bg_list, class_categories), class_categories)]
observed_class_counts_full <- setNames(numeric(length(class_categories)), class_categories)
if (length(observed_class_counts) > 0L) observed_class_counts_full[names(observed_class_counts)] <- as.numeric(observed_class_counts)
background_class_counts_full <- mean_named_table(repeat_class_bg_list, class_categories)

# IMPORTANT: define Top-N repeat names from the observed data, so the observed panel
# matches the original non-background script when repeat_summary_mode == "dominant".
# Then project the background onto exactly the same observed categories.
name_categories <- names(observed_name_counts)
name_categories <- name_categories[nzchar(name_categories)]
if (length(name_categories) > top_n_repeat_names) {
  name_categories <- name_categories[seq_len(top_n_repeat_names)]
}

# If the observed set has fewer than top_n categories, optionally append the most common
# background-only categories so the panel can still contain up to top_n rows.
if (length(name_categories) < top_n_repeat_names) {
  bg_name_pool <- unique(unlist(lapply(repeat_name_bg_list, names)))
  bg_name_pool <- bg_name_pool[nzchar(bg_name_pool)]
  bg_name_pool <- setdiff(bg_name_pool, name_categories)
  if (length(bg_name_pool) > 0L) {
    bg_name_pool <- bg_name_pool[order(-mean_named_table(repeat_name_bg_list, bg_name_pool), bg_name_pool)]
    n_to_add <- min(top_n_repeat_names - length(name_categories), length(bg_name_pool))
    name_categories <- c(name_categories, bg_name_pool[seq_len(n_to_add)])
  }
}

observed_name_counts_full <- setNames(numeric(length(name_categories)), name_categories)
if (length(observed_name_counts) > 0L) {
  nm <- intersect(names(observed_name_counts), name_categories)
  observed_name_counts_full[nm] <- as.numeric(observed_name_counts[nm])
}
background_name_counts_full <- mean_named_table(repeat_name_bg_list, name_categories)

class_total_observed <- sum(observed_class_counts_full)
class_total_background <- sum(background_class_counts_full)
# Use totals across ALL detected names (not just top N) so percentages are
# relative to the full detected repeat-name universe in each set.
name_total_observed <- observed_repeat_plot_data$name_total_all
name_total_background <- mean(unlist(repeat_name_bg_total_list))
if (is.na(name_total_background)) name_total_background <- 0

repeat_class_plot_df <- data.frame(
  category = rep(class_categories, 2L),
  source = rep(c("Observed", "Background"), each = length(class_categories)),
  n = c(as.numeric(observed_class_counts_full), as.numeric(background_class_counts_full)),
  total = c(rep(class_total_observed, length(class_categories)), rep(class_total_background, length(class_categories))),
  stringsAsFactors = FALSE
)
repeat_class_plot_df$percent <- ifelse(repeat_class_plot_df$total > 0, repeat_class_plot_df$n / repeat_class_plot_df$total, 0)

repeat_name_plot_df <- data.frame(
  category = rep(name_categories, 2L),
  source = rep(c("Observed", "Background"), each = length(name_categories)),
  n = c(as.numeric(observed_name_counts_full), as.numeric(background_name_counts_full)),
  total = c(rep(name_total_observed, length(name_categories)), rep(name_total_background, length(name_categories))),
  stringsAsFactors = FALSE
)
repeat_name_plot_df$percent <- ifelse(repeat_name_plot_df$total > 0, repeat_name_plot_df$n / repeat_name_plot_df$total, 0)

background_summary_df <- rbind(
  data.frame(
    plot = "feature_annotation",
    category = feature_levels,
    observed_n = feature_observed_counts,
    observed_percent = feature_observed_counts / length(regions),
    background_mean_n = feature_background_mean,
    background_mean_percent = feature_background_mean / length(regions),
    repeat_summary_mode = repeat_summary_mode,
    background_iterations = background_iterations,
    stringsAsFactors = FALSE
  ),
  data.frame(
    plot = "repeat_association",
    category = repeat_any_levels,
    observed_n = repeat_any_observed_counts,
    observed_percent = repeat_any_observed_counts / length(regions),
    background_mean_n = repeat_any_background_mean,
    background_mean_percent = repeat_any_background_mean / length(regions),
    repeat_summary_mode = repeat_summary_mode,
    background_iterations = background_iterations,
    stringsAsFactors = FALSE
  ),
  data.frame(
    plot = "repeat_class",
    category = class_categories,
    observed_n = as.numeric(observed_class_counts_full),
    observed_percent = if (class_total_observed > 0) as.numeric(observed_class_counts_full) / class_total_observed else 0,
    background_mean_n = as.numeric(background_class_counts_full),
    background_mean_percent = if (class_total_background > 0) as.numeric(background_class_counts_full) / class_total_background else 0,
    repeat_summary_mode = repeat_summary_mode,
    background_iterations = background_iterations,
    stringsAsFactors = FALSE
  ),
  data.frame(
    plot = "repeat_name_topN",
    category = name_categories,
    observed_n = as.numeric(observed_name_counts_full),
    observed_percent = if (name_total_observed > 0) as.numeric(observed_name_counts_full) / name_total_observed else 0,
    background_mean_n = as.numeric(background_name_counts_full),
    background_mean_percent = if (name_total_background > 0) as.numeric(background_name_counts_full) / name_total_background else 0,
    repeat_summary_mode = repeat_summary_mode,
    background_iterations = background_iterations,
    stringsAsFactors = FALSE
  )
)
background_summary_df$genome <- selected_genome
background_summary_df$excluded_repeat_classes <- paste(EXCLUDED_REPEAT_CLASSES, collapse = ",")

annotated_file <- paste0(out_prefix, ".annotated.tsv")
repeat_detail_file <- paste0(out_prefix, ".repeat_detail.tsv")
background_summary_file <- paste0(out_prefix, ".background_summary.tsv")
utils::write.table(out_df, file = annotated_file, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(repeat_detail_df, file = repeat_detail_file, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(background_summary_df, file = background_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote: ", annotated_file)
message("Wrote: ", repeat_detail_file)
message("Wrote: ", background_summary_file)

# Plots are sized as small paper-figure panels (a few cm on a side), not full-page
# figures: minimal text, no titles (captions belong in the manuscript legend), and a
# footprint that scales with the number of bars rather than a single fixed canvas.
panel_size_cm <- function(n_categories, flip) {
  if (flip) {
    list(width = 4.2, height = min(6, max(3.2, 0.42 * n_categories + 1.5)))
  } else {
    list(width = min(6, max(3.4, 0.75 * n_categories + 2.0)), height = 4.0)
  }
}

# At this panel size there isn't room for a "n=51.0"-style label without bars
# overlapping their neighbor, so bar labels are a bare rounded integer.
make_percent_compare_plot <- function(plot_df, x_label, y_label, out_file, flip = FALSE,
                                       width_cm = 4, height_cm = 4) {
  plot_df$source <- factor(plot_df$source, levels = c("Observed", "Background"))
  plot_df$category <- factor(plot_df$category, levels = unique(plot_df$category))
  ymax <- max(plot_df$percent, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 0.05

  p <- ggplot(plot_df, aes(x = category, y = percent, fill = source)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(
      aes(label = sprintf("%.0f", n)),
      position = position_dodge(width = 0.8),
      hjust = if (flip) -0.15 else 0.5,
      vjust = if (flip) 0.5 else -0.3,
      size = 1.6
    ) +
    scale_y_continuous(
      limits = c(0, ymax * 1.20),
      labels = function(x) sprintf("%.0f%%", x * 100)
    ) +
    labs(x = x_label, y = y_label, fill = NULL) +
    theme_bw(base_size = 6) +
    theme(
      plot.margin = margin(1, 2, 1, 1, unit = "mm"),
      legend.position = "bottom",
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 5),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0.5, "mm"),
      axis.text = element_text(size = 5),
      axis.text.x = if (flip) element_text(size = 5) else element_text(size = 5, angle = 40, hjust = 1, vjust = 1),
      axis.title = element_text(size = 6)
    )

  if (flip) {
    p <- p + coord_flip()
  }

  ggsave(paste0(out_file, ".png"), p, width = width_cm, height = height_cm, units = "cm", dpi = 600)
  ggsave(paste0(out_file, ".pdf"), p, width = width_cm, height = height_cm, units = "cm")
}

feature_plot <- paste0(out_prefix, ".feature_annotation_percent")
repeat_class_plot <- paste0(out_prefix, ".repeat_class_frequency")
repeat_name_plot <- paste0(out_prefix, ".repeat_name_top", top_n_repeat_names, "_frequency")

no_data_plot <- function(out_file, width_cm = 4, height_cm = 3.6) {
  png(paste0(out_file, ".png"), width = width_cm, height = height_cm, units = "cm", res = 600)
  plot.new()
  text(0.5, 0.5, "No repeat-overlapping\nregions found", cex = 0.5)
  dev.off()
  pdf(paste0(out_file, ".pdf"), width = width_cm / 2.54, height = height_cm / 2.54)
  plot.new()
  text(0.5, 0.5, "No repeat-overlapping\nregions found", cex = 0.5)
  dev.off()
}

message("Generating plots...")
feature_size <- panel_size_cm(length(unique(feature_plot_df$category)), flip = FALSE)
make_percent_compare_plot(
  feature_plot_df,
  x_label = "Region annotation",
  y_label = "% of regions",
  out_file = feature_plot,
  flip = FALSE,
  width_cm = feature_size$width,
  height_cm = feature_size$height
)

# Repeat association (percentage of regions overlapping any repeat) is not plotted:
# repeats cover most of the genome, so this comparison to background is uninformative.
# repeat_assoc_plot_df is still written to background_summary.tsv.

if (nrow(repeat_class_plot_df) > 0L && any(repeat_class_plot_df$n > 0)) {
  class_size <- panel_size_cm(length(unique(repeat_class_plot_df$category)), flip = FALSE)
  make_percent_compare_plot(
    repeat_class_plot_df,
    x_label = "Repeat class",
    y_label = if (repeat_summary_mode == "dominant") "Fraction of regions" else "Fraction of overlaps",
    out_file = repeat_class_plot,
    flip = FALSE,
    width_cm = class_size$width,
    height_cm = class_size$height
  )
} else {
  no_data_plot(repeat_class_plot)
}

if (nrow(repeat_name_plot_df) > 0L && any(repeat_name_plot_df$n > 0)) {
  name_size <- panel_size_cm(length(unique(repeat_name_plot_df$category)), flip = TRUE)
  make_percent_compare_plot(
    repeat_name_plot_df,
    x_label = "Repeat name",
    y_label = if (repeat_summary_mode == "dominant") "Fraction of regions" else "Fraction of overlaps",
    out_file = repeat_name_plot,
    flip = TRUE,
    width_cm = name_size$width,
    height_cm = name_size$height
  )
} else {
  no_data_plot(repeat_name_plot)
}

message("Done.")
message("Primary output table: ", annotated_file)
message("Detailed repeat overlap table: ", repeat_detail_file)
message("Background summary table: ", background_summary_file)
message("Repeat summary mode for class/name plots: ", repeat_summary_mode)
message("Background iterations: ", background_iterations)
message("Plots:")
message("  ", feature_plot)
message("  ", repeat_class_plot)
message("  ", repeat_name_plot)
