#!/usr/bin/env Rscript

# ==============================================================================
# 04 -- Merge per-pool STARsolo matrices into unified count tables
#       (PT x B73 cold BRB-seq, pools 1-4)
#
# STARsolo produces one Solo.out/Gene/raw/ directory per pool with:
#   umiDedup-1MM_Directional.mtx    UMI-collapsed counts (canonical)
#   umiDedup-NoDedup.mtx            raw read counts
#   features.tsv                    gene IDs (same across pools; sanity-checked)
#   barcodes.tsv                    the 96 whitelist barcodes (--soloCellFilter None)
#
# The SAME 96 barcodes are reused in every pool, so barcode -> sample_id must be
# mapped per pool (pool_N_barcode_map.tsv) BEFORE pools are combined.
#
# This script:
#   - reads all four pools' matrices
#   - maps barcode -> sample_id per pool, and checks the maps match the data
#   - checks gene sets match and sample IDs are unique across pools
#   - column-binds pools into single gene x sample matrices
#   - writes:
#       data/processed/Zea_mays_counts.txt        UMI-collapsed counts (canonical)
#       data/processed/Zea_mays_counts_raw.txt    no-dedup read counts (comparison)
#       data/processed/counts_umi.rds / counts_raw.rds   sparse matrices for R
#       data/processed/sample_metadata.csv        sample_id, pool, barcode
#       data/processed/sample_qc.csv              per-well UMIs, genes, flags
#       data/processed/starsolo_pool_qc.csv       per-pool trimming + STAR QC
#
# Run from the project root on the HPC:
#   cd /rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq
#   Rscript scripts/04_merge_pools.R
# Override locations with STARSOLO_ROOT / TRIM_ROOT / MAP_DIR env vars if needed.
# ==============================================================================

suppressPackageStartupMessages(library(Matrix))

project_dir <- Sys.getenv("PROJECT_DIR",
                          unset = "/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq")
solo_root <- Sys.getenv("STARSOLO_ROOT", unset = file.path(project_dir, "starsolo"))
trim_root <- Sys.getenv("TRIM_ROOT",     unset = file.path(project_dir, "trimmed"))
map_dir   <- Sys.getenv("MAP_DIR",       unset = file.path(project_dir, "data", "starsolo"))
out_dir   <- file.path(project_dir, "data", "processed")
pools     <- 1:4
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("STARsolo root: ", solo_root, "\n",
    "Trim logs:     ", trim_root, "\n",
    "Barcode maps:  ", map_dir,   "\n",
    "Output:        ", out_dir,   "\n", sep = "")

# --- read barcode maps once --------------------------------------------------
# Expected columns: sample_id, barcode (tab-separated, with header)
read_map <- function(p) {
  path <- file.path(map_dir, sprintf("pool_%d_barcode_map.tsv", p))
  if (!file.exists(path)) stop("Missing barcode map: ", path)
  m <- read.delim(path, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (!all(c("sample_id", "barcode") %in% names(m)))
    stop(path, " needs columns 'sample_id' and 'barcode'; found: ",
         paste(names(m), collapse = ", "))
  if (anyDuplicated(m$barcode))   stop(path, ": duplicated barcodes")
  if (anyDuplicated(m$sample_id)) stop(path, ": duplicated sample_ids")
  m$pool <- p
  m
}
maps <- lapply(pools, read_map)
names(maps) <- pools

# Sample IDs must be unique across ALL pools, or cbind gives duplicate columns
all_sids <- unlist(lapply(maps, `[[`, "sample_id"))
dup_sids <- unique(all_sids[duplicated(all_sids)])
if (length(dup_sids) > 0) {
  stop("sample_id values repeat across pools (e.g. parental checks?): ",
       paste(head(dup_sids, 10), collapse = ", "),
       "\nMake them unique in the barcode maps (e.g. add a _P<pool> suffix).")
}

# --- read one pool's matrix ---------------------------------------------------
read_pool <- function(p, matrix_name) {
  raw_dir  <- file.path(solo_root, sprintf("pool_%d", p), "Solo.out", "Gene", "raw")
  mtx_path <- file.path(raw_dir, matrix_name)
  bc_path  <- file.path(raw_dir, "barcodes.tsv")
  ft_path  <- file.path(raw_dir, "features.tsv")
  for (f in c(mtx_path, bc_path, ft_path)) if (!file.exists(f)) stop("Missing input: ", f)

  mat <- as(readMM(mtx_path), "CsparseMatrix")          # genes x barcodes
  bcs <- readLines(bc_path)
  fts <- read.delim(ft_path, header = FALSE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "")[[1]] # gene_id per row
  if (nrow(mat) != length(fts) || ncol(mat) != length(bcs))
    stop("Pool ", p, ": matrix dims don't match features/barcodes")
  rownames(mat) <- fts
  colnames(mat) <- bcs

  map <- maps[[as.character(p)]]
  missing_in_data <- setdiff(map$barcode, bcs)
  extra_in_data   <- setdiff(bcs, map$barcode)
  if (length(missing_in_data) > 0)
    warning(sprintf("Pool %d: %d mapped barcode(s) absent from STARsolo output: %s",
                    p, length(missing_in_data), paste(head(missing_in_data, 5), collapse = ", ")))
  if (length(extra_in_data) > 0)
    cat(sprintf("Pool %d: dropping %d barcode(s) not in the map\n", p, length(extra_in_data)))

  keep_bc <- intersect(bcs, map$barcode)
  mat <- mat[, keep_bc, drop = FALSE]
  colnames(mat) <- map$sample_id[match(colnames(mat), map$barcode)]
  cat(sprintf("Pool %d [%s]: %d genes x %d samples\n",
              p, matrix_name, nrow(mat), ncol(mat)))
  mat
}

merge_pools <- function(matrix_name, txt_path, rds_path) {
  cat("\n=== Merging pools for ", matrix_name, " ===\n", sep = "")
  mats <- lapply(pools, read_pool, matrix_name = matrix_name)

  gene_ok <- all(vapply(mats[-1], function(m)
    identical(rownames(m), rownames(mats[[1]])), logical(1)))
  if (!gene_ok) stop("Gene ID sets differ across pools -- check all pools used the same STAR index/GTF")

  combined <- do.call(cbind, mats)
  cat(sprintf("Combined: %d genes x %d samples\n", nrow(combined), ncol(combined)))

  saveRDS(combined, rds_path)                     # sparse; use this in R
  df <- as.data.frame(as.matrix(combined))        # TSV in legacy Zea_mays_counts.txt shape
  write.table(df, txt_path, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  cat("Wrote ", txt_path, " (", format(file.size(txt_path) / 1e6, digits = 3),
      " MB) and ", basename(rds_path), "\n", sep = "")
  invisible(combined)
}

umi <- merge_pools("umiDedup-1MM_Directional.mtx",
                   file.path(out_dir, "Zea_mays_counts.txt"),
                   file.path(out_dir, "counts_umi.rds"))
raw <- merge_pools("umiDedup-NoDedup.mtx",
                   file.path(out_dir, "Zea_mays_counts_raw.txt"),
                   file.path(out_dir, "counts_raw.rds"))
stopifnot(identical(colnames(umi), colnames(raw)))

# --- sample metadata -----------------------------------------------------------
meta <- do.call(rbind, lapply(maps, function(m) m[, c("sample_id", "pool", "barcode")]))
meta <- meta[match(colnames(umi), meta$sample_id), ]
write.csv(meta, file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)

# --- per-sample (per-well) QC --------------------------------------------------
umi_lib <- colSums(umi)
raw_lib <- colSums(raw)
genes_detected <- colSums(umi > 0)
dedup_rate <- ifelse(raw_lib > 0, 1 - umi_lib / raw_lib, NA_real_)

sample_qc <- data.frame(meta,
                        umi_total      = umi_lib,
                        reads_total    = raw_lib,
                        genes_detected = genes_detected,
                        dedup_rate     = round(dedup_rate, 4),
                        row.names = NULL)
# flag wells far below their own pool's median (failed or empty wells)
pool_med <- tapply(sample_qc$umi_total, sample_qc$pool, median)
sample_qc$frac_of_pool_median <- round(sample_qc$umi_total / pool_med[as.character(sample_qc$pool)], 3)
sample_qc$low_flag <- sample_qc$frac_of_pool_median < 0.10
write.csv(sample_qc, file.path(out_dir, "sample_qc.csv"), row.names = FALSE)

cat(sprintf("\nPer-sample PCR-duplicate rate: median %.1f%%, IQR [%.1f%%, %.1f%%]\n",
            100 * median(dedup_rate, na.rm = TRUE),
            100 * quantile(dedup_rate, 0.25, na.rm = TRUE),
            100 * quantile(dedup_rate, 0.75, na.rm = TRUE)))
cat("\nPer-pool well summary (UMIs and genes per well):\n")
print(do.call(rbind, lapply(split(sample_qc, sample_qc$pool), function(d) data.frame(
  pool = d$pool[1], wells = nrow(d),
  median_umi = median(d$umi_total), min_umi = min(d$umi_total),
  median_genes = median(d$genes_detected),
  low_wells = sum(d$low_flag)))), row.names = FALSE)
if (any(sample_qc$low_flag)) {
  cat("\nWells below 10% of their pool median:\n")
  print(sample_qc[sample_qc$low_flag, c("sample_id", "pool", "barcode", "umi_total", "genes_detected")],
        row.names = FALSE)
}

# --- per-pool QC (trimming + STAR + STARsolo) ---------------------------------
# helper: first number on the first line matching a regex
grep1 <- function(path, pattern) {
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  hit <- grep(pattern, lines, perl = TRUE, value = TRUE)
  if (length(hit) == 0) return(NA_real_)
  m <- regmatches(hit[1], regexpr(pattern, hit[1], perl = TRUE))
  as.numeric(sub(".*?([0-9]+(\\.[0-9]+)?)$", "\\1", m, perl = TRUE))
}
summary_val <- function(path, key) {               # Solo.out/Gene/Summary.csv
  if (!file.exists(path)) return(NA_real_)
  s <- read.csv(path, header = FALSE, stringsAsFactors = FALSE)
  v <- s$V2[s$V1 == key]
  if (length(v) == 0) NA_real_ else suppressWarnings(as.numeric(v))
}

qc <- do.call(rbind, lapply(pools, function(p) {
  trim_log <- file.path(trim_root, sprintf("pool_%d_trim.log", p))
  star_log <- file.path(solo_root, sprintf("pool_%d", p), "Log.final.out")
  solo_sum <- file.path(solo_root, sprintf("pool_%d", p), "Solo.out", "Gene", "Summary.csv")

  trim_in   <- grep1(trim_log, "Input Read Pairs:\\s*[0-9]+")
  trim_surv <- grep1(trim_log, "Both Surviving:\\s*[0-9]+")
  d <- sample_qc[sample_qc$pool == p, ]

  data.frame(
    pool                  = p,
    trim_input_pairs      = trim_in,
    trim_survived_pairs   = trim_surv,
    trim_survival_pct     = if (!is.na(trim_in) && trim_in > 0) 100 * trim_surv / trim_in else NA_real_,
    starsolo_input_reads  = grep1(star_log, "Number of input reads\\s*\\|\\s*[0-9]+"),
    uniquely_mapped_pct   = grep1(star_log, "Uniquely mapped reads %\\s*\\|\\s*[0-9.]+"),
    too_many_loci_pct     = grep1(star_log, "% of reads mapped to too many loci\\s*\\|\\s*[0-9.]+"),
    too_short_pct         = grep1(star_log, "% of reads unmapped: too short\\s*\\|\\s*[0-9.]+"),
    valid_barcodes        = summary_val(solo_sum, "Reads With Valid Barcodes"),
    mapped_to_gene        = summary_val(solo_sum, "Reads Mapped to Gene: Unique Gene"),
    saturation            = summary_val(solo_sum, "Sequencing Saturation"),
    umi_total             = sum(d$umi_total),
    reads_in_matrix       = sum(d$reads_total),
    dedup_rate_pct        = if (sum(d$reads_total) > 0) 100 * (1 - sum(d$umi_total) / sum(d$reads_total)) else NA_real_,
    median_umi_per_well   = median(d$umi_total),
    median_genes_per_well = median(d$genes_detected),
    low_wells             = sum(d$low_flag)
  )
}))
qc_path <- file.path(out_dir, "starsolo_pool_qc.csv")
write.csv(qc, qc_path, row.names = FALSE)
cat("\nPer-pool QC summary:\n")
print(qc, row.names = FALSE, digits = 4)
cat("Wrote ", qc_path, "\n\nDone.\n", sep = "")
