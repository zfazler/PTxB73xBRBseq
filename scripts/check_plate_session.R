
# ============================================================
# Check plate-session confounding and model estimability
# Project: PT x B73 BRB-Seq
# Run from the project root
# ============================================================

suppressPackageStartupMessages({
  library(limma)
})

# 1. Load metadata --------------------------------------------

meta <- read.csv(
  "data/metadata.csv",
  stringsAsFactors = FALSE,
  na.strings = c("NA", "NA ", "")
)

# Trim whitespace from character columns
meta[] <- lapply(
  meta,
  function(x) if (is.character(x)) trimws(x) else x
)

# Check required columns
required <- c(
  "sample_id", "plate", "run", "experimental_group",
  "genotype_class", "row", "column"
)

missing_cols <- setdiff(required, names(meta))

if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# 2. Apply the same sample exclusions as the DE analysis ------
# If sample_qc.csv exists, remove samples with low pool fraction.
# If counts_umi.rds exists, retain samples found in the count matrix.

if (file.exists("data/processed/sample_qc.csv")) {
  qc <- read.csv("data/processed/sample_qc.csv")

  if (all(c("sample_id", "frac_of_pool_median") %in% names(qc))) {
    low_quality <- qc$sample_id[
      !is.na(qc$frac_of_pool_median) &
        qc$frac_of_pool_median < 0.10
    ]

    meta <- meta[!meta$sample_id %in% low_quality, , drop = FALSE]
  }
}

if (file.exists("data/processed/counts_umi.rds")) {
  cnt <- readRDS("data/processed/counts_umi.rds")
  meta <- meta[meta$sample_id %in% colnames(cnt), , drop = FALSE]
}

meta <- meta[
  !is.na(meta$genotype_class) &
    !is.na(meta$plate) &
    !is.na(meta$run) &
    !is.na(meta$experimental_group) &
    !is.na(meta$row) &
    !is.na(meta$column),
  ,
  drop = FALSE
]

cat("\nNumber of samples analyzed:", nrow(meta), "\n")

# 3. Define model variables -----------------------------------

meta$plate <- factor(meta$plate)
meta$genotype_class <- factor(meta$genotype_class)

# Match the session definition from your DE script
meta$session <- factor(
  paste0("run", meta$run, "_grp", meta$experimental_group)
)

# IMPORTANT:
# This treats row position as a categorical effect.
# If your original DE model intentionally uses a numeric row trend,
# change this to: meta$door <- match(toupper(meta$row), LETTERS[1:8])
meta$door <- factor(match(toupper(meta$row), LETTERS[1:8]))

# Smooth trend across columns
ns_d <- splines::ns(as.numeric(meta$column), df = 3)
colnames(ns_d) <- paste0("depth_s", 1:3)
meta <- cbind(meta, ns_d)

# 4. Examine plate-session overlap ---------------------------

tab <- table(meta$plate, meta$session)

cat("\n--- Plate x session sample counts ---\n")
print(tab)

cat("\n--- Proportion of each plate in each session ---\n")
print(round(prop.table(tab, margin = 1), 2))

cat("\n--- Proportion of each session found on each plate ---\n")
print(round(prop.table(tab, margin = 2), 2))

# Number of plates and sessions represented
cat("\nNumber of plates:", nlevels(meta$plate), "\n")
cat("Number of sessions:", nlevels(meta$session), "\n")

# 5. Check design matrix rank ---------------------------------

formulas <- list(
  A_plate = ~ genotype_class + plate,

  B_plate_session =
    ~ genotype_class + plate + session,

  C_plate_session_position =
    ~ genotype_class + plate + session +
      door + depth_s1 + depth_s2 + depth_s3,

  D_session_position =
    ~ genotype_class + session +
      door + depth_s1 + depth_s2 + depth_s3
)

results <- list()

cat("\n--- Model rank and residual degrees of freedom ---\n")

for (nm in names(formulas)) {
  X <- model.matrix(formulas[[nm]], data = meta)

  rank_x <- qr(X)$rank
  residual_df <- nrow(X) - rank_x
  full_rank <- rank_x == ncol(X)

  results[[nm]] <- list(
    design = X,
    rank = rank_x,
    residual_df = residual_df,
    full_rank = full_rank
  )

  cat("\nModel:", nm, "\n")
  cat("  Samples:", nrow(X), "\n")
  cat("  Design columns:", ncol(X), "\n")
  cat("  Design rank:", rank_x, "\n")
  cat("  Full rank:", full_rank, "\n")
  cat("  Residual df:", residual_df, "\n")

  if (!full_rank) {
    cat("  Non-estimable coefficients:\n")
    print(limma::nonEstimable(X))
  }
}

# 6. Save outputs ---------------------------------------------

write.csv(
  as.data.frame.matrix(tab),
  "data/processed/plate_session_counts.csv"
)

summary_table <- data.frame(
  model = names(results),
  samples = vapply(results, function(x) nrow(x$design), numeric(1)),
  columns = vapply(results, function(x) ncol(x$design), numeric(1)),
  rank = vapply(results, function(x) x$rank, numeric(1)),
  full_rank = vapply(results, function(x) x$full_rank, logical(1)),
  residual_df = vapply(results, function(x) x$residual_df, numeric(1))
)

write.csv(
  summary_table,
  "data/processed/model_estimability_summary.csv",
  row.names = FALSE
)

cat("\n--- Summary table ---\n")
print(summary_table)

cat("\nFinished. Outputs saved in data/processed/.\n")
