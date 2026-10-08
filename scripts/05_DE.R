#!/usr/bin/env Rscript
# ==============================================================================
# 05 -- Expression analysis for PT x B73 cold BRB-seq (plates/pools 1-4)
#
# Design (from metadata.csv, 384 well):
#   ril  305 wells / 224 RILs   (77 RILs in 2 wells, 2 in 3 wells)
#   nil   51 wells /  50 NILs   (B73 background + PT segment)
#   b73   27 wells /  13 plants (B73_01..B73_13, most sampled 2-3x across plates)
#   1 empty well (PN4_SID331, plate 4 D07) -> negative control
#   No PT parent wells, so there is no PT-vs-B73 contrast.
#
# What this script does
#   0. Check barcode maps against metadata + whitelist (stops if they disagree)
#   1. Load UMI counts + merge QC + metadata; clean metadata; drop failed wells
#   2. Filter genes, TMM-normalize, logCPM
#   3. Technical QC: PCA, plate effect on B73 checks, replicate concordance
#   4. Class contrasts: RIL vs B73, NIL vs B73, RIL vs NIL
#      (limma-voom; replicate wells of one genotype handled with duplicateCorrelation)
#   5. Per-line screens: each NIL vs B73 and each RIL vs B73 (limma-voom, one joint fit)
#   6. Candidate genes: plots by class (+ allele test if a genotype file exists)
#   7. Export eQTL-ready expression (one value per genotype)
#
# Run from the project root on a compute node:
#   srun --cpus-per-task=4 --mem=32G --time=02:00:00 --pty bash
#   cd /rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq
#   module load R/4.5.0
#   Rscript scripts/05_expression_analysis.R
# Needs R packages: Matrix, edgeR, limma, ggplot2
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(edgeR); library(limma); library(ggplot2)
})

# ---- settings -------------------------------------------------------------------
proc_dir   <- file.path("data", "processed")
counts_rds <- file.path(proc_dir, "counts_umi.rds")      # from 04_merge_pools.R
qc_file    <- file.path(proc_dir, "sample_qc.csv")       # from 04_merge_pools.R
meta_file  <- file.path("data", "metadata.csv")
map_dir    <- file.path("data", "starsolo")              # pool_N_barcode_map.tsv + whitelist
cand_file  <- file.path("data", "candidate_genes.csv")   # OPTIONAL: gene_id, gene_name, category
gtf_path   <- file.path("data", "external", "Zea_mays.gtf")        # OPTIONAL (allele test)
GENO_FILE  <- file.path("data", "genotypes", "ril_genotypes.csv")  # OPTIONAL (allele test)
# GENO_FILE format: marker, chr, pos, then one column per RIL (LANPT16B001 ...)
# coded B73 / PT (or 0 / 2; 1 = het -> treated as missing)

RUN_PER_RIL_SCREEN <- TRUE   # FALSE = screen NILs only
min_frac_of_pool_median <- 0.10   # wells below this fraction of their pool median are dropped
fdr_cut <- 0.05

out_dir <- file.path("output", "expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 0. barcode map check --------------------------------------------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE, na.strings = c("NA", "NA ", ""))
meta[] <- lapply(meta, function(x) if (is.character(x)) trimws(x) else x)   # "nil " -> "nil"
meta$column <- suppressWarnings(as.integer(meta$column))

maps <- do.call(rbind, lapply(1:4, function(p)
  read.delim(file.path(map_dir, sprintf("pool_%d_barcode_map.tsv", p)), stringsAsFactors = FALSE)))
wl   <- readLines(file.path(map_dir, "barcode_whitelist.txt"))
mm   <- merge(meta[!is.na(meta$genotype), c("sample_id", "plate_pos")], maps, by = "sample_id")
per_well <- tapply(mm$barcode, mm$plate_pos, function(b) length(unique(b)))
n_off_wl <- sum(!mm$barcode %in% wl)
cat(sprintf("Barcode maps: %d of %d samples matched; %d wells with >1 barcode across plates; %d barcodes not in whitelist\n",
            nrow(mm), nrow(meta), sum(per_well > 1), n_off_wl))
if (nrow(mm) != sum(!is.na(meta$genotype)) || any(per_well > 1) || n_off_wl > 0)
stop("Barcode maps disagree with metadata/whitelist -- fix before analysis")

# ---- 1. load + clean -------------------------------------------------------------
counts <- readRDS(counts_rds)
qc     <- read.csv(qc_file, stringsAsFactors = FALSE)

# empty well = cross-talk check; near 0 also confirms barcode -> well orientation
for (s in meta$sample_id[is.na(meta$genotype)]) {
  if (s %in% colnames(counts)) {
    pm <- median(qc$umi_total[qc$pool == qc$pool[qc$sample_id == s]])
    cat(sprintf("Empty well %s: %d UMIs (%.1f%% of its pool median) -- expect near 0\n",
                s, round(sum(counts[, s])), 100 * sum(counts[, s]) / pm))
  }
}

missing <- setdiff(meta$sample_id[!is.na(meta$genotype)], colnames(counts))
if (length(missing)) warning(length(missing), " metadata samples not in counts, e.g. ",
                             paste(head(missing, 5), collapse = ", "))

low <- qc$sample_id[qc$frac_of_pool_median < min_frac_of_pool_median]
cat("Dropping", length(low), "low-UMI wells:", paste(low, collapse = ", "), "\n")

meta <- subset(meta, !is.na(genotype) & sample_id %in% colnames(counts) & !sample_id %in% low)
counts <- counts[, meta$sample_id]

meta$plate <- factor(meta$plate)
meta$run   <- factor(meta$run)
meta$experimental_group <- factor(meta$experimental_group)
meta$genotype_class <- factor(meta$genotype_class, levels = c("b73", "nil", "ril"))
# all B73 wells form one "line"; B73 plant IDs (B73_01..) are kept for replicate checks
meta$line <- ifelse(meta$genotype_class == "b73", "B73", meta$genotype)
cat("Samples kept:", ncol(counts), "\n"); print(table(meta$genotype_class, meta$plate))

# ---- 2. filter + normalize -------------------------------------------------------
y <- DGEList(counts = as.matrix(counts), samples = meta)
keep <- suppressWarnings(filterByExpr(y, min.count = 10, min.prop = 0.2))  # ~20% of samples
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
cat("Genes kept after filtering:", nrow(y), "\n")

logcpm <- cpm(y, log = TRUE, prior.count = 2)
saveRDS(logcpm, file.path(out_dir, "logCPM_TMM_filtered.rds"))

cand <- if (file.exists(cand_file)) read.csv(cand_file, stringsAsFactors = FALSE) else NULL
if (!is.null(cand)) cat(sum(cand$gene_id %in% rownames(y)), "of", nrow(cand),
                        "candidate genes pass expression filtering\n")

# ---- 3. technical QC -------------------------------------------------------------
# 3a. PCA on the 2,000 most variable genes
top <- head(order(apply(logcpm, 1, var), decreasing = TRUE), 2000)
pc  <- prcomp(t(logcpm[top, ]), center = TRUE, scale. = FALSE)
ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pca_df <- data.frame(meta, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = plate, shape = genotype_class)) +
  geom_point(alpha = 0.8) + theme_bw() +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
       title = "PCA of logCPM (top 2,000 variable genes)")
ggsave(file.path(out_dir, "pca_plate_class.png"), p_pca, width = 7, height = 5, dpi = 300)
pc_plate_r2 <- sapply(1:5, function(i) summary(lm(pc$x[, i] ~ meta$plate))$r.squared)
cat("R^2 of PC1-5 explained by plate:", round(pc_plate_r2, 3), "\n")

# 3b. plate effect on B73 checks only (genetically identical -> any DE is technical)
yb <- y[, meta$genotype_class == "b73"]
if (nlevels(droplevels(yb$samples$plate)) > 1) {
  db <- model.matrix(~ droplevels(plate), data = yb$samples)
  yb <- estimateDisp(yb, db)
  tb <- glmQLFTest(glmQLFit(yb, db), coef = 2:ncol(db))
  cat("B73 checks: genes with a plate effect (FDR < 0.05):",
      sum(p.adjust(tb$table$PValue, "BH") < fdr_cut), "of", nrow(yb), "\n")
  write.csv(topTags(tb, n = Inf)$table, file.path(out_dir, "b73_plate_effect.csv"))
}

# 3c. replicate concordance: same genotype (or same B73 plant) in 2+ wells
rep_geno <- names(which(table(meta$genotype) > 1))
rep_cor <- do.call(rbind, lapply(rep_geno, function(g) {
  s <- meta$sample_id[meta$genotype == g]
  cm <- cor(logcpm[, s])
  data.frame(genotype = g, class = meta$genotype_class[meta$genotype == g][1],
             n_wells = length(s), mean_r = mean(cm[upper.tri(cm)]),
             plates = paste(sort(unique(meta$plate[meta$genotype == g])), collapse = ","))
}))
write.csv(rep_cor, file.path(out_dir, "replicate_concordance.csv"), row.names = FALSE)
set.seed(1); rp <- replicate(500, { s <- sample(colnames(logcpm), 2); cor(logcpm[, s[1]], logcpm[, s[2]]) })
cat(sprintf("Replicate wells r = %.3f (median); random different-genotype pairs r = %.3f\n",
            median(rep_cor$mean_r), median(rp)))
cat("Lowest-concordance replicates (possible sample swaps):\n")
print(head(rep_cor[order(rep_cor$mean_r), ], 8), row.names = FALSE)

# ---- 4. class contrasts: RIL vs B73, NIL vs B73, RIL vs NIL -------------------------
# Wells of the same genotype are not independent, so genotype is a blocking factor
# (duplicateCorrelation). Interpretation:
#   RIL vs B73 : average effect of ~50% PT genome (genes with strong PT alleles shift the mean)
#   NIL vs B73 : average effect of small PT segments in B73 background
#   RIL vs NIL : genes responding to genome-wide PT ancestry beyond the NIL segments
dc <- model.matrix(~ 0 + genotype_class + plate, data = meta)
colnames(dc) <- sub("^genotype_class", "", colnames(dc))
# duplicateCorrelation is slow (~0.1 s/gene), so the consensus correlation is
# estimated on a random 2,000-gene subset; the fit itself uses all genes.
set.seed(42)
dc_idx <- sample(nrow(y), min(2000, nrow(y)))
v  <- voom(y, dc)
rho <- duplicateCorrelation(v[dc_idx, ], dc, block = meta$genotype)$consensus.correlation
v  <- voom(y, dc, block = meta$genotype, correlation = rho)
rho <- duplicateCorrelation(v[dc_idx, ], dc, block = meta$genotype)$consensus.correlation
cat(sprintf("Within-genotype correlation (duplicateCorrelation): %.3f\n", rho))
fc <- lmFit(v, dc, block = meta$genotype, correlation = rho)
cm <- makeContrasts(RIL_vs_B73 = ril - b73, NIL_vs_B73 = nil - b73, RIL_vs_NIL = ril - nil, levels = dc)
fc <- eBayes(contrasts.fit(fc, cm))

class_res <- do.call(rbind, lapply(colnames(cm), function(cn) {
  tt <- topTable(fc, coef = cn, number = Inf, sort.by = "none")
  data.frame(gene_id = rownames(tt), contrast = cn, logFC = tt$logFC, AveExpr = tt$AveExpr,
             t = tt$t, PValue = tt$P.Value, FDR = tt$adj.P.Val)
}))
write.csv(class_res, file.path(out_dir, "class_contrasts_all.csv"), row.names = FALSE)
class_sum <- do.call(rbind, lapply(split(class_res, class_res$contrast), function(d) data.frame(
  contrast = d$contrast[1], up = sum(d$FDR < fdr_cut & d$logFC > 0),
  down = sum(d$FDR < fdr_cut & d$logFC < 0))))
cat("Class contrasts, genes at FDR <", fdr_cut, ":\n"); print(class_sum, row.names = FALSE)
write.csv(class_sum, file.path(out_dir, "class_contrasts_summary.csv"), row.names = FALSE)

class_res$candidate <- if (!is.null(cand)) class_res$gene_id %in% cand$gene_id else FALSE
p_vol <- ggplot(class_res, aes(logFC, -log10(FDR))) +
  geom_point(data = subset(class_res, !candidate), colour = "grey65", size = 0.6, alpha = 0.6) +
  geom_point(data = subset(class_res, candidate), colour = "#d7301f", size = 1.6) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", colour = "grey40") +
  facet_wrap(~ contrast) + theme_bw() +
  labs(x = "log2 fold change", y = "-log10(FDR)",
       title = "Class contrasts (plate adjusted; candidate genes in red)")
ggsave(file.path(out_dir, "volcano_class_contrasts.png"), p_vol, width = 11, height = 4.5, dpi = 300)

# ---- 5. per-line screens: each NIL / RIL vs B73 -----------------------------------
# One limma-voom fit with every line as its own level (B73 = reference) + plate.
# Most lines are 1-2 wells, so variance is borrowed across genes and from the
# replicated B73 checks: treat these as screens for large effects, not final tests.
lines_in <- if (RUN_PER_RIL_SCREEN) c("b73", "nil", "ril") else c("b73", "nil")
ys <- y[, meta$genotype_class %in% lines_in]
ys$samples$line <- relevel(factor(ys$samples$line), ref = "B73")
ds <- model.matrix(~ line + plate, data = ys$samples)
ds <- ds[, !colnames(ds) %in% nonEstimable(ds), drop = FALSE]
cat(sprintf("Per-line screen: %d lines vs B73, %d samples, %d residual df\n",
            sum(grepl("^line", colnames(ds))), ncol(ys), ncol(ys) - ncol(ds)))
vs <- voom(ys, ds)
# single-well lines are noisy: robust variance moderation + a fold-change floor
# (treat, |log2FC| > 1) keep one odd well from producing many calls
fs <- treat(lmFit(vs, ds), lfc = 1, robust = TRUE)
line_coefs <- grep("^line", colnames(ds), value = TRUE)
class_of <- tapply(as.character(ys$samples$genotype_class), as.character(ys$samples$line), `[`, 1)

line_res <- do.call(rbind, lapply(line_coefs, function(cf) {
  tt <- topTreat(fs, coef = cf, number = Inf, sort.by = "none")
  ln <- sub("^line", "", cf)
  data.frame(gene_id = rownames(tt), line = ln, class = class_of[[ln]],
             logFC = tt$logFC, AveExpr = tt$AveExpr, PValue = tt$P.Value, FDR = tt$adj.P.Val,
             row.names = NULL)
}))
for (cl in intersect(c("nil", "ril"), unique(line_res$class))) {
  d <- line_res[line_res$class == cl, ]
  saveRDS(d, file.path(out_dir, sprintf("%s_vs_b73_per_line.rds", cl)))   # large: rds, not csv
  s <- aggregate(FDR ~ line, data = d, function(f) sum(f < fdr_cut))
  names(s)[2] <- "n_DE"
  s <- s[order(-s$n_DE), ]
  write.csv(s, file.path(out_dir, sprintf("%s_vs_b73_per_line_summary.csv", cl)), row.names = FALSE)
  cat(sprintf("\n%s lines with the most DE genes vs B73 (|log2FC| > 1, FDR < %.2f):\n", toupper(cl), fdr_cut))
  print(head(s, 10), row.names = FALSE)
}
if (!is.null(cand)) write.csv(line_res[line_res$gene_id %in% cand$gene_id, ],
                              file.path(out_dir, "candidate_per_line_vs_b73.csv"), row.names = FALSE)

# ---- 6. candidate genes ------------------------------------------------------------
if (!is.null(cand)) {
  cand_in <- cand[cand$gene_id %in% rownames(logcpm), ]
  cand_long <- do.call(rbind, lapply(seq_len(nrow(cand_in)), function(i)
    data.frame(meta[, c("sample_id", "genotype", "genotype_class", "plate")],
               gene = cand_in$gene_name[i], logCPM = logcpm[cand_in$gene_id[i], meta$sample_id])))
  p_cand <- ggplot(cand_long, aes(genotype_class, logCPM, fill = genotype_class)) +
    geom_boxplot(outlier.size = 0.6) + facet_wrap(~ gene, scales = "free_y") +
    theme_bw() + guides(fill = "none") + labs(x = NULL, title = "Candidate genes by genotype class")
  ggsave(file.path(out_dir, "candidate_genes_by_class.png"), p_cand, width = 10, height = 7, dpi = 300)
  write.csv(merge(cand_in, class_res, by = "gene_id"),
            file.path(out_dir, "candidate_class_contrasts.csv"), row.names = FALSE)
}

# Allele test (cis): plate-adjusted logCPM averaged per RIL ~ allele at nearest marker
if (!is.null(cand) && file.exists(GENO_FILE) && file.exists(gtf_path)) {
  gtf <- read.delim(gtf_path, header = FALSE, comment.char = "#", quote = "",
                    colClasses = c("character", "NULL", "character", "numeric", "numeric",
                                   "NULL", "NULL", "NULL", "character"))
  names(gtf) <- c("chr", "type", "start", "end", "attr")
  gtf <- gtf[gtf$type == "gene", ]
  gtf$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$attr)
  geno <- read.csv(GENO_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  ril  <- meta[meta$genotype_class == "ril" & meta$genotype %in% names(geno), ]
  ril_adj <- removeBatchEffect(logcpm[, ril$sample_id], batch = droplevels(ril$plate))
  cand_in <- cand[cand$gene_id %in% rownames(logcpm), ]

  allele_res <- do.call(rbind, lapply(seq_len(nrow(cand_in)), function(i) {
    g <- gtf[gtf$gene_id == cand_in$gene_id[i], ][1, ]
    if (is.na(g$chr)) return(NULL)
    mk <- geno[geno$chr == g$chr, ]
    if (!nrow(mk)) return(NULL)
    mk <- mk[which.min(abs(mk$pos - (g$start + g$end) / 2)), ]
    a  <- as.character(unlist(mk[1, ril$genotype]))
    a  <- factor(ifelse(a %in% c("0", "B73"), "B73", ifelse(a %in% c("2", "PT"), "PT", NA)),
                 levels = c("B73", "PT"))
    d  <- data.frame(logCPM = ril_adj[cand_in$gene_id[i], ril$sample_id], allele = a, genotype = ril$genotype)
    d  <- aggregate(logCPM ~ genotype + allele, data = d[!is.na(d$allele), ], FUN = mean)
    if (length(unique(d$allele)) < 2) return(NULL)
    co <- summary(lm(logCPM ~ allele, data = d))$coefficients["allelePT", ]
    data.frame(gene_id = cand_in$gene_id[i], gene_name = cand_in$gene_name[i], marker = mk$marker,
               marker_dist_bp = abs(mk$pos - (g$start + g$end) / 2),
               n_B73 = sum(d$allele == "B73"), n_PT = sum(d$allele == "PT"),
               PT_minus_B73_log2 = co[1], p = co[4])
  }))
  if (!is.null(allele_res)) {
    allele_res$FDR <- p.adjust(allele_res$p, "BH")
    write.csv(allele_res[order(allele_res$p), ], file.path(out_dir, "candidate_allele_tests.csv"), row.names = FALSE)
    cat("\nCandidate allele tests:\n"); print(allele_res[order(allele_res$p), ], row.names = FALSE)
  }
} else {
  cat("\nSkipping allele tests (optional): needs", cand_file, ",", GENO_FILE, "and", gtf_path, "\n")
}

# ---- 7. eQTL-ready export ---------------------------------------------------------
# One value per genotype (mean of its wells) after removing plate / run / group effects.
cov_design <- model.matrix(~ plate + run + experimental_group, data = meta)
cov_design <- cov_design[, !colnames(cov_design) %in% nonEstimable(cov_design), drop = FALSE]
resid <- removeBatchEffect(logcpm, covariates = cov_design[, -1, drop = FALSE])
geno_mean <- sapply(split(meta$sample_id, meta$genotype), function(s) rowMeans(resid[, s, drop = FALSE]))
saveRDS(geno_mean, file.path(out_dir, "expr_per_genotype_adjusted.rds"))
write.csv(data.frame(genotype = colnames(geno_mean),
                     class = meta$genotype_class[match(colnames(geno_mean), meta$genotype)],
                     n_wells = as.integer(table(meta$genotype)[colnames(geno_mean)])),
          file.path(out_dir, "expr_per_genotype_samples.csv"), row.names = FALSE)
cat("eQTL matrix:", nrow(geno_mean), "genes x", ncol(geno_mean), "genotypes\n\nDone.\n")
