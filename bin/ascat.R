#!/usr/bin/env Rscript
# ASCAT v3 on one tumour/normal WGS pair, in two steps so a refit does not redo the allele counting:
#
#   ascat.R prepare --tumor= --normal= --tumor_cram= --normal_cram= --sex=XX|XY --ref_fasta=
#                   --loci_prefix= --alleles_prefix= --threads=
#       -> <tumor>_{tumour,normal}{LogR,BAF}.txt           (alleleCounter; the slow step)
#
#   ascat.R fit     --tumor= --sex=XX|XY --gc= --rt= [--min_ploidy=1.5 --max_ploidy=5.5]
#                   [--purity= --ploidy=]                   (both set: manual fit)
#                   [--allele_refit=true]                   (the both-allele check)
#                   [--allele_cap=0.15]                     (cap on each allele's distance, in copies)
#                   [--min_aberrant=0.005]                  (refit guard)
#       -> <tumor>.ascat.{purity_ploidy,segments,segments_raw,metrics}.tsv, plots, .RData
#
# The fit follows the ASCAT HTS vignette: loadData -> correctLogR (GC + replication timing) -> aspcf ->
# runAscat(gamma = 1). Then:
#   - refit guard: purity can only be estimated from aberrant segments. At ASCAT's fit, the fraction of autosomal
#     het SNPs in segments > 0.15 copies from 1+1 is computed; below --min_aberrant (a flat genome) no refit is
#     tried and purity_status is "not identifiable" (the summary then reports shares of all cells);
#   - the both-allele check (allele_fit, fitlib.R): scores both alleles instead of ASCAT's minor allele, and
#     replaces ASCAT's fit when refit_decision says so (the refit is runAscat at the chosen purity/psi);
#   - QC of the pair: sex from the normal's X heterozygosity, tumour cells in the normal (fitlib.R).
# purity_ploidy.tsv holds the fit used, ASCAT's own and the both-allele fit, fit_used, purity_floor, the aberrant
# fraction, purity_status and the QC columns (sex_data, sex_het_ratio, sex_check, tin_fraction, tin_se,
# tin_n_het). segments_raw has n_het (het SNPs per segment). After a refit, ASCAT's own fit is also written as
# <tumor>.ascat.ownfit.{purity_ploidy,segments,segments_raw}.tsv, for the summary to fall back on when somatic
# SNVs do not confirm the refit (snv_purity.py).
script_dir <- function() dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))

# Het SNPs per segment (ASCAT's segmented BAF is defined at the germline het SNPs only).
het_counts <- function(bc, segs) {
    if (is.null(segs) || nrow(segs) == 0) return(numeric())
    ids <- rownames(bc$Tumor_BAF_segmented[[1]])
    chr <- as.character(bc$SNPpos[ids, 1]); pos <- as.numeric(bc$SNPpos[ids, 2])
    vapply(seq_len(nrow(segs)), function(i) {
        k <- chr == as.character(segs$chr[i])
        as.numeric(sum(pos[k] >= segs$startpos[i] & pos[k] <= segs$endpos[i]))
    }, 0)
}

# Fraction of autosomal het SNPs in segments whose unrounded copies are > dev from 1+1 (fitlib.R aberrant_share).
aberrant_fraction <- function(raw, n_het, dev = 0.15) {
    if (is.null(raw) || nrow(raw) == 0) return(0)
    aberrant_share(pmax(raw$nAraw, raw$nBraw), pmin(raw$nAraw, raw$nBraw), n_het,
                   !(as.character(raw$chr) %in% c("X", "Y")), dev)
}

# Sex from the normal's heterozygosity on X outside the PARs vs the autosomes (fitlib.R sex_call).
sex_check <- function(bc, given, xx_min = SEX_XX_MIN, xy_max = SEX_XY_MAX) {
    baf <- bc$Germline_BAF[, 1]; chr <- as.character(bc$SNPpos[, 1]); pos <- as.numeric(bc$SNPpos[, 2])
    ok  <- !is.na(baf)
    het <- ok & baf > 0.2 & baf < 0.8
    x   <- ok & chr == "X" & pos > PAR_X[1] & pos < PAR_X[2]
    au  <- ok & chr %in% as.character(1:22)
    if (sum(x) < 1000 || sum(au) == 0) return(SEX_NOT_ASSESSED)
    sex_call((sum(het & x) / sum(x)) / (sum(het & au) / sum(au)), given, xx_min, xy_max)
}

# Tumour cells in the normal (fitlib.R tin_core) on ASCAT's het SNPs and unrounded segments. This needs
# prepareHTS to orient the A/B alleles alike in both samples. NA without imbalanced segments or a purity.
tin_estimate <- function(bc, raw, purity, min_het = 200, min_dev = 0.12) {
    if (is.null(raw) || nrow(raw) == 0 || is.null(purity) || is.na(purity)) return(list(tin = NA, se = NA, n_het = 0))
    ids <- rownames(bc$Tumor_BAF_segmented[[1]])
    chr <- as.character(bc$SNPpos[ids, 1]); pos <- as.numeric(bc$SNPpos[ids, 2])
    tb  <- bc$Tumor_BAF[ids, 1]
    s   <- (bc$Germline_BAF[ids, 1] - 0.5) * sign(tb - 0.5)
    a <- pmax(raw$nAraw, raw$nBraw); b <- pmax(pmin(raw$nAraw, raw$nBraw), 0)
    ok  <- as.character(raw$chr) %in% as.character(1:22) & raw$endpos - raw$startpos + 1 >= 1e6
    in_seg <- function(i) chr == as.character(raw$chr[i]) & pos >= raw$startpos[i] & pos <= raw$endpos[i] & !is.na(s)
    tin_core(s, tb, in_seg, a, b, purity, ok, min_het, min_dev)
}

# The both-allele check on the segmented data runAscat scores (the runs of identical segmented LogR/BAF, split at
# chromosome ends; autosomes only, as the formula assumes two normal copies). At purity rho and psi, per run:
#   k = 2^(logR / gamma) ((1 - rho) 2 + rho psi),  minor = (rho - 1 - (BAF - 1) k) / rho,  major = (rho - 1 + BAF k) / rho
# Distance: both alleles' capped squared distance to an integer (fitlib.R allele_off), weighted by het SNPs.
# Grid and local minima (7x7 window) as in runAscat, infeasible points excluded; minima outside the ploidy bounds
# are dropped; prefer_undoubled picks the best. Confident runs: >= 200 het SNPs and >= 1 Mb.
allele_fit <- function(bc, min_ploidy, max_ploidy, gamma = 1, cap = 0.15, doubling_ratio = 1.5, tol = 0.1) {
    bs  <- bc$Tumor_BAF_segmented[[1]]
    ids <- rownames(bs)
    m   <- data.frame(r = bc$Tumor_LogR_segmented[ids, 1], b = bs[, 1], chr = as.character(bc$SNPpos[ids, 1]),
                      pos = as.numeric(bc$SNPpos[ids, 2]), stringsAsFactors = FALSE)
    m   <- m[!is.na(m$r) & !is.na(m$b), ]
    new <- c(TRUE, diff(m$r) != 0 | diff(m$b) != 0 | m$chr[-1] != m$chr[-nrow(m)])
    run <- cumsum(new)
    s   <- data.frame(r = m$r[new], b = pmax(m$b[new], 1 - m$b[new]), n = as.numeric(table(run)), chr = m$chr[new],
                      start = as.numeric(tapply(m$pos, run, min)), end = as.numeric(tapply(m$pos, run, max)))
    s   <- s[s$chr %in% as.character(1:22), ]
    conf <- s$n >= 200 & s$end - s$start + 1 >= 1e6
    copies <- function(rho, psi) {
        k <- 2^(s$r / gamma) * ((1 - rho) * 2 + rho * psi)
        list(minor = (rho - 1 - (s$b - 1) * k) / rho, major = (rho - 1 + s$b * k) / rho)
    }
    dist <- function(rho, psi) {
        cn <- copies(rho, psi)
        sum((allele_off(cn$minor, cap) + allele_off(cn$major, cap)) * s$n)
    }
    feasible <- function(rho, psi) {
        if (!any(conf)) return(TRUE)
        cn <- copies(rho, psi)
        alleles_feasible(cn$major, cn$minor, conf, tol)
    }
    psi_pos <- seq(1, 5.4, 0.05); rho_pos <- seq(0.1, 1, 0.01)
    floor_at <- function(psi) {
        ok <- vapply(rho_pos, function(r) feasible(r, psi), TRUE)
        if (any(ok)) rho_pos[which(ok)[1]] else 1
    }
    d <- outer(psi_pos, rho_pos, Vectorize(function(psi, rho) dist(rho, psi)))
    f <- outer(psi_pos, rho_pos, Vectorize(function(psi, rho) feasible(rho, psi)))
    dm <- ifelse(f, d, Inf)
    minima <- NULL
    for (i in seq_along(psi_pos)) for (j in seq_along(rho_pos)) {
        if (!is.finite(dm[i, j])) next
        w <- dm[max(1, i - 3):min(length(psi_pos), i + 3), max(1, j - 3):min(length(rho_pos), j + 3)]
        if (dm[i, j] <= min(w)) {
            cn <- copies(rho_pos[j], psi_pos[i])
            ploidy <- sum((cn$minor + cn$major) * s$n) / sum(s$n)
            if (ploidy > min_ploidy && ploidy < max_ploidy)
                minima <- rbind(minima, data.frame(dist = d[i, j], purity = rho_pos[j], psi = psi_pos[i], ploidy = ploidy))
        }
    }
    if (is.null(minima)) return(list(best = NULL, dist = dist, feasible = feasible, floor_at = floor_at))
    minima <- minima[order(minima$dist), ]
    c(prefer_undoubled(minima, doubling_ratio), list(dist = dist, feasible = feasible, floor_at = floor_at))
}

tsv_out <- function(x, file) write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE)

prepare <- function(get, tumor, sex, files) {
    # CRAM decoding needs the reference (no internet for the @SQ UR/M5 lookup): alleleCounter's -r is passed via
    # additional_allelecounter_flags (ASCAT 3.2.0 has no ref.fasta argument)
    ref_fasta <- normalizePath(get("ref_fasta"), mustWork = TRUE)
    hts_args <- names(formals(ascat.prepareHTS))
    ref_arg <- if ("additional_allelecounter_flags" %in% hts_args) {
        list(additional_allelecounter_flags = paste("-r", ref_fasta))
    } else if ("ref.fasta" %in% hts_args) {
        list(ref.fasta = ref_fasta)
    } else {
        stop("this ascat.prepareHTS cannot pass a reference to alleleCounter; CRAM input will not decode")
    }
    do.call(ascat.prepareHTS, c(list(
        tumourseqfile = get("tumor_cram"), normalseqfile = get("normal_cram"),
        tumourname = tumor, normalname = get("normal"), allelecounter_exe = Sys.which("alleleCounter"),
        alleles.prefix = get("alleles_prefix"), loci.prefix = get("loci_prefix"), gender = sex,
        genomeVersion = "hg38", nthreads = as.integer(get("threads")), chrom_names = c(1:22, "X"),
        tumourLogR_file = files[1], tumourBAF_file = files[2], normalLogR_file = files[3], normalBAF_file = files[4]
    ), ref_arg))
}

# The both-allele check and, if refit_decision says so, the refit. -> list(out, refit, af, purity_floor, note)
check_and_refit <- function(bc, out, run, get, min_ploidy, max_ploidy) {
    af <- allele_fit(bc, min_ploidy, max_ploidy, cap = as.numeric(get("allele_cap", "0.15")))
    a_rho <- out$purity[[1]]; a_psi <- if (!is.null(out$psi)) out$psi[[1]] else out$ploidy[[1]]
    floor <- af$floor_at(a_psi)
    b <- af$best
    if (is.null(b)) return(list(out = out, refit = FALSE, af = af, purity_floor = floor, note = sprintf(
        "allele check: no feasible minimum within the ploidy bounds (purity floor %.2f); ASCAT fit kept", floor)))
    dec <- refit_decision(a_rho, out$ploidy[[1]], af$dist(a_rho, a_psi), af$feasible(a_rho, a_psi), b, af$runner_up)
    vs_ascat  <- af$dist(a_rho, a_psi) / b$dist
    vs_runner <- if (is.null(af$runner_up)) Inf else af$runner_up$dist / b$dist
    msg <- sprintf("allele check: ASCAT %.2f/%.2f, both-allele %.2f/%.2f (ASCAT %.2fx worse, runner-up %.2fx; purity floor %.2f): %s",
                   a_rho, out$ploidy[[1]], b$purity, b$ploidy, vs_ascat, vs_runner, floor, dec$why)
    message(msg)
    if (!dec$refit) return(list(out = out, refit = FALSE, af = af, purity_floor = floor, note = paste0(msg, "; ASCAT fit kept")))
    new <- run(b$purity, b$psi)
    if (is.null(new$segments) || nrow(new$segments) == 0)
        return(list(out = out, refit = FALSE, af = af, purity_floor = floor, note = paste0(msg, "; refit failed, ASCAT fit kept")))
    list(out = new, refit = TRUE, af = af, purity_floor = floor, note = paste0(msg, "; refit used"))
}

fit <- function(get, tumor, sex, files) {
    purity <- get("purity", ""); ploidy <- get("ploidy", "")
    if (xor(nzchar(purity), nzchar(ploidy))) stop("set both --purity and --ploidy, or neither")
    manual <- nzchar(purity)
    # fail before the slow correction/segmentation if this ASCAT lacks an argument we pass
    missing_args <- setdiff(c("gamma", "write_segments", "min_ploidy", "max_ploidy", "rho_manual", "psi_manual"),
                            names(formals(ascat.runAscat)))
    if (length(missing_args)) stop("ascat.runAscat in ASCAT ", as.character(packageVersion("ASCAT")),
                                   " has no argument(s): ", paste(missing_args, collapse = ", "))

    bc <- ascat.loadData(Tumor_LogR_file = files[1], Tumor_BAF_file = files[2],
                         Germline_LogR_file = files[3], Germline_BAF_file = files[4], gender = sex, genomeVersion = "hg38")
    sx <- sex_check(bc, sex)
    message(sprintf("sex check: given %s, data %s (X/autosome het ratio %.3f)", sex, sx$data, sx$ratio))
    ascat.plotRawData(bc, img.prefix = "before_correction_")
    bc <- ascat.correctLogR(bc, GCcontentfile = get("gc"), replictimingfile = get("rt"))
    ascat.plotRawData(bc, img.prefix = "after_correction_")
    bc <- ascat.aspcf(bc)
    ascat.plotSegmentedData(bc)
    min_ploidy <- as.numeric(get("min_ploidy", "1.5")); max_ploidy <- as.numeric(get("max_ploidy", "5.5"))
    run <- function(rho, psi) ascat.runAscat(bc, gamma = 1, write_segments = TRUE, min_ploidy = min_ploidy,
                                             max_ploidy = max_ploidy, rho_manual = rho, psi_manual = psi)
    out <- ascat_fit <- run(if (manual) as.numeric(purity) else NA, if (manual) as.numeric(ploidy) else NA)

    # refit guard: enough aberrant genome to estimate purity at all?
    min_aberrant <- as.numeric(get("min_aberrant", "0.005"))
    has_fit <- !is.null(out$segments_raw) && nrow(out$segments_raw) > 0
    ab_frac <- if (has_fit) aberrant_fraction(out$segments_raw, het_counts(bc, out$segments_raw)) else 0
    identifiable <- ab_frac >= min_aberrant
    message(sprintf("aberrant het-SNP fraction at ASCAT's fit: %.4f (guard %.3f)", ab_frac, min_aberrant))

    chk <- list(out = out, refit = FALSE, af = NULL, purity_floor = NA, note = if (manual) "manual purity/ploidy" else "")
    if (!manual && !identifiable) {
        chk$note <- sprintf("purity not identifiable: %.1f%% of het SNPs in aberrant segments (< %.1f%%); ASCAT fit kept, no refit",
                            100 * ab_frac, 100 * min_aberrant)
    } else if (!manual && is_true(get("allele_refit", "true")) && !is.null(out$segments) && nrow(out$segments) > 0) {
        chk <- check_and_refit(bc, out, run, get, min_ploidy, max_ploidy)
    }
    out <- chk$out; af <- chk$af
    qc <- tryCatch(ascat.metrics(bc, out), error = function(e) data.frame(error = conditionMessage(e)))

    if (is.null(out$segments) || nrow(out$segments) == 0) {
        # ASCAT gives no fit for a non-aberrant sample: recorded, not an error
        pp <- data.frame(sample = tumor, purity = NA, ploidy = NA, goodness_of_fit = NA,
                         note = "no ASCAT fit (non-aberrant or low purity?)",
                         ascat_purity = NA, ascat_ploidy = NA, allele_purity = NA, allele_ploidy = NA, fit_used = "none",
                         aberrant_het_fraction = 0, purity_status = "not identifiable")
        segs <- data.frame(sample = character(), chr = character(), startpos = integer(),
                           endpos = integer(), nMajor = integer(), nMinor = integer())
        raw <- NULL
        tsv_out(segs, paste0(tumor, ".ascat.segments_raw.tsv"))
    } else {
        pp <- data.frame(sample = tumor, purity = out$purity[[1]], ploidy = out$ploidy[[1]],
                         goodness_of_fit = out$goodnessOfFit[[1]], note = chk$note,
                         ascat_purity = ascat_fit$purity[[1]], ascat_ploidy = ascat_fit$ploidy[[1]],
                         allele_purity = if (is.null(af$best)) NA else af$best$purity,
                         allele_ploidy = if (is.null(af$best)) NA else af$best$ploidy,
                         purity_floor = chk$purity_floor,
                         fit_used = if (manual) "manual" else if (chk$refit) "allele" else "ascat",
                         aberrant_het_fraction = round(ab_frac, 4),
                         purity_status = if (manual) "manual" else if (identifiable) "ok" else "not identifiable")
        segs <- out$segments
        raw <- out$segments_raw
        raw$n_het <- het_counts(bc, raw)   # for the summary's germline-homozygosity filter
        tsv_out(raw, paste0(tumor, ".ascat.segments_raw.tsv"))
    }
    tin <- tin_estimate(bc, raw, if (is.null(out$purity)) NA else out$purity[[1]])
    message(sprintf("tumour in normal: %s (se %s, %d het SNPs in imbalanced segments)",
                    format(round(tin$tin, 3)), format(round(tin$se, 3)), as.integer(tin$n_het)))
    pp <- cbind(pp, sex_data = sx$data, sex_het_ratio = round(sx$ratio, 3), sex_check = sx$status,
                tin_fraction = round(tin$tin, 3), tin_se = round(tin$se, 3), tin_n_het = tin$n_het)
    tsv_out(pp, paste0(tumor, ".ascat.purity_ploidy.tsv"))
    tsv_out(segs, paste0(tumor, ".ascat.segments.tsv"))
    tsv_out(qc, paste0(tumor, ".ascat.metrics.tsv"))
    if (chk$refit) {   # ASCAT's own fit too, for the SNV gate's fallback
        own <- pp
        own$purity <- ascat_fit$purity[[1]]; own$ploidy <- ascat_fit$ploidy[[1]]
        own$goodness_of_fit <- ascat_fit$goodnessOfFit[[1]]; own$fit_used <- "ascat"
        own$note <- paste0("ASCAT's own fit (the both-allele refit was not confirmed); ", chk$note)
        own_raw <- ascat_fit$segments_raw
        own_raw$n_het <- het_counts(bc, own_raw)
        tsv_out(own, paste0(tumor, ".ascat.ownfit.purity_ploidy.tsv"))
        tsv_out(ascat_fit$segments, paste0(tumor, ".ascat.ownfit.segments.tsv"))
        tsv_out(own_raw, paste0(tumor, ".ascat.ownfit.segments_raw.tsv"))
    }
    save(bc, out, ascat_fit, qc, file = paste0(tumor, ".ascat.RData"))
}

main <- function() {
    source(file.path(script_dir(), "fitlib.R"))
    suppressPackageStartupMessages(library(ASCAT))
    cmd  <- commandArgs(trailingOnly = TRUE)
    step <- cmd[1]
    get  <- cli_args(cmd[-1])
    tumor <- get("tumor"); sex <- get("sex")
    files <- paste0(tumor, c("_tumourLogR.txt", "_tumourBAF.txt", "_normalLogR.txt", "_normalBAF.txt"))
    message("ASCAT ", as.character(packageVersion("ASCAT")), " step=", step)
    if (identical(step, "prepare")) prepare(get, tumor, sex, files)
    else if (identical(step, "fit")) fit(get, tumor, sex, files)
    else stop("usage: ascat.R prepare|fit --key=value ...")
    message("ASCAT ", step, " done")
}

if (!interactive()) main()
